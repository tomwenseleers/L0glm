test_that('external clamping solves the unconstrained ridge before projection', {
  H <- matrix(c(1,.8,.8,1),2)
  X <- chol(H);y <- solve(t(X),c(-1,1));lambda <- rep(.1,2)
  expected <- pmax(as.numeric(solve(H+diag(lambda),crossprod(X,y))),0)
  for(solver in c('osqp','bigGlm')) {
    skip_if_not_installed(if(solver=='osqp')'osqp'else'glmnet')
    for(sparse in c(FALSE,TRUE)) {
      Z <- if(sparse)Matrix::Matrix(X,sparse=TRUE)else X
      exact <- lmRidge(Z,y,lambdas=lambda,nonnegative=TRUE,solver=solver)
      clipped <- lmRidge(Z,y,lambdas=lambda,nonnegative=TRUE,solver=solver,
        control.solver=controlSolverGen(solver,nonnegative.method="clamp"))
      expect_equal(clipped,expected,tolerance=1e-8)
      expect_equal(exact,c(0,1/1.1),tolerance=1e-8)
      expect_gt(abs(clipped[2]-exact[2]),1)
    }
  }
})

test_that('external projection preserves weights, priors, scaling and free bounds', {
  set.seed(418);X <- matrix(rnorm(90*7),90,7)
  X[,1] <- 1;y <- -3+X[,2]*2+rnorm(90)
  w <- seq(.3,2,length.out=90);w[1:5] <- 0
  penalty <- seq(.1,.7,length.out=7);prior <- seq(-.3,.3,length.out=7)
  lower <- c(-Inf,rep(0,6));upper <- c(Inf,.4,rep(Inf,5))
  target <- lmRidge(X,y,lambdas=penalty,weights=w,lower=lower,upper=upper,
                    prior.mean=prior,solver='eigen')
  for(solver in c('osqp','bigGlm')) {
    skip_if_not_installed(if(solver=='osqp')'osqp'else'glmnet')
    for(sparse in c(FALSE,TRUE)) {
      fit <- lmRidge(if(sparse)Matrix::Matrix(X,sparse=TRUE)else X,y,
        lambdas=penalty,weights=w,lower=lower,upper=upper,prior.mean=prior,
        scaling=seq(.5,2,length.out=7),solver=solver,
        control.solver=controlSolverGen(solver,nonnegative.method="clamp"))
      expect_equal(fit,target,tolerance=1e-7)
      expect_lt(fit[1],0)
      expect_equal(fit[2],.4,tolerance=1e-7)
    }
  }
})

test_that('native external bridges retain projected fit parity', {
  skip_if_not_installed('osqp');skip_if_not_installed('glmnet')
  for(noise in c('gaussian','poisson'))for(sparse in c(FALSE,TRUE)) {
    sim <- simul(method=if(noise=='gaussian')1L else 2L,n=100L,p=100L,
      s=.1,sparse=sparse,seed=2,options=list(sd_noise=1))
    family <- if(noise=='gaussian')gaussian()else poisson(identity)
    args <- list(X=sim$X,y=sim$y,family=family,nonnegative=TRUE,local.search=FALSE,
      ncores=1L,nthreads=1L,dispersion=1,start=1,maxit=50L)
    reference <- do.call(l0glmFit,c(args,list(solver='eigen',method=1L)))
    for(solver in c('osqp','bigGlm')) {
      fit <- do.call(l0glmFit,c(args,list(solver=solver,
        control.solver=controlSolverGen(solver,nonnegative.method="clamp"))))
      expect_equal(as.matrix(coef(fit)),as.matrix(coef(reference)),tolerance=1e-6)
      expect_identical(coef(fit)!=0,coef(reference)!=0)
      expect_equal(GIC(fit,mu_floor=1e-8),GIC(reference,mu_floor=1e-8),tolerance=1e-7)
    }
  }
})

test_that('external clamping retains a signed binomial intercept and starting exclusions', {
  skip_if_not_installed('osqp');skip_if_not_installed('glmnet')
  set.seed(419);X <- matrix(rnorm(120*2),120,2)
  offset <- seq(-.1,.1,length.out=120)
  y <- cbind(rbinom(120,1,plogis(-2+.7*X[,1]+offset)),
             rbinom(120,1,plogis(-1.7+.9*X[,1]+offset)))
  args <- list(X=X,y=y,family=binomial(),offset=offset,intercept=TRUE,
    nonnegative=TRUE,start=matrix(c(0,1,0),3,2),lambda=1,local.search=FALSE,
    ncores=1L,maxit=50L)
  reference <- do.call(l0glmFit,c(args,list(solver='eigen')))
  for(solver in c('osqp','bigGlm')) {
    fit <- do.call(l0glmFit,c(args,list(solver=solver,
      control.solver=controlSolverGen(solver,nonnegative.method="clamp"))))
    expect_equal(dim(coef(fit)),c(3L,2L))
    expect_true(all(coef(fit)[1,]<0))
    expect_true(all(coef(fit)[3,]==0))
    expect_equal(as.matrix(coef(fit)),as.matrix(coef(reference)),tolerance=1e-6)
  }
})
