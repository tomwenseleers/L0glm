test_that('rectangular and wide LSCG enforce the same bounded ridge quadratic', {
  set.seed(2116)
  for(shape in list(c(60L,20L),c(20L,60L)))for(sparse in c(FALSE,TRUE)) {
    n <- shape[1L];p <- shape[2L]
    X <- matrix(rnorm(n*p),n,p);X[,1] <- 1
    y <- rnorm(n);w <- runif(n,.2,2);w[1] <- 0
    d <- exp(seq(-3,2,length.out=p));d[1] <- 0
    prior <- seq(-.2,.2,length.out=p)
    lower <- c(-Inf,rep(0,p-1L));upper <- c(Inf,rep(.5,p-1L))
    G <- crossprod(X,w*X)+diag(d);h <- as.numeric(crossprod(X,w*y))+d*prior
    reference <- solve_bounded_gram(G,h,rep(0,p),lower,upper,batched=TRUE)
    expect_true(reference$success)
    for(precond in c(1L,2L)) {
      answer <- lmRidge(Matrix::Matrix(X,sparse=sparse),y,weights=w,lambdas=d,
        prior.mean=prior,lower=lower,upper=upper,method=4L,
        control.solver=controlSolverGen('eigen',4L,rowaug=TRUE,precond=precond,
          nonnegative.method='constrained'))
      expect_equal(as.numeric(answer),as.numeric(reference$coefficients),tolerance=1e-7)
      gradient <- as.numeric(G%*%answer-h)
      gradient[answer<=lower&gradient>0|answer>=upper&gradient<0] <- 0
      expect_lt(max(abs(gradient))/max(1,max(abs(h))),1e-8)
      expect_true(all(answer>=lower&answer<=upper))
    }
  }
})

test_that('constrained LSCG does not silently clamp when its Gram is disallowed', {
  set.seed(2117)
  X <- matrix(rnorm(50*8),50,8);y <- -X[,1]+rnorm(50)
  control <- controlSolverGen('eigen',4L,nonnegative.method='constrained',max.memory=0)
  expect_error(lmRidge(X,y,lambdas=rep(1,8),nonnegative=TRUE,method=4L,
    control.solver=control),'constrained LSCG polish')
  expect_error(l0glmFit(X,y,family=gaussian(),start=1,dispersion=1,
    nonnegative=TRUE,method=4L,ncores=1L,maxit=3L,local.search=FALSE,
    control.solver=control),'constrained LSCG polish')
})

test_that('native LSCG corrects both augmented and Woodbury proposals', {
  set.seed(2118)
  for(shape in list(c(50L,20L),c(20L,50L)))for(sparse in c(FALSE,TRUE)) {
    n <- shape[1L];p <- shape[2L]
    X <- matrix(rnorm(n*p),n,p);y <- rnorm(n);w <- runif(n,.3,2)
    prior <- seq(0,.2,length.out=p);d <- rep(.5,p)
    G <- crossprod(X,w*X)+diag(d);h <- as.numeric(crossprod(X,w*y))+d*prior
    reference <- solve_bounded_gram(G,h,rep(0,p),rep(0,p),rep(Inf,p),batched=TRUE)
    fit <- l0glmFit(Matrix::Matrix(X,sparse=sparse),y,family=gaussian(),weights=w,
      start=1,dispersion=1,lambda=2,normalize='none',nonnegative=TRUE,method=4L,
      maxit=1L,thresh=0,ncores=1L,prior.mean=prior,local.search=FALSE,
      control.solver=controlSolverGen('eigen',4L,nonnegative.method='constrained'))
    b <- as.numeric(coef(fit));gradient <- as.numeric(G%*%b-h)
    gradient[b<=0&gradient>0] <- 0
    expect_equal(b,as.numeric(reference$coefficients),tolerance=1e-7)
    expect_lt(max(abs(gradient))/max(1,max(abs(h))),1e-8)
    expect_true(all(b>=0))
  }
})
