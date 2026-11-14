test_that('bounded Gram polish agrees with an independent quadratic optimizer', {
  set.seed(511)
  X <- matrix(rnorm(80*6),80,6)
  X[,2] <- X[,1] + .03*X[,2]
  y <- as.numeric(X%*%c(6,0,1,-2,0,3)) + rnorm(80)
  lambda <- c(.02,.01,.2,.1,0,.1)
  prior <- c(1,-.2,.4,0,.1,1)
  G <- crossprod(X) + diag(lambda)
  h <- as.numeric(crossprod(X,y)) + lambda*prior
  lower <- c(0,0,-Inf,-Inf,0,.5)
  upper <- c(Inf,Inf,Inf,0,0,4)
  objective <- function(b)sum((y-X%*%b)^2)+sum(lambda*(b-prior)^2)
  gradient <- function(b)2*(as.numeric(G%*%b)-h)
  reference <- optim(pmin(upper,pmax(lower,rep(0,6))),objective,gradient,
    method='L-BFGS-B',lower=lower,upper=upper,
    control=list(factr=1e2,pgtol=1e-10,maxit=10000))
  for(sparse in c(FALSE,TRUE)) {
    A <- if(sparse)as(as(Matrix::Matrix(G,sparse=TRUE),'generalMatrix'),'CsparseMatrix')else G
    ans <- solve_bounded_gram(A,h,pmin(upper,pmax(lower,rep(0,6))),lower,upper)
    expect_true(ans$success)
    expect_lte(ans$error,1e-10)
    expect_equal(objective(as.numeric(ans$coefficients)),reference$value,tolerance=1e-9)
    expect_true(all(ans$coefficients>=lower & ans$coefficients<=upper))
    expect_equal(as.numeric(ans$coefficients)[5],0)
  }
})

test_that('bigGlm solves high dynamic range ridge steps to their box optimum', {
  skip_if_not_installed('glmnet')
  set.seed(913)
  X <- matrix(rnorm(300*12),300,12)
  X[,2] <- X[,1] + .015*X[,2]
  X[,4] <- X[,3] + .02*X[,4]
  truth <- c(10000,0,6000,0,0,500,rep(0,6))
  y <- as.numeric(X%*%truth)+rnorm(300)
  lambda <- rep(1e-6,12)
  G <- crossprod(X)+diag(lambda)
  h <- as.numeric(crossprod(X,y))
  reference <- solve_bounded_gram(G,h,truth,rep(0,12),rep(Inf,12))
  expect_true(reference$success)
  for(sparse in c(FALSE,TRUE)) {
    A <- if(sparse)Matrix::Matrix(X,sparse=TRUE)else X
    fit <- lmRidge(A,y,lambdas=lambda,nonnegative=TRUE,solver='bigGlm')
    expect_equal(fit,as.numeric(reference$coefficients),tolerance=1e-7)
  }
})

test_that('bigGlm preserves the native path and memory use for wide systems', {
  skip_if_not_installed('glmnet')
  local_mocked_bindings(solve_bounded_gram=function(...)stop('Unexpected Gram allocation'))
  set.seed(514)
  X <- matrix(rnorm(20*60),20,60)
  fit <- lmRidge(X,rnorm(20),lambdas=rep(1,60),nonnegative=TRUE,solver='bigGlm')
  expect_length(fit,60)
  expect_true(all(is.finite(fit) & fit>=0))
})

test_that('bigGlm recovers Gaussian PSF coefficients in dense and sparse storage', {
  skip_if_not_installed('glmnet')
  for(n in c(500L,1000L)) {
    d <- simul(method=1L,n=n,p=n,sparse=TRUE,seed=1L,options=list(sd_noise=1))
    lambda <- presetLambda('gic',p=n,n=n)
    start <- init_coefs(d$X,d$y,mu=d$y,family=gaussian(),lambda=lambda)
    coefs <- list()
    for(sparse in c(TRUE,FALSE)) {
      X <- if(sparse)d$X else as.matrix(d$X)
      fit <- l0glmFit(X,d$y,start=start,lambda=lambda,family=gaussian(),
        nonnegative=TRUE,solver='bigGlm',control.solver=utils::modifyList(list(),list(.native=FALSE)),ncores=1L,
        local.search=FALSE,maxit=25L,tol=1e-10,thresh=1e-5)
      met <- l0glmMetrics(fit,beta_true=d$beta_true)
      expect_gt(met$R2_coefs,.997)
      expect_gte(met$F1_score,if(n==1000L).99 else .94)
      coefs[[length(coefs)+1L]]<-as.numeric(coef(fit))
    }
    expect_equal(coefs[[1]],coefs[[2]],tolerance=1e-8)
  }
})
