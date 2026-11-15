test_that('endpoint certificates preserve independently fitted path models', {
  set.seed(420)
  X <- diag(10)
  for(family in list(poisson(identity),gaussian())) {
    y <- if(family$family=='poisson')c(400,100,30,5,1,0,0,0,0,0)else
      c(15,8,3,.1,-.1,.2,0,0,0,0)
    args <- list(X=X,y=y,family=family,dispersion=1,nonnegative=TRUE,
      ncores=1L,maxit=60L,plot=FALSE,printbest=FALSE,nlambda=5L)
    tuned <- do.call(l0glmTune,args)
    expect_true(tuned$endpoint_zero)
    for(i in seq_along(tuned$lambdas)) {
      fit <- do.call(l0glmFit,c(list(X=X,y=y,lambda=tuned$lambdas[i]),tuned$fit_args))
      expect_equal(coef(fit),tuned$coefficients[[i]],tolerance=1e-10)
      expect_equal(GIC(fit,mu_floor=1e-8),tuned$metrics$gic[i],tolerance=1e-10)
    }
  }
})

test_that('certified endpoint probes classify support without completing local moves', {
  X <- diag(8);y <- c(1000,100,30,5,0,0,0,0)
  args <- list(X=X,y=y,family=poisson(identity),nonnegative=TRUE,ncores=1L,
    lambda=500,maxit=50L,local.search=TRUE)
  full <- do.call(l0glmFit,args)
  probe <- do.call(l0glmFit,c(args,list(control.local=list(.endpoint=TRUE))))
  expect_identical(probe$local.search$stopping,'endpoint_nonzero_certificate')
  expect_gt(sum(coef(full)!=0),0)
  expect_gt(sum(coef(probe)!=0),0)
  expect_lte(probe$local.search$evaluations,full$local.search$evaluations)
  expect_true(all(diff(probe$local.search$history$value)<0))
})

test_that('endpoint certificates exclude estimated scales and compulsory terms', {
  X <- cbind(1,diag(8));y <- c(20,10,3,1,0,0,0,0)
  for(family in list(gaussian(),poisson(identity))) {
    args <- list(X=X,y=y,family=family,nonnegative=TRUE,always.include=1L,
      ncores=1L,lambda=50,maxit=50L,local.search=TRUE)
    reference <- do.call(l0glmFit,args)
    probe <- do.call(l0glmFit,c(args,list(control.local=list(.endpoint=TRUE))))
    expect_identical(coef(probe),coef(reference))
    expect_false(identical(probe$local.search$stopping,'endpoint_nonzero_certificate'))
  }
})
