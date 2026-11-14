test_that('equilibrated OSQP preserves weighted ridge units, priors and bounds', {
  skip_if_not_installed('osqp')
  local_mocked_bindings(solve_bounded_gram=function(...)list(success=FALSE))
  X <- diag(c(.001,1000,2,4))
  y <- as.numeric(X %*% c(-3,2,5,10))
  weights <- c(.4,2,1.3,.8)
  lambda <- c(.000002,100000,.3,.8)
  prior <- c(1,-.2,2,3)
  lower <- c(-Inf,0,0,0); upper <- c(Inf,0,4,8)
  expected <- (diag(X)*weights*y+lambda*prior)/(diag(X)^2*weights+lambda)
  expected <- pmin(upper,pmax(lower,expected))
  for(sparse in c(FALSE,TRUE)) {
    A <- if(sparse)Matrix::Matrix(X,sparse=TRUE)else X
    fit <- lmRidge(A,y,lambdas=lambda,weights=weights,prior.mean=prior,
      lower=lower,upper=upper,start=c(-1,0,2,3),solver='osqp')
    expect_equal(fit,expected,tolerance=1e-6)
    expect_true(all(fit>=lower & fit<=upper))
    expect_identical(fit[2],0)
  }
})

test_that('default OSQP matches polished bigGlm on correlated PSF fits', {
  skip_if_not_installed('osqp'); skip_if_not_installed('glmnet')
  for(fam in c('gaussian','poisson')) {
    n <- if(fam=='gaussian')500L else 1000L
    family <- if(fam=='gaussian')gaussian()else poisson(identity)
    d <- simul(method=if(fam=='gaussian')1L else 2L,n=n,p=n,
      sparse=TRUE,seed=1L,options=list(sd_noise=1))
    for(sparse in c(TRUE,FALSE)) {
      X <- if(sparse)d$X else as.matrix(d$X)
      args <- list(X=X,y=d$y,start=1,family=family,lambda='gic',
        nonnegative=TRUE,control.solver=utils::modifyList(list(),list(.native=FALSE)),ncores=1L,nthreads=1L,
        local.search=FALSE,maxit=25L,tol=1e-10,thresh=1e-5)
      reference <- do.call(l0glmFit,c(args,list(solver='bigGlm')))
      fit <- do.call(l0glmFit,c(args,list(solver='osqp')))
      a <- as.numeric(coef(fit)); b <- as.numeric(coef(reference))
      expect_identical(a!=0,b!=0)
      expect_lt(sqrt(sum((a-b)^2))/max(1,sqrt(sum(b^2))),1e-7)
      ma <- l0glmMetrics(fit,beta_true=d$beta_true)
      mb <- l0glmMetrics(reference,beta_true=d$beta_true)
      expect_equal(ma$gic,mb$gic,tolerance=1e-7)
      expect_equal(ma$R2_coefs,mb$R2_coefs,tolerance=1e-7)
      expect_equal(ma$F1_score,mb$F1_score)
    }
  }
})
