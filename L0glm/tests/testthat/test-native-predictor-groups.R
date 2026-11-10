group_R_reference <- function(args) {
  local_mocked_bindings(L0glm_group_native_pass=function(...)NULL,.package='L0glm')
  do.call(l0glmFit,args)
}

test_that('native predictor groups preserve the structured GLM solution', {
  set.seed(83);n <- 120L;p <- 12L
  X <- matrix(rnorm(n*p),n,p);X[,12] <- 0
  base <- 1+X[,1]-X[,2]+.6*X[,7]-.6*X[,8]
  families <- list(gaussian(),poisson(),binomial(),Gamma(log),inverse.gaussian(log),
    MASS::negative.binomial(4),quasipoisson())
  for(family in families) {
    mu <- if(family$family=='gaussian')base else family$linkinv(base/3)
    y <- switch(tolower(family$family),gaussian=mu+rnorm(n,sd=.3),
      poisson=rpois(n,mu),binomial=rbinom(n,1,mu),gamma=rgamma(n,shape=8,scale=mu/8),
      inverse.gaussian=pmax(.01,mu*exp(rnorm(n,sd=.15))),
      quasipoisson=rpois(n,mu),rnbinom(n,mu=mu,size=4))
    for(sparse in c(FALSE,TRUE))for(grouped in c(FALSE,TRUE)) {
      start <- matrix(1,p+1L,2L);start[4,1] <- 0
      args <- list(X=if(sparse)Matrix::Matrix(X,sparse=TRUE)else X,y=cbind(y,y),
        family=family,group=rep(1:4,each=3),always.include=11L,intercept=TRUE,
        offset=.1,weights=rep(c(0,1,2),length.out=n),prior.mean=.01,start=start,
        dispersion=if(family$family%in%c('gaussian','Gamma'))NULL else 1,
        grouped=grouped,relaxed=TRUE,local.search=FALSE,parallelize=FALSE,maxit=45L)
      reference <- group_R_reference(args);fit <- do.call(l0glmFit,args)
      expect_equal(coef(fit),coef(reference),tolerance=1e-7)
      expect_equal(GIC(fit),GIC(reference),tolerance=1e-9)
      expect_identical(as.matrix(coef(fit)!=0),as.matrix(coef(reference)!=0))
      expect_equal(fit$dispersion,reference$dispersion,tolerance=1e-9)
      expect_equal(fit$convergence,reference$convergence)
      expect_s4_class(coef(fit),'CsparseMatrix')
    }
  }
})

test_that('native group channel parallelism is deterministic and restores threads', {
  previous <- L0glm_set_thread_count(1L)
  withr::defer(L0glm_restore_thread_count(previous))
  set.seed(284);X <- matrix(rnorm(100*16),100,16)
  Y <- matrix(rpois(100*4,exp(.4+.7*X[,1]-.7*X[,2])),100,4)
  args <- list(X=X,y=Y,family=poisson(),group=rep(1:4,each=4),intercept=TRUE,
    relaxed=TRUE,parallelize=TRUE,local.search=FALSE)
  for(grouped in c(FALSE,TRUE)) {
    args$grouped <- grouped
    one <- do.call(l0glmFit,c(args,list(ncores=1L)))
    four <- do.call(l0glmFit,c(args,list(ncores=4L)))
    expect_identical(coef(one),coef(four))
    expect_identical(one$history,four$history)
    expect_identical(one$dispersion,four$dispersion)
  }
  if(requireNamespace('RhpcBLASctl',quietly=TRUE))
    expect_equal(RhpcBLASctl::omp_get_max_threads(),1L)
})

test_that('group routing preserves callbacks and validation checkpoints', {
  set.seed(42);X <- matrix(rnorm(80*6),80,6);y <- .8+X[,1]-X[,2]+rnorm(80)
  args <- list(X=X,y=y,group=rep(1:3,each=2),intercept=TRUE,dispersion=1,
    relaxed=TRUE,local.search=FALSE,parallelize=FALSE)
  native <- do.call(l0glmFit,args)
  local({
    local_mocked_bindings(glm_predictor_group_cpp=function(...)stop('native sentinel'),.package='L0glm')
    expect_error(do.call(l0glmFit,args),'native sentinel')
    custom <- gaussian();custom$variance <- function(mu)rep(1,length(mu))
    fallback <- do.call(l0glmFit,c(args,list(family=custom)))
    expect_equal(coef(native),coef(fallback),tolerance=1e-8)
  })
  validated <- do.call(l0glmFit,c(args,list(validation=list(X=X,y=y,starts='default'))))
  expect_true(is.finite(validated$validation$deviance))
})

test_that('native group relaxed refits handle aliased columns and wide designs', {
  set.seed(133);X <- matrix(rnorm(40*60),40,60);X[,60] <- X[,59]
  y <- X[,1]-X[,2]+rnorm(40)
  for(sparse in c(FALSE,TRUE)) {
    args <- list(X=if(sparse)Matrix::Matrix(X,sparse=TRUE)else X,y=y,
      group=rep(1:30,each=2),intercept=TRUE,start=1,lambda=0,dispersion=1,
      relaxed=TRUE,local.search=FALSE,parallelize=FALSE)
    reference <- group_R_reference(args);fit <- do.call(l0glmFit,args)
    expect_equal(as.numeric(predict(fit)),as.numeric(predict(reference)),tolerance=1e-8)
    expect_equal(GIC(fit),GIC(reference),tolerance=1e-8)
  }
})

test_that('unused initial link values do not warn in noncanonical Gaussian fits', {
  set.seed(51);X <- matrix(rnorm(100*6),100,6);y <- .4+.2*X[,1]+rnorm(100)
  args <- list(X=X,y=y,family=gaussian(log),group=rep(1:3,each=2),intercept=TRUE,
    dispersion=1,start=1,relaxed=TRUE,local.search=FALSE,parallelize=FALSE)
  reference <- group_R_reference(args)
  expect_silent(fit <- do.call(l0glmFit,args))
  expect_equal(coef(fit),coef(reference),tolerance=1e-7)
})
