test_that('ncores is the public native threading control', {
  expect_false(any(c('openmp','parallelize') %in% names(formals(l0glmFit))))
  expect_false('l0nnPois' %in% getNamespaceExports('L0glm'))
  expect_error(l0glmFit(diag(3),1:3,openmp=FALSE),'unused argument')
  expect_error(l0glmFit(diag(3),1:3,parallelize=FALSE),'unused argument')
  expect_false(exists('L0glm_native_workers',asNamespace('L0glm'),inherits=FALSE))
  expect_false(exists('L0glm_worker_init',asNamespace('L0glm'),inherits=FALSE))
})

test_that('family constructors remain supported by automatic routing', {
  X <- cbind(seq_len(30L)/30,cos(seq_len(30L)))
  y <- X[,1]+sin(seq_len(30L))
  a <- l0glmFit(X,y,family=gaussian,local.search=FALSE,ncores=1L)
  b <- l0glmFit(X,y,family=gaussian(),local.search=FALSE,ncores=1L)
  expect_identical(a$coefficients,b$coefficients)
})

test_that('one-thread and OpenMP built-in fits share the same native driver', {
  set.seed(20261114)
  X <- matrix(rnorm(100L*8L,sd=.2),100L,8L)
  beta <- matrix(0,8L,3L);beta[1:2,] <- c(.6,-.4)
  Y <- X %*% beta+matrix(rnorm(300L,sd=.15),100L,3L)
  original <- glm_native_cpp
  seen <- integer()
  local_mocked_bindings(glm_native_cpp=function(...,nthreads) {
    seen <<- c(seen,nthreads)
    original(...,nthreads=nthreads)
  })
  for(sparse in c(FALSE,TRUE)) {
    design <- if(sparse)Matrix::Matrix(X,sparse=TRUE)else X
    a <- l0glmFit(design,Y,lambda='gic',ncores=1L,local.search=FALSE)
    b <- l0glmFit(design,Y,lambda='gic',ncores=2L,local.search=FALSE)
    expect_identical(tail(seen,2L),c(1L,2L))
    for(name in c('coefficients','weights','lambdas','z','dispersion','convergence'))
      expect_identical(a[[name]],b[[name]])
  }
})

test_that('custom callbacks remain serial without starting R workers', {
  set.seed(613)
  X <- matrix(rnorm(70L*5L),70L,5L)
  Y <- cbind(X[,1]+rnorm(70L),-X[,1]+rnorm(70L))
  family <- gaussian()
  processes <- integer()
  family$variance <- function(mu) {
    processes <<- c(processes,Sys.getpid())
    rep(1,length(mu))
  }
  expect_null(L0glm_native_family(family))
  local_mocked_bindings(L0glm_native_model=function(...)stop('native callback route'))
  a <- l0glmFit(X,Y,family=family,dispersion=1,ncores=1L,local.search=FALSE)
  b <- l0glmFit(X,Y,family=family,dispersion=1,ncores=4L,local.search=FALSE)
  expect_identical(a$coefficients,b$coefficients)
  expect_identical(unique(processes),Sys.getpid())
})

test_that('custom Poisson callbacks bypass the specialized backend', {
  set.seed(615)
  X <- matrix(runif(40L*4L),40L,4L)
  Y <- matrix(rpois(80L,3),40L,2L)
  family <- poisson(identity)
  calls <- 0L
  family$variance <- function(mu) {calls <<- calls+1L;2*mu}
  local_mocked_bindings(L0nnpois_dense_cpp=function(...)stop('specialized callback route'),
    glm_native_cpp=function(...)stop('native callback route'))
  fit <- l0glmFit(X,Y,family=family,nonnegative=TRUE,dispersion=1,
    local.search=FALSE,ncores=2L)
  expect_gt(calls,0L)
  expect_true(all(is.finite(coef(fit))))
})

test_that('R solver bridges stay on the main thread', {
  skip_if_not_installed('glmnet')
  set.seed(614)
  X <- matrix(rnorm(80L*6L),80L,6L)
  Y <- cbind(X[,1]+rnorm(80L),-X[,1]+rnorm(80L))
  original <- glm_native_cpp;threads <- NULL
  local_mocked_bindings(glm_native_cpp=function(...,nthreads) {
    threads <<- nthreads;original(...,nthreads=nthreads)
  })
  fit <- l0glmFit(X,Y,solver='bigGlm',ncores=4L,local.search=FALSE)
  expect_identical(threads,1L)
  expect_true(all(is.finite(coef(fit))))
})

test_that('serial callbacks use the core budget for BLAS and restore settings', {
  X <- matrix(seq_len(120L)/120,40L,3L)
  Y <- cbind(X[,1]+1,X[,2]+2)
  family <- gaussian();family$variance <- function(mu)rep(2,length(mu))
  original <- L0glm_set_thread_count;seen <- integer()
  old <- Sys.getenv(L0glm_thread_env_vars(),unset=NA_character_)
  local_mocked_bindings(L0glm_set_thread_count=function(nthreads) {
    if(!is.null(nthreads))seen <<- c(seen,as.integer(nthreads))
    original(nthreads)
  })
  a <- l0glmFit(X,Y,family=family,dispersion=1,ncores=2L,local.search=FALSE)
  expect_true(2L %in% seen)
  expect_identical(Sys.getenv(L0glm_thread_env_vars(),unset=NA_character_),old)
  seen <- integer()
  b <- l0glmFit(X,Y,family=family,dispersion=1,ncores=4L,nthreads=1L,local.search=FALSE)
  expect_true(length(seen)>0L && all(seen==1L))
  expect_equal(a$coefficients,b$coefficients,tolerance=1e-12)
  expect_identical(Sys.getenv(L0glm_thread_env_vars(),unset=NA_character_),old)
  broken <- family;broken$variance <- function(mu)stop('callback failed')
  expect_error(l0glmFit(X,Y,family=broken,ncores=2L,local.search=FALSE),'callback failed')
  expect_identical(Sys.getenv(L0glm_thread_env_vars(),unset=NA_character_),old)
})

test_that('the specialized Poisson backend avoids nested BLAS threading', {
  X <- matrix(runif(60L*6L),60L,6L)
  Y <- matrix(rpois(60L*3L,3),60L,3L)
  original <- L0nnpois_dense_cpp;seen <- NULL
  before <- RhpcBLASctl::blas_get_num_procs()
  local_mocked_bindings(L0nnpois_dense_cpp=function(...) {
    seen <<- RhpcBLASctl::blas_get_num_procs();original(...)
  })
  fit <- l0glmFit(X,Y,family=poisson(identity),nonnegative=TRUE,
    ncores=2L,nthreads=4L,local.search=FALSE)
  expect_equal(seen,1L)
  expect_equal(RhpcBLASctl::blas_get_num_procs(),before)
  expect_true(all(is.finite(coef(fit))))
})
