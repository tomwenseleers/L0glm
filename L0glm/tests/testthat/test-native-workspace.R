test_that('a native working design preserves public storage and explicit routing', {
  X <- matrix(0,1000L,1000L);diag(X) <- 1
  a <- L0glm_native_design(X,list())
  expect_true(a$rowaug)
  expect_s4_class(a$X,'dgCMatrix')
  expect_identical(as.matrix(a$X),X)
  b <- L0glm_native_design(X,list(rowaug=FALSE))
  expect_false(b$rowaug)
  expect_identical(b$X,X)
  full <- matrix(1,1000L,1000L)
  expect_identical(L0glm_native_design(full,list())$X,full)
})

test_that('selected workspace reuse matches fresh solves with changing supports', {
  set.seed(114)
  X <- matrix(rnorm(120L*12L),120L,12L)
  X[abs(X)<.6] <- 0
  beta <- matrix(0,12L,3L);beta[1:3,] <- c(.8,-.6,.4)
  Y <- X%*%beta+matrix(rnorm(360L,sd=.3),120L,3L)
  for(sparse in c(FALSE,TRUE))for(method in c(1L,2L,3L,4L)) {
    design <- if(sparse)Matrix::Matrix(X,sparse=TRUE)else X
    control <- controlSolverGen('eigen',method)
    a <- l0glmFit(design,Y,lambda='gic',dispersion=.09,ncores=2L,
      local.search=FALSE,method=method,control.solver=control)
    control$.cache.selected <- FALSE
    b <- l0glmFit(design,Y,lambda='gic',dispersion=.09,ncores=2L,
      local.search=FALSE,method=method,control.solver=control)
    expect_identical(coef(a)!=0,coef(b)!=0)
    expect_equal(coef(a),coef(b),tolerance=1e-11)
    expect_equal(a$convergence,b$convergence,tolerance=1e-11)
  }
})

test_that('reused working systems respect zero observation weights, offsets and priors', {
  set.seed(115)
  X <- Matrix::Matrix(matrix(runif(90L*8L),90L,8L),sparse=TRUE)
  Y <- matrix(rpois(180L,as.numeric(X[,1])*4+2),90L,2L)
  weights <- rep(1,90L);weights[1:8] <- 0
  prior <- matrix(.1,8L,2L);start <- matrix(1,8L,2L);start[8,] <- 0
  for(method in 1:2) {
    control <- controlSolverGen('eigen',method)
    args <- list(X=X,y=Y,family=poisson(identity),nonnegative=TRUE,lambda='gic',
      ncores=2L,local.search=FALSE,method=method,control.solver=control,
      weights=weights,offset=1,prior.mean=prior,start=start,always.include=1L)
    a <- do.call(l0glmFit,args)
    args$control.solver$.native <- FALSE
    b <- do.call(l0glmFit,args)
    expect_true(all(coef(a)[8,]==0))
    expect_identical(coef(a)!=0,coef(b)!=0)
    expect_equal(coef(a),coef(b),tolerance=1e-7)
    expect_equal(GIC(a,mu_floor=1e-8),GIC(b,mu_floor=1e-8),tolerance=1e-8)
  }
})
