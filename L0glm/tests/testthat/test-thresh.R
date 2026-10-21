test_that("thresh is the fitting argument and is validated", {
  for (fun in list(L0glm.fit,L0nnpois)) {
    expect_identical(formals(fun)$thresh,1e-5)
    expect_false("delta" %in% names(formals(fun)))
    for (bad in list(-1,NA_real_,Inf,NaN,"1e-20",TRUE,numeric(),c(0,1)))
      expect_error(fun(diag(2),c(1,2),thresh=bad),"thresh must be")
    expect_error(do.call(fun,list(X=diag(2),y=c(1,2),delta=1e-5)),
                 "unused argument")
  }
})

test_that("adaptive weights are unsmoothed by default and honor explicit helper smoothing", {
  beta <- c(-2,-0.1,0,0.1,2)
  for (gamma in c(1,2,3))
    expect_equal(AR.weights(beta,gamma=gamma,delta=0.2),
                 1/(abs(beta)^gamma+0.2^gamma))
  expect_equal(AR.weights(c(0,1e-20)),c(Inf,1e40))
  expect_equal(AR.weights(c(0,2),delta=0),c(Inf,0.25))
  expect_equal(AR.weights(c(0,2),gamma=0),c(1,1))
})

test_that("selection uses a relative cutoff and preserves exact-zero exclusion", {
  for (sparse in c(FALSE,TRUE)) for (engine in c("native","R")) {
    fun <- if(engine=="native")L0nnpois else L0glm.fit
    args <- list(X=Matrix::Matrix(diag(2),sparse=sparse),y=c(1e6,0.02),family=poisson(identity),
                 start=c(1,1),lambda=0,thresh=0.01,
                 grouped=FALSE,relaxed=FALSE,maxit=2,parallelize=FALSE)
    if(engine=="R") args <- c(args,list(openmp=FALSE,nonnegative=TRUE))
    fit <- do.call(fun,args)
    expect_equal(as.vector(fit$coefficients),c(1e6,0),tolerance=1e-10)
    args$start <- c(0,1)
    fit <- do.call(fun,args)
    expect_identical(as.numeric(fit$coefficients[1]),0)
    expect_gt(as.numeric(fit$coefficients[2]),0)

    # Exact-boundary coefficients survive selection and the grouped handoff.
    args$X <- Matrix::Matrix(diag(2),sparse=sparse)
    args$y <- matrix(c(1,0.25),2,2); args$start <- matrix(1,2,2)
    args$thresh <- 0.25; args$grouped <- TRUE; args$group.maxit <- 2L
    fit <- do.call(fun,args)
    expect_equal(as.matrix(fit$coefficients),args$y,tolerance=1e-14)

    # thresh=0 disables only the relative cutoff, not the absolute floor.
    args$y <- c(1e-21,2e-20); args$start <- c(1,1)
    args$thresh <- 0; args$grouped <- FALSE
    fit <- do.call(fun,args)
    expect_identical(as.numeric(fit$coefficients[1]),0)
    expect_gt(as.numeric(fit$coefficients[2]),1e-20)
  }
})

test_that("relaxed fits retain coefficients below both cutoffs in all execution paths", {
  for (sparse in c(FALSE,TRUE)) for (engine in c("native","R")) {
    fun <- if(engine=="native")L0nnpois else L0glm.fit
    # Neither column is an intercept: both contribute to the relative cutoff.
    args <- list(X=Matrix::Matrix(matrix(c(1,0,1,1.01),2),sparse=sparse),
                 y=c(1.05,1),family=poisson(identity),start=c(1,1),lambda=1,thresh=0.1,
                 grouped=FALSE,relaxed=FALSE,maxit=1,parallelize=FALSE)
    if(engine=="R") args <- c(args,list(openmp=FALSE,nonnegative=TRUE))
    selection <- do.call(fun,args)
    expect_gt(min(selection$coefficients),max(selection$coefficients)*args$thresh)
    args$relaxed <- TRUE
    fit <- do.call(fun,args)
    expect_equal(as.vector(fit$coefficients),c(1.05 - 1 / 1.01, 1 / 1.01),tolerance=1e-12)

    # A retained row can have a coefficient below even the numerical floor.
    args$X <- Matrix::Matrix(diag(3),sparse=sparse)
    args$y <- matrix(c(1e-24,3,0,10,2,0),3,2)
    args$start <- matrix(1,3,2); args$lambda <- 0.01
    args$grouped <- TRUE; args$group.maxit <- 1L
    serial <- do.call(fun,args)
    expect_equal(as.matrix(serial$coefficients),args$y,tolerance=1e-12)
    expect_gt(as.numeric(serial$coefficients[1,1]),0)
    expect_lt(as.numeric(serial$coefficients[1,1]),1e-20)
    args$parallelize <- TRUE; args$ncores <- 2L
    parallel <- do.call(fun,args)
    expect_equal(parallel$coefficients,serial$coefficients,tolerance=1e-12)
  }
})

test_that("R and native penalties are unsmoothed and cutoffs respect coefficient signs", {
  args <- list(X=diag(3),y=c(2,3,4),family=poisson(identity),start=c(1,2,3),lambda=4,
               thresh=0.25,maxit=1,
               grouped=FALSE,relaxed=FALSE,parallelize=FALSE)
  native <- do.call(L0nnpois,args)
  general <- do.call(L0glm.fit,c(args,list(nonnegative=TRUE,openmp=FALSE)))
  expect_equal(as.numeric(native$lambdas),1/args$start^2)
  expect_equal(native$lambdas,general$lambdas)
  expect_equal(native$coefficients,general$coefficients,tolerance=1e-12)
  fit <- L0glm.fit(diag(3),c(-0.04,-0.5,0.4),family=gaussian(),
                  nonnegative=FALSE,openmp=FALSE,lambda=0,thresh=0.1,
                  start=rep(1,3),dispersion=1,maxit=1,
                  relaxed=FALSE,parallelize=FALSE)
  expect_equal(as.numeric(fit$coefficients),c(0,-0.5,0.4))
})
