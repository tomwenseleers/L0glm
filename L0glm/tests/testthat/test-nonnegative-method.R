test_that("nonnegative working methods preserve solver-specific defaults", {
  for (solver in c("eigen", "osqp", "bigGlm", "bigglm", "glmnet")) {
    expected <- if (solver == "eigen") "clamp" else "constrained"
    expect_identical(controlSolverGen(solver)$nonnegative.method, expected)
    expect_identical(controlSolverGen(solver, nonnegative.method=NULL)$nonnegative.method, expected)
    for (mode in c("clamp", "constrained")) {
      control <- controlSolverGen(solver, nonnegative.method=mode)
      expect_identical(control$nonnegative.method, mode)
      expect_false(any(c("nonnegative.polish", "nonnegative.clamp") %in% names(control)))
    }
  }
  for (method in 1:6) {
    expect_identical(controlSolverGen("eigen", method)$nonnegative.method, "clamp")
    expect_identical(controlSolverGen("eigen", method,
      nonnegative.method="constrained")$nonnegative.method, "constrained")
  }
  expect_false(any(c("nonnegative.polish", "nonnegative.clamp") %in% names(formals(controlSolverGen))))
  for (bad in list(NA, NA_character_, 1, TRUE, character(), c("clamp", "constrained"), "other"))
    expect_error(controlSolverGen("eigen", nonnegative.method=bad), 'must be "clamp" or "constrained"')
})

test_that("explicit default strategies give identical weighted ridge solutions", {
  set.seed(6291)
  X <- matrix(rnorm(100*8),100,8);X[,2] <- X[,1]+.02*X[,2]
  y <- as.numeric(X %*% c(2,-3,1,0,0,1,0,0))+rnorm(100)
  for (solver in c("eigen", "osqp", "bigGlm")) {
    skip_if_not_installed(if(solver=="osqp")"osqp"else"glmnet")
    mode <- if(solver=="eigen")"clamp"else"constrained"
    for (sparse in c(FALSE,TRUE)) {
      design <- if(sparse)Matrix::Matrix(X,sparse=TRUE)else X
      args <- list(X=design,y=y,lambdas=seq(.1,.8,.1),weights=seq(.5,2,length.out=100),
        prior.mean=seq(-.4,.3,.1),nonnegative=TRUE,solver=solver)
      implicit <- do.call(lmRidge,args)
      explicit <- do.call(lmRidge,c(args,list(control.solver=controlSolverGen(solver,nonnegative.method=mode))))
      expect_identical(explicit,implicit)
    }
  }
})
