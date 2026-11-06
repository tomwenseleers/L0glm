test_that("historical iteration defaults are retained by both fitting paths", {
  for (fun in list(l0glmFit, l0nnPois, L0glm.fit.internal)) {
    expect_equal(formals(fun)$maxit, 25)
    expect_equal(formals(fun)$tol, 1e-10)
  }
})

test_that("single-channel native fitting does not evaluate unused core counts", {
  X <- diag(3)
  single <- l0nnPois(X, 1:3, ncores = stop("unused core count was forced"))
  serial <- l0nnPois(X, cbind(1:3, 2:4), parallelize = FALSE,
                    ncores = stop("unused core count was forced"))
  expect_s4_class(single$coefficients, "dgCMatrix")
  expect_equal(dim(serial$coefficients), c(3, 2))
})
