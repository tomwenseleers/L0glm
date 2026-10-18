test_that("Gaussian Gram caching preserves weighted ridge selection", {
  set.seed(802)
  X <- matrix(rnorm(600 * 18), 600, 18)
  X <- sweep(X, 2L, seq(.4, 3, length.out = 18), `*`)
  y <- 2 + X[, 1] - .7 * X[, 4] + rnorm(600, sd = .4)
  weights <- rep(c(0, .5, 1, 2), length.out = 600)
  offset <- seq(-.2, .2, length.out = 600)
  for (method in 1:2) for (lambda in c(0, 12)) {
    args <- list(X = X, y = y, weights = weights, offset = offset,
      family = gaussian(), intercept = TRUE, lambda = lambda, method = method,
      prior.mean = rep(.1, 19), start = c(0, 1, 0, rep(1, 16)), always.include = 1L,
      maxit = 80L, tol = 1e-10, parallelize = FALSE)
    cached <- do.call(L0glm.fit, args)
    control <- control.solver.gen("eigen", method)
    control$rowaug <- TRUE
    reference <- do.call(L0glm.fit, c(args, list(control.solver = control)))
    expect_equal(cached$coefficients, reference$coefficients, tolerance = 1e-7)
    expect_identical(as.numeric(cached$coefficients)[3], 0)
  }
})

test_that("zero-weight rows do not turn an effectively wide fit into a Gram solve", {
  set.seed(803)
  X <- matrix(rnorm(80 * 10), 80, 10)
  y <- rnorm(80)
  calls <- 0L
  original <- eigen_ridge_dual
  testthat::local_mocked_bindings(eigen_ridge_dual = function(...) {
    calls <<- calls + 1L
    original(...)
  }, .package = "L0glm")
  fit <- L0glm.fit(X, y, weights = c(rep(1, 5), rep(0, 75)),
    family = gaussian(), intercept = TRUE, start = 1, dispersion = 1,
    lambda = 2, maxit = 3L, parallelize = FALSE)
  expect_gt(calls, 0L)
  expect_equal(dim(fit$coefficients), c(11L, 1L))
  expect_true(all(is.finite(fit$coefficients)))
})
