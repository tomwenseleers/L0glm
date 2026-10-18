test_that("default Eigen Gram routing follows numerical sparsity and honors overrides", {
  set.seed(1831)
  dense <- matrix(rnorm(80 * 12), 80, 12)
  zero_filled <- matrix(0, 1000, 1000)
  diag(zero_filled) <- 1
  original <- row_augX
  augmentations <- 0L
  local_mocked_bindings(row_augX = function(...) {
    augmentations <<- augmentations + 1L
    original(...)
  })
  for (method in 1:3) {
    for (spec in list(list(X = dense, sparse = FALSE),
                     list(X = Matrix::Matrix(dense, sparse = TRUE), sparse = TRUE),
                     list(X = zero_filled, sparse = TRUE))) {
      X <- spec$X
      y <- seq_len(nrow(X))/nrow(X)
      for (override in list(NULL, TRUE, FALSE)) {
        before <- augmentations
        control <- control.solver.gen("eigen", method, rowaug = override)
        beta <- lmridge(X, y, lambdas = rep(2, ncol(X)), control.solver = control)
        expected <- if (is.null(override)) spec$sparse else override
        expect_equal(augmentations - before, as.integer(expected))
        residual <- as.numeric(crossprod(as.matrix(X), as.numeric(X %*% beta) - y)) + 2*beta
        expect_lt(max(abs(residual)), 1e-6)
      }
    }
  }
})

test_that("dense row weighting preserves values, names and no-op inputs", {
  X <- matrix(seq_len(15), 5, dimnames = list(letters[1:5], LETTERS[1:3]))
  w <- c(.1, 2, .3, 1, 0)
  expect_equal(Xweights(X, w), sweep(X, 1L, sqrt(w), `*`))
  expect_identical(Xweights(X, rep(1, 5)), X)
  expect_identical(Xweights(X, NULL), X)
})

test_that("IRLS resolves the automatic Gram route once per fit", {
  set.seed(1832)
  X <- matrix(rnorm(100 * 12), 100, 12)
  y <- rbinom(100, 1, .4)
  original <- eigen_rowaug_default
  calls <- 0L
  local_mocked_bindings(eigen_rowaug_default = function(X) {
    calls <<- calls + 1L
    original(X)
  })
  fit <- L0glm.fit(X, y, family = binomial(), openmp = FALSE, start = 1,
    normalize = "none", maxit = 5L, relaxed = FALSE)
  expect_gt(max(fit$convergence$iterations), 1L)
  expect_equal(calls, 1L)
})
