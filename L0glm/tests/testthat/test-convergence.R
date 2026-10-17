test_that("a frozen estimated dispersion restarts after a later support loss", {
  calls <- 0L
  local_mocked_bindings(
    lmridge = function(X, iteration, ...) {
      value <- 2 + 1 / iteration
      if (ncol(X) == 1L) value else c(value, if (iteration >= 7L) 0 else 1)
    },
    L0glm_pearson_scale = function(...) {
      calls <<- calls + 1L
      if (calls <= 3L) 1 else 4
    })
  fit <- L0glm.fit.internal(matrix(1, 40, 2), rep(1, 40), family = gaussian(),
    start = c(2, 1), dispersion = -1, lambda = 1, maxit = 12L, tol = 0)
  expect_equal(fit$dispersion, 4)
  expect_equal(fit$convergence$dispersion_restarts, 1L)
  expect_equal(fit$convergence$dispersion_updates, 8L)
  expect_equal(fit$convergence$iterations, 12L)
  expect_false(fit$convergence$converged)
  before <- calls
  fixed <- L0glm.fit.internal(matrix(1, 40, 2), rep(1, 40), family = gaussian(),
    start = c(2, 1), dispersion = 1, lambda = 1, maxit = 12L, tol = 0)
  expect_identical(fixed$dispersion, 1)
  expect_identical(calls, before)
  expect_equal(fixed$convergence$dispersion_updates, 0L)
})

test_that("diagnostics distinguish iteration caps and converged coefficient updates", {
  capped <- L0glm.fit(diag(4), 1:4, family = gaussian(), dispersion = 1,
    start = 1, maxit = 1L, parallelize = FALSE, relaxed = FALSE)
  expect_equal(capped$convergence$iterations, 1L)
  expect_false(capped$convergence$converged)
  converged <- L0glm.fit(diag(4), 1:4, family = gaussian(), dispersion = 1,
    start = 1, lambda = 0, maxit = 100L, parallelize = FALSE, relaxed = TRUE)
  expect_true(all(converged$convergence$converged))
  expect_equal(converged$convergence$stage, c("selection", "relaxed"))
  expect_true(all(converged$convergence$iterations < 5L))
})

test_that("native and R diagnostics describe every channel and fitting pass", {
  set.seed(1640)
  X <- matrix(runif(30 * 8), 30, 8)
  y <- matrix(rpois(90, 3), 30, 3)
  for (sparse in c(FALSE, TRUE)) {
    args <- list(X = Matrix::Matrix(X, sparse = sparse), y = y, start = 1,
      lambda = 2, grouped = TRUE, group.maxit = 2L, relaxed = TRUE,
      maxit = 3L, parallelize = FALSE)
    a <- do.call(L0nnpois, args)
    b <- do.call(L0glm.fit, c(args, list(family = poisson(identity),
                                        nonnegative = TRUE, openmp = FALSE)))
    expect_equal(a$convergence, b$convergence)
    expect_equal(nrow(a$convergence), 9L)
    expect_equal(a$convergence$channel, rep(1:3, 3))
    expect_equal(a$convergence$stage, c(rep("selection", 6), rep("relaxed", 3)))
    expect_true(all(a$convergence$iterations <= 3L))
  }
  parallel <- L0nnpois(X, y, start = 1, lambda = 2, grouped = TRUE,
    group.maxit = 2L, relaxed = TRUE, maxit = 3L, parallelize = TRUE, ncores = 2L)
  serial <- L0nnpois(X, y, start = 1, lambda = 2, grouped = TRUE,
    group.maxit = 2L, relaxed = TRUE, maxit = 3L, parallelize = FALSE)
  expect_equal(parallel$convergence, serial$convergence)
})
