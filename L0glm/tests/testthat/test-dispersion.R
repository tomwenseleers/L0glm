test_that("positive dispersion is fixed and negative dispersion is an exact start", {
  set.seed(201)
  X <- matrix(rnorm(240 * 6), 240)
  y <- as.numeric(X %*% c(2, -1, 0, 0, 0, 0)) + rnorm(240, sd = 3)
  args <- list(X = X, y = y, family = gaussian(), lambda = "gic",
               start = rep(1, 6), parallelize = FALSE)
  fixed <- do.call(L0glm.fit, c(args, list(dispersion = 2, maxit = 1)))
  initial <- do.call(L0glm.fit, c(args, list(dispersion = -2, maxit = 1)))
  expect_equal(initial$coefficients, fixed$coefficients)
  expect_identical(initial$dispersion, 2)
  fixed <- do.call(L0glm.fit, c(args, list(dispersion = 2, maxit = 100, relaxed = TRUE)))
  estimated <- do.call(L0glm.fit, c(args, list(dispersion = -2, maxit = 100, relaxed = TRUE)))
  expect_identical(fixed$dispersion, 2)
  residual <- y - as.numeric(X %*% estimated$coefficients)
  expect_equal(estimated$dispersion,
               sum(residual^2) / (length(y) - sum(estimated$coefficients != 0)))
  expect_gt(estimated$dispersion, 5)
})

test_that("dispersion updates are bounded and wait for residual degrees of freedom", {
  set.seed(202)
  X <- matrix(rnorm(100 * 5), 100)
  y <- rnorm(100, sd = 20)
  args <- list(X = X, y = y, family = gaussian(), start = rep(1, 5),
               dispersion = -0.001, lambda = 1, thresh = 0)
  before <- do.call(L0glm.fit.internal, c(args, list(maxit = 2)))
  after <- do.call(L0glm.fit.internal, c(args, list(maxit = 3)))
  expect_equal(before$dispersion, 0.001)
  expect_equal(after$dispersion, 0.002)
  wide <- L0glm.fit.internal(diag(12), rep(1, 12), family = gaussian(),
    start = rep(1, 12), dispersion = -0.5, lambda = 1e-6, thresh = 0, maxit = 10)
  expect_equal(wide$dispersion, 0.5)
  expect_true(all(is.finite(wide$coefficients)))
})

test_that("automatic ridge scale respects response units and preserves the RNG", {
  set.seed(203)
  X <- matrix(rnorm(180 * 12), 180)
  y <- X[, 1] * 2 + rnorm(180)
  rng <- .Random.seed
  initial <- function(X, y) L0glm_initial_dispersion(
    X, y, y, rep(1, length(y)), rep(0, length(y)), gaussian())
  phi <- initial(X, y)
  expect_identical(.Random.seed, rng)
  expect_equal(initial(X, 10 * y), 100 * phi, tolerance = 1e-8)
  expect_equal(initial(Matrix::Matrix(X, sparse = TRUE), y), phi, tolerance = 1e-8)
  a <- L0glm.fit(X, y, lambda = "gic", maxit = 100, relaxed = TRUE,
                 family = gaussian(), parallelize = FALSE)
  b <- L0glm.fit(X, 10 * y, lambda = "gic", maxit = 100, relaxed = TRUE,
                 family = gaussian(), parallelize = FALSE)
  expect_equal(b$coefficients / 10, a$coefficients, tolerance = 1e-7)
  expect_equal(b$dispersion / 100, a$dispersion, tolerance = 1e-7)
})

test_that("a stable unpenalized fit does not wait for an unavailable scale update", {
  calls <- 0L
  original <- lmridge
  testthat::local_mocked_bindings(lmridge = function(...) {
    calls <<- calls + 1L
    original(...)
  })
  fit <- L0glm.fit.internal(diag(12), rep(1, 12), family = gaussian(),
    start = rep(1, 12), dispersion = -0.5, lambda = 0, maxit = 100)
  expect_equal(fit$coefficients, rep(1, 12))
  expect_lt(calls, 5L)
  expect_equal(fit$dispersion, 0.5)
})

test_that("weighted offset fits use the final mean for estimated scale", {
  set.seed(204)
  X <- matrix(rnorm(160 * 4), 160)
  offset <- seq(-2, 2, length.out = 160)
  w <- seq(0.5, 2, length.out = 160)
  w[1:8] <- 0
  y <- offset + X[, 1] + rnorm(160)
  fit <- L0glm.fit(X, y, weights = w, offset = offset, family = gaussian(),
    dispersion = NULL, lambda = 1, relaxed = TRUE, maxit = 100, parallelize = FALSE)
  mu <- offset + as.numeric(X %*% fit$coefficients)
  expect_equal(fit$dispersion, sum(w * (y - mu)^2) /
                 (sum(w > 0) - sum(fit$coefficients != 0)))
})

test_that("free non-Gaussian families report fresh Pearson scale after refitting", {
  set.seed(205)
  X <- cbind(1, rnorm(160))
  families <- list(Gamma("log"), inverse.gaussian("log"),
                   quasipoisson("log"), quasibinomial("logit"))
  for (family in families) {
    mu <- family$linkinv(as.numeric(X %*% c(0.5, 0.2)))
    w <- rep(if (family$family == "quasibinomial") 20 else 1, 160)
    y <- switch(family$family,
      Gamma = rgamma(160, shape = 5, scale = mu / 5),
      inverse.gaussian = statmod::rinvgauss(160, mean = mu, dispersion = 0.1),
      quasipoisson = rpois(160, mu),
      quasibinomial = rbinom(160, 20, mu) / 20)
    fit <- L0glm.fit(X, y, weights = w, family = family, intercept = FALSE,
      always.include = 1L, dispersion = NULL, lambda = 0.1,
      relaxed = TRUE, maxit = 100, parallelize = FALSE)
    fitted <- family$linkinv(as.numeric(X %*% fit$coefficients))
    expected <- sum(w * (y - fitted)^2 / family$variance(fitted)) /
      (160 - sum(fit$coefficients != 0))
    expect_equal(fit$dispersion, expected, info = family$family)
    expect_true(is.finite(fit$dispersion) && fit$dispersion > 0)
  }
})

test_that("fixed and estimated channels survive grouped and relaxed passes", {
  set.seed(206)
  X <- matrix(rnorm(120 * 6), 120)
  y <- cbind(2 * X[, 1] + rnorm(120), -X[, 1] + rnorm(120, sd = 2))
  args <- list(X = X, y = y, family = gaussian(), dispersion = c(4, -0.5),
    lambda = "gic", grouped = TRUE, group.maxit = 2, relaxed = TRUE, maxit = 100)
  serial <- do.call(L0glm.fit, c(args, list(parallelize = FALSE)))
  expect_identical(serial$dispersion[1], 4)
  residual <- y[, 2] - as.numeric(X %*% serial$coefficients[, 2])
  expect_equal(serial$dispersion[2], sum(residual^2) /
                 (120 - sum(serial$coefficients[, 2] != 0)))
  parallel <- do.call(L0glm.fit, c(args, list(parallelize = TRUE, ncores = 2)))
  expect_equal(parallel$coefficients, serial$coefficients)
  expect_equal(parallel$dispersion, serial$dispersion)
})

test_that("known-scale families keep their dispersion on both execution paths", {
  X <- matrix(seq(0.2, 2, length.out = 30), ncol = 1)
  y <- rep(1, 30)
  for (openmp in c(FALSE, TRUE)) {
    expect_warning(fit <- L0glm.fit(X, y, family = poisson(identity),
      nonnegative = TRUE, dispersion = -2, openmp = openmp,
      parallelize = FALSE), "Dispersion is fixed at 1")
    expect_equal(fit$dispersion, 1)
  }
  for (family in list(binomial(), MASS::negative.binomial(theta = 2))) {
    expect_warning(fit <- L0glm.fit(X, y, family = family, dispersion = 3,
      intercept = FALSE, parallelize = FALSE), "Dispersion is fixed at 1")
    expect_equal(fit$dispersion, 1)
  }
})

test_that("invalid dispersion inputs fail before fitting", {
  X <- matrix(seq_len(20), ncol = 1)
  for (bad in list(0, NA_real_, Inf, -Inf, NaN, "2", numeric(), c(1, 2), list(1))) {
    expect_error(L0glm.fit(X, seq_len(20), dispersion = bad,
                          parallelize = FALSE), "dispersion must be")
    expect_error(L0nnpois(X, seq_len(20), dispersion = bad,
                         parallelize = FALSE), "dispersion must be")
  }
})
