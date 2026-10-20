test_that("near-square Gaussian pilots retain enough noise to select sparse models", {
  set.seed(70201)
  X <- matrix(rnorm(50 * 50), 50, 50)
  beta <- numeric(50)
  beta[round(seq(1, 50, length.out = 5))] <- c(1, -1, 1, -1, 1)
  y <- as.numeric(X %*% beta) + rnorm(50, sd = sqrt(5 / 6))
  rng <- .Random.seed
  fit <- L0glm.fit(X, y, lambda = "gic", relaxed = TRUE, maxit = 100,
                   screen = "none", parallelize = FALSE)
  expect_identical(.Random.seed, rng)
  expect_lt(sum((as.numeric(coef(fit)) - beta)^2) / 5, 0.05)
  expect_equal(sum(coef(fit)[beta != 0, 1] != 0), 5)
  expect_lte(sum(coef(fit)[beta == 0, 1] != 0), 2)
  expect_gt(fit$dispersion, 0.4)
})

test_that("guarded normalized pilots preserve units, storage and nuisance profiling", {
  set.seed(2812)
  X <- cbind(Intercept = 1, matrix(rnorm(45 * 59), 45, 59))
  w <- seq(0.5, 2, length.out = 45)
  off <- seq(-1, 1, length.out = 45)
  y <- 2 + X[, 2] - X[, 4] + off + rnorm(45)
  initial <- function(X, y, offset = off) {
    prep <- L0glm_preprocessing(X, "rms", rep(-Inf, 60), 1L)
    L0glm_initial_dispersion(X, y, y, w, offset, gaussian(), 1L, prep)
  }
  rng <- .Random.seed
  phi <- initial(X, y)
  expect_identical(.Random.seed, rng)
  expect_true(is.finite(phi) && phi > 0)
  expect_equal(initial(X, 10 * y, 10 * off), 100 * phi, tolerance = 1e-7)
  expect_equal(initial(Matrix::Matrix(X, sparse = TRUE), y), phi, tolerance = 1e-7)
  shifted <- sweep(X, 2, c(0, seq(-2, 2, length.out = 59)), "+")
  expect_equal(initial(shifted, y), phi, tolerance = 1e-7)
})

test_that("a supplied dispersion bypasses the automatic pilot", {
  testthat::local_mocked_bindings(L0glm_canonical_dispersion = function(...) {
    stop("automatic pilot was called")
  })
  set.seed(2813)
  X <- matrix(rnorm(40 * 40), 40)
  y <- X[, 1] + rnorm(40)
  for (value in c(1, -1)) {
    fit <- L0glm.fit(X, y, dispersion = value, lambda = "gic",
                     maxit = 1, parallelize = FALSE)
    expect_identical(fit$dispersion, 1)
    expect_true(all(is.finite(coef(fit))))
  }
})
