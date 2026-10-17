test_that("a virtually interpolating pilot does not initialize at roundoff variance", {
  set.seed(1641)
  X0 <- matrix(rnorm(40 * 80), 40, 80)
  y0 <- 2 * X0[, 1] + rnorm(40)
  X <- rbind(X0, X0); y <- rep(y0, 2)
  preprocessing <- L0glm_preprocessing(X, "auto", rep(-Inf, 80), integer())
  initialize <- function(y) L0glm_initial_dispersion(X, y, y, rep(1, 80),
    rep(0, 80), gaussian(), preprocessing = preprocessing)
  rng <- .Random.seed
  phi <- initialize(y)
  expect_identical(.Random.seed, rng)
  expect_gt(phi, 1e-6 * mean(y^2) / preset.lambda("gic", 80, 80))
  expect_equal(initialize(10 * y), 100 * phi, tolerance = 1e-7)
  expect_true(is.finite(phi))
})
