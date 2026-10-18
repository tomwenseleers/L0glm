test_that("dense predictions retain offsets, weights, and excluded starting columns", {
  set.seed(1871)
  X <- matrix(rnorm(400*12), 400, 12)
  Z <- cbind(1, X)
  eligible <- setdiff(seq_len(ncol(Z)), c(4L, 8L))
  beta <- numeric(ncol(Z)); beta[eligible] <- seq(-.15, .15, length.out = length(eligible))
  offset <- seq(-.1, .1, length.out = nrow(X))
  weights <- sample(1:3, nrow(X), replace = TRUE)
  eta <- as.numeric(Z %*% beta) + offset
  for (family in list(gaussian(), poisson(), binomial())) {
    mu <- family$linkinv(eta)
    y <- switch(family$family, gaussian = rnorm(nrow(X), mu, .2),
      poisson = rpois(nrow(X), mu), binomial = rbinom(nrow(X), 1, mu))
    start <- matrix(as.numeric(seq_len(ncol(Z)) %in% eligible), ncol = 1L)
    actual <- L0glm.fit(X, y, family = family, weights = weights, offset = offset,
      intercept = TRUE, nonnegative = FALSE, start = start,
      lambda = 0, thresh = 0, maxit = 100L, tol = 1e-10)
    expected <- stats::glm.fit(Z[, eligible], y, family = family,
      weights = weights, offset = offset, control = list(maxit = 100L, epsilon = 1e-10))
    b <- as.numeric(actual$coefficients)
    expect_identical(b[-eligible], c(0, 0))
    expect_equal(b[eligible], unname(expected$coefficients), tolerance = 1e-6)
    expect_equal(as.numeric(Z %*% b), as.numeric(Z[, eligible] %*% expected$coefficients), tolerance = 1e-6)
  }
})
