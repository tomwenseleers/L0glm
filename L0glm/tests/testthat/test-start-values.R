test_that("coefficient mixing is no longer a fitting control", {
  for (fun in list(L0glm.fit, L0nnpois, L0glm.tune, L0glm.fit.internal))
    expect_false("warm.start" %in% names(formals(fun)))
  for (fun in list(L0glm.fit, L0nnpois))
    expect_error(do.call(fun, list(X = diag(3), y = 1:3, warm.start = 1)),
                 "unused argument")
})

test_that("explicit signed starts set the first ridge penalties without mixing", {
  start <- c(-1, 0, 2)
  y <- c(-2, 3, 4)
  for (sparse in c(FALSE, TRUE)) {
    fit <- L0glm.fit(Matrix::Matrix(diag(3), sparse = sparse), y,
                     family = gaussian(), nonnegative = FALSE,
                     start = start, lambda = 4, dispersion = 1,
                     maxit = 1, thresh = 0, relaxed = FALSE,
                     parallelize = FALSE, openmp = FALSE)
    expect_equal(as.numeric(fit$coefficients), c(-1, 0, 3.2))
    expect_equal(as.numeric(fit$lambdas)[c(1, 3)], c(1, 0.25))
  }
})

test_that("Gaussian initialization respects response units with fixed dispersion", {
  set.seed(16010)
  X <- matrix(rnorm(120 * 12), 120)
  beta <- c(-3, 2, rep(0, 10))
  y <- as.numeric(X %*% beta) + rnorm(120, sd = 0.5)
  fit <- function(scale) L0glm.fit(
    X, y * scale, family = gaussian(), nonnegative = FALSE,
    dispersion = 0.25 * scale^2, lambda = "gic", maxit = 100,
    relaxed = TRUE, parallelize = FALSE, openmp = FALSE)
  reference <- fit(1)$coefficients
  for (scale in c(0.01, 100))
    expect_equal(fit(scale)$coefficients / scale, reference, tolerance = 1e-8)
})

test_that("binomial initialization conditions slopes on the intercept", {
  x <- rep(seq(-1, 1, length.out = 20), 2)
  y <- round(1000 * plogis(-2 + 0.7 * x)) / 1000
  mu <- (1000 * y + 0.5) / 1001
  family <- binomial()
  eta <- family$linkfun(mu)
  W <- 1000 * family$mu.eta(eta)^2 / family$variance(mu)
  z <- eta + (y - mu) / family$mu.eta(eta)
  X <- cbind(1, x)
  lambda <- 0.3
  expected <- solve(crossprod(X, W * X) + diag(lambda, 2), crossprod(X, W * z))
  for (sparse in c(FALSE, TRUE)) {
    start <- init_coefs(Matrix::Matrix(X, sparse = sparse), y, mu,
                         weights = rep(1000, length(y)), family = family,
                         lambda = lambda)
    expect_gt(start[2], 0)
    expect_equal(start[2], as.numeric(expected[2]), tolerance = 1e-12)
  }
})
