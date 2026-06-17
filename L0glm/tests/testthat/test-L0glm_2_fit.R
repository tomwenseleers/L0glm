test_that("multichannel parallel fit matches serial fit", {
  ncores <- parallel::detectCores()
  skip_if(is.na(ncores) || ncores < 2L, "parallel test requires at least two cores")

  set.seed(42)
  n <- 18L
  p <- 10L
  m <- 4L

  X <- matrix(runif(n * p), n, p)
  beta <- matrix(0, p, m)
  beta[c(2L, 6L), ] <- c(1.5, 0.8)
  mu <- pmax(X %*% beta, 0.1)
  y <- matrix(rpois(n * m, as.vector(mu)), n, m)

  serial <- L0glm.fit(
    X, y,
    family = poisson(identity),
    lambda = "gic",
    nonnegative = TRUE,
    grouped = TRUE,
    relaxed = TRUE,
    parallelize = FALSE,
    solver = "eigen",
    method = 7L,
    maxit = 10L,
    group.maxit = 2L
  )

  parallel <- L0glm.fit(
    X, y,
    family = poisson(identity),
    lambda = "gic",
    nonnegative = TRUE,
    grouped = TRUE,
    relaxed = TRUE,
    parallelize = TRUE,
    ncores = 2L,
    solver = "eigen",
    method = 7L,
    maxit = 10L,
    group.maxit = 2L
  )

  expect_identical(parallel$coefficients, serial$coefficients)
  expect_identical(parallel$lambdas, serial$lambdas)
  expect_identical(parallel$weights, serial$weights)
  expect_identical(parallel$z, serial$z)
  expect_identical(parallel$dispersion, serial$dispersion)
})
