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
    method = 1L,
    maxit = 10L,
    group.maxit = 2L,
    warm.start = 0.5
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
    method = 1L,
    maxit = 10L,
    group.maxit = 2L,
    warm.start = 0.5
  )

  expect_identical(parallel$coefficients, serial$coefficients)
  expect_identical(parallel$lambdas, serial$lambdas)
  expect_identical(parallel$weights, serial$weights)
  expect_identical(parallel$z, serial$z)
  expect_identical(parallel$dispersion, serial$dispersion)
})

test_that("grouped relaxed fit keeps rows with cancelling channel coefficients", {
  X <- diag(3)
  y <- cbind(c(1, 0, 0), c(-1, 0, 0))

  fit <- L0glm.fit(
    X, y,
    family = gaussian(identity),
    lambda = 0.01,
    nonnegative = FALSE,
    grouped = TRUE,
    relaxed = TRUE,
    parallelize = FALSE,
    intercept = FALSE,
    normalize = FALSE,
    maxit = 10L,
    group.maxit = 2L,
    warm.start = 1
  )

  coefs <- as.matrix(fit$coefficients)
  lambdas <- as.matrix(fit$lambdas)
  expect_equal(coefs[1L, ], c(1, -1), tolerance = 1e-8)
  expect_true(all(coefs[-1L, ] == 0))
  expect_lt(max(abs(lambdas[1L, ])), 1e-25)
})

test_that("serial nthreads setting is restored after fit", {
  skip_if_not_installed("RhpcBLASctl")

  before_blas <- RhpcBLASctl::blas_get_num_procs()
  before_omp <- RhpcBLASctl::omp_get_max_threads()
  on.exit({
    utils::capture.output(RhpcBLASctl::blas_set_num_threads(before_blas))
    utils::capture.output(RhpcBLASctl::omp_set_num_threads(before_omp))
  }, add = TRUE)

  sim <- simul(method = 4, n = 40L, p = 40L, s = 0.1, seed = 2)
  fit <- L0glm.fit(
    sim$X, sim$y[, 1L, drop = FALSE],
    family = poisson(identity),
    lambda = "gic",
    nonnegative = TRUE,
    grouped = FALSE,
    relaxed = FALSE,
    parallelize = FALSE,
    nthreads = 1L,
    solver = "eigen",
    method = 1L,
    warm.start = 1
  )

  expect_s4_class(fit$coefficients, "dgCMatrix")
  expect_identical(RhpcBLASctl::blas_get_num_procs(), before_blas)
  expect_identical(RhpcBLASctl::omp_get_max_threads(), before_omp)
})
