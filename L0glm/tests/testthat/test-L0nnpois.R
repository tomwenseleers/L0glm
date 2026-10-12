test_that("L0nnpois shares fitting arguments but always enforces nonnegativity", {
  expect_identical(names(formals(L0nnpois)),
                   setdiff(names(formals(L0glm.fit)), c("nonnegative", "openmp")))
})

test_that("L0nnpois matches L0glm.fit for dense and sparse single-channel fits", {
  set.seed(1)
  n <- 40L
  p <- 20L
  X <- matrix(runif(n * p), n, p)
  beta <- rep(0, p)
  beta[c(3L, 8L)] <- c(2, 1)
  y <- rpois(n, pmax(as.numeric(X %*% beta), 0.1))

  fit_l0 <- L0glm.fit(
    X, y,
    openmp = FALSE,
    family = poisson(identity),
    lambda = "gic",
    nonnegative = TRUE,
    grouped = FALSE,
    relaxed = FALSE,
    parallelize = FALSE,
    solver = "eigen",
    method = 1L,
    warm.start = 1,
    maxit = 15L
  )
  fit_nn <- L0nnpois(
    X, y,
    lambda = "gic",
    grouped = FALSE,
    relaxed = FALSE,
    parallelize = FALSE,
    warm.start = 1,
    maxit = 15L
  )

  expect_equal(as.matrix(fit_nn$coefficients),
               as.matrix(fit_l0$coefficients),
               tolerance = 1e-12)
  expect_equal(fit_nn$lambdas, fit_l0$lambdas, tolerance = 1e-12)
  expect_equal(fit_nn$weights, fit_l0$weights, tolerance = 1e-10)
  expect_identical(fit_nn$z, fit_l0$z)
  expect_identical(fit_nn$dispersion, fit_l0$dispersion)

  set.seed(2)
  Xs <- Matrix::rsparsematrix(n, p, density = 0.15)
  Xs@x <- abs(Xs@x)
  beta[] <- 0
  beta[c(4L, 11L)] <- c(2.5, 1.5)
  ys <- rpois(n, pmax(as.numeric(Xs %*% beta), 0.1))

  fit_l0_s <- L0glm.fit(
    Xs, ys,
    openmp = FALSE,
    family = poisson(identity),
    lambda = "gic",
    nonnegative = TRUE,
    grouped = FALSE,
    relaxed = FALSE,
    parallelize = FALSE,
    solver = "eigen",
    method = 1L,
    warm.start = 1,
    maxit = 15L
  )
  fit_nn_s <- L0nnpois(
    Xs, ys,
    lambda = "gic",
    grouped = FALSE,
    relaxed = FALSE,
    parallelize = FALSE,
    warm.start = 1,
    maxit = 15L
  )

  expect_equal(as.matrix(fit_nn_s$coefficients),
               as.matrix(fit_l0_s$coefficients),
               tolerance = 1e-12)
  expect_equal(fit_nn_s$lambdas, fit_l0_s$lambdas, tolerance = 1e-8)
  expect_equal(fit_nn_s$weights, fit_l0_s$weights, tolerance = 1e-9)
  expect_identical(fit_nn_s$z, fit_l0_s$z)
  expect_identical(fit_nn_s$dispersion, fit_l0_s$dispersion)
})

test_that("L0nnpois ignores grouped outer passes for single-channel fits", {
  set.seed(4)
  n <- 50L
  p <- 25L
  X <- Matrix::rsparsematrix(n, p, density = 0.12)
  X@x <- abs(X@x)
  beta <- rep(0, p)
  beta[c(5L, 13L, 20L)] <- c(1.5, 2.0, 1.2)
  y <- rpois(n, pmax(as.numeric(X %*% beta), 0.1))

  grouped <- L0nnpois(
    X, y,
    lambda = "gic",
    grouped = TRUE,
    group.maxit = 4L,
    relaxed = FALSE,
    parallelize = FALSE,
    warm.start = 1,
    maxit = 10L
  )
  ungrouped <- L0nnpois(
    X, y,
    lambda = "gic",
    grouped = FALSE,
    group.maxit = 4L,
    relaxed = FALSE,
    parallelize = FALSE,
    warm.start = 1,
    maxit = 10L
  )

  expect_identical(grouped$coefficients, ungrouped$coefficients)
  expect_identical(grouped$lambdas, ungrouped$lambdas)
  expect_identical(grouped$weights, ungrouped$weights)
  expect_identical(grouped$z, ungrouped$z)
  expect_identical(grouped$dispersion, ungrouped$dispersion)
})

test_that("L0nnpois accepts NaturalOrdering control for sparse Cholesky", {
  set.seed(5)
  n <- 70L
  p <- 40L
  m <- 4L
  X <- Matrix::rsparsematrix(n, p, density = 0.1)
  X@x <- abs(X@x)
  beta <- matrix(0, p, m)
  beta[c(6L, 15L, 31L), ] <- c(2.1, 1.4, 1.8)
  y <- matrix(rpois(n * m, pmax(as.vector(X %*% beta), 0.1)), n, m)

  default <- L0nnpois(
    X, y,
    lambda = "gic",
    grouped = TRUE,
    relaxed = FALSE,
    parallelize = FALSE,
    warm.start = 1,
    maxit = 10L,
    group.maxit = 2L
  )
  natural <- L0nnpois(
    X, y,
    lambda = "gic",
    grouped = TRUE,
    relaxed = FALSE,
    parallelize = FALSE,
    warm.start = 1,
    maxit = 10L,
    group.maxit = 2L,
    control.solver = list(NaturalOrdering = TRUE)
  )
  amd <- L0nnpois(
    X, y,
    lambda = "gic",
    grouped = TRUE,
    relaxed = FALSE,
    parallelize = FALSE,
    warm.start = 1,
    maxit = 10L,
    group.maxit = 2L,
    control.solver = list(NaturalOrdering = FALSE)
  )

  expect_equal(as.matrix(natural$coefficients),
               as.matrix(amd$coefficients), tolerance = 1e-10)
  expect_equal(as.matrix(default$coefficients),
               as.matrix(amd$coefficients), tolerance = 1e-10)
  expect_equal(as.matrix(natural$lambdas),
               as.matrix(amd$lambdas), tolerance = 1e-8)
  expect_equal(as.matrix(natural$weights),
               as.matrix(amd$weights), tolerance = 1e-8)
  expect_equal(as.matrix(natural$z),
               as.matrix(amd$z), tolerance = 1e-12)
  expect_identical(natural$dispersion, amd$dispersion)
})

test_that("L0nnpois OpenMP multichannel fits match serial fits", {
  ncores <- parallel::detectCores()
  skip_if(is.na(ncores) || ncores < 2L, "parallel test requires at least two cores")

  set.seed(3)
  n <- 60L
  p <- 35L
  m <- 5L
  X <- Matrix::rsparsematrix(n, p, density = 0.1)
  X@x <- abs(X@x)
  beta <- matrix(0, p, m)
  beta[c(5L, 12L, 25L), ] <- c(2.5, 1.2, 2)
  y <- matrix(rpois(n * m, pmax(as.vector(X %*% beta), 0.1)), n, m)

  serial <- L0nnpois(
    X, y,
    lambda = "gic",
    grouped = TRUE,
    relaxed = TRUE,
    parallelize = FALSE,
    warm.start = 1,
    maxit = 12L,
    group.maxit = 3L
  )
  parallel <- L0nnpois(
    X, y,
    lambda = "gic",
    grouped = TRUE,
    relaxed = TRUE,
    parallelize = TRUE,
    ncores = 2L,
    warm.start = 1,
    maxit = 12L,
    group.maxit = 3L
  )

  expect_identical(parallel$coefficients, serial$coefficients)
  expect_identical(parallel$lambdas, serial$lambdas)
  expect_identical(parallel$weights, serial$weights)
  expect_identical(parallel$z, serial$z)
  expect_identical(parallel$dispersion, serial$dispersion)
})
