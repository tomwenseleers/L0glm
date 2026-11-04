test_that("one native driver supports every weighted-ridge solver", {
  set.seed(781)
  X <- matrix(rnorm(120 * 12, sd = .2), 120, 12)
  beta <- matrix(0, 12, 3); beta[1:3, ] <- c(.8, -.5, .4)
  y <- X %*% beta + matrix(rnorm(360, sd = .2), 120, 3)
  original <- glm_native_cpp
  calls <- 0L
  local_mocked_bindings(glm_native_cpp = function(...) {
    calls <<- calls + 1L
    original(...)
  })
  for (sparse in c(FALSE, TRUE)) {
    for (solver in c(paste0("eigen", 1:6), "osqp", "bigGlm")) {
      eigen <- startsWith(solver, "eigen")
      args <- list(X = Matrix::Matrix(X, sparse = sparse), y = y,
        family = gaussian(), lambda = "gic", intercept = TRUE,
        grouped = TRUE, group.maxit = 2L, relaxed = TRUE, maxit = 12L,
        parallelize = FALSE, solver = if (eigen) "eigen" else solver,
        method = if (eigen) as.integer(sub("eigen", "", solver)) else 1L)
      reference <- suppressWarnings(do.call(L0glm.fit, c(args, list(openmp = FALSE))))
      before <- calls
      native <- suppressWarnings(do.call(L0glm.fit, args))
      expect_equal(calls, before + 1L, info = solver)
      expect_equal(native$coefficients, reference$coefficients, tolerance = 1e-7)
      expect_identical(as.numeric(native$coefficients) != 0,
        as.numeric(reference$coefficients) != 0)
      expect_equal(native$dispersion, reference$dispersion, tolerance = 1e-7)
      expect_equal(GIC(native), GIC(reference), tolerance = 1e-8)
      expect_identical(native$convergence$iterations, reference$convergence$iterations)
    }
  }
})

test_that("native grouped fits retain masks, signed starts, priors and input metadata", {
  set.seed(784)
  X <- matrix(rnorm(140 * 14), 140, 14)
  colnames(X) <- paste0("x", 1:14)
  beta <- matrix(0, 14, 4)
  beta[1:3, ] <- outer(c(.8, -.5, .4), c(1, -1, 1, -1))
  y <- X %*% beta + matrix(rnorm(560, sd = .1), 140, 4)
  colnames(y) <- paste0("channel", 1:4)
  start <- matrix(1, 15, 4)
  start[1:2, ] <- 0 # Intercept and compulsory slope remain eligible.
  start[7:9, ] <- 0; start[5:6, 1] <- 0
  weights <- seq(.7, 1.3, length.out = 140)
  offset <- seq(0, .2, length.out = 140)
  for (normalize in c("none", "rms", "l2", "sd")) {
    args <- list(X = X, y = y, lambda = "gic", intercept = TRUE,
      always.include = 1L, start = start, prior.mean = rep(.01, 15),
      weights = weights, offset = offset, normalize = normalize,
      grouped = TRUE, relaxed = TRUE, group.maxit = 2L, parallelize = FALSE)
    reference <- suppressWarnings(do.call(L0glm.fit, c(args, list(openmp = FALSE))))
    native <- suppressWarnings(do.call(L0glm.fit, args))
    expect_equal(native$coefficients, reference$coefficients, tolerance = 1e-7)
    expect_equal(native$dispersion, reference$dispersion, tolerance = 1e-7)
    expect_true(all(native$coefficients[7:9, ] == 0))
    expect_true(all(native$coefficients[5:6, 1] == 0))
    expect_identical(dimnames(native$coefficients), dimnames(reference$coefficients))
    expect_identical(native$prior.weights, reference$prior.weights)
    expect_identical(native$offset, reference$offset)
    expect_identical(native$preprocessing, reference$preprocessing)
  }
})

test_that("native dispersion pilots preserve RNG and fixed/free scale conventions", {
  set.seed(785)
  X <- matrix(rnorm(180 * 12), 180, 12)
  y <- 1.2 * X[, 1] - .6 * X[, 2] + rnorm(180, sd = .4)
  for (dispersion in list(NULL, .16, -.1)) {
    args <- list(X = X, y = y, dispersion = dispersion, lambda = "gic",
      relaxed = TRUE, parallelize = FALSE)
    state <- .Random.seed
    native <- do.call(L0glm.fit, args)
    expect_identical(.Random.seed, state)
    reference <- do.call(L0glm.fit, c(args, list(openmp = FALSE)))
    expect_equal(native$coefficients, reference$coefficients, tolerance = 1e-7)
    expect_equal(native$dispersion, reference$dispersion, tolerance = 1e-7)
    if (identical(dispersion, .16)) expect_equal(native$dispersion, .16)
  }
})

test_that("Poisson identity retains its specialization and custom callbacks remain supported", {
  local_mocked_bindings(glm_native_cpp = function(...) stop("unexpected general native call"))
  X <- diag(8); y <- c(0, 5, 0, 9, 0, 3, 0, 0)
  fit <- L0glm.fit(X, y, family = poisson("identity"), nonnegative = TRUE,
    lambda = "gic", parallelize = FALSE)
  direct <- L0nnpois(X, y, lambda = "gic", parallelize = FALSE)
  expect_equal(fit[setdiff(names(fit), "family")], direct[setdiff(names(direct), "family")])
  expect_identical(fit$family$family, direct$family$family)
  expect_identical(fit$family$link, direct$family$link)
  family <- gaussian(); family$variance <- function(mu) rep(2, length(mu))
  expect_null(L0glm_native_family(family))
  custom <- L0glm.fit(X, y, family = family, parallelize = FALSE)
  expect_true(all(is.finite(custom$coefficients)))
  poisson_custom <- poisson("identity")
  poisson_custom$variance <- function(mu) mu + mu^2
  expect_null(L0glm_native_family(poisson_custom))
  expect_s3_class(L0glm.fit(X, y, family = poisson_custom,
    nonnegative = TRUE, parallelize = FALSE), "L0glm_fit")
  nb_custom <- MASS::negative.binomial(3)
  nb_custom$variance <- function(mu) mu + mu^2
  expect_null(L0glm_native_family(nb_custom))
})

test_that("native rectangular solver fallbacks remain safe in OpenMP workers", {
  set.seed(787)
  X <- matrix(rnorm(90 * 8), 90, 8)
  y <- cbind(X[, 1] + rnorm(90, sd = .2), X[, 2] + rnorm(90, sd = .2))
  control <- control.solver.gen("eigen", 6L)
  control$rowaug <- TRUE
  control$svd <- control$sparse.svd <- 1L
  control$precond <- 3L # NoQR cannot factor a rectangular augmented design.
  for (sparse in c(FALSE, TRUE)) {
    design <- Matrix::Matrix(X, sparse = sparse)
    reference <- L0glm.fit(design, y, lambda = "gic", grouped = FALSE,
      relaxed = FALSE, parallelize = FALSE)
    expect_warning(fit <- L0glm.fit(design, y, lambda = "gic", grouped = FALSE,
      relaxed = FALSE, method = 6L, control.solver = control,
      parallelize = TRUE, ncores = 2L), "Cholesky LLT was used as fallback")
    expect_equal(fit$coefficients, reference$coefficients, tolerance = 1e-7)
  }
})

test_that("native OpenMP restores caller thread settings", {
  set.seed(786)
  X <- matrix(rnorm(100 * 8), 100, 8)
  y <- cbind(X[, 1] + rnorm(100), -X[, 1] + rnorm(100))
  blas <- RhpcBLASctl::blas_get_num_procs()
  omp <- RhpcBLASctl::omp_get_max_threads()
  fit <- suppressWarnings(L0glm.fit(X, y, parallelize = TRUE, ncores = 2L, lambda = "gic"))
  expect_s3_class(fit, "L0glm_fit")
  expect_equal(RhpcBLASctl::blas_get_num_procs(), blas)
  expect_equal(RhpcBLASctl::omp_get_max_threads(), omp)
})

test_that("native Cholesky honors an ordering flag without an integer ordering", {
  set.seed(788)
  X <- Matrix::rsparsematrix(100, 12, density = .3)
  y <- as.numeric(X[, 1]) + rnorm(100, sd = .2)
  control <- control.solver.gen("eigen", 1L, NaturalOrdering = TRUE)
  control$ordering <- NULL
  seen <- NULL
  original <- glm_native_cpp
  local_mocked_bindings(glm_native_cpp = function(..., solver_control) {
    seen <<- solver_control$ordering
    original(..., solver_control = solver_control)
  })
  native <- L0glm.fit(X, y, lambda = "gic", control.solver = control, parallelize = FALSE)
  reference <- L0glm.fit(X, y, lambda = "gic", control.solver = control,
    openmp = FALSE, parallelize = FALSE)
  expect_identical(seen, 1L)
  expect_equal(native$coefficients, reference$coefficients, tolerance = 1e-7)
})
