test_that("eligible L0glm.fit calls forward all fitting controls to L0nnpois", {
  local_mocked_bindings(L0nnpois = function(...) list(args = list(...)))
  X <- matrix(1, 5, 3)
  y <- cbind(1:5, 2:6)
  args <- list(
    X = X, y = y, family = poisson(identity), nonnegative = TRUE,
    weights = rep(2, 5), dispersion = 1, start = matrix(1, 3, 2),
    offset = rep(0.1, 5), lambda = "bic", intercept = TRUE,
    always.include = 2L, normalize = TRUE, solver = "eigen", method = 1L,
    parallelize = TRUE, ncores = 2L, nthreads = 1L,
    maxit = 12L, grouped = FALSE, group.maxit = 3L, relaxed = TRUE,
    prior.mean = matrix(0.2, 3, 2), gamma = 1.8,
    thresh = 1e-6, tol = 1e-8,
    control.solver = list(NaturalOrdering = TRUE), verbose = TRUE
  )
  fit <- do.call(L0glm.fit, args)
  expect_equal(fit$args, args[names(fit$args)])
})

test_that("OpenMP dispatch handles vector responses and family functions", {
  local_mocked_bindings(L0nnpois = function(...) list(args = list(...)))
  fit <- L0glm.fit(diag(4), 1:4, nonnegative = TRUE,
                   family = function() poisson(identity))
  expect_false(fit$args$relaxed)
  expect_false(fit$args$intercept)
  expect_identical(fit$args$family$family, "poisson")
  # L0glm.fit retains its own public lambda default; dispatch does not change it.
  expect_identical(fit$args$lambda, 1)
})

test_that("ineligible or disabled OpenMP fits retain the general R path", {
  local_mocked_bindings(L0nnpois = function(...) stop("Unexpected dispatch"))
  base <- list(X = diag(4), y = c(4, 3, 2, 1), family = poisson(identity),
               nonnegative = TRUE, parallelize = FALSE, grouped = FALSE,
               relaxed = FALSE, lambda = 0.01)
  variants <- list(
    list(openmp = FALSE), list(nonnegative = FALSE),
    list(family = gaussian(identity)), list(family = poisson(log)),
    list(method = 2L), list(method = 3L), list(method = 4L),
    list(method = 5L), list(method = 6L), list(solver = "osqp"),
    list(control.solver = control.solver.gen("eigen", 2L)),
    list(control.solver = control.solver.gen("eigen", 5L)),
    list(control.solver = control.solver.gen("eigen", 1L, rowaug = FALSE))
  )
  for (variant in variants) {
    fit <- do.call(L0glm.fit, utils::modifyList(base, variant))
    expect_equal(dim(fit$coefficients), c(4L, 1L))
    expect_true(all(is.finite(fit$coefficients)))
  }
  for (invalid in list(NA, "true", 1, logical(), c(TRUE, FALSE))) {
    expect_error(L0glm.fit(diag(2), 1:2, openmp = invalid),
                 "openmp must be TRUE or FALSE", fixed = TRUE)
  }
})

test_that("dispatched grouped fits agree with direct and R fits", {
  set.seed(20261012)
  X <- matrix(runif(60 * 12), 60, 12)
  beta <- matrix(0, 12, 3)
  beta[c(2, 6, 10), ] <- c(8, 12, 6)
  y <- matrix(rpois(180, as.vector(X %*% beta)), 60, 3)
  for (sparse in c(FALSE, TRUE)) {
    for (grouped in c(FALSE, TRUE)) {
      args <- list(X = Matrix::Matrix(X, sparse = sparse), y = y,
                   weights = seq(0.7, 1.2, length.out = 60),
                   offset = Matrix::Matrix(0.2, 60, 3),
                   lambda = "gic", grouped = grouped, relaxed = TRUE,
                   parallelize = FALSE, group.maxit = 2L, maxit = 15L)
      direct <- do.call(L0nnpois, args)
      dispatched <- do.call(L0glm.fit, c(args, list(
        family = poisson(identity), nonnegative = TRUE)))
      general <- do.call(L0glm.fit, c(args, list(
        family = poisson(identity), nonnegative = TRUE, openmp = FALSE)))
      expect_equal(dispatched, direct)
      expect_equal(as.matrix(dispatched$coefficients), as.matrix(general$coefficients),
                   tolerance = 1e-8)
    }
  }
})

test_that("L0nnpois defaults to the GIC penalty", {
  expect_identical(formals(L0nnpois)$lambda, "gic")
  implicit <- L0nnpois(diag(6), c(0, 3, 0, 5, 0, 0), parallelize = FALSE)
  explicit <- L0nnpois(diag(6), c(0, 3, 0, 5, 0, 0), lambda = "gic", parallelize = FALSE)
  expect_equal(implicit, explicit)
  expect_equal(implicit$init.lambda, log(6) * log(log(6)))
})

test_that("Poisson dispatch rejects invalid responses", {
  for (bad in list(c(-1, 2), c(NA, 2), c(Inf, 2))) {
    expect_error(L0glm.fit(diag(2), bad, family = poisson(identity), nonnegative = TRUE),
                 "finite and nonnegative")
    expect_error(L0nnpois(diag(2), bad), "finite and nonnegative")
  }
})
