test_that("all-zero starts fit only mandatory terms and retain sparse dimensions", {
  X <- matrix(seq_len(30), 10, 3, dimnames = list(NULL, letters[1:3]))
  for (sparse in c(FALSE, TRUE)) for (m in c(1L, 2L)) {
    Y <- matrix(rep(seq_len(m) + 2, each = 10), 10, m)
    for (intercept in c(FALSE, TRUE)) for (backend in c("general", "native")) {
      args <- list(X = Matrix::Matrix(X, sparse = sparse), y = Y, start = 0,
                   intercept = intercept, lambda = "gic", maxit = 15L,
                   grouped = TRUE, relaxed = TRUE, parallelize = FALSE)
      fit <- if (backend == "native") do.call(L0nnpois, args) else
        do.call(L0glm.fit, c(args, list(family = poisson(identity),
                                      nonnegative = TRUE, openmp = FALSE)))
      expect_s4_class(fit$coefficients, "dgCMatrix")
      expect_equal(dim(fit$coefficients), c(3L + as.integer(intercept), m))
      if (intercept) {
        expect_equal(as.numeric(fit$coefficients[1, ]), seq_len(m) + 2, tolerance = 1e-8)
        expect_true(all(fit$coefficients[-1, ] == 0))
      } else expect_true(all(fit$coefficients == 0))
    }
  }
})

test_that("entry masks persist through grouped passes and relaxed refits", {
  set.seed(1650)
  X <- matrix(runif(30 * 6), 30, 6, dimnames = list(NULL, letters[1:6]))
  Y <- matrix(rpois(90, 4), 30, 3, dimnames = list(NULL, LETTERS[1:3]))
  start <- rbind(c(1, 0, 1), c(0, 1, 0), c(0, 0, 0),
                 c(1, 1, 0), c(0, 1, 1), c(0, 0, 0))
  for (sparse in c(FALSE, TRUE)) for (grouped in c(FALSE, TRUE)) {
    args <- list(X = Matrix::Matrix(X, sparse = sparse), y = Y, start = start,
                 lambda = 0.05, maxit = 15L, group.maxit = 3L,
                 grouped = grouped, relaxed = TRUE, parallelize = FALSE)
    native <- do.call(L0nnpois, args)
    general <- do.call(L0glm.fit, c(args, list(family = poisson(identity),
      nonnegative = TRUE, openmp = FALSE)))
    for (fit in list(native, general)) {
      expect_s4_class(fit$coefficients, "dgCMatrix")
      expect_equal(dimnames(fit$coefficients), list(letters[1:6], LETTERS[1:3]))
      expect_true(all(as.matrix(fit$coefficients)[start == 0] == 0))
      expect_true(all(as.matrix(fit$lambdas)[start == 0] == 1e11))
    }
    expect_equal(as.matrix(native$coefficients), as.matrix(general$coefficients), tolerance = 1e-7)
    parallel <- do.call(L0glm.fit, c(modifyList(args, list(parallelize = TRUE, ncores = 2L)),
                                   list(family = poisson(identity), nonnegative = TRUE, openmp = FALSE)))
    expect_equal(parallel$coefficients, general$coefficients, tolerance = 1e-8)
    threaded <- do.call(L0nnpois, modifyList(args, list(parallelize = TRUE, ncores = 2L)))
    expect_equal(threaded$coefficients, native$coefficients, tolerance = 1e-8)
  }
})

test_that("fixed zeros are removed before Gaussian dispersion initialization", {
  set.seed(1651)
  X <- matrix(rnorm(40 * 7), 40, 7)
  seen <- integer()
  original <- L0glm_initial_dispersion
  local_mocked_bindings(L0glm_initial_dispersion = function(X, ...) {
    seen <<- c(seen, ncol(X))
    original(X, ...)
  })
  fit <- L0glm.fit(X, rnorm(40), start = c(1, 0, 1, 0, 0, 0, 0),
                   family = gaussian(), parallelize = FALSE)
  expect_equal(seen, 2L)
  expect_equal(dim(fit$coefficients), c(7L, 1L))
  expect_true(all(fit$coefficients[-c(1, 3), ] == 0))
  fit <- L0glm.fit(X, rnorm(40), start = 0, family = gaussian(),
                   dispersion = -2, parallelize = FALSE)
  expect_true(is.finite(fit$dispersion) && fit$dispersion > 0)
  expect_true(all(fit$coefficients == 0))
})

test_that("zero-start mandatory coefficients remain eligible and unpenalized", {
  X <- cbind(0, seq_len(10), seq_len(10)^2)
  for (backend in c("general", "native")) {
    args <- list(X = X, y = 2 * X[, 2], start = 0, always.include = 2L,
                 thresh = 2, lambda = 1e5, relaxed = TRUE,
                 parallelize = FALSE, maxit = 10L)
    fit <- if (backend == "native") do.call(L0nnpois, args) else
      do.call(L0glm.fit, c(args, list(family = poisson(identity), nonnegative = TRUE, openmp = FALSE)))
    expect_equal(as.numeric(fit$coefficients), c(0, 2, 0), tolerance = 1e-8)
    expect_equal(as.numeric(fit$lambdas)[2], 0)
  }
  fit <- L0glm.fit(matrix(seq_len(10), 10, 1), rep(0.2, 10), weights = rep(20, 10),
                   family = binomial(), start = 0, intercept = TRUE,
                   nonnegative = TRUE, parallelize = FALSE)
  expect_equal(as.numeric(fit$coefficients), c(qlogis(0.2), 0), tolerance = 1e-7)
})
