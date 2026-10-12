test_that("the public constraint arguments are restricted to the supported API", {
  X <- matrix(seq_len(8), ncol = 1)
  y <- seq_len(8)
  expect_identical(formals(L0glm.fit)$nonnegative, FALSE)
  for (arg in c("lower", "upper")) {
    expect_false(arg %in% names(formals(L0glm.fit)))
    expect_false(arg %in% names(formals(L0nnpois)))
    expect_error(do.call(L0glm.fit, c(list(X = X, y = y), setNames(list(0), arg))),
                 "unused argument")
    expect_error(do.call(L0nnpois, c(list(X = X, y = y), setNames(list(0), arg))),
                 "unused argument")
  }
  expect_false("nonnegative" %in% names(formals(L0nnpois)))
  expect_error(L0nnpois(X, y, nonnegative = TRUE), "unused argument")
  expect_error(L0nnpois(X, y, nonnegative = FALSE), "unused argument")
  for (invalid in list(NA, NULL, logical(), c(TRUE, FALSE), 1, "TRUE")) {
    expect_error(L0glm.fit(X, y, nonnegative = invalid),
                 "nonnegative must be TRUE or FALSE", fixed = TRUE)
  }
})

test_that("all ridge backends respect optional nonnegativity including Gaussian intercepts", {
  X <- cbind(x = rep(c(-1, 1), 12), z = rep(c(-1, -1, 1, 1), 6))
  y <- -2 + 3 * X[, 1] - 1.5 * X[, 2]
  backends <- c(rep("eigen", 6), "osqp", "bigGlm", "nnls")
  methods <- c(1:6, 1L, 1L, 1L)
  for (sparse in c(FALSE, TRUE)) {
    for (normalize in c(FALSE, TRUE)) {
      for (i in seq_along(backends)) {
        for (nonnegative in c(FALSE, TRUE)) {
          fit <- L0glm.fit(
            Matrix::Matrix(X, sparse = sparse), y, family = gaussian(),
            intercept = TRUE, nonnegative = nonnegative, normalize = normalize,
            always.include = 2L,
            lambda = 1e-12, dispersion = 1, solver = backends[i], method = methods[i],
            warm.start = 1, relaxed = TRUE, parallelize = FALSE, maxit = 30
          )
          expected <- if (nonnegative) c(0, 3, 0) else c(-2, 3, -1.5)
          expect_equal(as.numeric(fit$coefficients), expected, tolerance = 1e-5,
                       info = paste(backends[i], methods[i], sparse, normalize, nonnegative))
        }
      }
    }
  }
})

test_that("binomial intercepts remain free across solvers and explicit intercept columns", {
  x <- rep(seq(-1, 1, length.out = 20), 2)
  y <- round(1000 * plogis(-2 + 0.7 * x)) / 1000
  expected <- unname(coef(glm(y ~ x, family = binomial(), weights = rep(1000, length(y)))))
  for (solver in c("eigen", "osqp", "bigGlm", "nnls")) {
    for (sparse in c(FALSE, TRUE)) {
      for (explicit in c(FALSE, TRUE)) {
        X <- if (explicit) cbind(x, 1) else matrix(x, ncol = 1)
        fit <- L0glm.fit(
          Matrix::Matrix(X, sparse = sparse), y, family = binomial,
          weights = rep(1000, length(y)), intercept = !explicit,
          nonnegative = TRUE, normalize = TRUE, lambda = 1e-12,
          solver = solver, warm.start = 0.9, relaxed = TRUE,
          parallelize = FALSE, maxit = 60
        )
        beta <- as.numeric(fit$coefficients)
        if (explicit) beta <- rev(beta)
        expect_true(beta[1] < 0, info = paste(solver, sparse, explicit))
        expect_gte(beta[2], 0)
        expect_equal(beta, expected, tolerance = 1e-4)
      }
    }
  }
})

test_that("binomial slopes remain constrained through grouped, relaxed and parallel fits", {
  x <- seq(-1, 1, length.out = 30)
  y <- cbind(round(1000 * plogis(-2 + x)), round(1000 * plogis(-1 - x))) / 1000
  for (grouped in c(FALSE, TRUE)) {
    for (relaxed in c(FALSE, TRUE)) {
      args <- list(
        X = Matrix::Matrix(matrix(x, ncol = 1), sparse = TRUE), y = y,
        family = binomial(), weights = rep(1000, length(x)),
        nonnegative = TRUE, lambda = 1e-8, warm.start = 0.9,
        grouped = grouped, relaxed = relaxed, group.maxit = 2L, maxit = 50L
      )
      serial <- do.call(L0glm.fit, c(args, list(parallelize = FALSE)))
      beta <- as.matrix(serial$coefficients)
      expect_true(all(beta[1, ] < 0))
      expect_true(all(beta[2, ] >= 0))
      expect_gt(beta[2, 1], 0.9)
      expect_equal(beta[2, 2], 0)
      parallel <- do.call(L0glm.fit, c(args, list(parallelize = TRUE, ncores = 2L)))
      expect_equal(parallel$coefficients, serial$coefficients, tolerance = 1e-12)
    }
  }
})

test_that("non-binomial intercepts follow nonnegative for other GLM families", {
  X <- matrix(rep(c(-1, 1), 20), ncol = 1)
  for (family in list(poisson("log"), Gamma("log"), inverse.gaussian("log"))) {
    y <- exp(-1 + 0.4 * X[, 1])
    # Poisson responses are counts; the other families accept positive real data.
    if (family$family == "poisson") y <- rep(c(0, 1), 20)
    for (nonnegative in c(FALSE, TRUE)) {
      fit <- L0glm.fit(
        X, y, family = family, intercept = TRUE, nonnegative = nonnegative,
        lambda = 1e-10, dispersion = 1, start = c(-1, 0.4), warm.start = 1,
        parallelize = FALSE, maxit = 60
      )
      beta <- as.numeric(fit$coefficients)
      expect_true(all(is.finite(beta)))
      if (nonnegative) expect_true(all(beta >= 0)) else expect_lt(beta[1], 0)
    }
  }
})

test_that("L0nnpois enforces nonnegativity on intercepts, starts and refits", {
  X <- matrix(rep(c(0, 1), 20), ncol = 1)
  y <- 2 + 3 * X[, 1]
  # The positive offset makes the unrestricted intercept negative while all means are valid.
  for (sparse in c(FALSE, TRUE)) {
    for (normalize in c(FALSE, TRUE)) {
      for (multichannel in c(FALSE, TRUE)) {
        for (relaxed in c(FALSE, TRUE)) {
          args <- list(
            X = Matrix::Matrix(X, sparse = sparse),
            y = if (multichannel) cbind(y, y + 1) else y,
            offset = rep(4, length(y)), intercept = TRUE, normalize = normalize,
            lambda = 1e-8, start = c(-1, 2), warm.start = 1,
            grouped = TRUE, group.maxit = 2, relaxed = relaxed, maxit = 40,
            parallelize = multichannel, ncores = 2L
          )
          fast <- do.call(L0nnpois, args)
          general <- do.call(L0glm.fit, c(args, list(family = poisson(identity),
                                                   nonnegative = TRUE, openmp = FALSE)))
          expect_true(all(is.finite(fast$coefficients)))
          expect_true(all(fast$coefficients >= 0))
          expect_true(all(as.matrix(fast$coefficients)[1, ] == 0))
          expect_equal(unname(as.matrix(fast$coefficients)),
                       unname(as.matrix(general$coefficients)), tolerance = 1e-8)
        }
      }
    }
  }
  unrestricted <- L0glm.fit(
    X, y, family = poisson(identity), offset = rep(4, length(y)),
    intercept = TRUE, nonnegative = FALSE, lambda = 1e-12, start = c(-1, 2),
    warm.start = 1, parallelize = FALSE, maxit = 50
  )
  expect_equal(as.numeric(unrestricted$coefficients), c(-2, 3), tolerance = 1e-6)
})
