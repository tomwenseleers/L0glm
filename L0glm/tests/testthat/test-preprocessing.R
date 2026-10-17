test_that("column information uses stable weighted dense and sparse moments", {
  set.seed(1630)
  X <- matrix(rnorm(40 * 5), 40, 5)
  X[abs(X) < 1] <- 0
  w <- runif(40); w[1:4] <- 0
  center <- colSums(w * X) / sum(w)
  for (sparse in c(FALSE, TRUE)) {
    Z <- Matrix::Matrix(X, sparse = sparse)
    moments <- L0glm_moments(Z, w)
    expect_equal(moments$mean, center)
    expect_equal(moments$ss, colSums(w * X^2))
    expect_equal(moments$centered_ss, colSums(w * sweep(X, 2, center)^2))
    expect_equal(L0glm_moments(Z, rep(0, 40))$centered_ss, rep(0, 5))
    expect_identical(inherits(Z, "sparseMatrix"), sparse)
  }
  X <- cbind(1e8 + seq(-1, 1, length.out = 40), rep(1e8, 40))
  expect_equal(L0glm_moments(X)$centered_ss,
               colSums(sweep(X, 2, colMeans(X))^2), tolerance = 1e-12)
  expect_equal(L0glm_moments(matrix(1:20, 10))$ss, colSums(matrix(1:20, 10)^2))
  constant <- matrix(0.123456789, 242, 1)
  for (sparse in c(FALSE, TRUE)) {
    moments <- L0glm_moments(Matrix::Matrix(constant, sparse = sparse), rep(0.3, 242))
    expect_identical(moments$mean, 0.123456789)
    expect_identical(moments$centered_ss, 0)
  }
})

test_that("a constant training predictor cannot create enormous held-out estimates", {
  set.seed(1634)
  X <- cbind(rep(0.123456789, 60), rnorm(60), rnorm(60))
  y <- 1 + 2 * X[, 2] + rnorm(60)
  for (sparse in c(FALSE, TRUE)) {
    fit <- L0glm.fit(Matrix::Matrix(X, sparse = sparse), y, family = gaussian(),
      intercept = TRUE, parallelize = FALSE, lambda = "gic")
    expect_equal(as.numeric(fit$coefficients)[2], 0)
    expect_equal(fit$preprocessing$scale[2], 1)
    validation <- X; validation[, 1] <- seq(-10, 10, length.out = 60)
    prediction <- as.numeric(cbind(1, validation) %*% fit$coefficients)
    expect_true(all(is.finite(prediction)))
    expect_lt(max(abs(prediction)), 20)
  }
})

test_that("nearly constant free-intercept columns have a stable working origin", {
  set.seed(1635)
  X <- cbind(0.004 + 1e-8 * rnorm(100), rnorm(100), rnorm(100))
  y <- 1 + 2e8 * (X[, 1] - 0.004) + rnorm(100, sd = 0.2)
  fits <- list()
  for (sparse in c(FALSE, TRUE)) {
    Z <- Matrix::Matrix(X, sparse = sparse)
    expect_warning(fit <- L0glm.fit(Z, y, family = gaussian(), intercept = TRUE,
      lambda = "gic", parallelize = FALSE), NA)
    expect_gt(abs(as.numeric(fit$coefficients)[2]), 1e8)
    expect_lt(mean((y - as.numeric(fit$X %*% fit$coefficients))^2), 0.1)
    expect_identical(inherits(fit$X, "sparseMatrix"), sparse)
    prep <- L0glm_preprocessing(cbind(Intercept = 1, Z), "auto", rep(-Inf, 4), 1L)
    working <- L0glm_center_working_design(cbind(Intercept = 1, Z), prep, rep(-Inf, 4))
    expect_identical(inherits(working$X, "sparseMatrix"), sparse)
    if (sparse) expect_equal(length(working$X@x), length(cbind(Intercept = 1, Z)@x))
    fits[[length(fits) + 1L]] <- fit
  }
  expect_equal(as.numeric(fits[[1]]$coefficients), as.numeric(fits[[2]]$coefficients),
               tolerance = 1e-6)
})

test_that("implicit centering requires a free intercept, not signed slopes", {
  X <- cbind(Intercept = 1, x = 1:20, z = (1:20)^2)
  free <- L0glm_preprocessing(X, "auto", c(-Inf, 0, 0), 1L)
  constrained <- L0glm_preprocessing(X, "rms", rep(0, 3), 1L)
  absent <- L0glm_preprocessing(X[, -1], "rms", rep(-Inf, 2), integer())
  expect_true(free$centered)
  expect_identical(free$mode, "rms")
  expect_equal(free$scale[-1], unname(sqrt(colMeans(scale(X[, -1], scale = FALSE)^2))))
  expect_equal(free$center[-1], unname(colMeans(X[, -1])))
  expect_false(constrained$centered)
  expect_false(absent$centered)
  expect_equal(constrained$scale[-1], unname(sqrt(colMeans(X[, -1]^2))))
  expect_equal(constrained$center, rep(0, 3))
  expect_identical(L0glm_preprocessing(X, "auto", rep(0, 3), 1L)$mode, "none")
  named_nuisance <- cbind(Intercept_like = 1:20, z = (1:20)^2)
  expect_false(L0glm_preprocessing(named_nuisance, "auto", rep(-Inf, 2), 1L)$centered)
})

test_that("free-intercept automatic Gaussian fits are invariant to predictor units and origin", {
  set.seed(1631)
  X <- matrix(rnorm(100 * 25), 100, 25)
  beta <- c(3, -2, 1.5, rep(0, 22))
  y <- 2 + as.numeric(X %*% beta) + rnorm(100, sd = 0.6)
  units <- exp(seq(-2, 2, length.out = ncol(X)))
  origin <- seq(-3, 3, length.out = ncol(X))
  args <- list(y = y, family = gaussian(), intercept = TRUE, maxit = 100,
               parallelize = FALSE, relaxed = FALSE, lambda = "gic")
  for (sparse in c(FALSE, TRUE)) {
    original <- do.call(L0glm.fit, c(list(X = Matrix::Matrix(X, sparse = sparse)), args))
    Z <- sweep(sweep(X, 2, origin, "-"), 2, units, "/")
    transformed <- do.call(L0glm.fit, c(list(X = Matrix::Matrix(Z, sparse = sparse)), args))
    recovered <- as.numeric(transformed$coefficients)
    recovered[-1] <- recovered[-1] / units
    recovered[1] <- recovered[1] - sum(origin * recovered[-1])
    expect_equal(recovered, as.numeric(original$coefficients), tolerance = 1e-7)
    expect_equal(transformed$dispersion, original$dispersion, tolerance = 1e-7)
    expect_equal(as.numeric(transformed$X %*% transformed$coefficients),
                 as.numeric(original$X %*% original$coefficients), tolerance = 1e-7)
    expect_identical(inherits(original$X, "sparseMatrix"), sparse)
    expect_s4_class(original$coefficients, "dgCMatrix")
    shifted <- do.call(L0glm.fit, c(list(X = Matrix::Matrix(X, sparse = sparse)),
                                   modifyList(args, list(y = y + 100))))
    expected <- as.numeric(original$coefficients); expected[1] <- expected[1] + 100
    expect_equal(as.numeric(shifted$coefficients), expected, tolerance = 1e-7)
  }
})

test_that("centering preserves nonnegative slopes with a free intercept only", {
  x <- 3:22; y <- 2 * x - 2
  X <- cbind(Intercept = 1, x = x)
  centered <- cbind(Intercept = 1, x = x - mean(x))
  for (solver in c("eigen", "osqp", "bigGlm")) {
    control <- control.solver.gen(solver, abstol = 1e-16, reltol = 1e-11)
    a <- lmridge(X, y, lambdas = c(0, 1), lower = c(-Inf, 0), solver = solver,
                 control.solver = control)
    b <- lmridge(centered, y, lambdas = c(0, 1), lower = c(-Inf, 0), solver = solver,
                 control.solver = control)
    b[1] <- b[1] - mean(x) * b[2]
    expect_equal(b, a, tolerance = 1e-6)
  }
  # Imposing the same nonnegative intercept bound after centering changes the
  # feasible set: this back-transformed fit has a forbidden negative intercept.
  a <- lmridge(X, y, lambdas = 0, lower = c(0, 0), solver = "osqp")
  b <- lmridge(centered, y, lambdas = 0, lower = c(0, 0), solver = "osqp")
  b[1] <- b[1] - mean(x) * b[2]
  expect_gte(a[1], -1e-7)
  expect_lt(b[1], -1)
  expect_gt(max(abs(X %*% a - X %*% b)), 0.1)
})

test_that("nonnegative binomial slopes permit centered initialization", {
  set.seed(1632)
  X <- matrix(rnorm(200 * 4), 200, 4)
  y <- rbinom(200, 1, plogis(-1 + 1.5 * X[, 1]))
  args <- list(y = y, family = binomial(), nonnegative = TRUE, intercept = TRUE,
               lambda = 2, parallelize = FALSE, maxit = 100, relaxed = FALSE)
  a <- do.call(L0glm.fit, c(list(X = X), args))
  b <- do.call(L0glm.fit, c(list(X = sweep(X, 2, c(2, -1, 3, -2), "-")), args))
  recovered <- as.numeric(b$coefficients)
  recovered[1] <- recovered[1] - sum(c(2, -1, 3, -2) * recovered[-1])
  expect_true(a$preprocessing$centered)
  expect_true(all(a$coefficients[-1, ] >= 0))
  expect_equal(recovered, as.numeric(a$coefficients), tolerance = 1e-6)
})

test_that("free intercepts are reprofiled after clipping or thresholding slopes", {
  set.seed(1636)
  X <- matrix(rbinom(120 * 5, 1, 0.3), 120, 5)
  y <- rbinom(120, 1, plogis(-1 + 3 * X[, 1] - 3 * X[, 2]))
  args <- list(y = y, family = binomial(), intercept = TRUE, nonnegative = TRUE,
    start = 1, lambda = 2, maxit = 1L, parallelize = FALSE, relaxed = FALSE)
  dense <- do.call(L0glm.fit, c(list(X = Matrix::Matrix(X, sparse = FALSE)), args))
  sparse <- do.call(L0glm.fit, c(list(X = Matrix::Matrix(X, sparse = TRUE)), args))
  expect_equal(as.numeric(dense$coefficients), as.numeric(sparse$coefficients), tolerance = 1e-7)
  expect_true(any(dense$coefficients[-1, ] == 0))
  expect_true(all(dense$coefficients[-1, ] >= 0))
  expect_identical(inherits(sparse$X, "sparseMatrix"), TRUE)
})

test_that("canonical Poisson initialization agrees between native and R paths", {
  set.seed(1633)
  X <- matrix(runif(50 * 12), 50, 12)
  y <- cbind(rpois(50, 2 + 3 * X[, 1]), rpois(50, 1 + 2 * X[, 2]))
  for (sparse in c(FALSE, TRUE)) for (mode in c("auto", "rms", "l2", "sd")) {
    args <- list(X = Matrix::Matrix(X, sparse = sparse), y = y, intercept = TRUE,
      family = poisson(identity),
      normalize = mode, relaxed = TRUE, grouped = TRUE, group.maxit = 2,
      parallelize = FALSE, maxit = 30, lambda = 2)
    a <- do.call(L0nnpois, args)
    b <- do.call(L0glm.fit, c(args, list(nonnegative = TRUE, openmp = FALSE)))
    expect_equal(a$coefficients, b$coefficients, tolerance = 1e-7)
    expect_equal(a$preprocessing, b$preprocessing)
    expect_false(a$preprocessing$centered)
    expect_identical(inherits(a$X, "sparseMatrix"), sparse)
  }
})
