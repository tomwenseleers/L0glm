test_that("column scaling preserves the weighted ridge objective for every solver", {
  set.seed(1620)
  X <- matrix(rnorm(40 * 8), 40, 8)
  X <- sweep(X, 2, exp(seq(-2, 2, length.out = 8)), "*")
  y <- rnorm(40); w <- runif(40, 0.2, 2)
  d <- seq(0.2, 3, length.out = 8); target <- rnorm(8) / 5
  reference <- as.numeric(solve(crossprod(X, w * X) + diag(d),
                                crossprod(X, w * y) + d * target))
  for (sparse in c(FALSE, TRUE)) for (method in 1:6) {
    estimate <- lmridge(Matrix::Matrix(X, sparse = sparse), y, weights = w,
                        lambdas = d, prior.mean = target, scaling = scaleX(X),
                        start = rep(0.1, 8), method = method,
                        control.solver = control.solver.gen("eigen", method,
                                                            tol = 1e-11))
    expect_equal(estimate, reference, tolerance = 1e-6)
  }
  for (solver in c("osqp", "bigGlm")) {
    estimate <- lmridge(X, y, weights = w, lambdas = d, prior.mean = target,
                        scaling = scaleX(X), solver = solver,
                        control.solver = control.solver.gen(solver, abstol = 1e-14,
                                                            reltol = 1e-11))
    expect_equal(estimate, reference, tolerance = 1e-6)
  }
})

test_that("scaled bounds and ridge targets describe the original coefficients", {
  X <- diag(c(0.5, 1, 4, 2)); y <- 1:4
  d <- c(0.1, 2, 1, 20); target <- c(0.1, -0.2, 0.3, 0)
  lower <- c(0, -0.3, 0, 0); upper <- c(0.8, 0.2, 1, Inf)
  expected <- pmin(upper, pmax(lower, (diag(X) * y + d * target) / (diag(X)^2 + d)))
  for (solver in c("eigen", "osqp", "bigGlm")) {
    beta <- lmridge(X, y, lambdas = d, prior.mean = target, scaling = c(0.5, 1.5, 1, 1),
                    lower = lower, upper = upper, solver = solver)
    expect_equal(beta, expected, tolerance = 1e-6)
  }
  beta <- lmridge(X, y, lambdas = d, prior.mean = target, scaling = scaleX(X),
                  nonnegative = TRUE, solver = "nnls")
  expect_equal(beta, pmax(0, (diag(X) * y + d * target) / (diag(X)^2 + d)), tolerance = 1e-8)
  expect_error(lmridge(X, y, scaling = c(1, 0, 1, 1)), "finite positive")
})

test_that("normalization choices have their stated column norms", {
  X <- cbind(1, 1:10, (-5:4)^2)
  expect_equal(L0glm_column_scale(X, FALSE), rep(1, 3))
  expect_equal(L0glm_column_scale(X, TRUE), sqrt(colSums(X^2)))
  expect_equal(L0glm_column_scale(X, "rms"), sqrt(colMeans(X^2)))
  expect_equal(L0glm_column_scale(X, "sd"), c(1, apply(X[, -1], 2, sd)), ignore_attr = TRUE)
  expect_error(L0glm_normalization("unsupported"), "normalize must")
})

test_that("all scaling choices preserve native and general fit coordinates", {
  set.seed(1621)
  X <- matrix(runif(30 * 45), 30, 45)
  X <- sweep(X, 2, exp(seq(-1, 1, length.out = 45)), "*")
  y <- rpois(30, 3)
  for (sparse in c(FALSE, TRUE)) {
    args <- list(X = Matrix::Matrix(X, sparse = sparse), y = y, lambda = 2,
                 start = 0.2, maxit = 5L, relaxed = FALSE, parallelize = FALSE,
                 prior.mean = rep(0.01, 45))
    reference <- do.call(L0nnpois, args)
    for (normalize in list(TRUE, "l2", "rms", "sd", "none")) {
      native <- do.call(L0nnpois, c(args, list(normalize = normalize)))
      general <- do.call(L0glm.fit, c(args, list(normalize = normalize,
        family = poisson(identity), nonnegative = TRUE, openmp = FALSE)))
      expect_equal(as.numeric(native$coefficients), as.numeric(reference$coefficients), tolerance = 1e-8)
      expect_equal(as.numeric(native$coefficients), as.numeric(general$coefficients), tolerance = 1e-8)
      expect_identical(inherits(native$X, "sparseMatrix"), sparse)
    }
  }
})
