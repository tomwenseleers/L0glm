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

test_that("internal scalar starts expand over the full design", {
  for (sparse in c(FALSE, TRUE)) {
    for (intercept in c(FALSE, TRUE)) {
      X <- rbind(diag(3), diag(3))
      if (intercept) X <- cbind(Intercept = 1, X)
      X <- Matrix::Matrix(X, sparse = sparse)
      args <- list(X = X, y = seq_len(6), family = gaussian(),
                   dispersion = 1, lambda = 4, maxit = 1, thresh = 0)
      for (start in c(1, -2)) {
        scalar <- do.call(L0glm.fit.internal, c(args, list(start = start)))
        vector <- do.call(L0glm.fit.internal,
                          c(args, list(start = rep(start, ncol(X)))))
        expect_identical(scalar, vector)
        expect_equal(as.numeric(scalar$lambdas), rep(1 / start^2, ncol(X)))
      }
    }
  }
})

test_that("scalar starts work for signed single and multichannel Gaussian fits", {
  set.seed(16011)
  X0 <- matrix(rnorm(36 * 4), 36, 4)
  for (sparse in c(FALSE, TRUE)) {
    X <- Matrix::Matrix(X0, sparse = sparse)
    for (intercept in c(FALSE, TRUE)) {
      p <- ncol(X) + as.integer(intercept)
      for (m in c(1L, 3L)) {
        y <- matrix(rnorm(36 * m), 36, m)
        args <- list(X = X, y = y, family = gaussian(), dispersion = 1,
                     intercept = intercept, nonnegative = FALSE,
                     lambda = 4, maxit = 1, thresh = 0,
                     grouped = FALSE, relaxed = FALSE,
                     parallelize = FALSE, openmp = FALSE)
        for (start in c(1, -2)) {
          scalar <- do.call(L0glm.fit, c(args, list(start = start)))
          explicit <- do.call(L0glm.fit,
                              c(args, list(start = matrix(start, p, m))))
          expect_identical(scalar$coefficients, explicit$coefficients)
          expect_identical(scalar$lambdas, explicit$lambdas)
        }
      }
    }
  }
})

test_that("scalar Poisson starts cover intercepts, grouping and both backends", {
  set.seed(16012)
  for (p0 in c(1L, 4L)) {
    X0 <- matrix(runif(30 * p0), 30, p0)
    y <- matrix(rpois(90, 3), 30, 3)
    for (sparse in c(FALSE, TRUE)) {
      X <- Matrix::Matrix(X0, sparse = sparse)
      for (intercept in c(FALSE, TRUE)) {
        p <- p0 + as.integer(intercept)
        for (grouped in c(FALSE, TRUE)) {
          args <- list(X = X, y = y, family = poisson(identity),
                       intercept = intercept, lambda = 0.5, maxit = 3,
                       group.maxit = 2, grouped = grouped, relaxed = grouped,
                       parallelize = FALSE)
          for (backend in c("general", "native", "dispatch")) {
            fitter <- if (backend == "native") L0nnpois else L0glm.fit
            controls <- if (backend == "native") list() else
              list(nonnegative = TRUE, openmp = backend == "dispatch")
            scalar <- do.call(fitter, c(args, controls, list(start = 1)))
            explicit <- do.call(fitter,
              c(args, controls, list(start = matrix(1, p, ncol(y)))))
            expect_identical(scalar$coefficients, explicit$coefficients)
            expect_identical(scalar$lambdas, explicit$lambdas)
          }
        }
      }
    }
  }
})

test_that("scalar starts survive multichannel parallel job splitting", {
  cores <- parallel::detectCores()
  skip_if(is.na(cores) || cores < 2L, "parallel test requires two cores")
  set.seed(16013)
  X <- matrix(runif(30 * 4), 30, 4)
  y <- matrix(rpois(90, 3), 30, 3)
  args <- list(X = X, y = y, family = poisson(identity),
               intercept = TRUE, lambda = 0.5, maxit = 3,
               grouped = TRUE, group.maxit = 2, relaxed = TRUE, ncores = 2L)
  for (backend in c("general", "native")) {
    fitter <- if (backend == "native") L0nnpois else L0glm.fit
    controls <- if (backend == "native") list() else
      list(nonnegative = TRUE, openmp = FALSE)
    serial <- do.call(fitter,
      c(args, controls, list(start = matrix(1, 5, 3), parallelize = FALSE)))
    parallel <- do.call(fitter,
      c(args, controls, list(start = 1, parallelize = TRUE)))
    expect_identical(parallel$coefficients, serial$coefficients)
    expect_identical(parallel$lambdas, serial$lambdas)
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
