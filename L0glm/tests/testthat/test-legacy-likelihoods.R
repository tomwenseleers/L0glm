test_that("legacy Gaussian path adapters preserve their weighted scale and metadata", {
  set.seed(1820)
  n <- 30L; p <- 5L
  X <- matrix(rnorm(n * p), n, p)
  beta <- cbind(null = 0, small = c(1, -0.5, 0, 0, 0), full = seq_len(p) / p)
  y <- rnorm(n)
  w <- seq(0.2, 2, length.out = n)
  for (sparse in c(FALSE, TRUE)) for (matrix_weights in c(FALSE, TRUE)) {
    design <- Matrix::Matrix(X, sparse = sparse)
    coefficients <- Matrix::Matrix(beta, sparse = sparse)
    weights <- if (matrix_weights) cbind(w, w / 2, rev(w)) else w
    df <- Matrix::colSums(coefficients != 0)
    prediction <- as.matrix(design %*% coefficients)
    expected <- vapply(seq_len(ncol(beta)), function(j) {
      residual <- sqrt(if (matrix_weights) weights[, j] else weights) * (y - prediction[, j])
      sigma <- sqrt(sum(residual^2) / (n - df[j]))
      sum(stats::dnorm(residual, sd = sigma, log = TRUE))
    }, numeric(1))
    glmnet <- structure(list(beta = coefficients, nobs = n, dim = c(n, p)), class = "glmnet")
    l0learn <- structure(list(beta = list(coefficients), n = n, p = p,
      suppSize = list(df), dim = c(n, p)), class = "L0Learn")
    for (object in list(glmnet, l0learn)) {
      value <- logLik(object, weights = weights, X = design, y = y)
      expect_equal(as.numeric(value), expected, tolerance = 1e-13)
      expect_identical(names(value), colnames(beta))
      expect_identical(attr(value, "df"), df)
      expect_identical(attr(value, "p"), p)
      expect_identical(attr(value, "nobs"), if (inherits(object, "glmnet")) as.numeric(n) else n)
      expect_identical(logLik(object, weights, design, y), value)
      expect_identical(logLik(object, X = design, weights, y = y), value)
    }
  }
})

test_that("the ordinis adapter retains its explicit intercept and coefficient count", {
  set.seed(1821)
  X <- matrix(rnorm(80), 20, 4)
  beta <- cbind(c(2, -1, 0, 0, 0), c(0, 0.5, -0.2, 0, 0))
  y <- rnorm(20); weights <- seq(0.1, 2, length.out = 20)
  for (sparse in c(FALSE, TRUE)) {
    object <- list(beta = Matrix::Matrix(beta, sparse = sparse), nobs = 20L, nvars = 4L)
    value <- logLik2(object, weights, Matrix::Matrix(X, sparse = sparse), y)
    prediction <- cbind(1, X) %*% beta
    df <- colSums(beta != 0)
    expected <- vapply(seq_len(2), function(j) {
      residual <- sqrt(weights) * (y - prediction[, j])
      sum(stats::dnorm(residual, sd = sqrt(sum(residual^2)/(20 - df[j])), log = TRUE))
    }, numeric(1))
    expect_equal(as.numeric(value), expected, tolerance = 1e-13)
    expect_identical(attr(value, "df"), Matrix::colSums(object$beta != 0))
    expect_identical(attr(value, "nobs"), 20)
    expect_identical(attr(value, "p"), 4L)
  }
})

test_that("legacy adapters retain missing-dimension and zero-residual behavior", {
  X <- matrix(1:12, 6, 2)
  beta <- matrix(c(1, 0), 2, 1)
  object <- structure(list(beta = beta, nobs = 6L), class = "glmnet")
  value <- logLik(object, weights = rep(1, 6), X = X, y = X[, 1])
  expect_identical(as.numeric(value), Inf)
  expect_null(attr(value, "p"))
  expect_error(logLik(object, X = X, y = X[, 1]), "require weights, X, and y")
})
