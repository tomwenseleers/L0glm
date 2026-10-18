test_that("multinomial row exclusions precede CLR centering and keep mandatory rows", {
  skip_if_not_installed("mgcv")
  set.seed(2101)
  X <- matrix(rnorm(180 * 5), 180, 5)
  y <- factor(rep(1:3, 60))
  initial <- matrix(1, 6, 3)
  initial[c(1, 4, 6), ] <- 0
  initial[2, 1] <- 0
  data <- L0glm_ext_prepare(X, y, mgcv::multinom(K = 2), NULL, NULL,
    TRUE, 6L, FALSE, "rms", initial, NULL, list())
  expect_true(all(!data$allowed[4, ]))
  expect_true(all(data$allowed[-4, ]))
  expect_true(all(data$B[4, ] == 0))
  expect_true(all(data$B[3, ] == 0)) # A nonzero common start centers to zero, but remains eligible.
  expect_true(all(data$allowed[2, ])) # An individual category zero is not a constraint.

  original <- L0glm_ext_multinomial_working
  calls <- list()
  fit <- testthat::with_mocked_bindings(
    L0glm.fit(X, y, family = mgcv::multinom(K = 2), start = initial,
      always.include = 6L, lambda = 0, relaxed = TRUE, maxit = 100L, tol = 1e-8),
    L0glm_ext_multinomial_working = function(data, B, D, state, args) {
      calls[[length(calls) + 1L]] <<- D
      original(data, B, D, state, args)
    }, .package = "L0glm")
  expect_true(all(vapply(calls, function(D) all(is.infinite(D[4, ])), logical(1))))
  expect_equal(dim(fit$coefficients), c(6L, 3L))
  expect_s4_class(fit$coefficients, "sparseMatrix")
  expect_true(all(fit$coefficients[4, ] == 0))
  expect_lt(max(abs(Matrix::rowSums(fit$coefficients))), 1e-12)
})

test_that("masked multinomial fits match fitting the retained predictors directly", {
  skip_if_not_installed("mgcv")
  set.seed(2102)
  X <- matrix(rnorm(260 * 5), 260, 5)
  colnames(X) <- paste0("x", 1:5)
  eta <- cbind(.2 + .7*X[, 1], -.1 - .6*X[, 1] + .5*X[, 3], -.1 - .5*X[, 3])
  P <- L0glm_ext_softmax(eta)$p
  y <- factor(vapply(seq_len(nrow(X)), function(i) sample.int(3, 1, prob = P[i, ]), integer(1)))
  keep <- c(1, 3, 5)
  initial <- matrix(1, 6, 3); initial[c(3, 5), ] <- 0
  for (sparse in c(FALSE, TRUE)) for (grouped in c(FALSE, TRUE)) for (relaxed in c(FALSE, TRUE)) {
    design <- if (sparse) Matrix::Matrix(X, sparse = TRUE) else X
    fit <- L0glm.fit(design, y, family = mgcv::multinom(K = 2), start = initial,
      lambda = 2, grouped = grouped, relaxed = relaxed, maxit = 200L, tol = 1e-8)
    subset <- L0glm.fit(design[, keep, drop = FALSE], y, family = mgcv::multinom(K = 2),
      start = 1, lambda = 2, grouped = grouped, relaxed = relaxed, maxit = 200L, tol = 1e-8)
    expect_true(all(fit$convergence$converged))
    expect_s4_class(fit$coefficients, "sparseMatrix")
    expect_equal(dim(fit$coefficients), c(6L, 3L))
    expect_true(all(fit$coefficients[c(3, 5), ] == 0))
    expect_equal(as.matrix(fit$coefficients[c(1, keep + 1L), ]),
      as.matrix(subset$coefficients), tolerance = 1e-7)
    expect_equal(predict(fit), predict(subset), tolerance = 1e-7)
    expect_equal(fit$objective, subset$objective, tolerance = 1e-7)
  }
})

test_that("zero multinomial starts restrict support but unit CLR starts remain eligible", {
  skip_if_not_installed("mgcv")
  set.seed(2103)
  X <- matrix(rnorm(240 * 4), 240, 4)
  P <- L0glm_ext_softmax(cbind(.8*X[, 1], -.7*X[, 1], -.1*X[, 1]))$p
  y <- factor(vapply(seq_len(nrow(X)), function(i) sample.int(3, 1, prob = P[i, ]), integer(1)))
  zero <- L0glm.fit(X, y, family = mgcv::multinom(K = 2), start = 0,
    lambda = 0, relaxed = TRUE, maxit = 100L, tol = 1e-9)
  expect_equal(dim(zero$coefficients), c(5L, 3L))
  expect_true(all(zero$coefficients[-1, ] == 0))
  expected <- matrix(rep(as.numeric(table(y))/length(y), each = nrow(X)), nrow(X))
  expect_equal(unname(predict(zero)), expected, tolerance = 1e-8)
  for (start in list(1, matrix(1, 5, 3))) {
    fit <- L0glm.fit(X, y, family = mgcv::multinom(K = 2), start = start,
      lambda = 0, relaxed = FALSE, maxit = 100L, tol = 1e-9)
    expect_true(any(fit$coefficients[-1, ] != 0))
    expect_lt(max(abs(Matrix::rowSums(fit$coefficients))), 1e-12)
  }
})

test_that("multinomial CV preserves excluded rows and full coefficient dimensions", {
  skip_if_not_installed("mgcv")
  set.seed(2104)
  X <- matrix(rnorm(180 * 4), 180, 4)
  P <- L0glm_ext_softmax(cbind(.9*X[, 1], -.7*X[, 1], -.2*X[, 1]))$p
  y <- factor(vapply(seq_len(nrow(X)), function(i) sample.int(3, 1, prob = P[i, ]), integer(1)))
  initial <- matrix(1, 5, 3); initial[4, ] <- 0
  original <- L0glm_ext_multinomial_working
  masked <- training <- logical()
  result <- testthat::with_mocked_bindings(
    L0glm.tune(X, y, family = mgcv::multinom(K = 2), start = initial,
      tune = "cv", foldid = rep(1:3, 60), nlambda = 3, grouped = TRUE,
      relaxed = TRUE, maxit = 120L, tol = 1e-7, plot = FALSE, printbest = FALSE),
    L0glm_ext_multinomial_working = function(data, B, D, state, args) {
      masked <<- c(masked, all(is.infinite(D[4, ])))
      training <<- c(training, data$n < nrow(X))
      original(data, B, D, state, args)
    }, .package = "L0glm")
  expect_true(any(training))
  expect_true(all(masked))
  coefficients <- c(result$coefficients,
    list(result$fit_best$coefficients, result$fit_1se$coefficients))
  expect_true(all(vapply(coefficients, function(B) inherits(B, "sparseMatrix"), logical(1))))
  expect_true(all(vapply(coefficients, function(B) identical(dim(B), c(5L, 3L)), logical(1))))
  expect_true(all(vapply(coefficients, function(B) all(B[4, ] == 0), logical(1))))
  expect_true(all(is.finite(result$cv$deviance)))
  replay <- do.call(L0glm.fit,
    c(list(X = X, y = y, lambda = result$lambda_best), result$fit_args))
  expect_equal(replay$coefficients, result$fit_best$coefficients, tolerance = 1e-10)
})
