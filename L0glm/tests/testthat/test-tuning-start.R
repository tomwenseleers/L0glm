test_that("a supplied common start replays every lambda independently", {
  X <- diag(6); y <- c(4, 1, 0.2, -0.1, 0.1, 0.2)
  args <- list(X = X, y = y, family = gaussian(), dispersion = 1,
    start = c(1, 0, 1, 1, 0, 1), parallelize = FALSE, relaxed = FALSE)
  tuned <- do.call(L0glm.tune, c(args, list(nlambda = 4, lambda_min = 0.01,
    plot = FALSE, printbest = FALSE)))
  expect_equal(tuned$start, args$start)
  expect_null(tuned$pilot)
  expect_equal(tuned$p, 6L)
  for (i in rev(seq_along(tuned$lambdas))) {
    direct <- do.call(L0glm.fit, c(list(X = X, y = y, lambda = tuned$lambdas[i]),
                                   tuned$fit_args))
    expect_equal(direct$coefficients, tuned$coefficients[[i]])
    expect_true(all(direct$coefficients[c(2, 5), ] == 0))
  }
  replay <- do.call(L0glm.fit, c(list(X = X, y = y, lambda = tuned$lambda_best),
                                 tuned$fit_args))
  expect_equal(replay, tuned$fit_best)
})

test_that("pilot pooling retains opposite signed features and structural zeros", {
  calls <- 0L
  local_mocked_bindings(L0glm.fit = function(..., lambda) {
    calls <<- calls + 1L
    beta <- matrix(c(if (calls %% 2L) 2 else -2, 0, 0, 1), 2, 2)
    list(coefficients = beta, always.include = integer())
  })
  result <- L0glm_pilot_start(list(), c(2, 20))
  expect_equal(result$start, matrix(c(2, 0, 0, 1), 2, 2))
  expect_equal(nrow(result$pilot), 10L)
  expect_identical(result$pilot$lambda[c(1, 10)], c(2, 20))
  expect_identical(calls, 10L)
  single <- L0glm_pilot_start(list(), c(2, 2))
  expect_equal(nrow(single$pilot), 1L)
  expect_equal(single$pilot$lambda, 2)
})

test_that("pilot starts match independent fits and reproduce the full tuning path", {
  set.seed(1650)
  X <- matrix(rnorm(70 * 10), 70, 10)
  y <- 1 + 2 * X[, 1] - X[, 3] + rnorm(70)
  args <- list(X = X, y = y, family = gaussian(), intercept = TRUE,
    parallelize = FALSE, relaxed = FALSE)
  tuned <- do.call(L0glm.tune, c(args, list(start = "pilot", nlambda = 4,
    plot = FALSE, printbest = FALSE)))
  coefficients <- vapply(tuned$pilot$lambda, function(lambda)
    as.numeric(do.call(L0glm.fit, c(args, list(lambda = lambda)))$coefficients),
    numeric(11))
  direction <- sign(rowMeans(coefficients))
  cancel <- which(direction == 0)
  if (length(cancel)) direction[cancel] <- vapply(cancel, function(j)
    sign(coefficients[j, which.max(abs(coefficients[j, ]))]), numeric(1))
  expect_equal(as.numeric(tuned$start), direction * rowMeans(abs(coefficients)))
  expect_equal(tuned$p, 10L)
  expect_equal(nrow(tuned$start), 11L)
  expect_true(tuned$endpoint_zero)
  for (i in rev(seq_along(tuned$lambdas))) {
    direct <- do.call(L0glm.fit, c(list(X = X, y = y, lambda = tuned$lambdas[i]),
                                   tuned$fit_args))
    expect_equal(direct$coefficients, tuned$coefficients[[i]])
  }
  expect_equal(tuned$preprocessing, tuned$fit_best$preprocessing)
})

test_that("CV builds pilot starts and scales using only each fold's training rows", {
  set.seed(1651)
  n <- 30L
  X <- matrix(rnorm(n * 8), n, 8, dimnames = list(as.character(seq_len(n)), NULL))
  y <- X[, 1] + rnorm(n)
  folds <- rep(1:3, length.out = n)
  observed <- list()
  original <- L0glm.fit
  local_mocked_bindings(L0glm.fit = function(X, y, ...) {
    args <- list(...)
    if (nrow(X) < n) observed[[length(observed) + 1L]] <<-
      list(rows = as.integer(rownames(X)), start = args$start, lambda = args$lambda)
    original(X, y, ...)
  })
  tuned <- L0glm.tune(X, y, family = gaussian(), start = "pilot", tune = "cv",
    foldid = folds, nlambda = 3L, parallelize = FALSE, plot = FALSE, printbest = FALSE)
  expect_length(observed, 3L * 13L)
  for (k in 1:3) {
    calls <- observed[(k - 1L) * 13L + 1:13]
    for (call in calls) expect_identical(call$rows, which(folds != k))
    for (call in calls[1:10]) expect_null(call$start)
    for (call in calls[11:13]) expect_equal(call$start, tuned$cv$starts[[k]])
    expect_equal(tuned$cv$pilots[[k]]$lambda,
      tuned$pilot$lambda * 2 / 3)
    scale <- L0glm_column_scale(X[folds != k, , drop = FALSE], "rms")
    expect_equal(tuned$cv$preprocessing[[k]]$scale, scale)
  }
  replay <- do.call(original, c(list(X = X, y = y, lambda = tuned$lambda_1se),
                                tuned$fit_args))
  expect_equal(replay, tuned$fit_1se)
})

test_that("multichannel common starts preserve each channel's restrictions", {
  X <- diag(6)
  y <- cbind(c(8, 0, 4, 1, 0, 0), c(0, 6, 2, 0, 0, 1))
  start <- matrix(1, 6, 2); start[c(2, 5), 1] <- 0; start[c(1, 5), 2] <- 0
  for (sparse in c(FALSE, TRUE)) {
    tuned <- L0glm.tune(Matrix::Matrix(X, sparse = sparse), y,
      family = poisson(identity), nonnegative = TRUE, grouped = TRUE,
      group.maxit = 2L, start = start, relaxed = TRUE, nlambda = 3L,
      plot = FALSE, printbest = FALSE, parallelize = FALSE)
    for (beta in tuned$coefficients) {
      expect_s4_class(beta, "dgCMatrix")
      expect_equal(dim(beta), c(6L, 2L))
      expect_true(all(beta[start == 0] == 0))
    }
  }
})
