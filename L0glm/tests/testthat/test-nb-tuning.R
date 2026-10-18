test_that("NB CV returns reproducible full-data best and one-SE fits", {
  skip_if_not_installed("mgcv")
  set.seed(1810)
  X <- matrix(rnorm(120 * 6), 120, 6)
  y <- rnbinom(120, mu = exp(.4 + X[, 1] * .7 - X[, 2] * .5), size = 3)
  folds <- rep(1:3, 40)
  args <- list(X = X, y = y, family = mgcv::nb(), intercept = TRUE,
    grouped = FALSE, relaxed = TRUE, maxit = 100, tol = 1e-7,
    nlambda = 4, plot = FALSE, printbest = FALSE, tune = "cv", foldid = folds)
  result <- do.call(L0glm.tune, args)
  expect_identical(result$cv$foldid, as.integer(folds))
  expect_true(result$endpoint_zero)
  expect_equal(result$fit_best$init.lambda, result$lambda_best)
  expect_equal(result$fit_1se$init.lambda, result$lambda_1se)
  expect_gte(result$lambda_1se, result$lambda_best)
  b <- which.min(result$cv$mean)
  expect_equal(result$lambda_1se, max(result$lambdas[result$cv$mean <= result$cv$mean[b] + result$cv$se[b]]))
  expect_equal(result$cv$lambda_train, outer(result$lambdas, rep(2/3, 3)))
  refit <- do.call(L0glm.fit, c(list(X = X, y = y, lambda = result$lambda_best), result$fit_args))
  expect_equal(refit$coefficients, result$fit_best$coefficients, tolerance = 1e-10)
  expect_equal(refit$theta, result$fit_best$theta, tolerance = 1e-10)
  train <- folds != 1
  train_args <- result$fit_args
  train_args$weights <- train_args$weights[train, , drop = FALSE]
  train_args$offset <- train_args$offset[train, , drop = FALSE]
  ref <- do.call(L0glm.fit, c(list(X = X[train, ], y = y[train],
    lambda = result$cv$lambda_train[1, 1]), train_args))
  expect_equal(ref$theta, result$cv$theta[[1]][[1]], tolerance = 1e-10)
  expect_equal(ref$preprocessing, result$cv$preprocessing[[1]])
})
