test_that("Cox CV uses full-risk-set partial-likelihood subtraction", {
  skip_if_not_installed("mgcv"); skip_if_not_installed("survival")
  set.seed(1812)
  X <- matrix(rnorm(90 * 4), 90, 4)
  y <- survival::Surv(rexp(90, exp(.5 * X[, 1] - .4 * X[, 2])), rep(c(1, 1, 0), 30))
  folds <- rep(1:3, each = 30)
  result <- L0glm.tune(X, y, family = mgcv::cox.ph(), tune = "cv", foldid = folds,
    maxit = 100, tol = 1e-7, nlambda = 3, plot = FALSE, printbest = FALSE)
  expect_true(all(is.finite(result$cv$deviance)))
  expect_true(result$endpoint_zero)
  train <- folds != 1
  train_args <- result$fit_args
  train_args$weights <- train_args$weights[train, , drop = FALSE]
  train_args$offset <- train_args$offset[train, , drop = FALSE]
  fit <- do.call(L0glm.fit, c(list(X = X[train, ], y = y[train, ],
    lambda = result$cv$lambda_train[1, 1]), train_args))
  full <- L0glm.metrics(fit, X, y)$loglik
  expected <- -2 * (full - as.numeric(logLik(fit)))/sum(y[!train, 2])
  expect_equal(result$cv$deviance[1, 1], expected, tolerance = 1e-10)
})
