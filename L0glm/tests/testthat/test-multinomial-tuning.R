test_that("multinomial row criterion, pilot seed, masks and plotting work together", {
  skip_if_not_installed("mgcv")
  set.seed(1811)
  X <- matrix(rnorm(150 * 5), 150, 5)
  P <- L0glm_ext_softmax(cbind(0, .8 * X[, 1] - .5 * X[, 2], -.7 * X[, 1] + .6 * X[, 3]))$p
  y <- factor(vapply(seq_len(150), function(i) sample.int(3, 1, prob = P[i, ]), integer(1)))
  result <- L0glm.tune(X, y, family = mgcv::multinom(K = 2), grouped = TRUE,
    relaxed = TRUE, start = "pilot", nlambda = 4, maxit = 100, tol = 1e-7,
    plot = FALSE, printbest = FALSE)
  expect_equal(result$criterion_support, "rows")
  expect_identical(result$coordinates, "clr")
  expect_true(result$endpoint_zero)
  expect_equal(as.numeric(rowSums(result$start)), rep(0, 6), tolerance = 1e-12)
  expect_equal(result$metrics$tuning_score,
    -2 * result$metrics$loglik + result$lambda_preset *
      (result$metrics$selected + result$metrics$mandatory + result$metrics$nuisance))
  replay <- do.call(L0glm.fit, c(list(X = X, y = y, lambda = result$lambda_best), result$fit_args))
  expect_equal(replay$coefficients, result$fit_best$coefficients, tolerance = 1e-10)
  plotfile <- tempfile(fileext = ".pdf")
  grDevices::pdf(plotfile)
  expect_silent(L0glm_plot_tuning(result))
  grDevices::dev.off()
  expect_gt(file.info(plotfile)$size, 1000)
  folds <- rep(1:3, 50)
  cv <- L0glm.tune(X, y, family = mgcv::multinom(K = 2), tune = "cv", foldid = folds,
    nlambda = 3, maxit = 100, tol = 1e-7, plot = FALSE, printbest = FALSE)
  expect_true(all(is.finite(cv$cv$deviance)))
  expect_s3_class(cv$fit_1se, "L0glm_extended_fit")
  expect_equal(rowSums(predict(cv$fit_1se)), rep(1, 150), tolerance = 1e-12)
})
