# Run from the repository root:
# Rscript benchmarks/tests/test-readme-deconvolution.R
options(L0glm.readme.run = FALSE)
source("benchmarks/readme_benchmarks_L0glm_deconvolution.R")
library(testthat)

test_that("manual logical flags accept logicals and unambiguous strings", {
  expect_true(readme_flag(TRUE, "run_glmnet"))
  expect_false(readme_flag(FALSE, "run_glmnet"))
  expect_true(readme_flag("true", "run_glmnet"))
  expect_false(readme_flag(" FALSE ", "openmp"))
  for (bad in list(NA, 1, NULL, "yes", c(TRUE, FALSE))) {
    expect_error(readme_flag(bad, "run_glmnet"), "must be TRUE or FALSE")
  }
})

test_that("benchmark support metrics use rows and agree with binary correlation", {
  truth <- cbind(c(0, 2, 1, 0, 0), c(3, 0, 0, 0, 0))
  estimate <- cbind(c(0, 0, 5, 2, 0), c(1, 0, 0, 0, 0))
  stats <- support_statistics(estimate, truth)
  expect_equal(unlist(stats[c("TP", "FP", "FN", "TN")], use.names = FALSE), c(2, 1, 1, 1))
  expect_equal(stats$F1, 2 / 3)
  expect_equal(stats$MCC, cor(as.numeric(stats$selected), as.numeric(stats$truth)))
  expect_equal(stats$MCC, L0glm:::support_mcc(2, 1, 1, 1))
  expect_equal(support_statistics(truth, truth)$MCC, 1)
  expect_true(is.na(support_statistics(0 * truth, truth)$MCC))
  expect_equal(support_statistics(0 * truth, truth)$F1, 0)
  # Several sub-threshold channels must not collectively select a row.
  expect_false(any(support_statistics(matrix(0.9e-5, 5, 2), truth)$selected))
})

test_that("coefficient R2 measures individual channel entries", {
  truth <- cbind(c(1, 2, 0), c(0, 0, 3))
  estimate <- truth[, 2:1]
  expect_equal(rowMeans(estimate), rowMeans(truth))
  expect_equal(coefficient_r2(estimate, truth), 0)
  expect_equal(coefficient_r2(2 * truth, truth), 1)
  expect_equal(coefficient_r2(estimate + truth / 3, truth),
               summary(lm(as.vector(truth) ~ 0 + as.vector(estimate + truth / 3)))$r.squared)
  expect_equal(coefficient_r2(0 * truth, truth), 0)
  expect_true(is.na(coefficient_r2(truth, 0 * truth)))
})

test_that("scores and plots handle single methods, zero fits and empty results", {
  truth <- cbind(c(0, 4, 1, 0, 0), c(3, 0, 0, 0, 0))
  sim <- list(X = Matrix::Diagonal(5), y = truth, y_true = truth, beta_true = truth)
  row <- score_fit("multi", "truth", 0, truth, sim)
  expect_equal(row$MCC, 1)
  expect_equal(row$F1, 1)
  expect_equal(row$R2_coefficients, 1)
  blank <- score_fit("multi", "zero", 0, 0 * truth, sim)
  expect_true(is.na(blank$MCC))
  expect_equal(blank$F1, 0)
  directory <- tempfile("readme-plot-tests-")
  dir.create(directory)
  files <- file.path(directory, paste0(c("support", "coefficients", "empty", "summary", "signal"), ".png"))
  plot_multi_support(sim, list(zero = 0 * truth), files[1])
  plot_coefficients("Multichannel", sim, list(zero = 0 * truth), files[2])
  plot_summary(data.frame(seconds = NA_real_), files[3])
  plot_summary(blank, files[4])
  single <- list(X = sim$X, y = truth[, 1], y_true = truth[, 1], beta_true = truth[, 1])
  plot_single_signal(single, list(zero = 0 * truth[, 1]), files[5])
  expect_true(all(file.exists(files)))
  expect_true(all(file.info(files)$size > 1000))
})

test_that("Poisson deviance uses the noise-free response without rounding", {
  truth <- c(0, 0.25, 1.5, 9)
  fitted <- c(0.2, 0.5, 1, 8)
  expected <- 2 * sum(fitted - truth + ifelse(truth == 0, 0, truth * log(truth / fitted)))
  expect_equal(poisson_deviance(truth, fitted), expected)
  expect_equal(poisson_deviance(truth[-1], truth[-1]), 0)
  expect_gt(poisson_deviance(truth, round(truth)), 0)
  expect_error(poisson_deviance(truth, -fitted), "nonnegative")
})

test_that("truth comparison uses the fixed floor without dropping zero predictions", {
  truth <- c(0, 0.25, 2)
  sim <- list(X = Matrix::Diagonal(3), y = c(0, 1, 2),
              y_true = truth, beta_true = truth)
  row <- score_fit("single", "zero", 0, rep(0, 3), sim)
  floor <- 1e-4
  expected <- 2 * (3 * floor - sum(truth) + sum(truth[-1] * log(truth[-1] / floor)))
  expect_equal(cfg$deviance_floor, floor)
  expect_equal(row$deviance_true, expected)
  expect_equal(poisson()$dev.resids(c(0, 1), c(0, 0), 1), c(0, Inf))
})

test_that("entry and row objectives use the fitter's lambda and observed likelihood", {
  X <- Matrix::Diagonal(3)
  beta <- rbind(c(1, 2), c(0, 3), c(0, 0))
  y <- rbind(c(1, 3), c(0, 4), c(0, 0))
  sim <- list(X = X, y = y, y_true = beta, beta_true = beta)
  row <- score_fit("multi", "test", 1, beta, sim)
  lambda <- log(3) * log(log(3))
  nll <- -2 * sum(dpois(as.vector(y), pmax(as.vector(beta), cfg$mu_floor), log = TRUE))
  expect_equal(row$lambda_L0, lambda)
  expect_equal(row$neg2loglik, nll)
  expect_equal(row$GIC, nll + 3 * lambda)
  expect_equal(row$GIC_rows, nll + 2 * lambda)
  expect_equal(row$nonzero_coefficients, 3)
  expect_equal(row$selected, 2)
  refs <- score_reference("multi", sim)
  expect_equal(refs$GIC_beta_true, row$GIC)
  expect_equal(refs$GIC_rows_beta_true, row$GIC_rows)
  expect_lte(refs$objective_lower_bound, min(row$GIC, row$GIC_rows))
  expect_equal(refs$deviance_true_lower_bound, 0)
})

test_that("comparator tuning selects the reported observed Poisson GIC", {
  skip_if_not_installed("glmnet")
  set.seed(91)
  X <- Matrix::Matrix(matrix(rexp(160), 40, 4), sparse = TRUE)
  y <- rpois(40, as.vector(X %*% c(5,0,2,0)))
  w <- 1 / (y + 0.1)
  sim <- list(X = X, y = y, y_true = X %*% c(5,0,2,0), beta_true = c(5,0,2,0))
  fit <- glmnet_single_gic(X,y,w,15L)
  row <- score_fit("single", "test", 1, fit$beta, sim)
  expect_equal(row$GIC, min(fit$GIC_path), tolerance = 1e-8)
  expect_true(all(fit$beta >= 0))
  expect_match(method_heading("L0glm", c(L0glm=1.23)), "Time = 1.23 s", fixed = TRUE)
})
