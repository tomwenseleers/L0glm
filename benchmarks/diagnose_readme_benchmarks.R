# Audit the saved README benchmark without rerunning the expensive fits.
# Run from the repository root or the L0glm/ package directory.
root_dir <- if (file.exists("L0glm/DESCRIPTION")) "." else ".."
root_dir <- normalizePath(root_dir, winslash = "/", mustWork = TRUE)
old_run <- getOption("L0glm.readme.run")
options(L0glm.readme.run = FALSE)
source(file.path(root_dir, "benchmarks", "readme_benchmarks_L0glm_deconvolution.R"))
options(L0glm.readme.run = old_run)

# Edit these paths directly to audit a different saved run.
result_file <- file.path(root_dir, "benchmarks", "readme_benchmarks", "readme_benchmark_results.rds")
diagnostic_dir <- file.path(root_dir, "benchmarks", "readme_benchmarks", "diagnostics")
dir.create(diagnostic_dir, recursive = TRUE, showWarnings = FALSE)
saved <- readRDS(result_file)
deviance_floor <- saved$config$deviance_floor
stopifnot(length(deviance_floor) == 1L, is.finite(deviance_floor), deviance_floor > 0)

manual_deviance_terms <- function(truth, fitted, floor) {
  mu <- pmax(as.vector(fitted), floor)
  truth <- as.vector(truth)
  value <- mu - truth
  positive <- truth > 0
  value[positive] <- value[positive] + truth[positive] * log(truth[positive] / mu[positive])
  2 * value
}

checks <- breakdown <- list()
for (case in c("single", "multichannel")) {
  sim <- if (case == "single") saved$single else saved$multi
  fits <- if (case == "single") saved$single_fits else saved$multi_fits
  truth <- as_matrix(sim$y_true)
  stopifnot(max(abs(sim$X %*% sim$beta_true - truth)) < 1e-8)
  bins <- cut(as.vector(truth), c(-Inf, 0, 1, 10, 100, Inf),
              labels = c("zero", "0 to 1", "1 to 10", "10 to 100", "over 100"))
  for (method in names(fits)) {
    beta <- as_matrix(fits[[method]])
    mu <- as.matrix(sim$X %*% beta)
    terms <- manual_deviance_terms(truth, mu, deviance_floor)
    row <- saved$results[saved$results$case == case & saved$results$method == method, ]
    stopifnot(abs(sum(terms) - row$deviance_true) < 1e-5)
    if (startsWith(method, "glmnet")) {
      path <- saved$tuning_paths[[if (case == "single") "glmnet_single" else "glmnet_multi"]]
      stopifnot(isTRUE(all.equal(row$GIC, min(path$GIC_path), tolerance = 1e-8)))
    }
    checks[[length(checks) + 1L]] <- data.frame(
      case, method, reported_deviance_true = row$deviance_true,
      independently_calculated_deviance_true = sum(terms),
      coefficient_relative_squared_error = sum((beta - sim$beta_true)^2) / sum(sim$beta_true^2),
      signal_relative_squared_error = sum((mu - truth)^2) / sum(truth^2),
      zero_predictions_with_positive_truth = sum(mu == 0 & truth > 0),
      floored_predictions_with_positive_truth = sum(mu <= deviance_floor & truth > 0))
    totals <- tapply(terms, bins, sum)
    breakdown[[length(breakdown) + 1L]] <- data.frame(
      case, method, true_mean_bin = names(totals), deviance_true = as.numeric(totals))
  }
}
checks <- do.call(rbind, checks)
write.csv(checks, file.path(diagnostic_dir, "independent_deviance_check.csv"), row.names = FALSE)
write.csv(do.call(rbind, breakdown), file.path(diagnostic_dir, "deviance_by_true_mean.csv"), row.names = FALSE)
print(checks)

# Minimal counterexample: OLS = (1, -0.3); clipping gives (1, 0), RSS = 3.
# Correct NNLS = (0.25, 0), RSS = 0.75. With both predictors available,
# beta.low=0 should solve the constrained problem rather than clip OLS.
X <- cbind(constant = rep(1, 4), increasing = 1:4)
y <- c(1, 0, 0, 0)
bound_checks <- list()
if (requireNamespace("abess", quietly = TRUE)) {
  for (sparse in c(FALSE, TRUE)) for (multichannel in c(FALSE, TRUE)) {
    design <- if (sparse) Matrix::Matrix(X, sparse = TRUE) else X
    response <- if (multichannel) cbind(y, 2*y) else y
    fit <- abess::abess(design, response,
                       family = if (multichannel) "mgaussian" else "gaussian",
                       support.size = 2, fit.intercept = FALSE, normalize = 0,
                       beta.low = 0, lambda = 0, num.threads = 1)
    beta <- if (multichannel) as.matrix(fit$beta[[1]]) else as.matrix(fit$beta)
    bound_checks[[length(bound_checks) + 1L]] <- data.frame(
      abess_version = as.character(packageVersion("abess")), sparse, multichannel,
      first_coefficient = beta[1, 1], second_coefficient = beta[2, 1],
      residual_sum_squares = sum((response - X %*% beta)^2),
      exact_NNLS_residual_sum_squares = if (multichannel) 3.75 else 0.75)
  }
  bound_checks <- do.call(rbind, bound_checks)
  write.csv(bound_checks, file.path(diagnostic_dir, "abess_bound_reproducer.csv"), row.names = FALSE)
  print(bound_checks)
}

if (requireNamespace("glmnet", quietly = TRUE)) {
  weight_message <- tryCatch({
    glmnet::glmnet(X, cbind(y, y), family = "mgaussian", weights = matrix(1, 4, 2))
    "Matrix weights accepted by this installed version."
  }, error = function(e) conditionMessage(e))
  writeLines(c(paste("glmnet version", packageVersion("glmnet")), weight_message),
             file.path(diagnostic_dir, "glmnet_matrix_weights.txt"))
  cat("glmnet matrix-weight check:", weight_message, "\n")
}
writeLines(capture.output(sessionInfo()), file.path(diagnostic_dir, "sessionInfo.txt"))
