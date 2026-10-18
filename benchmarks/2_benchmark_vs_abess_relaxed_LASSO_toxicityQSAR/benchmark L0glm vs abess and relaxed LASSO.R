# QSAR: practical L0glm, abess and relaxed LASSO comparison.
# Run from the repository root, package directory, or this script's directory.
root <- if (file.exists("L0glm/DESCRIPTION")) "." else
  if (file.exists("DESCRIPTION")) ".." else "../.."
data_dir <- file.path(root, "benchmarks", "2_benchmark_vs_abess_relaxed_LASSO_toxicityQSAR")
out_dir <- file.path(data_dir, "results_default_options")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pkgload::load_all(file.path(root, "L0glm"), compile = FALSE, quiet = TRUE)
# Use library(L0glm) instead of load_all() to benchmark the installed package.
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)
repetitions <- 3L

molecules <- read.csv(file.path(data_dir, "aquatictox_molecules_withsmiles.csv"))
X <- as.matrix(read.csv(file.path(data_dir, "aquatictox_covariates.csv")))
storage.mode(X) <- "double"
y <- molecules$toxicity

# Keep the original stratified split; filter descriptors using training rows only.
set.seed(1)
train <- as.vector(caret::createDataPartition(y, p = 0.75, list = FALSE))
test <- setdiff(seq_along(y), train)
remove <- caret::nearZeroVar(X[train, , drop = FALSE], freqCut = 99, uniqueCut = 1)
keep <- setdiff(seq_len(ncol(X)), remove)
X_train <- X[train, keep, drop = FALSE]
X_test <- X[test, keep, drop = FALSE]
y_train <- y[train]
y_test <- y[test]
n <- nrow(X_train)
p <- ncol(X_train)
lambda_gic <- log(p) * log(log(n))

# Common Gaussian GIC: profile sigma^2 as RSS/n; count slopes + intercept + sigma^2.
gic <- function(beta) {
  beta <- as.matrix(beta)  # coefficients always include the intercept first
  residuals <- sweep(cbind(1, X_train) %*% beta, 1L, y_train, "-")
  rss <- colSums(residuals^2)
  selected <- colSums(beta[-1L, , drop = FALSE] != 0)
  n * (log(2*pi) + 1 + log(rss/n)) + lambda_gic * (selected + 2L)
}
choose_gic <- function(beta, fit) {
  scores <- gic(beta)
  available <- which(is.finite(scores))
  if (!length(available)) stop("No candidate has finite Gaussian GIC.")
  best <- available[which.min(scores[available])]
  list(beta = beta[, best, drop = FALSE], fit = fit)
}
benchmark <- function(method, fit) {
  seconds <- numeric(repetitions)
  for (i in seq_len(repetitions)) {
    gc()
    seconds[i] <- system.time(model <- fit())[["elapsed"]]
  }
  beta <- as.numeric(model$beta)
  prediction <- as.numeric(beta[1L] + X_test %*% beta[-1L])
  metrics <- data.frame(method, selected = sum(beta[-1L] != 0),
    GIC = gic(model$beta), heldout_RMSE = sqrt(mean((prediction-y_test)^2)),
    seconds = median(seconds), seconds_min = min(seconds), seconds_max = max(seconds))
  list(metrics = metrics, model = model)
}

# auto normalization and NULL initialization are the defaults. They work well
# here; no manual centering, scaling, pilot or starting coefficients are needed.
results <- list()
results$L0glm <- benchmark("L0glm (GIC preset)", function() {
  fit <- L0glm::L0glm.fit(X_train, y_train, family = gaussian(),
    lambda = "gic", intercept = TRUE, nonnegative = FALSE, nthreads = 1L)
  list(beta = fit$coefficients, fit = fit)
})
results$abess <- benchmark("abess (GIC)", function() {
  fit <- abess::abess(X_train, y_train, family = "gaussian", normalize = 1L,
    fit.intercept = TRUE, tune.type = "gic", num.threads = 1L)
  choose_gic(coef(fit, support.size = fit$support.size), fit)
})
results$relaxed_LASSO <- benchmark("relaxed LASSO (GIC)", function() {
  fit <- glmnet::glmnet(X_train, y_train, family = "gaussian", alpha = 1,
    standardize = TRUE, intercept = TRUE, relax = TRUE)
  # gamma=0: refitted least squares; gamma=1: LASSO; 0.5: partial relaxation.
  beta <- do.call(cbind, lapply(c(0, 0.5, 1), function(gamma)
    coef(fit, s = fit$lambda, gamma = gamma)))
  choose_gic(beta, fit)
})

comparison <- do.call(rbind, lapply(results, `[[`, "metrics"))
rownames(comparison) <- NULL
print(comparison, row.names = FALSE)
write.csv(comparison, file.path(out_dir, "comparison.csv"), row.names = FALSE)
saveRDS(list(results = results, train = train, test = test, predictors = keep,
  lambda_gic = lambda_gic, repetitions = repetitions, session = sessionInfo()),
  file.path(out_dir, "models.rds"))

# Timings include fitting and training-GIC selection, excluding held-out scoring.
# L0glm uses its preset directly; the competitors search their regularization paths.
draw_comparison <- function() {
  old <- par(mfrow = c(1, 3), mar = c(7, 4.5, 3, 1), las = 2)
  on.exit(par(old))
  for (metric in c("GIC", "heldout_RMSE", "seconds")) {
    title <- switch(metric, GIC = "Training GIC", heldout_RMSE = "Held-out RMSE",
                    seconds = "Fit and selection time (s)")
    barplot(comparison[[metric]], names.arg = c("L0glm", "abess", "relaxed LASSO"),
      col = c("#238443", "#2166AC", "#CB181D"), border = NA,
      main = title, ylab = "", cex.names = 0.9)
  }
}
png(file.path(out_dir, "comparison.png"), width = 1500, height = 650, res = 130)
draw_comparison()
dev.off()
pdf(file.path(out_dir, "comparison.pdf"), width = 11.5, height = 5)
draw_comparison()
dev.off()
