# Full screening ranks and elapsed times on the same Gaussian design.
library(L0glm)
stopifnot(requireNamespace("SIS", quietly = TRUE))
RhpcBLASctl::blas_set_num_threads(1L)
n <- 500L
p <- 1000000L
nsis <- 100L
repetitions <- 3L

set.seed(123)
X <- matrix(rnorm(as.double(n) * p), n, p)
colnames(X) <- paste0("x", seq_len(p))
active <- sort(sample.int(p, 60L))
beta <- sample(c(-1, 1), 60L, TRUE) * runif(60L, .5, 1.5)
y <- as.numeric(X[, active, drop = FALSE] %*% beta) + rnorm(n)

results <- list()
reference_rank <- NULL
for (repetition in seq_len(repetitions)) {
  for (method in c("sis_glm_1", "sis_glm_auto", "SIS")) {
    gc()
    if (method == "sis_glm_1") seconds <- system.time(fit <-
      sis_glm(X, y, nsis = nsis, nthreads = 1L))[["elapsed"]]
    if (method == "sis_glm_auto") seconds <- system.time(fit <-
      sis_glm(X, y, nsis = nsis))[["elapsed"]]
    if (method == "SIS") seconds <- system.time(fit <-
      SIS::SIS(X, y, family = "gaussian", nsis = nsis, iter = FALSE,
        penalty = "lasso", tune = "bic", standardize = FALSE,
        parallel = FALSE, seed = 123L))[["elapsed"]]
    ranking <- as.integer(if (method == "SIS") fit$ix_list[[1L]] else fit$ranking)
    if (is.null(reference_rank)) reference_rank <- ranking
    row <- data.frame(method, repetition, seconds,
      threads = if (method == "SIS") NA_integer_ else fit$nthreads,
      identical_top_100 = identical(ranking[seq_len(nsis)], reference_rank[seq_len(nsis)]),
      identical_full_rank = identical(ranking, reference_rank),
      rank_correlation = cor(order(ranking), order(reference_rank)))
    stopifnot(row$identical_top_100, row$rank_correlation > 1 - 1e-10)
    results[[length(results) + 1L]] <- row
    print(row)
    rm(fit)
  }
}
results <- do.call(rbind, results)
aggregate(seconds ~ method, results, median)
root <- if (file.exists("L0glm/DESCRIPTION")) "." else ".."
output <- file.path(root, "benchmarks", "3_benchmarks_sis_vs_SIS_results")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
write.csv(results, file.path(output, "timings.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(output, "sessionInfo.txt"))

# SIS's initial Gaussian screen uses absolute correlations: standardization
# leaves these ranks unchanged. Its public call additionally fits a BIC-tuned
# LASSO on the 100 retained columns; compare ix_list[[1L]] before that selection.
