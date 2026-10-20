# High-dimensional Gaussian variable selection: full L0glm, SIS/ISIS,
# ABESS and relaxed LASSO. Use the current installed package, including SIS/ISIS.
library(L0glm)

# Editable settings. For a quick run, reduce p to e.g. 20000L.
# The original million-predictor design alone occupies about 4 GB as dense doubles.
n <- 500L
p <- 1000000L
n_nonzero <- 60L
beta_range <- c(-0.8, 0.8)
error_sd <- 1
seed_problem <- 123L
seed_noise <- 1L
nsis <- min(p, 2L * floor(n / log(n)))
nsis_isis <- min(p, max(1L, round(n / (log(p) * log(log(n))))))
isis_rounds <- 10L
nlambda <- 10L
abess_max_support <- min(p, n - 2L, 150L)

# Use one BLAS/OpenMP thread for comparable elapsed timings across methods.
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)

# Preserve the historical simulation, including positive and negative effects.
stopifnot(n > 2L, p >= n_nonzero, n_nonzero > 0L)
set.seed(seed_problem)
X <- matrix(rnorm(n * p), n, p)
true_support <- sample.int(p, n_nonzero)
beta_true <- numeric(p)
beta_true[true_support] <- runif(n_nonzero, beta_range[1L], beta_range[2L])
y_true <- as.numeric(X[, true_support, drop = FALSE] %*% beta_true[true_support])
set.seed(seed_noise)
y <- y_true + rnorm(n, sd = error_sd)

# Keep X dense: a fully populated sparse copy wastes memory. All methods omit
# the intercept and estimate Gaussian dispersion. GIC always uses the original p:
# GIC = -2*logLik + log(p)*log(log(n))*df, including one dispersion parameter.
lambda_gic <- log(p) * log(log(n))
fits_high_dimensional <- screening <- list()
seconds <- numeric()

message("Full L0glm")
seconds["Full L0glm"] <- system.time(
  fits_high_dimensional[["Full L0glm"]] <- L0glm.fit(X, y,
    family = gaussian(), lambda = lambda_gic, intercept = FALSE,
    nonnegative = FALSE, start = NULL, normalize = "auto", relaxed = TRUE,
    solver = "eigen", method = 1L, parallelize = FALSE, nthreads = 1L,
    screen = "none")
)[["elapsed"]]

# One call returns a full-dimension model for either screening strategy.
# ISIS starts afresh in each round and retains the best original-p GIC round.
# SIS uses a broad one-pass screen. ISIS uses the smaller historical batch size
# (20 here): admitting 160 per round gave much poorer recovery at p = 1 million.
# Neither screen guarantees signal retention.
for (name in c("SIS", "ISIS")) {
  label <- paste(name, "+ L0glm")
  message(label)
  seconds[label] <- system.time(
    fits_high_dimensional[[label]] <- L0glm.fit(X, y, family = gaussian(), screen = name,
      nsis = if (name == "SIS") nsis else nsis_isis,
      screen.maxit = isis_rounds, lambda = lambda_gic, intercept = FALSE,
      nonnegative = FALSE, start = NULL, normalize = "auto", relaxed = TRUE,
      solver = "eigen", method = 1L, parallelize = FALSE, nthreads = 1L)
  )[["elapsed"]]
  screening[[name]] <- fits_high_dimensional[[label]]$screening
}

# Competitors minimize the same GIC over their native paths. ABESS searches
# consecutive support sizes with warm starts, rather than its small default cap.
# Relaxed LASSO requests ten log-spaced native lambda values, refitting each
# support without shrinkage; the backend may stop a nearly saturated path early.
# The actual grid is retained in comparison$tuning.
# Timings include screening/fitting and, for competitors, path tuning/refitting.
message("ABESS and relaxed LASSO")
invisible(gc()) # Release large temporary design copies before native allocation.
comparison <- sparseglm(X, y, family = gaussian(), beta_true = beta_true,
  method = c("abess", "relaxed_lasso"), tune = "gic", nlambda = nlambda,
  intercept = FALSE, nonnegative = FALSE,
  support.size = 0:abess_max_support,
  options = list(abess = list(warm.start = TRUE)))

# Score every returned model on the original data with the package's common
# likelihood and support metrics. Screening and tuning never use beta_true.
columns <- c("method", "status", "seconds", "gic", "F1_score", "MCC",
             "support_size_estimated", "TP", "FP", "FN", "note")
scores <- lapply(names(fits_high_dimensional), function(name) {
  row <- cbind(data.frame(method = name, status = "ok", seconds = seconds[[name]],
                         note = ""),
               L0glm.metrics(fits_high_dimensional[[name]], beta_true = beta_true))
  row[columns]
})
competitor_scores <- L0glm.metrics(comparison)
competitor_scores$method <- c(abess = "ABESS", relaxed_lasso = "Relaxed LASSO")[
  competitor_scores$method]
benchmarks_high_dimensional <- rbind(do.call(rbind, scores), competitor_scores[columns])
rownames(benchmarks_high_dimensional) <- NULL
names(benchmarks_high_dimensional)[match(c("gic", "F1_score", "support_size_estimated"),
  names(benchmarks_high_dimensional))] <- c("GIC", "F1", "selected")

print(benchmarks_high_dimensional[setdiff(names(benchmarks_high_dimensional), "note")],
      row.names = FALSE, digits = 5)
if (any(nzchar(benchmarks_high_dimensional$note)))
  print(subset(benchmarks_high_dimensional, nzchar(note), select = c(method, note)),
        row.names = FALSE)

# Returned models retain the original dimensions; use coef(), predict(), GIC()
# or L0glm.metrics(). ISIS diagnostics: screening$ISIS$history / $best_iteration.
fits_high_dimensional <- c(fits_high_dimensional,
  setNames(comparison$fits, c("ABESS", "Relaxed LASSO")))
# Optional export (choose a path appropriate to the current working directory):
# write.csv(benchmarks_high_dimensional, "high_dimensional_metrics.csv", row.names = FALSE)
