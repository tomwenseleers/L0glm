# QSAR: normalization and initialization of L0glm versus abess / relaxed LASSO.
# Run from the repository root, L0glm package directory, or this directory.
# Settings are ordinary R values; no environment-variable configuration.
cfg <- list(
  seed = 1L,
  train_fraction = 0.75,
  normalize = c("none", "l2", "rms", "sd", "auto"),
  start = c("NULL", "1", "original_ones", "pilot"),
  relaxed = FALSE,                 # use c(FALSE, TRUE) to compare relaxed refits too
  repeats = 3L,                    # median elapsed seconds and their range
  nlambda = 10L,
  lambda_min = 2,                  # same lower bound for every L0 path
  nfolds = 5L,
  gamma = c(0, 0.5, 1),            # fully / partly relaxed, ordinary LASSO
  maxit = 25L,
  compare_maxit = 200L,            # also test longer runs for l2/rms scalar-one starts
  historical_reference = TRUE,     # reproduce the original all-row L2 transform separately
  tol = 1e-10,
  thresh = 1e-5,
  threads = 1L,
  load_checkout = TRUE,            # FALSE to benchmark installed L0glm
  plot = TRUE,
  output = "results_normalization_start"
)

repo_dir <- if (file.exists("L0glm/DESCRIPTION")) "." else
  if (file.exists("../L0glm/DESCRIPTION")) ".." else "../.."
repo_dir <- normalizePath(repo_dir, winslash = "/", mustWork = TRUE)
stopifnot(file.exists(file.path(repo_dir, "L0glm", "DESCRIPTION")))
data_dir <- file.path(repo_dir, "benchmarks",
                      "2_benchmark_vs_abess_relaxed_LASSO_toxicityQSAR")
out_dir <- file.path(data_dir, cfg$output)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
required <- c("Matrix", "caret", "glmnet", "abess", "RhpcBLASctl",
              if (cfg$load_checkout) "pkgload" else "L0glm")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install these packages first: ", paste(missing, collapse = ", "))
if (cfg$load_checkout) {
  pkgload::load_all(file.path(repo_dir, "L0glm"), compile = FALSE, quiet = TRUE)
} else library(L0glm)
RhpcBLASctl::blas_set_num_threads(cfg$threads)
RhpcBLASctl::omp_set_num_threads(cfg$threads)
stopifnot(cfg$repeats >= 1L, cfg$nlambda >= 3L, cfg$nfolds >= 3L,
          all(cfg$normalize %in% c("none", "l2", "rms", "sd", "auto")),
          all(cfg$start %in% c("NULL", "1", "original_ones", "pilot")))

# Preserve the original seed-1, outcome-stratified 242/80 split. All transforms
# and near-invariance filtering use training observations only. Test outcomes
# are used only for final scoring, never for tuning or pilot construction.
molecules <- read.csv(file.path(data_dir, "aquatictox_molecules_withsmiles.csv"))
X_raw <- as.matrix(read.csv(file.path(data_dir, "aquatictox_covariates.csv")))
y <- molecules$toxicity
stopifnot(nrow(X_raw) == length(y), all(is.finite(X_raw)), all(is.finite(y)))
storage.mode(X_raw) <- "double"
set.seed(cfg$seed)
train <- as.vector(caret::createDataPartition(y, p = cfg$train_fraction, list = FALSE))
test <- setdiff(seq_along(y), train)
remove <- caret::nearZeroVar(X_raw[train, , drop = FALSE], freqCut = 99, uniqueCut = 1)
keep <- setdiff(seq_len(ncol(X_raw)), remove)
X_train <- X_raw[train, keep, drop = FALSE]
X_test <- X_raw[test, keep, drop = FALSE]
y_train <- y[train]
y_test <- y[test]
n <- nrow(X_train)
p <- ncol(X_train)
set.seed(cfg$seed + 1L)
foldid <- sample(rep(seq_len(cfg$nfolds), length.out = n))
lambda_gic <- log(p) * log(log(n))
cat(sprintf("QSAR: %d training / %d test observations, %d predictors; GIC lambda %.6f\n",
            n, nrow(X_test), p, lambda_gic))

# Common Gaussian GIC, profiling sigma^2 = RSS/n, for ALL methods:
# -2 log L + log(p)*log(log(n)) * (selected slopes + intercept + sigma^2).
# Count actual nonzero slopes and keep the original p after support restriction.
# This is an IC score, not a ridge-effective-df correction or an L0 certificate.
training_scores <- function(beta, X, y) {
  beta <- as.matrix(beta)          # first row is always the free intercept
  stopifnot(nrow(beta) == ncol(X) + 1L, all(is.finite(beta)))
  predicted <- sweep(X %*% beta[-1L, , drop = FALSE], 2L, beta[1L, ], "+")
  rss <- colSums((predicted - y)^2)
  k <- colSums(beta[-1L, , drop = FALSE] != 0)
  min2LL <- length(y) * (log(2 * pi) + 1 + log(rss / length(y)))
  data.frame(selected = k, train_RMSE = sqrt(rss / length(y)),
             dispersion_ML = rss / length(y), logLik = -min2LL / 2,
             GIC = min2LL + log(ncol(X)) * log(log(length(y))) * (k + 2L))
}

# Include complete training/selection work; exclude test scoring and file I/O.
# Record failures/warnings explicitly instead of dropping unsuccessful variants.
timed <- function(fun) {
  warnings <- character()
  value <- NULL
  error <- NA_character_
  seconds <- system.time(tryCatch(
    value <- withCallingHandlers(fun(), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) error <<- conditionMessage(e)))[["elapsed"]]
  list(value = value, seconds = seconds, error = error,
       warnings = paste(unique(warnings), collapse = " | "))
}

# Exact interpolants have unbounded profiled Gaussian likelihood. Flag them
# instead of silently accepting -Inf as a meaningful model-selection score.
best_gic <- function(scores) {
  available <- which(is.finite(scores$GIC))
  if (length(available) < nrow(scores))
    warning("Excluded candidates with nonfinite Gaussian GIC (including exact interpolants).")
  if (!length(available)) stop("No candidate has finite Gaussian GIC.")
  available[which.min(scores$GIC[available])]
}

rows <- models <- paths <- list()
record <- function(run, id, method, selection, normalization = NA_character_,
                   start = NA_character_, relaxed = NA, repetition = 1L,
                   overhead = 0, maxit = NA_integer_) {
  value <- run$value
  metadata <- data.frame(id, method, selection, normalization, start, relaxed, maxit,
                         repetition, seconds = run$seconds + overhead,
                         pilot_seconds = overhead, lambda = NA_real_, gamma = NA_real_,
                         selected = NA_integer_, train_RMSE = NA_real_,
                         dispersion_ML = NA_real_, logLik = NA_real_, GIC = NA_real_,
                         heldout_RMSE = NA_real_, heldout_R2 = NA_real_,
                         converged = NA, error = run$error, warnings = run$warnings)
  if (!is.null(value)) {
    beta <- as.numeric(value$beta)
    score <- training_scores(beta, X_train, y_train)
    prediction <- as.numeric(beta[1L] + X_test %*% beta[-1L])
    metadata[names(score)] <- score
    metadata$heldout_RMSE <- sqrt(mean((prediction - y_test)^2))
    metadata$heldout_R2 <- 1 - sum((prediction - y_test)^2) / sum((y_test - mean(y_test))^2)
    metadata$lambda <- value$lambda
    metadata$gamma <- if (is.null(value$gamma)) NA_real_ else value$gamma
    metadata$converged <- if (is.null(value$converged)) NA else value$converged
    if (!is.null(value$pilot_seconds)) metadata$pilot_seconds <- value$pilot_seconds
    if (repetition == 1L) models[[id]] <<- value
  }
  rows[[length(rows) + 1L]] <<- metadata
  invisible(metadata)
}

fit_value <- function(fit, lambda) {
  list(beta = fit$coefficients, lambda = lambda, fit = fit,
       converged = all(fit$convergence$converged))
}

# Internal normalization: none = no scaling; l2 = centered L2 norm 1;
# rms = centered L2 norm sqrt(n); sd = centered sample SD 1. Here auto = rms.
# Scalar start=1 is in NORMALIZED working coordinates; original_ones passes an
# explicit matrix of ones in original descriptor units. NULL initializes each fit
# automatically. "pilot" uses L0glm.tune's signed-mean-absolute common seed from
# ten independent fits; zero entries restrict support. Save the seed for replay.
for (norm in cfg$normalize) for (seed in cfg$start) for (relax in cfg$relaxed)
  for (budget in unique(c(cfg$maxit,
      if (seed == "1" && norm %in% c("l2", "rms")) cfg$compare_maxit))) {
  tag <- paste("L0glm", norm, seed, if (relax) "relaxed" else "unrelaxed",
               paste0("iter", budget), sep = "_")
  cat("Fitting ", tag, "\n", sep = "")
  common <- list(X = X_train, y = y_train, family = gaussian(identity),
                 intercept = TRUE, nonnegative = FALSE, normalize = norm,
                 relaxed = relax, parallelize = FALSE, openmp = FALSE,
                 solver = "eigen", method = 1L, maxit = budget,
                 tol = cfg$tol, thresh = cfg$thresh,
                 start = if (seed == "NULL") NULL else if (seed == "1") 1 else
                   if (seed == "original_ones") matrix(1, p + 1L, 1L) else "pilot")
  for (repetition in seq_len(cfg$repeats)) {
    # Includes endpoint search, all pilots, lambda fits and GIC selection.
    tuned <- timed(function() do.call(L0glm.tune, c(common,
      list(tune = "gic", nlambda = cfg$nlambda, lambda_min = cfg$lambda_min,
           prop_range = 1, plot = FALSE, printbest = FALSE))))
    if (!is.null(tuned$value)) {
      result <- tuned$value
      tuned$value <- fit_value(result$fit_best, result$lambda_best)
      tuned$value$pilot_seconds <- sum(result$pilot$seconds)
      if (repetition == 1L) paths[[tag]] <- result
    }
    record(tuned, paste0(tag, "_GIC"), "L0glm", "GIC path", norm, seed, relax, repetition,
           maxit = budget)

    # Charge pilot construction to the preset fit too, but not endpoint search.
    if (seed == "pilot" && is.null(tuned$value)) {
      preset <- tuned
      pilot_seconds <- 0
    } else {
      preset_args <- common
      if (seed == "pilot") preset_args$start <- result$start
      pilot_seconds <- if (seed == "pilot") sum(result$pilot$seconds) else 0
      preset <- timed(function() {
        fit <- do.call(L0glm.fit, c(preset_args, list(lambda = "gic")))
        fit_value(fit, lambda_gic)
      })
    }
    record(preset, paste0(tag, "_preset"), "L0glm", "GIC preset", norm, seed, relax,
           repetition, overhead = pilot_seconds, maxit = budget)
  }
  write.csv(do.call(rbind, rows), file.path(out_dir, "runs.csv"), row.names = FALSE)
}

# Historical reference only: the original script normalized/centered all input
# rows before splitting. Retain that exact convention here for reproduction;
# it is separate from every training-only comparison above. No test outcomes
# enter these fits. Convert coefficients back to raw descriptor units to score.
if (cfg$historical_reference) {
  raw_scale <- sqrt(colSums(X_raw[, keep, drop = FALSE]^2))
  unit <- sweep(X_raw[, keep, drop = FALSE], 2, raw_scale, "/")
  unit_center <- colMeans(unit)
  unit <- sweep(unit, 2, unit_center, "-")
  unit_scale <- sqrt(colSums(unit^2))
  X_reference <- sweep(unit, 2, unit_scale, "/")[train, , drop = FALSE]
  for (budget in unique(c(cfg$maxit, cfg$compare_maxit))) for (repetition in seq_len(cfg$repeats)) {
    run <- timed(function() {
      fit <- L0glm.fit(X_reference, y_train, family = gaussian(), intercept = TRUE,
        start = matrix(1, p + 1L, 1L), normalize = "none", lambda = "gic",
        nonnegative = FALSE, relaxed = FALSE, maxit = budget, tol = cfg$tol,
        parallelize = FALSE, openmp = FALSE, solver = "eigen", method = 1L)
      value <- fit_value(fit, lambda_gic)
      b <- as.numeric(value$beta)
      b[-1L] <- b[-1L] / (raw_scale * unit_scale)
      b[1L] <- b[1L] - sum(raw_scale * unit_center * b[-1L])
      value$beta <- Matrix::Matrix(b, ncol = 1L, sparse = TRUE)
      value$reference_transform <- list(center = raw_scale * unit_center,
                                         scale = raw_scale * unit_scale)
      value
    })
    record(run, paste0("historical_l2_iter", budget), "L0glm", "historical preset",
      "historical_l2", "1", FALSE, repetition, maxit = budget)
  }
}

# abess: training-only normalization, signed slopes and free intercept.
# Keep the original broad support search with independent (non-warm) fits.
# Rescore its returned path using the same GIC as L0glm.
for (repetition in seq_len(cfg$repeats)) {
  cat("Fitting abess, repetition ", repetition, "\n", sep = "")
  run <- timed(function() {
    path <- abess::abess(X_train, y_train, family = "gaussian", normalize = 1L,
      tune.type = "gic", fit.intercept = TRUE, support.size = 0:(n - 1L),
      warm.start = FALSE, num.threads = cfg$threads)
    beta <- coef(path, support.size = path$support.size)
    scores <- training_scores(beta, X_train, y_train)
    best <- best_gic(scores)
    list(beta = beta[, best, drop = FALSE], lambda = NA_real_,
         support_size = path$support.size[best], scores = scores, path = path)
  })
  record(run, "abess_GIC", "abess", "GIC support path", repetition = repetition)
}

# Relaxed LASSO: separate training-GIC and training-only CV selection.
# gamma=0 is the unpenalized refit; gamma=1 is ordinary LASSO.
for (repetition in seq_len(cfg$repeats)) {
  cat("Fitting relaxed LASSO, repetition ", repetition, "\n", sep = "")
  run <- timed(function() {
    path <- glmnet::glmnet(X_train, y_train, family = "gaussian", alpha = 1,
      standardize = TRUE, intercept = TRUE, relax = TRUE)
    candidates <- lapply(cfg$gamma, function(gamma) {
      beta <- coef(path, s = path$lambda, gamma = gamma)
      scores <- training_scores(beta, X_train, y_train)
      list(beta = beta, scores = transform(scores, lambda = path$lambda, gamma = gamma))
    })
    scores <- do.call(rbind, lapply(candidates, `[[`, "scores"))
    beta <- do.call(cbind, lapply(candidates, `[[`, "beta"))
    best <- best_gic(scores)
    list(beta = beta[, best, drop = FALSE], lambda = scores$lambda[best],
         gamma = scores$gamma[best], scores = scores, path = path)
  })
  record(run, "relaxed_LASSO_GIC", "relaxed LASSO", "GIC lambda/gamma path",
         repetition = repetition)
  cv_run <- timed(function() {
    path <- glmnet::cv.glmnet(X_train, y_train, family = "gaussian", alpha = 1,
      standardize = TRUE, intercept = TRUE, relax = TRUE, gamma = cfg$gamma,
      foldid = foldid, nfolds = cfg$nfolds, type.measure = "mse", parallel = FALSE)
    lambda <- path$relaxed$lambda.min
    gamma <- path$relaxed$gamma.min
    list(beta = coef(path, s = "lambda.min", gamma = gamma),
         lambda = lambda, gamma = gamma, path = path)
  })
  record(cv_run, "relaxed_LASSO_CV", "relaxed LASSO", paste0(cfg$nfolds, "-fold CV"),
         repetition = repetition)
}

runs <- do.call(rbind, rows)
summary <- do.call(rbind, lapply(split(runs, factor(runs$id, levels = unique(runs$id))), function(z) {
  first <- z[1L, setdiff(names(z), "repetition"), drop = FALSE]
  first$seconds <- median(z$seconds)
  first$seconds_min <- min(z$seconds)
  first$seconds_max <- max(z$seconds)
  first$pilot_seconds <- median(z$pilot_seconds)
  first$completed_runs <- sum(is.na(z$error))
  first$clean_runs <- sum(is.na(z$error) & !nzchar(z$warnings) & !z$converged %in% FALSE)
  first$repeats <- nrow(z)
  first$GIC_range <- if (any(is.finite(z$GIC))) diff(range(z$GIC, finite = TRUE)) else NA_real_
  first$RMSE_range <- if (any(is.finite(z$heldout_RMSE)))
    diff(range(z$heldout_RMSE, finite = TRUE)) else NA_real_
  first$warnings <- paste(unique(z$warnings[nzchar(z$warnings)]), collapse = " | ")
  first$error <- paste(unique(z$error[!is.na(z$error)]), collapse = " | ")
  first
}))
rownames(summary) <- NULL
write.csv(runs, file.path(out_dir, "runs.csv"), row.names = FALSE)
write.csv(summary, file.path(out_dir, "comparison.csv"), row.names = FALSE)
write.csv(data.frame(row = seq_along(y), partition = ifelse(seq_along(y) %in% train, "train", "test"),
                     cv_fold = replace(rep(NA_integer_, length(y)), train, foldid)),
          file.path(out_dir, "split.csv"), row.names = FALSE)
write.csv(data.frame(column = keep, predictor = colnames(X_train)),
          file.path(out_dir, "predictors.csv"), row.names = FALSE)
saveRDS(list(config = cfg, models = models, paths = paths, summary = summary,
             runs = runs, train = train, test = test, foldid = foldid, keep = keep,
             lambda_gic = lambda_gic, session = sessionInfo()),
        file.path(out_dir, "models.rds"))
writeLines(trimws(capture.output(sessionInfo()), which = "right"),
           file.path(out_dir, "sessionInfo.txt"))

# Independently verify the common GIC and original-scale training predictions.
for (value in models) {
  beta <- value$beta
  met <- getFromNamespace("get.metrics", "L0glm")(
    beta, cbind(Intercept = 1, X_train), matrix(y_train),
    family = gaussian(identity), p = p, always.include = 1L)
  stopifnot(isTRUE(all.equal(met$gic, training_scores(beta, X_train, y_train)$GIC,
                            tolerance = 1e-8)))
  if (!is.null(value$fit)) {
    fit <- value$fit
    stopifnot(inherits(fit$coefficients, "dgCMatrix"), nrow(fit$coefficients) == p + 1L,
              fit$p == p,
              isTRUE(all.equal(as.numeric(fit$X %*% fit$coefficients),
                as.numeric(beta[1L] + X_train %*% beta[-1L]), tolerance = 1e-8)))
  }
}

if (cfg$plot) {
  draw_comparison <- function(zoom = TRUE) {
    colors <- c("L0glm" = "#238443", "abess" = "#2166AC", "relaxed LASSO" = "#CB181D")
    labels <- ifelse(summary$method == "L0glm",
      paste(summary$normalization, paste0("start=", summary$start), summary$selection,
            paste0("it=", summary$maxit),
            ifelse(summary$relaxed, "relaxed", "")),
      paste(summary$method, summary$selection))
    flagged <- nzchar(summary$warnings) | nzchar(summary$error) | summary$converged %in% FALSE
    labels <- paste0(labels, ifelse(flagged, " *", ""))
    at <- rev(seq_len(nrow(summary)))
    old <- par(mfrow = c(1, 3), mar = c(4.5, 1, 3, 1), oma = c(1.2, 25, 2, 0))
    on.exit(par(old))
    for (field in c("GIC", "heldout_RMSE", "seconds")) {
      x <- summary[[field]]
      finite <- is.finite(x) & (field != "seconds" | x > 0)
      shown <- x
      upper <- if (zoom && field != "seconds")
        min(max(x[finite]), median(x[finite]) + 3 * mad(x[finite])) else Inf
      clipped <- finite & x > upper
      shown[clipped] <- upper
      log_axis <- field == "seconds" || (!zoom && field == "heldout_RMSE")
      limits <- if (log_axis) range(shown[finite]) * c(0.8, 1.2) else
        extendrange(shown[finite], f = 0.08)
      plot(shown[finite], at[finite], xlim = limits,
        ylim = c(0.4, nrow(summary) + 0.6), pch = 19,
        col = colors[summary$method[finite]], yaxt = "n",
        log = if (log_axis) "x" else "", ylab = "",
        xlab = switch(field, GIC = "Training GIC", heldout_RMSE = "Held-out RMSE",
                       seconds = "Elapsed seconds (median; log axis)"),
        main = switch(field, GIC = "Lower GIC is better", heldout_RMSE = "Lower RMSE is better",
                       seconds = "Training and selection time"))
      abline(h = at, col = "grey92", lty = 3)
      if (field == "GIC") axis(2, at = at, labels = labels, las = 1, tick = FALSE, cex.axis = 0.8)
      if (field == "seconds") segments(pmax(summary$seconds_min, 1e-4), at,
        pmax(summary$seconds_max, 1e-4), at, col = colors[summary$method])
      points(shown[finite], at[finite], pch = ifelse(clipped[finite], 17, 19),
             col = colors[summary$method[finite]])
      if (any(clipped)) text(shown[clipped], at[clipped],
        labels = formatC(x[clipped], digits = 3, format = "g"), pos = 2, cex = 0.65)
    }
    mtext(sprintf("QSAR: n=%d, p=%d; %d held-out molecules; %d timing repetitions",
                  n, p, nrow(X_test), cfg$repeats), outer = TRUE, side = 3, line = 0.5)
    mtext(paste("* Warnings or incomplete convergence.",
      if (zoom) "Triangles: larger values beyond zoomed range (actual values labeled)." else
        "Full metric ranges; RMSE and time use log axes."),
      outer = TRUE, side = 1, line = 0.1, cex = 0.75)
  }
  png(file.path(out_dir, "comparison.png"), width = 2400,
      height = max(1300, 150 + nrow(summary) * 40), res = 150)
  draw_comparison()
  dev.off()
  png(file.path(out_dir, "comparison_all.png"), width = 2400,
      height = max(1300, 150 + nrow(summary) * 40), res = 150)
  draw_comparison(zoom = FALSE)
  dev.off()
  pdf(file.path(out_dir, "comparison.pdf"), width = 16, height = max(9, nrow(summary) * 0.27))
  draw_comparison()
  dev.off()
  full_summary <- summary
  focus <- summary$method != "L0glm" | summary$normalization == "historical_l2" |
    (summary$normalization %in% c("l2", "rms") & summary$start == "1") |
    (summary$normalization == "rms" & summary$start == "NULL" & summary$selection == "GIC preset")
  summary <- summary[which(focus), , drop = FALSE]
  png(file.path(out_dir, "comparison_unit_starts.png"), width = 2400, height = 1100, res = 150)
  draw_comparison()
  dev.off()
  pdf(file.path(out_dir, "comparison_unit_starts.pdf"), width = 16, height = 8)
  draw_comparison()
  dev.off()
  summary <- full_summary
}

print(summary[, c("method", "selection", "normalization", "start", "relaxed",
                  "maxit", "selected", "GIC", "heldout_RMSE", "seconds", "completed_runs", "clean_runs")], row.names = FALSE)
cat("\nResults: ", normalizePath(out_dir, winslash = "/"), "\n", sep = "")
if (any(runs$converged %in% FALSE) || any(!is.na(runs$error)) || any(nzchar(runs$warnings)))
  message("Some fits reported warnings, failures or incomplete convergence; inspect runs.csv.")
