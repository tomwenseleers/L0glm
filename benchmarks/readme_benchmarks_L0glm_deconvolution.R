# Reproduce the README benchmark examples for L0glm deconvolution.
#
# The examples use n = p = 10000 and npeaks = 1000 by default. Timings include
# regularisation selection for glmnet and abess. L0glm uses lambda = "gic",
# which applies the GIC preset directly inside L0glm.fit and therefore avoids a
# separate tuning pass.
#
# Outputs are written to benchmarks/readme_benchmarks/.

`%||%` <- function(x, y) if (is.null(x)) y else x

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file) && nzchar(script_file)) dirname(normalizePath(script_file)) else getwd()
root_dir <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
if (!file.exists(file.path(root_dir, "L0glm", "DESCRIPTION"))) {
  root_dir <- normalizePath(getwd(), mustWork = FALSE)
}

pkg_dir <- file.path(root_dir, "L0glm")
out_dir <- file.path(root_dir, "benchmarks", "readme_benchmarks")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (file.exists(file.path(pkg_dir, "DESCRIPTION")) &&
    requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(pkg_dir, quiet = TRUE)
} else {
  library(L0glm)
}

library(Matrix)

cfg <- list(
  n = as.integer(Sys.getenv("L0GLM_README_N", "10000")),
  p = as.integer(Sys.getenv("L0GLM_README_P", "10000")),
  npeaks = as.integer(Sys.getenv("L0GLM_README_NPEAKS", "1000")),
  m = as.integer(Sys.getenv("L0GLM_README_CHANNELS", "20")),
  seed_single = as.integer(Sys.getenv("L0GLM_README_SEED_SINGLE", "2")),
  seed_multi = as.integer(Sys.getenv("L0GLM_README_SEED_MULTI", "3")),
  nlambda = as.integer(Sys.getenv("L0GLM_README_NLAMBDA", "100")),
  nfolds = as.integer(Sys.getenv("L0GLM_README_NFOLDS", "5")),
  ncores = as.integer(Sys.getenv(
    "L0GLM_README_NCORES",
    as.character(max(1L, min(8L, parallel::detectCores(logical = FALSE) %||% 1L)))
  )),
  run_abess = tolower(Sys.getenv("L0GLM_README_RUN_ABESS", "true")) %in% c("1", "true", "yes"),
  run_glmnet = tolower(Sys.getenv("L0GLM_README_RUN_GLMNET", "true")) %in% c("1", "true", "yes")
)

support_grid <- sort(unique(as.integer(c(
  seq(100L, min(2500L, cfg$p), by = 100L),
  cfg$npeaks,
  round(cfg$npeaks * c(0.5, 0.75, 1.25, 1.5, 2))
))))
support_grid <- support_grid[support_grid >= 1L & support_grid < cfg$p]

gic_penalty <- function(n, p) log(p) * log(log(n))

time_it <- function(expr) {
  gc()
  t0 <- proc.time()[["elapsed"]]
  value <- force(expr)
  list(value = value, seconds = proc.time()[["elapsed"]] - t0)
}

safe_time <- function(label, expr) {
  message(sprintf("[%s] starting", label))
  ans <- tryCatch(time_it(expr), error = function(e) e)
  if (inherits(ans, "error")) {
    message(sprintf("[%s] failed: %s", label, conditionMessage(ans)))
    list(value = NULL, seconds = NA_real_, error = conditionMessage(ans))
  } else {
    message(sprintf("[%s] finished in %.2f s", label, ans$seconds))
    c(ans, list(error = NA_character_))
  }
}

as_matrix <- function(x) {
  if (inherits(x, "Matrix")) as.matrix(x) else if (is.null(dim(x))) matrix(x, ncol = 1L) else as.matrix(x)
}

as_sparse_matrix <- function(x) {
  if (inherits(x, "Matrix")) Matrix::Matrix(x, sparse = TRUE) else Matrix::Matrix(as_matrix(x), sparse = TRUE)
}

make_kernel_design <- function(n, p) {
  kernel <- convkernel("bigaussian", 51L, params = list(sigma1 = 1.25, sigma2 = 5), perc = 0)
  X <- convmatrix(kernel, n = n, n_cols = p)
  dimnames(X) <- list(paste0("t", seq_len(n)), paste0("x", seq_len(p)))
  X
}

simulate_single <- function(n, p, npeaks, seed) {
  sim <- simul(method = 2, n = n, p = p, s = npeaks / p, sparse = TRUE, seed = seed)
  sim$y <- as.vector(sim$y)
  sim$y_true <- as.vector(sim$y_true)
  dimnames(sim$X) <- list(paste0("t", seq_len(nrow(sim$X))), paste0("x", seq_len(ncol(sim$X))))
  names(sim$y) <- rownames(sim$X)
  names(sim$y_true) <- rownames(sim$X)
  rownames(sim$beta_true) <- colnames(sim$X)
  sim$support <- which(Matrix::rowSums(abs(sim$beta_true) > 1e-5) > 0)
  sim
}

simulate_multi <- function(n, p, npeaks, m, seed) {
  sim <- simulate_spike_train(
    n = n, p = p, k = npeaks,
    family = "poisson",
    multichannel = TRUE,
    sparse = TRUE,
    seed = seed,
    mean_beta = 10000
  )
  dimnames(sim$X) <- list(paste0("t", seq_len(nrow(sim$X))), paste0("x", seq_len(ncol(sim$X))))
  rownames(sim$beta_true) <- colnames(sim$X)

  if (!is.null(m) && is.finite(m) && m > 0L && m < ncol(sim$y)) {
    channel_score <- Matrix::colSums(sim$beta_true)
    keep <- order(channel_score, decreasing = TRUE)[seq_len(m)]
    sim$y <- sim$y[, keep, drop = FALSE]
    sim$y_true <- sim$y_true[, keep, drop = FALSE]
    sim$beta_true <- sim$beta_true[, keep, drop = FALSE]
  }
  colnames(sim$y) <- colnames(sim$y_true) <- colnames(sim$beta_true) <- paste0("ch", seq_len(ncol(sim$y)))
  rownames(sim$y) <- rownames(sim$y_true) <- rownames(sim$X)
  sim$support <- which(Matrix::rowSums(abs(sim$beta_true) > 1e-5) > 0)
  sim
}

score_fit <- function(case, method, seconds, beta_hat, sim, note = NA_character_, tol = 1e-5) {
  beta_hat <- as_sparse_matrix(beta_hat)
  beta_true <- as_sparse_matrix(sim$beta_true)
  X <- sim$X
  y <- as_matrix(sim$y)
  y_true <- as_matrix(sim$y_true)
  m <- ncol(y)

  selected <- Matrix::rowSums(abs(beta_hat) > tol) > 0
  truth <- Matrix::rowSums(abs(beta_true) > tol) > 0
  tp <- sum(selected & truth)
  fp <- sum(selected & !truth)
  fn <- sum(!selected & truth)
  tn <- sum(!selected & !truth)

  nz <- which(Matrix::rowSums(abs(beta_hat)) > tol)
  mu_hat <- matrix(0, nrow(X), m)
  if (length(nz)) {
    mu_hat <- as.matrix(X[, nz, drop = FALSE] %*% beta_hat[nz, , drop = FALSE])
  }
  mu_hat <- pmax(mu_hat, 1e-8)

  fam <- poisson(identity)
  dev_obs <- sum(fam$dev.resids(y = as.vector(y), mu = as.vector(mu_hat), wt = 1))
  dev_true <- sum(fam$dev.resids(y = as.vector(y_true), mu = as.vector(mu_hat), wt = 1))
  min2ll <- fam$aic(y = round(as.vector(y)), n = rep(1, length(y)),
                    mu = as.vector(mu_hat), wt = rep(1, length(y)),
                    dev = dev_obs)
  edf <- sum(abs(beta_hat) > tol)
  gic <- min2ll + gic_penalty(nrow(X) * m, ncol(X)) * edf

  r2 <- suppressWarnings(summary(lm(as.vector(beta_true) ~ 0 + as.vector(beta_hat)))$r.squared)
  precision <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  recall <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_

  data.frame(
    case = case,
    method = method,
    seconds = seconds,
    n = nrow(X),
    p = ncol(X),
    m = m,
    npeaks = sum(truth),
    selected = sum(selected),
    TP = tp,
    FP = fp,
    FN = fn,
    TN = tn,
    precision = precision,
    recall = recall,
    F1 = if ((2 * tp + fp + fn) > 0) 2 * tp / (2 * tp + fp + fn) else NA_real_,
    R2_coefficients = r2,
    deviance = dev_obs,
    deviance_true = dev_true,
    GIC = gic,
    note = note,
    check.names = FALSE
  )
}

glmnet_single_gic <- function(X, y, weights, nlambda) {
  fit <- glmnet::glmnet(
    x = X, y = y, weights = weights, family = "gaussian",
    alpha = 1, nlambda = nlambda, intercept = FALSE,
    standardize = FALSE, lower.limits = 0
  )
  pred <- X %*% fit$beta
  rss <- Matrix::colSums(weights * (matrix(y, nrow = length(y), ncol = ncol(pred)) - pred)^2)
  df <- Matrix::colSums(fit$beta != 0)
  rdf <- pmax(length(y) - df, 1)
  ll <- -0.5 * length(y) * (log(2 * pi) - log(rdf) + log(rss)) - 0.5 * rdf
  gic <- -2 * ll + gic_penalty(length(y), ncol(X)) * df
  j <- which.min(gic)
  list(beta = fit$beta[, j, drop = FALSE], lambda = fit$lambda[j], GIC_path = gic)
}

glmnet_multi_gic <- function(X, Y, row_weights, nlambda) {
  fit <- glmnet::glmnet(
    x = X, y = Y, weights = row_weights, family = "mgaussian",
    alpha = 1, nlambda = nlambda, intercept = FALSE,
    standardize = FALSE, standardize.response = FALSE, lower.limits = 0
  )
  nlambda_fit <- length(fit$lambda)
  rss <- df <- numeric(nlambda_fit)
  for (j in seq_len(nlambda_fit)) {
    beta_j <- do.call(cbind, lapply(fit$beta, function(B) B[, j]))
    pred <- as.matrix(X %*% beta_j)
    rss[j] <- sum(row_weights * rowSums((Y - pred)^2))
    df[j] <- sum(abs(beta_j) > 0)
  }
  nobs <- length(Y)
  rdf <- pmax(nobs - df, 1)
  ll <- -0.5 * nobs * (log(2 * pi) - log(rdf) + log(rss)) - 0.5 * rdf
  gic <- -2 * ll + gic_penalty(nobs, ncol(X)) * df
  j <- which.min(gic)
  beta <- do.call(cbind, lapply(fit$beta, function(B) B[, j]))
  list(beta = beta, lambda = fit$lambda[j], GIC_path = gic)
}

abess_single_gic <- function(X, y, weights, support_grid) {
  fit <- abess::abess(
    x = X, y = y, weight = weights, family = "gaussian",
    tune.path = "sequence", support.size = support_grid,
    fit.intercept = FALSE, beta.low = 0, lambda = 0,
    warm.start = TRUE, tune.type = "bic",
    ic.scale = gic_penalty(length(y), ncol(X)) / log(length(y))
  )
  best_size <- fit[["best.size"]][1]
  beta <- coef(fit, support.size = best_size, sparse = FALSE)[-1, , drop = FALSE]
  list(beta = beta, support.size = best_size)
}

abess_multi_gic <- function(X, Y, support_grid) {
  fit <- abess::abess(
    x = X, y = Y, family = "mgaussian",
    tune.path = "sequence", support.size = support_grid,
    fit.intercept = FALSE, beta.low = 0, lambda = 0,
    warm.start = TRUE, tune.type = "bic",
    ic.scale = gic_penalty(length(Y), ncol(X)) / log(length(Y))
  )
  best_size <- fit[["best.size"]][1]
  beta_raw <- coef(fit, support.size = best_size, sparse = FALSE)
  beta <- if (is.list(beta_raw)) beta_raw[[1]][-1, , drop = FALSE] else beta_raw[-1, , drop = FALSE]
  list(beta = beta, support.size = best_size)
}

plot_single_signal <- function(sim, fits, file) {
  png(file, width = 1400, height = 820, res = 150)
  on.exit(dev.off(), add = TRUE)
  idx <- seq_len(min(1200L, length(sim$y)))
  ymax <- as.numeric(stats::quantile(c(sim$y[idx], sim$y_true[idx]), 0.995, na.rm = TRUE)) * 1.15
  par(mar = c(4.2, 4.6, 3.5, 1.2), bty = "l")
  plot(idx, sim$y[idx], type = "l", col = "grey45", lwd = 1,
       xlab = "Time index", ylab = "Signal", ylim = c(0, ymax * 1.08),
       main = "Single-channel Poisson deconvolution (zoom; extreme fits capped)")
  lines(idx, pmin(sim$y_true[idx], ymax), col = "#1B9E77", lwd = 2)
  cols <- c(L0glm = "#0072B2", glmnet_lasso = "#D55E00", abess = "#CC79A7")
  for (nm in intersect(names(fits), names(cols))) {
    beta <- as_sparse_matrix(fits[[nm]])
    mu <- as.vector(sim$X %*% beta)
    lines(idx, pmin(mu[idx], ymax), col = cols[[nm]], lwd = 1.6)
  }
  true_beta <- as.vector(sim$beta_true)
  spike_idx <- which(true_beta > 0 & seq_along(true_beta) %in% idx)
  segments(spike_idx, 0, spike_idx, pmin(true_beta[spike_idx], ymax), col = "#E69F00")
  legend("topright", bty = "n", lwd = c(1, 2, 1.6, 1.6, 1.6, 1),
         col = c("grey45", "#1B9E77", cols), legend = c("observed", "no-noise truth", names(cols)))
}

plot_coefficients <- function(case, sim, fits, file) {
  png(file, width = 1400, height = 920, res = 150)
  on.exit(dev.off(), add = TRUE)
  methods <- names(fits)
  par(mfrow = c(1, length(methods)), mar = c(4.3, 4.3, 4.0, 1.0), oma = c(0, 0, 1, 0), bty = "l")
  truth <- as_matrix(sim$beta_true)
  xtruth <- if (ncol(truth) == 1L) truth[, 1] else rowMeans(truth)
  cap <- as.numeric(stats::quantile(xtruth[xtruth > 0], 0.995, na.rm = TRUE))
  for (nm in methods) {
    est <- as_matrix(fits[[nm]])
    xest <- if (ncol(est) == 1L) est[, 1] else rowMeans(est)
    keep <- xtruth != 0 | abs(xest) > 1e-5
    r2 <- suppressWarnings(summary(lm(xtruth ~ 0 + xest))$r.squared)
    plot(pmin(xest[keep], cap), pmin(xtruth[keep], cap), pch = 16, cex = 0.35, col = "#0072B280",
         xlim = c(0, cap), ylim = c(0, cap),
         xlab = "Estimated coefficient", ylab = "True coefficient",
         main = sprintf("%s\nR2 = %.3f", nm, r2))
    abline(0, 1, col = "#D55E00", lwd = 2)
  }
  mtext(sprintf("%s coefficient recovery (axes capped at 99th pct.)", case), outer = TRUE, cex = 1.1)
}

plot_multi_heatmap <- function(sim, fits, file) {
  png(file, width = 1400, height = 980, res = 150)
  on.exit(dev.off(), add = TRUE)
  methods <- c("truth", names(fits))
  panels <- length(methods)
  par(mfrow = c(panels, 1), mar = c(2.4, 4.2, 2.6, 1.0), oma = c(2, 0, 1, 0))
  cols <- colorRampPalette(c("white", "#56B4E9", "#0072B2", "#D55E00"))(256)
  truth_support <- which(Matrix::rowSums(abs(sim$beta_true) > 1e-5) > 0)
  for (nm in methods) {
    mat <- if (nm == "truth") as_matrix(sim$beta_true) else as_matrix(fits[[nm]])
    z <- log1p(mat)
    image(seq_len(nrow(mat)), seq_len(ncol(mat)), z,
          col = cols, useRaster = TRUE, axes = FALSE,
          xlab = "", ylab = "Channel", main = nm)
    axis(2, las = 1)
    if (nm != "truth") abline(v = truth_support, col = "#E69F0040", lwd = 0.4)
  }
  axis(1)
  mtext("Time index", side = 1, outer = TRUE)
  mtext("Multichannel grouped deconvolution coefficient maps", outer = TRUE, cex = 1.1)
}

plot_summary <- function(results, file) {
  png(file, width = 1400, height = 850, res = 150)
  on.exit(dev.off(), add = TRUE)
  ok <- results[!is.na(results$seconds), , drop = FALSE]
  par(mfrow = c(1, 3), mar = c(8, 4.5, 3, 1), bty = "l")
  cols <- c("#0072B2", "#D55E00", "#CC79A7", "#009E73", "#E69F00", "#999999")
  for (metric in c("seconds", "F1", "log10_deviance_true")) {
    vals <- if (metric == "log10_deviance_true") log10(ok$deviance_true) else ok[[metric]]
    names(vals) <- paste(ok$case, ok$method, sep = "\n")
    bp <- barplot(vals, las = 2, col = cols[seq_along(vals)], ylab = metric,
                  main = metric, cex.names = 0.72)
    if (metric == "seconds") text(bp, vals, labels = sprintf("%.1fs", vals), pos = 3, cex = 0.7)
  }
}

message("Simulating single-channel data")
single <- simulate_single(cfg$n, cfg$p, cfg$npeaks, cfg$seed_single)
message("Simulating multichannel data")
multi <- simulate_multi(cfg$n, cfg$p, cfg$npeaks, cfg$m, cfg$seed_multi)

results <- list()
single_fits <- list()
multi_fits <- list()

l0_single <- safe_time("single L0glm", L0glm.fit(
  single$X, single$y,
  family = poisson(identity), lambda = "gic",
  nonnegative = TRUE, relaxed = FALSE,
  intercept = FALSE, solver = "eigen", method = 7L,
  parallelize = FALSE, nthreads = 1L
))
if (!is.null(l0_single$value)) {
  single_fits$L0glm <- l0_single$value$coefficients
  results[[length(results) + 1L]] <- score_fit("single", "L0glm", l0_single$seconds,
                                               l0_single$value$coefficients, single)
} else {
  results[[length(results) + 1L]] <- data.frame(case = "single", method = "L0glm",
                                                seconds = NA_real_, note = l0_single$error)
}

if (cfg$run_glmnet && requireNamespace("glmnet", quietly = TRUE)) {
  w <- 1 / (single$y + 0.1)
  gl_single <- safe_time("single glmnet LASSO", glmnet_single_gic(single$X, single$y, w, cfg$nlambda))
  if (!is.null(gl_single$value)) {
    single_fits$glmnet_lasso <- gl_single$value$beta
    results[[length(results) + 1L]] <- score_fit("single", "glmnet_lasso", gl_single$seconds,
                                                 gl_single$value$beta, single,
                                                 note = sprintf("GIC-selected lambda %.4g", gl_single$value$lambda))
  } else {
    results[[length(results) + 1L]] <- data.frame(case = "single", method = "glmnet_lasso",
                                                  seconds = NA_real_, note = gl_single$error)
  }
}

if (cfg$run_abess && requireNamespace("abess", quietly = TRUE)) {
  w <- 1 / (single$y + 0.1)
  ab_single <- safe_time("single abess", abess_single_gic(single$X, single$y, w, support_grid))
  if (!is.null(ab_single$value)) {
    single_fits$abess <- ab_single$value$beta
    results[[length(results) + 1L]] <- score_fit("single", "abess", ab_single$seconds,
                                                 ab_single$value$beta, single,
                                                 note = sprintf("GIC-scaled BIC over %d support sizes; selected %d",
                                                                length(support_grid), ab_single$value$support.size))
  } else {
    results[[length(results) + 1L]] <- data.frame(case = "single", method = "abess",
                                                  seconds = NA_real_, note = ab_single$error)
  }
}

l0_multi <- safe_time("multichannel grouped L0glm", L0glm.fit(
  multi$X, multi$y,
  family = poisson(identity), lambda = "gic",
  nonnegative = TRUE, grouped = TRUE, relaxed = TRUE,
  intercept = FALSE, solver = "eigen", method = 7L,
  parallelize = cfg$ncores > 1L, ncores = cfg$ncores, nthreads = 1L,
  group.maxit = 6L
))
if (!is.null(l0_multi$value)) {
  multi_fits$L0glm_grouped <- l0_multi$value$coefficients
  results[[length(results) + 1L]] <- score_fit("multichannel", "L0glm_grouped", l0_multi$seconds,
                                               l0_multi$value$coefficients, multi,
                                               note = sprintf("parallelize=%s, ncores=%d",
                                                              cfg$ncores > 1L, cfg$ncores))
} else {
  results[[length(results) + 1L]] <- data.frame(case = "multichannel", method = "L0glm_grouped",
                                                seconds = NA_real_, note = l0_multi$error)
}

if (cfg$run_glmnet && requireNamespace("glmnet", quietly = TRUE)) {
  row_w <- 1 / (rowMeans(multi$y) + 0.1)
  gl_multi <- safe_time("multichannel glmnet group LASSO", glmnet_multi_gic(multi$X, multi$y, row_w, cfg$nlambda))
  if (!is.null(gl_multi$value)) {
    multi_fits$glmnet_group_lasso <- gl_multi$value$beta
    results[[length(results) + 1L]] <- score_fit("multichannel", "glmnet_group_lasso", gl_multi$seconds,
                                                 gl_multi$value$beta, multi,
                                                 note = sprintf("GIC-selected lambda %.4g", gl_multi$value$lambda))
  } else {
    results[[length(results) + 1L]] <- data.frame(case = "multichannel", method = "glmnet_group_lasso",
                                                  seconds = NA_real_, note = gl_multi$error)
  }
}

if (cfg$run_abess && requireNamespace("abess", quietly = TRUE)) {
  ab_multi <- safe_time("multichannel abess BSGS", abess_multi_gic(multi$X, multi$y, support_grid))
  if (!is.null(ab_multi$value)) {
    multi_fits$abess_BSGS <- ab_multi$value$beta
    results[[length(results) + 1L]] <- score_fit("multichannel", "abess_BSGS", ab_multi$seconds,
                                                 ab_multi$value$beta, multi,
                                                 note = sprintf("GIC-scaled BIC over %d support sizes; selected %d",
                                                                length(support_grid), ab_multi$value$support.size))
  } else {
    results[[length(results) + 1L]] <- data.frame(case = "multichannel", method = "abess_BSGS",
                                                  seconds = NA_real_, note = ab_multi$error)
  }
}

common_cols <- c("case", "method", "seconds", "n", "p", "m", "npeaks", "selected",
                 "TP", "FP", "FN", "TN", "precision", "recall", "F1",
                 "R2_coefficients", "deviance", "deviance_true", "GIC", "note")
results_df <- do.call(rbind, lapply(results, function(x) {
  missing <- setdiff(common_cols, names(x))
  for (nm in missing) x[[nm]] <- NA
  x[common_cols]
}))

results_df <- results_df[order(results_df$case, is.na(results_df$seconds), results_df$deviance_true), ]
write.csv(results_df, file.path(out_dir, "readme_benchmark_results.csv"), row.names = FALSE)
saveRDS(list(config = cfg, support_grid = support_grid, results = results_df,
             single = single, multi = multi,
             single_fits = single_fits, multi_fits = multi_fits),
        file.path(out_dir, "readme_benchmark_results.rds"))

if (length(single_fits)) {
  plot_single_signal(single, single_fits, file.path(out_dir, "single_channel_signal_zoom.png"))
  plot_coefficients("Single-channel", single, single_fits, file.path(out_dir, "single_channel_coefficients.png"))
}
if (length(multi_fits)) {
  plot_coefficients("Multichannel", multi, multi_fits, file.path(out_dir, "multichannel_coefficients.png"))
  plot_multi_heatmap(multi, multi_fits, file.path(out_dir, "multichannel_coefficient_maps.png"))
}
plot_summary(results_df, file.path(out_dir, "benchmark_summary.png"))

print(results_df)
message("Wrote outputs to: ", normalizePath(out_dir, winslash = "/"))
