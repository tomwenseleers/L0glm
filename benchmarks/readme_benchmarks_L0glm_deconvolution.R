# Reproduce the README benchmark examples for L0glm deconvolution.
#
# The examples use n = p = 10000, 200 m/z mass channels and npeaks = 1000 by default. Timings include
# regularisation selection for glmnet. All paths
# are scored with observed-data Poisson GIC. L0glm uses lambda = "gic",
# which applies the GIC preset directly inside L0glm.fit and therefore avoids a
# separate tuning pass. Eligible fits use the compiled L0nnpois backend through
# L0glm.fit(openmp = TRUE). Set cfg$openmp to FALSE to compare the R backend.
#
# Local outputs are written to codex_sandbox/readme_benchmarks/ (Git-ignored).
# Run interactively with source(), or execute sections after editing cfg below.
# In RStudio, with the working directory set to the package directory L0glm/:
# source("../benchmarks/readme_benchmarks_L0glm_deconvolution.R")
# From the repository root, use source("benchmarks/readme_benchmarks_L0glm_deconvolution.R").
# For a smaller run, edit cfg below: n = p = 600L, npeaks = 40L,
# m = 4L, nlambda = 30L, and choose a separate cfg$out_dir.

# Paths are relative to the working directory: package directory or repository root.
root_dir <- if (file.exists("L0glm/DESCRIPTION")) "." else ".."
root_dir <- normalizePath(root_dir, winslash = "/", mustWork = TRUE)
pkg_dir <- file.path(root_dir, "L0glm")
if (!file.exists(file.path(pkg_dir, "DESCRIPTION"))) {
  stop("Set the working directory to the L0glm package directory or repository root.")
}

if (file.exists(file.path(pkg_dir, "DESCRIPTION")) &&
    requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(pkg_dir, quiet = TRUE, export_all = FALSE)
} else {
  library(L0glm)
}

library(Matrix)

readme_flag <- function(x, name) {
  if (is.logical(x) && length(x) == 1L && !is.na(x)) return(x)
  if (is.character(x) && length(x) == 1L && !is.na(x)) {
    value <- tolower(trimws(x))
    if (value %in% c("true", "false")) return(value == "true")
  }
  stop(name, " must be TRUE or FALSE (or the string 'true' or 'false').")
}

cfg <- list(
  n = 10000L,
  p = 10000L,
  npeaks = 1000L,
  m = 200L,
  seed_single = 2L,
  seed_multi = 3L,
  nlambda = 100L,
  ncores = 8L,
  openmp = TRUE,
  thresh = 1e-5, # relative selection cutoff; absolute floor is 1e-20
  run_glmnet = TRUE,
  mu_floor = 1e-8, # numerical floor for observed likelihood/GIC and path selection
  deviance_floor = 1e-4, # prediction floor for comparison against noise-free truth
  plot_from = 1L,
  plot_n = 500L,
  out_dir = file.path(root_dir, "codex_sandbox", "readme_benchmarks")
)
for (flag in c("run_glmnet", "openmp")) cfg[[flag]] <- readme_flag(cfg[[flag]], flag)
for (field in c("n", "p", "npeaks", "m", "nlambda", "ncores", "plot_from", "plot_n")) {
  value <- cfg[[field]]
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value < 1 || value != floor(value)) stop("cfg$", field, " must be a positive integer.")
}
if (cfg$n < 3L || cfg$p < 2L || cfg$npeaks > min(cfg$n, cfg$p)) {
  stop("Use n >= 3, p >= 2 and npeaks <= min(n, p).")
}
stopifnot(length(cfg$mu_floor) == 1L, is.finite(cfg$mu_floor), cfg$mu_floor > 0)
stopifnot(length(cfg$deviance_floor) == 1L, is.finite(cfg$deviance_floor), cfg$deviance_floor > 0)
out_dir <- cfg$out_dir

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
  sim <- simulate_spike_train(n = n, p = p, k = npeaks, family = "poisson",
                              multichannel = FALSE, mean_beta = 2000,
                              sparse = TRUE, seed = seed)
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

support_statistics <- function(beta_hat, beta_true, tol = 1e-5) {
  selected <- Matrix::rowSums(abs(as_sparse_matrix(beta_hat)) > tol) > 0
  truth <- Matrix::rowSums(abs(as_sparse_matrix(beta_true)) > tol) > 0
  if (length(selected) != length(truth)) stop("Coefficient row counts differ.")
  tp <- as.double(sum(selected & truth))
  fp <- as.double(sum(selected & !truth))
  fn <- as.double(sum(!selected & truth))
  tn <- as.double(sum(!selected & !truth))
  margins <- c(tp + fp, tp + fn, tn + fp, tn + fn)
  list(selected = selected, truth = truth, TP = tp, FP = fp, FN = fn, TN = tn,
       precision = if (tp + fp > 0) tp / (tp + fp) else NA_real_,
       recall = if (tp + fn > 0) tp / (tp + fn) else NA_real_,
       F1 = if (2 * tp + fp + fn > 0) 2 * tp / (2 * tp + fp + fn) else NA_real_,
       # MCC = Pearson's phi for binary support; undefined for constant support.
       MCC = if (all(margins > 0)) (tp * tn - fp * fn) / prod(sqrt(margins)) else NA_real_)
}

# Pool individual channel entries within the specified rows (micro-averaging).
# TP-row conditioning measures spectral recovery after successful peak detection;
# all-true-row conditioning additionally counts spectra from missed peaks as FNs.
spectral_support_statistics <- function(beta_hat, beta_true, rows, tol = 1e-5) {
  selected <- as.matrix(beta_hat)[rows, , drop = FALSE]
  truth <- as.matrix(beta_true)[rows, , drop = FALSE]
  if (!identical(dim(selected), dim(truth))) stop("Coefficient dimensions differ.")
  selected <- abs(selected) > tol
  truth <- abs(truth) > tol
  tp <- as.double(sum(selected & truth)); fp <- as.double(sum(selected & !truth))
  fn <- as.double(sum(!selected & truth)); tn <- as.double(sum(!selected & !truth))
  margins <- c(tp + fp, tp + fn, tn + fp, tn + fn)
  list(TP = tp, FP = fp, FN = fn, TN = tn,
       F1 = if (2*tp + fp + fn > 0) 2*tp/(2*tp + fp + fn) else NA_real_,
       MCC = if (all(margins > 0)) (tp*tn-fp*fn)/prod(sqrt(margins)) else NA_real_)
}

# R-squared from regression through the origin, as used by get.metrics().
# Compute from all coefficient entries, including all channels (not row means).
coefficient_r2 <- function(beta_hat, beta_true) {
  est <- as_sparse_matrix(beta_hat)
  truth <- as_sparse_matrix(beta_true)
  ss_est <- sum(est^2)
  ss_true <- sum(truth^2)
  if (ss_true == 0) return(NA_real_)
  if (ss_est == 0) return(0)
  min(1, (sum(est * truth) / sqrt(ss_est) / sqrt(ss_true))^2)
}

# The floor is applied to predictions only; the reference response is unchanged.
# For noise-free y, this is 2 * KL(Poisson(y) || Poisson(pmax(mu, mu_floor))).
poisson_deviance <- function(y, mu, mu_floor = cfg$mu_floor) {
  if (length(y) != length(mu)) stop("Response and prediction sizes differ.")
  if (any(!is.finite(y)) || any(y < 0) || any(!is.finite(mu)) || any(mu < 0)) {
    stop("Poisson responses and predictions must be finite and nonnegative.")
  }
  sum(poisson()$dev.resids(as.vector(y), pmax(as.vector(mu), mu_floor), 1))
}

# Counts and lambda use the same n (rows of X) as L0glm's GIC preset.
# The entry criterion is used to select comparator paths; row GIC is also
# reported, but the current grouped L0 fitter is a support-sharing heuristic,
# not an exact optimizer of the row-count criterion.
poisson_gic <- function(X, y, beta, tol = 1e-5) {
  beta <- as_sparse_matrix(beta)
  mu <- as.matrix(X %*% beta)
  deviance <- poisson_deviance(y, mu)
  nonzero <- sum(abs(beta) > tol)
  selected <- sum(Matrix::rowSums(abs(beta) > tol) > 0)
  lambda <- gic_penalty(nrow(X), ncol(X))
  list(deviance = deviance, nonzero_coefficients = nonzero,
       selected = selected, lambda_L0 = lambda,
       # These omit the saturated-likelihood constant, which is common to a path.
       criterion = deviance + lambda * nonzero,
       criterion_rows = deviance + lambda * selected, mu = mu)
}

score_fit <- function(case, method, seconds, beta_hat, sim, note = NA_character_, tol = 1e-5) {
  beta_hat <- as_sparse_matrix(beta_hat)
  beta_true <- as_sparse_matrix(sim$beta_true)
  X <- sim$X
  y <- as_matrix(sim$y)
  y_true <- as_matrix(sim$y_true)
  m <- ncol(y)

  support <- support_statistics(beta_hat, beta_true, tol)
  spectra <- spectral_support_statistics(beta_hat, beta_true,
                                          support$selected & support$truth, tol)
  spectra_all <- spectral_support_statistics(beta_hat, beta_true, support$truth, tol)
  # A one-channel TP row has no within-row sparsity question to evaluate.
  if (m == 1L) spectra[] <- spectra_all[] <- NA_real_

  scores <- poisson_gic(X, y, beta_hat, tol)
  dev_true <- poisson_deviance(y_true, scores$mu, mu_floor = cfg$deviance_floor)
  min2ll <- -2 * sum(dpois(as.vector(y), pmax(as.vector(scores$mu), cfg$mu_floor), log = TRUE))

  r2 <- coefficient_r2(beta_hat, beta_true)

  data.frame(
    case = case,
    method = method,
    seconds = seconds,
    n = nrow(X),
    p = ncol(X),
    m = m,
    npeaks = sum(support$truth),
    selected = sum(support$selected),
    TP = support$TP,
    FP = support$FP,
    FN = support$FN,
    TN = support$TN,
    precision = support$precision,
    recall = support$recall,
    F1 = support$F1,
    MCC = support$MCC,
    spectra_TP_rows_TP = spectra$TP,
    spectra_TP_rows_FP = spectra$FP,
    spectra_TP_rows_FN = spectra$FN,
    spectra_TP_rows_TN = spectra$TN,
    spectra_TP_rows_F1 = spectra$F1,
    spectra_TP_rows_MCC = spectra$MCC,
    spectra_all_true_rows_F1 = spectra_all$F1,
    spectra_all_true_rows_MCC = spectra_all$MCC,
    R2_coefficients = r2,
    deviance = scores$deviance,
    deviance_true = dev_true,
    neg2loglik = min2ll,
    nonzero_coefficients = scores$nonzero_coefficients,
    lambda_L0 = scores$lambda_L0,
    GIC = min2ll + scores$lambda_L0 * scores$nonzero_coefficients,
    GIC_rows = min2ll + scores$lambda_L0 * scores$selected,
    note = note,
    check.names = FALSE
  )
}

# beta_true is an oracle reference, not the minimum of a noisy penalized fit.
# The unpenalized saturated likelihood is a rigorous lower bound for BOTH
# objectives (nonnegative penalties); it need not be attainable by X %*% beta.
score_reference <- function(case, sim) {
  oracle <- score_fit(case, "beta_true", NA_real_, sim$beta_true, sim)
  y <- as.vector(sim$y)
  data.frame(case = case, lambda_L0 = oracle$lambda_L0,
             true_rows = oracle$selected,
             true_nonzero_coefficients = oracle$nonzero_coefficients,
             neg2loglik_beta_true = oracle$neg2loglik,
             GIC_beta_true = oracle$GIC,
             GIC_rows_beta_true = oracle$GIC_rows,
             objective_lower_bound = -2 * sum(dpois(y, y, log = TRUE)),
             deviance_true_beta_true = oracle$deviance_true,
             deviance_true_lower_bound = 0)
}

glmnet_single_gic <- function(X, y, weights, nlambda) {
  fit <- glmnet::glmnet(
    x = X, y = y, weights = weights, family = "gaussian",
    alpha = 1, nlambda = nlambda, intercept = FALSE,
    standardize = FALSE, lower.limits = 0
  )
  gic <- vapply(seq_along(fit$lambda), function(j) {
    poisson_gic(X, y, fit$beta[, j, drop = FALSE])$criterion
  }, numeric(1))
  gic <- gic - 2 * sum(dpois(y, y, log = TRUE))
  j <- which.min(gic)
  list(beta = fit$beta[, j, drop = FALSE], lambda = fit$lambda[j], GIC_path = gic, lambda_path = fit$lambda)
}

glmnet_multi_gic <- function(X, Y, row_weights, nlambda) {
  # glmnet's mgaussian family requires length-n weights, shared by all channels.
  # These are a one-step Gaussian approximation to Poisson precision weights.
  fit <- glmnet::glmnet(
    x = X, y = Y, weights = row_weights, family = "mgaussian",
    alpha = 1, nlambda = nlambda, intercept = FALSE,
    standardize = FALSE, standardize.response = FALSE, lower.limits = 0
  )
  nlambda_fit <- length(fit$lambda)
  gic <- numeric(nlambda_fit)
  for (j in seq_len(nlambda_fit)) {
    beta_j <- do.call(cbind, lapply(fit$beta, function(B) B[, j]))
    gic[j] <- poisson_gic(X, Y, beta_j)$criterion
  }
  gic <- gic - 2 * sum(dpois(as.vector(Y), as.vector(Y), log = TRUE))
  j <- which.min(gic)
  beta <- do.call(cbind, lapply(fit$beta, function(B) B[, j]))
  list(beta = beta, lambda = fit$lambda[j], GIC_path = gic, lambda_path = fit$lambda)
}

method_label <- function(method) {
  labels <- c(L0glm = "Nonnegative L0 Poisson",
              L0glm_grouped = "Nonnegative relaxed group L0 Poisson",
              glmnet_lasso = "Nonnegative LASSO (glmnet)",
              glmnet_group_lasso = "Nonnegative group LASSO (glmnet)")
  if (method %in% names(labels)) labels[[method]] else method
}

method_heading <- function(method, timings = NULL) {
  label <- method_label(method)
  if (method %in% names(timings) && is.finite(timings[[method]])) {
    label <- paste0(label, sprintf("\nTime = %.2f s", timings[[method]]))
  }
  label
}

plot_indices <- function(n, from = 1L, size = n) {
  from <- min(max(1L, as.integer(from)), n)
  seq.int(from, min(n, from + size - 1L))
}

plot_single_signal <- function(sim, fits, file, from = 1L, size = 500L, timings = NULL) {
  if (!length(fits)) return(invisible(NULL))
  png(file, width = 1800, height = 450 * length(fits) + 150, res = 150)
  on.exit(dev.off(), add = TRUE)
  idx <- plot_indices(length(sim$y), from, size)
  true_beta <- as.vector(sim$beta_true)
  spike_idx <- intersect(idx, seq_along(true_beta))
  # Peak areas and the convolved signal have different units. Scale BOTH sets
  # of coefficient stems by the same PSF maximum, and state that factor.
  spike_scale <- max(sim$X)
  par(mfrow = c(length(fits), 1L), mar = c(3.1, 6.2, 4.0, 1.2),
      oma = c(1.7, 0, 3.1, 0), bty = "l", las = 1)
  for (nm in names(fits)) {
    beta <- as_sparse_matrix(fits[[nm]])
    est <- as.vector(beta)
    mu <- as.vector(sim$X[idx, , drop = FALSE] %*% beta)
    yrange <- range(0, sim$y[idx], true_beta[spike_idx] * spike_scale,
                    -sim$y_true[idx], -mu, -est[spike_idx] * spike_scale, finite = TRUE)
    if (diff(yrange) == 0) yrange <- c(-1, 1)
    plot(idx, sim$y[idx], type = "l", col = "grey40", lwd = 1,
         xlab = "Time index", ylab = "", ylim = yrange * 1.08,
         main = method_heading(nm, timings))
    mtext("Signal (fit mirrored)", side = 2, line = 4.5, las = 0)
    abline(h = 0, col = "grey70", lty = 3)
    segments(spike_idx, 0, spike_idx, true_beta[spike_idx] * spike_scale, col = "red", lwd = 1)
    segments(spike_idx, 0, spike_idx, -est[spike_idx] * spike_scale, col = "darkblue", lwd = 1)
    lines(idx, -sim$y_true[idx], col = "grey60", lwd = 1.1, lty = 2)
    lines(idx, -mu, col = "#E69F00", lwd = 1.5)
  }
  mtext("Single-channel Poisson deconvolution", outer = TRUE, line = 1.8, font = 2)
  mtext(sprintf("Spikes are peak areas multiplied by the PSF maximum (%.4g); negative values indicate mirroring only.",
                spike_scale), side = 1, outer = TRUE, line = 0.3, cex = 0.8)
  # Draw a shared legend in figure coordinates, above the first panel.
  par(fig = c(0, 1, 0, 1), mar = rep(0, 4), oma = rep(0, 4), new = TRUE)
  plot.new()
  legend("top", inset = c(0, 0.035), xpd = NA, bty = "n", horiz = TRUE,
         legend = c("Observed (+)", "True spikes (+)", "Estimated spikes (-)",
                    "No-noise signal (-)", "Fitted signal (-)"),
         col = c("grey40", "red", "darkblue", "grey60", "#E69F00"),
         lty = c(1, 1, 1, 2, 1), lwd = c(1, 1, 1, 1, 1.5), cex = 0.8)
}

plot_coefficients <- function(case, sim, fits, file, tol = 1e-5, timings = NULL) {
  if (!length(fits)) return(invisible(NULL))
  png(file, width = 650 * length(fits), height = 750, res = 150)
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(1, length(fits)), mar = c(4.5, 4.5, 6.8, 1.0),
      oma = c(2, 0, 2, 0), bty = "l")
  truth <- as.vector(sim$beta_true)
  for (nm in names(fits)) {
    est <- as.vector(fits[[nm]])
    keep <- abs(truth) > tol | abs(est) > tol
    r2 <- coefficient_r2(fits[[nm]], sim$beta_true)
    stats <- support_statistics(fits[[nm]], sim$beta_true, tol)
    xmax <- max(1, est, na.rm = TRUE)
    ymax <- max(1, truth, na.rm = TRUE)
    plot(est[keep], truth[keep], pch = 16, cex = 0.85, col = "#4682B4B0",
         xlim = c(0, xmax * 1.03), ylim = c(0, ymax * 1.03),
         xlab = "Estimated coefficients", ylab = "True coefficients",
         main = sprintf("%s\nR2 = %.4f\nF1 = %.3f; MCC = %.3f",
                        method_heading(nm, timings), r2, stats$F1, stats$MCC), cex.main = 0.95)
    abline(0, 1, col = "grey50", lty = 2)
    mtext(sprintf("TP = %d; FP = %d; FN = %d", stats$TP, stats$FP, stats$FN),
          side = 3, line = 0.2, cex = 0.75)
  }
  mtext(paste(case, "coefficient recovery"), outer = TRUE, font = 2, line = 0.5)
  mtext("All coefficient entries; axes are not capped. R2: regression through zero. F1, MCC and counts: row support.",
        side = 1, outer = TRUE, line = 0.6, cex = 0.75)
}

plot_multi_support <- function(sim, fits, file, from = 1L, size = nrow(sim$beta_true), tol = 1e-5, timings = NULL) {
  methods <- c("truth", names(fits))
  png(file, width = 1800, height = 300 * length(methods) + 160, res = 150)
  on.exit(dev.off(), add = TRUE)
  idx <- plot_indices(nrow(sim$beta_true), from, size)
  truth <- Matrix::rowSums(abs(sim$beta_true) > tol) > 0
  par(mfrow = c(length(methods), 1), mar = c(2.8, 4.4, 4.3, 1.2),
      oma = c(2, 0, 2.7, 0), bty = "l")
  for (nm in methods) {
    title <- if (nm == "truth") "True support" else method_heading(nm, timings)
    if (nm != "truth") {
      stats <- support_statistics(fits[[nm]], sim$beta_true, tol)
      title <- sprintf("%s\nTP %d; FP %d; FN %d  |  F1 %.3f; MCC %.3f",
                       title, stats$TP, stats$FP, stats$FN, stats$F1, stats$MCC)
    }
    plot(range(idx), c(0, 1), type = "n", xlim = range(idx) + c(-0.5, 0.5),
         ylim = c(0, 1), xlab = "Time index", ylab = "Support", yaxt = "n",
         main = title, cex.main = 0.95)
    abline(h = 0, col = "grey85")
    if (nm == "truth") {
      at <- idx[truth[idx]]
      segments(at, 0, at, 1, col = "darkblue", lwd = 1.3)
    } else {
      colours <- ifelse(truth & stats$selected, "darkblue",
                        ifelse(stats$selected, "red", "magenta"))
      at <- idx[(truth | stats$selected)[idx]]
      segments(at, 0, at, 1, col = colours[at], lwd = 1.3)
    }
  }
  mtext("Multichannel support recovery", outer = TRUE, line = 1.4, font = 2)
  mtext("Dark blue: true support / true positives     Red: false positives     Magenta: false negatives",
        outer = TRUE, line = 0.1, cex = 0.85)
  mtext("Counts and scores refer to the full support, including indices outside the displayed window.",
        side = 1, outer = TRUE, line = 0.4, cex = 0.75)
}

plot_summary <- function(results, file) {
  png(file, width = 1800, height = 850, res = 150)
  on.exit(dev.off(), add = TRUE)
  ok <- results[is.finite(results$seconds), , drop = FALSE]
  if (!nrow(ok)) {
    plot.new()
    text(0.5, 0.5, "No successful fits; see the result table for errors.")
    return(invisible(NULL))
  }
  par(mfrow = c(1, 4), mar = c(9, 4.3, 3, 1), bty = "l")
  cols <- rep(c("#0072B2", "#D55E00", "#CC79A7", "#009E73", "#E69F00", "#999999"), length.out = nrow(ok))
  for (metric in c("seconds", "F1", "MCC", "log10_deviance_true")) {
    vals <- if (metric == "log10_deviance_true") log10(pmax(ok$deviance_true, .Machine$double.eps)) else ok[[metric]]
    names(vals) <- paste(ok$case, ok$method, sep = "\n")
    limits <- if (metric == "MCC") c(-1, 1) else if (metric == "F1") c(0, 1) else range(0, vals, finite = TRUE)
    if (diff(limits) == 0) limits <- c(0, 1)
    if (metric == "seconds") limits[2] <- limits[2] * 1.12
    bp <- barplot(replace(vals, !is.finite(vals), 0), las = 2, col = cols,
                  ylab = metric, main = metric, cex.names = 0.72, ylim = limits)
    if (any(!is.finite(vals))) text(bp[!is.finite(vals)], 0, "NA", pos = 3, cex = 0.7)
    if (metric == "seconds") text(bp, vals, labels = sprintf("%.1fs", vals), pos = 3, cex = 0.7)
  }
}

# Set options(L0glm.readme.run = FALSE) to load helpers without running fits.
if (isTRUE(getOption("L0glm.readme.run", TRUE))) {
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Simulating single-channel data")
single <- simulate_single(cfg$n, cfg$p, cfg$npeaks, cfg$seed_single)
message("Simulating multichannel data")
multi <- simulate_multi(cfg$n, cfg$p, cfg$npeaks, cfg$m, cfg$seed_multi)

results <- list()
single_fits <- list()
multi_fits <- list()
tuning_paths <- list()

l0_single <- safe_time("single L0glm", L0glm.fit(
  single$X, single$y,
  family = poisson(identity), lambda = "gic", thresh = cfg$thresh,
  nonnegative = TRUE, relaxed = FALSE,
  intercept = FALSE, solver = "eigen", method = 1L,
  openmp = cfg$openmp, parallelize = FALSE, nthreads = 1L
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
    tuning_paths$glmnet_single <- gl_single$value[c("GIC_path", "lambda_path")]
    single_fits$glmnet_lasso <- gl_single$value$beta
    results[[length(results) + 1L]] <- score_fit("single", "glmnet_lasso", gl_single$seconds,
                                                 gl_single$value$beta, single,
                                                 note = sprintf("Poisson GIC-selected lambda %.4g", gl_single$value$lambda))
  } else {
    results[[length(results) + 1L]] <- data.frame(case = "single", method = "glmnet_lasso",
                                                  seconds = NA_real_, note = gl_single$error)
  }
}

l0_multi <- safe_time("multichannel grouped L0glm",
                      L0glm.fit(
  multi$X, multi$y,
  family = poisson(identity), lambda = "gic", thresh = cfg$thresh,
  nonnegative = TRUE, grouped = TRUE, relaxed = TRUE,
  intercept = FALSE, solver = "eigen", method = 1L,
  openmp = cfg$openmp,
  parallelize = cfg$ncores > 1L, ncores = cfg$ncores, nthreads = 1L,
  group.maxit = 6L
))
if (!is.null(l0_multi$value)) {
  multi_fits$L0glm_grouped <- l0_multi$value$coefficients
  results[[length(results) + 1L]] <- score_fit("multichannel", "L0glm_grouped", l0_multi$seconds,
                                               l0_multi$value$coefficients, multi,
                                               note = sprintf("openmp=%s, parallelize=%s, ncores=%d",
                                                              cfg$openmp, cfg$ncores > 1L, cfg$ncores))
} else {
  results[[length(results) + 1L]] <- data.frame(case = "multichannel", method = "L0glm_grouped",
                                                seconds = NA_real_, note = l0_multi$error)
}

if (cfg$run_glmnet && requireNamespace("glmnet", quietly = TRUE)) {
  row_w <- 1 / (rowMeans(multi$y) + 0.1)
  gl_multi <- safe_time("multichannel glmnet group LASSO", glmnet_multi_gic(multi$X, multi$y, row_w, cfg$nlambda))
  if (!is.null(gl_multi$value)) {
    tuning_paths$glmnet_multi <- gl_multi$value[c("GIC_path", "lambda_path")]
    multi_fits$glmnet_group_lasso <- gl_multi$value$beta
    results[[length(results) + 1L]] <- score_fit("multichannel", "glmnet_group_lasso", gl_multi$seconds,
                                                 gl_multi$value$beta, multi,
                                                 note = sprintf("Poisson GIC-selected lambda %.4g; shared row weights", gl_multi$value$lambda))
  } else {
    results[[length(results) + 1L]] <- data.frame(case = "multichannel", method = "glmnet_group_lasso",
                                                  seconds = NA_real_, note = gl_multi$error)
  }
}

common_cols <- c("case", "method", "seconds", "n", "p", "m", "npeaks", "selected",
                 "TP", "FP", "FN", "TN", "precision", "recall", "F1", "MCC",
                 "R2_coefficients", "deviance", "deviance_true", "neg2loglik",
                 "nonzero_coefficients", "lambda_L0", "GIC", "GIC_rows", "note")
results_df <- do.call(rbind, lapply(results, function(x) {
  missing <- setdiff(common_cols, names(x))
  for (nm in missing) x[[nm]] <- NA
  x[common_cols]
}))

results_df <- results_df[order(results_df$case, is.na(results_df$seconds), results_df$deviance_true), ]
references <- rbind(score_reference("single", single), score_reference("multichannel", multi))
for (nm in c("GIC_beta_true", "GIC_rows_beta_true", "objective_lower_bound")) {
  results_df[[nm]] <- references[[nm]][match(results_df$case, references$case)]
}
timings <- setNames(results_df$seconds, results_df$method)
write.csv(results_df, file.path(out_dir, "readme_benchmark_results.csv"), row.names = FALSE)
write.csv(references, file.path(out_dir, "readme_benchmark_references.csv"), row.names = FALSE)
saveRDS(list(config = cfg, session_info = sessionInfo(),
             results = results_df, references = references,
             tuning_paths = tuning_paths,
             single = single, multi = multi,
             single_fits = single_fits, multi_fits = multi_fits),
        file.path(out_dir, "readme_benchmark_results.rds"))

if (length(single_fits)) {
  plot_single_signal(single, single_fits, file.path(out_dir, "single_channel_signal_zoom.png"),
                     from = cfg$plot_from, size = cfg$plot_n, timings = timings)
  plot_coefficients("Single-channel", single, single_fits, file.path(out_dir, "single_channel_coefficients.png"), timings = timings)
}
if (length(multi_fits)) {
  plot_coefficients("Multichannel", multi, multi_fits, file.path(out_dir, "multichannel_coefficients.png"), timings = timings)
  plot_multi_support(multi, multi_fits, file.path(out_dir, "multichannel_support.png"), timings = timings)
  plot_multi_support(multi, multi_fits, file.path(out_dir, "multichannel_support_zoom.png"),
                    from = cfg$plot_from, size = cfg$plot_n, timings = timings)
}
plot_summary(results_df, file.path(out_dir, "benchmark_summary.png"))

print(results_df)
print(references)
message("Wrote outputs to: ", normalizePath(out_dir, winslash = "/"))
}
