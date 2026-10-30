#' @name estimate_psf_from_std
#' @title Estimate a time-varying PSF matrix from an alkane ladder TIC
#'
#' @description
#' Estimates the time-varying chromatographic point-spread function (PSF) from
#' the total ion chromatogram of an alkane ladder standard and returns a sparse
#' convolution design matrix.  The baseline is estimated as a low-degree
#' B-spline quantile regression underapproximation.  Major ladder peaks are
#' detected in the baseline-corrected TIC, each local peak profile is smoothed
#' as a non-negative, centroid-centred, log-concave unit-sum PSF, and the local
#' log-PSFs are interpolated or extrapolated over retention time.
#'
#' The returned sparse matrix \code{X} has one column per requested output
#' retention time.  Column \code{j} contains the PSF centred at
#' \code{output_timepoints[j]} and is normalised to sum to one.  It can be
#' passed as the \code{X} argument to \code{\link{L0glm_sparse_deco}} for
#' deconvolution with a time-varying PSF.
#'
#' @param y Numeric vector containing the alkane ladder total ion chromatogram.
#'   If \code{time} is omitted and \code{names(y)} are numeric, the names are
#'   used as retention times.
#' @param family GLM family object with an identity link, or one of
#'   \code{"poisson"}, \code{"gaussian"} or \code{"gamma"}.  The family is used
#'   to validate the intended non-negative identity-link noise model and to set
#'   family-aware weights for log-concave local PSF smoothing.
#' @param baseline_df Integer; degrees of freedom / number of spline knots for
#'   the low-dimensional B-spline baseline underapproximation.
#' @param output_timepoints Numeric vector of retention times for the rows and
#'   columns of the returned PSF matrix.  Defaults to the input \code{time}.
#' @param time Optional numeric vector of retention times for \code{y}.
#' @param baseline_tau Quantile used for the B-spline baseline
#'   underapproximation.  Small values such as \code{0.02} to \code{0.10} are
#'   useful for chromatograms with positive peaks.
#' @param n_peaks Optional integer number of major ladder peaks to keep.  For a
#'   C9--C40 ladder this is typically \code{32}.  If omitted, all sufficiently
#'   prominent local maxima after non-maximum suppression are used.
#' @param peak_min_distance Optional minimum distance, in retention-time units,
#'   between selected peaks.  If omitted, a conservative value is inferred from
#'   the chromatogram range.
#' @param peak_height Optional minimum baseline-corrected peak height.  If
#'   omitted, a robust threshold is inferred from the baseline residuals.
#' @param peak_snr Signal-to-noise multiplier used when \code{peak_height} is
#'   inferred automatically.
#' @param support_fraction Fraction of the peak apex height used to trim local
#'   peak windows.  The trimming is also bounded by neighbouring selected peaks.
#' @param psf_mass Fraction of each local peak's baseline-corrected signal mass
#'   used to infer the automatic PSF support.  The default uses the central
#'   99.9 percent of each detected peak and avoids very long low-intensity tails
#'   dominating \code{kernel_half_width}.
#' @param support_quantile Quantile of detected local peak half-widths used to
#'   infer \code{kernel_half_width} when it is not supplied.  Larger values keep
#'   more tail support but make the returned sparse matrix less sparse.
#' @param kernel_half_width Optional PSF half-width in retention-time units.
#'   If omitted, it is estimated from the detected peak windows.
#' @param psf_length Optional odd integer number of support points used for the
#'   local PSF grid.  If omitted, it is inferred from \code{kernel_half_width}
#'   and the median spacing of \code{output_timepoints}.
#' @param psf_nknots Optional number of knots for the log-concave local PSF
#'   smoothing.  If omitted, a small value is inferred from \code{psf_length}.
#' @param psf_tau Quantile used by the constrained log-spline smoother for each
#'   local PSF.  Values above \code{0.5} fit a more upper-envelope-like
#'   log-concave curve and tend to preserve heavier tails.
#' @param psf_weight_power Exponent applied to the local normalised peak profile
#'   when constructing log-spline weights.  The previous Poisson-like weighting
#'   corresponds to \code{1}; the default \code{0.5} downweights tails less and
#'   usually gives more realistic early alkane peak tails.  Use \code{0} for
#'   equal weights.
#' @param log_floor Positive relative floor used before interpolating local
#'   PSFs on the log scale.
#' @param drop_tol Relative tolerance for dropping tiny PSF entries from the
#'   sparse matrix.
#' @param return_details Logical; if \code{TRUE}, return a list with \code{X}
#'   and diagnostic details.  If \code{FALSE}, return \code{X} directly with the
#'   same diagnostics attached as attributes.
#' @param verbose Logical; print a short summary of the fitted baseline, peak
#'   detection, and PSF support.
#'
#' @return By default, a sparse \code{dgCMatrix} with dimensions
#'   \code{length(output_timepoints) x length(output_timepoints)}.  With
#'   \code{return_details = TRUE}, a list with components \code{X},
#'   \code{baseline}, \code{residual}, \code{peak_indices},
#'   \code{peak_timepoints}, \code{peak_heights}, \code{psf_offsets},
#'   \code{psf_kernels}, and \code{family}.
#'
#' @example examples/estimate_psf_from_std.R
#'
#' @seealso \code{\link{L0glm_sparse_deco}}, \code{\link{log_conc_smooth}},
#'   \code{\link{convmatrix}}
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom cobs cobs
#' @importFrom stats approx mad median na.pass predict quantile
#' @export
estimate_psf_from_std <- function(y,
                                  family = stats::poisson(link = "identity"),
                                  baseline_df = 10L,
                                  output_timepoints = NULL,
                                  time = NULL,
                                  baseline_tau = 0.05,
                                  n_peaks = NULL,
                                  peak_min_distance = NULL,
                                  peak_height = NULL,
                                  peak_snr = 5,
                                  support_fraction = 0.005,
                                  psf_mass = 0.999,
                                  support_quantile = 0.95,
                                  kernel_half_width = NULL,
                                  psf_length = NULL,
                                  psf_nknots = NULL,
                                  psf_tau = 0.75,
                                  psf_weight_power = 0.5,
                                  log_floor = 1e-12,
                                  drop_tol = 1e-10,
                                  return_details = FALSE,
                                  verbose = FALSE) {
  y_names <- names(y)
  y <- .psf_as_numeric_vector(y, "y")
  if (length(y) < 5L) stop("Need at least five TIC points.")

  if (is.null(time)) {
    if (!is.null(y_names) && length(y_names) == length(y)) {
      time <- suppressWarnings(as.numeric(y_names))
      if (any(!is.finite(time))) time <- NULL
    }
    if (is.null(time)) time <- seq_along(y)
  }
  time <- .psf_as_numeric_vector(time, "time")
  if (length(time) != length(y)) stop("'time' must have the same length as 'y'.")
  if (any(!is.finite(time)) || any(!is.finite(y))) {
    stop("'time' and 'y' must contain only finite values.")
  }

  ord <- order(time)
  time <- time[ord]
  y <- y[ord]
  if (any(diff(time) <= 0)) stop("'time' must contain unique, increasing values.")
  y <- pmax(y, 0)

  if (is.null(output_timepoints)) {
    output_timepoints <- time
  }
  output_timepoints <- .psf_as_numeric_vector(output_timepoints, "output_timepoints")
  if (length(output_timepoints) < 3L) {
    stop("'output_timepoints' must contain at least three values.")
  }
  if (any(!is.finite(output_timepoints)) || any(diff(output_timepoints) <= 0)) {
    stop("'output_timepoints' must contain finite, strictly increasing values.")
  }

  family <- .psf_identity_family(family)

  baseline_df <- as.integer(round(baseline_df))
  if (!is.finite(baseline_df) || baseline_df < 1L) {
    stop("'baseline_df' must be a positive integer.")
  }
  if (!is.numeric(baseline_tau) || length(baseline_tau) != 1L ||
      baseline_tau <= 0 || baseline_tau >= 1) {
    stop("'baseline_tau' must be a single number in (0, 1).")
  }
  if (!is.null(n_peaks)) {
    n_peaks <- as.integer(round(n_peaks))
    if (!is.finite(n_peaks) || n_peaks < 1L) {
      stop("'n_peaks' must be a positive integer when supplied.")
    }
  }
  if (!is.numeric(peak_snr) || length(peak_snr) != 1L || peak_snr <= 0) {
    stop("'peak_snr' must be a positive number.")
  }
  if (!is.numeric(support_fraction) || length(support_fraction) != 1L ||
      support_fraction <= 0 || support_fraction >= 1) {
    stop("'support_fraction' must be a single number in (0, 1).")
  }
  if (!is.numeric(psf_mass) || length(psf_mass) != 1L ||
      psf_mass <= 0 || psf_mass >= 1) {
    stop("'psf_mass' must be a single number in (0, 1).")
  }
  if (!is.numeric(support_quantile) || length(support_quantile) != 1L ||
      support_quantile <= 0 || support_quantile > 1) {
    stop("'support_quantile' must be a single number in (0, 1].")
  }
  if (!is.numeric(psf_tau) || length(psf_tau) != 1L ||
      psf_tau <= 0 || psf_tau >= 1) {
    stop("'psf_tau' must be a single number in (0, 1).")
  }
  if (!is.numeric(psf_weight_power) || length(psf_weight_power) != 1L ||
      psf_weight_power < 0 || psf_weight_power > 1) {
    stop("'psf_weight_power' must be a single number in [0, 1].")
  }
  if (!is.numeric(log_floor) || length(log_floor) != 1L ||
      log_floor <= 0 || log_floor >= 1) {
    stop("'log_floor' must be a single number in (0, 1).")
  }
  if (!is.numeric(drop_tol) || length(drop_tol) != 1L ||
      drop_tol < 0 || drop_tol >= 1) {
    stop("'drop_tol' must be a single number in [0, 1).")
  }

  baseline <- .psf_baseline_cobs(time, y, df = baseline_df, tau = baseline_tau)
  residual <- pmax(y - baseline, 0)
  noise <- .psf_noise_scale(y = y, baseline = baseline, family = family)

  peak_indices <- .psf_detect_major_peaks(
    time = time,
    residual = residual,
    noise = noise,
    n_peaks = n_peaks,
    peak_min_distance = peak_min_distance,
    peak_height = peak_height,
    peak_snr = peak_snr
  )
  if (length(peak_indices) < 1L) {
    stop("No major peaks could be detected in the baseline-corrected TIC.")
  }

  windows <- .psf_peak_windows(
    time = time,
    residual = residual,
    peak_indices = peak_indices,
    noise = noise,
    support_fraction = support_fraction,
    psf_mass = psf_mass
  )
  windows <- windows[vapply(windows, function(w) length(w$signal) >= 3L &&
                              sum(w$signal) > 0, logical(1))]
  if (length(windows) < 1L) {
    stop("Detected peaks did not contain enough positive signal to estimate a PSF.")
  }

  dt_out <- stats::median(diff(output_timepoints))
  if (!is.finite(dt_out) || dt_out <= 0) {
    stop("Could not infer a positive output-time spacing.")
  }

  if (is.null(kernel_half_width)) {
    half_widths <- vapply(windows, function(w) w$half_width, numeric(1))
    kernel_half_width <- as.numeric(stats::quantile(
      half_widths[is.finite(half_widths) & half_widths > 0],
      probs = support_quantile,
      names = FALSE,
      type = 8
    ))
    if (!is.finite(kernel_half_width) || kernel_half_width <= 0) {
      kernel_half_width <- 5 * dt_out
    }
  } else {
    kernel_half_width <- as.numeric(kernel_half_width)
    if (length(kernel_half_width) != 1L || !is.finite(kernel_half_width) ||
        kernel_half_width <= 0) {
      stop("'kernel_half_width' must be a positive number when supplied.")
    }
  }

  if (is.null(psf_length)) {
    psf_length <- 2L * ceiling(kernel_half_width / dt_out) + 1L
  } else {
    psf_length <- as.integer(round(psf_length))
    if (!is.finite(psf_length) || psf_length < 3L) {
      stop("'psf_length' must be an integer >= 3 when supplied.")
    }
    if (psf_length %% 2L == 0L) psf_length <- psf_length + 1L
  }

  if (is.null(psf_nknots)) {
    psf_nknots <- max(4L, min(25L, round(psf_length / 5)))
  } else {
    psf_nknots <- as.integer(round(psf_nknots))
    if (!is.finite(psf_nknots) || psf_nknots < 3L) {
      stop("'psf_nknots' must be an integer >= 3 when supplied.")
    }
  }

  psf_offsets <- seq(-kernel_half_width, kernel_half_width, length.out = psf_length)
  psf_kernels <- .psf_estimate_peak_kernels(
    windows = windows,
    offsets = psf_offsets,
    family = family,
    nknots = psf_nknots,
    tau = psf_tau,
    weight_power = psf_weight_power
  )
  keep <- Matrix::rowSums(Matrix::Matrix(psf_kernels, sparse = FALSE)) > 0
  psf_kernels <- psf_kernels[keep, , drop = FALSE]
  windows <- windows[keep]
  if (nrow(psf_kernels) < 1L) {
    stop("Log-concave smoothing failed for all detected peak windows.")
  }

  peak_timepoints <- vapply(windows, `[[`, numeric(1), "peak_time")
  peak_heights <- vapply(windows, `[[`, numeric(1), "peak_height")
  peak_indices <- vapply(windows, `[[`, integer(1), "peak_index")

  X <- .psf_build_sparse_matrix(
    output_timepoints = output_timepoints,
    peak_timepoints = peak_timepoints,
    psf_offsets = psf_offsets,
    psf_kernels = psf_kernels,
    log_floor = log_floor,
    drop_tol = drop_tol
  )

  details <- list(
    time = time,
    y = y,
    baseline = baseline,
    residual = residual,
    noise = noise,
    peak_indices = peak_indices,
    peak_timepoints = peak_timepoints,
    peak_heights = peak_heights,
    psf_offsets = psf_offsets,
    psf_kernels = psf_kernels,
    family = family,
    baseline_df = baseline_df,
    baseline_tau = baseline_tau,
    psf_mass = psf_mass,
    support_quantile = support_quantile,
    kernel_half_width = kernel_half_width,
    psf_length = psf_length,
    psf_nknots = psf_nknots,
    psf_tau = psf_tau,
    psf_weight_power = psf_weight_power
  )

  attr(X, "baseline") <- baseline
  attr(X, "residual") <- residual
  attr(X, "peak_indices") <- peak_indices
  attr(X, "peak_timepoints") <- peak_timepoints
  attr(X, "peak_heights") <- peak_heights
  attr(X, "psf_offsets") <- psf_offsets
  attr(X, "psf_kernels") <- psf_kernels
  attr(X, "family") <- family

  if (verbose) {
    message("Detected ", length(peak_timepoints), " major peaks; PSF length = ",
            psf_length, ", half-width = ", signif(kernel_half_width, 4), ".")
  }

  if (return_details) {
    details$X <- X
    return(details[c("X", setdiff(names(details), "X"))])
  }
  X
}

.psf_as_numeric_vector <- function(x, name) {
  x_names <- names(x)
  if (is.factor(x)) x <- as.character(x)
  if (is.matrix(x) || inherits(x, "Matrix")) {
    if (length(dim(x)) != 2L || min(dim(x)) != 1L) {
      stop("'", name, "' must be a numeric vector or a one-column/one-row matrix.")
    }
    x <- as.vector(x)
  }
  if (!is.numeric(x)) stop("'", name, "' must be numeric.")
  out <- as.numeric(x)
  if (!is.null(x_names) && length(x_names) == length(out)) names(out) <- x_names
  out
}

.psf_identity_family <- function(family) {
  if (is.character(family)) {
    family <- match.arg(tolower(family), c("poisson", "gaussian", "gamma"))
    family <- switch(
      family,
      poisson = stats::poisson(link = "identity"),
      gaussian = stats::gaussian(link = "identity"),
      gamma = stats::Gamma(link = "identity")
    )
  } else if (is.function(family) && is.null(family$family)) {
    family <- family(link = "identity")
  }

  if (!is.list(family) || is.null(family$family) || is.null(family$link) ||
      is.null(family$variance)) {
    stop("'family' must be a GLM family object or one of 'poisson', 'gaussian', 'gamma'.")
  }
  if (!identical(family$link, "identity")) {
    stop("'family' must use an identity link.")
  }
  supported <- c("poisson", "gaussian", "Gamma")
  if (!family$family %in% supported) {
    warning("Family '", family$family, "' is not one of poisson, gaussian or Gamma; ",
            "continuing with its identity-link variance function.", call. = FALSE)
  }
  family
}

.psf_baseline_cobs <- function(time, y, df, tau) {
  if (df <= 1L || length(unique(y)) < 2L) {
    return(rep(max(0, stats::quantile(y, tau, names = FALSE)), length(y)))
  }

  nknots <- max(4L, min(length(y) - 1L, df))
  knots <- seq(min(time), max(time), length.out = nknots)
  fit <- try(
    cobs::cobs(
      x = time,
      y = y,
      tau = tau,
      degree = 2,
      knots = knots,
      nknots = nknots,
      knots.add = FALSE,
      repeat.delete.add = FALSE,
      lambda = 0
    ),
    silent = TRUE
  )

  if (inherits(fit, "try-error")) {
    warning("B-spline quantile baseline fit failed; using a constant lower quantile baseline.",
            call. = FALSE)
    return(rep(max(0, stats::quantile(y, tau, names = FALSE)), length(y)))
  }

  baseline <- try(
    stats::predict(fit, z = time, na.action = stats::na.pass)[, "fit"],
    silent = TRUE
  )
  if (inherits(baseline, "try-error") || any(!is.finite(baseline))) {
    warning("B-spline quantile baseline prediction failed; using a constant lower quantile baseline.",
            call. = FALSE)
    baseline <- rep(max(0, stats::quantile(y, tau, names = FALSE)), length(y))
  }
  pmax(as.numeric(baseline), 0)
}

.psf_noise_scale <- function(y, baseline, family) {
  resid <- y - baseline
  q50 <- stats::quantile(resid, 0.50, names = FALSE)
  lower <- resid[resid <= q50]
  noise <- stats::mad(lower, center = stats::median(lower), constant = 1.4826,
                      na.rm = TRUE)
  if (!is.finite(noise) || noise <= 0) {
    mu <- pmax(baseline, .Machine$double.eps)
    var_mu <- family$variance(mu)
    noise <- sqrt(stats::median(var_mu[is.finite(var_mu) & var_mu > 0]))
  }
  if (!is.finite(noise) || noise <= 0) noise <- .Machine$double.eps
  noise
}

.psf_detect_major_peaks <- function(time,
                                    residual,
                                    noise,
                                    n_peaks,
                                    peak_min_distance,
                                    peak_height,
                                    peak_snr) {
  n <- length(residual)
  smooth_resid <- .psf_running_mean(residual, width = 3L)
  candidates <- which(smooth_resid[2L:(n - 1L)] > smooth_resid[1L:(n - 2L)] &
                        smooth_resid[2L:(n - 1L)] >= smooth_resid[3L:n]) + 1L
  if (!length(candidates)) candidates <- which.max(residual)

  if (is.null(peak_min_distance)) {
    peak_min_distance <- max(3 * stats::median(diff(time)), diff(range(time)) / 150)
  }
  if (!is.numeric(peak_min_distance) || length(peak_min_distance) != 1L ||
      !is.finite(peak_min_distance) || peak_min_distance <= 0) {
    stop("'peak_min_distance' must be a positive number when supplied.")
  }

  if (is.null(peak_height)) {
    nz <- residual[residual > 0]
    q_height <- if (length(nz)) stats::quantile(nz, 0.75, names = FALSE) else 0
    peak_height <- max(peak_snr * noise, q_height)
  } else {
    peak_height <- as.numeric(peak_height)
    if (length(peak_height) != 1L || !is.finite(peak_height) || peak_height < 0) {
      stop("'peak_height' must be a non-negative number when supplied.")
    }
  }

  candidates <- candidates[residual[candidates] >= peak_height]
  if (!length(candidates)) {
    n_fallback <- if (is.null(n_peaks)) 10L else n_peaks
    candidates <- order(residual, decreasing = TRUE)[seq_len(min(n, n_fallback))]
  }

  candidates <- candidates[order(residual[candidates], decreasing = TRUE)]
  selected <- integer(0)
  for (idx in candidates) {
    if (!length(selected) ||
        all(abs(time[idx] - time[selected]) >= peak_min_distance)) {
      selected <- c(selected, idx)
      if (!is.null(n_peaks) && length(selected) >= n_peaks) break
    }
  }

  selected[order(time[selected])]
}

.psf_peak_windows <- function(time,
                              residual,
                              peak_indices,
                              noise,
                              support_fraction,
                              psf_mass) {
  n <- length(time)
  lapply(seq_along(peak_indices), function(k) {
    idx <- peak_indices[k]
    left_limit_time <- if (k == 1L) min(time) else mean(time[peak_indices[k - 1L]], time[idx])
    right_limit_time <- if (k == length(peak_indices)) {
      max(time)
    } else {
      mean(time[idx], time[peak_indices[k + 1L]])
    }
    left_limit <- .psf_lower_bound(time, left_limit_time)
    right_limit <- .psf_upper_bound(time, right_limit_time)
    left_limit <- max(1L, min(left_limit, idx))
    right_limit <- min(n, max(right_limit, idx))

    edge <- max(noise, support_fraction * residual[idx])
    left <- idx
    while (left > left_limit && residual[left] > edge) left <- left - 1L
    right <- idx
    while (right < right_limit && residual[right] > edge) right <- right + 1L

    min_half_points <- 2L
    left <- max(left_limit, min(left, idx - min_half_points))
    right <- min(right_limit, max(right, idx + min_half_points))
    x <- time[left:right]
    signal <- residual[left:right]
    local_floor <- stats::quantile(signal, 0.05, names = FALSE)
    signal <- pmax(signal - local_floor, 0)
    if (sum(signal) <= 0) signal <- pmax(residual[left:right], 0)
    centroid <- sum(x * signal) / sum(signal)
    q <- .psf_weighted_quantile(
      x = x,
      w = signal,
      probs = c((1 - psf_mass) / 2, 1 - (1 - psf_mass) / 2)
    )
    half_width <- max(abs(q - centroid))
    if (!is.finite(half_width) || half_width <= 0) {
      half_width <- max(abs(range(x) - centroid))
    }

    list(
      x = x,
      signal = signal,
      centroid = centroid,
      half_width = half_width,
      peak_index = idx,
      peak_time = time[idx],
      peak_height = residual[idx]
    )
  })
}

.psf_weighted_quantile <- function(x, w, probs) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  x <- x[keep]
  w <- w[keep]
  if (!length(x) || sum(w) <= 0) return(rep(NA_real_, length(probs)))
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)
  stats::approx(cw, x, xout = probs, ties = "ordered", rule = 2)$y
}

.psf_estimate_peak_kernels <- function(windows,
                                       offsets,
                                       family,
                                       nknots,
                                       tau,
                                       weight_power) {
  kernels <- lapply(windows, function(w) {
    x_centered <- w$x - w$centroid
    raw <- stats::approx(
      x = x_centered,
      y = w$signal,
      xout = offsets,
      yleft = 0,
      yright = 0,
      ties = "ordered"
    )$y
    if (sum(raw) <= 0 || all(!is.finite(raw))) return(rep(0, length(offsets)))
    raw <- pmax(raw, 0)
    raw <- raw / sum(raw)
    smoothed <- .psf_log_concave_smooth(raw,
                                        family = family,
                                        nknots = nknots,
                                        tau = tau,
                                        weight_power = weight_power)
    smoothed <- shift_kernel(smoothed, perc = "mean", method = "linear")
    if (sum(smoothed) > 0) smoothed / sum(smoothed) else smoothed
  })
  do.call(rbind, kernels)
}

.psf_log_concave_smooth <- function(y,
                                    family,
                                    nknots,
                                    tau,
                                    weight_power) {
  y <- pmax(as.numeric(y), 0)
  if (length(y) < 3L || sum(y) <= 0) return(y)

  n <- length(y)
  nknots <- max(3L, min(as.integer(nknots), n))
  x <- seq_len(n)
  offs <- max(y) * 1e-6
  if (!is.finite(offs) || offs <= 0) offs <- .Machine$double.eps

  weights <- (y + offs)^weight_power
  weights <- pmax(weights, .Machine$double.eps)
  weights <- weights / sum(weights) * length(weights)

  fit <- try(
    cobs::cobs(
      x = x,
      y = log10(y + offs),
      constraint = "concave",
      w = weights,
      knots = seq(min(x), max(x), length.out = nknots),
      nknots = nknots,
      knots.add = FALSE,
      repeat.delete.add = FALSE,
      tau = tau,
      degree = 2,
      lambda = 0
    ),
    silent = TRUE
  )

  if (inherits(fit, "try-error")) {
    out <- log_conc_smooth(y, tau = tau, nknots = nknots)
    return(out / sum(out))
  }

  out <- try(10^stats::predict(fit, z = x, na.action = stats::na.pass)[, "fit"] + offs,
             silent = TRUE)
  if (inherits(out, "try-error") || any(!is.finite(out))) {
    out <- log_conc_smooth(y, tau = tau, nknots = nknots)
    return(out / sum(out))
  }
  out <- pmax(as.numeric(out) - min(out, na.rm = TRUE), 0)
  if (sum(out) > 0) out / sum(out) else y / sum(y)
}

.psf_build_sparse_matrix <- function(output_timepoints,
                                     peak_timepoints,
                                     psf_offsets,
                                     psf_kernels,
                                     log_floor,
                                     drop_tol,
                                     normalize = TRUE,
                                     smooth = FALSE) {
  n <- length(output_timepoints)
  half_width <- max(abs(psf_offsets))
  if (log_floor == 0) {
    # Linear log-tail continuation preserves concavity when supports differ.
    log_kernels <- t(apply(psf_kernels, 1L, .psf_log_tail))
  } else {
    floor_abs <- pmax(psf_kernels, apply(psf_kernels, 1L, max) * log_floor)
    log_kernels <- log(floor_abs)
  }
  common <- which(colSums(psf_kernels) > 0)

  i_list <- vector("list", n)
  j_list <- vector("list", n)
  x_list <- vector("list", n)

  for (j in seq_len(n)) {
    center <- output_timepoints[j]
    log_kernel <- .psf_interp_log_kernel(center, peak_timepoints, log_kernels, smooth)
    kernel_grid <- exp(log_kernel - max(log_kernel))
    if (log_floor == 0) kernel_grid[seq_along(kernel_grid) < min(common) |
                                     seq_along(kernel_grid) > max(common)] <- 0
    if (sum(kernel_grid) <= 0 || any(!is.finite(kernel_grid))) next
    kernel_grid <- kernel_grid / sum(kernel_grid)

    lo <- .psf_lower_bound(output_timepoints, center - half_width)
    hi <- .psf_upper_bound(output_timepoints, center + half_width)
    if (lo > hi || lo > n || hi < 1L) next
    rows <- lo:hi
    vals <- stats::approx(
      x = psf_offsets,
      y = kernel_grid,
      xout = output_timepoints[rows] - center,
      yleft = 0,
      yright = 0,
      ties = "ordered"
    )$y
    vals <- pmax(as.numeric(vals), 0)
    if (length(vals) && max(vals) > 0 && drop_tol > 0) {
      vals[vals < max(vals) * drop_tol] <- 0
    }
    val_sum <- sum(vals)
    if (!is.finite(val_sum) || val_sum <= 0) {
      nearest <- which.min(abs(output_timepoints - center))
      rows <- nearest
      vals <- 1
    } else if (normalize) {
      vals <- vals / val_sum
    }
    keep <- vals > 0 & is.finite(vals)
    i_list[[j]] <- rows[keep]
    j_list[[j]] <- rep.int(j, sum(keep))
    x_list[[j]] <- vals[keep]
  }

  Matrix::sparseMatrix(
    i = unlist(i_list, use.names = FALSE),
    j = unlist(j_list, use.names = FALSE),
    x = unlist(x_list, use.names = FALSE),
    dims = c(n, n)
  )
}

.psf_interp_log_kernel <- function(center, peak_timepoints, log_kernels, smooth = FALSE) {
  n_peaks <- length(peak_timepoints)
  if (n_peaks == 1L || center <= peak_timepoints[1L]) {
    return(log_kernels[1L, ])
  }
  if (center >= peak_timepoints[n_peaks]) {
    return(log_kernels[n_peaks, ])
  }
  left <- findInterval(center, peak_timepoints)
  right <- left + 1L
  w <- (center - peak_timepoints[left]) /
    (peak_timepoints[right] - peak_timepoints[left])
  if (smooth) w <- w * w * (3 - 2 * w)
  (1 - w) * log_kernels[left, ] + w * log_kernels[right, ]
}

.psf_log_tail <- function(kernel) {
  at <- which(kernel > 0)
  if (length(at) == 1L) return(-36 * abs(seq_along(kernel) - at))
  z <- stats::approx(at, log(kernel[at]), xout = seq_along(kernel), rule = 2)$y
  left <- seq_len(min(at) - 1L)
  right <- seq.int(max(at), length(kernel))
  z[left] <- z[min(at)] + (left - min(at)) * (z[at[2L]] - z[at[1L]])
  z[right] <- z[max(at)] + (right - max(at)) * (z[at[length(at)]] - z[at[length(at) - 1L]])
  z
}

.psf_running_mean <- function(x, width = 3L) {
  width <- as.integer(width)
  if (width <= 1L || length(x) < width) return(x)
  filt <- rep(1 / width, width)
  y <- as.numeric(stats::filter(x, filt, sides = 2L))
  missing <- !is.finite(y)
  y[missing] <- x[missing]
  y
}

.psf_lower_bound <- function(x, value) {
  idx <- findInterval(value, x)
  if (idx >= 1L && idx <= length(x) &&
      abs(x[idx] - value) <= sqrt(.Machine$double.eps) * max(1, abs(value))) {
    return(idx)
  }
  max(1L, idx + 1L)
}

.psf_upper_bound <- function(x, value) {
  min(length(x), findInterval(value, x))
}
