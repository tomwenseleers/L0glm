#' @name L0glm_blind_deco
#' @title Alternating sparse L0 glm powered blind deconvolution
#'
#' @description
#' Runs a simple block-coordinate (alternating) optimisation:
#' repeatedly
#' \enumerate{
#'   \item fix the PSF and estimate spikes with \code{L0glm_sparse_deco};
#'   \item fix the spikes and estimate the PSF with
#'         \code{L0glm_deco_estkernel};
#' }
#' until \code{maxiter} iterations have elapsed or the relative L1 change in 
#' the PSF drops below \code{tol}. This allows one to simultaneously determine 
#' the point spread function and spike locations of an observed blurred single-
#' or multichannel spike train corrupted by exponential-family noise 
#' (Poisson, Gaussian, Gamma, Inverse Gaussian, etc.).
#' The centering specified by \code{perc} sets the time origin of the estimated
#' PSF and spikes. Use the same convention for simulated reference kernels when
#' comparing spike locations or coefficients directly.
#' Set \code{perc = 0} to retain the estimated time origin without recentering.
#' 
#' @details
#' The algorithm is inspired by the L0-penalised blind deconvolution approach
#' of de Rooi & Eilers (2011) for chromatographic pulse trains, but it
#' leverages the general machinery of \pkg{L0glm}. This results in the following
#' improvements:
#' \itemize{
#'   \item No manual tuning of the regularisation parameter is required; instead
#'         it can be preset to approximate information criteria such as GIC, AIC,
#'         BIC, EBIC or MBIC
#'   \item Both single channel data with L0 best subset selection or 
#'         multichannel data with grouped L0 best subset of group selection 
#'         are allowed for  
#'    \item We allow the estimated convolution kernel \eqn{g} to be smoothed by
#'          a custom function (e.g. log-concave splines) 
#'   \item Any exponential-family noise model with an \emph{identity} link
#'         can be used (Poisson, Gaussian, Gamma, ...)
#'   \item Estimates of the PSF are shrunk towards an optional prior guess or
#'         the estimate from the previous iteration and warm starts are used, 
#'         which alongside PSF centroid centering results in improved 
#'         stability and better convergence behaviour  
#'   \item The underlying sparse designs are handled by Eigen's efficient sparse
#'         Cholesky solver, enabling large problems to be tackled quickly
#' }
#' 
#' In the non-negative-matrix-factorisation (NMF) literature this model
#' would be referred to as \emph{single-atom convolutional NMF} or \emph{sparse 
#' convolutional dictionary learning}, whereby a one-dimensional signal (or
#' multichannel matrix) \eqn{y} is approximated as \eqn{W(g)\,H}, where \eqn{W(g)}
#' is a banded Toeplitz matrix that holds time-shifted copies of
#' an estimated optimal dictionary atom/PSF \eqn{g}, and \eqn{H} is a
#' (potentially multi-column) matrix of sparse, non-negative activations
#' \eqn{a}.
#' However, relative to standard single-atom convolutional NMF we 
#' offer several extensions:
#' \itemize{
#'   \item \emph{Optimal sparsity}: we do not fix the sparsity level of the 
#'       activations \eqn{a} in advance; instead we let \pkg{L0glm} perform
#'       best-subset or best subset of group selection through an
#'       L\eqn{_{0}} or grouped L\eqn{_{0}} penalty,
#'   \item \emph{Smoothness and shape-constraints on \eqn{g}}: we allow 
#'       the estimated convolution kernel atom \eqn{g} to be smoothed by
#'       a custom function (e.g. log-concave splines)    
#'   \item \emph{General noise models}: any exponential-family likelihood
#'         with an \emph{identity} link is supported (Poisson, Gaussian,
#'         Gamma, inverse-Gaussian, ...), instead of only Gaussian
#'         (Euclidean-distance NMF) or Poisson (NMF with KL divergence)
#'   \item \emph{Speed}: the Toeplitz designs in the alternating updates are
#'         solved by sparse-Cholesky linear algebra rather than multiplicative 
#'         updates, giving orders-of-magnitude speed-ups on long signals
#' }
#' Methodologically this could be referred to as
#' \emph{generalised / deviance matrix factorisation}
#' (Kidzinski et al. 2022; Castiglione et al. 2024; Wang & Carvalho 2023), 
#' specialised to the case where \eqn{W} is a Toeplitz single-atom convolutional 
#' dictionary matrix and where the matrix factorisation is solved via
#' alternating block-coordinate descent, with each update being powered by 
#' a L\eqn{_{0}} or grouped L\eqn{_{0}}-penalised generalised linear model (GLM)
#' fit via adaptive ridge regularised IRLS.
#'
#' A classical baseline for Poisson deconvolution is the
#' Richardson-Lucy algorithm (Richardson 1972; Lucy 1974). In modern
#' terms, the blind deconvolution version of Richardson-Lucy 
#' equals Kullback-Leibler NMF with
#' multiplicative updates (Lee & Seung 1999) and 
#' factor matrix \eqn{W} constrained to be a convolution matrix.  
#' Empirically, \code{L0glm_blind_deco()} has much better support recovery,
#' attains much lower reconstruction error and is dramatically faster thanks to 
#' (i) closed-form sparse linear solves, (ii) the automatic selection of the penalty
#' parameter through information-criterion presets, (iii) semi-parametric shape
#' constraints on the estimated convolution kernel and (iv) the possibility
#' to exploit grouped-sparsity when several channels are deconvolved jointly.
#'
#' @inheritParams L0glm_sparse_deco
#' @param kernel_est Initial PSF (L1-normalised internally). For
#'   \code{nwindows > 1}, a matrix with one row per window can supply local
#'   starting PSFs of equal length.
#' @param nwindows Number of overlapping windows used to estimate a
#'   time-varying PSF. The default \code{1L} estimates one shared PSF. Larger
#'   values fit local PSFs, interpolate their log densities smoothly across
#'   scans, and refit the original outcomes using the resulting sparse design.
#'   For multichannel Poisson data, local PSF calibration uses the summed
#'   channels; the final fit retains every original channel.
#' @param perc       Passed to \code{\link{shift_kernel}} when the PSF is
#'                   updated. Either '"mean"' (shift the centroid, default) or a 
#'                   numeric percentile in [0, 1] to shift a given percentile in
#'                   the CDF of the estimated kernel.
#' @param kernel_smoothfun Optional smoother applied to the raw PSF in
#'                   each iteration (default: \code{\link{log_conc_smooth}}).
#' @param kernel_true Optional vector of length L with convolution kernel / PSF 
#'                    ground truth for cases where simulated data is used. 
#'                    Can be plotted alongside the estimated PSF in each 
#'                    iteration when \code{plot=TRUE}.
#' @param plot Boolean specifying whether to plot the estimated PSF in each iteration
#'             (alongside the true PSF if \code{kernel_true} is given).                                        
#' @param maxiter  Integer >= 1.  Maximum alternating steps. Default 15.
#' @param tol      Relative L1 change in the PSF that triggers early stop (default 1e-3).
#'
#' @return A list with components
#'   \describe{
#'     \item{\code{a_hat}}{Final spike estimate.}
#'     \item{\code{g_hat}}{Final unit-sum PSF, or a matrix with one PSF per
#'       window when \code{nwindows > 1}.}
#'     \item{\code{deviance_trace}}{Deviance value after each iteration.}
#'     \item{\code{fitted_values}, \code{fitted_baseline}, \code{coef_baseline}}{
#'       Fitted signal, background and background coefficients on the original data scale.}
#'     \item{\code{X}, \code{X_spikes}, \code{X_baseline}}{Final sparse PSF design
#'       (\code{X} and \code{X_spikes}) and background basis.}
#'     \item{\code{psf_timepoints}, \code{windows}}{For time-varying fits,
#'       calibration centers in scan coordinates and local diagnostics.}
#'     \item{\code{deviance}, \code{iterations}, \code{converged}}{
#'       Final deviance, number of alternating steps and PSF convergence flag.}
#'     \item{\code{best_iteration}, \code{selection_score}, \code{selection_trace}}{
#'       For a shared PSF, the selected update, its information criterion and
#'       the criterion at each completed candidate fit. Local diagnostics for
#'       time-varying PSFs are stored in \code{windows}.}
#'   }
#'   Background coefficients are fitted with the spikes in each iteration.
#'   The kernel update uses the fitted background as a channel-specific offset.
#'   Original observations are used in every fit. A final spike/background fit
#'   ensures that the returned coefficients and deviance use the returned PSF.
#'   Time-varying fits blend local log PSFs using nonnegative smoothstep weights,
#'   with linear log-tail continuation when their supports differ. This preserves
#'   log-concavity when the local PSFs are log-concave. Edge columns are truncated
#'   using the same convention as \code{convmatrix()}.
#'   The best completed spike/background fit according to \code{IC} is retained
#'   across PSF updates. \code{best_iteration} identifies its PSF update; zero
#'   denotes the supplied starting kernel. \code{selection_trace} records the
#'   criterion for every completed candidate fit.
#'   
#' @references
#'   \itemize{
#'     \item Richardson, W. H. (1972)\cr
#'           Bayesian-based iterative method of image restoration.\cr
#'           \emph{J. Opt. Soc. Am.} 62 (1): 55-59.
#'     \item Lucy, L. B. (1974)\cr
#'           An iterative technique for the rectification of observed
#'           distributions.\cr
#'           \emph{Astron. J.} 79: 745-754.
#'     \item Smaragdis, P. (2004)\cr
#'           Non-negative matrix factor
#'           deconvolution; extraction of multiple sound sources from
#'           monophonic inputs.\cr
#'           \emph{Proceedings of ICA 2004}.
#'     \item Lee, D. D., & Seung, H. S. (1999)\cr
#'           Learning the parts of objects by non-negative matrix
#'           factorization.\cr
#'           \emph{Nature} 401 (6755): 788-791.
#'     \item de Rooi, J., & Eilers, P. H. C. (2011)\cr
#'           Deconvolution of pulse trains with the L0 penalty.\cr
#'           \emph{Analytica Chimica Acta} 705 (1-2): 218-226.
#'     \item Kidzinski, L., Stumpe, M. C., & Hastie, T. (2022)\cr
#'           Generalized matrix factorization.\cr
#'           arXiv:2208.05987.
#'     \item Castiglione, F., van Kolk, E., & Chiappori, G. (2024)\cr
#'           Generalized non-negative matrix factorization for biological
#'           count data.\cr
#'           \emph{Bioinformatics} 40 (2): btae024.
#'     \item Wang, W., & Carvalho, C. M. (2023)\cr
#'           Deviance matrix factorization.\cr
#'           arXiv:2303.01517.
#'   }
#'
#' @seealso \code{\link{L0glm_sparse_deco}}, \code{\link{L0glm_deco_estkernel}}
#'
#' @example examples/L0glm_blind_deco.R
#' @importFrom stats setNames
#' @export
L0glm_blind_deco <- function(y,
                          family = poisson(identity),
                          kernel_est = NULL,
                          L = if (is.matrix(kernel_est)) ncol(kernel_est) else
                            if (!is.null(kernel_est)) length(kernel_est) else 41,
                          baseline = NULL,
                          nknots = max(4L, floor(nrow(as.matrix(y))/500)),
                          kernelfunc = function(x) convkernel("gaussian", L, list(sigma = 6)),
                          kernel_smoothfun = function (y) log_conc_smooth(y),
                          perc = "mean",
                          IC = "gic",
                          grouped = TRUE,
                          kernel_true = NULL,
                          plot = FALSE,
                          maxiter = 15,
                          tol = 1e-3,
                          nwindows = 1L, ...)
{
  response <- L0glm_deco_response(y, family)
  y <- response$y; family <- response$family
  n <- nrow(y)
  L0glm_deco_length(nwindows, "nwindows", 1L)
  if (!is.null(kernel_est)) L <- if (is.matrix(kernel_est)) ncol(kernel_est) else length(kernel_est)
  L0glm_deco_length(maxiter, "maxiter", 1L)
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0)
    stop("tol must be a finite nonnegative number.", call. = FALSE)
  if (nwindows > 1L) {
    return(L0glm_blind_windows(y, family, kernel_est, L, baseline, nknots,
      kernelfunc, kernel_smoothfun, perc, IC, grouped, kernel_true, plot,
      maxiter, tol, nwindows, ...))
  }
  if (is.matrix(kernel_est)) {
    if (nrow(kernel_est) != 1L) stop("kernel_est must have one row per window.", call. = FALSE)
    kernel_est <- as.numeric(kernel_est)
  }
  
  ## --- initial PSF -------------------------------------------------------
  g_hat <- if (is.null(kernel_est)) {
    L0glm_deco_length(L, "L", 3L)
    kernelfunc(seq_len(L))
  } else kernel_est
  g_hat <- L0glm_deco_kernel(g_hat, normalize = TRUE)
  L <- length(g_hat)
  
  dev_trace <- numeric()
  selection_trace <- numeric()
  best <- NULL
  criterion <- if (is.character(IC)) L0glm_ic_name(IC) else "gic"
  checkpoint <- function(fit, kernel, iteration) {
    score <- L0glm_model_criterion(fit$fit, criterion, mu_floor = 1e-8)
    selection_trace[iteration + 1L] <<- score
    if (is.finite(score) && (is.null(best) || score < best$score))
      best <<- list(fit = fit, kernel = kernel, iteration = iteration, score = score)
  }
  converged <- FALSE
  for (it in seq_len(maxiter)) {
    if (plot) {
      if (is.null(kernel_true)) { plot(g_hat, type = "l", 
                                       ylab = "Density",
                                       main = paste0("ESTIMATED PSF\niteration ", it-1), 
                                       col = "red", lwd = 3) } else {   
        plot(kernel_true, type="l", 
             ylab = "Density",
             main = paste0("TRUE (grey) & ESTIMATED (red) PSF\niteration ", it-1),
             col = "grey",lwd = 7)
        lines(g_hat, col = "red", lwd = 3) 
    }
    }
    
    ## 1. spike update -----------------------------------------------------
    Xg <- convmatrix(g_hat, n)              # banded design
    fit <- L0glm_sparse_deco(y = y,
                             family = family, 
                             X = Xg,
                             baseline = baseline,
                             nknots = nknots,
                             IC = IC, 
                             grouped = grouped, ...)
    a_hat <- fit$coef_spikes                          # estimated spike train
    checkpoint(fit, g_hat, it - 1L)
    
    ## 2. PSF update -------------------------------------------------------
    g_new <- L0glm_deco_estkernel(y = y,
                                  y_deconv = a_hat,
                                  offset = fit$fitted_baseline,
                                  L = L,
                                  kernel_est = g_hat,
                                  family = family,
                                  kernel_smoothfun = kernel_smoothfun,
                                  perc = perc,
                                  IC = IC)$kernel
    
    ## 3. convergence ------------------------------------------------------
    delta <- sum(abs(g_new - g_hat)) / sum(abs(g_hat))
    g_hat <- g_new
    y_hat <- convmatrix(g_hat, n) %*% a_hat + fit$fitted_baseline
    dev_trace[it] <- L0glm_deco_deviance(y, y_hat, family)
    if (delta < tol) { converged <- TRUE; break }
  }
  # Synchronize the final spike and background estimates with the returned PSF.
  fit <- L0glm_sparse_deco(y, family = family, kernel = g_hat, baseline = baseline,
                          nknots = nknots, IC = IC, grouped = grouped, ...)
  checkpoint(fit, g_hat, it)
  converged <- converged && (is.null(best) || sum(abs(best$kernel - g_hat)) < tol)
  if (!is.null(best)) { fit <- best$fit; g_hat <- best$kernel }
  dev_trace[it] <- fit$deviance
  
  list(a_hat = fit$coef_spikes,
       g_hat = g_hat,
       deviance_trace = dev_trace[seq_len(it)],
       fitted_values = fit$fitted_values,
       fitted_baseline = fit$fitted_baseline,
       coef_baseline = fit$coef_baseline,
       X = fit$X_spikes, X_spikes = fit$X_spikes, X_baseline = fit$X_baseline,
       deviance = fit$deviance, iterations = it,
       best_iteration = if (is.null(best)) it else best$iteration,
       criterion = criterion, selection_score = if (is.null(best)) NA_real_ else best$score,
       selection_trace = selection_trace,
       converged = converged)
}

L0glm_blind_windows <- function(y, family, kernel_est, L, baseline, nknots,
                                kernelfunc, kernel_smoothfun, perc, IC, grouped,
                                kernel_true, plot, maxiter, tol, nwindows, ...) {
  n <- nrow(y)
  L0glm_deco_length(L, "L", 3L)
  if (n / nwindows < L)
    stop("Each blind-deconvolution window must span at least L scans.", call. = FALSE)
  if (is.matrix(kernel_est) && (nrow(kernel_est) != nwindows || ncol(kernel_est) != L))
    stop("kernel_est must have one row per window and L columns.", call. = FALSE)
  centers <- .5 + (seq_len(nwindows) - .5) * n / nwindows
  half_width <- n / (2 * nwindows) + L / 2
  local_y <- if (ncol(y) > 1L && identical(family$family, "poisson"))
    matrix(rowSums(y), ncol = 1L) else y
  kernels <- matrix(0, nwindows, L)
  windows <- vector("list", nwindows)
  for (j in seq_len(nwindows)) {
    rows <- seq.int(max(1L, ceiling(centers[j] - half_width)),
                    min(n, floor(centers[j] + half_width)))
    initial <- if (is.matrix(kernel_est)) kernel_est[j, ] else kernel_est
    fit <- L0glm_blind_deco(local_y[rows, , drop = FALSE], family = family,
      kernel_est = initial, L = L, baseline = baseline,
      nknots = max(4L, ceiling(nknots / nwindows) + 3L), kernelfunc = kernelfunc,
      kernel_smoothfun = kernel_smoothfun, perc = perc, IC = IC, grouped = grouped,
      kernel_true = if (is.matrix(kernel_true)) kernel_true[j, ] else kernel_true,
      plot = plot, maxiter = maxiter, tol = tol, ...)
    kernels[j, ] <- fit$g_hat
    windows[[j]] <- list(scans = range(rows), iterations = fit$iterations,
      best_iteration = fit$best_iteration, selection_trace = fit$selection_trace,
      converged = fit$converged, deviance_trace = fit$deviance_trace)
  }
  rownames(kernels) <- format(centers, trim = TRUE)
  X <- .psf_build_sparse_matrix(seq_len(n), centers, seq_len(L) - (L + 1L) %/% 2L,
    kernels, log_floor = 0, drop_tol = 1e-12, normalize = FALSE, smooth = TRUE)
  fit <- L0glm_sparse_deco(y, family = family, X = X, baseline = baseline,
                           nknots = nknots, IC = IC, grouped = grouped, ...)
  list(a_hat = fit$coef_spikes, g_hat = kernels, psf_timepoints = centers,
    windows = windows, deviance_trace = fit$deviance,
    fitted_values = fit$fitted_values, fitted_baseline = fit$fitted_baseline,
    coef_baseline = fit$coef_baseline, X = fit$X_spikes, X_spikes = fit$X_spikes,
    X_baseline = fit$X_baseline, deviance = fit$deviance,
    iterations = vapply(windows, `[[`, 0L, "iterations"),
    converged = all(vapply(windows, `[[`, FALSE, "converged")))
}
