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
#' 
#' @details
#' The algorithm is inspired by the L0-penalised blind deconvolution approach
#' of de Rooi & Eilers (2011) for chromatographic pulse trains, but it
#' leverages the general machinery of \pkg{L0glm}. This results in the following
#' improvements:
#' \itemize{
#'   \item No manual tuning of the regularisation parameter is required; instead
#'         it can be preset to approximate information criteria such as GIC, MEBIC, AIC,
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
#' @param kernel_est Initial starting estimate of PSF (L1-normalised internally).
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
#'     \item{\code{g_hat}}{Final PSF estimate (L1-normalised).}
#'     \item{\code{deviance_trace}}{Deviance value after each iteration.}
#'   }
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
                          L = if (!is.null(kernel_est)) length(kernel_est) else 41,
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
                          tol = 1e-3)
{
  if (is.vector(y)) y <- as.matrix(y)
  n  <- nrow(y)
  offset <- NULL
  
  ## --- initial PSF -------------------------------------------------------
  g_hat <- if (is.null(kernel_est)) {
    kernelfunc(1:L)
  } else kernel_est
  g_hat <- g_hat / sum(g_hat)
  
  dev_trace <- numeric()
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
    est_baseline <- (it==1)&(!is.null(baseline))     # should baseline be estimated?
    fit <- L0glm_sparse_deco(y = if (est_baseline|is.null(baseline)) y else pmax(y-offset, 0), 
                             family = family, 
                             X = Xg,
                             baseline = if (est_baseline) baseline else NULL,
                             IC = IC, 
                             grouped = grouped)
    a_hat <- fit$coef_spikes                          # estimated spike train
    a_baseline <- fit$coef_baseline                   # estimated baseline coefs
    if (est_baseline) { X_baseline <- fit$X_baseline         # design matrix for baseline 
                        if (!is.null(X_baseline)) { offset <- X_baseline %*% a_baseline } else { 
                           offset <- NULL }
    }
    
    ## 2. PSF update -------------------------------------------------------
    g_new <- L0glm_deco_estkernel(y = if (est_baseline|is.null(baseline)) y else pmax(y-offset, 0), 
                                  y_deconv = a_hat,
                                  offset = if (est_baseline) offset else NULL,
                                  L = L,
                                  kernel_est = g_hat,
                                  family = family,
                                  kernel_smoothfun = kernel_smoothfun,
                                  perc = perc,
                                  IC = IC)$kernel
    
    ## 3. convergence ------------------------------------------------------
    delta <- sum(abs(g_new - g_hat)) / sum(abs(g_hat))
    g_hat <- g_new
    y_hat <- Xg %*% a_hat
    if (!is.null(offset)) y_hat <- y_hat + offset
    dev   <- sum(family$dev.resids(y = as.vector(y), 
                                   mu = pmax(as.vector(y_hat), .Machine$double.eps), # to avoid Inf values with poisson
                                   wt = 1))
    dev_trace[it] <- dev
    if (delta < tol) break
  }
  
  list(a_hat = a_hat,
       g_hat = g_hat,
       deviance_trace = dev_trace[seq_len(it)])
}
