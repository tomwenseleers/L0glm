#' @name L0_blind_deco
#' @title Alternating sparse L0 glm powered blind deconvolution
#'
#' @description
#' Runs a simple block-coordinate (alternating) optimisation:
#' repeatedly
#' \enumerate{
#'   \item fix the PSF and estimate spikes with \code{L0_sparse_deco};
#'   \item fix the spikes and estimate the PSF with
#'         \code{L0_deco_estkernel};
#' }
#' until \code{maxiter} iterations have elapsed or the PSF change
#' drops below \code{tol}.
#'
#' @inheritParams L0_sparse_deco
#' @param kernel_est Initial starting estimate of PSF (L1-normalised internally).
#' @param p          Passed to \code{\link{shift_kernel}} when the PSF is
#'                   updated; `"mean"` recentres the centroid.
#' @param kernel_smoothfun Optional smoother applied to the raw PSF in
#'                   each iteration (default: \code{\link{log_conc_smooth}}).
#' @param kernel_true Optional vector with convolution kernel / PSF ground truth 
#'                    for cases where simulated data is used. Can be plotted 
#'                    alongside the estimated PSF in each iteration when 
#'                    \code{plot=TRUE}.
#' @param plot Boolean specifying whether to plot the estimated PSF in each iteration
#'             (alongside the true PSF if \code{kernel_true} is given).                                        
#' @param maxiter  Integer ≥ 1.  Maximum alternating steps. Default 30.
#' @param tol      Relative L1 change in the PSF that triggers early stop (default 1e-9).
#'
#' @return A list with components
#'   \describe{
#'     \item{\code{a_hat}}{Final spike estimate.}
#'     \item{\code{g_hat}}{Final PSF estimate (L1-normalised).}
#'     \item{\code{deviance_trace}}{Deviance value after each iteration.}
#'   }
#'
#' @seealso \code{\link{L0_sparse_deco}}, \code{\link{L0_deco_estkernel}}
#'
#' @examples
#' if (requireNamespace("L0glm", quietly = TRUE)) {
#'   set.seed(123)
#'   n       <- 500
#'   ks_true <- convkernel("bigaussian", 41,
#'                         list(sigma1 = 2, sigma2 = 5))
#'   a_true  <- rpois(n, 0.1)
#'   y       <- rpois(n, as.vector(convmatrix(ks_true, n) %*% a_true))
#'
#'   fit <- L0_blind_deco(y, L = 41, maxiter = 10)
#'   matplot(cbind(ks_true, fit$g_hat), type = "l", lty = 1,
#'           col = c("grey", "red"), main = "PSF: true vs. blind estimate")
#'   legend("topright", c("true", "blind-deconv"),
#'          col = c("grey", "red"), lty = 1, bty = "n")
#' }
#'
#' @importFrom stats setNames
#' @export
L0_blind_deco <- function(y,
                          family = poisson(identity),
                          kernel_est = NULL,
                          L = if (!is.null(kernel_est)) length(kernel_est) else 41,
                          kernelfunc = function(x) convkernel("gaussian", L, list(sigma = 6)),
                          kernel_smoothfun = function (y) log_conc_smooth(y),
                          p = "mean",
                          IC = "gic",
                          grouped = FALSE,
                          kernel_true = NULL,
                          plot = FALSE,
                          maxiter = 30,
                          tol = 1e-9)
{
  if (is.vector(y)) y <- as.matrix(y)
  n  <- nrow(y)
  
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
    Xg  <- convmatrix(g_hat, n)              # banded design
    a_hat <- L0_sparse_deco(y, family, X = Xg,
                            IC = IC, grouped = grouped)$y_deconv
    
    ## 2. PSF update -------------------------------------------------------
    g_new <- L0_deco_estkernel(y, a_hat,
                               L = L,
                               kernel_est = g_hat,
                               family = family,
                               kernel_smoothfun = kernel_smoothfun,
                               perc = p,
                               IC = IC)$kernel
    
    ## 3. convergence ------------------------------------------------------
    delta <- sum(abs(g_new - g_hat)) / sum(abs(g_hat))
    g_hat <- g_new
    dev   <- sum(family$dev.resids(y, as.vector(Xg %*% a_hat) + 1e-10, 1))
    dev_trace[it] <- dev
    if (delta < tol) break
  }
  
  list(a_hat = a_hat,
       g_hat = g_hat,
       deviance_trace = dev_trace[seq_len(it)])
}