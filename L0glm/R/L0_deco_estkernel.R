#' @name L0_deco_estkernel
#' @title Estimate the PSF given an (approximate) spike train
#'
#' @description
#' Solves the *transpose* problem:
#' given an (estimated) latent signal \eqn{a} and the observation \eqn{y},
#' find the PSF \eqn{g} that minimises the deviance in
#' \eqn{y \approx X(g)\,a}.  Regularisation is again via \pkg{L0glm};
#' optional post-processing can log-concave-smooth and/or recentre the PSF.
#'
#' @inheritParams L0_sparse_deco
#' @param y_deconv Matrix or vector with the (already estimated) latent
#'                 spikes.
#' @param L        Length of the PSF to estimate (ignored when
#'                 \code{kernel_est} is supplied).
#' @param kernel_est Optional prior guess for the PSF.  It is
#'                   L1-normalised internally and used as the shrinkage
#'                   centre in \pkg{L0glm}.
#' @param kernel_smoothfun Optional function applied to the raw estimate
#'                         (default: \code{\link{log_conc_smooth}}).
#' @param perc   Passed to \code{\link{shift_kernel}} – use `"mean"` to
#'               recenter the centroid, or a value in \eqn{[0,1]} for a
#'               percentile.
#'
#' @return A list with
#' \describe{
#'   \item{\code{kernel}}{The final (average) PSF estimate.}
#'   \item{\code{kernels}}{Channel-wise estimates before averaging.}
#' }
#'
#' @seealso
#'   \code{\link{L0_sparse_deco}}, \code{\link{shift_kernel}},
#'   \code{\link{log_conc_smooth}}
#'
#' @examples
#' if (requireNamespace("L0glm", quietly = TRUE)) {
#'   ## toy example continuing the previous one
#'   n   <- 200
#'   ks0 <- convkernel("gaussian", 25, list(sigma = 3))  # true PSF
#'   a   <- rpois(n, 0.1)
#'   y   <- rpois(n, as.vector(convmatrix(ks0, n) %*% a))
#'
#'   a_est <- L0_sparse_deco(y, kernel = ks0)$y_deconv
#'   g_est <- L0_deco_estkernel(y, a_est, L = 25)$kernel
#'   matplot(cbind(ks0, g_est), type = "l", lty = 1,
#'           col = c("grey", "blue"), main = "PSF: true vs. estimate")
#'   legend("topright", c("true", "estimated"),
#'          col = c("grey", "blue"), lty = 1, bty = "n")
#' }
#'
#' @importFrom Matrix Matrix Diagonal
#' @importFrom L0glm L0glm.fit
#' @export
L0_deco_estkernel <- function(y,     # vector or matrix (for multivariate case) with response(s) to be deconvolved
                              y_deconv, # vector or matrix (for multivariate case) with estimated deconvolved unblurred signal(s)
                              kernel_est = NULL, # vector with prior estimate of kernel, if available - is always normalised to unit L1 norm before use
                              L = ifelse(!is.null(kernel_est), length(kernel_est), 41), # desired length of  convolution kernel / PSF to be estimated
                              family = poisson(identity),   # noise model, we only support identity link families here
                              kernel_smoothfun = function (y) log_conc_smooth(y), # NULL or function to smooth estimated kernel with (here we use log-concave spline)
                              perc = "mean", # shift PSF to centroid
                              IC = "gic") { # preset L0glm regularisation to approximately optimise this IC, choices aic, bic, gic, mbic or ebic)  
  if (is.vector(y)) y <- as.matrix(y)
  if (is.vector(y_deconv)) y_deconv <- Matrix(y_deconv, sparse=TRUE)
  n <- nrow(y)
  
  kernel_est_scaling <- sum(kernel_est) # scale kernel_est / PSF_est to unit L1 norm (standard for PSF)
  kernel_est <- kernel_est/kernel_est_scaling
  
  # rs <- rowSums(y_deconv) # L1 norm normalise rows of y_deconv so that with Poisson noise estimated kernel would also have expected Poisson noise
  # y_deconv_scalings <- ifelse(rs > 0, 1/rs, 1)                  # 1/row-sum   (rows that are all-zero keep factor 1)
  # y_deconv <- Diagonal(x = y_deconv_scalings) %*% y_deconv      # left-multiply by a sparse diagonal
  
  # L1 norm normalise rows of y_deconv so that with Poisson noise kernel_raw would also have expected Poisson noise
  cs <- colSums(y_deconv)
  y_deconv_scalings <- ifelse(cs > 0, 1/cs, 1)
  y_deconv <- y_deconv %*% Diagonal(x = y_deconv_scalings)
  
  kernels <- do.call(rbind, lapply(1:ncol(y), # loop over outcome channels for multivariate case
                                   function (col) {
                                     X <- convmatrix(kernel = y_deconv[, col],
                                                     n = n,
                                                     shape = "same",
                                                     n_cols = L,
                                                     origin = (L + 1L) %/% 2L)
                                     kernel <- L0glm.fit(y = y[, col, drop = FALSE],  # nonnegative L0-norm penalised GLM, but now shrinking towards shrink towards prior kernel estimate 
                                                         X = X, 
                                                         family = family, 
                                                         intercept = FALSE,
                                                         # no.pen = seq(1, ncol(X), by=2), # not penalising every other scan line could be similar to penalising the difference between consecutive coefficients, i.e. to using a fused ridge penalty, but with pre-set lambda that will approx optimise GIC
                                                         nonnegative = TRUE,
                                                         # start = XXX,
                                                         warm.start = 1,
                                                         # prior.mean = kernel_est*y_deconv_scalings[col], # shrink towards prior kernel estimate (more stable)
                                                         grouped = FALSE, 
                                                         lambda = IC,
                                                         method = 7L, # using sparse Cholesky LLT solver in adaptive ridge IRLS iterations
                                                         solver = "eigen")$coefficients
                                     # optional: smooth estimated psf (NOTE: maybe also consider total variation denoiser tvR::denoise1?)
                                     if (!is.null(kernel_smoothfun)) { kernel <- kernel_smoothfun(y = kernel) }
                                     # center PSF at given percentile p or with "mean" put centroid in middle
                                     kernel <- shift_kernel(kernel, 
                                                            perc = perc, # put centroid in middle 
                                                            method = "linear")
                                     return(kernel/sum(kernel)) # normalise PSF to unit L1 norm
                                   } ))
  # average obtained kernels for different channels - we assume they should be identical, but also return the separate kernels
  kernel <- colMeans(kernels) 
  # plot(kernel, type="l", lwd=5)    
  # matplot(t(kernels), type="l", add = TRUE)    
  
  return(list(kernels = kernels, kernel = kernel))
}