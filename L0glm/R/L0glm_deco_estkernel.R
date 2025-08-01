#' @name L0glm_deco_estkernel
#' @title Estimate the PSF given a spike train estimate via L0‐penalised GLMs
#'
#' @description Solves the \emph{transpose} problem: given an estimated latent
#' single- or multichannel spike train signal \eqn{a} and the vector (or matrix)
#' of observations \eqn{y}, find the PSF \eqn{g} that minimises the deviance in
#' \eqn{y \approx X(g)\,a}, where \eqn{X(g)} is a Toeplitz design matrix that
#' contains time-shifted copies of point-spread function \eqn{g} (PSF, a.k.a.
#' convolution kernel) and errors \eqn{\varepsilon} are assumed to be drawn from
#' an exponential family error distribution (e.g. Poisson, Gaussian, Gamma,
#' Inverse Gaussian). Regularisation is again done via \pkg{L0glm}. An optional
#' prior guess for the PSF can be given, in which case it is used as the
#' shrinkage centre in \pkg{L0glm}. Optional post-processing applies a
#' log-concave spline smooth or custom smoothing function and/or recentres the
#' PSF.
#'
#' @details The algorithm is inspired by the L0-penalised blind deconvolution
#' approach of de Rooi & Eilers (2011) for chromatographic pulse trains, but it
#' leverages the general machinery of \pkg{L0glm}.  This results in the
#' following improvements:
#' \itemize{
#'   \item No manual tuning of the regularisation parameter is required; instead
#'         it can be preset to optimise information criteria such as GIC, AIC,
#'         BIC, EBIC or MBIC
#'   \item Both single channel or multichannel data are allowed for
#'    \item We allow the estimated convolution kernel \eqn{g} to be smoothed by
#'          a custom function (e.g. log-concave splines)
#'   \item Any exponential-family noise model with an \emph{identity} link
#'         can be used (Poisson, Gaussian, Gamma, …)
#'   \item Shrinkage towards an optional prior guess for the PSF
#'         is allowed for, resulting in enhanced stability and improved convergence
#'         behaviour
#'   \item The underlying sparse design is handled by Eigen's efficient sparse
#'         Cholesky solver, enabling large problems to be tackled quickly
#' }
#'
#' @inheritParams L0glm_sparse_deco
#' @param y_deconv Matrix or vector with the (already estimated) latent spikes.
#' @param L        Length of the PSF to estimate (ignored when \code{kernel_est}
#'   is supplied).
#' @param kernel_est Optional prior guess for the PSF.  It is L1-normalised
#'   internally and used as the shrinkage centre in \pkg{L0glm}.
#' @param kernel_smoothfun Optional function applied to the raw estimate
#'   (default: \code{\link{log_conc_smooth}}).
#' @param perc   Passed to \code{\link{shift_kernel}} – use `"mean"` to recenter
#'   the centroid, or a value in \eqn{[0,1]} for a percentile.
#'
#' @return A list with
#' \describe{
#'   \item{\code{kernel}}{The final (average) PSF estimate.}
#'   \item{\code{kernels}}{Channel-wise estimates before averaging.}
#' }
#'
#' @references de Rooi, J. J. & Eilers, P. H. C. (2011). \emph{Deconvolution of
#' pulse trains with the L0 penalty}. \strong{Analytica Chimica Acta},
#' \emph{705}(1-2), 218-226. DOI: 10.1016/j.aca.2011.07.012
#'
#' @seealso \code{\link{L0glm_sparse_deco}}, \code{\link{shift_kernel}},
#' \code{\link{log_conc_smooth}}
#'
#' @examples
#' if (requireNamespace("L0glm", quietly = TRUE)) {
#'
#'   ## ------------------------------------------------------------------------
#'   ## 1  Simulate a blurred spike-train with a BIGAUSSIAN PSF & Poisson noise
#'   ## ------------------------------------------------------------------------
#'   set.seed(123)
#'   library(Matrix)
#'
#'   ## --- true point–spread function (PSF) -----------------------------
#'   L       <- 51                       # kernel length
#'   sigma1  <- 1.25                     # left-hand σ
#'   sigma2  <- 5                        # right-hand σ
#'   g_true  <- convkernel("bigaussian", L,
#'                         params = list(sigma1 = sigma1,
#'                                       sigma2 = sigma2),
#'                         perc   = 0)   # already centred
#'
#'   ## --- latent sparse spike-train ------------------------------------
#'   n        <- 10000                   # timepoints
#'   npeaks   <- round(n / 10)           # 1000, ~10 % non-zero
#'   idx      <- sample(n, npeaks)
#'   amps     <- sample(1000:100000, npeaks)
#'   a_true   <- numeric(n);  a_true[idx] <- amps   # true spikes
#'
#'   ## --- generate blurred, noisy observations -------------------------
#'   X_sim  <- convmatrix(g_true, n)
#'   y_true <- as.vector(X_sim %*% a_true)     # noiseless convolution
#'   y      <- rpois(n, y_true)                # with Poisson noise
#'
#'   ## ------------------------------------------------------------------------
#'   ## 2  Deconvolution with an initial PSF guess
#'   ## ------------------------------------------------------------------------
#'
#'   g0 <- convkernel("gaussian", L, list(sigma = 6))  # initial guess for PSF
#'   fit <- L0glm_sparse_deco(y, kernel = g0,
#'                            family = poisson(identity), IC = "gic")
#'   a_est <- as.vector(fit$y_deconv)                  # estimated spike train
#'
#'   ## ------------------------------------------------------------------------
#'   ## 3  Update PSF estimate given initial spike train estimate
#'   ## ------------------------------------------------------------------------
#'
#'   fit <- L0glm_deco_estkernel(y,
#'                               y_deconv = a_est,  # estimated spike train
#'                               kernel_est = g0,   # for more stable & faster convergence
#'                               L = L,
#'                               family = poisson(identity),
#'                               IC = "gic")
#'   g_est <- fit$kernel
#'
#'   plot(g_true, col="red",
#'        lwd=7, type="l",
#'        ylab = "Density",
#'        main = "True (red) vs. estimated (blue) bigaussian PSF\nafter 1 iteration of Poisson L0 glm blind deconvolution")
#'   lines(g_est, col="blue", lwd=3)
#' }
#'
#' @importFrom Matrix Matrix Diagonal colSums
#' @importFrom L0glm L0glm.fit
#' @export
L0glm_deco_estkernel <- function(y,     # vector or matrix (for multivariate case) with response(s) to be deconvolved
                              y_deconv, # vector or matrix (for multivariate case) with estimated deconvolved unblurred signal(s)
                              kernel_est = NULL, # vector with prior estimate of kernel, if available - is always normalised to unit L1 norm before use
                              L = ifelse(!is.null(kernel_est), length(kernel_est), 51), # desired length of  convolution kernel / PSF to be estimated
                              family = poisson(identity),   # noise model, we only support identity link families here
                              kernel_smoothfun = function (y) log_conc_smooth(y), # NULL or function to smooth estimated kernel with (here we use log-concave spline)
                              perc = "mean", # shift PSF to centroid
                              IC = "gic") { # preset L0glm regularisation to approximately optimise this IC, choices aic, bic, gic, mbic or ebic)  
  if (is.vector(y)) y <- as.matrix(y)
  if (is.vector(y_deconv)) y_deconv <- Matrix(y_deconv, sparse=TRUE)
  n <- nrow(y)
  
  if (!is.null(kernel_est)) {
    kernel_est <- kernel_est/sum(kernel_est) # scale kernel_est / PSF_est to unit L1 norm (standard for PSF)
  } else { prior.mean <- NULL }
    
  # L1 norm normalise each column of y_deconv so that with Poisson noise estimated kernel would have expected Poisson noise
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
                                     # shrink towards prior kernel estimate if kernel_est given (more stable)
                                     if (!is.null(kernel_est)) { prior.mean = kernel_est/y_deconv_scalings[col] } 
                                     kernel <- L0glm.fit(y = y[, col, drop = FALSE],  # nonnegative L0-norm penalised GLM, but now shrinking towards shrink towards prior kernel estimate 
                                                         X = X, 
                                                         family = family, 
                                                         intercept = FALSE,
                                                         # no.pen = seq(1, ncol(X), by=2), # not penalising every other scan line could be similar to penalising the difference between consecutive coefficients, i.e. to using a fused ridge penalty, but with pre-set lambda that will approx optimise GIC
                                                         nonnegative = TRUE,
                                                         start = prior.mean,
                                                         warm.start = 1,
                                                         prior.mean = prior.mean, 
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