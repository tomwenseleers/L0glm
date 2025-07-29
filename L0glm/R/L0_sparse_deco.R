#' @name L0_sparse_deco
#' @title Sparse (L0-penalised) deconvolution given a fixed PSF
#'
#' @description
#' Estimates the *latent*, un-blurred spike train \eqn{a} from an observed
#' signal \eqn{y} that is assumed to follow the convolution model
#' \deqn{y \;=\; X(\theta)\,a + \varepsilon,}
#' where \eqn{X(\theta)} is a Toeplitz design matrix that contains
#' time-shifted copies of a known point-spread function (PSF, a.k.a. kernel).
#' The coefficient vector (or matrix) \eqn{a} is found with **best-subset
#' (L0)** regularisation via \pkg{L0glm}.
#'
#' @param y       Numeric vector or matrix; the observed (blurred) signal.
#'                A matrix is treated as \eqn{n \times p} parallel channels.
#' @param family  A GLM \link[stats]{family} object with an *identity*
#'                link (e.g. `poisson(identity)`, `gaussian(identity)`
#'                or `Gamma(identity)`).
#' @param kernel  Either
#'   * a numeric vector that is the fixed PSF **in forward time order**, or
#'   * `NULL` together with `L` and `kernelfunc`, in which case the PSF is
#'     generated on the fly (useful for simulations).
#' @param kernelfunc A function that returns a length-\eqn{L} numeric vector
#'                   when called as `kernelfunc(1:L)`. Default a Gaussian with 
#'                   sigma=6.
#' @param L       Desired kernel length when `kernel` is `NULL`.
#' @param X       **Optional.**  A sparse design matrix that overrides
#'                the automatic call to \code{\link{convmatrix}}.  This is
#'                useful when you need time-varying PSFs.
#' @param IC      Information criterion that sets the penalty parameter in
#'                \pkg{L0glm}. One of `"gic"` (default), `"aic"`, `"bic"`,
#'                `"mbic"` or `"ebic"`.
#' @param grouped Logical, passed unchanged to \code{L0glm::L0glm.fit}, specifying
#'                whether to perform grouped L0 feature selection, i.e. best subset
#'                of group selection (BSGS).
#'
#' @return A list with
#' \describe{
#'   \item{\code{y_deconv}}{Estimated un-blurred spikes (same shape as \code{y}).}
#'   \item{\code{fitted.values}}{\eqn{X a}\, – the fitted (re-blurred) signal.}
#'   \item{\code{deviance}}{Residual deviance of the fit.}
#' }
#'
#' @seealso
#'   \code{\link{L0_deco_estkernel}} for the complementary
#'   “estimate-the-kernel” step,  
#'   \code{\link{L0_blind_deco}} for the full alternating blind
#'   deconvolution,  
#'   \code{\link{convmatrix}} (design matrix) and
#'   \code{\link{convkernel}} (PSF generator).
#'
#' @examples
#' if (requireNamespace("L0glm", quietly = TRUE)) {
#'   set.seed(1)
#'   n  <- 200
#'   psf <- convkernel("gaussian", 25, list(sigma = 4)) # PSF
#'   X  <- convmatrix(psf, n)                # design matrix
#'   a  <- rpois(n, 0.1)                     # sparse ground-truth spikes
#'   y  <- as.vector(X %*% a)                # noiseless convolution
#'   y  <- rpois(n, y)                       # add Poisson noise
#'
#'   fit <- L0_sparse_deco(y, kernel = psf, family = poisson(identity))
#'   a_est <- as.vector(fit$y_deconv)        # estimated spikes
#'   plot(a, type = "h", lty = 1,
#'        col = "red", ylim = c(-max(a, a_est), max(a, a_est)),
#'        xlab = "Time index", ylab = "Coefficients",
#'        main = "True (red) vs. estimated (blue) spikes")
#'   lines(-a_est, type = "h", lty = 1, col = "blue")
#' }
#'
#' @importFrom Matrix Matrix Diagonal
#' @importFrom L0glm L0glm.fit
#' @export
L0_sparse_deco <- function(y,                            # response to be deconvolved (vector or matrix, i.e. univariate or multivariate)
                           family = poisson(identity),   # noise model, e.g. poisson/gaussian/Gamma, we only support identity link families here
                           kernel = NULL, # kernel: either vector with blur kernel or length of blur kernel L & kernelfunc
                           kernelfunc  = function(x)
                             convkernel("gaussian",
                                        L        = length(x),   # <- same length
                                        params   = list(sigma = 6),
                                        perc     = 0),          # no recentring
                           # kernelfunc = function(x) gauspeak(x, u = (1+L)/2, w = 6, h = 1),
                           L = 41,
                           X = NULL, # (sparse) deconvolution kernel - can be used to specify matrix with time-varying PSFs & which would then override kernel
                           IC = "gic", # preset L0glm regularisation to approximately optimise this IC, choices aic, bic, gic, mbic or ebic 
                           grouped = TRUE) {  # use grouped feature selection (best subset of group selection), only used if outcome y is multivariate
  dev_res <- family$dev.resids # deviance residuals function
  if (is.vector(y)) y <- as.matrix(y)
  n <- nrow(y)
  p <- ncol(y)
  
  if (is.null(kernel)) {
    kernel <- kernelfunc(1L:L)
    L      <- length(kernel)
  }  
  
  if (is.null(X)) { 
    X <- convmatrix(kernel = kernel, n = n)  # covariate matrix with time-shifted versions of kernel / psf 
  } else { X <- Matrix(X, sparse = TRUE) }
  y_deconv <- L0glm.fit(y = y,   # nonnegative L0-norm penalised GLM
                        X = X,  
                        family = family, 
                        intercept = FALSE,
                        # start = XXX,
                        warm.start = 1,
                        nonnegative = TRUE, # imposing nonnegativity constraints on coefficients
                        lambda = IC, 
                        grouped = grouped,
                        method = 7L, # using sparse Cholesky LLT solver
                        solver = "eigen")$coefficients
  y_hat <- X %*% y_deconv # fitted values # cbind(1, )
  deviance <- sum(dev_res(y = as.numeric(y), 
                          mu = as.numeric(y_hat+1E-10), # to avoid Inf values with poisson
                          wt = rep(1, prod(dim(y)))))
  # if (is.vector(y_deconv)) { y_deconv <- y_deconv[-1] } else { y_deconv <- y_deconv[-1,,drop=FALSE] } # remove intercept again
  return(list(y_deconv = y_deconv, fitted.values = y_hat, deviance = deviance))
}
