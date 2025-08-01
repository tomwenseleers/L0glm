#' @name L0glm_sparse_deco
#' @title Sparse deconvolution with a fixed PSF via L0‐penalised GLMs
#'
#' @description
#' Estimates the \emph{latent}, un-blurred single or multichannel spike train 
#' \eqn{a} from an observed signal \eqn{y} that is assumed to follow the 
#' convolution model
#' \deqn{y \;=\; X(g)\,a + \varepsilon,}
#' where \eqn{X(g)} is a Toeplitz design matrix that contains
#' time-shifted copies of a known point-spread function \eqn{g} (PSF, a.k.a. convolution 
#' kernel) and errors \eqn{\varepsilon} are assumed to be drawn from an exponential 
#' family error distribution (e.g. Poisson, Gaussian, Gamma, Inverse Gaussian). 
#' The coefficient vector (or matrix) \eqn{a} is found with \emph{best-subset} 
#' or \emph{best-subset of group selection} via \pkg{L0glm}.
#' 
#' @details
#' The algorithm is inspired by the L0-penalised deconvolution approach
#' of de Rooi & Eilers (2011) for chromatographic pulse trains, but it
#' leverages the general machinery of \pkg{L0glm}.  This results in the following
#' improvements:
#' \itemize{
#'   \item No manual tuning of the regularisation parameter is required; instead
#'         it can be preset to optimise information criteria such as GIC, AIC, 
#'         BIC, EBIC or MBIC
#'   \item Both single channel data with L0 best subset selection or 
#'         multichannel data with grouped L0 best subset of group selection 
#'         are allowed for         
#'   \item Any exponential-family noise model with an \emph{identity} link
#'         can be used (Poisson, Gaussian, Gamma, …)
#'   \item The underlying sparse design is handled by Eigen's efficient sparse
#'         Cholesky solver, enabling large problems to be tackled quickly
#' }
#' 
#' In the non-negative-matrix-factorisation (NMF) literature this model
#' would be referred to as \emph{single-atom convolutional NMF} or \emph{sparse 
#' convolutional dictionary learning}, whereby a one–dimensional signal (or
#' multichannel matrix) \eqn{y} is approximated as \eqn{W(g)\,H}, where \eqn{W(g)}
#' is a banded Toeplitz matrix that holds time–shifted copies of
#' a fixed dictionary atom/PSF \eqn{g}, and \eqn{H} is a
#' (potentially multi-column) matrix of sparse, non-negative activations
#' \eqn{a}.
#' However, relative to standard single-atom convolutional NMF we 
#' offer several extensions:
#' \itemize{
#'   \item \emph{Optimal sparsity}: we do not fix the sparsity level of the 
#'       activations \eqn{a} in advance; instead we let \pkg{L0glm} perform
#'       best-subset or best subset of group selection through an
#'       L\eqn{_{0}} or grouped L\eqn{_{0}} penalty,
#'   \item \emph{General noise models}: any exponential-family likelihood
#'         with an \emph{identity} link is supported (Poisson, Gaussian,
#'         Gamma, inverse-Gaussian, ...), instead of only Gaussian
#'         (Euclidean-distance NMF) or Poisson (NMF with KL divergence)
#'   \item \emph{Speed}: the Toeplitz design is solved by sparse-Cholesky 
#'         linear algebra rather than multiplicative updates, giving 
#'         orders-of-magnitude speed-ups on long signals
#' }
#' Methodologically this could be referred to as
#' \emph{generalised / deviance matrix factorisation}
#' (Kidzinski et al. 2022; Castiglione et al. 2024; Wang & Carvalho 2023), 
#' specialised to the case where \eqn{W} is a Toeplitz single-atom convolutional 
#' dictionary matrix and where the matrix factorisation is solved with an 
#' L\eqn{_{0}} or grouped L\eqn{_{0}}-penalised generalised linear model (GLM)
#' fit via adaptive ridge regularised iteratively reweighted least squares 
#' (IRLS).
#'
#' A classical baseline for Poisson deconvolution is the
#' Richardson–Lucy algorithm (Richardson 1972; Lucy 1974). In modern
#' terms, Richardson–Lucy equals Kullback–Leibler NMF with
#' multiplicative updates (Lee & Seung 1999) and  
#' factor matrix \eqn{W} constrained to be a convolution matrix.  
#' Empirically, \code{L0glm_sparse_deco()} has much better support recovery,
#' attains much lower reconstruction error and is dramatically faster thanks to 
#' (i) closed-form sparse linear solves, (ii) the automatic selection of the penalty
#' parameter through information-criterion presets, and (iii) the possibility
#' to exploit grouped-sparsity when several channels are deconvolved jointly.
#'
#' @param y       Numeric vector or matrix; the observed (blurred) signal.
#'                A matrix is treated as \eqn{n \times m} outcome channels.
#' @param family  A GLM \link[stats]{family} object with an \emph{identity}
#'                link specifying the desired error distribution 
#'                (e.g. `poisson(identity)`, `gaussian(identity)`
#'                or `Gamma(identity)`).
#' @param kernel  Either a numeric vector that is the fixed PSF in forward time 
#'                order, or `NULL` together with `L` and `kernelfunc`, in which 
#'                case the PSF is generated on the fly (useful for simulations).
#' @param kernelfunc A function that returns a length-\eqn{L} numeric vector
#'                   when called as `kernelfunc(1:L)`. Defaults to a Gaussian with 
#'                   sigma=6.
#' @param L       Desired kernel length when `kernel` is `NULL`.
#' @param X       Optional.  A sparse design matrix that overrides
#'                the automatic call to \code{\link{convmatrix}}.  This is
#'                useful when you need time-varying PSFs.
#' @param IC      Information criterion that sets the penalty parameter in
#'                \pkg{L0glm}. One of `"gic"` (default), `"aic"`, `"bic"`,
#'                `"mbic"` or `"ebic"`.
#' @param grouped Logical, passed unchanged to \code{L0glm.fit}, specifying
#'                whether to perform grouped L0 feature selection, i.e. best subset
#'                of group selection (BSGS), in case there are multiple outcome
#'                channels. Defaults to TRUE.
#'
#' @return A list with
#' \describe{
#'   \item{\code{y_deconv}}{Estimated un-blurred spikes (same shape as \code{y}).}
#'   \item{\code{fitted.values}}{\eqn{X a}, the fitted (re-blurred) signal.}
#'   \item{\code{deviance}}{Residual deviance of the fit.}
#' }
#'
#' @references
#'   \itemize{
#'     \item Richardson, W. H. (1972)\cr
#'           Bayesian‐based iterative method of image restoration.\cr
#'           \emph{J. Opt. Soc. Am.} 62 (1): 55–59.
#'     \item Lucy, L. B. (1974)\cr
#'           An iterative technique for the rectification of observed
#'           distributions.\cr
#'           \emph{Astron. J.} 79: 745–754.
#'     \item Smaragdis, P. (2004)\cr
#'           Non‐negative matrix factor
#'           deconvolution; extraction of multiple sound sources from
#'           monophonic inputs.\cr
#'           \emph{Proceedings of ICA 2004}.
#'     \item Lee, D. D., & Seung, H. S. (1999)\cr
#'           Learning the parts of objects by non‐negative matrix
#'           factorization.\cr
#'           \emph{Nature} 401 (6755): 788–791.
#'     \item de Rooi, J., & Eilers, P. H. C. (2011)\cr
#'           Deconvolution of pulse trains with the L0 penalty.\cr
#'           \emph{Analytica Chimica Acta} 705 (1-2): 218–226.
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
#' @seealso
#'   \code{\link{L0glm_deco_estkernel}} for the complementary
#'   “estimate-the-kernel” step,  
#'   \code{\link{L0glm_blind_deco}} for the full alternating blind
#'   deconvolution,  
#'   \code{\link{convmatrix}} (design matrix) and
#'   \code{\link{convkernel}} (PSF generator).
#'
#' @examples
#' if (requireNamespace("L0glm", quietly = TRUE)) {
#'
#'   ## ------------------------------------------------------------------------
#'   ## 1  Simulate a blurred spike-train with a BIGAUSSIAN PSF & Poisson noise
#'   ## ------------------------------------------------------------------------
#'   set.seed(123)
#'
#'   ## --- true point–spread function (PSF) -----------------------------
#'   L       <- 51                       # kernel length
#'   sigma1  <- 1.25                     # left-hand σ
#'   sigma2  <- 5                        # right-hand σ
#'   g_true  <- convkernel("bigaussian", L,
#'                         params = list(sigma1 = sigma1,
#'                                        sigma2 = sigma2),
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
#'   ## 2  Deconvolution with the known PSF
#'   ## ------------------------------------------------------------------------
#'   fit <- L0glm_sparse_deco(y, kernel = g_true,
#'                            family = poisson(identity), IC = "gic")
#'   a_est <- as.vector(fit$y_deconv)
#'
#'   ## ---------- summary statistics ------------------------------------
#'   n_true <- sum(a_true != 0)
#'   n_est  <- sum(a_est  != 0)
#'   tp     <- sum(a_true != 0 & a_est != 0)
#'   precision <- tp / n_est
#'   recall    <- tp / n_true
#'   F1        <- 2 * precision * recall / (precision + recall)
#'   R2p <- 100 *
#'     summary(lm(a_true ~ 0 + a_est))$r.squared
#'
#'   ## ---------- scatter-plot with true vs. estimated peak areas -------
#'   main_expr <- bquote(
#'     atop("Poisson L0-glm deconvolution with known PSF",
#'          atop("true vs. estimated peak areas",
#'               italic(R)^2 == .(sprintf("%.2f", R2p)) * "%" ~ "|" ~
#'               "peaks:" ~ .(n_true) ~ "(true)," ~ .(n_est) ~ "(est.)" ~ "|" ~
#'               "F1 =" ~ .(sprintf("%.2f", 100 * F1)) * "%")))
#'
#'   plot(a_est, a_true,
#'        xlim = range(a_est, a_true),
#'        ylim = range(a_est, a_true),
#'        pch  = 16, col = "blue",
#'        xlab = "estimated value",
#'        ylab = "true value",
#'        main = main_expr)
#' }
#'
#' @importFrom Matrix Matrix Diagonal
#' @importFrom L0glm L0glm.fit
#' @export
L0glm_sparse_deco <- function(y,                         # response to be deconvolved (vector or matrix, i.e. univariate or multivariate)
                           family = poisson(identity),   # noise model, e.g. poisson/gaussian/Gamma, we only support identity link families here
                           kernel = NULL, # kernel: either vector with blur kernel or length of blur kernel L & kernelfunc
                           kernelfunc  = function(x)
                             convkernel("gaussian",
                                        L        = length(x),   # <- same length
                                        params   = list(sigma = 6),
                                        perc     = 0),          # no recentering as this PSF is already centered
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
                        warm.start = 1,
                        nonnegative = TRUE, # imposing nonnegativity constraints on coefficients
                        lambda = IC, 
                        grouped = grouped,
                        method = 7L, # using sparse Cholesky LLT solver
                        solver = "eigen")$coefficients
  y_hat <- X %*% y_deconv # fitted values
  deviance <- sum(dev_res(y = as.numeric(y), 
                          mu = as.numeric(y_hat+1E-10), # to avoid Inf values with poisson
                          wt = rep(1, prod(dim(y)))))
  return(list(y_deconv = y_deconv, fitted.values = y_hat, deviance = deviance))
}
