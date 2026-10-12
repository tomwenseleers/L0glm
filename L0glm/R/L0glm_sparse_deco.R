#' @name L0glm_sparse_deco
#' @title Sparse deconvolution with a fixed PSF via L0-penalised GLMs
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
#'         it can be preset to approximate information criteria such as GIC, MEBIC, AIC,
#'         BIC, EBIC or MBIC
#'   \item Both single channel data with L0 best subset selection or 
#'         multichannel data with grouped L0 best subset of group selection 
#'         are allowed for         
#'   \item Any exponential-family noise model with an \emph{identity} link
#'         can be used (Poisson, Gaussian, Gamma, ...)
#'   \item The underlying sparse design is handled by Eigen's efficient sparse
#'         Cholesky solver, enabling large problems to be tackled quickly
#' }
#' 
#' In the non-negative-matrix-factorisation (NMF) literature this model
#' would be referred to as \emph{single-atom convolutional NMF} or \emph{sparse 
#' convolutional dictionary learning}, whereby a one-dimensional signal (or
#' multichannel matrix) \eqn{y} is approximated as \eqn{W(g)\,H}, where \eqn{W(g)}
#' is a banded Toeplitz matrix that holds time-shifted copies of
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
#' Richardson-Lucy algorithm (Richardson 1972; Lucy 1974). In modern
#' terms, Richardson-Lucy equals Kullback-Leibler NMF with
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
#' @param baseline  Character.  One of
#'                  \code{NULL} (no baseline, default),
#'                  \code{"intercept"} (a single constant column),
#'                  or \code{"bspline"} (cubic B-spline basis).
#'                  The baseline columns are prepended to the Toeplitz
#'                  design and are \emph{always included} in the model.
#'                  Their fitted contribution is used as an \code{offset}
#'                  in the \code{L0glm_deco_estkernel} kernel-update step.
#' @param nknots    Integer >= 4.  Number of spline knots
#'                  to use when \code{baseline = "bspline"}. With \code{nknots < 4}
#'                  and \code{baseline="bspline"} \code{baseline} is changed
#'                  to \code{intercept}. 
#' @inheritParams L0glm.tune
#' @param grouped Logical, passed unchanged to \code{L0glm.fit}, specifying
#'                whether to perform grouped L0 feature selection, i.e. best subset
#'                of group selection (BSGS), in case there are multiple outcome
#'                channels. Defaults to TRUE.
#' @param control.solver Optional solver-control list passed to \code{L0glm.fit}.
#'                By default this uses Eigen sparse Cholesky LLT with AMD
#'                ordering. For strictly banded PSF designs, users can pass
#'                \code{control.solver.gen("eigen", method=1L, NaturalOrdering=TRUE)}
#'                to request Eigen natural ordering.
#'
#' @return A list with
#' \describe{
#'   \item{\code{coef_spikes}}{Estimated un-blurred spikes, i.e. the coefficients
#'                             for design matrix \code{X_spikes} (same shape as \code{y}).}
#'   \item{\code{coef_baseline}}{Estimated baseline coefficients for baseline
#'                               design matrix \code{X_baseline}.}
#'   \item{\code{fitted_values}}{\eqn{X a}, the fitted signal.}
#'   \item{\code{fitted_baseline}}{\eqn{X a}, the fitted baseline}
#'   \item{\code{deviance}}{Residual deviance of the fit.}
#'   \item{\code{X_spikes}}{Toeplitz design matrix for the spikes.}
#'   \item{\code{X_baseline}}{Design matrix for the baseline.}
#' }
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
#' @seealso
#'   \code{\link{L0glm_deco_estkernel}} for the complementary
#'   "estimate-the-kernel" step,  
#'   \code{\link{L0glm_blind_deco}} for the full alternating blind
#'   deconvolution,  
#'   \code{\link{convmatrix}} (design matrix) and
#'   \code{\link{convkernel}} (PSF generator).
#'
#' @example examples/L0glm_sparse_deco.R
#' @importFrom Matrix Matrix Diagonal
#' @importFrom splines bs
#' @export
L0glm_sparse_deco <- function(y,                         # response to be deconvolved (vector or matrix, i.e. univariate or multivariate)
                           family = poisson(identity),   # noise model, e.g. poisson/gaussian/Gamma, we only support identity link families here
                           kernel = NULL, # kernel: either vector with blur kernel or length of blur kernel L & kernelfunc
                           kernelfunc  = function(x)
                             convkernel("gaussian",
                                        L        = length(x),   # <- same length
                                        params   = list(sigma = 6),
                                        perc     = 0),          # no recentering as this PSF is already centered
                           L = if (!is.null(kernel)) length(kernel) else 51,
                           X = NULL, # (sparse) deconvolution kernel - can be used to specify matrix with time-varying PSFs & which would then override kernel
                           baseline    = NULL,
                           nknots      = max(4L, floor(nrow(as.matrix(y))/500)),
                           IC = "gic", # preset L0glm regularisation to approximately optimise this IC, choices aic, bic, gic, mebic, ebic or mbic
                           grouped = TRUE,  # use grouped feature selection (best subset of group selection), only used if outcome y is multivariate
                           control.solver = do.call("control.solver.gen",
                                                    list(solver = "eigen", method = 1L))) {
  dev_res <- family$dev.resids # deviance residuals function
  if (is.vector(y)) y <- as.matrix(y)
  n <- nrow(y)
  
  ## ---------- PSF design ----------
  if (is.null(kernel)) kernel <- kernelfunc(1L:L)
  if (is.null(X)) X <- convmatrix(kernel = kernel, n = n)  # covariate matrix with time-shifted versions of kernel / psf 

  ## ---------- baseline design ----------
  if (!is.null(baseline) & nknots < 4) baseline <- "intercept"
  if (!is.null(baseline)) {
    if (baseline == "intercept") {
      X_base <- Matrix::Matrix(1, n, 1, sparse = TRUE)
    } else if (baseline == "bspline") {
      X_base <- splines::bs(seq_len(n), df = nknots, degree = 3, intercept = TRUE)
      X_base <- Matrix::Matrix(X_base, sparse = TRUE)
  }} else {
    X_base <- NULL
  }
  
  ## ---------- final model matrix ----------
  X_full <- if (is.null(X_base)) X else cbind(X_base, X)
  ## ---------- L0glm fit -------------------
  fit <- L0glm.fit(y = y,   # nonnegative L0-norm penalised GLM
                   X = X_full, 
                   family = family, 
                   intercept = FALSE,
                   start = rep(1, ncol(X_full)), # WORKED WHILE start=NULL DID NOT
                   warm.start = 1,
                   nonnegative = TRUE, # nonnegative coefficients, including any intercept
                   lambda = IC, 
                   grouped = grouped,
                   prior.mean = if (!is.null(X_base)) c(rep(min(y), ncol(X_base)), rep(0,ncol(X))) else NULL, # NEEDED FOR STABLE CONVERGENCE WHEN THERE IS A BASELINE
                   # using a quantile regression underapproximation with tau=0.05 (=sum(y<a_baseline)/n) 
                   # instead of rep(min(y), ncol(X_base))
                   # would also be an option (likely better solution)
                   always.include = if (!is.null(X_base)) 1:ncol(X_base) else NULL,
                   normalize = !is.null(X_base), # required
                   maxit = 100,
                   method = 1L, # sparse Cholesky LLT
                   solver = "eigen",
                   control.solver = control.solver,
                   # verbose = TRUE,
                   gamma = 2,
                   thresh = 1e-10,
                   tol = 1E-4)
  # NOTE: we include X_baseline in design matrix & include those terms in always.include
  # (those terms get their prior.mean updated in each IRLS iteration as opposed
  # to having their estimates shrunk to zero).
  # + use start=rep(1, ncol(X_baseline)+ncol(X_sim)) & normalize = TRUE
  # & prior.mean = c(rep(min(y), ncol(X_baseline)), rep(0,ncol(X_sim)))
  # to shrink baseline coefficients towards expected baseline.
  # This works much better than presubtracting baseline, e.g. using
  # X = X_sim & y = pmax(y-quantile(y, 0.05), 0).
  
  beta      <- fit$coefficients
  
  k_base    <- if (is.null(X_base)) 0 else ncol(X_base)
  beta_base <- if (k_base) beta[seq_len(k_base), , drop = FALSE] else NULL
  beta_spk  <- if (k_base!=0) beta[-seq_len(k_base), , drop = FALSE] else beta
  
  y_hat  <- X_full %*% beta
  y_base <- if (!is.null(X_base)) X_base %*% beta_base else Matrix(0, ncol = ncol(y_hat), nrow = nrow(y_hat), sparse = TRUE)
  dev    <- sum(dev_res(y = as.vector(y), 
                        mu = pmax(as.vector(y_hat), .Machine$double.eps), # to avoid Inf values with poisson
                        wt = 1))
  
  return(list(coef_spikes     = beta_spk,
              coef_baseline   = beta_base,
              fitted_values   = y_hat,
              fitted_baseline = y_base,
              deviance        = dev,
              X_spikes        = X,
              X_baseline      = X_base))
}
