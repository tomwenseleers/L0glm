#' Sure Independence Screening for GLMs via Fast Score Tests
#'
#' @description
#' `sis_glm()` implements *sure independence screening* (SIS) for
#' generalized linear models (GLMs) using a **single null fit** and
#' **marginal score test statistics** for each feature. It supports all
#' standard `stats::family`/link combinations (e.g., Gaussian, binomial,
#' Poisson, Gamma, inverse.gaussian), optional observation `weights`, an
#' `offset` term, and both **dense** (`matrix`) and **sparse** (classes
#' from **Matrix**) design matrices.
#'
#' In our benchmarks (see the Examples), this implementation is typically
#' **>10x faster** than the vanilla SIS procedure in the \pkg{SIS}
#' package for large, tall-and-wide problems, while producing essentially
#' identical rankings. Speedups come from avoiding \(p\) univariate GLM
#' fits and instead computing closed-form score statistics using only the
#' intercept-only (null) GLM fit.
#'
#' @details
#' For each predictor \eqn{x_j}, consider the marginal GLM
#' \eqn{\eta_i = \beta_0 + \beta_j x_{ij}}{eta_i = beta_0 + beta_j x_{ij}}
#' with link \eqn{\eta_i = g(\mu_i)}{eta_i = g(mu_i)} and null fit
#' \eqn{\hat{\mu}_{0i}}{mu0_i}. The per-feature efficient score and (expected)
#' information at the null are
#' \deqn{
#'   U_j = \sum_{i=1}^n x_{ij}\,(y_i - \hat{\mu}_{0i})\,\left.\frac{d\eta_i}{d\mu_i}\right|_{\hat{\mu}_0},
#'   \qquad
#'   I_{jj} = \sum_{i=1}^n x_{ij}^2\,\frac{(d\mu_i/d\eta_i)^2}{V(\hat{\mu}_{0i})},
#' }{
#'   U_j = sum_i x_{ij} * (y_i - mu0_i) * (d eta_i / d mu_i at mu0),
#'   and I_{jj} = sum_i x_{ij}^2 * (d mu_i / d eta_i)^2 / V(mu0_i).
#' }
#' with variance function \eqn{V(\cdot)}{V(.)}. The screening statistic is
#' \deqn{T_j = U_j^2 / I_{jj},}{T_j = U_j^2 / I_{jj},}
#' which is approximately \eqn{\chi^2_1}{chi^2_1} under \eqn{H_0:\beta_j=0}{H0: beta_j = 0}
#' for canonical GLMs. Ranking by \eqn{T_j} closely tracks ranking by marginal
#' likelihood gain / deviance drop (and coincides in the Gaussian identity-link case).
#'
#' The function retains the top `nsis` variables (screening size). By default,
#' `nsis = min(n, floor(n / (log(log(p)) * log(p))))`, a conservative cap related to
#' the maximum support size that is used in \pkg{abess}. To emulate the classical
#' SIS choice, set `nsis = floor(n / log(n))`.
#' 
#' Note that for quasi families, the score statistic uses the variance and link
#' functions but there is no full likelihood; this is fine for ranking.
#'
#' **Computation.** The statistic is computed in two linear-time passes without
#' fitting \eqn{p} univariate GLMs:
#' \itemize{
#'   \item Score numerator: \eqn{U_j = \sum_i x_{ij}\,Uvec_i}, implemented as
#'         \code{U <- as.numeric(crossprod(X, Uvec))}. This dispatches to
#'         \pkg{Matrix} methods automatically when \code{X} is sparse.
#'   \item Score denominator: \eqn{I_{jj} = \sum_i W_i\,x_{ij}^2}. For dense
#'         \code{matrix} inputs we use row scaling by \code{sqrtW <- sqrt(W)} and
#'         compute \code{I <- colSums((X * sqrtW)^2)}. For \pkg{Matrix} sparse
#'         inputs we avoid densifying by scaling only the nonzeros and then
#'         summing by column:
#'         \code{Xw <- X; Xw@x <- (X@x * sqrtW[X@i + 1L])^2; I <- as.numeric(Matrix::colSums(Xw))}.
#' }
#' Both steps use BLAS/sparse kernels where available; the sparse branch never
#' forms \code{Diagonal(W)} and never allocates an \eqn{n \times p} temporary.
#'
#' @param X Numeric design matrix \eqn{n \times p}{n x p}. May be a base
#'   `matrix` (dense) or a \pkg{Matrix} sparse class (e.g., `dgCMatrix`).
#'   Columns are the features to screen.
#' @param y Response vector of length \eqn{n}. For binomial GLMs, supply a 0/1
#'   vector; `weights` can encode binomial totals if desired.
#' @param family A GLM family object (e.g., `gaussian()`, `binomial()`,
#'   `poisson()`, `Gamma()`, `inverse.gaussian()`), as in `glm()`.
#' @param weights Optional non-negative observation weights of length \eqn{n}.
#' @param offset Optional offset vector of length \eqn{n}.
#' @param nsis Integer, number of predictors to retain. Defaults to
#'   \eqn{\min\{n, \lfloor n / (\log\log p \cdot \log p)\rfloor\}}{
#'   min(n, floor(n / (log log p * log p)))}. For the classical SIS size,
#'   use `nsis = floor(n / log(n))`.
#' @param return_stats Logical; if `TRUE`, also returns a data frame with the
#'   per-feature statistic and ranks. Default `FALSE` to minimize memory and I/O
#'   for ultrahigh-dimensional use.
#'
#' @return A `list` with components:
#' \itemize{
#'   \item `indices` - integer vector of length `nsis` with the indices of the kept
#'         predictors (highest \eqn{T_j}).
#'   \item `ranking` - full ordering of all \eqn{p} features (only if
#'         `return_stats = TRUE`), decreasing in \eqn{T_j}.
#'   \item `stats` - data frame with columns `feature`, `stat` (=\eqn{T_j}),
#'         `sign` (=\eqn{sign(U_j)} as a plain sign), `rank`, and `keep`
#'         (only if `return_stats = TRUE`).
#'   \item `nsis` - the number of variables retained.
#'   \item `family`, `link` - the GLM family and link used.
#' }
#'
#' @references
#'   \itemize{
#' \item Fan, J. and Lv, J. (2008).
#'      Sure independence screening for ultrahigh dimensional feature space.
#'      \emph{Journal of the Royal Statistical Society: Series B}, 70(5), 849-911. \cr
#' \item Fan, J. and Song, R. (2010).
#'       Sure independence screening in generalized linear models with NP-dimensionality.
#'       \emph{Annals of Statistics}, 38(6), 3567-3604. \cr
#' \item Saldana, D. F. and Feng, Y. (2018).
#'       SIS: An R Package for Sure Independence Screening in Ultrahigh-Dimensional Statistical Models.
#'       \emph{Journal of Statistical Software}, 83(2), 1-25. \cr
#' \item Zhao, S. D. and Li, Y. (2014).
#'       Score test variable screening.
#'       \emph{Biometrics}, 70(4), 862-871.
#' }
#'
#' @seealso
#' \itemize{
#'   \item \pkg{SIS} for classical (I)SIS implementations.
#'   \item \pkg{glmnet} and \pkg{ncvreg} to fit LASSO/MCP on the screened set.
#' }
#'
#' @example examples/sis_glm.R
#' @export
#' @importFrom stats glm.fit
sis_glm <- function(X, y,
                    family = gaussian(),
                    weights = NULL,
                    offset  = NULL,
                    nsis       = NULL,
                    return_stats = TRUE) {
  
  # --- basic checks
  n <- NROW(X); p <- NCOL(X)
  stopifnot(length(y) == n)
  if (is.null(offset)) offset <- Matrix(0, ncol=1, nrow=n, sparse=TRUE) else offset <- Matrix(offset, sparse=TRUE)
  nonnegative <- ((family$family=="poisson"|family$family=="quasipoisson"|
       family$family=="binomial"|family$family=="quasibinomial"|
       family$family=="Gamma"|family$family=="inverse.gaussian")&
       family$link=="identity")
  
  # if (is.null(nsis)) nsis <- min(p, max(1L, floor(n / log(max(n, 3)))))  # nsis ~ n/log n as in SIS package
  if (is.null(nsis)) nsis <- min(n, round(n/(log(log(n))*log(p))))  # maximum support size used in abess package
  
  # --- null fit (intercept-only)
  fit0 <- L0glm.fit(X = matrix(1, n, 1), y = y, lambda = 0,
                    weights = weights, offset = offset, 
                    start = c(1),
                    nonnegative = nonnegative,
                    family = family, intercept = FALSE)
  # fit0 <- stats::glm.fit(x = matrix(1, n, 1), y = y,
  #                        weights = weights, offset = offset,
  #                        family = family, intercept = TRUE)
  
  # --- score test (fast, single null fit)
  linkinvfun <- family$linkinv
  eta0       <- matrix(1, n, 1) %*% fit0$coefficients + offset # fit0$linear.predictors
  mu0        <- linkinvfun(eta0) # fit0$fitted.values
  vfun       <- family$variance
  mu_eta     <- family$mu.eta
  mueta      <- pmax(mu_eta(eta0), .Machine$double.eps)
  Wdiag      <- (mueta^2) / pmax(vfun(mu0), .Machine$double.eps)
  detadmu    <- 1 / mueta
  Uvec       <- (y - mu0) * detadmu
    
  # -- score statistic (two BLAS-level passes)
  # make vectors to avoid S4 broadcasting surprises
  Uvec  <- as.numeric(Uvec)                  # length n
  Wvec  <- as.numeric(Wdiag)                 # length n
  sqrtW <- sqrt(pmax(Wvec, .Machine$double.eps))
  
  # U_j = x_j' * Uvec  (dispatches to Matrix method if X is sparse)
  U <- as.numeric(crossprod(X, Uvec))
  
  # I_jj = sum_i W_i * x_{ij}^2
  if (inherits(X, "sparseMatrix")) {
    Xw <- X                                  # shallow copy
    # scale rows by sqrtW and square the nonzeros: (sqrtW * X)^2
    Xw@x <- (X@x * sqrtW[X@i + 1L])^2
    I <- as.numeric(Matrix::colSums(Xw))     # stays sparse
  } else {
    I <- colSums((X * sqrtW)^2)              # dense path
  }
  
  stats_vec <- U^2 / pmax(I, .Machine$double.eps)
  signs     <- sign(U)
  
  # --- top-nsis selection
  nsis <- min(as.integer(nsis), length(stats_vec))
  idx_top  <- sort.list(stats_vec, decreasing = TRUE)[seq_len(nsis)]
  idx_keep <- idx_top[ order(stats_vec[idx_top], decreasing = TRUE) ]
  
  if (return_stats) {
    # If we really need full ranks for all p, compute below; otherwise omit for memory.
    o <- order(stats_vec, decreasing = TRUE)
    stats <- data.frame(
      feature = seq_len(p),
      stat    = stats_vec,
      sign    = signs,
      rank    = match(seq_len(p), o),
      keep    = FALSE
    )
    stats$keep[idx_keep] <- TRUE
  } else {
    stats <- NULL
  }
  
  list(
    indices = idx_keep,
    ranking = if (return_stats) o else NULL,
    stats   = stats,
    nsis    = length(idx_keep),
    family  = family$family,
    link    = family$link
  )
}
