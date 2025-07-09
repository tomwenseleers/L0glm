#' @name glm_dist
#' @title Calculate dissimilarity matrix for multivariate data conforming to any GLM 
#' exponential-family distribution
#'
#' @description
#' `glm_dist()` computes a column–by–column dissimilarity matrix for a
#' dense or sparse data matrix \eqn{Y} whose columns are variables
#' and rows are observations and whose distribution conforms to
#' any standard GLM (generalized linear model) exponential family (e.g.
#' Poisson, Binomial, Gaussian, Gamma, Inverse-Gaussian).
#'
#' Two flavours are implemented:
#'
#' * **Symmetrised unit deviance** (`method = "deviance"`, default) —
#'   Jensen–Shannon divergence / Bregman information, calculated from the 
#'   symmetrised unit deviance of any exponential-family GLM distribution. 
#'   Reduces to squared Euclidean distance for the Gaussian case but more general
#'   as it can be used as a squared distance metric also for any of the other
#'   exponential-family distributions. Accepts observation-level precision 
#'   weights, e.g. total count weights in case of Binomial data.
#' * **Poisson–size-factor distance**
#'   (`method = "poissondist"`) — the Witten (2011)
#'   distance that removes a column-specific size factor *before*
#'   computing the symmetrised Poisson deviance.  This distance
#'   assumes that scaling differences across outcomes are nuisance). Equivalent 
#'   to `PoiClaClu::PoissonDistance(t(Y), type = "mle", transform = FALSE)$dd` but
#'   slightly faster.
#'
#' @param Y  A numeric matrix (`nrow(Y)` observations × `ncol(Y)` variables).
#'           Sparse `Matrix::dgCMatrix` is accepted. Proportion data should
#'           be provided in case of binomial data.
#' @param weights  Optional `nrow(Y) × ncol(Y)` matrix of precision
#'   weights (e.g. total counts for Binomial, `1/var` for Gaussian).
#'   Only used with `method = "deviance"`.
#' @param family  A **stats** family object (default
#'   `poisson(identity)`).  Ignored when `method = "poissondist"`.
#' @param method  `"deviance"` or `"poissondist"` (default for Poisson).
#'
#' @return An object of class **dist**.
#'
#' @references
#' Witten, D. M. (2011). *Classification and Clustering of Sequencing Data 
#' using a Poisson Model.* Ann. Appl. Stat. 5(4): 2493–2518.
#'
#' @examples
#' ## example: calculating distance among different m/z mass traces in GC/MS datasets
#' data(simulated_spectra)
#' Y <- simulate_GCMS(n = 5000,
#'                    k = 200,
#'                    S = simulated_spectra,
#'                    family = "poisson")$D
#' dim(Y) # 5000 x 258
#' D1 <- glm_dist(Y = Y, family  = poisson(identity), method="poissondist")
#' # identical to D1b <- PoiClaClu::PoissonDistance(t(Y), type = "mle", transform = FALSE)$dd
#' 
#' @export
glm_dist <- function(Y,
                     weights = NULL,
                     family  = poisson(identity),
                     method  = c("deviance", "poissondist")) {
  
  method <- match.arg(method)
  Y <- as.matrix(Y)                       # ensure matrix
  n <- nrow(Y); p <- ncol(Y)
  
  ## ----------------------------------------------------------------------------
  ## 1.  Jensen–Shannon divergence calculated from the symmetrised unit deviance 
  ## ----------------------------------------------------------------------------
  if (method == "deviance") {
    D <- matrix(0, p, p)
    devfun <- family$dev.resids
    
    if (is.null(weights)) {
      for (j in 2:p) for (i in 1:(j - 1)) {
        mu <- 0.5 * (Y[, i] + Y[, j])
        D[i, j] <- D[j, i] <-
          sum(devfun(Y[, i], mu, 1) +
                devfun(Y[, j], mu, 1))
      }
    } else {
      W <- as.matrix(weights)
      for (j in 2:p) for (i in 1:(j - 1)) {
        wi <- W[, i]; wj <- W[, j]
        yi <- Y[, i]; yj <- Y[, j]
        mu <- (wi * yi + wj * yj) / (wi + wj)
        D[i, j] <- D[j, i] <-
          sum(devfun(yi, mu, wi) +
                devfun(yj, mu, wj))
      }
    }
    return(as.dist(D))
  }
  
  ## ----------------------------------------------------------------------------
  ## 2.  Poisson size-factor distance (Witten 2011)                     
  ## ----------------------------------------------------------------------------
  if (method == "poissondist") {
    size_fac <- colSums(Y)
    beta <- 1L   # default of PoiClaClu::PoissonDistance , here hard coded
    D <- matrix(0, p, p)
    
    for (j in 2:p) {
      yj <- Y[, j]; sj <- size_fac[j]
      for (i in 1:(j - 1)) {
        yi <- Y[, i]; si <- size_fac[i]
        
        cvec <- yi + yj
        if (all(cvec == 0)) { D[i, j] <- 0; next }
        
        n_i <- si / (si + sj) * cvec
        n_j <- sj / (si + sj) * cvec
        
        di <- (yi + beta) / (n_i + beta)
        dj <- (yj + beta) / (n_j + beta)
        
        D[i, j] <- D[j, i] <-
          sum(cvec - n_i * di - n_j * dj +
                yi * log(di) + yj * log(dj))
      }
    }
    return(as.dist(D))
  }
}