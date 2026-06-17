#' @name glm_dist
#' @title Calculate dissimilarity matrix for multivariate data conforming to any GLM 
#' exponential-family distribution
#'
#' @description
#' `glm_dist()` computes a column-by-column dissimilarity matrix for a
#' dense or sparse data matrix \eqn{Y} whose columns are variables
#' and rows are observations and whose distribution conforms to
#' any standard GLM (generalized linear model) exponential family (e.g.
#' Poisson, Binomial, Gaussian, Gamma, Inverse-Gaussian).
#'
#' Two flavours are implemented:
#' \itemize{
#'   \item \strong{Symmetrised unit deviance}
#'         (\code{method = "deviance"}, default) -
#'         equal to \emph{two} times the Jensen-Shannon divergence for any
#'         exponential-family distributed data.  
#'         For the Gaussian(identity) family this becomes
#'         one half of the weighted squared Euclidean distance.
#'         Observation-level precision weights (e.g. Binomial totals,
#'         1/variance for Gaussian) can be taken into account.
#'
#'   \item \strong{Poisson size-factor distance}
#'         (\code{method = "poissondist"}) -
#'         Witten (2011) distance that divides each column by its
#'         size factor before computing the symmetrised Poisson deviance,
#'         treating global scale differences as nuisance.
#'         Numerically identical to
#'         \code{PoiClaClu::PoissonDistance(t(Y), type = "mle",
#'         transform = FALSE)$dd} but computed much faster.
#' }
#' @details
#' \itemize{
#'        \item \strong{Symmetrised unit deviance}
#'                     (\code{method = "deviance"}, default)\cr
#' For two columns *i* and *j* let
#' \deqn{\mu_t^{(i,j)} \;=\;
#'       \frac{w_{ti} y_{ti} + w_{tj} y_{tj}}
#'            {w_{ti} + w_{tj}}}{
#'       \mu = (w_i y_i + w_j y_j)/(w_i + w_j)}
#' be the precision-weighted mean at row *t*.  With
#' \eqn{d(y,\mu,w)} the unit deviance returned by
#' \code{family\$dev.resids}, the distance is
#' \deqn{D_{ij}^{(\mathrm{dev})}
#'       = \sum_t \bigl[d(y_{ti},\mu_t^{(i,j)},w_{ti})
#'                     + d(y_{tj},\mu_t^{(i,j)},w_{tj})\bigr].}
#' For the Gaussian(identity) family
#' \eqn{d(y,\mu,w)=w(y-\mu)^2}, hence
#' \eqn{D_{ij}^{(\mathrm{dev})}= \tfrac12\;
#'       \Vert y_{\bullet i}-y_{\bullet j}\Vert_2^2}, i.e. **one half of
#' the squared Euclidean distance** (so its square-root equals Euclidean
#' distance up to \eqn{\sqrt{1/2}}).
#'
#'      \item \strong{Poisson size-factor distance (Witten 2011)}
#'                   (\code{method = "poissondist"})\cr
#' With column size factors
#' \eqn{s_i = \mathrm{size\_fac}[i]}, row totals
#' \eqn{c_t = y_{ti}+y_{tj}}, and a Gamma-Poisson smoothing constant
#' \eqn{\beta} (default 1), define
#' \deqn{n_{ti} = \frac{s_i c_t}{s_i+s_j},\qquad
#'       d_i = \frac{y_{ti}+\beta}{n_{ti}+\beta}}
#' and analogously for *j*.  The Witten (2011) distance is
#' \deqn{D_{ij}^{(\mathrm{pois})}
#'       = \sum_t\!\bigl[n_{ti}+n_{tj}
#'       - n_{ti}d_i-n_{tj}d_j
#'       + y_{ti}\log d_i + y_{tj}\log d_j\bigr].}
#' Setting \code{size_fac = rowSums(Y)} yields the
#' row-wise variant identical to
#' \code{PoiClaClu::PoissonDistance(x,\dots)} without transposing *x*.
#' }
#'
#' Both distances are twice the Jensen-Shannon divergence of the
#' corresponding exponential-family distributions; for Gaussian this
#' again equals \eqn{\tfrac12\Vert\Delta\Vert_2^2}.
#'
#' @param Y  A numeric matrix (`nrow(Y)` observations x `ncol(Y)` variables).
#'           Sparse `Matrix::dgCMatrix` is accepted. Proportion data should
#'           be provided in case of binomial data.
#' @param weights  Optional `nrow(Y) x ncol(Y)` matrix of precision
#'   weights (e.g. total counts for Binomial, `1/var` for Gaussian).
#'   Only used with `method = "deviance"`.
#' @param family  A **stats** family object (default
#'   `poisson(identity)`).
#' @param method  `"deviance"` or `"poissondist"` (default for family `poisson(identity)`).
#' @param size_fac  Size factors used for method `"poissondist"` (default colSums(Y),
#'                  should be changed to rowSums(Y) if t(Y) is analysed to cluster
#'                  rows rather than columns).
#'
#' @return An object of class **dist**.
#'
#' @references
#' Witten, D. M. (2011). 
#' \emph{Classification and Clustering of Sequencing Data 
#' using a Poisson Model.} 
#' Annals of Applied Statistics 5(4): 2493-2518.
#'
#' @example examples/glm_dist.R
#' @export
#' @author
#' Tom Wenseleers \email{tom.wenseleers@kuleuven.be} \\
#' Laboratory of Socioecology and Social Evolution, KU Leuven \\
#' \url{https://bio.kuleuven.be/eeb/tw/research} \\
#' ORCID: \href{https://orcid.org/0000-0002-1434-861X}{0000-0002-1434-861X}
glm_dist <- function(Y,
                     weights = NULL,
                     family  = poisson(identity),
                     method  = c("deviance", "poissondist"),
                     size_fac = colSums(Y)) {
  
  method <- match.arg(method)
  Y <- as.matrix(Y)                       # ensure matrix
  n <- nrow(Y); p <- ncol(Y)
  
  ## ----------------------------------------------------------------------------
  ## 1.  Jensen-Shannon divergence calculated from the symmetrised unit deviance 
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
