#' @name glm_cluster
#' @title Hierarchical clustering and aggregation for multivariate GLM 
#' exponential-family data
#'
#' @description
#' \code{glm_cluster()} performs Ward hierarchical clustering on a
#' deviance (Jensen-Shannon divergence) dissimilarity matrix from
#' \code{\link{glm_dist}}, adaptively cuts the dendrogram using 
#' (\pkg{dynamicTreeCut}), and then calls \code{glm_aggregate()} to 
#' aggregate columns belonging to the same cluster.
#'
#' \code{glm_aggregate()} can also be called directly when you already
#' have a vector of cluster memberships and want to reuse it across
#' multiple datasets.
#' 
#' This function can be used for clustering or as a feature-grouping 
#' dimensionality-reduction technique for
#' multivariate GLM data conforming to any exponential-family distribution 
#' (e.g. Poisson, Binomial, Gaussian, Gamma, Inverse-Gaussian). 
#' For aggregation, we use summation for count families and weighted or 
#' unweighted means for the other distributions. 
#' 
#' Applications, among others, include grouping correlated mass traces in 
#' GC/MS datasets or clustering genes with correlated gene-expression 
#' patterns in RNA-seq datasets. 
#' 
#' @details
#' Clustering via \emph{Ward's minimisation of within-cluster deviance},
#' using the \eqn{\sqrt{\mathrm{JS}}}{sqrt(JS)} divergence, is a special
#' case of \strong{agglomerative Bregman clustering}
#' (Telgarsky & Dasgupta 2012); the Jensen-Shannon divergence is itself a
#' Bregman divergence for exponential-family log-likelihoods
#' (Banerjee et al. 2005).  Because \eqn{\sqrt{\mathrm{JS}}} is a true
#' metric (Endres & Schindelin 2003), Ward linkage yields monotone
#' dendrograms (see Nielsen 2016, Example 3). The
#' dendrogram height equals the cumulative increase in cluster
#' Bregman information.
#' 
#' \strong{Distance formulas used internally} (see \code{\link{glm_dist}}) \cr 
#' For two columns \emph{i} and \emph{j} we write
#' \eqn{y_{ti}}, \eqn{y_{tj}} for their values at row \emph{t} and
#' \eqn{w_{ti}}, \eqn{w_{tj}} for optional precision weights.
#'
#' \itemize{
#' \item \strong{Symmetrised unit deviance}
#'                     (\code{method = "deviance"}, default)\cr
#' \deqn{
#'   D_{ij}^{\text{dev}} \;=\;
#'   \sum_t
#'     \Bigl[ d\!\bigl(y_{ti},\;\bar\mu_t,\;w_{ti}\bigr) \;+\;
#'            d\!\bigl(y_{tj},\;\bar\mu_t,\;w_{tj}\bigr) \Bigr], \qquad
#'   \bar\mu_t \;=\;
#'   \frac{w_{ti}y_{ti}+w_{tj}y_{tj}}{w_{ti}+w_{tj}},
#' }{
#'   D = Sigma_t [ d(y_i,mu_bar,w_i)+d(y_j,mu_bar,w_j) ],  mu_bar = (w_i y_i + w_j y_j)/(w_i+w_j)
#' }
#' where \eqn{d(y,\mu,w)} is the unit deviance returned by
#' \code{family$dev.resids(y, mu, w)}. For the Gaussian(identity) family
#' \eqn{d(y,\mu,w)=w(y-\mu)^2}, hence
#' \eqn{D_{ij}^{\text{dev}}=\tfrac12\sum_t w_t(y_{ti}-y_{tj})^2} -
#' exactly \emph{one-half} the (weighted) squared Euclidean distance, so
#' \eqn{\sqrt{\mathrm{JS}}} is Euclidean up to the constant
#' \eqn{\sqrt{1/2}}.
#'
#' \item \strong{Poisson size-factor distance (Witten 2011)}
#'                   (\code{method = "poissondist"})\cr
#' With column size factors \eqn{s_i}, row totals
#' \eqn{c_t=y_{ti}+y_{tj}} and smoothing \eqn{\beta} (default 1):
#' \deqn{
#'   n_{ti}= \frac{s_i\,c_t}{s_i+s_j}, \qquad
#'   d_i   = \frac{y_{ti}+\beta}{n_{ti}+\beta}
#' }
#' \deqn{
#'   D_{ij}^{\mathrm{pois}} =
#'     \sum_t \Bigl[ n_{ti}+n_{tj} - n_{ti}d_i-n_{tj}d_j +
#'                     y_{ti}\log d_i + y_{tj}\log d_j \Bigr].
#' }
#' }
#' 
#' Both expressions equal \emph{twice} the Jensen-Shannon divergence
#' between the empirical distributions of columns \emph{i} and \emph{j};
#' taking the square-root gives a proper metric.  Ward's linkage then minimises 
#' the increase in within-cluster Bregman information (see
#' Telgarsky & Dasgupta 2012).
#'
#' Note that hierarchical clustering is feasible only for moderately large data
#' sets; for very large matrices more scalable Bregman
#' k-means/X-means alternatives could be considered (Banerjee et al. 2005).
#' 
#' @inheritParams glm_dist
#' @param rowbins Positive integer.  If \code{rowbins > 1} the
#'   function first aggregates every \code{rowbins} successive rows
#'   (time points, samples) to speed up distance computation and
#'   clustering; clusters are inferred on the binned matrix, but the
#'   final \code{Y_agg} is computed from the original (unbinned) data.
#'   Default \code{rowbins = 1} skips this step.
#' @param clusterids Integer vector of cluster labels (length =
#'   \code{ncol(Y)}); required by \code{glm_aggregate()}, returned by
#'   \code{glm_cluster()}.
#'
#' @return A list   
#' \describe{
#'   \item{D}{\code{dist} object returned by \code{glm_dist()}}
#'   \item{hc}{\code{hclust} object (using method Ward.D)}
#'   \item{clusterids}{integer vector with cluster ID memberships returned by 
#'                     \code{cutreeDynamic()}}
#'   \item{Y_agg}{matrix with per-cluster aggregated columns (summed for poisson 
#'                and binomial, averaged within clusters for the other families)}
#'   \item{weights_agg}{aggregated weights (NULL if \code{weights=NULL})}
#' }
#' 
#' @references
#' Telgarsky, M. & Dasgupta, S. (2012).  
#' \emph{Agglomerative Bregman Clustering}.  
#' Proceedings of the 29th International Conference on Machine Learning, 59-66. \cr
#' Banerjee, A., Merugu, S., Dhillon, I., & Ghosh, J. (2005).  
#' \emph{Clustering with Bregman Divergences}.  
#' Journal of Machine Learning Research, 6, 1705-1749. \cr
#' Endres, D. M. & Schindelin, J. E. (2003).  
#' \emph{A New Metric for Probability Distributions}.  
#' IEEE Transactions on Information Theory, 49(7), 1858-1860. \cr
#' \emph{Hierarchical Clustering}.
#' In: Nielsen, F., Barbaresco, F. (eds.)
#' \emph{Geometric Theory of Information},
#' Lecture Notes in Computer Science 10589,
#' pp. 369-395. Springer, Cham.
#' \doi{10.1007/978-3-319-21903-5_8} \cr
#' Witten, D. M. (2011). 
#' \emph{Classification and Clustering of Sequencing Data 
#' using a Poisson Model.} 
#' Annals of Applied Statistics 5(4): 2493-2518.
#' 
#' @example examples/glm_cluster.R
#' @importFrom fastcluster hclust
#' @importFrom dynamicTreeCut cutreeDynamic
#' @export
#' @author
#' Tom Wenseleers \email{tom.wenseleers@kuleuven.be} \\
#' Laboratory of Socioecology and Social Evolution, KU Leuven \\
#' \url{https://bio.kuleuven.be/eeb/tw/research} \\
#' ORCID: \href{https://orcid.org/0000-0002-1434-861X}{0000-0002-1434-861X}
glm_cluster <- function(Y,
                        weights = NULL,
                        family  = poisson(identity),
                        method  = if (family$family == "poisson") "poissondist" else "deviance",
                        size_fac = colSums(Y),
                        rowbins   = 1L) {
  
  if (rowbins < 1 || rowbins != as.integer(rowbins))
    stop("`rowbins` must be a positive integer.")
  
  ## 0. (optional) bin rows for speed -------------------------------------
  Y_work  <- Y                  # default: original data
  W_work  <- weights
  if (rowbins > 1L) {
    row_clust <- rep(seq_len(ceiling(nrow(Y)/rowbins)),
                     each = rowbins, length.out = nrow(Y))
    
    tmp <- glm_aggregate(t(Y),
                         if (!is.null(weights)) t(weights) else NULL,
                         family      = family,
                         clusterids  = row_clust)
    
    Y_work <- t(tmp$Y_agg)
    if (!is.null(weights))
      W_work <- if (!is.null(tmp$weights_agg)) t(tmp$weights_agg) else NULL
  }
  
  
  ## 1. distance & clustering on (possibly) row-binned matrix -------------
  ### 1a. calculate symmetrised deviance (2 x Jensen-Shannon divergence) distance matrix 
  D <- glm_dist(Y_work, weights = W_work, family = family, 
                method = method, 
                size_fac = if (rowbins == 1L) size_fac else colSums(Y_work))
  
  ### 1b. Ward clustering (uses squared distances)
  hc <- fastcluster::hclust(D, method = "ward.D")  # Ward on squared JS divergences
  
  ### 1c. adaptive tree cut
  clusterids <- dynamicTreeCut::cutreeDynamic(hc,
                                      distM = as.matrix(D),
                                      deepSplit = 2,
                                      minClusterSize = 1)
  
  ## 2. aggregate ORIGINAL matrix with those column cluster ---------------
  agg <- glm_aggregate(Y, weights, family, clusterids)
  
  list(D = D,
       hc = hc,
       clusterids = clusterids,
       Y_agg = agg$Y_agg,
       weights_agg = agg$weights_agg)
}



#' @describeIn glm_cluster  Aggregate columns of \eqn{Y} over given cluster IDs.
#' @export
glm_aggregate <- function(Y,
                          weights     = NULL,
                          family      = poisson(identity),
                          clusterids) {
  
  if (length(clusterids) != ncol(Y))
    stop("clusterids must have length ncol(Y)")
  
  P <- Matrix::sparseMatrix(i = seq_along(clusterids),
                            j = clusterids,
                            x = 1L,
                            dims = c(length(clusterids),
                                     max(clusterids)))
  
  weights_agg <- NULL
  
  if (family$family %in% c("gaussian", "Gamma", "inverse.gaussian")) {
    
    if (is.null(weights)) {
      Y_agg <- (Y %*% P) / rep(colSums(P), each = nrow(Y))
    } else {
      W      <- as.matrix(weights)
      num    <- (Y * W) %*% P
      den    <-  W       %*% P
      Y_agg  <- num / den
      weights_agg <- den
    }
    
  } else if (family$family == "binomial") {
    
    if (is.null(weights))
      stop("For Binomial data please supply 'weights' (total counts).")
    
    num    <- (Y * weights) %*% P
    den    <-  weights       %*% P
    Y_agg  <- num / den
    weights_agg <- den
    
  } else {                      # Poisson or other count family
    Y_agg <- Y %*% P
    if (!is.null(weights))
      weights_agg <- weights %*% P
  }
  
  list(Y_agg       = as.matrix(Y_agg),
       weights_agg = if (!is.null(weights_agg)) as.matrix(weights_agg)
       else NULL)
}
