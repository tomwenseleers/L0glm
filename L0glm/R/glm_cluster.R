#' @name glm_cluster
#' @title Hierarchical clustering and aggregation for multivariate GLM 
#' exponential-family data
#'
#' @description
#' \code{glm_cluster()} performs Ward hierarchical clustering on a
#' Jensen–Shannon / deviance dissimilarity matrix from
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
#' Applications include grouping correlated mass traces in GC/MS datasets or 
#' clustering genes with correlated gene-expression patterns in RNA-seq datasets. 
#' 
#' @details
#' Clustering via \emph{Ward’s minimisation of within-cluster deviance},
#' using the \eqn{\sqrt{\mathrm{JS}}}{sqrt(JS)} divergence, is a special
#' case of \strong{agglomerative Bregman clustering}
#' (Telgarsky & Dasgupta 2012); the Jensen–Shannon divergence is itself a
#' Bregman divergence for exponential-family log-likelihoods
#' (Banerjee et al. 2005).  Because \eqn{\sqrt{\mathrm{JS}}} is a true
#' metric (Endres & Schindelin 2003), Ward linkage yields monotone
#' dendrograms (see Nielsen 2016, Example 3).
#'
#' Note that hierarchical clustering is feasible only for moderately large data
#' sets; for very large matrices more scalable Bregman
#' k-means/X-means alternatives could be considered (Banerjee et al. 2005).
#' 
#' @inheritParams glm_dist
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
#' Proceedings of the 29th International Conference on Machine Learning, 59–66. \cr
#' Banerjee, A., Merugu, S., Dhillon, I., & Ghosh, J. (2005).  
#' \emph{Clustering with Bregman Divergences}.  
#' Journal of Machine Learning Research, 6, 1705–1749. \cr
#' Endres, D. M. & Schindelin, J. E. (2003).  
#' \emph{A New Metric for Probability Distributions}.  
#' IEEE Transactions on Information Theory, 49(7), 1858–1860. \cr
#' \emph{Hierarchical Clustering}.
#' In: Nielsen, F., Barbaresco, F. (eds.)
#' \emph{Geometric Theory of Information},
#' Lecture Notes in Computer Science 10589,
#' pp. 369–395. Springer, Cham.
#' \url{https://doi.org/10.1007/978-3-319-21903-5_8}
#' 
#' @examples
#' ## example: clustering and aggregating m/z mass traces in GC/MS datasets
#' data(simulated_spectra)
#' Y <- simulate_GCMS(n = 5000,
#'                    k = 200,
#'                    S = simulated_spectra,
#'                    family = "poisson")$D
#' dim(Y) # 5000 x 258
#' matplot(Y ^ 0.2, type="l", lty=1, col=rgb(0,0,0,0.1),
#'         xlab = "Time index",
#'         ylab = "Mass trace ion intensity ^ 0.2",
#'         main="258 m/z mass channels in original GC/MS dataset")
#' 
#' # aggregate correlated m/z mass channels by clustering using Witten 2011's Poisson size-factor distance
#' Yagg1 <- glm_cluster(Y = Y, family  = poisson(identity), method="poissondist")$Y_agg
#' dim(Yagg1) # 5000 x 10
#' matplot(Yagg1 ^ 0.2, type="l", lty=1,
#'         xlab = "Time index",
#'         ylab = "Mass cluster ion intensity ^ 0.2",
#'         main="258 m/z mass channels clustered into 10 aggregate traces\nusing Witten 2011's Poisson size-factor distance")
#' # NOTE: this uses distance matrix D1 <- glm_dist(Y = Y, family  = poisson(identity), method="poissondist")
#' # which is identical to D1b <- PoiClaClu::PoissonDistance(t(Y), type = "mle", transform = FALSE)$dd
#' 
#' # aggregate correlated m/z mass channels by clustering using Jensen–Shannon divergence / symmetrised deviance
#' Yagg2 <- glm_cluster(Y = Y, family  = poisson(identity), method="deviance")$Y_agg
#' dim(Yagg2) # 5000 x 4
#' matplot(Yagg2 ^ 0.2, type="l", lty=1,
#'         xlab = "Time index",
#'         ylab = "Mass cluster ion intensity ^ 0.2",
#'         main="258 m/z mass channels clustered into 4 aggregate traces\nusing Jensen–Shannon deviance distance")
#'
#' @importFrom fastcluster hclust
#' @importFrom dynamicTreeCut cutreeDynamic
#' @export
glm_cluster <- function(Y,
                        weights = NULL,
                        family  = poisson(identity),
                        method  = if (family$family == "poisson") "poissondist" else "deviance") {
  
  ## 1. calculate Jensen-Shannon divergence (= squared Jensen-Shannon distance matrix)
  D <- glm_dist(Y, weights = weights, family = family, method = method)
  
  ## 2. Ward clustering (uses squared distances)
  hc <- fastcluster::hclust(D, method = "ward.D")  # Ward on squared JS divergences
  
  ## 3. adaptive tree cut
  clusterids <- dynamicTreeCut::cutreeDynamic(hc,
                                      distM = as.matrix(D),
                                      deepSplit = 2,
                                      minClusterSize = 1)
  
  ## 4. aggregation ----------------------------------------------------
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