#' Dissimilarities between scalar GLM response channels
#' @md
#'
#' @description
#' Compute pairwise dissimilarities between columns of a response matrix.
#' `method="poissondist"` uses Witten's size-adjusted, smoothed Poisson
#' likelihood gain with a native OpenMP kernel. `method="deviance"` pools
#' each pair at its precision-weighted mean and sums their unit deviances.
#'
#' @param Y Finite numeric observations-by-channels matrix. Sparse matrices
#'   are accepted; response calculations use a dense matrix internally.
#' @param weights Nonnegative case/precision weights: a scalar, a vector with
#'   one entry per observation, or an observations-by-channels matrix.
#'   `NULL` uses unit weights. Binomial observations are proportions and
#'   weights give the trial totals.
#' @param family Scalar GLM family object, default `poisson(identity)`.
#'   Response values must lie in the family's mean/deviance domain.
#' @param method `"deviance"` (default) or `"poissondist"`.
#' @param size_fac Nonnegative channel size factors for `"poissondist"`,
#'   default `colSums(Y)`. A zero factor is allowed for an uninformative
#'   zero-count channel. Counts remain on their original scale.
#' @param nthreads Number of native worker threads. `NULL` uses the detected
#'   physical cores; `1L` requests serial execution. Results have the same
#'   within-pair summation order at every thread count.
#'
#' @details
#' For channels \eqn{a,b}, the pooled-deviance dissimilarity is
#' \deqn{D_{ab}=\sum_i\{d(y_{ia},\bar\mu_i,w_{ia})+
#' d(y_{ib},\bar\mu_i,w_{ib})\},\qquad
#' \bar\mu_i=\frac{w_{ia}y_{ia}+w_{ib}y_{ib}}{w_{ia}+w_{ib}}.}
#' Native pooled-deviance kernels cover Gaussian, Poisson, binomial, Gamma,
#' inverse Gaussian and fixed-shape NB families; custom scalar families use
#' their R `dev.resids` function. NB shape is taken from the family object
#' and held fixed. Rows with zero total weight contribute zero. For a regular scalar
#' exponential family with convex mean potential \eqn{F},
#' \eqn{d(y,\mu,w)=2w B_F(y,\mu)}. The pooled loss is twice a weighted
#' Jensen gap. This connects the dissimilarity to exponential-family
#' likelihood and Bregman clustering (Banerjee et al., 2005).
#' For Gaussian observations with unit weights it is
#' \eqn{\tfrac12\sum_i(y_{ia}-y_{ib})^2}.
#' A generic pooled unit deviance is a dissimilarity; distributional
#' Jensen--Shannon terminology and a universal metric claim do not apply.
#'
#' For the Poisson method let \eqn{s_a,s_b} be the size factors and define
#' \deqn{c_i=\frac{w_{ia}y_{ia}+w_{ib}y_{ib}}
#' {w_{ia}s_a+w_{ib}s_b},\quad m_{ia}=s_a c_i,\quad
#' r_{ia}=\frac{w_{ia}y_{ia}+1}{w_{ia}m_{ia}+1}.}
#' The contribution for channel \eqn{a} is
#' \eqn{w_{ia}\{m_{ia}-m_{ia}r_{ia}+y_{ia}\log r_{ia}\}},
#' with the analogous contribution for \eqn{b}. Unit weights recover
#' Witten's Gamma(1,1)-smoothed likelihood-gain formula. The weighted version
#' applies the same smoothing to the weighted count and expected count.
#' Zero-weight and zero-exposure terms contribute zero. Common rescaling of
#' all size factors leaves the model unchanged.
#'
#' @return A labelled `dist` object in packed lower-triangular storage.
#' @references
#' Witten, D. M. (2011). Classification and Clustering of Sequencing Data
#' using a Poisson Model. Annals of Applied Statistics 5(4), 2493--2518.
#' \doi{10.1214/11-AOAS493}.
#'
#' Banerjee, A., Merugu, S., Dhillon, I. S. and Ghosh, J. (2005).
#' Clustering with Bregman Divergences. JMLR 6, 1705--1749.
#' \url{https://jmlr.org/papers/v6/banerjee05b.html}.
#' @seealso [glm_cluster()]
#' @example examples/glm_dist.R
#' @export
 glm_dist <- function(Y, weights=NULL, family=poisson(identity),
                      method=c("deviance","poissondist"),
                      size_fac=colSums(Y), nthreads=NULL) {
  method<-match.arg(method)
  Y<-L0glm_cluster_matrix(Y)
  w<-L0glm_cluster_weights(weights,nrow(Y),ncol(Y))
  threads<-L0glm_sis_threads(nthreads)
  if (method=="poissondist") {
    s<-L0glm_cluster_sizes(size_fac,Y,w)
    values<-L0glm_poissondist_cpp(Y,w,s,threads)
    return(structure(values,Size=ncol(Y),Labels=colnames(Y),Diag=FALSE,
      Upper=FALSE,method="poissondist",class="dist"))
  }
  if (!inherits(family,"family")) family<-family()
  if (!is.function(family$dev.resids)) stop("family needs scalar unit deviances.",call.=FALSE)
  info<-L0glm_cluster_family(family)
  if(info$id>0L) {
    values<-L0glm_deviance_dist_cpp(Y,w,info$id,info$theta,threads)
    return(structure(values,Size=ncol(Y),Labels=colnames(Y),Diag=FALSE,
      Upper=FALSE,method="deviance",class="dist"))
  }
  W<-matrix(w,nrow(Y),ncol(Y))
  D<-matrix(0,ncol(Y),ncol(Y),dimnames=list(colnames(Y),colnames(Y)))
  if(ncol(Y)>1L) for(j in 2:ncol(Y)) for(i in seq_len(j-1L)) {
    total<-W[,i]+W[,j];good<-total>0
    if(!any(good)) next
    mu<-(W[good,i]*Y[good,i]+W[good,j]*Y[good,j])/total[good]
    value<-sum(family$dev.resids(Y[good,i],mu,W[good,i])+
      family$dev.resids(Y[good,j],mu,W[good,j]))
    if(!is.finite(value) || value< -1e-10) stop("The family returned invalid unit deviances.",call.=FALSE)
    D[i,j]<-D[j,i]<-max(0,value)
  }
  as.dist(D)
}

L0glm_cluster_matrix<-function(Y) {
  Y<-as.matrix(Y)
  if(!(is.numeric(Y)||is.logical(Y)) || !nrow(Y) || !ncol(Y) || any(!is.finite(Y)))
    stop("Y must be a finite, nonempty numeric matrix.",call.=FALSE)
  if(!is.double(Y)) storage.mode(Y)<-"double"
  Y
}
L0glm_cluster_weights<-function(weights,n,m) {
  if(is.null(weights)) return(1)
  if(!(is.numeric(weights)||inherits(weights,"Matrix"))) stop("weights must be numeric.",call.=FALSE)
  if(!is.null(dim(weights)) && !identical(as.integer(dim(weights)),c(n,m)) &&
      !identical(as.integer(dim(weights)),c(n,1L))) stop("weights have incompatible dimensions.",call.=FALSE)
  w<-as.numeric(weights)
  if(!length(w)%in%c(1L,n,n*m) || any(!is.finite(w)|w<0))
    stop("weights must be finite, nonnegative and match Y.",call.=FALSE)
  w
}
L0glm_cluster_sizes<-function(size_fac,Y,w) {
  if(!is.numeric(size_fac) || length(size_fac)!=ncol(Y) || any(!is.finite(size_fac)|size_fac<0))
    stop("size_fac needs one finite nonnegative value per channel.",call.=FALSE)
  if(any(size_fac==0) && any(size_fac==0 & colSums(Y*matrix(w,nrow(Y),ncol(Y)))>0))
    stop("A positive-count channel needs a positive size factor.",call.=FALSE)
  as.numeric(size_fac)
}

L0glm_cluster_family<-function(family) {
  name<-tolower(family$family)
  if(grepl("cox|multinom",name)) stop("Clustering requires independent scalar response families.",call.=FALSE)
  id<-match(name,c("poisson","gaussian","binomial","gamma","inverse.gaussian"),nomatch=0L)
  if(name=="quasipoisson") id<-1L
  if(name=="quasibinomial") id<-3L
  theta<-1
  if(grepl("^negative binomial",name)) {
    theta<-if(is.function(family$getTheta)) family$getTheta(TRUE) else
      get0(".Theta",envir=environment(family$variance),inherits=FALSE)
    if(!is.numeric(theta)||length(theta)!=1L||!is.finite(theta)||theta<=0)
      stop("Negative-binomial clustering needs a fixed positive shape in family.",call.=FALSE)
    id<-6L
  }
  list(id=id,theta=theta,likelihood=id>0L&&!grepl("^quasi",name))
}
