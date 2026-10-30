#' Likelihood-based clustering and aggregation of GLM response channels
#' @md
#'
#' @description
#' Cluster columns of a scalar-family response matrix, select a tree cut by an
#' information criterion, and aggregate the original, unbinned channels.
#' Native OpenMP kernels calculate distances, Bregman merges and profile losses.
#'
#' @inheritParams glm_dist
#' @param rowbins Positive integer, default `1L`. Values above one sum adjacent
#'   rows for clustering only; use this acceleration for ordered Poisson scans
#'   with unit weights. Returned aggregates retain every original row.
#' @param linkage `"ward"` or `"bregman"`. The size-adjusted Poisson distance
#'   defaults to the Ward hierarchy. Regular scalar GLM likelihoods default
#'   to Bregman merges; [signal_families] use Ward likelihood-ratio distances.
#' @param criterion `"gic"` (default), `"aic"`, `"bic"`, `"ebic"`, `"mbic"`,
#'   or `"dynamic"` for dynamicTreeCut on a Ward hierarchy. Quasibinomial and
#'   custom unit-deviance families default to dynamic cutting. For quasipoisson,
#'   supply a common dispersion to obtain quasi-information criteria. Supplying `k`
#'   requests that fixed number of groups.
#' @param k Optional fixed number of groups, between one and `ncol(Y)`.
#' @param k.max Largest number of groups considered by criterion selection,
#'   default `32L`, capped at `ncol(Y)`.
#' @param dispersion Known positive scalar or per-channel dispersion for
#'   Gaussian, Gamma and inverse Gaussian clustering. Default `1` gives unit
#'   dispersion. It scales clustering precision and is held fixed across cuts.
#'   Poisson and binomial use dispersion one; NB uses its fixed family shape.
#' @param clusterids One positive integer label per channel. Labels need not be
#'   consecutive; `glm_aggregate()` puts groups in increasing label order.
#'
#' @details
#' For a partition \eqn{\mathcal G}, the raw-mean clustering model assigns a
#' free mean \eqn{\mu_{ig}} at each row \eqn{i} of group \eqn{g}. Its profile
#' distortion is
#' \deqn{L(\mathcal G)=\sum_g\sum_{j\in G_g}\sum_i
#'      \frac{w_{ij}}{\phi_j}d(y_{ij},\widehat\mu_{ig},1),\qquad
#' \widehat\mu_{ig}=\frac{\sum_{j\in G_g}w_{ij}y_{ij}/\phi_j}
#' {\sum_{j\in G_g}w_{ij}/\phi_j}.}
#' Bregman clustering greedily merges the pair minimising
#' \eqn{L(A\cup B)-L(A)-L(B)}. The implementation retains sufficient
#' statistics and a priority queue, recalculating the true likelihood-merge
#' cost after each merge. This is agglomerative Bregman information loss
#' (Banerjee et al., 2005; Telgarsky and Dasgupta, 2012).
#'
#' For `method="poissondist"`, clustering profiles global channel intensity
#' using fixed size factors: \eqn{\mu_{ij}=s_j c_{ig}} and
#' \eqn{\widehat c_{ig}=\sum_j w_{ij}y_{ij}/\sum_j w_{ij}s_j}.
#' The Ward route uses Witten's smoothed pairwise distance followed by
#' `fastcluster::hclust(method="ward.D")`. The Bregman route uses the
#' unsmoothed Poisson likelihood-merge cost for this size-factor model.
#' Ward's recurrence is exact for squared-Euclidean sum-of-squares geometry;
#' general deviance dissimilarities give a computational Ward hierarchy.
#' Both routes score their candidate cuts under the specified profile model.
#' Thus the same size factors adjust both linkages. The factors are nuisance
#' exposures used to compare profile shapes; the regression stages continue
#' to use raw counts, and the returned aggregates sum those raw counts.
#'
#' For each candidate cut the criterion is
#' \deqn{L(\mathcal G)-2\sum_g m_g\log(m_g/m)+c(n,n)q_{\mathcal G}.}
#' Here \eqn{m_g} is the number of channels in group \eqn{g}, the second term
#' is the fitted classification-mixture allocation loss, and
#' \eqn{q_{\mathcal G}} counts informative row/group means plus \eqn{K-1}
#' mixing proportions. \eqn{c} is the package IC penalty and \eqn{n} is the
#' original number of rows (at least three for scoring). The saturated
#' likelihood constant is common to all cuts and omitted. For GIC,
#' \eqn{c(n,n)=\log(n)\log\log(n)}. This selects a cut within the fitted
#' hierarchy; the regression likelihood and predictor count used by
#' [L0glm.fit()] are separate. Dispersion and NB shape are fixed here rather
#' than re-estimated at every cut. Gaussian/Gamma/inverse Gaussian dispersion
#' should be supplied when its known scale differs from one.
#'
#' Native scalar kernels support Gaussian, Poisson, binomial proportions,
#' Gamma, inverse Gaussian and fixed-shape negative binomial families. Joint
#' multinomial and Cox responses have different likelihood structures.
#' Quasipoisson and quasibinomial have native working-deviance kernels and
#' support fixed cuts or Ward/dynamic cutting.
#' Custom scalar unit deviances remain available through [glm_dist()] and
#' the Ward/dynamic route.
#' NB1, common-scale Gamma and censored Gamma use the fixed-parameter exact
#' likelihoods in [signal_families]. Their pooled means are optimized with a
#' native scalar score solve, and the saturated reference uses separate
#' optimized means. They support Ward clustering with IC, fixed or dynamic cuts.
#'
#' Aggregation summarises the original columns. Poisson counts and their
#' means add. Shared Poisson case weights stay unchanged; differing weights
#' within a group cannot be combined into one exact weighted Poisson model.
#' Gaussian means use precision-weighted averaging and summed precisions,
#' including the supplied dispersion in those precisions.
#' Binomial proportions use trial-weighted proportions. For Gamma, inverse
#' Gaussian and NB, the returned summaries do not generally have the same
#' fitted GLM family/link/dispersion. Clustering support therefore does not
#' imply that automatic response aggregation is valid for every GLM.
#'
#' @return `glm_cluster()` returns a list with `D` (packed distances for Ward),
#'   `hc` (hierarchy; `NULL` for one channel), `clusterids`, `Y_agg`,
#'   `weights_agg`, `history` (cut losses, degrees of freedom and criterion),
#'   `centres` on the clustering row grid, `size_fac`, `linkage`, `criterion`
#'   and `rowbins`. Incompatible count weights give `Y_agg=NULL` and an
#'   explanatory `aggregation_note`; clustering itself remains available.
#'   `glm_aggregate()` returns `Y_agg` and `weights_agg`.
#'
#' @references
#' Banerjee, A. et al. (2005). Clustering with Bregman Divergences.
#' JMLR 6, 1705--1749. \url{https://jmlr.org/papers/v6/banerjee05b.html}.
#'
#' Telgarsky, M. and Dasgupta, S. (2012). Agglomerative Bregman Clustering.
#' ICML. \url{https://icml.cc/2012/papers/758.pdf}.
#'
#' Witten, D. M. (2011). Classification and Clustering of Sequencing Data
#' using a Poisson Model. AOAS 5(4), 2493--2518. \doi{10.1214/11-AOAS493}.
#' @example examples/glm_cluster.R
#' @seealso [signal_families] for NB1, common-scale Gamma and censored Gamma
#'   clustering with exact likelihood-ratio distances and a Ward hierarchy.
#' @importFrom fastcluster hclust
#' @importFrom dynamicTreeCut cutreeDynamic
#' @export
 glm_cluster<-function(Y,weights=NULL,family=poisson(identity),
                       method=if(family$family=="poisson") "poissondist" else "deviance",
                       size_fac=colSums(Y),rowbins=1L,nthreads=NULL,
                       linkage=NULL,criterion=NULL,k=NULL,k.max=32L,dispersion=1) {
  Y<-L0glm_cluster_matrix(Y);n<-nrow(Y);m<-ncol(Y)
  if(!inherits(family,"family")) family<-family()
  method<-match.arg(method,c("deviance","poissondist"));info<-L0glm_cluster_family(family)
  if(is.null(linkage)) linkage<-if(method=="poissondist"||!info$likelihood||!is.null(info$scalar)||grepl("^quasi",family$family)) "ward" else "bregman"
  if(is.null(criterion)) criterion<-if(info$likelihood) "gic" else "dynamic"
  linkage<-match.arg(linkage,c("ward","bregman"))
  criterion<-match.arg(tolower(criterion),c("gic","bic","aic","ebic","mbic","dynamic"))
  if(!info$likelihood&&criterion!="dynamic"&&is.null(k))
    stop("This family has no clustering likelihood IC; use fixed k or a dynamic cut.",call.=FALSE)
  valid_integer<-function(x) is.numeric(x)&&length(x)==1L&&is.finite(x)&&x>=1&&x==floor(x)&&x<=.Machine$integer.max
  if(!valid_integer(rowbins)||!valid_integer(k.max)||(!is.null(k)&&(!valid_integer(k)||k>m)))
    stop("rowbins, k.max and k must be valid positive integers.",call.=FALSE)
  w<-L0glm_cluster_weights(weights,n,m);threads<-L0glm_sis_threads(nthreads)
  if(method=="poissondist" && info$id!=1L) stop("poissondist clustering requires a Poisson family.",call.=FALSE)
  if(!is.numeric(dispersion)||!length(dispersion)%in%c(1L,m)||any(!is.finite(dispersion)|dispersion<=0))
    stop("dispersion must be positive and scalar or per-channel.",call.=FALSE)
  if(info$likelihood&&info$id%in%c(1L,3L,6L)&&family$family!="quasipoisson"&&any(dispersion!=1)) stop("This family has fixed dispersion one.",call.=FALSE)
  if(info$id==0L && is.null(info$scalar) && (linkage!="ward"||criterion!="dynamic") && is.null(k))
    stop("This scalar family supports Ward clustering with a fixed or dynamic cut.",call.=FALSE)
  if(info$id==0L && linkage=="bregman") stop("No native Bregman kernel for this scalar family.",call.=FALSE)
  if(!is.null(info$scalar)&&any(dispersion!=1)) stop("Pass the fixed noise parameter in the family constructor.",call.=FALSE)
  profile<-function(labels) {
    if(is.null(info$scalar)) return(L0glm_cluster_profile_cpp(work,precision,sizes,labels,info$id,info$theta,threads))
    s<-info$scalar
    L0glm_scalar_profile_cpp(work,precision,labels,s$kind,s$parameter,s$censor,threads)
  }
  sizes<-if(method=="poissondist") L0glm_cluster_sizes(size_fac,Y,w) else rep(1,m)
  work<-Y;precision<-if(length(dispersion)==1L) w/dispersion else
    as.numeric(sweep(matrix(w,n,m),2L,dispersion,"/"))
  if(rowbins>1L) {
    if(info$id!=1L||any(w!=1)) stop("Row binning requires Poisson data with unit weights.",call.=FALSE)
    if(length(dispersion)!=1L) stop("Row binning requires a common dispersion.",call.=FALSE)
    bins<-rep(seq_len(ceiling(n/rowbins)),each=rowbins,length.out=n)
    work<-rowsum(Y,bins,reorder=FALSE);precision<-1/dispersion
  }
  D<-hc<-NULL
  if(m>1L) {
    if(linkage=="ward") {
      D<-glm_dist(work,weights=precision,family=family,method=method,size_fac=sizes,nthreads=threads)
      hc<-fastcluster::hclust(D,method="ward.D")
    } else {
      tree<-L0glm_bregman_tree_cpp(work,precision,sizes,info$id,info$theta,threads)
      hc<-structure(c(tree,list(labels=colnames(Y),method="Bregman likelihood merge",call=match.call(),dist.method=method)),class="hclust")
    }
  }
  history<-NULL;centres<-NULL
  if(criterion=="dynamic" && is.null(k) && m>1L) {
    if(linkage!="ward") stop("Dynamic cuts use a Ward hierarchy; choose an IC or fixed k for Bregman.",call.=FALSE)
    group<-dynamicTreeCut::cutreeDynamic(hc,distM=as.matrix(D),deepSplit=2,minClusterSize=1,verbose=0)
    missing<-which(group==0L);group[missing]<-max(group)+seq_along(missing)
    group<-match(group,sort(unique(group)))
  } else {
    ks<-if(m==1L) 1L else if(!is.null(k)) as.integer(k) else seq_len(min(k.max,m))
    labels<-if(m==1L) matrix(1L,1L,1L) else matrix(stats::cutree(hc,k=ks),m,length(ks))
    if(info$id>0L || !is.null(info$scalar)) {
      loss<-profile(labels)
      mixing<-vapply(seq_along(ks),function(i) {z<-tabulate(labels[,i]);-2*sum(z*log(z/m))},numeric(1))
      df<-loss$df+ks-1L;penalty<-if(criterion=="dynamic") "gic" else criterion
      score<-if(info$likelihood) get.IC(-(loss$deviance+mixing)/2,max(3,n),max(2,n),df)[[penalty]] else rep(NA_real_,length(ks))
      history<-data.frame(k=ks,deviance=loss$deviance,allocation_loss=mixing,df=df,IC=score)
      best<-if(info$likelihood) which.min(score) else 1L
    } else best<-1L
    group<-labels[,best]
  }
  if(info$id>0L || !is.null(info$scalar)) centres<-profile(matrix(as.integer(group),m,1L))$centres
  compatible<-TRUE
  if((info$id%in%c(1L,6L)||!is.null(info$scalar))&&length(w)>n) {
    W<-matrix(w,n,m)
    compatible<-all(vapply(split(seq_len(m),group),function(js) all(W[,js,drop=FALSE]==W[,js[1L]]),logical(1)))
  }
  agg<-if(compatible) glm_aggregate(Y,if(info$id==2L) precision else weights,family,group) else
    list(Y_agg=NULL,weights_agg=NULL,aggregation_note="Count weights differ within groups; use the clustering labels with the original weighted responses.")
  c(list(D=D,hc=hc,clusterids=as.integer(group),history=history,centres=centres,
    size_fac=if(method=="poissondist") sizes else NULL,linkage=linkage,
    criterion=if(!is.null(k)) "fixed" else criterion,rowbins=as.integer(rowbins)),agg)
}

#' @describeIn glm_cluster Summarise the original channels using supplied groups.
#' @export
 glm_aggregate<-function(Y,weights=NULL,family=poisson(identity),clusterids) {
  Y<-L0glm_cluster_matrix(Y);n<-nrow(Y);m<-ncol(Y)
  if(!inherits(family,"family")) family<-family()
  if(!is.numeric(clusterids)||length(clusterids)!=m||any(!is.finite(clusterids)|clusterids<1|clusterids!=floor(clusterids)))
    stop("clusterids needs one positive integer per channel.",call.=FALSE)
  group<-match(clusterids,sort(unique(clusterids)));k<-max(group)
  P<-Matrix::sparseMatrix(i=seq_len(m),j=group,x=1,dims=c(m,k))
  w<-L0glm_cluster_weights(weights,n,m);kind<-tolower(family$family)
  count<-kind%in%c("poisson","quasipoisson","nb1","gamma_scale")||grepl("^negative binomial",kind)
  if(count) {
    Y_agg<-as.matrix(Y%*%P);weights_agg<-NULL
    if(!is.null(weights)) {
      W<-matrix(w,n,m);weights_agg<-matrix(0,n,k)
      for(g in seq_len(k)) {
        js<-which(group==g);wg<-W[,js[1L]]
        if(any(W[,js,drop=FALSE]!=wg)) stop("Count aggregation requires shared row weights within each group.",call.=FALSE)
        weights_agg[,g]<-wg
      }
    }
  } else {
    W<-matrix(w,n,m);den<-as.matrix(W%*%P)
    Y_agg<-as.matrix((Y*W)%*%P)/pmax(den,.Machine$double.xmin)
    Y_agg[den==0]<-0;weights_agg<-den
  }
  colnames(Y_agg)<-paste0("cluster_",seq_len(k));rownames(Y_agg)<-rownames(Y)
  if(!is.null(weights_agg)) dimnames(weights_agg)<-dimnames(Y_agg)
  list(Y_agg=Y_agg,weights_agg=weights_agg)
}
