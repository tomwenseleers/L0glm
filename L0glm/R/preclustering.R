#' Preclustering response channels before sparse Poisson fitting
#' @name preclustering
#' @md
#' @description
#' `L0glm.fit(precluster=TRUE)` and `L0nnpois(precluster=TRUE)` cluster and sum
#' response channels, screen predictors on the reduced responses, then refine
#' the sparse fit on every original channel. The workflow supports
#' nonnegative identity-link Poisson models. Standalone [glm_cluster()] also
#' supports scalar-family Bregman/deviance clustering.
#'
#' @details
#' For conditionally independent Poisson channels with an identity link,
#' \deqn{Y^A_{ig}=\sum_{j\in G_g}Y_{ij}\sim
#' \mathrm{Poisson}(x_i^T\Gamma_g+o^A_{ig}),\qquad
#' \Gamma_g=\sum_{j\in G_g}\beta_j,\quad o^A_{ig}=\sum_{j\in G_g}o_{ij}.}
#' Nonnegativity preserves row support under summation. The cluster screen
#' uses raw counts and a permissive `screen.factor * lambda` (default `0.1`).
#' Original-channel refinement uses the requested penalty, calculated once
#' from the original predictor universe. Its default is two grouped passes;
#' an explicitly supplied `group.maxit` takes precedence unless
#' `control.precluster$refit.passes` is supplied. `relaxed` retains its usual
#' selected-support coefficient-refitting meaning.
#'
#' With `start=NULL`, original-channel refinement broadcasts each candidate's
#' summed aggregate coefficient divided by the original channel count.
#' Explicit starts retain their original coefficient-coordinate meaning;
#' scalar `1` is expanded in normalized coordinates before summation.
#' Starting zeroes remain excluded in their original channels; compulsory
#' terms remain eligible. Coefficients always retain the original predictor
#' dimensions, with excluded rows zero. Shared observation weights remain
#' unchanged; offsets and coefficient prior centres are summed across groups.
#' This workflow requires weights shared across channels.
#'
#' With `cluster.refit=FALSE`, the returned model fits summed cluster outcomes
#' at the requested `lambda`. Its coefficients/predictions have one column per
#' cluster and its metrics describe those aggregated outcomes. This differs
#' from the permissive screening model retained inside an original-channel fit.
#' Aggregate likelihood criteria and original-channel criteria describe
#' different data and should be compared on their own response scales.
#'
#' `control.precluster` is a named list. Supported entries are `k.max=32L`,
#' `k=NULL`, `linkage="ward"`, `criterion="gic"`, `rowbins=1L`,
#' `size_fac=NULL` (use training column totals), `nthreads=NULL` (physical
#' cores), `screen.factor=0.1`, `refit.passes=2L` and `clusterids=NULL`.
#' Supplied cluster labels reuse a partition; responses are always aggregated
#' afresh. Clustering threads are serial when `parallelize=FALSE`.
#' Five-row binning is useful for ordered chromatographic scans; the general
#' default of one row avoids dependence on arbitrary observation order.
#'
#' Clustering uses training response profiles. Observation offsets are carried
#' into the aggregated and original regression stages. Validation uses the
#' original channels for original-channel refinement; aggregate-only validation
#' uses the same training-derived membership with summed responses/offsets.
#' [L0glm.tune()] supports original-channel refinement, learning a partition
#' once on the full data and separately in each CV training fold. Training
#' column totals determine size factors in that fold. Partitions are reused
#' across lambda values; selection coefficients are initialized independently.
#' Externally supplied `clusterids` or `size_fac` are reused as provided.
#' To tune an aggregate-only model, call [glm_aggregate()] first and tune its
#' returned responses directly.
#'
#' The workflow does not add a background basis or change compulsory terms.
#' Include a low-df B-spline basis through `always.include` for a smooth
#' background. Its nonnegative coefficients are refined separately from the
#' selected peaks on the original channels; see [compulsory_terms].
#'
#' The returned ordinary `L0glm_fit` gains `preclustering`: memberships,
#' channel names, candidates, screening coefficients, actual penalties,
#' initialization, clustering settings, cut history and phase timings.
#' One-channel original fits, zero-penalty original fits and starts excluding
#' all selectable rows use the ordinary fitter and record the bypass reason.
#'
#' @seealso [glm_cluster()], [glm_aggregate()], [L0glm.fit()], [L0nnpois()]
#' @examples
#' set.seed(8)
#' X <- diag(80)
#' B <- matrix(0,80,8)
#' B[c(12,35,64),] <- matrix(runif(24,2,12),3,8)
#' Y <- matrix(rpois(80*8, X %*% B + .1),80,8)
#' fit <- L0nnpois(X,Y,lambda="gic",offset=.1,precluster=TRUE,
#'   parallelize=FALSE,control.precluster=list(k.max=4L))
#' dim(coef(fit))
#' fit$preclustering$candidates
#' agg <- L0nnpois(X,Y,lambda="gic",offset=.1,precluster=TRUE,
#'   cluster.refit=FALSE,parallelize=FALSE,control.precluster=list(k=4L))
#' dim(coef(agg))
NULL

L0glm_precluster_flags<-function(precluster,cluster.refit) {
  for(value in list(precluster,cluster.refit)) if(!is.logical(value)||length(value)!=1L||is.na(value))
    stop("precluster and cluster.refit must be TRUE or FALSE.",call.=FALSE)
}

L0glm_precluster_coef<-function(value,name,X,added,mandatory,normalize,m) {
  if(is.null(value)) return(NULL)
  p<-ncol(X)
  if(!(is.numeric(value)||inherits(value,"dMatrix"))) stop(name," must contain numeric coefficients.",call.=FALSE)
  if(name=="start"&&is.null(dim(value))&&length(value)==1L&&isTRUE(value==1)) {
    prep<-L0glm_preprocessing(X,normalize,rep(0,p),mandatory)
    value<-L0glm_unit_start(value,X,prep,1L)
  } else if(length(value)==1L) value<-matrix(as.numeric(value),p,1L) else {
    if(is.null(dim(value))) value<-matrix(value,ncol=1L)
    if(added&&nrow(value)==p-1L) value<-rbind(0,value)
  }
  if(nrow(value)!=p||!ncol(value)%in%c(1L,m)||any(!is.finite(value)))
    stop(name," must match the original coefficient dimensions and be finite.",call.=FALSE)
  if(name=="start") value[value<0]<-0
  value
}

L0glm_precluster_sumcoef<-function(value,P) {
  if(is.null(value)) return(NULL)
  if(ncol(value)==1L) return(outer(as.numeric(value),as.numeric(Matrix::colSums(P))))
  as.matrix(value%*%P)
}

L0glm_precluster_validation<-function(validation,group=NULL,rows=NULL,added=FALSE,predictor_names=NULL) {
  if(is.null(validation)) return(NULL)
  if(is.null(group)) {
    if(inherits(validation,"L0glm_validation")) return(L0glm_subset_validation(validation,rows))
    if(!is.null(predictor_names)&&!is.null(colnames(validation$X))) {
      if(anyDuplicated(predictor_names)||anyDuplicated(colnames(validation$X))||
          !setequal(predictor_names,colnames(validation$X))) stop("Validation predictor names must match training predictors.",call.=FALSE)
      validation$X<-validation$X[,predictor_names,drop=FALSE]
    }
    if(added) validation$X<-cbind(Intercept=1,validation$X)
    validation$X<-validation$X[,rows,drop=FALSE]
    colnames(validation$X)<-NULL
    return(validation)
  }
  n<-nrow(validation$y);m<-length(group)
  agg<-glm_aggregate(validation$y,validation$weights,poisson(identity),group)
  P<-Matrix::sparseMatrix(i=seq_len(m),j=group,x=1,dims=c(m,max(group)))
  validation$y<-agg$Y_agg;validation$weights<-agg$weights_agg
  if(!is.null(validation$offset)) validation$offset<-as.matrix(L0glm_ext_matrix(validation$offset,n,m,0,"Validation offset")%*%P)
  if(inherits(validation,"L0glm_validation")) {
    validation$m<-max(group);validation$mass<-rowSums(validation$weights)
    validation$response<-list(y=validation$y,m=validation$m)
    for(name in c("baseline","refit_baseline")) if(!is.null(validation[[name]])) validation[[name]]<-as.matrix(validation[[name]]%*%P)
    validation$baseline_dispersion<-rep(1,validation$m)
  }
  validation
}

L0glm_precluster_fit<-function(X,y,options,fit.fun,supplied=character()) {
  begun<-proc.time()[[3L]];family<-options$family
  if(!identical(family$family,"poisson")||!identical(family$link,"identity")||
      !isTRUE(options$nonnegative))
    stop("precluster currently supports nonnegative identity-link Poisson fits. Use glm_cluster for other scalar families.",call.=FALSE)
  Y<-L0glm_cluster_matrix(y);n<-nrow(Y);m<-ncol(Y)
  if(nrow(X)!=n||any(Y<0)) stop("X and nonnegative Poisson responses must have compatible rows.",call.=FALSE)
  control<-options$control.precluster
  defaults<-list(k.max=32L,k=NULL,linkage="ward",criterion="gic",rowbins=1L,
    size_fac=NULL,nthreads=NULL,screen.factor=.1,refit.passes=if("group.maxit"%in%supplied) options$group.maxit else 2L,clusterids=NULL)
  if(!is.list(control)||length(control)&&(is.null(names(control))||any(!nzchar(names(control)))||
      anyDuplicated(names(control))||any(!names(control)%in%names(defaults))))
    stop("Use named, supported control.precluster entries; see ?preclustering.",call.=FALSE)
  control<-utils::modifyList(defaults,control,keep.null=TRUE)
  if(length(control$screen.factor)!=1L||!is.numeric(control$screen.factor)||!is.finite(control$screen.factor)||
      control$screen.factor<=0||control$screen.factor>1||!is.numeric(control$refit.passes)||
      length(control$refit.passes)!=1L||!is.finite(control$refit.passes)||control$refit.passes<1||
      control$refit.passes!=floor(control$refit.passes)) stop("Invalid screening factor or refinement pass count.",call.=FALSE)
  refit<-options$cluster.refit;added<-isTRUE(options$intercept)
  Z<-if(added) cbind(Intercept=1,X) else X;p<-ncol(Z)
  mandatory<-L0glm_mandatory(Z,if(added) c(1L,options$always.include+1L) else options$always.include)
  p_total<-p-length(mandatory)
  lambda<-L0glm_resolve_lambda(options$lambda,max(1L,p_total),n)
  phi<-L0glm_dispersion(options$dispersion,family,m)
  if(!is.null(options$dispersion)) options$dispersion<-phi
  original_start<-L0glm_precluster_coef(options$start,"start",Z,added,mandatory,options$normalize,m)
  prior<-L0glm_precluster_coef(options$prior.mean,"prior.mean",Z,added,mandatory,options$normalize,m)
  allowed<-if(is.null(original_start)) rep(TRUE,p) else as.numeric(Matrix::rowSums(original_start!=0))>0
  allowed[mandatory]<-TRUE
  call_fit<-function(design,response,opt) {
    opt$precluster<-FALSE
    do.call(fit.fun,c(list(X=design,y=response),opt[intersect(names(opt),names(formals(fit.fun)))]))
  }
  bypass<-if(m==1L) "one response channel" else if(refit&&lambda==0) "zero penalty" else
    if(refit&&!any(allowed[setdiff(seq_len(p),mandatory)])) "starting zeroes exclude all selectable rows" else NULL
  if(!is.null(bypass)) {
    fit<-call_fit(X,Y,options)
    fit$preclustering<-list(used=FALSE,reason=bypass,original_channels=m,seconds=c(total=proc.time()[[3L]]-begun))
    return(fit)
  }
  w<-L0glm_cluster_weights(options$weights,n,m)
  if(!any(w>0)) stop("Preclustering needs positive total observation weight.",call.=FALSE)
  if(length(w)>n) {
    W<-matrix(w,n,m)
    if(any(W!=W[,1L])) stop("Preclustering requires observation weights shared across channels.",call.=FALSE)
    w<-W[,1L]
  }
  threads<-if(isTRUE(options$parallelize)) L0glm_sis_threads(control$nthreads) else 1L
  cluster_start<-proc.time()[[3L]]
  if(is.null(control$clusterids)) {
    cl<-glm_cluster(Y,weights=w,family=family,method="poissondist",
      size_fac=if(is.null(control$size_fac)) colSums(Y) else control$size_fac,
      rowbins=control$rowbins,nthreads=threads,linkage=control$linkage,
      criterion=control$criterion,k=control$k,k.max=control$k.max)
  } else {
    labels<-control$clusterids
    if(!is.numeric(labels)||length(labels)!=m||any(!is.finite(labels)|labels<1|labels!=floor(labels)))
      stop("clusterids needs one positive integer per original channel.",call.=FALSE)
    cl<-list(clusterids=match(labels,sort(unique(labels))),history=NULL,linkage="supplied")
  }
  cluster_seconds<-proc.time()[[3L]]-cluster_start;group<-cl$clusterids;k<-max(group)
  P<-Matrix::sparseMatrix(i=seq_len(m),j=group,x=1,dims=c(m,k))
  response<-if(is.null(cl$Y_agg)) glm_aggregate(Y,options$weights,family,group)$Y_agg else cl$Y_agg
  opt<-options;opt$lambda<-if(refit) control$screen.factor*lambda else lambda
  opt$dispersion<-if(is.null(options$dispersion)) NULL else 1
  opt$weights<-if(is.null(options$weights)) NULL else w
  opt$offset<-if(is.null(options$offset)) NULL else as.matrix(L0glm_ext_matrix(options$offset,n,m,0,"offset")%*%P)
  opt$start<-L0glm_precluster_sumcoef(original_start,P);opt$prior.mean<-L0glm_precluster_sumcoef(prior,P)
  opt$validation<-if(refit) NULL else L0glm_precluster_validation(options$validation,group=group)
  selected_start<-proc.time()[[3L]];reduced<-call_fit(X,response,opt)
  selection_seconds<-proc.time()[[3L]]-selected_start
  C<-coef(reduced);candidates<-sort(unique(c(mandatory,which(Matrix::rowSums(abs(C))>0))))
  initialization<-NULL;refit_seconds<-0
  if(refit) {
    opt<-options;opt$lambda<-lambda;opt$intercept<-FALSE;opt$always.include<-match(mandatory,candidates)
    opt$group.maxit<-control$refit.passes;opt$screen<-"none"
    if(is.null(original_start)) {
      mean_start<-as.numeric(Matrix::rowSums(C[candidates,,drop=FALSE]))/m
      initialization<-list(method="aggregate_row_mean",rows=candidates,values=mean_start)
      opt$start<-matrix(mean_start,length(candidates),m)
    } else {
      opt$start<-original_start[candidates,,drop=FALSE]
      initialization<-list(method="supplied",rows=candidates,values=opt$start)
    }
    opt$prior.mean<-if(is.null(prior)) NULL else prior[candidates,,drop=FALSE]
    opt$validation<-L0glm_precluster_validation(options$validation,rows=candidates,added=added,predictor_names=colnames(X))
    if(!length(candidates)) {
      opt<-options;opt$lambda<-lambda;opt$start<-0;opt$screen<-"none"
      refit_start<-proc.time()[[3L]];fit<-call_fit(X,Y,opt)
    } else {
      refit_start<-proc.time()[[3L]];fit<-call_fit(Z[,candidates,drop=FALSE],Y,opt)
      entry<-summary(methods::as(methods::as(coef(fit),"generalMatrix"),"CsparseMatrix"))
      fit$coefficients<-Matrix::sparseMatrix(i=candidates[entry$i],j=entry$j,x=entry$x,
        dims=c(p,m),dimnames=list(colnames(Z),colnames(Y)))
      if(!is.null(fit$lambdas)) {
        penalties<-matrix(Inf,p,m);penalties[candidates,]<-as.numeric(fit$lambdas)
        fit$lambdas<-penalties
      }
      fit$X<-reduced$X;fit$p<-p_total;fit$always.include<-mandatory;fit$added_intercept<-added
      fit$preprocessing<-reduced$preprocessing
    }
    refit_seconds<-proc.time()[[3L]]-refit_start
  } else fit<-reduced
  # Refresh working quantities at the returned mean after all coefficient updates.
  mu<-fitted(fit);W<-as.matrix(fit$prior.weights)/pmax(mu,.Machine$double.eps)
  fit$weights<-if(ncol(mu)==1L) as.numeric(W) else Matrix::Matrix(W,sparse=TRUE)
  used_sizes<-if(is.null(cl$size_fac)) if(is.null(control$size_fac)) colSums(Y) else control$size_fac else cl$size_fac
  control$clusterids<-NULL;control$size_fac<-NULL
  fit$preclustering<-list(used=TRUE,refitted=refit,clusterids=group,
    original_channels=m,channel_names=colnames(Y),clusters=k,candidates=candidates,
    screen.coefficients=if(refit) C else NULL,lambda.screen=if(refit) control$screen.factor*lambda else NULL,
    lambda.target=lambda,initialization=initialization,history=cl$history,
    size_fac=used_sizes,screening=reduced$screening,
    control=control,nthreads=threads,
    seconds=c(clustering=cluster_seconds,selection=selection_seconds,refit=refit_seconds,total=proc.time()[[3L]]-begun))
  fit
}

# A tuning path shares a training-derived partition, never fitted coefficients.
L0glm_precluster_reuse <- function(args, fit) {
  if (isTRUE(fit$preclustering$used) && is.null(args$control.precluster$clusterids))
    args$control.precluster$clusterids <- fit$preclustering$clusterids
  args
}
