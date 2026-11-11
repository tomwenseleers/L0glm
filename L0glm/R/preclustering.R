#' Preclustering response channels before sparse GLM fitting
#' @name preclustering
#' @md
#' @description
#' `l0glmFit(precluster=TRUE)` clusters and aggregates
#' response channels, screen predictors on the reduced responses, then refine
#' the sparse fit on every original channel. The identity-link workflow
#' supports Gaussian models with or without nonnegativity, and nonnegative
#' Poisson, quasipoisson, NB1, NB2, Gamma and inverse Gaussian models. Gamma
#' includes [gammaScale()] and [gammaCensored()]. `l0nnPois(precluster=TRUE)`
#' provides the same workflow for nonnegative identity-link Poisson models.
#'
#' @details
#' For conditionally independent Poisson channels with an identity link,
#' \deqn{Y^A_{ig}=\sum_{j\in G_g}Y_{ij}\sim
#' \mathrm{Poisson}(x_i^T\Gamma_g+o^A_{ig}),\qquad
#' \Gamma_g=\sum_{j\in G_g}\beta_j,\quad o^A_{ig}=\sum_{j\in G_g}o_{ij}.}
#' Nonnegativity preserves row support under summation. Gaussian channels
#' use arithmetic means and multiply each group's case weights by its size.
#' Signed effects can cancel within a Gaussian group, so group fits may miss
#' predictors with opposing channel effects. Sums preserve the count law for
#' independent Poisson and common-dispersion NB1, and the continuous law for
#' uncensored common-scale Gamma. In the other count and positive families
#' the summed mean is retained and the variance is estimated or approximated.
#' For supplied scales, typical channel intensity shares `r` are computed
#' from weighted column totals. The working scale of a sum uses
#' \eqn{\sum_j\phi_j r_j^q}, with \eqn{q=1,2,3} for quasipoisson, Gamma and
#' inverse Gaussian variance, respectively. NB2 uses the corresponding
#' moment-matched starting or fixed shape, taking a median across groups
#' because the aggregate family specifies one shape. Censored Gamma sums
#' approximate the observed censoring process.
#'
#' The cluster screen uses a permissive `screen.factor * lambda` (default `0.1`).
#' Original-channel refinement uses the requested penalty, calculated once
#' from the original predictor universe. Its default is two grouped passes;
#' an explicitly supplied `group.maxit` takes precedence unless
#' `control.precluster$refit.passes` is supplied. Automatic multichannel support
#' refinement accepts an unpenalized refit only when it improves the target.
#'
#' With `start=NULL`, Poisson refinement broadcasts each candidate's summed
#' aggregate coefficient divided by the original channel count. Other families
#' use their usual initialization on the reduced predictor set.
#' Explicit starts retain their original coefficient-coordinate meaning;
#' scalar `1` is expanded in normalized coordinates before summation.
#' Starting zeroes remain excluded in their original channels; compulsory
#' terms remain eligible. Coefficients always retain the original predictor
#' dimensions, with excluded rows zero. This workflow requires training
#' weights shared across channels. Offsets and coefficient prior centres use
#' the same aggregation as the response.
#'
#' With `cluster.refit=FALSE`, the returned model fits aggregated cluster outcomes
#' at the requested `lambda`. Its coefficients/predictions have one column per
#' cluster and its metrics describe those aggregated outcomes. This differs
#' from the permissive screening model retained inside an original-channel fit.
#' Aggregate likelihood criteria and original-channel criteria describe
#' different data and should be compared on their own response scales.
#'
#' `control.precluster` is a named list. Supported entries are `k.max=32L`,
#' `k=NULL`, `linkage="ward"`, `criterion="gic"`, `rowbins=1L`,
#' `size_fac=NULL` (use training column totals for count responses), `nthreads=NULL` (physical
#' cores), `screen.factor=0.1`, `refit.passes=2L` and `clusterids=NULL`.
#' Supplied cluster labels reuse a partition; responses are always aggregated
#' afresh. Clustering threads are serial when `parallelize=FALSE`.
#' Five-row binning is useful for ordered chromatographic scans; the general
#' default of one row avoids dependence on arbitrary observation order. Count
#' responses use size-adjusted Poisson clustering. Gaussian and positive
#' continuous responses use their own likelihood or deviance distances.
#' The specialized Gamma adapters use Ward clustering.
#'
#' Clustering uses training response profiles. Observation offsets are carried
#' into the aggregated and original regression stages. Validation uses the
#' original channels for original-channel refinement; aggregate-only validation
#' uses the same training-derived membership and aggregation.
#' [l0glmTune()] supports original-channel refinement, learning a partition
#' once on the full data and separately in each CV training fold. Training
#' column totals determine size factors in that fold. Partitions are reused
#' across lambda values; selection coefficients are initialized independently.
#' Externally supplied `clusterids` or `size_fac` are reused as provided.
#' To tune an aggregate-only model, call [glmAggregate()] first and tune its
#' returned responses directly.
#'
#' The workflow does not add a background basis or change compulsory terms.
#' Include a low-df B-spline basis through `always.include` for a smooth
#' background. Its nonnegative coefficients are refined separately from the
#' selected peaks on the original channels; see [compulsory_terms].
#'
#' The returned fit gains `preclustering`: memberships,
#' channel names, candidates, screening coefficients, actual penalties,
#' initialization, clustering settings, aggregate family and approximation,
#' cut history and phase timings.
#' One-channel original fits, zero-penalty original fits and starts excluding
#' all selectable rows use the ordinary fitter and record the bypass reason.
#'
#' @seealso [glmCluster()], [glmAggregate()], [l0glmFit()], [l0nnPois()]
#' @examples
#' set.seed(8)
#' X <- diag(80)
#' B <- matrix(0,80,8)
#' B[c(12,35,64),] <- matrix(runif(24,2,12),3,8)
#' Y <- matrix(rpois(80*8, X %*% B + .1),80,8)
#' fit <- l0nnPois(X,Y,lambda="gic",offset=.1,precluster=TRUE,
#'   parallelize=FALSE,control.precluster=list(k.max=4L))
#' dim(coef(fit))
#' fit$preclustering$candidates
#' agg <- l0nnPois(X,Y,lambda="gic",offset=.1,precluster=TRUE,
#'   cluster.refit=FALSE,parallelize=FALSE,control.precluster=list(k=4L))
#' dim(coef(agg))
NULL

L0glm_precluster_flags<-function(precluster,cluster.refit) {
  for(value in list(precluster,cluster.refit)) if(!is.logical(value)||length(value)!=1L||is.na(value))
    stop("precluster and cluster.refit must be TRUE or FALSE.",call.=FALSE)
}

L0glm_precluster_kind <- function(family, nonnegative) {
  kind <- tolower(family$family)
  if (startsWith(kind, "negative binomial")) kind <- "nb2"
  supported <- c("gaussian", "poisson", "quasipoisson", "nb1", "nb2",
    "gamma", "gamma_scale", "gamma_censored", "inverse.gaussian")
  if (!identical(family$link,"identity") || !kind %in% supported ||
      (kind != "gaussian" && !isTRUE(nonnegative)))
    stop("precluster requires identity-link Gaussian, or nonnegative identity-link Poisson, quasipoisson, negative binomial, Gamma or inverse Gaussian fits.",
      call. = FALSE)
  kind
}

L0glm_precluster_coef<-function(value,name,X,added,mandatory,normalize,m,nonnegative) {
  if(is.null(value)) return(NULL)
  p<-ncol(X)
  if(!(is.numeric(value)||inherits(value,"dMatrix"))) stop(name," must contain numeric coefficients.",call.=FALSE)
  if(name=="start"&&is.null(dim(value))&&length(value)==1L&&isTRUE(value==1)) {
    lower<-rep(if(nonnegative) 0 else -.Machine$double.xmax,p)
    prep<-L0glm_preprocessing(X,normalize,lower,mandatory)
    value<-L0glm_unit_start(value,X,prep,1L)
  } else if(length(value)==1L) value<-matrix(as.numeric(value),p,1L) else {
    if(is.null(dim(value))) value<-matrix(value,ncol=1L)
    if(added&&nrow(value)==p-1L) value<-rbind(0,value)
  }
  if(nrow(value)!=p||!ncol(value)%in%c(1L,m)||any(!is.finite(value)))
    stop(name," must match the original coefficient dimensions and be finite.",call.=FALSE)
  if(name=="start"&&nonnegative) value[value<0]<-0
  value
}

L0glm_precluster_sumcoef<-function(value,P) {
  if(is.null(value)) return(NULL)
  if(ncol(value)==1L) return(outer(as.numeric(value),as.numeric(Matrix::colSums(P))))
  as.matrix(value%*%P)
}

L0glm_precluster_validation<-function(validation,group=NULL,groups=NULL,rows=NULL,added=FALSE,
                                     predictor_names=NULL,family=NULL,dispersion=NULL) {
  if(is.null(validation)) return(NULL)
  if(is.null(groups)) {
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
  n<-nrow(validation$y);m<-nrow(groups);k<-ncol(groups)
  # Sums use the same case weights as the original responses.
  agg<-glmAggregate(validation$y,validation$weights,poisson(identity),group)
  average<-family$family=="gaussian"
  validation$y<-if(average) as.matrix(as.matrix(validation$y)%*%groups) else agg$Y_agg
  validation$weights<-if(average) {
    w<-if(is.null(agg$weights_agg)) matrix(1,n,k) else agg$weights_agg
    sweep(w,2L,tabulate(group,k),`*`)
  } else agg$weights_agg
  if(!is.null(validation$offset)) validation$offset<-L0glm_precluster_offset(validation$offset,n,m,groups)
  if(inherits(validation,"L0glm_validation")) {
    validation$m<-k;validation$mass<-rowSums(validation$weights)
    validation$family<-family
    validation$response<-if(identical(validation$kind,"nb"))
      L0glm_ext_nb_response(validation$y,family,n) else list(y=validation$y,m=k)
    for(name in c("baseline","refit_baseline")) if(!is.null(validation[[name]])) validation[[name]]<-as.matrix(validation[[name]]%*%groups)
    validation$baseline_dispersion<-rep_len(if(is.null(dispersion)) 1 else abs(dispersion),k)
    if(identical(validation$kind,"nb"))
      validation$baseline_state<-list(theta=validation$response$theta)
  }
  validation
}

# Match the variance of an independent sum at the typical channel proportions.
L0glm_precluster_noise <- function(kind, family, phi, Y, w, P, group) {
  fit <- list(family=family, dispersion=phi, approximation=NULL)
  if (kind %in% c("poisson", "nb1", "gamma_scale")) {
    if (kind == "gamma_scale" && family$scalar$censor > 0)
      fit$approximation <- "censored recorded intensities"
    return(fit)
  }
  degree <- switch(kind, quasipoisson=1, gamma=2, gamma_censored=2, nb2=2,
    inverse.gaussian=3, gaussian=0)
  k <- ncol(P)
  if (kind != "gaussian") fit$approximation <- "aggregate variance model"
  if (is.null(phi) && !kind %in% c("nb2", "gamma_censored")) return(fit)
  if (degree == 0) {
    share <- 1/tabulate(group,k)[group]
  } else {
    totals <- as.numeric(crossprod(rep_len(w,nrow(Y)),Y))
    mass <- as.numeric(totals %*% P)[group]
    share <- ifelse(mass > 0, totals/pmax(mass,.Machine$double.xmin),
      1/tabulate(group,k)[group])
  }
  group_sum <- function(value) as.numeric(value %*% P)
  if (kind == "nb2") {
    theta <- L0glm_cluster_family(family)$theta / stats::median(group_sum(share^2))
    fit$family <- if (identical(L0glm_ext_kind(family),"nb"))
      mgcv::nb(theta=if(isTRUE(family$n.theta>0)) -theta else theta, link="identity") else
      MASS::negative.binomial(theta,link="identity")
  } else if (kind == "gamma_censored") {
    scale <- family$scalar$parameter * stats::median(group_sum(share^2))
    fit$family <- gammaCensored(scale, link="identity", censor=family$scalar$censor)
  } else if (!is.null(phi)) {
    fit$dispersion <- group_sum(abs(phi) * share^max(1,degree))
    estimate <- group_sum(phi < 0) > 0
    fit$dispersion[estimate] <- -fit$dispersion[estimate]
  }
  fit
}

L0glm_precluster_fit<-function(X,y,options,fit.fun,supplied=character()) {
  begun<-proc.time()[[3L]];family<-options$family
  kind<-L0glm_precluster_kind(family,options$nonnegative)
  Y<-L0glm_cluster_matrix(y);n<-nrow(Y);m<-ncol(Y)
  if(nrow(X)!=n) stop("X and the responses must have compatible rows.",call.=FALSE)
  if(kind!="gaussian"&&any(Y<0)) stop("This family requires nonnegative responses.",call.=FALSE)
  if(kind%in%c("gamma","inverse.gaussian")&&any(Y<=0))
    stop("This family requires positive responses.",call.=FALSE)
  if(is.function(family$validate_response)&&!family$validate_response(Y))
    stop("Responses lie outside the family domain.",call.=FALSE)
  control<-options$control.precluster
  size_adjusted<-kind%in%c("poisson","quasipoisson","nb1","nb2")
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
  original_start<-L0glm_precluster_coef(options$start,"start",Z,added,mandatory,options$normalize,m,options$nonnegative)
  prior<-L0glm_precluster_coef(options$prior.mean,"prior.mean",Z,added,mandatory,options$normalize,m,options$nonnegative)
  allowed<-if(is.null(original_start)) rep(TRUE,p) else as.numeric(Matrix::rowSums(original_start!=0))>0
  allowed[mandatory]<-TRUE
  call_fit<-function(design,response,opt) {
    opt$precluster<-FALSE
    if (kind %in% c("gaussian","nb2")) opt$control.family$dispersion.calibration <- FALSE
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
  sizes<-if(size_adjusted) if(is.null(control$size_fac)) colSums(Y) else control$size_fac
  if(!size_adjusted&&!is.null(control$size_fac))
    stop("Preclustering size factors apply to count responses.",call.=FALSE)
  if(is.null(control$clusterids)) {
    cl<-glmCluster(Y,weights=w,family=if(size_adjusted) poisson(identity) else family,
      method=if(size_adjusted) "poissondist" else "deviance",size_fac=sizes,
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
  noise<-L0glm_precluster_noise(kind,family,if(is.null(options$dispersion)) NULL else phi,Y,w,P,group)
  sizes_by_group<-tabulate(group,k)
  if(kind=="gaussian") P<-P%*%Matrix::Diagonal(x=1/sizes_by_group)
  response<-if(size_adjusted&&!is.null(cl$Y_agg)) cl$Y_agg else as.matrix(Y%*%P)
  dimnames(response)<-list(rownames(Y),paste0("cluster_",seq_len(k)))
  opt<-options;opt$lambda<-if(refit) control$screen.factor*lambda else lambda
  opt$family<-noise$family;opt$dispersion<-noise$dispersion
  if(kind%in%c("poisson","nb1","nb2","gamma_scale","gamma_censored")&&!is.null(phi))
    opt$dispersion<-1
  opt$weights<-if(kind=="gaussian") outer(rep_len(w,n),sizes_by_group) else if(!is.null(options$weights)) w
  opt$offset<-if(is.null(options$offset)) NULL else L0glm_precluster_offset(options$offset,n,m,P)
  opt$start<-L0glm_precluster_sumcoef(original_start,P);opt$prior.mean<-L0glm_precluster_sumcoef(prior,P)
  opt$validation<-if(refit) NULL else L0glm_precluster_validation(options$validation,group=group,groups=P,
    family=opt$family,dispersion=opt$dispersion)
  selected_start<-proc.time()[[3L]];reduced<-call_fit(X,response,opt)
  selection_seconds<-proc.time()[[3L]]-selected_start
  C<-coef(reduced);candidates<-sort(unique(c(mandatory,which(Matrix::rowSums(abs(C))>0))))
  initialization<-NULL;refit_seconds<-0
  if(refit) {
    opt<-options;opt$lambda<-lambda;opt$intercept<-FALSE;opt$always.include<-match(mandatory,candidates)
    opt$group.maxit<-control$refit.passes;opt$screen<-"none"
    if(is.null(original_start)&&kind=="poisson") {
      mean_start<-as.numeric(Matrix::rowSums(C[candidates,,drop=FALSE]))/m
      initialization<-list(method="aggregate_row_mean",rows=candidates,values=mean_start)
      opt$start<-matrix(mean_start,length(candidates),m)
    } else if(!is.null(original_start)) {
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
      if(!is.null(fit$metrics)) {
        fit$metrics<-l0glmMetrics(fit,mu_floor=1e-8)
        fit$objective<-fit$metrics$objective
      }
    }
    refit_seconds<-proc.time()[[3L]]-refit_start
  } else fit<-reduced
  # Refresh working quantities at the returned mean after all coefficient updates.
  if(kind=="poisson") {
    mu<-fitted(fit);W<-L0glm_observation_broadcast(fit$prior.weights)/pmax(mu,.Machine$double.eps)
    fit$weights<-if(ncol(mu)==1L) as.numeric(W) else Matrix::Matrix(W,sparse=TRUE)
  }
  used_sizes<-if(is.null(cl$size_fac)) sizes else cl$size_fac
  control$clusterids<-NULL;control$size_fac<-NULL
  fit$preclustering<-list(used=TRUE,refitted=refit,clusterids=group,
    original_channels=m,channel_names=colnames(Y),clusters=k,candidates=candidates,
    screen.coefficients=if(refit) C else NULL,lambda.screen=if(refit) control$screen.factor*lambda else NULL,
    lambda.target=lambda,initialization=initialization,history=cl$history,
    size_fac=used_sizes,screening=reduced$screening,
    aggregate.family=noise$family$family,aggregate.approximation=noise$approximation,
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

L0glm_precluster_offset <- function(value,n,m,groups) {
  offset<-L0glm_observation_input(value,n,m,0,"offset")
  if(ncol(offset)==1L)as.numeric(offset) %o% as.numeric(Matrix::colSums(groups)) else
    as.matrix(offset%*%groups)
}
