# Offset screening is an outer support search; every candidate fit starts afresh.
L0glm_isis <- function(X,y,family,weights,offset,nsis,full,nonnegative,niter,options,nthreads=1L) {
  if (!is.numeric(niter) || length(niter)!=1L || !is.finite(niter) ||
      niter<1 || niter!=floor(niter)) stop("niter must be a positive integer.",call.=FALSE)
  if (length(options) && (is.null(names(options)) || any(!nzchar(names(options))) ||
      anyDuplicated(names(options)) || any(!names(options)%in%names(formals(L0glm.fit)))))
    stop("Use named L0glm.fit arguments for ISIS fitting options.",call.=FALSE)
  started <- proc.time()[["elapsed"]]
  nthreads <- L0glm_sis_threads(nthreads);threads_used <- 1L
  kind <- L0glm_ext_kind(family);n <- nrow(X);p <- ncol(X)
  m <- if (identical(kind,"cox")) 1L else if (identical(kind,"multinomial")) family$nlp+1L else NCOL(y)
  args <- utils::modifyList(list(lambda="gic",start=NULL,relaxed=TRUE,parallelize=FALSE,
    intercept=identical(kind,"multinomial") || identical(family$family,"binomial")),options,keep.null=TRUE)
  args$screen <- "none" # The outer search owns screening; candidate fits never recurse.
  if (!is.logical(args$intercept) || length(args$intercept)!=1L || is.na(args$intercept))
    stop("intercept must be TRUE or FALSE.",call.=FALSE)
  if (identical(kind,"cox") && args$intercept) {
    warning("A Cox intercept is not identifiable; intercept was omitted.",call.=FALSE)
    args$intercept <- FALSE
  }
  added <- as.integer(args$intercept)
  design <- if (added) cbind(Intercept=1,X) else X
  included <- L0glm_mandatory(design,if (added) c(1L,args$always.include+1L) else args$always.include)
  mandatory <- included[included>added]-added
  if (identical(kind,"cox") && length(L0glm_intercepts(design)))
    stop("Remove constant/intercept columns from a Cox design.",call.=FALSE)
  p_total <- ncol(design)-length(included)
  joint <- identical(kind,"cox") || identical(kind,"multinomial")
  w <- L0glm_ext_matrix(weights,n,if (joint) 1L else m,1,"weights")
  off <- L0glm_ext_matrix(offset,n,m,0,"offset")
  if (any(!is.finite(w)|w<0) || any(colSums(w)==0) || any(!is.finite(off)))
    stop("Each response needs positive finite weights; offsets must be finite.",call.=FALSE)
  n_ic <- if (identical(kind,"cox")) sum(w[,1L]*L0glm_ext_cox_response(y,n)$status) else
    if (is.null(kind)) n else sum(rowSums(w)>0)
  args$lambda <- L0glm_resolve_lambda(args$lambda,max(1L,p_total),max(2,n_ic))
  # Explicit coefficient arrays are in original coordinates; scalar one remains
  # a scalar so the selected fit can initialize in its normalized coordinates.
  expand <- function(value,name) {
    if (is.null(value) || name=="start" && length(value)==1L) return(value)
    if (length(value)==1L) return(matrix(value,ncol(design),m))
    B <- as.matrix(value)
    if (added && nrow(B)==p) B <- rbind(0,B)
    if (ncol(B)==1L && m>1L) B <- matrix(B,nrow(B),m)
    if (!identical(dim(B),c(ncol(design),as.integer(m))) || any(!is.finite(B)))
      stop(name," must have the original coefficient dimensions and finite values.",call.=FALSE)
    B
  }
  args$start <- expand(args$start,"start");args$prior.mean <- expand(args$prior.mean,"prior.mean")
  if (!is.null(args$start) && any(!is.finite(args$start))) stop("start must be finite.",call.=FALSE)
  allowed <- if (is.null(args$start)) rep(TRUE,p) else if (length(args$start)==1L)
    rep(args$start!=0,p) else rowSums(args$start[seq_len(p)+added,,drop=FALSE]!=0)>0
  allowed[mandatory] <- TRUE;eligible <- which(allowed)
  screen_X <- if (length(eligible)==p) X else X[,eligible,drop=FALSE]
  B <- Matrix::Matrix(0,ncol(design),m,sparse=TRUE)
  selected <- mandatory;previous <- NULL;best <- NULL;history <- list();converged <- FALSE
  for (iteration in seq_len(niter)) {
    if (length(eligible)) {
      active <- which(Matrix::rowSums(B!=0)>0)
      eta <- as.matrix(design[,active,drop=FALSE]%*%B[active,,drop=FALSE])+off
      screen <- sis_glm(screen_X,y,family,weights=w,offset=eta,nsis=min(nsis,length(eligible)),
        return_stats=full,nonnegative=nonnegative,nthreads=nthreads)
      threads_used <- max(threads_used,screen$nthreads)
      candidates <- sort(unique(c(selected,mandatory,eligible[screen$indices])))
    } else {screen <- NULL;candidates <- seq_len(p)}
    if (identical(candidates,previous)) {converged <- TRUE;break}
    rows <- if (added) c(1L,candidates+1L) else candidates
    fit_args <- args
    fit_args$always.include <- match(mandatory,candidates)
    for (name in c("start","prior.mean")) if (!is.null(args[[name]]) && length(args[[name]])>1L)
      fit_args[[name]] <- args[[name]][rows,,drop=FALSE]
    fit <- do.call(L0glm.fit,c(list(X=X[,candidates,drop=FALSE],y=y,family=family,
      weights=w,offset=off,nonnegative=nonnegative),fit_args))
    score <- L0glm.metrics(fit,mu_floor=1e-8)
    gic <- get.IC(score$loglik,max(2,n_ic),max(1L,p_total),score$df)$gic
    entry <- summary(methods::as(methods::as(coef(fit),"generalMatrix"),"CsparseMatrix"))
    B <- Matrix::sparseMatrix(i=rows[entry$i],j=entry$j,x=entry$x,dims=c(ncol(design),m),
      dimnames=list(colnames(design),colnames(coef(fit))))
    selected <- which(Matrix::rowSums(B[seq_len(p)+added,,drop=FALSE]!=0)>0)
    history[[iteration]] <- data.frame(iteration,candidates=length(candidates),selected=length(setdiff(selected,mandatory)),
      df=score$df,gic,deviance=score$deviance)
    if (is.null(best) || is.na(gic) && is.na(best$gic) ||
        !is.na(gic) && (is.na(best$gic) || gic<best$gic))
      best <- list(fit=fit,coefficients=B,candidates=candidates,rows=rows,screen=screen,gic=gic,iteration=iteration)
    previous <- candidates
    if (!length(eligible)) {converged <- TRUE;break}
  }
  fit <- best$fit;fit$coefficients <- best$coefficients;fit$X <- design
  fit$p <- p_total;fit$always.include <- included
  if (!is.null(fit$lambdas)) {
    penalties <- matrix(Inf,ncol(design),m)
    penalties[best$rows,] <- as.numeric(fit$lambdas)
    fit$lambdas <- if (is.null(dim(fit$lambdas))) as.numeric(penalties) else penalties
  }
  for (name in c("center","scale")) {
    value <- rep(if (name=="center") 0 else 1,ncol(design))
    value[best$rows] <- fit$preprocessing[[name]];fit$preprocessing[[name]] <- value
  }
  fit$preprocessing$p <- p_total
  fit$metrics <- L0glm.metrics(fit,mu_floor=1e-8)
  if (!is.null(fit$objective)) fit$objective <- fit$metrics$objective
  indices <- sort(unique(c(mandatory,which(Matrix::rowSums(fit$coefficients[seq_len(p)+added,,drop=FALSE]!=0)>0))))
  ranking <- stats <- NULL
  if (full && !is.null(best$screen)) {
    ranking <- eligible[best$screen$ranking]
    stats <- best$screen$stats;stats$feature <- eligible[stats$feature];stats$keep <- stats$feature%in%indices
  }
  list(indices=indices,top_statistics=NULL,ranking=ranking,stats=stats,nsis=length(indices),
    screen_size=nsis,family=family$family,link=family$link,theta=fit$theta,
    fit=fit,coefficients=fit$coefficients,candidates=best$candidates,lambda=args$lambda,
    history=do.call(rbind,history),best_iteration=best$iteration,iterations=length(history),
    converged=converged,seconds=proc.time()[["elapsed"]]-started,nthreads=threads_used)
}
