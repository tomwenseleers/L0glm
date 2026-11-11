# Add cached candidate solutions while retaining the ordinary global ISIS path.
# Coefficient fits and original-coordinate reconstruction use the same engine.
L0glm_isis_cache_size <- function(X,family,args,nonnegative,full,niter,nsis,mandatory,p_total,n_ic) {
  if(!inherits(X,'L0glm_data') || X$format!='bigmemory' || X$missing!='error' ||
     !identical(X$descriptor@description$sharedType,'FileBacked') ||
     niter<3L || full || nonnegative || family$family!='gaussian' || family$link!='identity' ||
      !is.null(args$prior.mean) || !is.null(args$validation) ||
     !is.null(args$maxit) && args$maxit<=1L ||
     !is.null(args$solver) && args$solver!='eigen' ||
     !is.null(args$method) && args$method!=1L ||
     !is.null(args$start) && (length(args$start)!=1L || args$start==0) ||
     p_total<1L || as.double(nrow(X))*ncol(X)<=X$control$max.cells || ncol(X)<=10*nrow(X))return(0L)
  charge <- presetLambda('gic',p_total,n_ic)
  if(!is.finite(charge)||charge<=0 || !isTRUE(all.equal(args$lambda,charge,tolerance=1e-12)))return(0L)
  limit <- X$control$working.columns
  if(is.null(limit))limit <- 10000L
  capacity <- floor(X$control$max.cells/nrow(X))-length(mandatory)-as.integer(args$intercept)
  size <- as.integer(max(0,min(limit,capacity,ncol(X))))
  if(size>nsis)size else 0L
}

L0glm_isis_cached <- function(X,y,family,nsis,niter,nthreads,size,state) {
  args <- state$args;added <- state$added;mandatory <- state$mandatory
  first <- sisGLM(X,y,family,weights=state$w,offset=state$off,
    nsis=size,return_stats=FALSE,nthreads=nthreads)
  pool <- sort(unique(c(mandatory,first$indices)))
  A <- L0glm_data_matrix(X[,pool,drop=FALSE])
  inner_args <- args;inner_args$always.include <- match(mandatory,pool)
  inner <- L0glm_isis(A,y,family,state$w,state$off,nsis,FALSE,FALSE,
    max(1L,niter%/%2L),inner_args,nthreads,p_reference=state$p_total,working=FALSE)
  candidates <- sort(unique(c(mandatory,head(first$indices,nsis))))
  result <- L0glm_isis(X,y,family,state$w,state$off,nsis,FALSE,FALSE,
    1L,args,nthreads,fixed_candidates=candidates,working=FALSE)
  history <- result$history;best_iteration <- 1L
  current <- result;previous <- NULL;global_scans <- 0L;converged <- FALSE
  for(iteration in seq_len(niter-1L)) {
    B <- coef(current$fit);active <- which(Matrix::rowSums(B!=0)>0)
    eta <- as.matrix(state$design[,active,drop=FALSE]%*%B[active,,drop=FALSE])+state$off
    global <- sisGLM(X,y,family,weights=state$w,offset=eta,nsis=nsis,
      return_stats=FALSE,nthreads=nthreads)
    global_scans <- global_scans+1L
    candidates <- sort(unique(c(mandatory,current$indices,global$indices)))
    if(identical(candidates,previous)){converged<-TRUE;break}
    proposal <- L0glm_isis(X,y,family,state$w,state$off,nsis,FALSE,FALSE,
      1L,args,nthreads,fixed_candidates=candidates,working=FALSE)
    history <- rbind(history,transform(proposal$history,iteration=nrow(history)+1L))
    if(proposal$history$gic[1L]<min(result$history$gic)) {
      result <- proposal;best_iteration <- nrow(history)
    }
    if(identical(coef(proposal$fit),coef(current$fit))){converged<-TRUE;break}
    current <- proposal;previous <- candidates
  }
  if(min(inner$history$gic)<min(result$history$gic)) {
    candidates <- pool[inner$candidates]
    result <- L0glm_isis(X,y,family,state$w,state$off,nsis,FALSE,FALSE,
      1L,args,nthreads,fixed_candidates=candidates,working=FALSE)
    history <- rbind(history,transform(result$history,iteration=nrow(history)+1L))
    best_iteration <- nrow(history)
  }
  retained <- result$candidates
  pruned <- L0glm_isis_gaussian_prune(result,X,y,family,state,nsis,nthreads)
  if(!is.null(pruned)) {
    result <- pruned
    history <- rbind(history,transform(pruned$history,iteration=nrow(history)+1L))
    best_iteration <- nrow(history)
  }
  retained <- sort(unique(c(result$candidates,retained)))
  extra <- setdiff(pool[inner$candidates],retained)
  capacity <- floor(X$control$max.cells/nrow(X))-added
  result$candidates <- sort(c(retained,head(extra,max(0,capacity-length(retained)))))
  result$history <- history;result$iterations <- nrow(history)
  result$best_iteration <- best_iteration;result$converged <- converged
  result$nthreads <- max(first$nthreads,inner$nthreads)
  result$working.set <- list(columns=length(pool),global.scans=global_scans+1L,
    pool.iterations=inner$iterations,pool.history=inner$history,pruned=!is.null(pruned))
  result$seconds <- proc.time()[['elapsed']]-state$started
  result
}

# Exact original-p GIC comparisons after unconstrained Gaussian refitting.
# Compulsory terms remain eligible, and positive supplied variance stays fixed.
L0glm_isis_gaussian_prune <- function(result,X,y,family,state,nsis,nthreads) {
  B <- coef(result$fit);rows <- sort(unique(c(state$included,
    which(Matrix::rowSums(B!=0)>0))))
  if(!length(setdiff(rows,state$included)))return(NULL)
  if(length(rows)>nsis+length(state$included))return(NULL)
  Z <- as.matrix(state$design[,rows,drop=FALSE]);z <- as.numeric(y-state$off)
  w <- as.numeric(state$w);positive <- sum(w>0)
  charge <- presetLambda('gic',state$p_total,state$n_ic)
  original <- rows
  repeat {
    R <- tryCatch(chol(crossprod(Z,w*Z)),error=function(e)NULL)
    if(is.null(R))return(NULL)
    b <- as.numeric(backsolve(R,forwardsolve(t(R),crossprod(Z,w*z))))
    rss <- sum(w*(z-as.numeric(Z%*%b))^2)
    delta <- b^2/diag(chol2inv(R))
    phi <- state$args$dispersion
    change <- if(is.null(phi)||phi<0)positive*log1p(delta/rss)-charge else delta/phi-charge
    change[rows%in%state$included | !is.finite(change)] <- Inf
    j <- which.min(change)
    if(!length(j) || change[j]>=-1e-10*max(1,abs(min(result$history$gic))))break
    rows <- rows[-j];Z <- Z[,-j,drop=FALSE]
    if(!length(rows)){b <- numeric();break}
  }
  if(identical(rows,original))return(NULL)
  start <- Matrix::sparseMatrix(i=rows,j=rep(1L,length(rows)),x=b,
    dims=dim(B),dimnames=dimnames(B))
  args <- state$args;args$start <- start
  survivors <- rows[rows>state$added]-state$added
  proposal <- L0glm_isis(X,y,family,state$w,state$off,nsis,FALSE,FALSE,
    1L,args,nthreads,fixed_candidates=survivors,working=FALSE)
  if(proposal$history$gic[1L]<min(result$history$gic))proposal else NULL
}
