L0glm_local_pattern <- function(model) {
  S<-methods::as(methods::as(coef(model),"generalMatrix"),"CsparseMatrix")
  list(matrix=S,rows=S@i[S@x!=0]+1L)
}

L0glm_local_better <- function(value,current,tol) {
  if (!is.finite(value)) return(FALSE)
  if (is.infinite(current)) return(current>0)
  is.na(current) || value<current-tol*max(1,abs(current))
}

L0glm_local_eta <- function(B,d) {
  if(isTRUE(d$incremental) && !is.null(d$reference.B)) {
    if(identical(B,d$reference.B))return(d$reference.eta)
    delta<-B-d$reference.B
    rows<-which(rowSums(delta!=0)>0);cols<-which(colSums(delta!=0)>0)
    eta<-d$reference.eta
    if(length(rows)&&length(cols))eta[,cols]<-eta[,cols,drop=FALSE]+
      as.matrix(d$D[,rows,drop=FALSE]%*%delta[rows,cols,drop=FALSE])
    if(isTRUE(d$clip.eta) && all(B>=0))eta<-pmax(eta,0)
    return(eta)
  }
  selected <- which(rowSums(B!=0)>0)
  as.matrix(d$D[,selected,drop=FALSE]%*%Matrix::Matrix(B[selected,,drop=FALSE],sparse=TRUE))+d$O
}

L0glm_local_phi <- function(model,d) {
  if(!is.null(d$kind)) return(rep(1,d$m))
  phi <- rep_len(abs(model$dispersion),d$m)
  if(length(d$common.phi)) return(d$common.phi)
  if(!isTRUE(d$family$consistent.scale) && d$family$family=="gaussian" && any(d$estimated)) {
    prediction<-matrix(d$family$linkinv(L0glm_local_eta(as.matrix(coef(model)),d)),d$n,d$m)
    residual <- d$Y-prediction
    phi[d$estimated] <- colSums(d$W*residual^2)[d$estimated]/colSums(d$W>0)[d$estimated]
  }
  pmax(phi,1e-12)
}

L0glm_local_score <- function(model,d,control) {
  pattern<-L0glm_local_pattern(model);B<-as.matrix(pattern$matrix);eta <- L0glm_local_eta(B,d)
  ids<-pattern$rows[!pattern$rows%in%d$mandatory]
  rows<-length(unique(ids));entries<-length(ids)
  mandatory <- length(d$mandatory)*(d$m-as.integer(identical(d$kind,"multinomial")))
  nuisance <- if(identical(d$kind,"nb")) sum(model$theta.estimated) else
    if(is.null(d$kind)) sum(d$estimated) else 0L
  df <- entries+mandatory+nuisance
  quasi <- FALSE
  if(identical(d$kind,"multinomial")) {
    lp <- L0glm_ext_softmax(eta)$logp
    ll <- sum(d$W[,1L]*(lgamma(rowSums(d$Y)+1)-rowSums(lgamma(d$Y+1))+
      rowSums(ifelse(d$Y==0,0,d$Y*lp))))
    df <- if(d$grouped) rows*(d$m-1L)+mandatory else entries-rows+mandatory
  } else if(identical(d$kind,"cox")) {
    ll <- L0glm_ext_cox_eta(as.numeric(eta),d$response,d$W[,1L],d$cox.cache)$loglik
  } else {
    mu <- matrix(d$family$linkinv(eta),d$n,d$m)
    if(any(!is.finite(mu)) ||
       ((d$family$family%in%c("poisson","quasipoisson","Gamma","inverse.gaussian") ||
          identical(d$kind,"nb") || !is.null(d$family$scalar)) && any(mu<0)) ||
       (d$family$family%in%c("binomial","quasibinomial") && any(mu<0|mu>1)))
      return(list(value=Inf,loglik=-Inf,df=df,selected=rows,nonzero=entries,
        GIC=Inf,row.objective=Inf,objective=Inf,quasi=FALSE))
    positive <- d$family$family%in%c("poisson","quasipoisson","Gamma","inverse.gaussian") ||
      identical(d$kind,"nb")
    if(positive) mu <- pmax(mu,control$mu.floor)
    if(length(d$poisson.constant) && isTRUE(control$native.score)) {
      dev<-as.numeric(nuisance_poisson_deviance_cpp(d$poisson.Y,mu,d$W,
        as.integer(control$nthreads)))
      scale<-if(d$family$family=="quasipoisson") L0glm_local_phi(model,d) else rep(1,d$m)
      ll<--.5*sum((dev+d$poisson.constant)/scale)
      quasi<-d$family$family=="quasipoisson"
    } else if(identical(d$kind,"nb")) {
      ll <- sum(vapply(seq_len(d$m),function(j) {
        good <- d$W[,j]>0
        sum(d$W[good,j]*dnbinom(d$Y[good,j],mu=mu[good,j],size=model$theta[j],log=TRUE))
      },numeric(1)))
    } else {
      phi <- L0glm_local_phi(model,d)
      if(!length(d$common.phi)) phi[d$estimated] <- -phi[d$estimated]
      ll <- suppressWarnings(L0glm_loglik_mean(d$Y,mu,d$W,d$family,phi))
      if(!is.finite(ll) && startsWith(d$family$family,"quasi")) {
        ll <- -.5*sum(vapply(seq_len(d$m),function(j)
          sum(d$family$dev.resids(d$Y[,j],mu[,j],d$W[,j]))/abs(phi[j]),numeric(1)))
        quasi <- TRUE
      }
    }
  }
  if(is.na(ll)) ll <- -Inf
  units <- if(d$grouped && d$m>1L) rows else entries
  charge <- if(control$target=="IC") df else units
  criterion <- -2*ll+d$lambda*charge
  validation <- if(is.null(d$validation)) NULL else
    L0glm_validation_loss(d$validation,B,state=list(theta=model$theta))
  list(value=if(is.null(validation))criterion else validation,criterion=criterion,
    validation=validation,loglik=ll,df=df,selected=rows,nonzero=entries,
    GIC=-2*ll+d$gic*df,row.objective=-2*ll+d$lambda*rows,
    objective=-2*ll+d$lambda*units,quasi=quasi)
}

L0glm_local_remap <- function(model,rows,d,template) {
  E <- methods::as(methods::as(coef(model),"generalMatrix"),"TsparseMatrix")
  model$coefficients <- Matrix::sparseMatrix(i=rows[E@i+1L],j=E@j+1L,x=E@x,dims=c(d$p,d$m),
    dimnames=dimnames(template$coefficients))
  model$X <- d$D;model$p <- template$p;model$always.include <- d$mandatory
  model$offset <- template$offset
  if(is.null(d$kind))model$prior.weights <- template$prior.weights else model$weights <- template$weights
  model$added_intercept <- template$added_intercept;model$init.lambda <- template$init.lambda
  model$grouped <- d$grouped
  model$preprocessing <- template$preprocessing
  if(!is.null(model$lambdas)) {
    penalty <- matrix(Inf,d$p,d$m);penalty[rows,] <- as.numeric(model$lambdas)
    model$lambdas <- penalty
  }
  if(is.null(d$kind)) {
    model$dispersion.estimated <- d$estimated
    if(length(d$common.phi)) model$dispersion <- d$common.phi
  }
  model$screening <- template$screening
  model$preclustering <- template$preclustering
  if (!is.null(template$nuisance)) {
    model$nuisance <- template$nuisance
    model$nuisance$kkt <- NULL
    model$nuisance$local.search <- TRUE
  }
  if (!is.null(template$selected)) model$selected <- sort(unique(L0glm_local_pattern(model)$rows))
  model$metrics <- model$objective <- NULL
  for(name in c("dispersionCalibration","initialDispersion","initialTheta"))
    model[[name]] <- template[[name]]
  model
}

L0glm_local_refit <- function(mask,current,d,control) {
  mask <- mask & d$allowed;mask[d$mandatory,] <- TRUE
  if(identical(d$kind,"multinomial")) mask[rowSums(mask)<2L,] <- FALSE
  if(isTRUE(control$channel.refit) && d$m>1L &&
      (is.null(d$kind)||identical(d$kind,"nb"))) {
    previous<-as.matrix(coef(current))!=0;previous[d$mandatory,]<-TRUE
    changed<-which(colSums(mask!=previous)>0)
    if(length(changed)>0L && length(changed)<d$m)
      return(L0glm_local_refit_channels(mask,current,d,control,changed))
  }
  if(is.finite(control$block.max)) {
    partial<-L0glm_local_refit_block(mask,current,d,control)
    if(!is.null(partial))return(partial)
  }
  rows <- which(rowSums(mask)>0)
  A <- d$D[,rows,drop=FALSE];B <- as.matrix(coef(current))[rows,,drop=FALSE]
  included <- match(d$mandatory,rows);submask <- mask[rows,,drop=FALSE]
  if(identical(d$family$family,"gaussian") && d$family$link=="identity" && !d$nonnegative) {
    C <- matrix(0,length(rows),d$m)
    for(j in seq_len(d$m)) {
      keep <- which(submask[,j]);Z <- as.matrix(A[,keep,drop=FALSE])
      N <- sum(d$W[,j]>0)
      nuisance <- if(is.null(d$control.solver$.nuisance.df))0L else d$control.solver$.nuisance.df
      if(length(keep)>N || d$estimated[j] && length(keep)+nuisance>=N)
        stop('Selected support has insufficient information for an unpenalized refit.',call.=FALSE)
      if(length(keep)) {
        G <- crossprod(Z,d$W[,j]*Z);rhs <- crossprod(Z,d$W[,j]*(d$Y[,j]-d$O[,j]))
        R <- tryCatch(chol(G),error=function(e)NULL)
        if(is.null(R)||!L0glm_relaxed_gram(G,R))
          stop('Selected support has insufficient information for an unpenalized refit.',call.=FALSE)
        b <- as.numeric(backsolve(R,forwardsolve(t(R),rhs)))
        if(any(!is.finite(b))) stop("Nonfinite Gaussian refit.")
        C[keep,j] <- b
      }
    }
    model <- current;model$local.refit<-NULL;model$coefficients <- Matrix::Matrix(C,sparse=TRUE)
    model$dispersion <- pmax(colSums(d$W*(d$Y-as.matrix(A%*%C)-d$O)^2)/
      pmax(1,colSums(d$W>0)-colSums(submask)),1e-12)
    model$dispersion[!d$estimated] <- current$dispersion[!d$estimated]
    model$weights <- d$W;model$z <- d$Y-d$O
    model$lambdas <- matrix(ifelse(submask,0,Inf),length(rows),d$m)
    model$support.refit <- TRUE
    model$convergence <- data.frame(stage="local_refit",pass=1L,channel=seq_len(d$m),
      iterations=1L,converged=TRUE)
    return(L0glm_local_remap(model,rows,d,current))
  }
  # Zero starting entries lock coefficients out, so seed only eligible zeros.
  zero <- submask & B==0
  B[zero] <- if(d$nonnegative) 1e-6 else 1e-6
  B[!submask] <- 0
  if(identical(d$kind,"multinomial")) {
    for(j in seq_len(nrow(B))) {
      keep <- which(submask[j,]);B[j,keep] <- B[j,keep]-mean(B[j,keep])
    }
    # Explicitly impose CLR-entry masks on the existing exact joint IRLS engine.
    data <- L0glm_ext_prepare(A,d$Y,d$family,d$W,d$O,FALSE,included,
      FALSE,d$normalize,NULL,NULL,d$control.family)
    data$allowed <- submask
    start <- L0glm_ext_working_coefficients(data,B)
    args <- list(solver=d$solver,method=d$method,control.solver=d$control.solver,validation=NULL)
    result <- L0glm_ext_pass(data,start,submask,list(),0,FALSE,0,
      control$refit.maxit,control$refit.tol,args,TRUE)
    if(any(result$relaxation_skipped))
      stop('Selected support has insufficient information for an unpenalized refit.',call.=FALSE)
    model <- L0glm_ext_result(data,result,list(local_refit=result),0,d$grouped,TRUE,args)
  } else {
    supported <- if(identical(d$kind,'cox')) {
      data <- list(kind='cox',Z=A[d$cox.cache$order,,drop=FALSE],m=d$m,
        weights=d$W,offset=d$O,response=d$response,cox.cache=d$cox.cache)
      L0glm_relaxed_check(data,B,submask)
    } else L0glm_relaxed_standard(A,d$W,submask,
      if(identical(d$kind,'nb'))current$theta.estimated else d$estimated,
      if(is.null(d$control.solver$.nuisance.df))0L else d$control.solver$.nuisance.df)
    if(any(!supported))
      stop('Selected support has insufficient information for an unpenalized refit.',call.=FALSE)
    args <- list(X=A,y=d$original.y,family=d$family,weights=d$W,offset=d$O,
      intercept=FALSE,always.include=included,lambda=0,start=B,
      nonnegative=d$nonnegative,normalize=d$normalize,solver=d$solver,method=d$method,
      parallelize=FALSE,ncores=control$ncores,nthreads=1L,openmp=TRUE,
      maxit=control$refit.maxit,tol=control$refit.tol,thresh=0,
      grouped=d$grouped,group.maxit=1L,screen="none",local.search=FALSE,
      control.solver=d$control.solver,control.family=d$control.family)
    if(d$nonnegative && identical(d$family$family,'poisson') && identical(d$family$link,'identity') &&
       identical(d$solver,'eigen') && isTRUE(d$method==1L))args$control.solver$.refit.only <- TRUE
    # The constrained native backend uses ncores for its outcome threads.
    args$parallelize <- control$ncores>1L
    if(is.null(d$kind)) {
      phi <- L0glm_local_phi(current,d)
      if(!length(d$common.phi)) phi[d$estimated] <- -phi[d$estimated]
      args$dispersion <- phi
    } else if(identical(d$kind,"nb")) args$control.solver$.theta.start <- current$theta
    model <- do.call(L0glm_internal_fit,c(args,list(relaxed=TRUE)))
  }
  model$support.refit <- TRUE
  L0glm_local_remap(model,rows,d,current)
}

L0glm_local_refine <- function(fit,args,control=list(),support.only=FALSE) {
  clock <- proc.time()[[3L]]
  control <- L0glm_local_control(control,args)
  family <- fit$family;kind <- L0glm_ext_kind(family)
  B <- as.matrix(coef(fit));p <- nrow(B);m <- ncol(B);n <- nrow(fit$X)
  added <- as.integer(isTRUE(fit$added_intercept))
  mandatory <- fit$always.include
  allowed <- matrix(TRUE,p,m)
  if(!is.null(args$start)) {
    S <- args$start
    if(length(S)==1L) allowed[] <- S!=0 else {
      S <- as.matrix(S)
      if(nrow(S)==p-added && added) S <- rbind(0,S)
      if(ncol(S)==1L && m>1L) S <- matrix(S,nrow(S),m)
      if(!identical(dim(S),dim(B))) stop("Starting mask dimensions differ.")
      allowed <- S!=0
      if(identical(kind,"multinomial")) allowed <- matrix(rowSums(allowed)>0,p,m)
    }
  }
  allowed[mandatory,] <- TRUE
  if(any(B[!allowed]!=0)) stop("Initial fit violates an explicit starting mask.")
  estimated <- if(is.null(kind)) rep_len(fit$dispersion.estimated,m) else rep(FALSE,m)
  quasi <- startsWith(family$family,"quasi")
  calibrated <- identical(fit$dispersionCalibration$method,"cross-channel") &&
    isTRUE(fit$dispersionCalibration$applied)
  W <- L0glm_ext_matrix(if(is.null(kind))fit$prior.weights else fit$weights,n,
    if(!is.null(kind) && kind%in%c("cox","multinomial"))1L else m,1,"weights")
  O <- L0glm_ext_matrix(fit$offset,n,m,0,"offset")
  pool <- which(rowSums(allowed)>0)
  scope <- if(control$scope=="auto") if(length(fit$screening$candidates)) "screened" else "all" else control$scope
  if(scope=="screened" && length(fit$screening$candidates))
    pool <- intersect(pool,unique(c(mandatory,which(rowSums(B!=0)>0),fit$screening$candidates+added)))
  d <- list(D=fit$X,X=args$X,Y=fit$y,original.y=fit$y,O=O,W=W,family=family,
    kind=kind,p=p,n=n,m=m,added=added,mandatory=mandatory,allowed=allowed,pool=setdiff(pool,mandatory),
    nonnegative=isTRUE(args$nonnegative),grouped=isTRUE(fit$grouped),estimated=estimated,
    common.phi=if(is.null(kind) && (quasi||calibrated)) pmax(abs(fit$dispersion),1e-8) else numeric(),
    normalize=if(is.null(args$normalize)) "auto" else args$normalize,
    solver=if(is.null(args$solver)) "eigen" else args$solver,
    method=if(is.null(args$method)) 1L else args$method,
    refit.maxit=if(is.null(args$maxit))25L else args$maxit,
    refit.tol=if(is.null(args$tol))1e-10 else args$tol,
    control.solver=if(is.null(args$control.solver))controlSolverGen(
      solver=if(is.null(args$solver))"eigen" else args$solver,
      method=if(is.null(args$method))1L else args$method) else args$control.solver,
    control.family=if(is.null(args$control.family))list() else args$control.family,
    response=fit$response,cox.cache=fit$cox.cache,
    gic=presetLambda("gic",max(1L,fit$p),max(2,stats::nobs(fit))))
  if(calibrated)d$control.solver$.dispersion.fixed <- TRUE
  d$validation <- if(is.null(args$validation)) NULL else
    if(inherits(args$validation,"L0glm_validation")) args$validation else
      L0glm_validation_context(args$validation,args$X,fit$y,family,isTRUE(fit$added_intercept))
  d$lambda <- if(is.null(control$criterion)) fit$init.lambda else
    L0glm_resolve_lambda(control$criterion,max(1L,fit$p),max(2,stats::nobs(fit)))
  d$clip.eta<-d$nonnegative && identical(family$link,"identity") && all(O>=0) &&
    all(if(inherits(d$D,"sparseMatrix"))d$D@x>=0 else d$D>=0)
  d$incremental<-isTRUE(control$incremental) &&
    (is.null(family$scalar)||!identical(family$link,"identity")||d$clip.eta)
  d$in.pool<-logical(p);d$in.pool[d$pool]<-TRUE
  if(is.null(kind) && family$family%in%c("poisson","quasipoisson")) {
    d$poisson.Y<-if(family$family=="poisson")round(d$Y)else d$Y
    storage.mode(d$poisson.Y)<-"double"
    Y<-d$poisson.Y
    constant<-ifelse(Y==0,0,Y*log(pmax(Y,1e-300)))-Y-lgamma(Y+1)
    d$poisson.constant<--2*colSums(W*constant)
  }
  d$reference.eta<-L0glm_local_eta(B,d);d$reference.B<-B
  if(support.only)return(L0glm_local_support_refit(fit,d,control))
  if(!is.null(fit[["group"]]))return(L0glm_group_refine(fit,args,d,control))
  score <- L0glm_local_score(fit,d,control)
  initial <- score;current <- fit;evaluations <- 0L;failures <- character()
  screen.seconds <- refit.seconds <- 0
  block.evaluations<-channel.evaluations<-cache.hits<-0L
  history <- list(data.frame(iteration=0L,move="initial",value=score$value,GIC=score$GIC,
    selected=score$selected,nonzero=score$nonzero,evaluations=0L))
  proposals <- list();stop.reason <- "iteration_limit"
  unbounded <- is.infinite(score$value) && score$value<0
  if(unbounded) stop.reason <- "unbounded_criterion"
  for(iteration in if(unbounded)integer()else seq_len(control$maxit)) {
    cache <- new.env(parent=emptyenv());best <- list(model=current,score=score,move="none")
    base <- as.matrix(coef(current))!=0;base[d$mandatory,] <- TRUE
    active <- which(rowSums(base)>0)
    evaluate <- function(mask,move) {
      if(evaluations>=control$max.evals) return(NULL)
      mask <- mask & allowed;mask[d$mandatory,] <- TRUE
      if(identical(kind,"multinomial"))mask[rowSums(mask)<2L,] <- FALSE
      signature <- which(mask)
      key <- paste0("s:",digest::digest(signature,algo="xxhash64"))
      while(exists(key,cache,inherits=FALSE)) {
        cached <- get(key,cache,inherits=FALSE)
        if(identical(cached$signature,signature)) {cache.hits<<-cache.hits+1L;return(cached)}
        key <- paste0(key,":collision")
      }
      evaluations <<- evaluations+1L;t <- proc.time()[[3L]]
      candidate <- tryCatch(withCallingHandlers(L0glm_local_refit(mask,current,d,control),
        warning=function(w) invokeRestart("muffleWarning")),error=function(e)e)
      refit.seconds <<- refit.seconds+proc.time()[[3L]]-t
      if(inherits(candidate,"error")) {
        if(control$trace) cat("Refit failed:",conditionMessage(candidate),"at",deparse(conditionCall(candidate)),"\n")
        failures <<- c(failures,conditionMessage(candidate));return(NULL)
      }
      if(length(candidate$local.refit$rows))block.evaluations<<-block.evaluations+1L
      if(length(candidate$local.refit$channels))channel.evaluations<<-channel.evaluations+1L
      value <- L0glm_local_score(candidate,d,control)
      trial <- list(model=candidate,score=value,move=move,signature=signature)
      assign(key,trial,cache)
      proposals[[length(proposals)+1L]] <<- data.frame(iteration,move,value=value$value,
        improves=is.finite(value$value) && value$value<score$value-control$tol,
        selected=value$selected,nonzero=value$nonzero)
      if(L0glm_local_better(value$value,best$score$value,control$tol)) best <<- trial
      trial
    }
    if(iteration==1L || control$refit.every.round) {
      evaluate(base,"refit")
      if(best$move=="refit") {
        current<-best$model;score<-best$score
        d$reference.eta<-L0glm_local_eta(as.matrix(coef(current)),d)
        d$reference.B<-as.matrix(coef(current))
        history[[length(history)+1L]]<-data.frame(iteration=iteration-1L,move="refit",
          value=score$value,GIC=score$GIC,selected=score$selected,nonzero=score$nonzero,evaluations)
        base<-as.matrix(coef(current))!=0;base[d$mandatory,]<-TRUE
        active<-which(rowSums(base)>0)
        best<-list(model=current,score=score,move="none")
      }
    }
    scan <- L0glm_local_scan(current,d,control,active);screen.seconds <- screen.seconds+scan$seconds
    add <- scan$indices
    drops <- head(L0glm_local_drop_rank(current,d,scan$context),control$drop)
    for(j in add) {mask<-base;mask[j,]<-allowed[j,];evaluate(mask,paste0("add:",j))}
    for(k in seq.int(2L,max(2L,min(control$max.swap,length(add))))) if(k<=length(add) && k<=control$max.swap) {
      mask<-base;mask[head(add,k),]<-allowed[head(add,k),];evaluate(mask,paste0("add_batch:",k))
    }
    for(j in drops) {
      mask<-base;mask[j,]<-FALSE
      dropfit <- evaluate(mask,paste0("drop:",j))
      replacement <- add
      if(control$conditional && !is.null(dropfit)) {
        scan2 <- L0glm_local_scan(dropfit$model,d,control,union(which(rowSums(mask)>0),j))
        screen.seconds <- screen.seconds+scan2$seconds
        replacement <- unique(c(scan2$indices,add))
      }
      for(k in head(replacement,control$top)) {
        swapped<-mask;swapped[k,]<-allowed[k,];evaluate(swapped,paste0("swap:",j,":",k))
      }
    }
    for(k in seq.int(2L,max(2L,min(control$max.swap,length(drops),length(add)))))
      if(k<=length(drops) && k<=length(add) && k<=control$max.swap) {
        mask<-base;mask[head(drops,k),]<-FALSE
        reduced <- evaluate(mask,paste0("drop_batch:",k));incoming<-head(add,k)
        if(control$conditional && !is.null(reduced)) {
          scan2<-L0glm_local_scan(reduced$model,d,control,union(which(rowSums(mask)>0),head(drops,k)))
          screen.seconds<-screen.seconds+scan2$seconds
          incoming<-head(scan2$indices,k)
        }
        if(length(incoming)==k) {mask[incoming,]<-allowed[incoming,];evaluate(mask,paste0("splice:",k))}
      }
    # Entry-level improvements preserve the existing within-row spectral sparsity.
    # In CLR space an isolated category cannot form a nonzero zero-sum row.
    if(control$entry.moves && m>1L) {
      selectable <- setdiff(which(rowSums(base)>0),mandatory)
      add.rows <- sort(unique(c(selectable,add)))
      if(length(add.rows)) {
        Z<-d$D[,add.rows,drop=FALSE]
        if(identical(kind,"multinomial")) {
          P<-scan$context$context$probability;mass<-d$W[,1L]*rowSums(d$Y)
          s<-(d$Y-rowSums(d$Y)*P)*d$W[,1L];h<-mass*P*(1-P)
        } else {s<-scan$context$context$score;h<-scan$context$context$weight}
        U<-as.matrix(Matrix::crossprod(Z,s));I<-as.matrix(Matrix::crossprod(Z^2,h))
        gain<-(if(d$nonnegative)pmax(U,0)else U)^2/pmax(I,1e-20)
        possible<-allowed[add.rows,,drop=FALSE]&!base[add.rows,,drop=FALSE]
        if(identical(kind,"multinomial"))possible[rowSums(base[add.rows,,drop=FALSE])<2L,]<-FALSE
        gain[!possible]<--Inf
        for(u in head(order(gain,decreasing=TRUE),control$top)) if(is.finite(gain[u])) {
          pos<-arrayInd(u,dim(gain));j<-add.rows[pos[1L]];k<-pos[2L]
          mask<-base;mask[j,k]<-TRUE;evaluate(mask,paste0("add_entry:",j,":",k))
        }
      }
      indices <- which(base & allowed,arr.ind=TRUE)
      indices <- indices[indices[,1L]%in%selectable,,drop=FALSE]
      if(nrow(indices)) {
        magnitude <- abs(as.matrix(coef(current))[indices])
        for(u in head(order(magnitude),control$drop)) {
          j<-indices[u,1L];k<-indices[u,2L];mask<-base;mask[j,k]<-FALSE
          if(!identical(kind,"multinomial")||sum(mask[j,])>=2L) evaluate(mask,paste0("drop_entry:",j,":",k))
        }
      }
    }
    if(best$move=="none") {stop.reason<-if(evaluations>=control$max.evals) "evaluation_limit" else "no_improving_tested_move";break}
    current<-best$model;score<-best$score
    d$reference.eta<-L0glm_local_eta(as.matrix(coef(current)),d)
    d$reference.B<-as.matrix(coef(current))
    history[[length(history)+1L]] <- data.frame(iteration,move=best$move,value=score$value,
      GIC=score$GIC,selected=score$selected,nonzero=score$nonzero,evaluations)
    if(control$trace) cat(iteration,best$move,format(score$value,digits=12),"\n")
    if(evaluations>=control$max.evals) {stop.reason<-"evaluation_limit";break}
  }
  if(!unbounded && is.finite(control$block.max) && control$final.refit) {
    final.control<-control;final.control$block.max<-Inf;final.control$channel.refit<-FALSE
    mask<-as.matrix(coef(current))!=0;mask[d$mandatory,]<-TRUE
    t<-proc.time()[[3L]]
    refined<-tryCatch(withCallingHandlers(L0glm_local_refit(mask,current,d,final.control),
      warning=function(w)invokeRestart("muffleWarning")),error=function(e)e)
    refit.seconds<-refit.seconds+proc.time()[[3L]]-t;evaluations<-evaluations+1L
    if(!inherits(refined,"error")) {
      next.score<-L0glm_local_score(refined,d,control)
      if(L0glm_local_better(next.score$value,score$value,control$tol)) {
        current<-refined;score<-next.score
        history[[length(history)+1L]]<-data.frame(iteration=length(history),move="final_refit",
          value=score$value,GIC=score$GIC,selected=score$selected,nonzero=score$nonzero,evaluations)
      }
    } else failures<-c(failures,conditionMessage(refined))
  }
  current$local.search <- list(initial=initial,final=score,history=do.call(rbind,history),
    proposals=if(length(proposals))do.call(rbind,proposals)else data.frame(),
    seconds=proc.time()[[3L]]-clock,screen.seconds=screen.seconds,refit.seconds=refit.seconds,
    evaluations=evaluations,failures=failures,scope=scope,pool=pool,
    target=if(is.null(d$validation))control$target else "validation",
    block.evaluations=block.evaluations,channel.evaluations=channel.evaluations,cache.hits=cache.hits,
    lambda=d$lambda,scoring.floor=control$mu.floor,common.dispersion=d$common.phi,
    stopping=stop.reason,control=control,after.screening=!is.null(fit$screening))
  if (!is.null(d$validation)) {
    current$validation <- fit$validation
    current$validation$deviance <- score$validation
    current$validation$local.search <- c(initial=initial$validation,final=score$validation)
  }
  if (identical(kind,"cox")) current$baseline <- L0glm_ext_cox_baseline(current)
  if(inherits(current,"L0glm_extended_fit")) {
    current$metrics <- l0glmMetrics(current,mu_floor=control$mu.floor)
    current$objective <- current$metrics$objective
  }
  stopifnot(all(diff(current$local.search$history$value)<0),all(as.matrix(coef(current))[!allowed]==0))
  current
}

# Above the combinatorial-search limit, compare just one selected-support refit.
# This retains useful spectral boundary corrections without add/drop/swap scans.
L0glm_local_support_refit <- function(fit,d,control) {
  started <- proc.time()[[3L]]
  control$block.max <- Inf;control$channel.refit <- FALSE
  control$mu.floor <- .Machine$double.eps
  control$refit.maxit <- d$refit.maxit;control$refit.tol <- d$refit.tol
  mask <- as.matrix(coef(fit))!=0
  if(d$grouped)mask <- matrix(rowSums(mask)>0,d$p,d$m)&d$allowed
  if(!is.null(fit[['group']]))mask <- fit$group.active
  mask[d$mandatory,] <- TRUE
  decorate <- function(model) {
    if(!is.null(fit[['group']])) {
      for(name in c('group','group.multiplier','group.allowed','group.cost.type'))model[[name]] <- fit[[name]]
      model$group.active <- L0glm_group_active(as.matrix(coef(model)),fit$group.allowed,d$mandatory,
        d$grouped,list(id=match(fit[['group']],unique(fit[['group']]))))
    }
    model
  }
  score <- function(model) {
    value <- L0glm_local_score(model,d,control)
    if(!is.null(fit[['group']]) && control$target=='L0') {
      model <- decorate(model)
      value$criterion <- l0glmMetrics(model,mu_floor=0)$group_objective
      value$value <- if(is.null(d$validation))value$criterion else value$validation
    }
    value
  }
  initial <- score(fit)
  candidate <- tryCatch(L0glm_local_refit(mask,fit,d,control),error=identity)
  final <- initial;accepted <- FALSE;failures <- character()
  if(inherits(candidate,'error'))failures <- conditionMessage(candidate)else {
    next.score <- score(candidate)
    accepted <- L0glm_local_better(next.score$value,initial$value,control$tol)
    if(accepted) {fit <- decorate(candidate);final <- next.score}
  }
  row <- function(s,move,iteration)data.frame(iteration,move,value=s$value,GIC=s$GIC,
    selected=s$selected,nonzero=s$nonzero,evaluations=as.integer(iteration>0L))
  history <- row(initial,'initial',0L)
  if(accepted)history <- rbind(history,row(final,'support_refit',1L))
  fit$local.search <- list(initial=initial,final=final,history=history,proposals=data.frame(),
    evaluations=1L,failures=failures,seconds=proc.time()[[3L]]-started,screen.seconds=0,
    refit.seconds=proc.time()[[3L]]-started,scoring.floor=control$mu.floor,lambda=d$lambda,pool=d$pool,
    stopping='support_refit_only',target=if(is.null(d$validation))control$target else 'validation',
    control=control,after.screening=!is.null(fit$screening))
  if(!is.null(d$validation))fit$validation$deviance <- final$validation
  if(inherits(fit,'L0glm_extended_fit')) {
    fit$metrics <- l0glmMetrics(fit,mu_floor=0);fit$objective <- fit$metrics$objective
  }
  fit
}
