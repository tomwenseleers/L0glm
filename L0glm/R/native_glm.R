# Native scalar-family driver. R handles public input conventions and each
# built-in family initialization expression; numerical preparation and all
# grouped/relaxed passes share one C++ state machine.
L0glm_native_family <- function(family) {
  links <- c('identity','log','sqrt','inverse','logit','probit','cloglog','cauchit','1/mu^2')
  link <- match(family$link,links)-1L
  if(length(link)!=1L || is.na(link)) return(NULL)
  scalar <- family$scalar;theta <- 1
  if(!is.null(scalar)) {
    if(!scalar$kind %in% c(4L,5L,6L)) return(NULL)
    reference <- L0glm_scalar_family(scalar$kind,scalar$parameter,family$link,scalar$censor,family$family)
    descriptor <- list(family=0L,link=link,theta=1,scalar=scalar$kind,
      parameter=scalar$parameter,censor=scalar$censor)
  } else {
    name <- tolower(family$family)
    kind <- match(name,c('gaussian','poisson','binomial','gamma','inverse.gaussian'))-1L
    if(name=='quasipoisson') kind <- 1L
    if(name=='quasibinomial') kind <- 2L
    if(name=='quasi') kind <- match(family$varfun,c('constant','mu','mu(1-mu)','mu^2','mu^3'))-1L
    if(startsWith(name,'negative binomial(')) {
      kind <- 5L
      theta <- tryCatch(get('.Theta',environment(family$variance),inherits=FALSE),error=function(e)NA_real_)
      if(length(theta)!=1L || !is.finite(theta) || theta<=0)return(NULL)
    }
    if(length(kind)!=1L || is.na(kind)) return(NULL)
    reference <- tryCatch(if(name=='quasi') do.call(stats::quasi,list(link=family$link,variance=family$varfun)) else
      if(kind==5L) MASS::negative.binomial(theta,family$link) else
        get(if(name=='gamma')'Gamma'else name,asNamespace('stats'))(family$link),error=function(e)NULL)
    if(!is.null(family$score)) return(NULL)
    descriptor <- list(family=kind,link=link,theta=theta,scalar=-1L,parameter=1,censor=0)
  }
  callbacks <- c('variance','linkfun','linkinv','mu.eta','dev.resids','validmu','valideta','score','log_density')
  body_or_value <- function(x)if(is.function(x))body(x)else x
  if(is.null(reference) || any(vapply(callbacks,function(nm)
    !identical(body_or_value(family[[nm]]),body_or_value(reference[[nm]])),TRUE)) ||
    !identical(deparse(family$initialize),deparse(reference$initialize))) return(NULL)
  descriptor
}

L0glm_native_probes <- function(n) {
  had <- exists('.Random.seed',.GlobalEnv,inherits=FALSE)
  old <- if(had)get('.Random.seed',.GlobalEnv)else NULL
  on.exit(if(had)assign('.Random.seed',old,.GlobalEnv)else rm('.Random.seed',envir=.GlobalEnv))
  set.seed(1701)
  sample(c(-1,1),4L*n,replace=TRUE)
}

L0glm_native_context <- function(env,family) {
  fields <- c('y','weights','offset','eta','mu','prior.mean','scaling','lower','upper',
    'lambda','tol','thresh','maxit','offset_df','estimate_dispersion','scale_frozen',
    'hard_threshold','always.include','free_intercept')
  context <- mget(fields,env,inherits=FALSE)
  for(name in fields[seq_len(9L)]) context[[name]] <- as.numeric(context[[name]])
  context$auto_beta <- is.null(env$start);context$auto_phi <- is.null(env$dispersion)
  context$beta <- if(is.null(env$start)) numeric(ncol(env$X)) else as.numeric(env$start)
  context$dispersion <- if(is.null(env$dispersion)) NA_real_ else as.numeric(env$dispersion)
  context$gaussian_name <- family$family=='gaussian'
  context$binomial_name <- family$family %in% c('binomial','quasibinomial')
  context$signal_family <- L0glm_signal_family(family)
  context$reference_p <- env$preprocessing$p
  context$identity_boundary <- family$link=='identity' && family$family %in% c('gaussian','poisson','quasipoisson')
  eigen_chol <- env$solver=='eigen' && identical(env$control.solver$method,'chol') &&
    (is.null(env$control.solver$pivot) || env$control.solver$pivot==0L) &&
    (is.null(env$control.solver$ordering) || env$control.solver$ordering==0L)
  context$cache_gaussian <- eigen_chol && context$gaussian_name && family$link=='identity' &&
    env$preprocessing$mode!='none' && !inherits(env$X,'sparseMatrix') &&
    !isTRUE(env$control.solver$rowaug) && ncol(env$X)<=sum(context$weights>0) && ncol(env$X)<=2000L
  context$cache_fixed_weights <- eigen_chol && is.null(family$scalar) &&
    (family$family=='Gamma' && family$link=='log' || family$family %in% c('poisson','quasipoisson') && family$link=='sqrt')
  context$normalized <- env$preprocessing$mode!='none'
  context$always.include <- as.integer(context$always.include)
  context$free_intercept <- as.integer(context$free_intercept)
  context
}

L0glm_native_solver_bridge <- function(solver,control) {
  if(solver=='eigen')return(NULL)
  force(solver);force(control)
  function(Xw,yw,lambdas,prior.mean,start,lower,upper,iteration) {
    if(solver=='osqp')lmridge.osqp(Xw,as.numeric(yw),lambdas,lower,upper,as.numeric(start),
      control,as.numeric(prior.mean),iteration)else
      lmridge.glmnet(Xw,as.numeric(yw),lambdas,lower,upper,control,as.numeric(prior.mean),iteration)
  }
}

# Dependency wrappers use R's API. Parallelize their complete native channel
# passes in independent R workers; grouped support and dispersion remain in C++.
L0glm_native_workers <- function(design,descriptor,probes,control,bridge,threads,nthreads) {
  type <- if(Sys.info()[['sysname']]=='Windows')'PSOCK'else'FORK'
  cluster <- if(type=='PSOCK')parallel::makeCluster(threads,type=type,useXDR=FALSE,
    master='127.0.0.1',rscript_args='--vanilla')else parallel::makeCluster(threads,type=type)
  ready <- FALSE
  on.exit(if(!ready)parallel::stopCluster(cluster))
  init <- L0glm_worker_init;environment(init) <- baseenv()
  dll <- getLoadedDLLs()[['L0glm']]
  dev <- requireNamespace('pkgload',quietly=TRUE) && pkgload::is_dev_package('L0glm')
  parallel::clusterCall(cluster,init,.libPaths(),if(is.null(dll))NULL else dll[['path']],
    nthreads,getNamespaceInfo('L0glm','path'),dev)
  parallel::clusterCall(cluster,function(data){assign('.L0glm_native_worker',data,.GlobalEnv);NULL},
    list(design=design,descriptor=descriptor,probes=probes,control=control,bridge=bridge))
  run <- function(contexts) {
    chunks <- split(seq_along(contexts),cut(seq_along(contexts),threads,labels=FALSE))
    jobs <- lapply(chunks,function(i)contexts[i])
    results <- parallel::clusterApply(cluster,jobs,L0glm_native_worker_pass)
    messages <- unique(unlist(lapply(results,attr,'warnings'),use.names=FALSE))
    for(message in messages)warning(message,call.=FALSE)
    do.call(c,results)
  }
  ready <- TRUE
  list(run=run,stop=function()parallel::stopCluster(cluster))
}
L0glm_native_worker_pass <- function(contexts) {
  data <- get('.L0glm_native_worker',.GlobalEnv,inherits=FALSE)
  messages <- character()
  result <- withCallingHandlers(do.call(glm_native_pass_cpp,
    c(list(design=data$design,contexts=contexts),data$descriptor,
      list(probes=data$probes,solver_control=data$control,solver_bridge=data$bridge))),
    warning=function(w){messages<<-c(messages,conditionMessage(w));invokeRestart('muffleWarning')})
  attr(result,'warnings') <- unique(messages)
  result
}

L0glm_native_model <- function(frame,descriptor,ncores=1L) {
  take <- function(name)get(name,frame,inherits=TRUE)
  X <- take('X');Y <- take('y');p <- ncol(X);m <- ncol(Y);family <- take('family')
  control <- take('control.solver')
  if(take('solver')=='eigen' && identical(control$method,'chol') && is.null(control$ordering))
    control$ordering <- as.integer(isTRUE(control$NaturalOrdering))
  control$.native.context <- TRUE
  control$.native.centered <- L0glm_center_working_design(X,take('preprocessing'),take('lower'))
  control$.native.intercept <- intersect(take('always.include.start'),
    intersect(L0glm_unit_intercepts(control$.native.centered$X),which(take('lower')<0)))
  common <- list(X=X,family=family,lambda=take('lambda'),hard_threshold=TRUE,
    normalize=take('normalize'),preprocessing=take('preprocessing'),
    lower=take('lower'),upper=take('upper'),maxit=take('maxit'),thresh=take('thresh'),
    solver=take('solver'),method=take('method'),tol=take('tol'),
    always.include=take('always.include.start'),control.solver=control,verbose=take('verbose'))
  start <- take('start');allowed <- take('allowed');contexts <- vector('list',m)
  solver_control <- NULL
  for(col in seq_len(m)) {
    prior <- as.numeric(take('prior.mean')[,col])
    args <- c(common,list(y=Y[,col],weights=L0glm_observation_column(take('weights'),col),
      offset=L0glm_observation_column(take('offset'),col),prior.mean=prior,
      dispersion=take('disp_list')[[col]],start=if(is.null(start)||length(start)==1L)start else start[,col,drop=FALSE]))
    env <- do.call(L0glm.fit.internal,args)
    if(is.null(solver_control))solver_control <- env$control.solver
    a <- L0glm_native_context(env,family)
    a$allowed <- if(is.null(allowed))rep(TRUE,p)else as.logical(allowed[,col])
    a$center_shift <- as.numeric(env$centered$shift)
    a$center_intercept <- if(length(env$centered$intercept)) as.integer(env$centered$intercept[1L]-1L) else -1L
    a$original_prior <- prior
    contexts[[col]] <- a
  }
  design <- control$.native.centered$X
  if(inherits(design,'sparseMatrix')) design <- as(as(design,'generalMatrix'),'CsparseMatrix') else {
    design <- as.matrix(design);storage.mode(design) <- 'double'
  }
  threads <- as.integer(min(m,ncores))
  solver_control[c('.native.context','.native.centered','.native.intercept')] <- NULL
  bridge <- L0glm_native_solver_bridge(take('solver'),solver_control)
  if(threads>1L) {
    thread_state <- L0glm_set_thread_count(take('worker_nthreads'))
    on.exit(L0glm_restore_thread_count(thread_state),add=TRUE)
  }
  probes <- if(any(vapply(contexts,function(a)a$auto_phi,TRUE))) L0glm_native_probes(nrow(X)) else numeric()
  executor <- NULL
  if(!is.null(bridge) && threads>1L) {
    pool <- L0glm_native_workers(design,descriptor,probes,solver_control,bridge,threads,take('worker_nthreads'))
    on.exit(pool$stop(),add=TRUE)
    executor <- pool$run
  }
  result <- do.call(glm_native_cpp,c(list(design=design,contexts=contexts),descriptor,
    list(nthreads=threads,probes=probes,max_passes=if(m>1L && take('grouped'))as.integer(take('group.maxit'))else 1L,
      relaxed=take('relaxed'),nonnegative=take('nonnegative'),
      solver_control=solver_control,solver_bridge=bridge,executor=executor)))
  if(isTRUE(take('verbose')))print(result$convergence,row.names=FALSE)
  status <- result$native_stats;result$native_stats <- NULL
  if(status$failed_steps) warning('Native IRLS retained previous coefficients after line search in ',
    status$failed_steps,' channel passes.',call.=FALSE)
  if(status$failed_solves) warning('Native ridge solves retained previous coefficients in ',
    status$failed_solves,' channel passes.',call.=FALSE)
  if(status$inexact_solves)warning('Iterative solver did not converge to the requested accuracy in ',
    status$inexact_solves,' ridge solves; finite approximate solutions were retained.',call.=FALSE)
  if(status$solver_fallbacks)warning('The requested Eigen solver failed in ',status$solver_fallbacks,
    ' ridge solves; Cholesky LLT was used as fallback.',call.=FALSE)
  for(name in c('coefficients','lambdas','weights','z'))
    result[[name]] <- if(m==1L)as.numeric(result[[name]])else Matrix(result[[name]],sparse=TRUE)
  result
}
