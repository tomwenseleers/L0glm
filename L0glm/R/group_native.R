# Built-in scalar GLMs share the native family and ridge implementations.
# Public preparation and custom/special-family callbacks retain the R driver.
L0glm_group_native_pass <- function(data,B,allowed,state,lambda,grouped,
                                   thresh,maxit,tol,args,relaxed) {
  native <- data$native.group
  if(is.null(native) || !is.null(args$validation))return(NULL)
  control <- args$control.solver
  auto_rowaug <- is.null(control$rowaug)
  if(is.null(control$ordering) && !is.null(control$NaturalOrdering))
    control$ordering <- as.integer(isTRUE(control$NaturalOrdering))
  threads <- if(is.null(args$group.nthreads))1L else min(data$m,args$group.nthreads)
  bridge <- L0glm_native_solver_bridge(args$solver,control)
  if(!is.null(bridge))threads <- 1L
  if(threads>1L) {
    previous_threads <- L0glm_set_thread_count(1L)
    on.exit(L0glm_restore_thread_count(previous_threads),add=TRUE)
  }
  context <- list(B=B,y=data$y,weights=as.matrix(data$weights),offset=as.matrix(data$offset),
    prior=data$prior,lower=data$lower,allowed=allowed,included=as.integer(data$included),
    phi=state$phi,estimated=state$estimated,stable=as.integer(state$stable),
    initial_mu=data$initial.mu,
    initial_eta=matrix(suppressWarnings(data$family$linkfun(data$initial.mu)),data$n,data$m),
    ids=as.integer(data$groups$id),multiplier=data$groups$multiplier,
    lambda=lambda,grouped=grouped,thresh=thresh,maxit=as.integer(maxit),tol=tol,
    theta_tol=data$control$theta.tol,relaxed=relaxed,loglik_kind=native$loglik,auto_rowaug=auto_rowaug)
  design <- data$Z
  if(inherits(design,'sparseMatrix'))design <- methods::as(methods::as(design,'generalMatrix'),'CsparseMatrix')
  else {design <- as.matrix(design);storage.mode(design) <- 'double'}
  result <- do.call(glm_predictor_group_cpp,c(list(design=design,context=context),native$descriptor,
    list(solver_control=control,solver_bridge=bridge,nthreads=as.integer(threads),cache=native$cache$pointer)))
  native$cache$pointer <- result$cache
  if(result$failed_iteration)warning('No finite descent step in glm IRLS iteration ',result$failed_iteration,'.',call.=FALSE)
  if(result$failed_solves)warning('Native group ridge solves retained previous coefficients in ',result$failed_solves,' channels.',call.=FALSE)
  if(result$inexact_solves)warning('Iterative solver did not reach its requested accuracy in ',result$inexact_solves,' native group solves.',call.=FALSE)
  if(result$solver_fallbacks)warning('Cholesky fallback was used in ',result$solver_fallbacks,' native group solves.',call.=FALSE)
  result[c('cache','failed_iteration','failed_solves','inexact_solves','solver_fallbacks')] <- NULL
  if(!nrow(result$history))result$history <- data.frame()
  result
}
