# R boundary for native extended-family coordination and working solves.
L0glm_extended_native_spec <- function(data,args) {
  if(isFALSE(args$control.solver$.native)||data$kind!='nb')return(NULL)
  links <- c('identity','log','sqrt','inverse','logit','probit','cloglog','cauchit','1/mu^2')
  link <- match(data$family$link,links)-1L
  if(is.na(link)||!link%in%c(0:3,8L))return(NULL)
  reference <- stats::make.link(data$family$link)
  same <- vapply(c('linkfun','linkinv','mu.eta'),function(nm)
    identical(body(data$family[[nm]]),body(reference[[nm]])),TRUE)
  if(!all(same))return(NULL)
  list(family=5L,link=link,theta=1,scalar=-1L,parameter=1,censor=0)
}

L0glm_extended_native <- function(data,B,allowed,state,lambda,grouped,thresh,maxit,tol,args,
                                  refit=FALSE,D=NULL) {
  spec <- L0glm_extended_native_spec(data,args)
  if(is.null(spec)||is.null(D)&&!is.null(args$validation))return(NULL)
  control <- L0glm_solver_memory(args$control.solver)
  cores <- if(!is.null(args$ncores))args$ncores else if(!is.null(args$group.nthreads))args$group.nthreads else 1L
  cores <- max(1L,min(data$m,as.integer(cores)))
  bridge <- L0glm_native_solver_bridge(args$solver,control)
  if(!is.null(bridge))cores <- 1L
  if(cores>1L) {
    old <- L0glm_set_thread_count(1L)
    on.exit(L0glm_restore_thread_count(old),add=TRUE)
  }
  if(is.null(control$ordering))control$ordering <- as.integer(isTRUE(control$NaturalOrdering))
  cache <- data$native.cache;signature <- list(solver=args$solver,rowaug=control$rowaug)
  working <- if(!is.null(cache)&&identical(cache$Z,data$Z)&&identical(cache$signature,signature))cache$working else NULL
  if(is.null(working)) {
    working <- L0glm_native_design(data$Z,control,compress=args$solver=='eigen')
    if(!is.null(cache)){cache$Z <- data$Z;cache$signature <- signature;cache$working <- working}
  }
  control$rowaug <- working$rowaug
  groups <- data$groups
  if(is.null(groups))groups <- list(id=seq_len(data$p),multiplier=rep(1,data$p))
  numeric_matrix <- function(x) {
    x <- as.matrix(x);if(typeof(x)!='double')storage.mode(x) <- 'double';x
  }
  context <- list(kind=0L,B=B,y=numeric_matrix(data$y),weights=numeric_matrix(data$weights),offset=numeric_matrix(data$offset),
    prior=numeric_matrix(data$prior),lower=data$lower,allowed=allowed,included=as.integer(data$included),
    theta=state$theta,estimated=state$estimated,stable=if(is.null(state$stable))integer(data$m)else as.integer(state$stable),
    noise_df=if(is.null(data$noise.df))rep(-1,data$m)else as.numeric(data$noise.df),
    ids=as.integer(groups$id),multiplier=as.numeric(groups$multiplier),predictor_groups=!is.null(data$groups),
    theta_min=data$control$theta.min,theta_max=data$control$theta.max,
    fixed_noise=isTRUE(data$fixed.dispersion),max_working=as.integer(data$control$max.working),
    lambda=lambda,grouped=grouped,thresh=thresh,maxit=as.integer(maxit),tol=tol,
    theta_tol=data$control$theta.tol,refit=refit,pilot=!is.null(D),D=D)
  result <- do.call(glm_extended_cpp,c(list(design=working$X,context=context),spec,
    list(solver_control=control,solver_bridge=bridge,nthreads=cores)))
  dimnames(result$B) <- dimnames(B)
  if(!is.null(D))return(result$B)
  state$theta <- as.numeric(result$state$theta)
  if(!is.null(state$stable)||any(state$estimated))state$stable <- result$state$stable
  result$state <- state
  if(result$failed_iteration)warning('No finite descent step in ',data$kind,' IRLS iteration ',result$failed_iteration,'.',call.=FALSE)
  result$failed_iteration <- NULL
  result
}
