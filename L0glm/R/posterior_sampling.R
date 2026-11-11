# Use explicit seeds locally, following the bootstrap's RNG-state convention.
L0glm_inference_seed <- function(seed) {
  if(is.null(seed))return(function()invisible(NULL))
  previous <- if(exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE))
    get('.Random.seed',envir=.GlobalEnv)else NULL
  set.seed(seed)
  function() {
    if(!is.null(previous))assign('.Random.seed',previous,envir=.GlobalEnv)else
      if(exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE))rm('.Random.seed',envir=.GlobalEnv)
    invisible(NULL)
  }
}

L0glm_posterior_draws <- function(cache,nsim,control) {
  if(length(nsim)!=1L || !is.numeric(nsim) || !is.finite(nsim) || nsim<1 ||
     nsim!=floor(nsim) || nsim>.Machine$integer.max)
    stop('nsim must be a positive integer.',call.=FALSE)
  k <- nrow(cache$map)
  if(as.double(k)*nsim>control$max.cells || prod(as.double(cache$dimensions))>.Machine$integer.max ||
     as.double(control$burnin)+nsim*control$thin>.Machine$integer.max)
    stop('Posterior draws exceed the memory or sampler iteration limit.',call.=FALSE)
  restore <- L0glm_inference_seed(control$seed);on.exit(restore(),add=TRUE)
  sample <- L0glm_posterior_sample(cache,cache$quadratic,nsim,control$burnin,control$thin)
  entries <- which(sample!=0,arr.ind=TRUE)
  draws <- Matrix::sparseMatrix(i=cache$index[entries[,1L]],j=entries[,2L],x=sample[entries],
    dims=c(prod(cache$dimensions),nsim))
  attr(draws,'coefficient.dim') <- cache$dimensions
  attr(draws,'coefficient.dimnames') <- cache$dimnames
  attr(draws,'coefficient.index') <- cache$index
  attr(draws,'seed') <- control$seed
  attr(draws,'sampler') <- attr(sample,'sampler')
  draws
}

L0glm_posterior_inference <- function(object,control) {
  context <- L0glm_posterior_setup(object,control)
  mode <- L0glm_posterior_mode(context,control)
  prepared <- context$prepared;B <- prepared$B;m <- ncol(B)
  position <- arrayInd(context$index,dim(B));row <- prepared$rows[position[,1L]]
  column <- position[,2L];dimensions <- dim(coef(object));dim_names <- dimnames(coef(object))
  index <- row+(column-1L)*dimensions[1L]
  cache <- list(map=context$map,constrained=context$constrained,quadratic=mode,index=index,
    dimensions=dimensions,dimnames=dim_names)
  draws <- L0glm_posterior_draws(cache,control$nsim,control)
  sample <- as.matrix(draws[index,,drop=FALSE]);k <- length(index)
  V <- as.matrix(context$map%*%mode$covariance%*%Matrix::t(context$map))
  center <- as.numeric(context$map%*%mode$mean)
  constrained <- any(context$constrained)
  covariance <- if(constrained && k)stats::cov(t(sample))else V
  estimate <- if(constrained)rowMeans(sample)else center
  se <- sqrt(pmax(0,diag(covariance)))
  probability <- c((1-control$level)/2,(1+control$level)/2)
  interval <- if(constrained && k)t(apply(sample,1L,stats::quantile,probs=probability,names=FALSE))else
    cbind(estimate+stats::qnorm(probability[1L])*se,estimate+stats::qnorm(probability[2L])*se)
  terms <- dim_names[[1L]];if(is.null(terms))terms <- paste0('x',seq_len(dimensions[1L]))
  missing <- is.na(terms)|!nzchar(terms);terms[missing] <- paste0('x',which(missing))
  channels <- dim_names[[2L]];if(is.null(channels))channels <- as.character(seq_len(m))
  table <- data.frame(term=terms[row],channel=channels[column],
    selected.estimate=as.numeric(coef(object)[cbind(row,column)]),
    mode=as.numeric(context$map%*%mode$mode),estimate,std.error=se,
    cred.low=interval[,1L],cred.high=interval[,2L],row,column)
  labels <- paste(table$term,table$channel,sep=':')
  dimnames(covariance) <- dimnames(V) <- list(labels,labels)
  structure(list(method='posterior',coefficients=table,draws=draws,covariance=covariance,
    gaussian.covariance=V,prior.precision=context$prior$precision[context$index],
    prior.mean=context$prior$center[context$index],prior.source=context$prior$source,
    level=control$level,nsim=control$nsim,
    seed=control$seed,converged=mode$converged,iterations=mode$iterations,score=mode$score,
    dispersion=object$dispersion,theta=object$theta,conditional.on.support=TRUE,
    constrained=constrained,sampler=attr(draws,'sampler'),
    control=control,sampling=cache),class='L0glm_posterior')
}

#' @rdname coefficient_inference
#' @param nsim Number of posterior coefficient draws. The simulation methods
#'   return coefficient draws, with excluded entries fixed at zero.
#' @param seed Optional random seed; preserves the caller's random-number state.
#' @param ... Reserved; unused arguments are rejected.
#' @method simulate L0glm_fit
#' @importFrom stats simulate
#' @export
simulate.L0glm_fit <- function(object,nsim=1L,seed=NULL,control.inference=list(),...) {
  if(length(list(...)))stop('Unused posterior simulation arguments.',call.=FALSE)
  control <- L0glm_inference_control(control.inference)
  control$seed <- seed;control$nsim <- nsim
  # Validate before preparing a selected design or dense information matrix.
  if(length(nsim)!=1L || !is.numeric(nsim) || !is.finite(nsim) || nsim<1 || nsim!=floor(nsim))
    stop('nsim must be a positive integer.',call.=FALSE)
  L0glm_inference_control(utils::modifyList(control,list(nsim=max(2,nsim))))
  context <- L0glm_posterior_setup(object,control)
  mode <- L0glm_posterior_mode(context,control)
  position <- arrayInd(context$index,dim(context$prepared$B))
  dimensions <- dim(coef(object))
  cache <- list(map=context$map,constrained=context$constrained,quadratic=mode,
    index=context$prepared$rows[position[,1L]]+(position[,2L]-1L)*dimensions[1L],
    dimensions=dimensions,dimnames=dimnames(coef(object)))
  L0glm_posterior_draws(cache,nsim,control)
}

#' @rdname coefficient_inference
#' @method simulate L0glm_posterior
#' @export
simulate.L0glm_posterior <- function(object,nsim=1L,seed=NULL,control.inference=list(),...) {
  if(length(list(...)))stop('Unused posterior simulation arguments.',call.=FALSE)
  control <- L0glm_inference_control(utils::modifyList(object$control,control.inference))
  control$seed <- seed
  L0glm_inference_control(control)
  L0glm_posterior_draws(object$sampling,nsim,control)
}

#' @rdname coefficient_inference
#' @param x A conditional posterior object.
#' @param digits Significant digits for printed posterior summaries.
#' @method print L0glm_posterior
#' @export
print.L0glm_posterior <- function(x,digits=4L,...) {
  cat('Conditional coefficient posterior\nPrior:',x$prior.source,'\n',
    x$nsim,'draws;',x$sampler,'\n')
  print(utils::head(x$coefficients,20L),digits=digits,row.names=FALSE)
  if(nrow(x$coefficients)>20L)cat('Showing 20 of',nrow(x$coefficients),'coefficients.\n')
  invisible(x)
}
