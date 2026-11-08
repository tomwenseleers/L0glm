L0glm_bootstrap_inference <- function(object,control) {
  if(is.null(object$specification) || is.null(object$engine))
    stop('Bootstrap inference requires a saved fitting specification.',call.=FALSE)
  if(isTRUE(object$specification$precluster) && isFALSE(object$specification$cluster.refit))
    stop('Bootstrap original-channel fits or aggregate responses before fitting.',call.=FALSE)
  n <- nrow(object$X);added <- as.integer(object$added_intercept)
  original <- if(is.null(object$interactions))ncol(object$X)-added else object$interactions$p.original
  if(as.double(n)*original>control$max.cells)
    stop('Bootstrap row resampling exceeds control.inference$max.cells.',call.=FALSE)
  X <- if(!is.null(object$interactions))object$X$data else if(added)
    object$X[,-1L,drop=FALSE]else object$X
  y <- object$y;options <- object$specification
  options <- utils::modifyList(options,control$refit.args,keep.null=TRUE)
  if(is.null(control$refit) && any(!names(options)%in%names(formals(get(object$engine)))))
    stop('refit.args contains an unknown fitting argument.',call.=FALSE)
  if(!is.null(control$seed)) {
    previous <- if(exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE))
      get('.Random.seed',envir=.GlobalEnv)else NULL
    on.exit(if(is.null(previous)) {
      if(exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE))rm('.Random.seed',envir=.GlobalEnv)
    } else assign('.Random.seed',previous,envir=.GlobalEnv),add=TRUE)
    set.seed(control$seed)
  }
  terms <- object$interactions$terms
  if(is.null(terms))terms <- list()
  term_keys <- vapply(terms,paste,'',collapse=':')
  snapshots <- vector('list',control$nboot);status <- rep('',control$nboot)
  seconds <- numeric(control$nboot);converged <- rep(NA,control$nboot);warnings <- rep('',control$nboot)
  for(b in seq_len(control$nboot)) {
    started <- proc.time()[['elapsed']]
    rows <- sample.int(n,n,replace=TRUE);args <- options
    for(name in c('weights','offset'))if(!is.null(args[[name]]) && length(args[[name]])>1L)
      args[[name]] <- if(is.null(dim(args[[name]])))args[[name]][rows]else args[[name]][rows,,drop=FALSE]
    if(!is.null(args$family))args$family <- unserialize(serialize(args$family,NULL))
    args$X <- X[rows,,drop=FALSE];args$y <- y[rows,,drop=FALSE]
    snapshots[b] <- list(tryCatch(withCallingHandlers({
      fit <- if(is.null(control$refit))do.call(object$engine,args)else control$refit(args)
      if(!inherits(fit,'L0glm_fit') && !is.null(fit$fit_best))fit <- fit$fit_best
      if(!inherits(fit,'L0glm_fit') || ncol(coef(fit))!=ncol(coef(object)) ||
         as.integer(fit$added_intercept)!=added)
        stop('The bootstrap refit must return a model with matching response and intercept formats.')
      products <- fit$interactions$terms
      keys <- vapply(products,paste,'',collapse=':')
      for(i in which(!keys%in%term_keys)) {
        terms[[length(terms)+1L]] <- products[[i]];term_keys <- c(term_keys,keys[i])
      }
      entries <- summary(methods::as(methods::as(coef(fit),'generalMatrix'),'CsparseMatrix'))
      index <- entries$i
      product <- which(index>original+added)
      index[product] <- original+added+match(keys[index[product]-original-added],term_keys)
      status[b] <- 'ok'
      convergence <- L0glm_summary_convergence(fit$convergence)
      if(!is.null(convergence) && !any(convergence$unavailable>0))
        converged[b] <- all(convergence$not_converged==0)
      list(row=index,column=entries$j,x=entries$x)
    },warning=function(w) {
      warnings[b] <<- paste(c(warnings[b][nzchar(warnings[b])],conditionMessage(w)),collapse='; ')
      invokeRestart('muffleWarning')
    }),error=function(e) {status[b] <<- conditionMessage(e);NULL}))
    seconds[b] <- proc.time()[['elapsed']]-started
  }
  success <- which(status=='ok')
  if(length(success)<2L)stop('Fewer than two bootstrap fits succeeded; inspect the refit specification.',call.=FALSE)
  if(length(success)<control$nboot)
    warning(sprintf('%d of %d bootstrap fits failed; intervals use successful fits. See inference$diagnostics.',
      control$nboot-length(success),control$nboot),call.=FALSE)
  p <- original+added+length(terms);m <- ncol(coef(object))
  if(as.double(p)*m>.Machine$integer.max)stop('Bootstrap coefficient index exceeds sparse-matrix limits.',call.=FALSE)
  values <- lapply(snapshots[success],function(s)s$row+(s$column-1L)*p)
  draws <- Matrix::sparseMatrix(i=as.integer(unlist(values)),j=rep(seq_along(success),lengths(values)),
    x=as.numeric(unlist(lapply(snapshots[success],`[[`,'x'))),dims=c(p*m,length(success)))
  point <- summary(methods::as(methods::as(coef(object),'generalMatrix'),'CsparseMatrix'))
  point_index <- point$i+(point$j-1L)*p
  mandatory <- as.vector(outer(object$always.include,(seq_len(m)-1L)*p,'+'))
  index <- sort(unique(c(unlist(values),point_index,mandatory)))
  if(as.double(length(index))*length(success)>control$max.cells)
    stop('Bootstrap summaries exceed control.inference$max.cells.',call.=FALSE)
  samples <- as.matrix(draws[index,,drop=FALSE])
  row <- (index-1L)%%p+1L;column <- (index-1L)%/%p+1L
  estimate <- numeric(length(index));estimate[match(point_index,index)] <- point$x
  interval <- if(length(index))t(apply(samples,1L,stats::quantile,
    probs=c((1-control$level)/2,(1+control$level)/2),names=FALSE))else matrix(numeric(),0L,2L)
  probability <- rowMeans(samples!=0)
  row_names <- if(!is.null(object$interactions))object$X$columns else colnames(X)
  if(is.null(row_names))row_names <- paste0('x',seq_len(original))
  if(added)row_names <- c('Intercept',row_names)
  if(length(row_names)<p) {
    original_names <- row_names[seq_len(original)+added]
    row_names <- c(row_names[seq_len(original+added)],vapply(terms,function(j)
      paste(original_names[j],collapse=':'),''))
  }
  missing <- is.na(row_names)|!nzchar(row_names);row_names[missing] <- paste0('x',which(missing)-added)
  channel_names <- colnames(coef(object));if(is.null(channel_names))channel_names <- as.character(seq_len(m))
  table <- data.frame(term=row_names[row],channel=channel_names[column],estimate,
    bootstrap.mean=rowMeans(samples),std.error=apply(samples,1L,stats::sd),
    selection.frequency=probability,conf.low=interval[,1L],conf.high=interval[,2L],row,column)
  frequency <- Matrix::sparseMatrix(i=row,j=column,x=probability,dims=c(p,m),
    dimnames=list(row_names,channel_names))
  list(method='bootstrap',coefficients=table,selection.frequency=frequency,draws=draws,
    terms=terms,replicates=success,diagnostics=data.frame(replicate=seq_along(status),status,converged,warnings,seconds),
    nboot=control$nboot,successful=length(success),level=control$level,seed=control$seed,
    conditional.on.support=FALSE)
}
