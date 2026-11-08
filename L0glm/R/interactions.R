L0glm_interaction_order <- function(interactions) {
  if(length(interactions)!=1L || !is.numeric(interactions) || !is.finite(interactions) ||
     interactions<0 || interactions>=.Machine$integer.max || interactions!=floor(interactions))
    stop('interactions must be a nonnegative integer smaller than .Machine$integer.max.',call.=FALSE)
  as.integer(interactions)
}

L0glm_interaction_control <- function(control,args) {
  defaults <- list(max.terms=10000L,max.cells=5e7,tol=1e-8,criterion=NULL)
  if(!is.list(control) || length(control) && (is.null(names(control)) ||
     any(!names(control)%in%names(defaults)) || anyDuplicated(names(control))))
    stop('Use named control.interactions entries.',call.=FALSE)
  defaults[names(control)] <- control
  for(name in c('max.terms','max.cells','tol'))if(length(defaults[[name]])!=1L ||
     !is.numeric(defaults[[name]]) || !is.finite(defaults[[name]]) || defaults[[name]]<=0)
    stop('Interaction limits and tolerance must be positive finite values.',call.=FALSE)
  if(is.null(defaults$criterion))defaults$criterion <-
    if(is.character(args$lambda))args$lambda else 'gic'
  defaults$criterion <- L0glm_ic_name(defaults$criterion)
  defaults
}

L0glm_interaction_groups <- function(group,terms,base) {
  if(is.null(group))return(NULL)
  original <- unique(group[group!='0'])
  parents <- lapply(terms,function(j)sort(group[j]))
  key <- vapply(parents,function(g)paste(nchar(g),g,sep='=',collapse='|'),'')
  unique_key <- unique(key)
  labels <- paste0('interaction[',vapply(parents[match(unique_key,key)],paste,'',collapse=':'),']')
  for(i in seq_along(labels))while(labels[i]%in%c(original,head(labels,i-1L)))
    labels[i] <- paste0('.',labels[i])
  label <- c(group,labels[match(key,unique_key)])
  m <- ncol(coef(base));kind <- L0glm_ext_kind(base$family)
  cost <- L0glm_group_cost(label,m,base$grouped,kind,identical(base$group.cost.type,'dimensions'))
  cost[names(base$group.multiplier)] <- base$group.multiplier
  list(label=label,cost=cost)
}

L0glm_interaction_restore <- function(fit,X,columns,p,groups=NULL) {
  added <- as.integer(isTRUE(fit$added_intercept));X$intercept <- added==1L
  map <- if(added)c(1L,columns+1L)else columns
  m <- ncol(coef(fit));q <- ncol(X)
  old <- fit$coefficients
  entries <- summary(methods::as(methods::as(old,'generalMatrix'),'CsparseMatrix'))
  fit$coefficients <- Matrix::sparseMatrix(i=map[entries$i],j=entries$j,x=entries$x,
    dims=c(q,m),dimnames=list(colnames(X),colnames(old)))
  fit$X <- X;fit$p <- p;fit$always.include <- map[fit$always.include]
  for(name in c('center','scale')) {
    value <- rep(if(name=='center')0 else 1,q)
    value[map] <- fit$preprocessing[[name]];fit$preprocessing[[name]] <- value
  }
  fit$preprocessing$p <- q-added
  if(!is.null(fit$lambdas)) {
    value <- matrix(Inf,q,m);value[map,] <- as.numeric(fit$lambdas)
    fit$lambdas <- value
  }
  if(!is.null(groups)) {
    label <- if(added)c('0',groups$label)else groups$label
    label[fit$always.include] <- '0'
    fit[["group"]] <- label;fit$group.multiplier <- groups$cost
    attr(fit$group.multiplier,'dimension_based') <- fit$group.cost.type=='dimensions'
    for(name in c('group.allowed','group.active')) {
      value <- matrix(FALSE,q,m);value[map,] <- fit[[name]];fit[[name]] <- value
    }
  }
  if(!is.null(fit$selected))fit$selected <- map[fit$selected]
  if(!is.null(fit$local.search))fit$local.search$coordinate_rows <- map
  fit$metrics <- fit$objective <- NULL
  fit
}

L0glm_interaction_array <- function(value,rows,p,m,products,A,args,included,kind,name) {
  if(is.null(value) || name=='start' && length(value)==1L)return(value)
  added <- length(rows)>0L && rows[1L]==0L
  B <- L0glm_ext_matrix(value,p+as.integer(added),m,0,name)
  original <- rows+as.integer(added)
  out <- B[original,,drop=FALSE]
  if(length(products)) {
    extra <- matrix(0,length(products),m)
    if(name=='start') {
      lower <- rep(if(isTRUE(args$nonnegative))0 else -.Machine$double.xmax,ncol(A))
      family <- args$family
      if(!inherits(family,'family'))family <- family()
      if(identical(kind,'multinomial'))lower[] <- -.Machine$double.xmax
      if(family$family%in%c('binomial','quasibinomial'))
        lower[L0glm_intercepts(A)] <- -.Machine$double.xmax
      scale <- L0glm_preprocessing(A,if(is.null(args$normalize))'auto'else args$normalize,
        lower,included)$scale
      extra[] <- 1/tail(scale,length(products))
      eligible <- B!=0
      if(identical(kind,'multinomial'))eligible <- matrix(rowSums(eligible)>0,nrow(B),m)
      for(j in seq_along(products))extra[j,!apply(eligible[products[[j]]+as.integer(added),,drop=FALSE],2,all)] <- 0
      if(identical(kind,'multinomial')) {
        contrast <- c(1,rep(-1/(m-1L),m-1L))
        extra <- extra*rep(contrast,each=nrow(extra))
      }
    }
    out <- rbind(out,extra)
  }
  out
}

L0glm_interaction_fit <- function(args,interactions,control) {
  control <- L0glm_interaction_control(control,args)
  if(isTRUE(args$precluster) && isFALSE(args$cluster.refit))
    stop('Interaction generation requires cluster.refit=TRUE; aggregate responses before fitting an aggregate-only model.',call.=FALSE)
  args$interactions <- 0L;args$control.interactions <- NULL
  X <- args$X
  if(is.data.frame(X))X <- as.matrix(X)
  base <- do.call(l0glmFit,args)
  p <- ncol(X);n <- nrow(X);m <- ncol(coef(base));added <- as.integer(base$added_intercept)
  kind <- L0glm_ext_kind(base$family)
  mandatory <- base$always.include[base$always.include>added]-added
  if(is.null(args$lambda))args$lambda <- if(!is.null(base$screening))'gic'else base$init.lambda
  args$intercept <- base$added_intercept
  args$family <- if(is.null(args$family))gaussian()else args$family
  if(is.null(args$control.solver))args$control.solver <- controlSolverGen(
    if(is.null(args$solver))'eigen'else args$solver,if(is.null(args$method))1L else args$method)
  group <- if(is.null(args[["group"]]))NULL else L0glm_group_labels(args[["group"]],p)
  if(!is.null(group))group[mandatory] <- '0'
  products <- list();view <- L0glm_design_view(X)
  groups <- L0glm_interaction_groups(group,products,base)
  current <- L0glm_interaction_restore(base,view,seq_len(p),base$p,groups)
  models <- list(current);orders <- 0L;lambdas <- base$init.lambda
  validation <- if(is.null(args$validation))NULL else base$validation$deviance
  best_validation <- if(is.null(validation))Inf else validation
  for(degree in seq.int(2L,interactions+1L)) {
    selected <- setdiff(which(Matrix::rowSums(coef(current)[seq_len(p)+added,,drop=FALSE]!=0)>0),mandatory)
    keep <- sort(unique(c(mandatory,selected)))
    if(!is.null(group)) {
      selected_groups <- unique(group[selected])
      keep <- sort(unique(c(keep,which(group%in%selected_groups))))
    }
    if(length(selected)<degree)break
    count <- choose(length(selected),degree)
    if(!is.finite(count) || count+length(products)>control$max.terms)
      stop('Interaction candidates exceed control.interactions$max.terms.',call.=FALSE)
    if(as.double(n)*(length(keep)+count+length(products)+added)>control$max.cells)
      stop('Interaction working design exceeds control.interactions$max.cells.',call.=FALSE)
    products <- c(products,combn(selected,degree,simplify=FALSE))
    view <- L0glm_design_view(X,products)
    columns <- c(keep,p+seq_along(products))
    A <- view[,columns,drop=FALSE]
    options <- args;options$X <- A;options$screen <- 'none'
    options$always.include <- match(mandatory,keep)
    options$control.solver$.ic.p <- base$p+length(products)
    design <- if(added)cbind(Intercept=1,A)else A
    rows <- if(added)c(0L,keep)else keep
    for(name in c('start','prior.mean'))options[[name]] <- L0glm_interaction_array(
      args[[name]],rows,p,m,products,design,args,
      if(added)c(1L,options$always.include+1L)else options$always.include,kind,name)
    groups <- L0glm_interaction_groups(group,products,base)
    if(!is.null(groups)) {
      options[["group"]] <- groups$label[columns]
      labels <- unique(options[["group"]]);labels <- labels[labels!='0']
      options$group.multiplier <- groups$cost[labels]
      if(identical(base$group.cost.type,'dimensions'))attr(options$group.multiplier,'dimension_based') <- TRUE
      if(is.character(args$lambda) && is.null(options$control.local$target))options$control.local$target <- 'IC'
    }
    if(!is.null(args$validation)) {
      options$validation <- args$validation
      options$validation$X <- L0glm_design_view(args$validation$X,products)[,columns,drop=FALSE]
    }
    if(startsWith(base$family$family,'quasi'))options$dispersion <- base$dispersion
    fitted <- do.call(l0glmFit,options)
    if(startsWith(base$family$family,'quasi'))fitted$dispersion.estimated <- base$dispersion.estimated
    current <- L0glm_interaction_restore(fitted,view,columns,base$p+length(products),groups)
    models[[length(models)+1L]] <- current
    orders <- c(orders,degree-1L);lambdas <- c(lambdas,current$init.lambda)
    if(!is.null(validation)) {
      loss <- fitted$validation$deviance;validation <- c(validation,loss)
      if(loss>best_validation+control$tol*max(1,abs(best_validation)))break
      best_validation <- min(best_validation,loss)
    }
  }
  pool <- base$p+length(products)
  metrics <- lapply(models,function(model) {
    model$p <- pool
    l0glmMetrics(model,mu_floor=1e-8)
  })
  score <- if(is.null(validation))vapply(seq_along(models),function(i) {
    value <- metrics[[i]][[control$criterion]]
    if(is.na(value) && identical(base$family$family,'quasibinomial')) {
      model <- models[[i]];mu <- predict(model)
      loss <- sum(vapply(seq_len(m),function(j)sum(base$family$dev.resids(model$y[,j],mu[,j],
        L0glm_observation_column(model$prior.weights,j)))/base$dispersion[j],1))
      value <- get.IC(-loss/2,n,max(1,pool),metrics[[i]]$df)[[control$criterion]]
    }
    value
  },1)else validation
  best <- which.min(score)
  if(!length(best) || all(is.na(score)|score==Inf))
    stop('No interaction stage has an available selection score.',call.=FALSE)
  chosen <- models[[best]]
  columns <- seq_len(nrow(coef(chosen))-added)
  fit <- L0glm_interaction_restore(chosen,L0glm_design_view(X,products),columns,pool,groups)
  fit$screening <- base$screening
  fit$interactions <- list(order=interactions,selected_order=orders[best],terms=products,
    input.names=colnames(X),columns=colnames(view),p.original=p,control=control,
    specification=utils::modifyList(args[setdiff(names(args),c('X','y'))],
      list(interactions=interactions,control.interactions=control)),
    history=data.frame(order=orders,lambda=lambdas,
      df=vapply(metrics,function(s)s$df,1),
      score=score,selected=seq_along(models)==best))
  fit$metrics <- l0glmMetrics(fit,mu_floor=1e-8)
  fit$objective <- fit$metrics$objective
  fit
}
