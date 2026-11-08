# Whole-block moves reuse the family-specific restricted refits and score contexts.
L0glm_group_refine <- function(fit,args,d,control) {
  d$allowed <- fit$group.allowed
  control$block.max <- Inf;control$channel.refit <- FALSE
  B <- as.matrix(coef(fit));label <- fit[["group"]]
  groups <- list(id=match(label,unique(label)))
  blocks <- L0glm_group_blocks(B,groups$id,d$grouped)
  eligible <- d$allowed;eligible[d$mandatory,] <- FALSE
  units <- split(which(eligible),blocks[eligible])
  ids <- as.integer(names(units))
  cost <- fit$group.multiplier[unique(label)[(ids-1L)%%max(groups$id)+1L]]
  score <- function(model) {
    s <- L0glm_local_score(model,d,control)
    charge <- if(control$target=='IC')s$df else
      sum(cost[vapply(units,function(i)any(coef(model)[i]!=0),TRUE)])
    s$objective <- -2*s$loglik+d$lambda*charge
    s$criterion <- s$objective
    s$value <- if(is.null(s$validation))s$objective else s$validation
    s
  }
  decorate <- function(model) {
    for(name in c('group','group.multiplier','group.allowed','group.cost.type'))model[[name]] <- fit[[name]]
    model$group.active <- L0glm_group_active(as.matrix(coef(model)),d$allowed,d$mandatory,d$grouped,groups)
    model$metrics <- model$objective <- NULL
    model
  }
  current <- fit;value <- initial <- score(current)
  evaluations <- 0L;failures <- character();screen.seconds <- refit.seconds <- 0
  history <- list(data.frame(iteration=0L,move='initial',value=value$value,GIC=value$GIC,
    selected=value$selected,nonzero=value$nonzero,evaluations=0L))
  stopping <- 'iteration_limit'
  for(iteration in seq_len(control$maxit)) {
    B <- as.matrix(coef(current));base <- current$group.active
    active <- vapply(units,function(i)any(B[i]!=0),TRUE)
    best <- list(model=current,score=value,move='none')
    cache <- new.env(parent=emptyenv())
    evaluate <- function(mask,move) {
      if(evaluations>=control$max.evals)return(NULL)
      key <- paste(which(mask),collapse=',')
      if(exists(key,cache,inherits=FALSE))return(get(key,cache,inherits=FALSE))
      evaluations <<- evaluations+1L;started <- proc.time()[[3L]]
      model <- tryCatch(L0glm_local_refit(mask,current,d,control),error=identity)
      refit.seconds <<- refit.seconds+proc.time()[[3L]]-started
      if(inherits(model,'error')) {
        failures <<- c(failures,conditionMessage(model));return(NULL)
      }
      model <- decorate(model);s <- score(model)
      trial <- list(model=model,score=s,move=move);assign(key,trial,cache)
      if(L0glm_local_better(s$value,best$score$value,control$tol))best <<- trial
      trial
    }
    if(iteration==1L)evaluate(base,'refit')
    started <- proc.time()[[3L]]
    context <- L0glm_local_context(current,d,control$adjust.max)
    statistics <- sis_score_cpp(L0glm_sis_design(d$D),context$context,as.integer(d$p),TRUE,
      as.integer(control$nthreads))$statistics
    statistic <- vapply(units,function(i)sum(statistics[unique((i-1L)%%d$p+1L)]),1)/cost
    add <- head(which(!active)[order(statistic[!active],decreasing=TRUE)],control$top)
    scaled <- B*fit$preprocessing$scale
    norm <- vapply(units,function(i)sum(scaled[i]^2),1)
    screen.seconds <- screen.seconds+proc.time()[[3L]]-started
    drops <- head(which(active)[order(norm[active]/cost[active])],control$drop)
    for(g in add) {mask <- base;mask[units[[g]]] <- TRUE;evaluate(mask,paste0('add_group:',ids[g]))}
    if(length(add)>1L && control$max.swap>1L) {
      mask <- base;mask[unlist(units[head(add,control$max.swap)])] <- TRUE;evaluate(mask,'add_groups')
    }
    for(g in drops) {
      mask <- base;mask[units[[g]]] <- FALSE;evaluate(mask,paste0('drop_group:',ids[g]))
      for(h in add) {
        exchange <- mask;exchange[units[[h]]] <- TRUE
        evaluate(exchange,paste0('swap_groups:',ids[g],':',ids[h]))
      }
    }
    if(best$move=='none') {stopping <- 'neighbourhood';break}
    current <- best$model;value <- best$score
    d$reference.eta <- L0glm_local_eta(as.matrix(coef(current)),d)
    d$reference.B <- as.matrix(coef(current))
    history[[length(history)+1L]] <- data.frame(iteration,move=best$move,value=value$value,
      GIC=value$GIC,selected=value$selected,nonzero=value$nonzero,evaluations)
    if(control$trace)cat('Group refinement:',best$move,value$value,'\n')
    if(evaluations>=control$max.evals) {stopping <- 'evaluation_limit';break}
  }
  current$metrics <- l0glmMetrics(current,mu_floor=control$mu.floor)
  current$objective <- current$metrics$group_objective
  current$local.search <- list(initial=initial,final=value,history=do.call(rbind,history),
    evaluations=evaluations,failures=unique(failures),stopping=stopping,
    screen.seconds=screen.seconds,refit.seconds=refit.seconds,scoring.floor=control$mu.floor,
    units='predictor_groups',pool=d$pool,lambda=d$lambda,control=control)
  if(!is.null(d$validation)) {
    current$validation <- fit$validation
    current$validation$deviance <- value$validation
    current$validation$local.search <- c(initial=initial$validation,final=value$validation)
  }
  if(identical(d$kind,'cox'))current$baseline <- L0glm_ext_cox_baseline(current)
  current
}
