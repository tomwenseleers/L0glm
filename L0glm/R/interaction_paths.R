L0glm_interaction_keys <- function(fit) {
  recipe <- fit$interactions
  c(if(fit$added_intercept)'intercept',paste0('m:',seq_len(recipe$p.original)),
    vapply(recipe$terms,function(i)paste0('i:',paste(i,collapse=':')),''))
}

L0glm_interaction_truth <- function(truth,fit) {
  if(is.null(truth) || is.null(fit$interactions))return(truth)
  beta <- coef(fit);recipe <- fit$interactions
  B <- as.matrix(truth)
  if(ncol(B)!=ncol(beta))stop('coefs.true has the wrong number of response columns.',call.=FALSE)
  original <- recipe$p.original+as.integer(fit$added_intercept)
  if(is.null(rownames(B)) && nrow(B)==original)
    return(rbind(B,matrix(0,nrow(beta)-original,ncol(beta))))
  if(identical(dim(B),dim(beta)) && is.null(rownames(B)))return(B)
  if(is.null(rownames(B)) || anyDuplicated(rownames(B)))
    stop('Supply original main-effect truth or uniquely named interaction truth.',call.=FALSE)
  index <- match(rownames(B),rownames(beta))
  if(anyNA(index))stop('Named truth includes terms absent from this interaction path.',call.=FALSE)
  result <- matrix(0,nrow(beta),ncol(beta),dimnames=dimnames(beta));result[index,] <- B
  result
}

L0glm_interaction_finish_path <- function(result,models,X,p,truth,mu_floor) {
  keys <- lapply(models,L0glm_interaction_keys)
  terms <- list()
  for(model in models)for(term in model$interactions$terms) {
    key <- paste(term,collapse=':')
    if(!key%in%names(terms))terms[[key]] <- term
  }
  names(terms) <- NULL
  original <- models[[1]]$interactions$p.original
  full <- p+length(terms)
  for(i in seq_along(models)) {
    model <- models[[i]];recipe <- model$interactions
    group <- model[["group"]]
    if(!is.null(group)) {
      added <- as.integer(model$added_intercept)
      group <- group[seq_len(original)+added]
    }
    groups <- L0glm_interaction_groups(group,terms,model)
    target <- L0glm_design_view(X,terms)
    new_keys <- c(if(model$added_intercept)'intercept',paste0('m:',seq_len(original)),
      vapply(terms,function(j)paste0('i:',paste(j,collapse=':')),''))
    map <- match(keys[[i]],new_keys)
    if(model$added_intercept)map <- map[-1L]-1L
    model <- L0glm_interaction_restore(model,target,map,full,groups)
    recipe$terms <- terms;recipe$columns <- target$columns
    model$interactions <- recipe
    models[[i]] <- model
  }
  scores <- lapply(models,function(model) {
    s <- l0glmMetrics(model,mu_floor=mu_floor,beta_true=L0glm_interaction_truth(truth,model))
    if(is.null(s$support_size_estimated))s$support_size_estimated <- s$selected
    s
  })
  iterations <- result$metrics$maxit
  result$metrics <- do.call(rbind,lapply(seq_along(models),function(i)
    cbind(lambda=result$lambdas[i],scores[[i]])))
  if(!is.null(iterations))result$metrics$maxit <- iterations
  if(result$tune=='cv') {
    result$metrics$cv_deviance <- result$cv$mean
    result$metrics$cv_se <- result$cv$se
  }
  selection <- if(result$tune=='cv')result$metrics$cv_deviance else result$metrics[[result$tune]]
  result$metrics$tuning_score <- selection
  best <- which.min(selection)
  result$lambda_best <- result$lambdas[best];result$fit_best <- models[[best]]
  if(!is.null(result$lambda_1se)) {
    one <- match(result$lambda_1se,result$lambdas)
    if(!is.null(iterations))one <- which(result$lambdas==result$lambda_1se &
      iterations==result$maxit_1se)[1L]
    result$fit_1se <- models[[one]]
  }
  result$coefficients <- lapply(models,coef);result$p <- full
  result$preprocessing <- result$fit_best$preprocessing
  result$interaction_terms <- terms
  result$interaction_keys <- L0glm_interaction_keys(result$fit_best)
  list(result=result,truth=L0glm_interaction_truth(truth,result$fit_best))
}
