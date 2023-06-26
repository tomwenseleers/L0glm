#' @importFrom qlcMatrix colMax
#' 
#' @describeIn L0glm
#' Selects optimal value of regularization parameter lambda
#' @export
L0glm.tune <- function(X, y, family, 
                       weights = matrix(1, nrow = nrow(y), ncol = ncol(y)), 
                       offset = Matrix(0, nrow = nrow(X), ncol=ncol(y), sparse=T), 
                       dispersion = NULL, 
                       IC = "aic", 
                       warm.start.ratio = 0.9, 
                       coefs.true = NULL,
                       printbest = TRUE, # Want the best empirical lambda printed out based on the provided IC?
                       plot = TRUE, # Want plots produced?
                       ...){
  
  # Checking if the IC is one of the valid criteria
  IC = tolower(IC)
  if(!(IC %in% c("aic", "bic", "ebic", "gic", "mbic"))){
    stop("Invalid information criteria.")
  }
  
  if(!is.null(dispersion)){
    lambda.preset <- preset.lambda(IC = IC, dispersion = dispersion, family = family)
  } else {
    lambda.preset = preset.lambda(IC = IC, y = y[,which.max(colMeans(abs(y)))],
                                  X = X, weights=rep(1, nrow(y)),
                                  family = family)
  }
  n <- nobs <- nrow(X)
  p <- ncol(X)
  
  
  # for the lambda sequence below we use lambda_preset_aic 
  # (round(log10(lambda_preset_aic))-3) as the lower bound and the theoretical
  # lambda that would penalize all variables to zero with an L0 penalty as maximum
  # see https://stats.stackexchange.com/questions/416144/minimum-and-maximum-regularization-in-l0-pseudonorm-penalized-regression
  maxlambda = max((1/colSums(X^2)) * (crossprod(X, Y)^2)) # TODO make this Xw and Yw
  lambdas = 10^seq((round(log10(lambda.preset))-3), round(log10(maxlambda))+1, by=1)
  
  # Get warm start
  warmup <- L0glm.fit(X = X, y = y,
                          lambda = min(lambdas), 
                          family = family,
                          offset = offset, 
                          weights = weights,
                          ...) # using provided / default params

  warm_start <- warmup$coefficients
  Z <- warmup$response
  W <- warmup$weights
  w.col <- which.max(colMeans(abs(y))) 
  Xw <- Xweights(X, W[,w.col]) # TODO: talk with Tom about this (what column should we use?)
  Zw <- sapply(seq_len(ncol(Z)), function(col) {
    Zw <- Yweights(as.matrix(Z[, col, drop = FALSE]), as.matrix(W[, col, drop=F]))
    return(Zw)
  })

  # Recalculate maxlambda and lambdas with obtained weights 
  maxlambda = max((1/colSums(Xw^2)) * (crossprod(Xw, Zw)^2)) # TODO make this Xw and Yw
  lambdas = 10^seq((round(log10(lambda.preset))-3), round(log10(maxlambda))+1, by=1)
  no.warm.start = matrix(1, nrow=nrow(warm_start), ncol=ncol(warm_start))
  coefficients <- lapply(lambdas, function (lambda) {
    new_coefs <- L0glm.fit(X = X, y = y, 
                           family = family,
                           start = warm.start.ratio*warm_start+(1-warm.start.ratio)*no.warm.start, # warm.start.ratio= 0 => no warm starts,warm.start.ratio= 1 => full warm starts
                           lambda = lambda, # needs to be larger with larger sd_noise for family=gaussian, e.g. 0.5 with sd_noise=1 & 5 with sd_noise=5 
                           weights = weights,
                           offset = offset,
                           ...)$coefficients
    warm_start <<- new_coefs
  }) # 82s, TODO check how warm starts would affect solution quality & stability
  metrics <- data.frame(do.call(rbind, lapply(coefficients, function (coefs) get.metrics(beta = coefs, X = X, y = y, 
                                                                                         weights = weights, 
                                                                                         offset = offset, 
                                                                                         family = family, 
                                                                                         beta_true = coefs.true))))
  
  # metrics <- structure(metrics, class = "tune.L0glm") # Deformats the data frame. TODO: should I make a constructor for the tune.L0glm class?
  if(printbest) {
    bestlambda = lambdas[which.min(metrics[,IC])]
    message <- paste("Best empirical lambda based on", IC, "over lambda sequence is", bestlambda)
    print(message)
  }
  if(plot){
    plot.tune(metrics = metrics, IC = IC, coefs.true = coefs.true, lambdas = lambdas)
  }
  return(metrics)
}

abess.fit <- function(X, y, ...){
  fit <- abess(x = X, 
               y = y, 
               ...)
  
  
  plot(fit)
  return(fit)
}

glmnet.fit <- function(X, y, ...){
  
}

tune.plots <- function(metrics, lambdas, IC = "aic", coefs.true = NULL){
  IC = tolower(IC)
  if(!(IC %in% c("aic", "bic", "ebic", "gic", "mbic"))){
    stop("Invalid information criteria.")
  }
  
  bestlambda = lambdas[which.min(metrics[,IC])]
  bestlambda = lambdas[which.min(abs(lambdas-bestlambda))] 
  beta_relaxedgroupL0 = as.matrix(coefficients[[which(lambdas==bestlambda)]]) 
  IC = toupper(IC) # better for plot titles
  
  # PLOT TRUE VS ESTIMATED COEFFICIENTS 
  if(!is.null(coefs.true)){
    R2 = round(cor(as.vector(coefs.true), as.vector(beta_relaxedgroupL0)), digits = 4)^2 # Calculate R2
    title = paste("True vs estimated coefficients\nnonnegative relaxed L0 group GLM\nR2=", R2)
    plot(as.vector(beta_relaxedgroupL0), as.vector(coefs.true), pch=16, 
         xlab="estimated coefficients nonnegative relaxed L0 group GLM",
         ylab="true coefficients", col="steelblue",
         main=title)
  }
  
  # PLOT ESTIMATED SUPPORT SIZE IN FUNCTION OF REGULARISATION
  title = paste("Estimated support size in function of regularisation\n(red = estimate at best", IC, ", blue = true support size)")
  plot(x = log10(lambdas), y = metrics$support_size_estimated,
       ylab="estimated support size",
       main = title,
       pch=16, type="o")
  abline(v=log10(bestlambda), col="red")
  abline(h=metrics$support_size_true, col="blue")
  
  # PLOT DEVIANCE
  title = paste("Theoretical deviance relative to true underlying signal without noise\n(red=estimate at best", IC,", blue=best deviance)")
  plot(x=log10(lambdas), y=metrics$deviance, pch=16, type="o", 
       ylim=c(0,max(unlist(metrics$deviance))), 
       ylab="deviance",
       main = title) # theoretical deviance relative to true underlying signal without noise
  lambda_bestdeviance = lambdas[which.min(metrics$deviance)]
  abline(v=log10(bestlambda), col="red")
  abline(v=log10(lambda_bestdeviance), col="blue")
  par(oma=c(0,2,0,0))
  
  # PLOT LOG LIKELIHOOD IN FUNCTION OF REGULARISATION
  title = paste("Log likelihood in function of regularisation\n(red=estimate at best", IC,")")
  plot(x=log10(lambdas), y=metrics$loglik,
       main=title,
       ylab="log likelihood",
       pch=16, type="o") # log likelihood
  lambda_bestloglik = lambdas[which.max(metrics$loglik)]
  abline(v=log10(bestlambda), col="red")
  
  # PLOT MAD
  title = paste("Mean absolute deviation of coefficients vs true coefs\n(red=estimate at best", IC,")")
  plot(x=log10(lambdas), y=metrics$MAD_coefs, 
       main=title,
       ylab="Mean absolute deviation coefficients",
       pch=16, type="o")
  lambdas[which.min(metrics$MAD_coefs)] 
  abline(v=log10(bestlambda), col="red")
  
  # PLOT ACCURACY
  title = paste("Accuracy to recover true support\n(red=estimate at best", IC,")")
  plot(x=log10(lambdas), y=metrics$accuracy,
       main=title,
       ylab="accuracy",
       pch=16, type="o") 
  lambda_best_accuracy = lambdas[which.max(metrics$accuracy)] 
  abline(v=log10(bestlambda), col="red")
  
  # TODO add more plots?
}
