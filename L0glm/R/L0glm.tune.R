#' @importFrom qlcMatrix colMax
#' 
#' @description
#' Tunes the lambda parameter for the \code{\link{L0glm}} and \code{\link{L0glm.fit}} 
#' functions. It uses a list of penalty parameter values (\code{lambda}) to fit 
#' L0-penalized GLMs, and selects the best lambda value according to a specified 
#' information criterion. It can handle multiple output variables (channels), 
#' and can produce a plots of the regularization path,true vs. estimated coefficients,
#' estimated support size, deviance, log likelihood, mean, absolute deviation, 
#' accuracy, precision, recall and the IC of interest, all in function of regularization 
#' (some plots are only produced if the user passes a vector of true coefficients
#' \code{coefs.true}).
#' The function produces these plots by calling another function \code{\link{tune.plots}}.
#' This function can also be called separately with the object \code{L0glm.tune$metrics}.
#' 
#' @param X 
#' A matrix of predictor variables.
#' @param y 
#' A vector or matrix of response variables.
#' @param family 
#' A description of the error distribution and link function to be
#' used in the model. This should be a call to a family function. (See
#' \code{stats::family} for details of family functions.)
#' @param weights 
#' A matrix of observation weights.
#' @param offset 
#' A matrix to be added to the linear predictor.
#' @param dispersion 
#' A dispersion parameter for the family function. If not NULL, it is used to 
#' set the lambda sequence.
#' @param IC A character string specifying the information criterion to use for
#' model selection. It should be one of \code{"aic"}, 
#' \code{"bic"}, \code{"gic"}, \code{"ebic"}, \code{"mbic"}.
#' @param warm.start.ratio A numeric value indicating the ratio of warm start 
#' coefficients to use.
#' @param coefs.true A vector or matrix of true coefficient values, for comparison 
#' with the estimated coefficients.
#' @param printbest A logical value indicating whether to print the best lambda 
#' value according to the selected IC (default is TRUE).
#' @param plot A logical value indicating whether to produce the plots mentioned
#' in the description, all in function of regularization (lambda).
#' @param ... Other arguments to pass to the `L0glm.fit` function except for \code{start}.
#' 
#' @return
#' \code{L0glm.tune} returns a list containing the following elements:
#' \item{\code{metrics}} {a data frame containing metrics evaluating the models fit
#' with different regularization parameters.}
#' \item{\code{best.lambda}} {best lambda according to the \code{IC} 
#' criterion specified.}
#' \item{\code{lambdas}} {a vector containing all the lambdas that were used.}
#' 
#' @export
L0glm.tune <- function(X, y, family, 
                       weights = matrix(1, nrow = nrow(as.matrix(y)), ncol = ncol(as.matrix(y))), 
                       offset = Matrix(0, nrow = nrow(X), ncol = ncol(as.matrix(y)), sparse=T),
                       dispersion = NULL, 
                       IC = "aic", 
                       warm.start.ratio = 0.01, 
                       coefs.true = NULL,
                       printbest = TRUE, # Want the best empirical lambda printed out based on the provided IC?
                       plot = TRUE, # Want plots produced?
                       ...){
  require(scales)
  # Checking if the IC is one of the valid criteria
  IC = tolower(IC)
  if(!(IC %in% c("aic", "bic", "ebic", "gic", "mbic"))){
    stop("Invalid information criteria.")
  }

  y <- as.matrix(y)
  n <- nrow(X)
  p <- ncol(X)
  m <- ncol(y)
  if(!is.null(dispersion)){
    lambda.preset <- preset.lambda(IC = IC, dispersion = dispersion, family = family)
  } else if (m > 1){
    lambda.preset = preset.lambda(IC = IC, y = y[,which.max(colMeans(abs(y)))],
                                  X = X, weights=rep(1, nrow(y)),
                                  family = family)
  } else {
    lambda.preset = preset.lambda(IC = IC, y = y, X = X, family = family)
  }
  
  # Get warm start
  warmup <- L0glm.fit(X = X, y = y,
                      lambda = lambda.preset, 
                      family = family,
                      offset = offset, 
                      weights = weights,
                      ...) # using provided / default params
  warm_start <- as.matrix(warmup$coefficients)
  if(sum(warm_start) == 0) stop("Zero coefficients produced.")
  Z <- as.matrix(warmup$z)
  W <- as.matrix(warmup$weights)
  
  # Recalculate maxlambda
  maxlambda <- maxlambda(X = X, weights = W, z = warmup$z, 
                         coefs = warm_start, family = family)

  lambdas <- 10^seq(round(log10(lambda.preset))-3, log10(maxlambda), by=0.2)
                    
  if(length(lambdas) > 30)  lambdas <- 10^seq(round(log10(lambda.preset))-3, log10(maxlambda), by = 0.5)
  if(length(lambdas) > 30)  lambdas <- 10^seq(round(log10(lambda.preset))-3, log10(maxlambda), by = 1)
  
  if(warm.start.ratio > 0){
    no.warm.start = matrix(1, nrow = nrow(warm_start), ncol = ncol(warm_start))
    start <- warm_start*warm.start.ratio + no.warm.start*(1-warm.start.ratio)
    start[warm_start == 0] <- 0 # slightly faster & stable with this
  } else {
    start = NULL
  }

  coefmat <- Matrix(0, nrow=length(lambdas), ncol = p)
  
  metrics <- data.frame(do.call(rbind, lapply(1:(length(lambdas)), function (i) { 
    lambda <- lambdas[[i]]
    L0glm.out <- L0glm.fit(X = X, y = y, 
                           family = family, 
                           weights = weights,
                           lambda = lambda, 
                           start = start,
                           offset = offset,
                           ...)$coefficients 
    if(m > 1){
      coefmat[i,] <<- rowMeans(L0glm.out) # avg coefficients across channels
    } else {
      coefmat[i,] <<- as.vector(L0glm.out)
    }
    return(unlist(c(lambda = lambda, get.metrics(beta = L0glm.out, 
                                                X = X, y = y, 
                                                family = poisson(identity), 
                                                beta_true = coefs.true))))
  })))
  
  if(plot){
    title = ifelse(m > 1, "Regularisation path group L0 Poisson GLM", "Regularisation path L0 Poisson GLM")
    matplot(x = log10(metrics$lambda), y=as.matrix(coefmat), type="l", lty=1,
            xlab = "Log10(lambda)", ylab="Coefficient estimates", main = title)
    abline(v = log10(maxlambda), col="red") # max lambda where all coefs would be penalized to 0
    if(!is.null(coefs.true)) abline(h = as.vector(coefs.true), col = alpha(1:6,0.2)) # TODO alpha gives error
  }
  
  if(printbest) {
    bestlambda = lambdas[which.min(metrics[,IC])]
    message <- paste("Best empirical lambda based on", IC, "over lambda sequence is", bestlambda)
    print(message)
  }
  if(plot){ 
    tune.plots(metrics = metrics, IC = IC, coefficients = coefmat, coefs.true = coefs.true, lambdas = lambdas)
  }
  return(list(metrics = metrics, best.lambda = bestlambda, lambdas = lambdas))
  
  # TODO: add more returns
  # - the optimal values of lambda & the index in lambdas of that optimal value based on given metrics - so this could also 
  # be a dataframe with columns lambda, lambda_index, metric (e.g. different information criteria plus accuracy & F1 score for 
  # support recovery & metrics for bias in coefficients  if beta_true was given
  # - a sparse matrix with the coefficient estimates for those lambda values (rowMeans of the coefficients for multichannel) (with the right optimal threshold applied for each)
  # - a list with the originally estimated sparse coefficient matrices for the different lambda values (but with the right optimal threshold applied for each)
  
}

tune.plots <- function(metrics, coefficients, lambdas, IC = "aic", coefs.true = NULL){
  IC = tolower(IC)
  if(!(IC %in% c("aic", "bic", "ebic", "gic", "mbic"))){
    stop("Invalid information criteria.")
  }
  bestlambda = lambdas[which.min(metrics[,IC])]
  bestlambda = lambdas[which.min(abs(lambdas-bestlambda))] 
  beta_relaxedgroupL0 = as.matrix(coefficients[which(metrics$lambda==bestlambda),]) 
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
  if(!is.null(coefs.true)){
    title = paste("Estimated support size in function of regularisation\n(red = estimate at best", IC, ", cyan = true support size)")
  } else {
    title = paste("Estimated support size in function of regularisation\n(red = estimate at best", IC, ")")
  }
  
  
  plot(x = log10(lambdas), y = metrics$support_size_estimated,
       ylab="estimated support size",
       main = title,
       pch=16, type="o")
  abline(v = log10(bestlambda), col="red")
  if(!is.null(coefs.true)) abline(h = metrics$support_size_true, col="cyan")
  
  # PLOT DEVIANCE
  title = paste("Theoretical deviance relative to true underlying signal without noise\n(red=estimate at best", IC,", cyan=best deviance)")
  plot(x=log10(lambdas), y=metrics$deviance, pch=16, type="o", 
       ylim=c(0,max(unlist(metrics$deviance))), 
       ylab="deviance",
       main = title) # theoretical deviance relative to true underlying signal without noise
  lambda_bestdeviance = lambdas[which.min(metrics$deviance)]
  abline(v=log10(bestlambda), col="red")
  abline(v=log10(lambda_bestdeviance), col="cyan")
  par(oma=c(0,2,0,0))
  
  # PLOT LOG LIKELIHOOD IN FUNCTION OF REGULARISATION
  title = paste("Log likelihood in function of regularisation\n(red = estimate at best", IC,", cyan = estimate at best log likelihood)")
  plot(x=log10(lambdas), y=metrics$loglik,
       main=title,
       ylab="log likelihood",
       pch=16, type="o") 
  lambda_bestloglik = lambdas[which.max(metrics$loglik)]
  abline(v=log10(bestlambda), col="red")
  abline(v=log10(lambda_bestloglik), col="cyan")
  
  # PLOT MAD
  # TODO plot MD too?
  title = paste("Mean absolute deviation of coefficients vs true coefs\n(red = estimate at best", IC,")")
  plot(x=log10(lambdas), y=metrics$MAD_coefs, 
       main=title,
       ylab="Mean absolute deviation coefficients",
       pch=16, type="o")
  abline(v=log10(bestlambda), col="red")
  
  # PLOT ACCURACY
  title = paste("Accuracy to recover true support\n(red = estimate at best", IC,", cyan = estimate at best accuracy)")
  plot(x=log10(lambdas), y=metrics$accuracy,
       main=title,
       ylab="Accuracy",
       pch=16, type="o") 
  lambda_best_accuracy = lambdas[which.max(metrics$accuracy)] 
  abline(v=log10(bestlambda), col="red")
  abline(v=log10(lambda_best_accuracy), col="cyan")
  
  # PLOT PRECISION
  title = paste("Precision in function of regularisation\n(red = estimate at best", IC, ", blue = estimate at best precision)")
  plot(x=log10(lambdas), y=metrics$precision,
       main=title,
       ylab="Precision",
       pch=16, type="o") 
  lambda_best_precision = lambdas[which.max(metrics$precision)] 
  abline(v=log10(bestlambda), col="red")
  abline(v=log10(lambda_best_precision), col="cyan")
  
  # PLOT RECALL
  # TODO should it be precision and recall separately or F1 score?
  title = paste("Recall in function of regularisation\n(red = estimate at best", IC, ", blue = estimate at best recall)")
  plot(x=log10(lambdas), y=metrics$recall,
       main=title,
       ylab="Recall",
       pch=16, type="o") 
  lambda_best_recall = lambdas[which.max(metrics$recall)] 
  abline(v=log10(bestlambda), col="red")
  abline(v=log10(lambda_best_recall), col="cyan")
  
  # PLOT SPECIFICITY
  title = paste("Specificity in function of regularisation\n(red = estimate at best", IC, ", blue = estimate at best specificity)")
  plot(x=log10(lambdas), y=metrics$specificity,
       main=title,
       ylab="Specificity",
       pch=16, type="o") 
  lambda_best_specficity = lambdas[which.max(metrics$specificity)] 
  abline(v=log10(bestlambda), col="red")
  abline(v=log10(lambda_best_specficity), col="cyan")
  
  
  # R2
  title = paste("R2 value in function of regularisation\n(red = estimate at best", IC,")")
  plot(x=log10(lambdas), y = metrics$R2_coefs,
       main = title, 
       ylab = "R2", 
       pch = 16, type="o")
  abline(v=log10(bestlambda), col="red")
  
  # TP, FP, TN, FN
  color = c("steelblue", "purple")
  title = paste("Effect of regularisation on True and False Positives\n(red = estimate at best", IC, ")")
  plot(x=log10(lambdas), y = metrics$TP,
       main = title,
       ylab = "Value",
       col = color[1],
       pch = 16,
       type = "o")
  lines(log10(metrics$lambda), metrics$FP, type="o", pch=16, col=color[2])
  # lines(log10(metrics$lambda), metrics$TN, type="o", pch=16, col=color[1])
  # lines(log10(metrics$lambda), metrics$FN, type="o", pch=16, col=color[4])
  
  legend("topright", legend=c("TP", "FP"), col=color, pch=16, lty=1)
  abline(v=log10(bestlambda), col="red")
  
  # PLOT IC OF INTEREST
  # TODO should it also plot other ICs?
  title = paste("Effect of regularisation on", IC, "\n(red = estimate at best", IC, ")")
  plot(x = log10(lambdas), y = metrics[[paste0(tolower(IC))]],
       main = title,
       ylab = IC, 
       pch = 16,
       type = "o")
  abline(v=log10(bestlambda), col="red")
}

optimize.threshold <- function(thresh, IC = "bic", coefficients, Xmat, y, family){
  coefficients[abs(coefficients) < thresh[1]] = 0
  IC.val <- get.metrics(beta = coefficients, 
                        X = Xmat, y = y, 
                        family = family)[[paste0(IC)]]
  return(IC.val)
}


