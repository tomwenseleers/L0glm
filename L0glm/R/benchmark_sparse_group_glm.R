#' @name benchmark_sparse_group_glm
#' @title Benchmark several sparse group penalized GLM methods
#' 
#' @description
#' Fit simulated multichannel data created by the [simulate()] function
#' by one of several sparse group penalized GLM methods.
#' Various aspects of the data generation can be controlled using the \code{options} parameter.
#' @param sim Simulated data generated using [simulate()] function. 
#' @param family XXX
#' @param method A numeric value indicating the method that needs to be used to fit the GLM. 
#' The possible values are:
#' \enumerate{
#'   \item 1: ...
#' }
#' @param options Here specify the options that each method accepts the options that the package
#' that method calls accepts (and to refer to the respective packages and functions)
#' @param plot XXX
#' @param IC XXX
#' @param internal For internal calls (like method = 3) where we need a starting fit, the function self calls
#' and requires the whole resulting fit (not only the coefficients, which is what is returned)
#' @return A vector or matrix of coefficients.
#' 
#' @importFrom glmnet glmnet cv.glmnet coef.glmnet
#' @importFrom abess abess extract
#' 
#' @export
benchmark_sparse_group_glm <- function(sim, 
                                       family = NULL, method = 1, 
                                       options = NULL, plot = TRUE, IC = "gic", 
                                       internal = FALSE){

  X <- sim$X
  y <- sim$y
  p <- ncol(sim$X)
  n <- nrow(sim$X)
  m <- ncol(sim$y)
  
  if(method == 1){ # best subset multitask learning using family="mgaussian" using abess
    if(is.null(options)){
      fit <- abess:abess(x = X, 
                   y = y, 
                   weights = ifelse(family$family == "poisson", 1/(y+0.1), NULL),  # approx 1/variance weights to approximate family="poisson", which is not supported
                   family = "mgaussian",
                   tune.path = "sequence", 
                   support.size = c(1:(p-1)),
                   fit.intercept = FALSE,
                   beta.low = 0, # nonnegativity constraints
                   lambda=0, 
                   warm.start = TRUE,
                   tune.type = IC)
    } else {
      fit <- abess:abess(x = X, 
                   y = y, 
                   weights = options$weights,  # approx 1/variance weights to approximate family="poisson", which is not supported
                   family = options$family,
                   tune.path = options$tune.path, 
                   support.size = c(1:(p-1)),
                   fit.intercept = FALSE,
                   beta.low = 0, # nonnegativity constraints
                   lambda=options$lambda, 
                   warm.start = options$warm.start,
                   tune.type = IC)
    }
    
    if(plot){
      plot(fit, type="tune")  # optimal support size 
      coefs = as.matrix(abess::extract(fit)$beta) # coefficient matrix for best subset
      coefs[abs(coefs)<1E-5] = 0 # apply same thresholding as in L0glm
      
      image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=coefs^0.1, col = topo.colors(255),
            useRaster=TRUE,
            xlab="Time", ylab="Channel", main="abess multitask weighted mgaussian\n(red = true peaks, cyan = estimated peaks)")
      abline(v=(1:nrow(sim$X))[as.vector(qlcMatrix::rowMax(sim$beta_true)!=0)], col="red")
      abline(v=(1:nrow(sim$X))[as.vector(qlcMatrix::rowMin(coefs)!=0)], col="cyan")
      
      R2 = round(stats::cor(as.vector(rowMeans(sim$beta_true)), 
                            as.vector(rowMeans(coefs))), digits = 4)^2 # Calculate R2
      title = paste("True vs estimated coefficients\nabess multitask weighted mgaussian\nR2=", R2)
      plot(as.vector(rowMeans(coefs)), as.vector(rowMeans(sim$beta_true)), pch=16, 
           xlab="estimated coefficients nonnegative L0 group GLM",
           ylab="true coefficients", col="steelblue",
           main=title)
    }
    
  }else if(method ==2){ # glmnet nonnegative mgaussian multivariate group LASSO fit. Only for Gaussian
    if(is.null(options)){
      cvfit <- glmnet::cv.glmnet(x = X, 
                         y = y, 
                         alpha = 1,  # lasso
                         family = "mgaussian", 
                         nlambda = 100,
                         nfolds = 5, 
                         standardize = FALSE,
                         standardize.response = FALSE,
                         intercept = FALSE, 
                         relax = FALSE,
                         lower.limits = 0) # for nonnegativity constraints
      
      # Fit the nonnegative group LASSO glmnet model over a sequence of lambda values
      fit <- glmnet::glmnet(x = X, 
                    y = y,
                    alpha = 1, # lasso
                    family = "mgaussian", 
                    nlambda = 100,
                    standardize = FALSE,
                    standardize.response = FALSE,
                    intercept = FALSE, 
                    relax = FALSE, 
                    lower.limits = rep(0, ncol(X)))
    } else {
      cvfit <- glmnet::cv.glmnet(x = X, 
                         y = y, 
                         alpha = 1,  # lasso
                         family = "mgaussian", 
                         nlambda = options$nlambda,
                         nfolds = 5, 
                         standardize = options$standardize,
                         standardize.response = options$standardize.response,
                         intercept = options$intercept, 
                         relax = options$relax,
                         lower.limits = 0) # for nonnegativity constraints
      
      # Fit the nonnegative group LASSO glmnet model over a sequence of lambda values
      fit <- glmnet::glmnet(x = X, 
                    y = y,
                    alpha = 1, # lasso
                    family = "mgaussian", 
                    nlambda = options$nlambda,
                    standardize = options$standardize,
                    standardize.response = options$standardize.response,
                    intercept = options$intercept, 
                    relax = options$relax, 
                    lower.limits = 0)
    }
    best_lambda <- cvfit$lambda.1se
    # get the coefficients for each task for the best lambda value
    coefs <- glmnet::coef.glmnet(fit, s = best_lambda)
    beta <- do.call(cbind, lapply(seq_len(m), 
                                  function (channel) as.matrix(coefs[[channel]][-1,,drop=F])))
    
    colnames(beta) = 1:ncol(beta)
    beta[abs(beta)<1E-5] = 0 # slight amount of thresholding of small coefficients using same thresh as used in L0glm
    
    coefs <- beta
    if(plot){
      # plot(cvfit) # TODO remove?
      # plot(fit, xvar="lambda")
      
      # image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=coefs^0.1, col = topo.colors(255),
      #      useRaster=TRUE,
      #      xlab="Time", ylab="Channel", main="nonnegative group Lasso glmnet\n(red = true peaks, cyan = estimated peaks)")
      # abline(v=(1:nrow(sim$X))[as.vector(qlcMatrix::rowMax(sim$beta_true)!=0)], col="red")
      # abline(v=(1:nrow(sim$X))[as.vector(qlcMatrix::rowMin(coefs)!=0)], col="cyan")
      
      
      R2 = round(stats::cor(as.vector(rowMeans(sim$beta_true)), 
                            as.vector(rowMeans(coefs))), digits = 4)^2 # Calculate R2
      title = paste("True vs estimated coefficients\nnonnegative group Lasso glmnetn\nR2=", R2)
      plot(as.vector(rowMeans(coefs)), as.vector(rowMeans(sim$beta_true)), pch=16, 
           xlab="estimated coefficients nonnegative L0 group GLM",
           ylab="true coefficients", col="steelblue",
           main=title)
      
    }
    
    
  }
  
  
  return(coefs)
}


