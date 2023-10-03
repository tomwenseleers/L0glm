# Functions for benchmarking L0glm algo against other algorithms

#' Fit a simulated object created with \code{\link{simulate}} function.
#' This function simulates different types of data based on the method specified. 
#' Various aspects of the data generation can be controlled using the \code{options} parameter.
#' @param sim
#' Simulated data generated using \code\link{simulate} function. 
#' @param method 
#' A numeric value indicating the method that needs to be used to fit the GLM. 
#' The possible values are:
#' 1: ...
#' @param options 
#' Here specify the options that each method accepts the options that the package
#' that method calls accepts (and to refer to the respective packages and functions)
#' @param internal
#' For internal calls (like method = 3) where we need a starting fit, the function self calls
#' and requires the whole resulting fit (not only the coefficients, which is what is returned)
#'
#' @return A vector or matrix of coefficients.
#'
#' 
#' @export
fit.GLM <- function(sim, family = NULL, method = 1, options = NULL, plot = TRUE, 
                    IC = "gic", internal = FALSE){
  require(scales)
  X <- sim$X
  y <- sim$y
  if(family$family == "gaussian"){
    W <- as.matrix(rep(1,ncol(X)))
  } else if(family$family == "poisson"){
    W <- 1/(y+1E-10)
  }
  
  if(method == 1) { # L0glm
    if(is.null(family)){
      fam.fun <- get(sim$family)
      family <- fam.fun(link = "identity")
    }
    
    if(is.null(options)){
      result <- L0glm.fit(X, y, 
                          family = family, 
                          lambda = 1, 
                          prior.mean = NULL,
                          nonnegative = TRUE,
                          solver = "eigen", 
                          method = 7L)
    } else {
      args <- c(list(X=X, y=y, family=family), options)
      result <- do.call(L0glm.fit, args)
    }
    if(plot){
      plot(x = result$coefficients, y=sim$beta_true, pch=16, col='blue',
           main="L0 glm identity link Poisson with lambda preset to optimise BIC",
           xlab="Estimated coefficients", ylab="True coefficients")
    }
    
    coefs <- result$coefficients
    
    
  } else if(method == 2) { # glmnet
    if(is.null(options)){
      result <- glmnet(x = X, 
                       y = y, 
                       weights = W,
                       family = "gaussian", # poisson(identity) should be supported but gives an error
                       nlambda = 20L,
                       intercept = F, 
                       standardize = F, 
                       alpha = 1, 
                       lower.limits = 0)
    } else {
      args <- c(list(x = X, y = y, family = "gaussian"), options)
      result <- do.call(glmnet(args))
    }
    LLs <- logLik(result, weights = W, X = X, y = y)
    edfs = attr(LLs, "df")
    IC.val = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,IC]
    
    bestlambda_lasso_glmnet = result$lambda[which.min(IC.val)]
    coefs = coef(result, s = bestlambda_lasso_glmnet)
    
    if(plot){
      plot(result, xvar='lambda',
           main = "Regularisation path negative weighted LASSO with Glmnet")
      abline(h=as.vector(sim$beta_true[sim$beta_true!=0]), col=alpha(1:6,0.2))
      
      plot(x = coefs[-1], y=sim$beta_true, pch=16, col='blue',
           main="Glmnet nonnegative weighted LASSO",
           xlab="Estimated coefficients", ylab="True coefficients")
    }
    
  } else if(method == 3) { # adaotive LASSO using glmnet 
    # since glmnet is so much faster than ordinis & supports sparse
    # covariate matrices, let's use that to fit an adaptive LASSO model, cf
    # https://rpubs.com/kaz_yos/alasso
    # note that for the first step we use LASSO rather than ridge coefficients, i.e. we use
    # glmnet.fit and then in a second step use a glmnet LASSO fit with
    # adaptive penalty weights penalty.factor = 1 / (abs(glmnet.fit))+delta)
    # with delta=1E-5
    # NOTE: If this would be iterated further this would be a way to approximate L0 norm
    # penalized regression (cf. Nuel 2017 talk)
    
    # Start with regular LASSO fit
    glmnet.fit <- fit.GLM(sim = sim, method = 2, family = family, plot = F, internal = TRUE)
    
    LLs <- logLik(glmnet.fit, weights=W, X=X, y=y)
    edfs = attr(LLs, "df")
    IC.val = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,IC]
    bestlambda_lasso_glmnet = glmnet.fit$lambda[which.min(IC.val)]
    coefs.glmnet.fit = coef(glmnet.fit, s = bestlambda_lasso_glmnet)
    
    # in second step we do another LASSO fit with adaptive penalty weights
    delta = 1E-5
    if(is.null(options)){
      result <- glmnet(x = X, 
                       y = y, 
                       weights = W,
                       family = "gaussian", # poisson(identity) should be supported but gave an error
                       nlambda = 20L,
                       penalty.factor = 1 / (abs(coefs.glmnet.fit[-1])+delta), # adaptive penalty weights
                       intercept = F, 
                       standardize = F, 
                       alpha = 1, # lasso (lower value, ie elastic net, made solution worse)
                       lower.limits = 0) # for nonnegativity constraints
    } else {
      args <- c(list(x = X, y = y, weights = W, family = "gaussian", penalty.factor = 1/(abs(coefs.glmnet.fit[-1])+delta)), options)
      do.call(glmnet, args)
    }
    LLs <- logLik(result, weights=W, X=X, y=y)
    edfs = attr(LLs, "df")
    IC.val = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,IC]
    bestlambda_alasso_glmnet = result$lambda[which.min(IC.val)]
    coefs = coef(result, s = bestlambda_lasso_glmnet)
    
    if(plot){
      plot(result, xvar='lambda',
           main = "Regularisation path adaptive LASSO with Glmnet")
      abline(h=as.vector(sim$beta_true[sim$beta_true!=0]), col=alpha(1:6,0.2))
      
      plot(x=coefs[-1], y=sim$beta_true, pch=16, col='blue',
           main="True vs. estimated coefficients \n glmnet nonnegative weighted adaptive LASSO",
           xlab="Estimated coefficients", ylab="True coefficients")
    }
    
  } else if(method == 4) { # nonnegative L0 learn
    if(sim$family == "poisson"){
      X = as.matrix(Diagonal(x=sqrt(W)) %*% X)
      y = as.vector(y*sqrt(W))
    } else {
      X = as.matrix(X)
      y = as.vector(y)
    }
    if(is.null(options)){
      result <- L0Learn.fit(x = X,  # note: no support for sparse matrices
                            y = y,
                            loss = "SquaredError",
                            penalty = "L0", algorithm="CD",
                            maxSuppSize = round(p/2),
                            intercept = FALSE,
                            nLambda = 100, # NOTE: to make timings fair we should probably use 100 or 20 lambda values everywhere?
                            autoLambda=TRUE)
    } else {
      args = c(list(x = X, y = y), options)
      result <- do.call(L0learn.fit, args)
    }
    
    LLs = logLik(result, weights = W, X, y)
    edfs = attr(LLs, "df")
    LLs = as.vector(LLs)
    IC.val = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,IC]
    coefs = coef(result, lambda = result$lambda[[1]][which.min(IC.val)], gamma = 0)
    
    if(plot){
      plot(result)
      
      plot(x=coefs, y=sim$beta_true, pch=16, col='blue',
           main="True vs. estimated coefficients\n nnL0Learn weighted least squares best subset",
           xlab="Estimated coefficients", ylab="True coefficients")
    }
    
  } else if(method == 5) { # nonnegative L0glm after glmnet LASSO
    # Start with regular LASSO fit
    glmnet.fit <- fit.GLM(sim = sim, method = 2, family = family, plot = F, internal = TRUE)
    LLs <- logLik(glmnet.fit, weights = W, X = X, y = y)
    edfs = attr(LLs, "df")
    IC.val = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,IC]
    
    bestlambda = glmnet.fit$lambda[which.min(IC.val)]
    glmnet.coefs = coef(glmnet.fit, s = bestlambda)
    
    selcoefs = which(glmnet.coefs[-1]!=0)
    p_selcoefs = length(selcoefs)
    
    if(sim$family == "poisson"){
      X = as.matrix(Diagonal(x=sqrt(W)) %*% X)
      y = as.vector(y*sqrt(W))
    } else {
      X = as.matrix(X)
      y = as.vector(y)
    }
    
    if(is.null(options)){
      fit <- L0Learn.fit(x = X[,selcoefs,drop=F],  # note: no support for sparse matrices
                         y = y,
                         loss = "SquaredError",
                         penalty = "L0", algorithm = "CD",
                         maxSuppSize = (p_selcoefs-1), 
                         intercept = FALSE, 
                         nLambda = 100,
                         autoLambda=TRUE)
      
    } else {
      args <- c(list(x = X[,selcoefs,drop=F], y = y, maxSuppSize = (p_selcoefs-1)), options)
      fit <- do.call(L0learn.fit, args)
    }
    if(plot) plot(fit)
    
    LLs = logLik(fit, weights=W, X[,selcoefs,drop=F], y)
    edfs = attr(LLs, "df")
    LLs = as.vector(LLs)
    IC.val = get.IC(loglik = LLs, n = n, p = p_selcoefs, edf = edfs)[,IC]
    coefs.sel = coef(fit,
                     lambda = fit$lambda[[1]][which.min(IC.val)],
                     gamma = 0)
    coefs = Matrix(0, p)
    coefs[selcoefs] = coefs.sel
    if(plot){
      plot(fit)
      plot(x=coefs, y=sim$beta_true, pch=16, col='blue',
           main="True vs. estimated coefficients\nweighted nnL0Learn after nonnegative LASSO",
           xlab="Estimated coefficients", ylab="True coefficients")
    } 
    
  } else if(method == 6) { # abess
    # fast for small problems, but becomes slow for n=p=10000
    if(IC != "ebic") warning("Ebic is recommended with method == 6 for better performance")
    if(is.null(options)){
      fit <- abess(x = X, 
                   y = y,
                   weights = W, 
                   family = sim$family,
                   tune.path = "sequence",
                   support.size = c(1:round(p*0.5)),
                   fit.intercept = FALSE,
                   beta.low = 0, # nonnegativity constraints
                   lambda=0,
                   warm.start = TRUE,
                   tune.type = IC)
    } else {
      args <- c(list(x = X, y = y, weights = W, tune.type = IC, family = sim$family), options)
      fit <- do.call(abess, args)
    }
    
    coefs = as.matrix(extract(fit)$beta) # coefficient matrix for best subset
    
    if(plot){
      plot(fit, type="tune")
      
      plot(x=coefs, y=sim$beta_true, pch=16, col='blue',
           main="True vs. estimated coefficients abess best\nsubset weighted least squares",
           xlab="Estimated coefficients", ylab="True coefficients")
    }
    
  } else if(method == 7) { # nonnegative abess after nonnegative LASSO
    # very good solution, and fast for small problems, but
    # becomes very slow for large problems, e.g. n=p=10000
    glmnet.fit <- fit.GLM(sim = sim, method = 2, family = family, plot = F, internal = TRUE)
    LLs <- logLik(glmnet.fit, weights = W, X = X, y = y)
    edfs = attr(LLs, "df")
    IC.val = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,IC]
    
    bestlambda = glmnet.fit$lambda[which.min(IC.val)]
    glmnet.coefs = coef(glmnet.fit, s = bestlambda)
    
    selcoefs = which(glmnet.coefs[-1]!=0)
    p_selcoefs = length(selcoefs)
    
    if(is.null(options)){
      fit <- abess(x = X[,selcoefs,drop=F], 
                   y = y, 
                   weights = W,
                   family = "gaussian", # poisson not supported
                   tune.path = "sequence", 
                   support.size = c(1:(p_selcoefs-1)),
                   fit.intercept = FALSE,
                   beta.low = 0, # nonnegativity constraints
                   lambda = 0, 
                   warm.start = TRUE, # FALSE worked better here, but much slower
                   tune.type = IC)
    } else {
      args <- c(list(x = X[,selcoefs,drop=F], y = y, tune.type = IC, support.size = c(1:(p_selcoefs-1))), options)
      fit <- do.call(abess, args)
    }
    
    plot(fit, type="tune")  # optimal support size 
    coefs.selcoefs = as.matrix(extract(fit)$beta) # coefficient matrix for best subset
    coefs = Matrix(rep(0, p), sparse=T)
    coefs[selcoefs] = coefs.selcoefs
    
    if(plot){
      plot(fit, type="tune")  # optimal support size 
      
      plot(x=coefs, y=sim$beta_true, pch=16, col='blue',
           main="True vs estimated coefficients abess best subset\nweighted least squares after nonnegative LASSO",
           xlab="Estimated coefficients", ylab="True coefficients")
    }
  } else if(method == 8) { # nonnegative weighted MCP penalized regression in ordinis 
    if(is.null(options)){
      fit <- ordinis(X, y, 
                     weights = W,
                     nlambda = 100L, 
                     penalty = "mcp", # using MCP penalty
                     intercept = FALSE,
                     lower.limits = rep(0, p), # force all coefficients to be positive
                     standardize = FALSE,
                     alpha = 1) 
    } else {
      args <- c(list(X = X, y = y, weights = W), options)
      fit <- do.call(ordinis, args)
    }
    
    LLs = logLik2(fit, weights= as.vector(W), X, as.vector(y))
    edfs = attr(LLs, "df")
    LLs = as.vector(LLs)
    IC.val = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,IC]
    coefs = fit$beta[,which.min(IC.val)]
    
    if(plot){
      plot(fit)
      abline(h=as.vector(sim$beta_true[sim$beta_true!=0]), col=alpha(1:6,0.2))
      
      plot(x=coefs[-1], y=sim$beta_true, pch=16, col='blue',
           main="True vs. estimated coefficients\nordinis nonnegative MCP",
           xlab="Estimated coefficients", ylab="True coefficients")
    }
  } else if(method == 9) { # ordinis scad
    if(is.null(options)){
      fit <- ordinis(X, y, 
                     weights = W,
                     nlambda = 100L, 
                     penalty = "scad", # using scad penalty
                     intercept = FALSE,
                     lower.limits = rep(0, p), # force all coefficients to be positive
                     standardize = FALSE,
                     alpha = 1) 
    } else {
      args <- c(list(X = X, y = y, weights = W), options)
      fit <- do.call(ordinis, args)
    }
    
    LLs = logLik2(fit, weights=W, X, y)
    edfs = attr(LLs, "df")
    LLs = as.vector(LLs)
    IC.val = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,IC]
    coefs = fit$beta[,which.min(IC.val)]
    
    
    
    if(plot){
      plot(fit)
      abline(h=as.vector(sim$beta_true[sim$beta_true!=0]), col=alpha(1:6,0.2))
      
      plot(x=coefs[-1], y=sim$beta_true, pch=16, col='blue',
           main="True vs. estimated coefficients\nordinis nonnegative SCAD",
           xlab="Estimated coefficients", ylab="True coefficients")
    }
  }
  if(internal) return(result)
  return(coefs)
}





#' Fit a multichannel simulated object created with \code{\link{simulate}} function.
#' This function simulates different types of data based on the method specified. 
#' Various aspects of the data generation can be controlled using the \code{options} parameter.
#' @param sim
#' Simulated data generated using \code\link{simulate} function. 
#' @param method 
#' A numeric value indicating the method that needs to be used to fit the GLM. 
#' The possible values are:
#' 1: ...
#' @param options 
#' Here specify the options that each method accepts the options that the package
#' that method calls accepts (and to refer to the respective packages and functions)
#' @param internal
#' For internal calls (like method = 3) where we need a starting fit, the function self calls
#' and requires the whole resulting fit (not only the coefficients, which is what is returned)
#'
#' @return A vector or matrix of coefficients.
#'
#' 
#' @export
fit.GLM.mc <- function(sim, family = NULL, method = 1, options = NULL, plot = TRUE, IC = "gic", internal = FALSE){
  require(scales)
  X <- sim$X
  y <- sim$y
  p <- ncol(X)
  if(method == 1){ # best subset multitask learning using family="mgaussian" using abess
    if(is.null(options)){
      fit <- abess(x = X, 
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
      fit <- abess(x = X, 
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
      coefs = as.matrix(extract(fit)$beta) # coefficient matrix for best subset
      coefs[abs(coefs)<1E-5] = 0 # apply same thresholding as in L0glm
      library(qlcMatrix)
      image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=coefs^0.1, col = topo.colors(255),
            useRaster=TRUE,
            xlab="Time", ylab="Channel", main="abess multitask weighted mgaussian\n(red = true peaks, cyan = estimated peaks)")
      abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")
      abline(v=(1:nrow(sim$X))[as.vector(rowMin(coefs)!=0)], col="cyan")
      
      R2 = round(cor(as.vector(rowMeans(sim$beta_true)), as.vector(rowMeans(coefs))), digits = 4)^2 # Calculate R2
      title = paste("True vs estimated coefficients\nabess multitask weighted mgaussian\nR2=", R2)
      plot(as.vector(rowMeans(coefs)), as.vector(rowMeans(sim$beta_true)), pch=16, 
           xlab="estimated coefficients nonnegative L0 group GLM",
           ylab="true coefficients", col="steelblue",
           main=title)
    }
    
  }else if(method ==2){ # glmnet nonnegative mgaussian multivariate group LASSO fit. Only for Gaussian
    if(is.null(options)){
      cvfit <- cv.glmnet(x = X, 
                         y = y, 
                         alpha = 1,  # lasso
                         family = "mgaussian", 
                         nlambda = 100,
                         nfolds = 5, 
                         standardize = FALSE,
                         standardize.response = FALSE,
                         intercept = FALSE, 
                         relax = FALSE,
                         lower.limits = rep(0, ncol(X)+1)) # for nonnegativity constraints
      
      # Fit the nonnegative group LASSO glmnet model over a sequence of lambda values
      fit <- glmnet(x = X, 
                    y = y,
                    alpha = 1, # lasso
                    family = "mgaussian", 
                    nlambda = 100,
                    standardize = FALSE,
                    standardize.response = FALSE,
                    intercept = FALSE, 
                    relax = FALSE, 
                    lower.limits = rep(0, ncol(X)+1))
    } else {
      cvfit <- cv.glmnet(x = X, 
                         y = y, 
                         alpha = 1,  # lasso
                         family = "mgaussian", 
                         nlambda = options$nlambda,
                         nfolds = 5, 
                         standardize = options$standardize,
                         standardize.response = options$standardize.response,
                         intercept = options$intercept, 
                         relax = options$relax,
                         lower.limits = rep(0, ncol(X)+1)) # for nonnegativity constraints
      
      # Fit the nonnegative group LASSO glmnet model over a sequence of lambda values
      fit <- glmnet(x = X, 
                    y = y,
                    alpha = 1, # lasso
                    family = "mgaussian", 
                    nlambda = options$nlambda,
                    standardize = options$standardize,
                    standardize.response = options$standardize.response,
                    intercept = options$intercept, 
                    relax = options$relax, 
                    lower.limits = rep(0, ncol(X)+1))
    }
    best_lambda <- cvfit$lambda.1se
    # get the coefficients for each task for the best lambda value
    coefs <- coef(fit, s = best_lambda)
    beta <- do.call(cbind, lapply(seq_len(m), function (channel) as.matrix(coefs[[channel]][-1,,drop=F])))
    
    colnames(beta) = 1:ncol(beta)
    beta[abs(beta)<1E-5] = 0 # slight amount of thresholding of small coefficients using same thresh as used in L0glm
    
    coefs <- beta
    if(plot){
      # plot(cvfit) # TODO remove?
      # plot(fit, xvar="lambda")
      
      # image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=coefs^0.1, col = topo.colors(255),
      #      useRaster=TRUE,
      #      xlab="Time", ylab="Channel", main="nonnegative group Lasso glmnet\n(red = true peaks, cyan = estimated peaks)")
      # abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")
      # abline(v=(1:nrow(sim$X))[as.vector(rowMin(coefs)!=0)], col="cyan")
      
      
      R2 = round(cor(as.vector(rowMeans(sim$beta_true)), as.vector(rowMeans(coefs))), digits = 4)^2 # Calculate R2
      title = paste("True vs estimated coefficients\nnonnegative group Lasso glmnetn\nR2=", R2)
      plot(as.vector(rowMeans(coefs)), as.vector(rowMeans(sim$beta_true)), pch=16, 
           xlab="estimated coefficients nonnegative L0 group GLM",
           ylab="true coefficients", col="steelblue",
           main=title)
      
    }
    
    
  }
  
  
  return(coefs)
}



# UTILITY FUNCTIONS ####

# logLik method for gaussian glmnet
# cf logLik.ordinis in https://github.com/jaredhuling/ordinis/blob/master/R/ordinis_methods.R
logLik.glmnet = function(object, weights, X, y) {
  n  = as.numeric(object$nobs)
  df = Matrix::colSums(object$beta!=0)
  # df = nr of nonzero coefficients (incl intercept) (NOTE: more correct for the LASSO would be to calculate this as the trace of the hat matrix)

  # Compute y_pred for all columns of object$beta at once
  y_pred = X %*% object$beta

  # Compute squared residuals, multiply by weights, and sum over rows to get resid.ss
  resid.ss = Matrix::colSums(weights * (matrix(rep(y, ncol(y_pred)), ncol = ncol(y_pred) ) - y_pred) ^ 2)

  rdf = n - df
  logL = -0.5 * n * (log(2 * pi) - log(rdf) + log(resid.ss)) - 0.5 * rdf
  attr(logL, "df") = df  # nr of nonzero coefficients
  attr(logL, "nobs") = n
  attr(logL, "p") = object$dim[[2]]
  return(logL)
}



# log-likelihood method for gaussian L0Learn fit,
# cf logLik.ordinis in https://github.com/jaredhuling/ordinis/blob/master/R/ordinis_methods.R
logLik.L0Learn = function(object, weights, X, y) {
  n  = object$n
  p = object$p
  df = object$suppSize[[1]]
  # df = nr of nonzero coefficients (incl intercept)
  
  # Compute y_pred for all columns of object$beta at once
  y_pred = X %*% object$beta[[1]]
  
  # Compute squared residuals, multiply by weights, and sum over rows to get resid.ss
  resid.ss = colSums(weights * (matrix(rep(y, ncol(y_pred)), ncol = ncol(y_pred) ) -
                                  y_pred) ^ 2)
  
  rdf = n - df
  logL = -0.5 * n * (log(2 * pi) - log(rdf) + log(resid.ss)) - 0.5 * rdf
  attr(logL, "df") = df  # nr of nonzero coefficients
  attr(logL, "nobs") = n
  attr(logL, "p") = object$dim[[2]]
  return(logL)
}



# log-likelihood method for gaussian ordinis fit, 
# cf logLik.ordinis in https://github.com/jaredhuling/ordinis/blob/master/R/ordinis_methods.R
# but correctly taking into account weights
# NOTE: ordinis is currently not taking into account the weights in
# the calculation of the log likelihood as loss=residual SS do not take into 
# account weights here https://github.com/jaredhuling/ordinis/blob/master/R/ordinis_call.R
# hence the redefinition of the log likelihood here
logLik2 = function(object, weights, X, y) {  
  n  = as.numeric(object$nobs)
  p = object$nvars
  df = Matrix::colSums(object$beta!=0)  
  # df = nr of nonzero coefficients (incl intercept) (NOTE: more correct for the LASSO would be to calculate this as the trace of the hat matrix)
  
  # Compute y_pred for all columns of object$beta at once
  y_pred =  cbind(1, X) %*% object$beta
  
  # Compute squared residuals, multiply by weights, and sum over rows to get resid.ss
  resid.ss = colSums(weights * (matrix(rep(y, ncol(y_pred)), ncol = ncol(y_pred) ) - 
                                  y_pred) ^ 2)
  
  rdf = n - df 
  logL = -0.5 * n * (log(2 * pi) - log(rdf) + log(resid.ss)) - 0.5 * rdf
  attr(logL, "df") = df  # nr of nonzero coefficients
  attr(logL, "nobs") = n
  attr(logL, "p") = p
  return(logL)
}

