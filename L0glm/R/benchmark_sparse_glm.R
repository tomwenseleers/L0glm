#' @name benchmark_sparse_glm
#' @title Benchmark several sparse GLM methods
#' 
#' @description
#' Fit simulated data created by the [simulate()] function with one
#' of several sparse GLM methods. 
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
#' @importFrom L0Learn L0Learn.fit
#' @importFrom ordinis ordinis
#' @importFrom scales alpha
#' 
#' @export
benchmark_sparse_glm <- function(sim, family = NULL, method = 1, options = NULL, 
                                 plot = TRUE, IC = "gic", internal = FALSE){
  
  X <- sim$X
  y <- sim$y
  p <- ncol(sim$X)
  n <- nrow(sim$X)

  if(family$family == "gaussian"){
    W <- rep(1,p)
  } else if(family$family == "poisson"){
    W <- as.vector(1/(y+0.1))
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
           xlab="estimated coefficients", ylab="true coefficients")
    }
    
    coefs <- result$coefficients
    
    
  } else if(method == 2) { # glmnet
    if(is.null(options)){
      result <- glmnet::glmnet(x = X, 
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
      result <- do.call(glmnet::glmnet(args))
    }
    LLs <- logLik(result, weights = W, X = X, y = y)
    edfs = attr(LLs, "df")
    IC.val = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,IC]
    
    bestlambda_lasso_glmnet = result$lambda[which.min(IC.val)]
    coefs = glmnet::coef.glmnet(result, s = bestlambda_lasso_glmnet)
    
    if(plot){
      plot(result, xvar='lambda',
           main = "Regularisation path negative weighted LASSO with Glmnet")
      abline(h=as.vector(sim$beta_true[sim$beta_true!=0]), col=scales::alpha(1:6,0.2))
      
      plot(x = coefs[-1], y=sim$beta_true, pch=16, col='blue',
           main="Glmnet nonnegative weighted LASSO",
           xlab="estimated coefficients", ylab="true coefficients")
    }
    
  } else if(method == 3) { # adaptive LASSO using glmnet 
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
    glmnet.fit <- benchmark_sparse_glm(sim = sim, method = 2, family = family, 
                                       plot = F, internal = TRUE)
    
    LLs <- logLik(glmnet.fit, weights=W, X=X, y=y)
    edfs = attr(LLs, "df")
    IC.val = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,IC]
    bestlambda_lasso_glmnet = glmnet.fit$lambda[which.min(IC.val)]
    coefs.glmnet.fit = glmnet::coef.glmnet(glmnet.fit, s = bestlambda_lasso_glmnet)
    
    # in second step we do another LASSO fit with adaptive penalty weights
    delta = 1E-5
    if(is.null(options)){
      result <- glmnet::glmnet(x = X, 
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
      do.call(glmnet::glmnet, args)
    }
    LLs <- logLik(result, weights=W, X=X, y=y)
    edfs = attr(LLs, "df")
    IC.val = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,IC]
    bestlambda_alasso_glmnet = result$lambda[which.min(IC.val)]
    coefs = glmnet::coef.glmnet(result, s = bestlambda_lasso_glmnet)
    
    if(plot){
      plot(result, xvar='lambda',
           main = "Regularisation path adaptive LASSO with Glmnet")
      abline(h=as.vector(sim$beta_true[sim$beta_true!=0]), col=alpha(1:6,0.2))
      
      plot(x=coefs[-1], y=sim$beta_true, pch=16, col='blue',
           main="True vs. estimated coefficients \n glmnet nonnegative weighted adaptive LASSO",
           xlab="estimated coefficients", ylab="true coefficients")
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
      result <- L0Learn::L0Learn.fit(x = X,  # note: no support for sparse matrices
                            y = y,
                            loss = "SquaredError",
                            penalty = "L0", algorithm="CD",
                            maxSuppSize = round(p/2),
                            intercept = FALSE,
                            nLambda = 100, # NOTE: to make timings fair we should probably use 100 or 20 lambda values everywhere?
                            autoLambda=TRUE)
    } else {
      args = c(list(x = X, y = y), options)
      result <- do.call(L0Learn::L0Learn.fit, args)
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
           xlab="estimated coefficients", ylab="true coefficients")
    }
    
  } else if(method == 5) { # nonnegative L0glm after glmnet LASSO
    # Start with regular LASSO fit
    glmnet.fit <- benchmark_sparse_glm(sim = sim, method = 2, family = family, 
                                       plot = F, internal = TRUE)
    LLs <- logLik(glmnet.fit, weights = W, X = X, y = y)
    edfs = attr(LLs, "df")
    IC.val = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,IC]
    
    bestlambda = glmnet.fit$lambda[which.min(IC.val)]
    glmnet.coefs = glmnet::coef.glmnet(glmnet.fit, s = bestlambda)
    
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
      fit <- L0Learn::L0Learn.fit(x = X[,selcoefs,drop=F],  # note: no support for sparse matrices
                         y = y,
                         loss = "SquaredError",
                         penalty = "L0", algorithm = "CD",
                         maxSuppSize = (p_selcoefs-1), 
                         intercept = FALSE, 
                         nLambda = 100,
                         autoLambda=TRUE)
      
    } else {
      args <- c(list(x = X[,selcoefs,drop=F], y = y, maxSuppSize = (p_selcoefs-1)), options)
      fit <- do.call(L0Learn::L0Learn.fit, args)
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
           xlab="estimated coefficients", ylab="true coefficients")
    } 
    
  } else if(method == 6) { # abess
    # fast for small problems, but becomes slow for n=p=10000
    # if(IC != "ebic") warning("Ebic is recommended with method == 6 for better performance")
    if(is.null(options)){
      fit <- abess::abess(x = X, 
                   y = y,
                   weight = W, 
                   family = "gaussian",
                   tune.path = "sequence",
                   support.size = c(1:round(p*0.5)),
                   fit.intercept = FALSE,
                   beta.low = 0, # nonnegativity constraints
                   lambda=0,
                   warm.start = TRUE,
                   tune.type = "bic", # trick to support any IC that L0glm supports
                   ic.scale = preset.lambda(IC = IC, # options aic, bic, gic, ebic, mbic 
                                            p = p, 
                                            n = n)/log(n)) # trick to support any IC that L0glm supports
    } else {
      args <- c(list(x = X, y = y, weight = W, 
                     family = "gaussian",
                     tune.type = "bic", # trick to support any IC that L0glm supports
                     ic.scale = preset.lambda(IC = IC, # options aic, bic, gic, ebic, mbic 
                                              p = p, 
                                              n = n)/log(n)), # trick to support any IC that L0glm supports, 
                options)
      fit <- do.call(abess::abess, args)
    }
    
    # coefs = as.matrix(abess::extract(fit)$beta) # coefficient matrix for best subset
    coefs = coef(fit,  
                 support.size = fit[["best.size"]], 
                 sparse = FALSE)[-1,] # coefficient matrix for best subset
    
    if(plot){
      plot(fit, type="tune")
      
      plot(x=coefs, y=sim$beta_true, pch=16, col='blue',
           main="True vs. estimated coefficients abess best\nsubset weighted least squares",
           xlab="estimated coefficients", ylab="true coefficients")
    }
    
  } else if(method == 7) { # nonnegative abess after nonnegative LASSO
    # very good solution, and fast for small problems, but
    # becomes very slow for large problems, e.g. n=p=10000
    glmnet.fit <- benchmark_sparse_glm(sim = sim, method = 2, family = family, 
                                       plot = F, internal = TRUE)
    LLs <- logLik(glmnet.fit, weights = W, X = X, y = y)
    edfs = attr(LLs, "df")
    IC.val = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,IC]
    
    bestlambda = glmnet.fit$lambda[which.min(IC.val)]
    glmnet.coefs = glmnet::coef.glmnet(glmnet.fit, s = bestlambda)
    
    selcoefs = which(glmnet.coefs[-1]!=0)
    p_selcoefs = length(selcoefs)
    
    if(is.null(options)){
      fit <- abess::abess(x = X[,selcoefs,drop=F], 
                   y = y, 
                   weight = W,
                   family = "gaussian", # identity link poisson not supported
                   tune.path = "sequence", 
                   support.size = c(1:(p_selcoefs-1)),
                   fit.intercept = FALSE,
                   beta.low = 0, # nonnegativity constraints
                   lambda = 0, 
                   warm.start = TRUE, # FALSE worked better here, but much slower
                   tune.type = "bic", # trick to support any IC that L0glm supports
                   ic.scale = preset.lambda(IC = IC, # options aic, bic, gic, ebic, mbic 
                                            p = p, 
                                            n = n)/log(n)) # trick to support any IC that L0glm supports
    } else {
      args <- c(list(x = X[,selcoefs,drop=F], y = y, 
                     weight = W,
                     family = "gaussian",
                     tune.type = IC, 
                     support.size = c(1:(p_selcoefs-1))), options)
      fit <- do.call(abess::abess, args)
    }
    
    plot(fit, type="tune")  # optimal support size 
    # coefs.selcoefs = as.matrix(abess::extract(fit)$beta) # coefficient matrix for best subset
    coefs.selcoefs = coef(fit,  
                          support.size = fit[["best.size"]], 
                          sparse = FALSE)[-1,] # coefficient matrix for best subset
    coefs = Matrix(rep(0, p), sparse=T)
    coefs[selcoefs] = coefs.selcoefs
    
    if(plot){
      plot(fit, type="tune")  # optimal support size 
      
      plot(x=coefs, y=sim$beta_true, pch=16, col='blue',
           main="True vs estimated coefficients abess best subset\nweighted least squares after nonnegative LASSO",
           xlab="estimated coefficients", ylab="true coefficients")
    }
  } else if(method == 8) { # nonnegative weighted MCP penalized regression in ordinis 
    if(is.null(options)){
      fit <- ordinis::ordinis(X, y, 
                     weights = W,
                     nlambda = 100L, 
                     penalty = "mcp", # using MCP penalty
                     intercept = FALSE,
                     lower.limits = rep(0, p), # force all coefficients to be positive
                     standardize = FALSE,
                     alpha = 1) 
    } else {
      args <- c(list(X = X, y = y, weights = W), options)
      fit <- do.call(ordinis::ordinis, args)
    }
    
    LLs = logLik2(fit, weights= as.vector(W), X, as.vector(y))
    edfs = attr(LLs, "df")
    LLs = as.vector(LLs)
    IC.val = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,IC]
    coefs = fit$beta[,which.min(IC.val[-1])+1]
    
    if(plot){
      plot(fit)
      abline(h=as.vector(sim$beta_true[sim$beta_true!=0]), col=scales::alpha(1:6,0.2))
      
      plot(x=coefs[-1], y=sim$beta_true, pch=16, col='blue',
           main="True vs. estimated coefficients\nordinis nonnegative MCP",
           xlab="estimated coefficients", ylab="true coefficients")
    }
  } else if(method == 9) { # ordinis scad
    if(is.null(options)){
      fit <- ordinis::ordinis(X, y, 
                     weights = W,
                     nlambda = 100L, 
                     penalty = "scad", # using scad penalty
                     intercept = FALSE,
                     lower.limits = rep(0, p), # force all coefficients to be positive
                     standardize = FALSE,
                     alpha = 1) 
    } else {
      args <- c(list(X = X, y = y, weights = W), options)
      fit <- do.call(ordinis::ordinis, args)
    }
    
    LLs = logLik2(fit, weights=W, X, y)
    edfs = attr(LLs, "df")
    LLs = as.vector(LLs)
    IC.val = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,IC]
    coefs = fit$beta[,which.min(IC.val[-1])+1]
    
    if(plot){
      plot(fit)
      abline(h=as.vector(sim$beta_true[sim$beta_true!=0]), col=alpha(1:6,0.2))
      
      plot(x=coefs[-1], y=sim$beta_true, pch=16, col='blue',
           main="True vs. estimated coefficients\nordinis nonnegative SCAD",
           xlab="estimated coefficients", ylab="true coefficients")
    }
  }
  if(internal) return(result)
  return(coefs)
}




