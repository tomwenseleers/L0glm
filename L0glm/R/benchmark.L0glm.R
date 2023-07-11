#' Simulate Data
#'
#' This function simulates different types of data based on the method specified. 
#' Various aspects of the data generation can be controlled using the \code{options} parameter.
#'
#' @param method 
#' A numeric value indicating the type of data to simulate. The possible values are:
#' 1: Single channel (single variable) data with Gaussian noise.
#' 2: Single channel data with Poisson noise.
#' 3: Multichannel data with Gaussian noise.
#' 4: Multichannel data with Poisson noise.
#' 5: Orthogonal data.
#' 6: High dimensional data generation (p > 40000 suggested).
#' @param n 
#' A numeric value indicating the number of observations to generate.
#' @param p 
#' A numeric value indicating the number of predictors (variables) to generate.
#' @param options 
#' A list of optional parameters to control the data generation process. This can include:
#'   For methods 1 to 4: \code{sd_noise} (standard deviation of the noise), 
#'   \code{mean_beta} (mean of the beta coefficients), \code{sd_logbeta} (standard deviation of the log-beta).
#'   For method 5: \code{n_nonzero} (number of nonzero coefficients), \code{beta_sd} (standard deviation of the beta coefficients), 
#'   \code{error_sd} (standard deviation of the error term).
#'   For method 6: \code{sd_noise} (standard deviation of the noise).
#'   If options are \code{NULL}, default values will be used.
#'
#' @return A variable containing the simulated data.
#'
#' @examples
#' \dontrun{
#' sim_data <- simulate(method=1, n=500, p=500)
#' }
#' 
#' @export
simulate <- function(method=1, n = 500, p = 500, options = NULL){
  if(method == 1){
    family = "gausian"
    multichannel = F
    sd_noise <- ifelse(is.null(options), 5, options$sd_noise)
    mean_beta <- ifelse(is.null(options), 2000, options$mean_beta)
  } else if(method == 2){
    family <- "poisson"
    multichannel =  F
    sd_noise <- ifelse(is.null(options), 10, options$sd_noise)
    mean_beta <- ifelse(is.null(options), 2000, options$mean_beta)
  } else if(method == 3){
    family = "gaussian"
    multichannel = T
    sd_noise = ifelse(is.null(options), 5, options$sd_noise)
    mean_beta <- ifelse(is.null(options), 10000, options$mean_beta)
  } else if(method == 4){
    family = "poisson"
    multichannel = T
    sd_noise = ifelse(is.null(options), 10, options$sd_noise)
    mean_beta <- ifelse(is.null(options), 10000, options$mean_beta)
  } else if(method == 6){
    if(p < 40000) warning("This method is suggested for data with higher dimension.")
    m <- ifelse(is.null(options), 50, m)
  }
  
  if(method < 5){
    sim <- simulate_spike_train(n = n, p = p, family = family, 
                                k = round(0.1*n), 
                                multichannel = multichannel,
                                mean_beta = mean_beta,
                                sd_logbeta = ifelse(is.null(options), 1, options$sd_logbeta),
                                sd_noise = sd_noise)
  } else if(method == 5){
    sim <- simulate_orthogonal(n = n, p = p, 
                               n_nonzero = ifelse(is.null(options), 24, options$n_nonzero),
                               beta_sd = ifelse(is.null(options), 0.5, options$beta_sd), 
                               error_sd = ifelse(is.null(options), 1, options$error_sd))
  } else if(method == 6){
    sim <- simulate_ordinis(n = n, p = p, m = m, sd_noise = ifelse(is.null(options), 3, options$sd_noise))
  }
  
  return(sim)
}

#' Fit a simulated object created with \code{\link{simulate}} function.
#' This function simulates different types of data based on the method specified. 
#' Various aspects of the data generation can be controlled using the \code{options} parameter.
#' @param sim
#' @param method 
#' A numeric value indicating the type of data to simulate. The possible values are:
#' 1: Single channel (single variable) data with Gaussian noise.
#' 2: Single channel data with Poisson noise.
#' 3: Multichannel data with Gaussian noise.
#' 4: Multichannel data with Poisson noise.
#' 5: Orthogonal data.
#' 6: High dimensional data generation (p > 40000 suggested).
#' @param options 
#' A list of optional parameters to control the data generation process. This can include:
#'   For methods 1 to 4: \code{sd_noise} (standard deviation of the noise), 
#'   \code{mean_beta} (mean of the beta coefficients), \code{sd_logbeta} (standard deviation of the log-beta).
#'   For method 5: \code{n_nonzero} (number of nonzero coefficients), \code{beta_sd} (standard deviation of the beta coefficients), 
#'   \code{error_sd} (standard deviation of the error term).
#'   For method 6: \code{sd_noise} (standard deviation of the noise).
#'   If options are \code{NULL}, default values will be used.
#'
#' @return A variable containing the simulated data.
#'
#' @examples
#' \dontrun{
#' sim_data <- simulate(method=1, n=500, p=500)
#' }
#' 
#' @export
fit.GLM <- function(sim, family = NULL, method = 1, options = NULL, plot = TRUE, IC = "gic"){
  X <- sim$X
  y <- sim$y
  W <- 1/(y+1E-10)
  
  if(method == 1) { # L0glm
    if(is.null(family)){
      fam.fun <- get(sim$family)
      family <- fam.fun(link = "identity")
    }
    
    if(is.null(options)){
      result <- L0glm.fit(X, y, 
                          family = family, 
                          lambda = 1, 
                          ridge.center = NULL,
                          nonnegative = TRUE,
                          warmupiters = 1,
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
      plot(result, xvar='lambda')
      abline(h=as.vector(sim$beta_true[sim$beta_true!=0]), col=alpha(1:6,0.2))
      
      plot(x = coefs[-1], y=sim$beta_true, pch=16, col='blue',
           main="glmnet nonnegative weighted LASSO",
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
    glmnet.fit <- fit.GLM(sim = sim, method = 2, plot = F)
    
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
      plot(result, xvar='lambda')
      abline(h=as.vector(sim$beta_true[sim$beta_true!=0]), col=alpha(1:6,0.2))
      
      plot(x=coefs[-1], y=sim$beta_true, pch=16, col='blue',
           main="glmnet nonnegative weighted adaptive LASSO",
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
      args = c(list(x = X, y = y, options)
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
           main="nnL0Learn weighted least squares best subset",
           xlab="Estimated coefficients", ylab="True coefficients")
    }
    
  } else if(method == 5) { # nonnegative L0glm after glmnet LASSO
    # Start with regular LASSO fit
    glmnet.fit <- fit.GLM(sim = sim, method = 2, plot = F)
    selcoefs <- which(glmnet.fit[-1]!=0)
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
    
    
    LLs = logLik(fit, weights=W, X[,selcoefs,drop=F], y)
    edfs = attr(LLs, "df")
    LLs = as.vector(LLs)
    IC.val = get.IC(loglik = LLs, n = n, p = p_selcoefs, edf = edfs)[,IC]
    coefs.sel = coef(fit,
                     lambda = fit$lambda[[1]][which.min(IC)],
                     gamma = 0)
    coefs = Matrix(0, p)
    coefs[selcoefs] = coefs.sel
    if(plot){
      plot(fit)
      plot(x=coefs_nnL0Learn_after_lasso, y=sim$beta_true, pch=16, col='blue',
           main="weighted nnL0Learn after nonnegative LASSO",
           xlab="Estimated coefficients", ylab="True coefficients")
    } 
    
  } else if(method == 6) { # abess
    # note: currently no support for nonnegativity constraints
    # fast for small problems, but becomes slow for n=p=10000
    if(IC != "ebic") warning("Ebic is recommended with method == 6 for better performance")
    if(is.null(options)){
      fit <- abess(x = X, 
                   y = y,
                   weights = W, 
                   family = sim$family,
                   tune.path = "sequence",
                   support.size = c(1:round(p*0.5)),
                   lambda=0,
                   warm.start = TRUE,
                   tune.type = IC)
    } else {
      args <- c(list(x = X, y = y, weights = W, tune.type = IC, family = sim$family), options)
      fit <- do.call(abess, args)
    }
    
    coefs = as.matrix(extract(fit_abess)$beta) # coefficient matrix for best subset
    
    if(plot){
      plot(fit, type="tune")
      
      plot(x=coefs_abess, y=sim$beta_true, pch=16, col='blue',
           main="abess best subset weighted least squares",
           xlab="Estimated coefficients", ylab="True coefficients")
    }
    
  } else if(method == 7) { # abess after LASSO
    # very good solution, and fast for small problems, but
    # becomes very slow for large problems, e.g. n=p=10000
    glmnet.fit <- fit.GLM(sim = sim, method = 2, plot = F)
    selcoefs <- which(glmnet.fit[-1]!=0)
    p_selcoefs = length(selcoefs)
    
    if(is.null(options)){
      fit <- abess(x = X[,selcoefs,drop=F], 
                   y = y, 
                   weights = W,
                   family = "gaussian", # poisson not supported
                   tune.path = "sequence", 
                   support.size = c(1:(p_selcoefs-1)),
                   lambda = 0, 
                   warm.start = TRUE, # FALSE worked better here, but much slower
                   tune.type = IC)
    } else {
      args <- c(list(x = X[,selcoefs,drop=F], y = y, tune.type = IC, support.size = c(1:(p_selcoefs-1))), options)
      fit <- do.call(abess, args)
    }
    
    plot(fit_abessafterlasso_selcoefs, type="tune")  # optimal support size 
    coefs.selcoefs = as.matrix(extract(fit)$beta) # coefficient matrix for best subset
    coefs = Matrix(rep(0, p), sparse=T)
    coefs[selcoefs] = coefs.selcoefs
    
    if(plot){
      plot(fit, type="tune")  # optimal support size 
      
      plot(x=coefs, y=sim$beta_true, pch=16, col='blue',
           main="abess best subset weighted least squares\nafter nonnegative LASSO",
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
    
    LLs = logLik2(fit, weights=W, X, y)
    edfs = attr(LLs, "df")
    LLs = as.vector(LLs)
    IC.val = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,IC]
    coefs = fit$beta[,which.min(IC.val)]
    
    
    
    if(plot){
      plot(fit_MCP)
      abline(h=as.vector(sim$beta_true[sim$beta_true!=0]), col=alpha(1:6,0.2))
      
      plot(x=coefs_MCP_ordinis[-1], y=sim$beta_true, pch=16, col='blue',
           main="ordinis nonnegative MCP",
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
      plot(fit_MCP)
      abline(h=as.vector(sim$beta_true[sim$beta_true!=0]), col=alpha(1:6,0.2))
      
      plot(x=coefs_MCP_ordinis[-1], y=sim$beta_true, pch=16, col='blue',
           main="ordinis nonnegative MCP",
           xlab="Estimated coefficients", ylab="True coefficients")
    }
  }
  return(coefs)
}
