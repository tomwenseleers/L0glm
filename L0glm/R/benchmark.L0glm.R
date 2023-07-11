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
      
      if(plot){
        plot(result, xvar='lambda')
        abline(h=as.vector(sim$beta_true[sim$beta_true!=0]), col=alpha(1:6,0.2))
      }
    }
    
  } else if(method == 4) { # loglik glmnet, cf logLik.ordinis in https://github.com/jaredhuling/ordinis/blob/master/R/ordinis_methods.R
    glmnet.fit <- fit.GLM(sim = sim, method = 2, options = options, plot = F)
    LLs <- logLik(glmnet.fit, weights = W, X = X, y = y)
    edfs = attr(LLs, "df")
    IC = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,IC]
    which.min(IC) # 
    
    bestlambda_lasso_glmnet = glmnet.fit$lambda[which.min(IC)]
    coefs_lasso_glmnet = coef(glmnet.fit, s = bestlambda_lasso_glmnet)
    
    if(plot){
      plot(x=coefs_lasso_glmnet[-1], y=sim$beta_true, pch=16, col='blue',
           main="glmnet nonnegative weighted LASSO",
           xlab="Estimated coefficients", ylab="True coefficients")
    }
  } else if(method == 5) { # adaotive LASSO using glmnet 
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
    glmnet.fit <- fit.GLM(sim = sim, method = 2, options = options, plot = F)
    
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
    
    if(plot){
      plot(result, xvar='lambda')
      abline(h=as.vector(sim$beta_true[sim$beta_true!=0]), col=alpha(1:6,0.2))
      
      coefs.result = coef(result, s = bestlambda_lasso_glmnet)
      plot(x=coefs.result[-1], y=sim$beta_true, pch=16, col='blue',
           main="glmnet nonnegative weighted adaptive LASSO",
           xlab="Estimated coefficients", ylab="True coefficients")
    }
  } else if(method == 6) { # L0 learn
    
  } else if(method == 7) { # abess
    
  } else if(method == 8) { # abess after LASSO
    
  } else if(method == 9) { # ordinis mcp
    
  } else if(method == 10) { # ordinis scad
    
  }
  
}

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
