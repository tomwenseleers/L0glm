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
                                = ifelse(is.null(options), 1, options$error_sd))
  } else if(method == 6){
    sim <- simulate_ordinis(n = n, p = p, m = m, sd_noise = ifelse(is.null(options), 3, options$sd_noise))
  }

  return(sim)
}


fit.GLM <- function(X, y, method = 1, options){
  if(method == 1) { # glmnet
    
  } else if(method == 2) { # adaotive lasso glmnet
    
  } else if(method == 3) { # L0 learn
    
  } else if(method == 4){ # abess
    
  } else if(method == 5){ # abess after LASSO
    
  } else if(method == 6){ # ordinis mcp
    
  } else if(method == 7) # ordinis SCAD
  
}

fit.glmnet <- function(){
  
}


