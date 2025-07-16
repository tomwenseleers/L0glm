#' @name simul
#' @aliases simul
#' @title Simulate data.
#' 
#' @description
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
#' 6: High dimensional data (p > 40000 suggested).
#' @param n 
#' A numeric value indicating the number of observations to generate.
#' @param p 
#' A numeric value indicating the number of predictors (variables) to generate.
#' @param S 
#' Sparsity of true solution (proportion of total coefficients p that are nonzero).
#' @param sparse
#' Boolean indicating if the returned design matrix should be sparse (of Matrix class dsCMatrix) or dense.
#' @param options 
#' A list of optional parameters to control the data generation process. This can include:
#'   For methods 1 to 4: \code{sd_noise} (standard deviation of the noise), 
#'   \code{mean_beta} (mean of the beta coefficients), \code{sd_logbeta} (standard deviation of the log-beta).
#'   For method 5: \code{n_nonzero} (number of nonzero coefficients), \code{beta_sd} (standard deviation of the beta coefficients), 
#'   \code{error_sd} (standard deviation of the error term).
#'   For method 6: \code{sd_noise} (standard deviation of the noise).
#'   If options are \code{NULL}, default values will be used.
#' @param seed
#' Random seed.
#' 
#' @return A variable containing the simulated data.
#'
#' @importFrom statmod rinvgauss
#' @importFrom stats gaussian plogis poisson qlogis rbinom rgamma rpois runif
#' @importFrom methods as is
#' @importFrom utils flush.console
#'
#' @examples
#' \dontrun{
#' sim_data <- simul(method=1, n=500, p=500)
#' }
#' 
#' @export
simul <- function(method = 1, n = 500, p = 500, s = 0.1, sparse = T, options = NULL, seed = 1){
  if(method == 1){
    family = "gaussian"
    multichannel = FALSE
    k <- round(s*n)
    sd_noise <- ifelse(is.null(options), 5, options$sd_noise)
    mean_beta <- ifelse(is.null(options$mean_beta), 2000, options$mean_beta)
  } else if(method == 2){
    family = "poisson"
    k <- round(s*n)
    multichannel =  FALSE
    mean_beta <- ifelse(is.null(options$mean_beta), 2000, options$mean_beta)
  } else if(method == 3){
    family = "gaussian"
    multichannel = TRUE
    k <- round(s*n)
    if(k > 100) k <- 100
    sd_noise = ifelse(is.null(options$sd_noise), 5, options$sd_noise)
    mean_beta <- ifelse(is.null(options$mean_beta), 10000, options$mean_beta)
  } else if(method == 4){
    family = "poisson"
    multichannel = TRUE
    k <- round(s*n)
    if(k > 100) k <- 100
    sd_noise = ifelse(is.null(options$sd_noise), 10, options$sd_noise)
    mean_beta <- ifelse(is.null(options$mean_beta), 10000, options$mean_beta)
  } else if(method == 6){
    if(p < 40000) warning("This method is suggested for data with higher dimension.")
    m <- ifelse(is.null(options), 50, m)
  }
  
  if(method == 1 | method == 3){
    sim <- simulate_spike_train(n = n, p = p, family = family, 
                                seed = seed,
                                k = k, 
                                multichannel = multichannel,
                                mean_beta = mean_beta,
                                sparse = sparse,
                                sd_logbeta = ifelse(is.null(options$sd_logbeta), 1, options$sd_logbeta),
                                sd_noise = sd_noise)
  } else if(method == 2 | method == 4){
    sim <- simulate_spike_train(n = n, p = p, family = family, 
                                seed = seed,
                                k = k, 
                                multichannel = multichannel,
                                mean_beta = mean_beta,
                                sparse = sparse,
                                sd_logbeta = ifelse(is.null(options$sd_logbeta), 1, options$sd_logbeta))
  } else if(method == 5){
    sim <- simulate_standardnormal(n = n, p = p, 
                               n_nonzero = ifelse(is.null(options$n_nonzero), 24, options$n_nonzero),
                               beta_sd = ifelse(is.null(options$beta_sd), 0.5, options$beta_sd), 
                               error_sd = ifelse(is.null(options$error_sd), 1, options$error_sd))
  } else if(method == 6){
    sim <- simulate_ordinis(n = n, p = p, m = m, sd_noise = ifelse(is.null(options$sd_noise), 3, options$sd_noise))
  }
  
  return(sim)
}




# simulate single channel data with L2 norm renormalized standard normal variables & gaussian noise with SD of 1 ####
# cf. Frommlet & Nuel 2016 p. 9

simulate_standardnormal <- function(n=100, p=100, n_nonzero=24, beta_sd=0.5, error_sd=1, 
                                    seed.problem=1, seed.noise=sample.int(1E4, 1, replace=T)) {
  # Generate X matrix
  set.seed(seed.problem)
  X <- matrix(rnorm(n*p), n, p)
  
  # # we L2 norm normalize the columns of X first
  # X.n <- sqrt(Matrix::colSums(X^2)) 
  # X <- as.matrix(X %*% Matrix::Diagonal(x = 1 / X.n))
  
  # Generate beta vector with all zeros initially
  beta <- rep(0, p)
  
  # Randomly select n_nonzero indices to assign non-zero values
  non_zero_indices <- sample.int(n=p, size=n_nonzero)
  
  # Generate beta vector
  beta[non_zero_indices] <- rnorm(n_nonzero, 0, beta_sd)
  
  # Generate errors
  set.seed(seed.noise)
  epsilon <- rnorm(n, 0, error_sd)
  
  # Generate y vector
  y <- X %*% beta + epsilon
  
  return(list(X=X, y=y, beta_true=beta))
}



simulate_ordinis <- function(n = 500, p = 50000, m = 50, sd_noise = 3){
  beta_true <- matrix(c(runif(m), rep(0, p - m))) 
  X <- matrix(rnorm(n * p, sd = sd_noise), n, p)
  y <- drop(X %*% beta_true) + rnorm(n)
  
  return(list(X = X, y = y, beta_true = beta_true))
}




# # simulate some univariate data with sparse solution for any GLM family ####
# simulate_glm <- function(n = 10000, 
#                          p = 100, 
#                          n_nonzero = 50, 
#                          intercept = 1, 
#                          family = "gaussian",
#                          error_sd = 1,
#                          ntot = 10) {
#   # n: number of observations
#   # p: number of variables (excluding intercept)
#   # n_nonzero: nr of variables (excluding intercept term) with nonzero coefficients
#   # intercept: coefficient for intercept term
#   # family: response distribution ("gaussian", "poisson", "binomial", "gamma", "inverse.gaussian")
#   # error_sd: sd of error for Gaussian case
#   # ntot: total sample size for the binomial case
#   
#   # Create matrix of predictor variables
#   X <- cbind(1, # for intercept 
#              matrix(runif(n * p), nrow = n, ncol = p))
#   
#   # Set true coefficients 
#   beta <- c(intercept, rep(0, p))
#   # if (family=="gamma") { beta[sample(2:(p+1), n_nonzero)] <- abs(rnorm(n_nonzero)) } else { 
#   # sample from half normal distribution for gamma to make sure that true coefficients are always positive 
#   beta[sample(2:(p+1), n_nonzero)] <- rnorm(n_nonzero) 
#   # }
#   
#   # Set link function and inverse link function based on family
#   if (family == "gaussian") {
#     linkinv <- identity
#     linkf <- identity
#   } else if (family == "poisson") {
#     linkinv <- exp
#     linkf <- log
#   } else if (family == "binomial") {
#     linkinv <- plogis
#     linkf <- qlogis
#   } else if (family == "gamma") { # log link Gamma
#     linkinv <- exp
#     linkf <- log
#   } else if (family == "inverse.gaussian") { # log link inverse Gaussian
#     linkinv <- exp # function(eta) 1 / eta^2
#     linkf <- log # function(mu) sqrt(1 / mu)
#   } else {
#     stop("Invalid family")
#   }
#   
#   # Generate response variable
#   eta <- X %*% beta
#   mu <- linkinv(eta)
#   if (family == "gaussian") {
#     y <- rnorm(n, mean = mu, sd = error_sd)
#   } else if (family == "poisson") {
#     y <- rpois(n, lambda = mu)
#   } else if (family == "binomial") {
#     y <- rbinom(n, size = ntot, prob = mu)/ntot
#   } else if (family == "gamma") {
#     # shape = 1 and scale = 1/rate. scale here is mu & rate 1/mu to get the mean mu
#     y <-  rgamma(n, shape = 1, scale = mu)
#   } else if (family == "inverse.gaussian") {
#     y <- statmod::rinvgauss(n, mean = mu, shape = 1)
#   }
#   
#   # Set prior weights based on family
#   if (family == "binomial") {
#     weights <- rep(ntot, n)
#   } else {
#     weights <- rep(1, n)
#   }
#   
#   # L2 norm normalize columns of X
#   # X.n <- sqrt(Matrix::colSums(X^2)) # Calculate L2 norm of columns
#   # X <- as.matrix(X %*% Matrix::Diagonal(x = 1 / X.n))
#   # beta <- beta * X.n
#   
#   list(X = X, y = y, weights = weights, beta_true = beta, eta_true = eta, mu_true = mu)
# }
# 
# # data = simulate_glm(n = 100, p = 1000000, n_nonzero = 50, intercept = 1, family = "gaussian")






