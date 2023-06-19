# some data simulation functions ####

# simulate single channel data with orthogonal variables & gaussian noise with SD of 1 ####
# cf. Frommlet & Nuel 2016 p. 9

simulate_orthogonal <- function(n=100, p=100, k_star=24, beta_sd=0.5, error_sd=1) {
  # Generate X matrix
  X <- matrix(rnorm(n*p), n, p)
  
  # Generate beta vector
  beta <- c(rnorm(k_star, 0, beta_sd), rep(0, p-k_star))
  
  # Generate errors
  epsilon <- rnorm(n, 0, error_sd)
  
  # Generate y vector
  y <- X %*% beta + epsilon
  
  return(list(X=X, y=y, beta=beta))
}

# # Settings
# n <- 100
# p_values <- c(100, 250, 500, 1000, 2500, 5000, 10000)
# k_star <- 24
# beta_sd <- 0.5
# error_sd <- 1
# 
# # Generate datasets
# set.seed(42)  # for reproducibility
# datasets <- lapply(p_values, function(p) {
#   replicate(1000, simulate_orthogonal(n, p, k_star, beta_sd, error_sd), simplify=FALSE)
# })



# simulate some univariate data with sparse solution for any GLM family & canonical link ####
simulate_glm <- function(n, p, n_nonzero, family = "gaussian") {
  # n: number of observations
  # p: number of variables
  # family: response distribution ("gaussian", "poisson", "binomial", "gamma")
  
  # Create matrix of predictor variables
  X <- matrix(runif(n * p), nrow = n, ncol = p)
  
  # Set true coefficients 
  beta <- rep(0, p)
  beta[sample(1:p, n_nonzero)] <- rnorm(n_nonzero)
  
  # Set link function and inverse link function based on family
  if (family == "gaussian") {
    linkinv <- identity
    linkf <- identity
  } else if (family == "poisson") {
    linkinv <- exp
    linkf <- log
  } else if (family == "binomial") {
    linkinv <- plogis
    linkf <- qlogis
  } else if (family == "gamma") {
    linkinv <- function(mu) 1 / mu
    linkf <- function(y) 1 / y
  } else {
    stop("Invalid family")
  }
  
  # Generate response variable
  mu <- X %*% beta
  eta <- linkf(mu)
  if (family == "gaussian") {
    y <- rnorm(n, mean = linkinv(eta), sd = 1)
  } else if (family == "poisson") {
    y <- rpois(n, lambda = linkinv(eta))
  } else if (family == "binomial") {
    y <- rbinom(n, size = 1, prob = linkinv(eta))
  } else if (family == "gamma") {
    y <- rgamma(n, shape = 1, rate = linkinv(eta))
  }
  
  list(X = X, y = y, beta_true = beta)
}

# data = simulate_glm(n = 100, p = 1000000, n_nonzero = 50, family = "gaussian")

