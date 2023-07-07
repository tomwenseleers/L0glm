# some data simulation functions ####

# simulate single channel data with orthogonal variables & gaussian noise with SD of 1 ####
# cf. Frommlet & Nuel 2016 p. 9

simulate_orthogonal <- function(n=100, p=100, k_star=24, beta_sd=0.5, error_sd=1) {
  # Generate X matrix
  X <- matrix(rnorm(n*p), n, p)
  
  # # we L2 norm normalize the columns of X first
  # X.n <- sqrt(Matrix::colSums(X^2)) 
  # X <- as.matrix(X %*% Matrix::Diagonal(x = 1 / X.n))

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



# simulate some univariate data with sparse solution for any GLM family ####
simulate_glm <- function(n = 10000, 
                         p = 100, 
                         n_nonzero = 50, 
                         intercept = 1, 
                         family = "gaussian",
                         ntot = 10) {
  # n: number of observations
  # p: number of variables (excluding intercept)
  # n_nonzero: nr of variables (excluding intercept term) with nonzero coefficients
  # intercept: coefficient for intercept term
  # family: response distribution ("gaussian", "poisson", "binomial", "gamma", "inverse.gaussian")
  # ntot: total sample size for the binomial case
  
  require(statmod)
  
  # Create matrix of predictor variables
  X <- cbind(1, # for intercept 
             matrix(runif(n * p), nrow = n, ncol = p))
  
  # Set true coefficients 
  beta <- c(intercept, rep(0, p))
  # if (family=="gamma") { beta[sample(2:(p+1), n_nonzero)] <- abs(rnorm(n_nonzero)) } else { 
    # sample from half normal distribution for gamma to make sure that true coefficients are always positive 
  beta[sample(2:(p+1), n_nonzero)] <- rnorm(n_nonzero) 
    # }
  
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
  } else if (family == "gamma") { # log link Gamma
    linkinv <- exp
    linkf <- log
  } else if (family == "inverse.gaussian") { # log link inverse Gaussian
    linkinv <- exp # function(eta) 1 / eta^2
    linkf <- log # function(mu) sqrt(1 / mu)
  } else {
    stop("Invalid family")
  }
  
  # Generate response variable
  eta <- X %*% beta
  mu <- linkinv(eta)
  if (family == "gaussian") {
    y <- rnorm(n, mean = mu, sd = 1)
  } else if (family == "poisson") {
    y <- rpois(n, lambda = mu)
  } else if (family == "binomial") {
    y <- rbinom(n, size = ntot, prob = mu)/ntot
  } else if (family == "gamma") {
    # shape = 1 and scale = 1/rate. scale here is mu & rate 1/mu to get the mean mu
    y <-  rgamma(n, shape = 1, scale = mu)
  } else if (family == "inverse.gaussian") {
    y <- rinvgauss(n, mean = mu, shape = 1)
  }
  
  # Set prior weights based on family
  if (family == "binomial") {
    weights <- rep(ntot, n)
  } else {
    weights <- rep(1, n)
  }
  
  # L2 norm normalize columns of X
  # X.n <- sqrt(Matrix::colSums(X^2)) # Calculate L2 norm of columns
  # X <- as.matrix(X %*% Matrix::Diagonal(x = 1 / X.n))
  # beta <- beta * X.n
  
  list(X = X, y = y, weights = weights, beta_true = beta, eta_true = eta, mu_true = mu)
}

# data = simulate_glm(n = 100, p = 1000000, n_nonzero = 50, intercept = 1, family = "gaussian")

