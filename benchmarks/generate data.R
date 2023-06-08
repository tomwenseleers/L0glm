# simulate data in high dimensional setting as described in Frommlet & Nuel 2016 p. 9
# (single channel data with orthogonal variables & gaussian noise with SD of 1)

# Function to generate one dataset
generate_dataset <- function(n=100, p=100, k_star=24, beta_sd=0.5, error_sd=1) {
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

# Settings
n <- 100
p_values <- c(100, 250, 500, 1000, 2500, 5000, 10000)
k_star <- 24
beta_sd <- 0.5
error_sd <- 1

# Generate datasets
set.seed(42)  # for reproducibility
datasets <- lapply(p_values, function(p) {
  replicate(1000, generate_dataset(n, p, k_star, beta_sd, error_sd), simplify=FALSE)
})

