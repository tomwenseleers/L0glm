library(NMF)
library(ggplot2)

# Load L0glm-based NMF function
source("R/NNMF.R")
source("R/L0glm_2_fit.R")
source("R/L0glm_1_fit_internal.R")
source("R/utils.R")
source("R/Lmridge.R")
source("R/control.solver.gen.R")
source("R/RcppExports.R")

# Generate synthetic high-dimensional non-negative data
set.seed(42)
m <- 500  # Number of rows (samples)
n <- 300  # Number of columns (features)
k <- 10   # Rank of factorization

X <- abs(matrix(rnorm(m * n), nrow = m, ncol = n))  # High-dimensional dataset

# Initialize results storage
results <- list()

# Run Multiplicative Update (MU) NMF
mu_time <- system.time({
  nmf_mu <- nmf(X, rank = k, method = "brunet", nrun = 1)
})
results$MU <- list(W = nmf_mu@fit@W, H = nmf_mu@fit@H, time = mu_time["elapsed"])

# Run Alternating Least Squares (ALS) NMF
als_time <- system.time({
  nmf_als <- nmf(X, rank = k, method = "lee", nrun = 1)
})
results$ALS <- list(W = nmf_als@fit@W, H = nmf_als@fit@H, time = als_time["elapsed"])


# Run L0glm-based NMF
l0_time <- system.time({
  nmf_l0 <- L0glm.NMF(X, k = k, family = gaussian(), start = matrix(0.1, k, ncol(X)))
})
results$L0glm <- list(W = nmf_l0$W, H = nmf_l0$H, time = l0_time["elapsed"])

# Compute reconstruction error
calc_reconstruction_error <- function(W, H, X) {
  return(norm(X - W %*% H, "F"))
}

for (method in names(results)) {
  results[[method]]$error <- calc_reconstruction_error(results[[method]]$W, results[[method]]$H, X)
  results[[method]]$sparsity <- sum(results[[method]]$W > 0) + sum(results[[method]]$H > 0)
}

# Convert results to a data frame for visualization
df_results <- data.frame(
  Method = rep(names(results), each = 3),
  Metric = rep(c("Reconstruction Error", "Sparsity", "Time (s)"), times = length(results)),
  Value = unlist(lapply(results, function(x) c(x$error, x$sparsity, x$time)))
)

# Plot results
ggplot(df_results, aes(x = Method, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "NMF Benchmark: L0glm vs. MU vs. ALS", y = "Metric Value", x = "Method") +
  theme_minimal()

