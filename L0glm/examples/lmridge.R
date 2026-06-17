\dontrun{
# Simulate some data
set.seed(123)
n <- 100
p <- 10
X <- matrix(rnorm(n * p), n, p)
beta_true <- rnorm(p)
y <- X %*% beta_true + rnorm(n)

# Fit ridge regression model with default settings
beta_est <- lmridge(X, y, lambdas = 1)

# Fit ridge regression with nonnegativity constraints using the osqp solver
beta_est_nn <- lmridge(X, y, lambdas = 1, nonnegative = TRUE, solver = "osqp")

# Compare estimated coefficients to true coefficients
plot(beta_true, beta_est, xlab = "True Coefficients", ylab = "Estimated Coefficients",
     main = "Ridge Regression Coefficients")
abline(0, 1, col = "red")
}
