#' L0-Regularized Non-negative Matrix Factorization
#'
#'@param X Input matrix (m x n) with non-negative values.
#'@param k Rank of the factorization (number of components).
#'@param lambda_W L0 penalty for W.
#'@param lambda_H L0 penalty for H.
#'@param max_iter Maximum number of iterations.
#'@param tol Convergence tolerance.
#'@param family GLM family (e.g., "gaussian", "poisson").
#'@return A list with W and H matrices.
#'@examples
#' X <- matrix(abs(rnorm(100)), 10, 10) 
#' result <- L0glm.NMF(X, k = 3)
#'@export

library(Matrix)
library(parallel)

L0glm.NMF <- function(X, k, lambda_W = 0.1, lambda_H = 0.1, max_iter = 100, tol = 1e-4, 
                      family = gaussian(), start = NULL) {
  m <- nrow(X)
  n <- ncol(X)
  
  # Convert X to sparse if applicable
  X <- Matrix(X, sparse = TRUE)
  
  # Initialize W and H
  W <- matrix(runif(m * k), nrow = m, ncol = k)
  H <- matrix(runif(k * n), nrow = k, ncol = n)
  
  # Normalize to prevent divergence
  W <- W / rowSums(W + 1e-10)
  H <- H / colSums(H + 1e-10)
  
  if (is.null(start)) {
    start <- matrix(0.1, k, n)
  }
  
  cl <- makeCluster(detectCores() - 1)  # Parallel computing setup
  clusterExport(cl, c("L0glm.fit", "W", "X", "lambda_H", "family"))
  
  for (iter in 1:max_iter) {
    W_old <- W
    H_old <- H
    
    # Parallel update for H
    H <- t(parSapply(cl, 1:ncol(X), function(j) {
      fit <- L0glm.fit(W, X[, j], family = family, lambda = lambda_H, nonnegative = TRUE)
      if (length(fit$beta) == k) fit$beta else rep(0.1, k)
    }))
    
    # Parallel update for W
    W <- t(parSapply(cl, 1:nrow(X), function(i) {
      fit <- L0glm.fit(t(H), X[i, ], family = family, lambda = lambda_W, nonnegative = TRUE)
      if (length(fit$beta) == k) fit$beta else rep(0.1, k)
    }))
    
    # Check convergence using matrix norms
    diff_W <- norm(W - W_old, "F") / (norm(W_old, "F") + 1e-10)
    diff_H <- norm(H - H_old, "F") / (norm(H_old, "F") + 1e-10)
    
    if (iter %% 10 == 0) {  
      print(paste("Iteration:", iter, "diff_W:", diff_W, "diff_H:", diff_H))
    }
    
    if (max(diff_W, diff_H) < tol) {
      message("Converged at iteration ", iter)
      break
    }
  }
  
  stopCluster(cl)  # Close parallel workers
  
  return(list(W = W, H = H))
}
