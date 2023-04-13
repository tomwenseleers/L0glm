multichannel_l0glm.fit <- function(X, Y, family,
                                   lambda = 0,
                                   start = NULL,
                                   lower = rep(-.Machine$double.xmax,ncol(X)),
                                   upper = rep(.Machine$double.xmax,ncol(X)),
                                   solver = "osqp",
                                   method = NULL,
                                   nonnegative = FALSE,
                                   block = FALSE,
                                   control.solver = NULL,
                                   post.filter.fn = function(u) return(u),
                                   control.l0 = list(maxit = 100, reltol = 1E-4, delta = 1E-5, gamma = 2, warn = FALSE),
                                   control.iwls = list(maxit = 1, reltol = 1E-4, thresh = 1E-5, warn = FALSE),
                                   control.fit = list(maxit = 10, block.size = NULL, tol = 1E-7),
                                   verbose = FALSE){
  p <- ncol(X)
  m <- ncol(Y)
  # Create objects to store betas, weights, fitted values, eta, iterations, convergence
  beta.matrix <- matrix(0, nrow = p, ncol = m)
  weights.matrix <- matrix(0, nrow = p, ncol = m)
  prior.matrix <- matrix(0, nrow = p, ncol = m)
  fitted.matrix <- matrix(0, nrow = p, ncol = m)
  eta.matrix <- matrix(0, nrow = p, ncol = m)
  converged.iwls <- vector("numeric", length = 0)
  iter.iwls <- vector("list", m)
  iter.l0 <- vector("list", m)
  lambda.matrix <- matrix(0, nrow = p, ncol = m)
  lambda.w.matrix <- matrix(0, nrow = p, ncol = m)
  converged.l0 <- vector("numeric", length = 0)
  
  
  # # Order columns of y decreasingly based on sparsity (decreasing)
  # colnames(Y) <- seq(ncol(Y))
  # Y_dec <- Y[, order(colSums(Y != 0), decreasing = TRUE)] # Order by decreasing level of sparsity
  
  # Fit model for each channel
  for (channel in 1:m){
    # Calculate the average of previously fitted coefficients as the starting weights
    if (channel > 1) {
      if(!current_fit$zero.coefs){
        start <- rowMeans(beta.matrix[, 1:(channel - 1), drop = FALSE])
      } else {
        # start <- start
        start <- rowMeans(beta.matrix[, 1:(channel - 1), drop = FALSE])
        print(channel)
      }
    } else {
      start <- start
    }
    
    # Fit the model for the current channel
    current_fit <-  L0glm.fit2(X = X,
                               y = as.vector(Y[, channel]),
                               lambda = lambda,
                               weights = 1 / (as.vector(Y[, channel]) + 0.1),
                               family = family,
                               lower = lower,
                               upper = upper,
                               nonnegative = nonnegative,
                               solver = solver,
                               method = method,
                               start = start,
                               block = block,
                               control.solver = control.solver, 
                               control.l0 = control.l0,
                               control.iwls = control.iwls,
                               control.fit = control.fit,
                               post.filter.fn = post.filter.fn,
                               verbose = verbose)
    
    # Store channel values into respective objects
    beta.matrix[, channel] <- current_fit$coefficients
    weights.matrix[, channel] <- current_fit$weights
    prior.matrix[, channel] <- current_fit$prior.weights
    fitted.matrix[, channel] <- current_fit$fitted.values
    eta.matrix[, channel] <- as.vector(current_fit$eta)
    converged.iwls[channel] <- c(converged.iwls, current_fit$converged_iwls)
    iter.iwls[[channel]] <- current_fit$iter.iwls
    names(iter.iwls)[channel] <- paste("channel", channel)
    iter.l0[[channel]] <- current_fit$iter.l0
    names(iter.l0)[channel] <- paste("channel", channel)
    lambda.matrix <- current_fit$lambda
    lambda.w.matrix <- current_fit$lambda.w
    converged.l0 <- c(converged.l0, current_fit$converged.l0)
    
  }
  
  # # Reorder channels
  # colnames(beta.matrix) <- colnames(y_dec) 
  # true_dec <- true_coefs[, order(colSums(sim$y != 0), decreasing = TRUE)]
  # beta.matrix <- beta_matrix[, order(as.numeric(colnames(beta_matrix)))]
  
  return(list(coefficients = beta.matrix,
              weights = weights.matrix,
              prior.weights = prior.matrix,
              fitted.values = fitted.matrix,
              eta = eta.matrix,
              family = family,
              converged.iwls = converged.iwls,
              iter.iwls = iter.iwls,
              iter.l0 = iter.l0,
              lambda = lambda.matrix,
              lambda.w = lambda.w.matrix,
              converged.l0 = converged.l0,
              post.filter.fn = post.filter.fn))
}