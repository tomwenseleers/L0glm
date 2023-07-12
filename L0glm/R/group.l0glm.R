# TODO: delete this file?

#' @describeIn L0glm
#'
#' Implementation of multichannel fitting
#' It is called by the upper code \code{L0glm.fit2}.
#'
#' @export

group.cpglm <- function(X, Y, family,
                        lambda = 0,
                        start = NULL,
                        lower = rep(-.Machine$double.xmax,ncol(X)),
                        upper = rep(.Machine$double.xmax,ncol(X)),
                        solver = "osqp",
                        method = NULL,
                        nonnegative = FALSE,
                        refit.coefs = 10,
                        block.size = NULL,
                        control.block = control.block,
                        control.solver = NULL,
                        post.filter.fn = function(u) return(u),
                        control.l0 = list(maxit = 100, reltol = 1E-4, delta = 1E-5, gamma = 2, warn = FALSE),
                        control.iwls = list(maxit = 1, reltol = 1E-4, thresh = 1E-5, warn = FALSE),
                         verbose = FALSE){
  p <- ncol(X) # Number of covariates
  n <- nrow(X) # Number of observations
  m <- ncol(Y) # Number of correlated channels
  
  if(solver == "eigen" && method %in% c(1, 2, 7, 8) && verbose) print("The selected solver does not use warm starts, which can result in suboptimal solutions. If the result is not satisfactory, try using a different solver or method with warm start capability.")
  
  
  # Order columns of y decreasingly based on sparsity (decreasing)
  colnames(Y) <- seq(ncol(Y))
  Y.dec <- Y[, order(colSums(Y != 0), decreasing = TRUE)] # Order by decreasing level of sparsity


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



  # Fit model for each channel
  for (channel in 1:m){
    # Calculate the average of previously fitted coefficients as the starting weights
    if (channel > 1) {
      start <- rowMeans(beta.matrix[, 1:(channel - 1), drop = FALSE])
    } else {
      start <- start
    }


    # Using the sum of all coefficients and scaling by channel maximum ==> makes no difference
    # if (channel > 1) {
    #   start <- rowSums(beta.matrix[, 1:(channel - 1), drop = FALSE])
    #   max.val <- max(abs(Y[, channel]))
    #   start <- start / max.val
    # } else {
    #   start <- start
    # }

    # Fit the model for the current channel
    current.fit <-  L0glm.fit(X = X,
                               y = as.vector(Y.dec[, channel]),
                               lambda = lambda,
                               family = family,
                               lower = lower,
                               upper = upper,
                               nonnegative = nonnegative,
                               solver = solver,
                               method = method,
                               start = start,
                               block.size = block.size,
                               control.solver = control.solver,
                               control.l0 = control.l0,
                               control.iwls = control.iwls,
                               post.filter.fn = post.filter.fn,
                               verbose = verbose)

    # Store channel values into respective objects
    beta.matrix[, channel] <- current.fit$coefficients
    weights.matrix[, channel] <- current.fit$weights
    prior.matrix[, channel] <- current.fit$prior.weights
    fitted.matrix[, channel] <- current.fit$fitted.values
    eta.matrix[, channel] <- as.vector(current.fit$eta)
    converged.iwls <- c(converged.iwls, current.fit$converged.iwls)
    iter.iwls[[channel]] <- current.fit$iter.iwls
    names(iter.iwls)[channel] <- paste("channel", channel)
    iter.l0[[channel]] <- current.fit$iter.l0
    names(iter.l0)[channel] <- paste("channel", channel)
    lambda.matrix <- current.fit$lambda
    lambda.w.matrix <- current.fit$lambda.w
    converged.l0 <- c(converged.l0, current.fit$converged.l0)

  }
  
  # Reorder channels
  colnames(beta.matrix) <- colnames(Y.dec)
  beta.matrix <- beta.matrix[, order(as.numeric(colnames(beta.matrix)))]
  
  
  



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
