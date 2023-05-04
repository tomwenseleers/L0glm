# TODO: documentation
#' @describeIn L0glm
#'
#' Implementation of multichannel fitting
#' It is called by the upper code \code{L0glm.fit2}.
#'
#' @export

group.cpglm2 <- function(X, Y, family,
                         lambda = 0,
                         start = NULL,
                         lower = rep(-.Machine$double.xmax,ncol(X)),
                         upper = rep(.Machine$double.xmax,ncol(X)),
                         solver = "osqp",
                         method = NULL,
                         nonnegative = FALSE,
                         refit.coefs = 3,
                         block = FALSE,
                         control.solver = NULL,
                         post.filter.fn = function(u) return(u),
                         control.l0 = list(maxit = 100, reltol = 1E-4, delta = 1E-5, gamma = 2, warn = FALSE),
                         control.iwls = list(maxit = 1, reltol = 1E-4, thresh = 1E-5, warn = FALSE),
                         control.fit = list(maxit = 10, block.size = NULL, tol = 1E-7),
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
    current.fit <-  L0glm.fit2(X = X,
                               y = as.vector(Y.dec[, channel]),
                               lambda = lambda,
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
  
  # # Refit first few coefficients
  # for(channel in 1:refit.coefs){
  #     # Calculate the average of previously fitted coefficients as the starting weights
  #     if (channel > 1) {
  #       start <- rowMeans(beta.matrix[, 1:(channel - 1), drop = FALSE])
  #     
  #     # Fit the model for the current channel
  #     current.fit <-  L0glm.fit2(X = X,
  #                                y = as.vector(Y.dec[, channel]),
  #                                lambda = lambda,
  #                                family = family,
  #                                lower = lower,
  #                                upper = upper,
  #                                nonnegative = nonnegative,
  #                                solver = solver,
  #                                method = method,
  #                                start = start,
  #                                block = block,
  #                                control.solver = control.solver,
  #                                control.l0 = control.l0,
  #                                control.iwls = control.iwls,
  #                                control.fit = control.fit,
  #                                post.filter.fn = post.filter.fn,
  #                                verbose = verbose)
  #     
  #     # Store channel values into respective objects
  #     beta.matrix[, channel] <- current.fit$coefficients
  #     weights.matrix[, channel] <- current.fit$weights
  #     prior.matrix[, channel] <- current.fit$prior.weights
  #     fitted.matrix[, channel] <- current.fit$fitted.values
  #     eta.matrix[, channel] <- as.vector(current.fit$eta)
  #     converged.iwls <- c(converged.iwls, current.fit$converged.iwls)
  #     iter.iwls[[channel]] <- current.fit$iter.iwls
  #     names(iter.iwls)[channel] <- paste("channel", channel)
  #     iter.l0[[channel]] <- current.fit$iter.l0
  #     names(iter.l0)[channel] <- paste("channel", channel)
  #     lambda.matrix <- current.fit$lambda
  #     lambda.w.matrix <- current.fit$lambda.w
  #     converged.l0 <- c(converged.l0, current.fit$converged.l0)
  # }

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


# group.cpglm2 <- function(X = X, Y = Y,
#                          weights = matrix(1, nrow = nrow(Y), ncol = ncol(Y)),
#                          family = poisson(identity),
#                          intercept = FALSE,
#                          start = NULL,
#                          nonnegative = FALSE,
#                          lower = rep(-.Machine$double.xmax,ncol(X)),
#                          upper = rep(.Machine$double.xmax,ncol(X)),
#                          lambdas = 0,
#                          solver = "osqp",
#                          method = NULL,
#                          group.pen.fn = function(mat) sqrt(rowSums(mat^2)), # default = L2 norm of the estimated coefficients per component
#                          no.pen = 0,
#                          control.l0 = list(maxit = 100, reltol = 1E-4, delta = 1E-5, gamma = 2, warn = FALSE),
#                          control.iwls = list(maxit = 1, reltol = 1E-4, thresh = 1E-5, warn = FALSE),
#                          control.fit = list(maxit = 10, block.size = NULL, tol = 1E-7),
#                          control.solver = NULL,
#                          cv = "none",
#                          verbose = FALSE
# ){
#   n <- nrow(Y) # nr of obersvations
#   m <- ncol(Y) # nr of correlated channels
#   p <- ncol(X) # nr of covariates
# 
#   # Check grouping penalty function
#   test.fun <- group.pen.fn(matrix(rep(c(-1, 0, 1), nrow = 3, ncol = 5)))
#   if(any(is.na(test.fun) | is.infinite(test.fun)) || !is.vector(test.fun) || length(test.fun) != 3)
#     stop("Invalid function supplied. \nThe function should take a matrix and return a vector of length \'ncol(X)\' with finite and non-NA elements")
# 
#   # Check family argument
#   if (is.character(family)) family <- get(family, mode = "function", envir = parent.frame())
#   if (is.function(family)) family <- family()
#   if (family$link!="identity" | !(family$family %in% c("gaussian","poisson","Gamma","inverse.gaussian","quasipoisson")) )
#     stop("family(link) must be one of: gaussian(identity), poisson(identity), Gamma(identity), inverse.gaussian(identity), quasipoisson(identity)")
#   # TODO: support any link function? also support binomial & quasibinomial
# 
#   # Add an intercept if required
#   if(intercept){
#     X <- cbind(Intercept = 1, X)
#     no.pen <- unique(c(1, no.pen + 1)) # No penalization of the intercept
#   }
# 
#   # Get penalization factors lambdas
# 
#   if(identical(lambdas, 0)) lambdas <- rep(0, m)
#   if(any(lambdas < 0)) stop("Negative lambdas not allowed")
#   if(cv == "none"){
#     if(length(lambdas) != m) stop("You must provide a lambda for every channel, ie length(lambdas) must equal ncol(Y), or choose a cross-validation method")
#   } else {
#     # TODO
#     stop("Not implemented")
#   }
#   lambdas <- rep(lambdas, p + intercept)
#   lambdas[no.pen] <- 0
# 
#   # Fit the L0 penalized nonnegative GLM
#   fit <- group.cpglm.fit(X = X, Y = Y, weights = weights, family = family,
#                          lambdas = lambdas, group.pen.fn = group.pen.fn,
#                          solver = solver, method = method,
#                          start = start, lower = lower, upper = upper, nonnegative = nonnegative,
#                          control.l0 = control.l0,
#                          control.iwls = control.iwls,
#                          control.fit = control.fit,
#                          control.solver = control.solver,
#                          post.filter.fn = post.filter.fn,
#                          verbose = verbose)
# 
#   # TODO compute glm slots: residuals, deviance, etc.
#   # TODO compute effective degrees of freedom
#   # TODO perform inference on estimated coefficients
# 
#   names(fit)[names(fit) == "mu"] <- "fitted.values"
#   return(structure(fit, class = c("cpglm"))) # TODO glm and lm compatibility? Maybe not possible because of multivariate Y
# }
# 
# 
# 
# group.cpglm.fit2 <- function(X = X, Y = Y,
#                              weights = matrix(1, nrow = nrow(Y), ncol = ncol(Y)),
#                              family = poisson(identity),
#                              lambdas = rep(0, ncol(X)),
#                              group.pen.fn = function(mat) sqrt(Matrix::rowSums(mat^2)), # default = L2 norm of the estimated coefficients per component
#                              #group.pen.fn = function(mat) Matrix::rowSums(abs(mat)),# try L1 norm
#                              start = NULL,
#                              lower = rep(-.Machine$double.xmax,ncol(X)),
#                              upper = rep(.Machine$double.xmax,ncol(X)),
#                              solver = "osqp",
#                              method = NULL,
#                              nonnegative = FALSE,
#                              control.l0 = list(maxit = 100, reltol = 1E-4, delta = 1E-5, gamma = 2, warn = FALSE),
#                              control.iwls = list(maxit = 1, reltol = 1E-4, thresh = 1E-5, warn = FALSE),
#                              control.fit = list(maxit = 10, block.size = NULL, tol = 1E-7),
#                              control.solver = NULL,
#                              post.filter.fn = function(u) return(u),
#                              verbose = FALSE){
# 
#   n <- nrow(Y) # nr of obersvations
#   m <- ncol(Y) # nr of correlated channels
#   p <- ncol(X) # nr of covariates
# 
#   if(is.null(group.pen.fn)) group.pen.fn <- function(mat) sqrt(Matrix::rowSums(mat^2))
# 
#   # Check if starting values are valid
#   if(is.null(start)){
#     start <- matrix(1, nrow = p, ncol = m)
#   } else if(is.vector(start) || (is.matrix(start) && ncol(start) == 1)){ #TODO test this case
#     start <- matrix(rep(start, m), nrow = p, ncol = m, byrow = FALSE)
#   } else if(is.matrix(start) && (nrow(start) != p || ncol(start) != m)){
#     warning("Invalid dimensions for starting values. Starting values set to 1.")
#   } else {
#     warning("Invalid input for starting values. Starting values set to 1.")
#   }
# 
# 
#   # Create objects to store betas, weights, fitted values, eta, iterations, convergence
#   beta.matrix <- matrix(1, nrow = p, ncol = m)
#   weights.matrix <- matrix(0, nrow = p, ncol = m)
#   prior.matrix <- matrix(0, nrow = p, ncol = m)
#   fitted.matrix <- matrix(0, nrow = p, ncol = m)
#   eta.matrix <- matrix(0, nrow = p, ncol = m)
#   converged.iwls <- vector("numeric", length = m)
#   iter.l0 <- vector("numeric", length = m)
#   lambda.matrix <- matrix(0, nrow = p, ncol = m)
#   lambda.w.matrix <- matrix(0, nrow = p, ncol = m)
#   converged.l0 <- vector("numeric", length = m)
#   zero.coefs.vec <- vector("numeric", length = m)
# 
# 
#   # Adaptive ridge loop
#   mask <- matrix(0, nrow = p, ncol = m)
#   reldiff <- matrix(0, nrow = p, ncol = m)
#   iter.iwls <- vector("list", length = m)
# 
#   converged <- logical(m) # Store convergence for each channel
# 
#   for (i in 1:control.l0$maxit){
#     beta0 <- beta.matrix
#     beta.group <- group.pen.fn(beta.matrix)
#     lambda.adap <- lambdas%o% AR.weights(beta = beta.group, gamma = control.l0$gamma, delta = control.l0$delta)
#     # Fit the model with current adaptive weights for all channels
#     for(channel in 1:m){
#       if(i == 1){ # Set channel names in first iteration
#         names(iter.iwls)[[channel]] <- paste("channel", channel)
#         names(iter.l0)[channel] <- paste("channel", channel)
#       }
# 
#       if(!converged[channel]){
#         current.fit <- glm.iwls2(X = X,
#                                  y = as.vector(Y[, channel]),
#                                  weights = weights[, channel],
#                                  family = family,
#                                  solver = solver,
#                                  method = method,
#                                  start = start[, channel],
#                                  lambda = as.vector(lambda.adap), #[, channel],
#                                  nonnegative = nonnegative,
#                                  # block = block, # TODO remove block?
#                                  lower = lower,
#                                  upper = upper,
#                                  control.iwls = control.iwls,
#                                  control.fit = control.fit, # TODO remove ?
#                                  control.solver = control.solver,
#                                  converged.set = rep(FALSE, ncol(X)),
#                                  verbose = verbose)
# 
#         # Store channel results for current iteration
#         iter.iwls[[channel]] <- c(iter.iwls[[channel]], current.fit$iter.iwls) # Store iwls iterations
#         beta.matrix[, channel] <- start[,channel] <- current.fit$coefficients # Store coefficients
#         weights.matrix[, channel] <- current.fit$weights
#         prior.matrix[, channel] <- current.fit$prior.weights
#         fitted.matrix[, channel] <- current.fit$fitted.values
#         eta.matrix[, channel] <- as.vector(current.fit$eta)
#         converged.iwls[channel] <- current.fit$converged.iwls
#         iter.l0[channel] <- i
#         lambda.matrix[, channel] <- lambdas
#         lambda.w.matrix[, channel] <- lambda.adap/lambdas # Converged lambda weights
#         lambda.w.matrix[is.na(lambda.w.matrix)] <- 0 # happens when dividing by 0
#         zero.coefs.vec[channel] <- current.fit$zero.coefs
#       }
#     }
# 
#     # Control iteration
#     reldiff <- 2*abs((beta.matrix-beta0)/(beta.matrix+beta0))
#     mask <- reldiff <= control.l0$reltol
#     mask[is.na(mask)] <- TRUE
# 
#     for(j in 1:m){
#       if (sum(mask[,j]) >= p-1) {
#         converged[j] <- TRUE
#       } else {
#         converged[j] <- FALSE
#       }
#     }
#     if(all(converged)) break
#   }
#   if(control.l0$warn && control.l0$maxit > 1 && !all(converged)){
#     warning("Best subset algorithm did not converge, maybe increase control.l0$maxit")
#   }
# 
#   # Perform post filtering of the coefficients using user supplied function
#   beta.matrix <- apply(beta.matrix, MARGIN = 2, FUN = post.filter.fn)
# 
#   return(list(coefficients = beta.matrix,
#               weights = weights.matrix,
#               prior.weights = prior.matrix,
#               fitted.values = fitted.matrix,
#               eta = eta.matrix,
#               family = family,
#               converged.iwls = converged.iwls,
#               iter.iwls = iter.iwls,
#               iter.l0 = iter.l0,
#               converged.l0 = converged,
#               lambda = lambdas,
#               lambda.w = lambda.w.matrix,
#               zero.coefs = zero.coefs.vec,
#               post.filter.fn = post.filter.fn))
# }