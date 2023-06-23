#' @describeIn L0glm
#'
#'  Fits a generalized linear model with L0 regularization and option for nonnegativity and box constraints. It takes as input parameters the design matrix X, the response 
#'  variable y, regularization parameters, solver options, and control parameters. It supports parallel computation for independent channels if specified. The function 
#'  applies L0 regularization iteratively, potentially including relaxed group L0 feature selection. The result is a fitted model with coefficients estimated using 
#'  the specified regularization method and solver. The function contains minimal code for fitting the GLM, as it calls \code{L0glm.fit.internal}, and it is called
#'  by the higher level \code{L0glm}. It supports parallellization when fitting a multichannel response variable.
#'
#' @export
L0glm.fit_new = function(X, 
                         y, 
                         weights = NULL,
                         family = gaussian(identity),
                         offset = NULL,
                         start = NULL,
                         lower = rep(-Inf, ncol(X)),
                         upper = rep(Inf, ncol(X)),
                         lambda = 1, 
                         nonnegative = FALSE,
                         parallellize = TRUE,
                         no.pen = 0, # Vector of unpenalized indices. Include intercept if not to be penalized.
                         grouped = TRUE, # TODO: bit unclear what this one does in cobmination with relaxed, maybe single argument group_penalty="independent", "relaxed independent", "group" or "relaxed group" might be clearer?
                         group.maxit = 6, # note: increased default to 6 to increase performance
                         relaxed = TRUE,
                         gamma = 2,
                         nonzero.centered = FALSE,
                         ridge.center = NULL, # Should be a vector of length ncol(X)
                         thresh = 1, # TODO: think about sensible default that should be a function of the scale of the data
                         solver = "eigen", # or "osqp" or "glmnet"
                         method = 7L, # = Cholesky # 3(+4) = LSCG, 5(+6) = CG, 7(+8) = Cholesky  
                         maxit = 25,
                         tol = 1E-10,
                         warmupiters = 4, # TODO change to 1 in gaussian
                         control.solver = do.call("control.solver.gen", list(solver = solver, method = method)),
                         verbose = FALSE){
  
  require(parallel)
  require(qlcMatrix)
  y <- as.matrix(y)
  # Make weights, offset and start into a matrix if they are provided as a vector and y is multichannel
  if(!is.null(weights)){
    if(ncol(y) > 1 & ncol(as.matrix(weights)) == 1){
      weights <- matrix(rep(weights, times = ncol(y)), nrow = length(weights), ncol = ncol(y))
    }
  } else { # If no weights are passed
    if(ncol(y) > 1){ # Multichannel
      weights = matrix(1, nrow = nrow(y), ncol = ncol(y))
    } else { # Single channel
      weights = rep(1, nrow(X))
    }
  }
  
  if(!is.null(offset)){
    if(ncol(y) > 1 & ncol(as.matrix(offset)) == 1){
      offset <- matrix(rep(offset, times = ncol(y)), nrow = length(offset), ncol = ncol(y))
    }
  } else { # If no offset is passed
    if(ncol(y) > 1){ # Multichannel
      offset <- Matrix(0, nrow = nrow(X), ncol=ncol(y), sparse=T)
    } else { # Single channel
      offset = rep(0, nrow(X))
    }
  }
  
  if(!is.null(start)){
    if(ncol(y) > 1 & dim(as.matrix(start))[2] == 1){
      start <- matrix(rep(start, times = ncol(y)), nrow = length(start), ncol = ncol(y))
    }
  }
  
  if(nonnegative) lower[lower < 0] <- 0
  
  
  if(ncol(y) > 1){
    if(!grouped) group.maxit <- 1 # Set group max iteration to 1 if group selection is negative -> channels will be fit independently
    # Define the number of cores to use (if it is 1, don't parallellize)
    n_cores <- round(detectCores()/2) 
    
    if(parallellize & n_cores > 1){
      # Determine the type of cluster based on the operating system
      cluster_type <- if (Sys.info()["sysname"] == "Windows") "PSOCK" else "FORK"
      
      # Create a cluster using the number of cores
      cl <- makeCluster(n_cores, type = cluster_type) # on Mac & Ubuntu type="FORK" is faster
      
      # Export the functions to the workers
      clusterExport(cl, list("fun", "L0glm.fit.internal"))
      
      # Load the required libraries and functions on each worker
      # (this only needs to be done once)
      cleval <- clusterEvalQ(cl, {
        require(L0glm)
        devtools::load_all(".") # TODO remove this? this should not be necessary once all callable functions are correctly exported by the L0glm package
      })
      
      # Create a list of column index chunks
      cols <- 1:ncol(y)
      
      # chunks <- split(cols, cut(cols, n_cores))
      
      chunk.size <- ceiling(length(cols) / n_cores)
      extra <- length(cols) - (chunk.size * (n_cores-1))
      chunk_sizes <- c(rep(chunk.size, n_cores-1), extra)
      split_points <- c(0, cumsum(chunk_sizes))
      split_vector <- rep(1:n_cores, chunk_sizes)
      chunks <- split(cols, split_vector)
      
      for(i in 1L:group.maxit){
        if(i > 1){
          avgcoefs = rowMeans(coefficients)
          variable.indices = which(avgcoefs > 0)
          start <- matrix(rep(avgcoefs, times = ncol(y)), nrow = length(avgcoefs), ncol = ncol(y))
        } else {
          variable.indices <- 1:ncol(X)
        }
        
        coefficients <- pen.weights.mat <- lambdas.mat <- Matrix(0, nrow = ncol(X), ncol = ncol(y), sparse = TRUE)
        weights.mat <- Matrix(0, nrow = nrow(X), ncol = ncol(y), sparse = TRUE)
        
        # Use clusterApply to apply the function in parallel & combine results
        fit <- clusterApply(cl, chunks, fun,
                            X = X[,variable.indices, drop = F], 
                            y = y, 
                            family = family, 
                            lambda = lambda, 
                            gamma = gamma,
                            maxit = maxit,
                            tol = tol,
                            no.pen = no.pen,
                            thresh = thresh, 
                            solver = solver, 
                            nonzero.centered = nonzero.centered,
                            ridge.center = ridge.center[variable.indices],
                            weights = weights, 
                            start = start[variable.indices, , drop = F], 
                            lower = lower[variable.indices], 
                            upper = upper[variable.indices], 
                            offset = offset,
                            method = method, 
                            nonnegative = nonnegative, 
                            control.solver = control.solver, 
                            warmupiters = warmupiters, 
                            verbose = verbose)
        # Extract results
        coefficients[variable.indices,] <- do.call(cbind, lapply(fit, function(x) sapply(x, function(y) y$coefficients)))
        weights.mat <- do.call(cbind, lapply(fit, function(x) sapply(x, function(y) y$weights)))
        pen.weights.mat[variable.indices,] <- do.call(cbind, lapply(fit, function(x) sapply(x, function(y) y$pen.weights)))
        lambdas.mat[variable.indices,] <- do.call(cbind, lapply(fit, function(x) sapply(x, function(y) y$lambdas)))
      }
      
      if(relaxed){
        # Apply relaxed group L0 feature selection
        avgcoefs = rowMeans(coefficients)
        variable.indices = which(avgcoefs > 0)
        start <- matrix(rep(avgcoefs, times = ncol(y)), nrow = length(avgcoefs), ncol = ncol(y))
        
        coefficients <- pen.weights.mat <- lambdas.mat <- Matrix(0, nrow = ncol(X), ncol = ncol(y), sparse = TRUE)
        weights.mat <- Matrix(0, nrow = nrow(X), ncol = ncol(y), sparse = TRUE)
        
        fit <- clusterApply(cl, chunks, fun,
                            X = X[,variable.indices, drop = F], 
                            y = y, 
                            family = family, 
                            lambda = lambda, 
                            gamma = gamma,
                            maxit = maxit,
                            tol = tol,
                            no.pen = no.pen,
                            thresh = thresh, 
                            solver = solver, 
                            nonzero.centered = nonzero.centered,
                            ridge.center = ridge.center[variable.indices],
                            weights = weights, 
                            start = start[variable.indices, , drop = F], 
                            lower = lower[variable.indices], 
                            upper = upper[variable.indices], 
                            offset = offset,
                            method = method, 
                            nonnegative = nonnegative, 
                            control.solver = control.solver, 
                            warmupiters = warmupiters, 
                            verbose = verbose)
        # Extract results
        coefficients[variable.indices,] <- do.call(cbind, lapply(fit, function(x) sapply(x, function(y) y$coefficients)))
        weights.mat <- do.call(cbind, lapply(fit, function(x) sapply(x, function(y) y$weights)))
        pen.weights.mat[variable.indices,] <- do.call(cbind, lapply(fit, function(x) sapply(x, function(y) y$pen.weights)))
        lambdas.mat[variable.indices,] <- do.call(cbind, lapply(fit, function(x) sapply(x, function(y) y$lambdas)))
      }
      
      stopCluster(cl)
      
    } else { # If parallellize = FALSE
      for(i in 1L:group.maxit){
        if(i > 1){
          avgcoefs <- rowMeans(coefficients)
          variable.indices = which(avgcoefs > 0)
          start <- matrix(rep(avgcoefs, times = ncol(y)), nrow = length(avgcoefs), ncol = ncol(y))
        } else {
          variable.indices <- 1:ncol(X)
        }
        
        
        coefficients <- pen.weights.mat <- lambdas.mat <- Matrix(0, nrow = ncol(X), ncol = ncol(y), sparse = TRUE)
        weights.mat <- Matrix(0, nrow = nrow(X), ncol = ncol(y), sparse = TRUE)
        
        fit <- do.call(cbind, lapply(1:ncol(y), function(col) {
          L0glm.fit.internal(X = X[,variable.indices,drop=F],
                             y = y[,col,drop=F],
                             family=family,
                             lambda = lambda,
                             gamma = gamma,
                             nonnegative = nonnegative,
                             offset = as.vector(offset[col, drop = F]),
                             weights = as.vector(weights[,col,drop=F]),
                             nonzero.centered = nonzero.centered,
                             ridge.center = ridge.center[variable.indices],
                             start = as.vector(start[variable.indices, col, drop = F]),
                             lower = lower[variable.indices],
                             upper = upper[variable.indices],
                             maxit = maxit,
                             thresh = thresh,
                             solver = solver,  
                             method = method,  
                             tol = tol,
                             no.pen = no.pen,
                             control.solver = control.solver,
                             warmupiters = warmupiters,
                             verbose = verbose)
        }))
        
        coefficients[variable.indices,] <-  Matrix(do.call(cbind, fit[1, ]), sparse = TRUE)
        weights.mat <-  Matrix(do.call(cbind, fit[2, ]), sparse = TRUE)
        pen.weights.mat[variable.indices,] <-  Matrix(do.call(cbind, fit[3, ]), sparse = TRUE)
        lambdas.mat[variable.indices,] <-  Matrix(do.call(cbind, fit[4, ]), sparse = TRUE)
      }
      
      if(relaxed){
        # Apply relaxed group L0 feature selection
        # using 1/ avg coefficients over channels from previous fit as starting values for another L0 GLM fit but now using lambda=0, 
        # i.e. using a nonnegative Poisson GLM, using variables selected earlier via iterated L0 GLMs
        avgcoefs = rowMeans(coefficients)
        variable.indices = which(avgcoefs>0)
        start <- matrix(rep(avgcoefs, times = ncol(y)), nrow = length(avgcoefs), ncol = ncol(y))
        
        coefficients <- pen.weights.mat <- lambdas.mat <- Matrix(0, nrow = ncol(X), ncol = ncol(y), sparse = TRUE)
        weights.mat <- Matrix(0, nrow = nrow(X), ncol = ncol(y), sparse = TRUE)
        
        fit <- do.call(cbind, lapply(1:ncol(y), function(col) {
          L0glm.fit.internal(X = X[,variable.indices,drop=F],
                             y = as.vector(y[,col,drop=F]),
                             family = family,
                             lambda = 0,
                             gamma = gamma,
                             nonnegative = nonnegative,
                             weights = as.vector(weights[,col,drop=F]),
                             start = as.vector(start[variable.indices, col, drop = F]),
                             lower = lower[variable.indices],
                             upper = upper[variable.indices],
                             maxit = maxit,
                             nonzero.centered = nonzero.centered,
                             ridge.center = ridge.center[variable.indices],
                             thresh = thresh,
                             solver = solver,
                             method = method,
                             tol = tol,
                             no.pen = no.pen,
                             control.solver = control.solver,
                             warmupiters = warmupiters,
                             verbose = verbose)}))
        
        coefficients[variable.indices,] <-  Matrix(do.call(cbind, fit[1, ]), sparse = TRUE)
        weights.mat <-  Matrix(do.call(cbind, fit[2, ]), sparse = TRUE)
        pen.weights.mat[variable.indices,] <-  Matrix(do.call(cbind, fit[3, ]), sparse = TRUE)
        lambdas.mat[variable.indices,] <-  Matrix(do.call(cbind, fit[4, ]), sparse = TRUE)
      }
    }
    weights <- weights.mat
    pen.weights <- pen.weights.mat
    lambdas <- lambdas.mat
  } else { # Single channel
    fit <- L0glm.fit.internal(X = X, 
                              y = y,
                              weights = weights,
                              family = family,
                              nonnegative = nonnegative,
                              offset = offset,
                              start = start,
                              lower = lower,
                              upper = upper,
                              lambda = lambda, 
                              nonzero.centered = nonzero.centered,
                              ridge.center = ridge.center,
                              gamma = gamma,
                              thresh = thresh,
                              solver = solver, 
                              method = method, 
                              maxit = maxit,
                              tol = tol,
                              no.pen = no.pen,
                              control.solver = control.solver,
                              warmupiters = warmupiters,
                              verbose = verbose)
    
    coefficients <- fit$coefficients
    pen.weights <- fit$pen.weights
    lambdas <- fit$lambdas
    weights <- fit$weights
  }
  
  weights[is.infinite(weights)] <- 0
  pen.weights[is.infinite(pen.weights)] <- 0
  lambdas[is.infinite(lambdas)] <- 0
  
  return(list(coefficients = coefficients, weights = weights, pen.weights = pen.weights, lambdas = lambdas))
}

#' @describeIn L0glm
#'
#'  Fits a generalized linear model with L0 regularization and option for nonnegativity and box constraints. It takes as input parameters the design matrix X, the response 
#'  variable y, regularization parameters, solver options, and control parameters. The result is a fitted model with coefficients estimated using 
#'  the specified regularization method and solver. The function it is called by the higher level \code{L0glm.fit_new}.
#'
#' @export
L0glm.fit.internal = function(X, 
                              y, 
                              weights = rep(1, ncol(X)),
                              family = poisson(identity),
                              offset = rep(0, nrow(X)),
                              start = NULL,
                              no.pen = 0,
                              nonnegative = FALSE,
                              lower = rep(0, ncol(X)),
                              upper = rep(Inf, ncol(X)),
                              lambda = 1, 
                              gamma = 2,
                              nonzero.centered = FALSE,
                              ridge.center = NULL,
                              thresh = 1,
                              solver = "eigen", # or "osqp" or "glmnet"
                              method = 7L, # = Cholesky # 3(+4) = LSCG, 5(+6) = CG, 7(+8) = Cholesky  
                              maxit = 25,
                              tol = 1E-10,
                              control.solver = do.call("control.solver.gen", list(solver = solver, method = method)),
                              warmupiters = 4,
                              verbose = FALSE) { # nr of iterations before obs weights are updated in addition to penalty weights in each iteration, 3 or 4 gave best performance
  
  if (!is(family, "family")) family = family()
  variance = family$variance # variance function for given family
  linkinv = family$linkinv # inverse link function
  if (!is.function(variance) || !is.function(linkinv) )
    stop("'family' argument seems not to be a valid family object", call. = FALSE)
  mu.eta = family$mu.eta # function of derivative of the expected values
  etastart = NULL
  n <- nobs <- nrow(X)    # needed by the initialize expression below
  p <- nvars <- ncol(X)   # needed by the initialize expression below
  
  mustart = NULL
  eval(family$initialize) # initializes n and fitted values mustart based on y (mu=expected value on response scale)
  #if (is.null(mustart)) mustart <- y
  
  if (is.null(start)) {  
    beta <- beta_old <- rep(1,p) # initial model coefficients
    eta = family$linkfun(mustart) 
  } else {  # eta = expected value on link scale
    beta <- beta_old <- start # when starting coefficient values were provided
    eta = offset + as.vector(if (p == 1L) X * start else X %*% start) # eta = expected value on link scale
  }
  names(beta) <- colnames(X)
  
  mu = linkinv(eta) # expected value on response scale
  dev.resids = family$dev.resids # deviance residuals function
  dev = sum(dev.resids(y, mu, weights)) # deviance 
  dev_old = 0
  nonzero = as.vector(beta!=0) # nonzero coefficients
  mask <- rep(FALSE, p)
  
  delta = 0
  
  for (i in 1L:maxit) {
    varg = variance(mu) # variance
    gprime = mu.eta(eta) # derivative of the expected values
    
    good = !is.na(weights) & !(weights==0) & !is.infinite(weights) & (gprime != 0) # valid rows
    
    
    if (i==1|i>warmupiters) {
      if (family$link == "identity"){ # simplified formula and faster than formula in else statement
        W = as.vector(weights/variance(mu)) # 1/variance (* weights) with identity link, where variance = 1 with Gaussian noise, variance = mu with Poisson noise, variance = mu^2 with Gamma noise
        good = (!is.na(W) & !(W==0) & !is.infinite(W) & (gprime != 0)) # valid rows
        z = (y - offset)[good] # response
      } else {
        mu.prime <- family$mu.eta(eta) # derivative of the expected values
        good = (weights != 0) & (mu.prime != 0) & !is.infinite(weights) & !is.na(weights)# valid rows
        z <- (eta - offset)[good] + (y - mu)[good] / mu.prime[good] # working (adjusted) response
        W <- weights[good] * as.vector(mu.prime[good]^2 / variance(mu)[good]) # updated weights 
      }
    } else {
      z <- (y - offset)[good]
    }
    
    if(is.null(ridge.center)){ # If it's not null we use user provided
      if(i > warmupiters & nonzero.centered){
        center = start # If start is NULL, this is also set to zero. 
      } else if (nonzero.centered){
        center = beta_old
      } else {
        center = NULL
      }
    } else {
      center = ridge.center
    }
    
    
    if (all(!good)) {
      conv <- FALSE
      if(verbose) print(gettextf("no observations informative at iteration %d",
                                 i), domain = NA)
      break
    }
    
    pen.weights = AR.weights(beta = beta, delta = delta, gamma = gamma) # adaptive ridge pentalty weights
    lambdas = lambda*pen.weights # ridge penalties on individual coefficients
    
    beta_old = beta
    lambdas[no.pen] = 0
    res <- lmridge(X = X[good, nonzero, drop=F], y = z,
                             weights = W[good],
                             pen.weights = pen.weights[good], # TODO: remove this argument
                             lambdas = lambdas[nonzero],
                             nonnegative = nonnegative,
                             start = start[nonzero], # as.vector(crossprod(Xwaug, ywaug) / colSums(Xwaug^2)),
                             lower = lower[nonzero], 
                             upper = upper[nonzero],
                             ridge.center = center[nonzero],
                             solver = solver,
                             method = method)
    if(sum(res) != 0){
      beta[nonzero] <- res
    } else {
      if(verbose) print("Zero coefficients produced. Method is switched to 7.")
      method <- 7L
      lmridge(X = X[good, nonzero, drop=F], y = z,
              weights = W[good],
              pen.weights = pen.weights[good], # TODO: remove this argument
              lambdas = lambdas[nonzero],
              nonnegative = nonnegative,
              start = start[nonzero], # as.vector(crossprod(Xwaug, ywaug) / colSums(Xwaug^2)),
              lower = lower[nonzero], 
              upper = upper[nonzero],
              ridge.center = center[nonzero],
              solver = solver,
              method = method)
    }
    beta[nonzero][abs(beta[nonzero]) < thresh] = 0
    nonzero = as.vector(beta!=0) # nonzero coefficients
    if(sum(nonzero) == 0){
      conv <- FALSE
      if(verbose) warning("Zero coefficients produced.")
      break
    }
    if (any(!is.finite(beta))) {
      conv <- FALSE
      warning(gettextf("non-finite coefficients at iteration %d", i), domain = NA)
      break
    }
    
    ## calculate updated values of eta and mu with the new coefs:
    eta = X[,nonzero,drop=F] %*% beta[nonzero] + offset # expected value / predictions on link scale
    mu = linkinv(eta) # expected value / predictions on response scale
    dev_old = dev
    dev = sum(dev.resids(y, mu, weights))
    
    
    reldiff <- 2*abs((beta-beta_old)/(beta+beta_old))
    mask <- reldiff <= tol # keep coefficients fixed in subsequent iterations if they converged
    mask[is.na(mask)] <- TRUE
    stop.loop <- sum(mask) >= (p-1)
    if(stop.loop & i > warmupiters) break # Only break if at glm.iwls was called at least once.
    
    
    # if (abs(dev - dev_old) / (0.1 + abs(dev)) < tol) break # This gives error because dev goes to infinity
    # or if (norm(beta-beta_old, "2") < tol) break
    dev_old = dev
    beta_old = beta # TODO remove dev since converged is calculated differently?
  }
  
  return(list(coefficients = as.vector(beta), weights = as.vector(W), pen.weights = as.vector(pen.weights), lambdas = as.vector(lambdas)))
  # return(as.vector(beta))
}

#' @describeIn L0glm
#'
#' Helper function used in the \code{L0glm.fit_new} function for parallel computation. It takes a subset of column indices cols along with 
#' other parameters such as regularization parameters, solver options, and data matrices X and y. It performs L0 regularization using 
#' the \code{L0glm.fit.internal} function on each column independently and returns the results as a matrix. The function supports weighting 
#' of observations and nonnegative constraints if specified. 
#'
#'
#' @export
fun <- function(cols, lambda, gamma, thresh, solver, method, warmupiters, nonnegative, X, y, family, nonzero.centered = FALSE, ridge.center = NULL,
                weights = NULL, start = NULL, lower, upper, verbose = FALSE, maxit, tol, offset, no.pen = 0,
                control.solver = do.call("control.solver.gen", list(solver = solver, method = method))) { 
  # require(L0glm)
  # devtools::load_all(".")
  fit <- lapply(cols, function(col) {
    y = as.vector(y[,col,drop=F])
    
    # Call L0glm.fit.internal
    L0glm.fit.internal(X = X,
                       y = y,
                       family = family,
                       lambda = lambda,
                       gamma = gamma,
                       weights = as.vector(weights[,col,drop=F]),
                       start = as.vector(start[, col, drop = F]),
                       maxit = maxit,
                       tol = tol,
                       no.pen = no.pen,
                       lower = lower,
                       offset = as.vector(offset[,col, drop ? F]),
                       upper = upper,
                       nonzero.centered = nonzero.centered,
                       ridge.center = ridge.center,
                       thresh = thresh,
                       solver = solver,
                       nonnegative = nonnegative,
                       method = method,
                       warmupiters = warmupiters,
                       verbose = verbose)
  })
  
  return(fit)
}

