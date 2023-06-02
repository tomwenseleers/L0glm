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
                         weights = rep(1, ncol(X)),
                         family = gaussian(identity),
                         offset = rep(0, nrow(X)),
                         start = NULL,
                         lower = rep(0, ncol(X)),
                         upper = rep(Inf, ncol(X)),
                         lambda = 1, 
                         nonnegative = FALSE,
                         parallellize = TRUE,
                         relaxedGroupL0 = TRUE,
                         delta = 1E-10, 
                         gamma = 2,
                         thresh = 1,
                         solver = "eigen", # or "osqp" or "glmnet"
                         method = 7L, # = Cholesky # 3(+4) = LSCG, 5(+6) = CG, 7(+8) = Cholesky  
                         maxit = 25,
                         tol = 1E-10,
                         warmupiters = 4,
                         control.solver = do.call("control.solver.gen", list(solver = solver, method = method)),
                         verbose = FALSE){
  
  y <- as.matrix(y)
  
  if(ncol(y) > 1){
    if(parallellize){
      # Fit parallelly all channels (independent L0)
      # Define the number of cores to use
      n_cores <- round(detectCores()/2) # 8 for me
      
      # Determine the type of cluster based on the operating system
      cluster_type <- if (Sys.info()["sysname"] == "Windows") "PSOCK" else "FORK"
      
      # Create a cluster using the number of cores
      cl <- makeCluster(n_cores, type = cluster_type) # on Mac & Ubuntu type="FORK" is faster
      
      # Export the sim data and the function to the workers
      clusterExport(cl, list("sim", "fun", "L0glm.fit.internal"))
      
      # Load the required libraries and functions on each worker
      # (this only needs to be done once)
      cleval <- clusterEvalQ(cl, {
        require(L0glm)
        devtools::load_all(".") # TODO remove this? this should not be necessary once all callable functions are correctly exported by the L0glm package
      })
      
      # Create a list of column index chunks
      cols <- 1:ncol(y)
      chunks <- split(cols, cut(cols, n_cores))
      # Use clusterApply to apply the function in parallel & combine results
      fit1 <- do.call(cbind, clusterApply(cl, chunks, fun, 
                                          X = X, y = y, lambda=lambda, delta=delta, gamma=gamma, thresh=thresh, solver=solver, 
                                          weights = weights, start = start, lower = lower, upper = upper,
                                          method=method, nonnegative = nonnegative, control.solver = control.solver, warmupiters = warmupiters, verbose = verbose))
      
      fit <- fit1
      if(relaxedGroupL0){
        # Apply relaxed group L0 feature selection
        avgcoefs = rowMeans(fit1)
        variable.indices = which(avgcoefs>0) # TODO remove?
        fit2 = Matrix(0, nrow=nrow(fit1), ncol=ncol(fit1), sparse=TRUE)
        cols <- 1:ncol(y)
        chunks <- split(cols, cut(cols, n_cores))
        fit2[variable.indices,] <- do.call(cbind, clusterApply(cl, chunks, fun, 
                                                               X = X[,variable.indices, drop = F], y = y, lambda = lambda, delta = delta, gamma = gamma, thresh = thresh, solver = solver, 
                                                               weights = NULL, start = as.vector(avgcoefs[variable.indices]), lower = lower, upper = upper,
                                                               method = method, nonnegative = nonnegative, control.solver = control.solver, warmupiters = warmupiters, verbose = verbose))
        fit <- fit2
        
        avgcoefs = rowMeans(fit2)
        variable.indices = which(avgcoefs>0)
        fit3 = Matrix(0, nrow=nrow(fit1), ncol=ncol(fit1), sparse=TRUE)
        cols <- 1:ncol(y)
        chunks <- split(cols, cut(cols, n_cores))
        
        fit3[variable.indices,] <- do.call(cbind, clusterApply(cl, chunks, fun,
                                                               X = X[,variable.indices, drop = F], y = y, lambda = lambda, delta = delta, gamma = gamma, thresh = thresh, solver = solver,
                                                               weights = NULL, start = as.vector(avgcoefs[variable.indices]), lower = lower, upper = upper,
                                                               method = method, nonnegative = nonnegative, control.solver = control.solver, warmupiters = warmupiters, verbose = verbose))
        
        fit <- fit3
        
        # avgcoefs = rowMeans(fit3)
        # variable.indices = which(avgcoefs>0)
        # fit4 = Matrix(0, nrow=nrow(fit1), ncol=ncol(fit1), sparse=TRUE)
        # cols <- 1:ncol(y)
        # chunks <- split(cols, cut(cols, n_cores))
        # 
        # fit4[variable.indices,] <- do.call(cbind, clusterApply(cl, chunks, fun,
        #                                                        X = X[,variable.indices, drop = F], y = y, lambda = 0, delta = delta, gamma = gamma, thresh = thresh, solver = solver,
        #                                                        weights = NULL, start = as.vector(avgcoefs[variable.indices]), lower = lower, upper = upper,
        #                                                        method = method, nonnegative = nonnegative, control.solver = control.solver, warmupiters = warmupiters, verbose = verbose))
        # fit <- fit4
        
        
      } else {
        fit <- fit1
      }
      
      
      
      stopCluster(cl)
    } else { # If parallellize = FALSE
      # Fit all channels sequentially
      fit1 <- do.call(cbind, lapply(1:ncol(y), function(col) {
        L0glm.fit.internal(X = X,
                           y = y[,col,drop=F],
                           lambda = lambda,
                           delta = delta,
                           weights = weights,
                           start = start,
                           lower = lower,
                           upper = upper,
                           gamma = gamma,
                           thresh = thresh,
                           solver = solver,  # TODO check why glmnet & osqp crash or give error. Osqp crashes R, glmnet gives error because sum(Xwaug) == 0 before calling bigGLM in lmridge
                           method = method,  # TODO add the rest of the parameters
                           control.solver = control.solver,
                           warmupiters = warmupiters,
                           verbose = verbose)
      }))
      
      if(relaxedGroupL0){
        # Apply relaxed group L0 feature selection
        # 2nd step using 1/ avg coefficients over channels from previous fit as starting values for another L0 GLM fit
        avgcoefs = rowMeans(fit1)
        variable.indices = which(avgcoefs>0)
        fit2 = Matrix(0, nrow=nrow(fit1), ncol=ncol(fit1), sparse=TRUE)
        
        fit2[variable.indices,] <- do.call(cbind, lapply(1:ncol(y), function(col) {
          L0glm.fit.internal(X=X[,variable.indices,drop=F],
                             y=y[,col,drop=F],
                             weights=1/(as.vector(y[,col,drop=F])+0.1),
                             start=as.vector(avgcoefs[variable.indices]),
                             family = family,
                             solver = solver,
                             method = method,
                             lambda = lambda,
                             control.solver = control.solver,
                             nonnegative = nonnegative,
                             verbose = verbose)
        }))
        
        # 3rd step using 1/ avg coefficients over channels from previous fit as starting values for another L0 GLM fit
        avgcoefs = rowMeans(fit2)
        variable.indices = which(avgcoefs>0)
        fit3 = Matrix(0, nrow=nrow(fit1), ncol=ncol(fit1), sparse=TRUE)
        
        fit3[variable.indices,] <- do.call(cbind, lapply(1:ncol(y), function(col) {
          L0glm.fit.internal(X = X[,variable.indices,drop=F], 
                             y = as.vector(y[,col,drop=F]),
                             weights = 1/(as.vector(y[,col,drop=F])+0.1),
                             start = as.vector(avgcoefs[variable.indices]),
                             family = family, 
                             solver = solver, 
                             method = method, 
                             lambda = lambda, 
                             control.solver = control.solver,
                             nonnegative = nonnegative,
                             verbose = verbose)}))
        fit <- fit3
        
        # 4th step using 1/ avg coefficients over channels from previous fit as starting values for another L0 GLM fit but now using lambda=0, 
        # i.e. using a nonnegative Poisson GLM, using variables selected earlier via iterated L0 GLMs
        avgcoefs = rowMeans(fit3)
        variable.indices = which(avgcoefs>0)
        fit4 = Matrix(0, nrow=nrow(fit1), ncol=ncol(fit1), sparse=TRUE)
        
        fit4[variable.indices,] <- do.call(cbind, lapply(1:ncol(y), function(col) {
          L0glm.fit.internal(X = X[,variable.indices,drop=F],
                             y = as.vector(y[,col,drop=F]),
                             weights = 1/(as.vector(y[,col,drop=F])+0.1),
                             start = as.vector(avgcoefs[variable.indices]),
                             family = family,
                             solver = solver,
                             method = method,
                             lambda = 0,
                             control.solver = control.solver,
                             nonnegative = nonnegative,
                             verbose = verbose)}))
        fit <- fit4
        
      } else {
        fit <- fit1
      }
      
    }
  } else {
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
                              delta = delta, 
                              gamma = gamma,
                              thresh = thresh,
                              solver = solver, 
                              method = method, 
                              maxit = maxit,
                              tol = tol,
                              control.solver = control.solver,
                              warmupiters = warmupiters,
                              verbose = verbose)
  }
  return(fit)
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
                              nonnegative = FALSE,
                              lower = rep(0, ncol(X)),
                              upper = rep(Inf, ncol(X)),
                              lambda = 1, 
                              delta = 1E-10, 
                              gamma = 2,
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
  # if (nvars == 0) { # In case a NULL model is fit, code taken from 'glm'
  #   return(list(coefficients = numeric(),
  #               weights = numeric(),
  #               prior.weights = weights,
  #               fitted.values = rep(0, n),
  #               family = family,
  #               converged.iwls = TRUE,
  #               iter.iwls = 0))
  # }
  
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
      }}
    
    if (all(!good)) {
      conv <- FALSE
      if(verbose) print(gettextf("no observations informative at iteration %d",
                                 i), domain = NA)
      break
    }
    
    pen.weights = AR.weights(beta = beta, delta = delta, gamma = gamma) # adaptive ridge pentalty weights
    lambdas = lambda*pen.weights # ridge penalties on individual coefficients
    
    # weighted least squares ridge regression
    # Xw = Diagonal(x=sqrt(W[good])) %*% X[good,,drop=F] # = X * sqrt(W)
    # Xwaug = rbind(Xw, Diagonal(x=sqrt(lambdas))) # row augmented covariate matrix
    # require(microbenchmark)
    # beta = solve(crossprod(Xw), crossprod(Xw,z*sqrt(W)), tol=2*.Machine$double.eps)
    beta_old = beta

    beta[nonzero] <- lmridge(X = X[good, nonzero, drop=F], y = z,
                             weights = W[good],
                             pen.weights = pen.weights[good], # TODO: remove this argument
                             lambdas = lambdas[nonzero],
                             nonnegative = nonnegative,
                             start = start[nonzero], # as.vector(crossprod(Xwaug, ywaug) / colSums(Xwaug^2)),
                             lower = lower[nonzero], 
                             upper = upper[nonzero],
                             solver = solver,
                             method = method)
    beta[nonzero][abs(beta[nonzero])<thresh] = 0
    nonzero = as.vector(beta!=0) # nonzero coefficients
    if(sum(nonzero) == 0){
      conv <- FALSE
      if(verbose) print("Zero coefficients produced.")
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
  
  return(beta) # TODO add more returns? Or should the function only return the coefs?
  
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
fun <- function(cols, lambda, delta, gamma, thresh, solver, method, warmupiters, nonnegative, X, y, 
                weights = NULL, start = NULL, lower, upper, verbose = FALSE, control.solver = do.call("control.solver.gen", list(solver = solver, method = method))) { 
  # require(L0glm)
  # devtools::load_all(".")
  results <- lapply(cols, function(col) {
    y = as.vector(y[,col,drop=F])
    
    if (is.null(weights) && !is.null(y)) {
      weights = 1/as.vector(y+0.1)
    }
    
    # Call L0glm.fit.internal
    L0glm.fit.internal(X = X,
                       y = y,
                       lambda = lambda,
                       weights = weights,
                       start = start,
                       lower = lower,
                       upper = upper,
                       delta = delta,
                       gamma = gamma,
                       thresh = thresh,
                       solver = solver,
                       nonnegative = nonnegative,
                       method = method,
                       warmupiters = warmupiters,
                       verbose = verbose)
  })
  
  # Return the results as a matrix
  do.call(cbind, results)
}