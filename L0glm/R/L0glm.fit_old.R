# TODO: delete this file?


#' @describeIn L0glm
#'
#' fits an L0 penalized regression using the adaptive ridge algorithm (Frommlet
#' et Nuel 2016). The function contains minimal code for fitting the adaptive
#' ridge and is called by the higher level \code{L0glm}.
L0glm.fit.old <- function(X, y,
                      weights = rep(1, length(y)),
                      family = poisson(identity),
                      lambda = 0,
                      start = NULL,
                      lower = rep(-.Machine$double.xmax,ncol(X)), #TODO document
                      upper = rep(.Machine$double.xmax,ncol(X)), #TODO document
                      solver = "eigen",
                      method = ifelse(ncol(X)>=20000, 3L, 5L),
                      group.pen.fn = NULL,
                      block.size = NULL,
                      control.solver = NULL,
                      control.block = NULL,
                      nonnegative = FALSE,
                      post.filter.fn = function(u) return(u),
                      control.l0 = list(maxit = 100, reltol = 1E-4, delta = 1E-5, gamma = 2, warn = FALSE),
                      control.iwls = list(maxit = 1, reltol = 1E-4, thresh = 1E-5, warn = FALSE),
                      verbose = FALSE){
  n <- nrow(X)
  p <- ncol(X)
  
  
  control.l0 <- do.call("control.l0.gen", control.l0)
  control.iwls <- do.call("control.iwls.gen", control.iwls)
  if(is.null(solver) & is.null(control.solver)) stop("Please choose a solver. The available solvers are osqp, glmnet, nnls and eigen.")
  
  # Create limit vector if box constraints are provided as integers, throw warning if provided constraints are not in the correct format
  if(is.null(lower)){
    lower <- rep(-Inf, p)
  } else if(length(lower) == 1 & class(lower) %in% c("numeric", "integer")){
    lower <- rep(lower, p)
  } else if(length(lower) == p & class(lower) %in% c("numeric", "integer")) {
    lower <- lower
  } else {
    stop("Lower limit should be numeric value or a numeric vector of length p.")
  }
  
  if(is.null(upper)){
    upper <- rep(Inf, p)
  } else if(length(upper) == 1 & class(lower) %in% c("numeric", "integer")){
    upper <- rep(upper, p)
  } else if(length(upper) == p & class(lower) %in% c("numeric", "integer")){
    upper <- upper
  } else {
    stop("Upper limit should be numeric value or a numeric vector of length p.")
  }
  
  if(is.null(control.block)) control.block <- do.call("control.block.gen", list(block.size = block.size))
  
  if(length(lambda) == 1) lambda <- rep(lambda, p)
  
  if("matrix" %in% class(y)){ 
    if((ncol(y))>1){
      return(multi.fit <- group.cpglm(X = X,
                                      Y = y,
                                      family = family,
                                      lambda = lambda,
                                      lower = lower,
                                      upper = upper,
                                      nonnegative = nonnegative,
                                      solver = solver,
                                      method = method,
                                      start = start,
                                      block.size = block.size,
                                      control.block = control.block,
                                      control.solver = control.solver,
                                      control.l0 = control.l0,
                                      control.iwls = control.iwls,
                                      post.filter.fn = post.filter.fn,
                                      verbose = verbose))
    } else {
      y <- as.vector(y)
    }
    
  } else {
    y <- as.vector(y)
  }
  
  if (family$family == "gaussian") { # TODO for which families should i do this?
    #fixed.weights <- rep(1, length(y))
    fixed.weights <- 1/(y + 0.1) # TODO talk with Tom
  } else if (family$family == "binomial") {
    fixed.weights <- 1/(y*(1-y))
    fixed.weights[y == 0 | y == 1] <- 1  # avoid division by zero
  } else if (family$family == "poisson") {
    fixed.weights <- 1/(y + 0.1)
  } else if (family$family == "Gamma") {
    fixed.weights <- 1/(y^2 + 0.1)
  } else {
    stop("Unsupported family!")
  }
  
  pen.weights = rep(1, p) # ridge regression penalty weights on individual coefficients (sum to p)
  #lambdas = lambda*pen.weights # ridge penalties on individual coefficients
  good = (!is.na(fixed.weights)&!(fixed.weights==0)&!is.infinite(fixed.weights))
  
  
  # The adaptive ridge loop
  mask <- rep(FALSE, p)
  iter.iwls <- vector()
  beta <- if(is.null(start)) rep(1,p) else start
  nonzero = as.vector(beta!=0) # nonzero coefficients
  
  for (i in 1:control.l0$maxit) {
    # print(i)
    beta0 <- beta # coefs from initialization / previous iteration
    pen.weights <- AR.weights(beta = beta, delta = control.l0$delta, gamma = control.l0$gamma)
    lambda.adap <- lambda * pen.weights # ridge penalty with adaptive weights
    if (i <= control.l0$warmupiters & control.l0$maxit > 1){ #TODO: do warmupiters-1? since the initialization could also be counted as an iteration
      # Fit the model with fixed weights
      beta[nonzero] <- lmridge(X = X[good, nonzero, drop = F], y = y[good],
                               pen.weights = pen.weights[good],
                               weights = fixed.weights[good],
                               lambdas = lambda.adap[nonzero],
                               start = beta[nonzero], # as.vector(crossprod(Xwaug, ywaug) / colSums(Xwaug^2)),
                               lower = lower[nonzero],
                               upper = upper[nonzero],
                               nonnegative = nonnegative,
                               solver = solver,
                               method = method)
      
      beta[abs(beta) < control.iwls$thresh] <- 0
      nonzero = as.vector(beta!=0)
      iter.iwls[i] <- 1
    } else {
      # Fit the model with current adaptive pen.weights and update obs weights (done within iwls loop)
      fit <- glm.iwls(X = X, y = y, i = i, weights = weights, 
                      lower = lower, upper = upper, solver = solver, method = method,
                      family = family, start = beta0, # set to NULL to have LSCG & CG being initialised by diagonal / Jacobi initialisation
                      lambda = lambda.adap, pen.weights = pen.weights,
                      converged.set = mask, nonnegative = nonnegative,
                      control.iwls = control.iwls, control.solver = control.solver, control.block = control.block, verbose = verbose)
      beta <- start <- fit$coefficients  # Note: beta and start is updated because this allows that if only 1 iteration is required, the loop boils down to a non penalized glm
      iter.iwls[i] <- fit$iter.iwls
    }
    
    
    # Control iteration
    reldiff <- 2*abs((beta-beta0)/(beta+beta0))
    mask <- reldiff <= control.l0$reltol # keep coefficients fixed in subsequent iterations if they converged
    mask[is.na(mask)] <- TRUE
    stop.loop <- sum(mask) >= (p-1)
    if(stop.loop & i > control.l0$warmupiters) break # Only break if at glm.iwls was called at least once.
  }
  
  if(control.l0$warn && control.l0$maxit > 1 && !stop.loop)
    warning("Best subset algorithm did not converge, maybe increase control.l0$maxit")
  
  fit$iter.iwls <- iter.iwls
  fit$iter.l0 <- i
  fit$lambda <- lambda
  fit$lambda.w <- lambda.adap/lambda # Converged lambda weights
  fit$lambda.w[is.na(fit$lambda.w)] <- 0 # happens when dividing by 0
  fit$converged.l0 <- stop.loop || control.l0$maxit == 1
  
  # Perform post filtering of the coefficients using user supplied function
  fit$coefficients <- post.filter.fn(fit$coefficients)
  fit$post.filter.fn <- post.filter.fn
  
  return(fit)
}


#' @describeIn L0glm
#'
#' fits a GLM using iteratively weighted least squares (McCullagh and Nelder 1989). Note that
#' \code{lambda} should be a vector of length \code{p} and no lambda tuning is
#' performed.
#'
glm.iwls <- function(X, y, weights, family, i = NULL,
                     solver = "eigen",
                     pen.weights = rep(1, ncol(X)),
                     method = ifelse(ncol(X)>=20000, 3L, 5L),
                     start = NULL,
                     lambda = rep(0, ncol(X)),
                     nonnegative = FALSE,
                     lower = rep(-.Machine$double.xmax,ncol(X)),
                     upper = rep(.Machine$double.xmax,ncol(X)),
                     control.iwls = list(maxit = 1, reltol = 1E-4, thresh = 1E-5, warn = FALSE),
                     control.solver = NULL,
                     control.block = NULL,
                     converged.set = rep(FALSE, ncol(X)),
                     verbose = FALSE){ # TODO match
  
  n <- N <- nrow(X)
  p <- ncol(X)
  
  if (p == 0) { # IN case a NULL model is fit, code taken from 'glm'
    return(list(coefficients = numeric(),
                weights = numeric(),
                prior.weights = weights,
                fitted.values = rep(0, n),
                family = family,
                converged.iwls = TRUE,
                iter.iwls = 0))
  }
  
  # Initialize coefficients
  beta <- if(!is.null(start)){
    start
    betaold <- start
  } else {
    rep(1, ncol(X))
  }
  names(beta) <- colnames(X)
  
  # Initialize fit
  if (!is.null(start)){
    eta <- X %*% beta
    mu <- family$linkinv(eta)
  } else {
    # Initial estimate
    nobs <- n # needed for evaluating the initialization
    etastart <- mustart <- NULL # needed for evaluating the initialization
    eval(family$initialize) # generates mustart
    n <- N
    mu <- mustart
    eta <- family$linkfun(mustart)
  }
  
  # Initialize iteration control variables
  niters.nochange <- converged.set * control.iwls$maxit # nr of iterations over which coefficients showed less than reltol relative change.
  nonzero <- rep(TRUE, ncol(X)) # nonzero coefficients
  
  # Start loop
  for(j in 1:control.iwls$maxit) {
    # Linearize the response and adapt weights
    # The updates use Fisher scoring, which is Newton-Raphson with the expected Hessian, see
    # https://stats.stackexchange.com/questions/344309/why-using-newtons-method-for-logistic-regression-optimization-is-called-iterati
    if(family$link == "identity"){ # simplified formula and faster than formula in else statement
      W <- W0 <- as.vector(weights/family$variance(mu)) # 1/variance (* weights) with identity link, where variance = 1 with Gaussian noise, variance = mu with Poisson noise, variance = mu^2 with Gamma noise
      z <- y
    } else {
      mu.prime <- family$mu.eta(eta)
      z <- as.vector(eta + (y - mu) / mu.prime)
      W <- W0 <- as.vector(weights * mu.prime^2 / family$variance(mu))
    }
    
    # Update partial residuals
    # We work on partial residuals on transformed scale z, ie the scaled model
    # predictions minus the predictions for covariates with nonzero coefficients
    # that already converged. This allows for a dramatic speed up after the 1st iteration.
    z <- z - as.vector(X[, converged.set & nonzero, drop=F] %*% beta[converged.set & nonzero])
    
    # Update the observation subset
    good <- is.finite(W) & W>0
    betaold <- beta
    
    if(is.null(control.block$block.size)) control.block$block.size = p
    if(control.block$block.size < p){
      block.ids <- ceiling((1:p)/control.block$block.size)
      beta <- block.lmridge(X = X, y=y, block.ids = block.ids, lambdas = lambda, weights = W, control.block = control.block,
                            pen.weights = pen.weights, lower = lower, upper = upper, nonnegative = nonnegative, solver = solver, method = method,
                            start = betaold, control.solver = control.solver, verbose = verbose)
    } else {
      if(length(beta[beta != 0]) > 1){
        beta[beta!=0] <- lmridge(X[good,beta != 0, drop = F], y = y[good], lambdas = lambda[beta != 0], weights = W[good],
                                 pen.weights = pen.weights[good], lower = lower[beta != 0], upper = upper[beta != 0], 
                                 nonnegative = nonnegative, solver = solver, method = method, 
                                 start = betaold[beta != 0], 
                                 control.solver = control.solver, verbose = verbose)
      } else {
        beta <- lmridge(X = X[good,,drop=F], y = y[good], pen.weights = pen.weights[good],
                        lambdas = lambda, weights = W[good],lower = lower, upper = upper,
                        nonnegative = nonnegative, solver = solver, method = method, 
                        start = betaold,
                        control.solver = control.solver, verbose = verbose)
        
      }
    }
    
    beta[abs(beta) < control.iwls$thresh] <- 0
    nonzero <- beta != 0
    
    # Update fitted data
    eta <- X[,nonzero,drop=F] %*% beta[nonzero] # + offset if required
    mu <- family$linkinv(eta)
    
    # Update iteration control
    reldiff <- 2*abs((beta-betaold)/(beta+betaold)) # Vector of relative differences in coefficient estimates
    reldiff[is.na(reldiff)] <- 0
    niters.nochange <- niters.nochange + (reldiff <= control.iwls$reltol)*1
    # if coefficients stay within reltol in terms of relative difference over
    # 2 subsequent iterations we consider they converged. If initial estimates
    # are supplied, we look only at previous iteration
    # This is perfromed only when working with identity link
    converged.set <- niters.nochange >= (2 - !is.null(start))
    if (anyNA(beta)) {
      if(control.iwls$warn) warning("Error: fit failed - NaN coefficients produced")
      break
    }
    
    if(sum(beta) == 0) {
      beta <- betaold
      if(control.iwls$warn) warning("Error: fit failed - 0 coefficients produced") # TODO warning or print?
      #if(verbose) print("Error: fit failed - 0 coefficients produced")
      zero.coefs <- TRUE 
      break
    } else {
      zero.coefs <- FALSE
    }
    
    if(all(converged.set)) break
  }
  if(control.iwls$warn && !all(converged.set) && control.iwls$maxit > 1)
    warning("IWLS algorithm did not converge, try increasing control.iwls$maxit") # TODO change warning?
  
  # Return output
  return(list(coefficients = beta,
              weights = W0,
              prior.weights = weights,
              fitted.values = as.vector(mu),
              eta = eta,
              family = family,
              converged.iwls = all(converged.set) || control.iwls$maxit == 1,
              iter.iwls = j,
              zero.coefs = zero.coefs))
}

# TODO Explore standardization of y and X and influence on lambda selection vs theoretical lambda


