# Internal helper function with the core modified IRLS L0glm algorithm
#
# The weighted least square step of the standard IRLS Fisher-scoring based
# GLM algorithm is replaced with a weighted least square regression with
# adaptive ridge penalty (as in Frommlet & Nuel 2016, but also estimating
# the correct dispersion multiplier in the ridge penalty, and using
# Fisher scoring instead of Newton-Raphson, plus adding provisions for
# nonnegativity or box constraints).
#
#  Fits a generalized linear model with L0 regularization and option for
#  nonnegativity and box constraints. It takes as input parameters the design
#  matrix X, the response variable y, regularization parameters, solver options,
#  and control parameters. The result is a fitted model with coefficients
#  estimated using the specified regularization method and solver. The function
#  it is called by the higher level \code{L0glm.fit}.
#
L0glm.fit.internal = function(X, 
                              y, 
                              family = poisson(identity),
                              weights = rep(1, nrow(X)),
                              start = NULL,
                              dispersion = NULL,
                              offset = rep(0, nrow(X)),
                              lambda = log(ncol(X))*log(log(nrow(X))), # = GIC criterion 
                              no.pen = 0,
                              nonnegative = FALSE,
                              normalize = TRUE,
                              solver = "eigen", 
                              method = 7L, 
                              lower = rep(-.Machine$double.xmax, ncol(X)),
                              upper = rep(.Machine$double.xmax, ncol(X)),
                              maxit = 25,
                              warm.start = 0.9,
                              prior.mean = NULL,
                              gamma = 2,
                              thresh = 1e-5,
                              tol = 1E-10,
                              control.solver = do.call("control.solver.gen", list(solver = solver, method = method)),
                              verbose = FALSE) { 

  
  # NOTE currently I haven't implemented any SCREENING RULES like
  # http://proceedings.mlr.press/v119/atamturk20a/atamturk20a.pdf
  # or sure independent screening prescreening, as implemented in abess
  # https://cran.r-project.org/web/packages/abess/abess.pdf
  # Implementing these might further speed up the algo.
  
  if (grepl("_noreg",solver)) { solver = gsub("_noreg", "", solver) # use given solver without regularisation
                                lambda = 0
                                maxit = 1 }
  
  # Normalize columns of X to unit L2 norm? If so, calculate scaling
  if(normalize) { scaling = scaleX(X) } else { # L2 norm of columns of X
                  scaling = rep(1, ncol(X)) }
  
  # interpret all zero starting coefficient values as start values still requiring initialisation
  if (!is.null(start)) { if (sum(start==0)==length(start)) start = NULL } 

  dispersion_init = dispersion 
  # if NULL dispersion will be estimated (=variance for gaussian), 
  # except for binomial & poisson where it is 1
  if (family$family %in% c("binomial", "poisson")) { dispersion_init <- 1 }
  if (is.null(dispersion)) { dispersion <- 1 }

  if (!is(family, "family")) family = family()
  variance = family$variance # variance in function of the mean µ
  linkinv = family$linkinv # inverse link function, to get predictions on response scale µ from η
  link = family$linkfun # link function g(µ), to get predictions on link scale η from µ
  init = family$initialize # initialisation function to initialise mustart (predictions on response scale µ) and n (nr of observations) and m (total nr of trials) for binomial data
  dev.resids = family$dev.resids # deviance residuals function in function of y, mu and prior weights wt 
  mu.eta = family$mu.eta # function of derivative of the predicted values on response scale µ with respect to predicted values on link scale η = d(g-1)/dη=dμ/dη
  
  if (!is.function(variance) || !is.function(linkinv) )
    stop("'family' argument seems not to be a valid family object", call. = FALSE)

  etastart = NULL
  mustart = NULL
  n = NULL
  nobs <- nrow(X)    # needed by the initialize expression below
  p <- nvars <- ncol(X)   # needed by the initialize expression below
  
  if(!is.null(prior.mean)){ 
    center = prior.mean # center of nonzero centered ridge regression if prior means specified = center of Gaussian priors
  } else {
    prior.mean = rep(0, p)
    center = NULL # zero centered ridge regression, i.e. assume zero centered Gaussian prior
  }
  
  eval(init) # initialise mustart (predictions on response scale µ) and n (nr of observations) and m (total nr of trials) for binomial data
  
  if (is.null(start)) {  
    eta = as.matrix(link(mustart)) # initialize η = g(µ) = predictions on link scale
    # to initialise using an orthogonal initialisation
    # beta <-  beta_old <- start <- init_coefs(X = X, weights = weights, offset = offset, mu = mu, family = family, lambda = lambda)
    # eta = as.matrix(offset + as.vector(if (p == 1L) X * start else X %*% start)) # eta = expected value on link scale
    # to initialise with constant starting values
    mu = linkinv(eta) # predictions on response scale µ
    
    # option1 for start: equal coefficients
    # start <- start_coefs(X=X, family=family) 
    # option2 for start: coordinate-wise coefficients 
    # = crossprod(Xwaug, zwaug) / colSums(Xwaug ^ 2), 
    # which would be exact if variables were orthogonal
    # gave a more stable solution & converged faster
    start <- init_coefs(X=X, y=y, mu=mustart, weights=weights, offset=offset, prior.mean = prior.mean, family=family, lambda=lambda) 
    start <- clip(coefs=start, lower, upper)
    # start <- warm.start*start + (1-warm.start)*rep(sd(start)^(1/gamma), ncol(X))
    
    beta <-  beta_old <- start
  } else {  
    # eta = as.matrix(link(mustart))
    beta <- beta_old <- start # when starting coefficient values were provided
    # WE ASSUME THAT IF STARTING VALUES ARE GIVEN THEY ARE ALREADY QUITE GOOD AND SO USE
    eta = as.matrix(offset + as.vector(if (p == 1L) X * start else X %*% start)) # initialize η = g(µ) = predictions on link scale
    mu = linkinv(eta) # predictions on response scale µ
  }
  names(beta) = colnames(X)
  
  dev = sum(dev.resids(y, mu, weights)) # deviance = sum of deviance residuals
  residuals = rep(0, nobs) # worker residuals, here initialised at 0 for 1st iteration
  
  # nonzero coefficients 
  # TODO: if warm.start=1 we should set nonzero from the start & kick
  # out variables with start=0, which would get infinite penalties, from the very beginning
  if (is.null(start)) { nonzero = as.vector(beta!=0) } else { nonzero = as.vector(start!=0) } 
  # mask <- rep(FALSE, p)
  
  for (i in 1L:maxit) {
    # print(i)
    # print(dispersion)
    
    var = variance(mu) # variance in function of the mean µ
    
      if (family$link == "identity"){ # simplified formula and faster than formula in else statement
        W.ret <- as.vector(weights / var) # updated 1/variance (* weights) weights with identity link, where variance = 1 with Gaussian noise, variance = mu, mu^2 and mu^3 with Poisson, Gamma & inverse Gaussian noise
        good = which(!is.na(W.ret) & (W.ret!=0) & !is.infinite(W.ret)) # valid rows # (gprime != 0)
        z.ret <- (y - offset) # working (adjusted) response (returned non-subsetted version)
        z <- z.ret[good] # subset used for fit
        W <- W.ret[good] # subset of weights used for fit
      } else {
        gprime = mu.eta(eta) # derivative of the expected values µ with respect to predictions on link scale η = d(g-1)/dη=dμ/dη
        good = which((weights != 0) & (gprime != 0) & !is.infinite(weights) & !is.na(weights)) # valid rows
        gradient = y - mu # gradient of log-likelihood with respect to η = ∂ℓ/∂η
        z.ret <- eta - offset + gradient / gprime # working (adjusted) response = linearised quadratic approximation of log-likelihood function ℓ around η
        z <- z.ret[good] 
        W.ret <- weights * as.vector(gprime^2 / var) # updated working weights = diagonal of expected Hessian = reliability of observations
        W <- W.ret[good] # subset of weights used for fit
      }
    
    
    # UPDATE RIDGE PENALTIES ####
    
    # penalty weights & penalties for the L0 norm penalty ####
    pen.weights <- AR.weights(beta = beta, gamma = gamma, delta = 0, q = 0) # adaptive ridge pentalty weights, cf. Frommlet & Nuel (2016)
    lambdas = (1/4)*lambda*dispersion*pen.weights
    # ridge penalties on individual coefficients, we scale these penalties 
    # by the estimated dispersion & by 1/4 (cf Frommlet & Nuel 2016) so that we 
    # can directly use the lambda penalty values implied by the AIC & BIC criteria 
    # (e.g. 2 & log(n) for AIC & BIC) (stricly speaking this was derived for 
    # the Gaussian case only, but it works in practice for other noise
    # families too)
    
    # NOTE: I think that with an Lq penalty and for any q between 0 and 2 the factor would be
    # L0overLqnorm*4*(1/2^(2-q)) as for q=0 this would give (1/4) as a factor,
    # and for q=2, where the penalty weights would become 1 (ie nonadaptive)
    # the correct factor should be the L0overL2norm to get lambda on the equivalent
    # scale of the L0 penalty - the L0overL2norm (for the solution approximating
    # the true solution) could be estimated from an L0glmfit with
    # q=0
      
    # NOTE: later we could consider also allowing for fused Lq ridge penalties,
    # useful e.g. to enforce smoothness among the recovered coefficients
    # beta_1st_order_diff = c(0, diff(beta, 1))
    # 
    # # lambda_fused_ridge = 100
    # pen.weights.fused_ridge <- AR.weights(beta = beta_1st_order_diff, 
    #                                       delta = max(max(abs(beta_1st_order_diff))/10, 1E-10), 
    #                                       gamma = 2, 
    #                                       q = 1) # ridge pentalty weights, 1 with q=2
    # pen.weights.fused_ridge <- length(pen.weights.fused_ridge) * pen.weights.fused_ridge/sum(abs(pen.weights.fused_ridge))
    # lambdas_fused_ridge = lambda*pen.weights.fused_ridge*sum((beta_1st_order_diff)^2)

    # total ridge penalties ####
    # lambdas = lambdas # + lambdas_fused_ridge
    lambdas[no.pen] = 0 # terms that should be left unpenalized (e.g. intercept)
    
    if (length(good)==0) {
      conv <- FALSE
      if(verbose) print(gettextf("No observations informative at iteration %d",
                                 i), domain = NA)
      break
    }
    
    beta_old = beta

    # browser()
    
    # weighted least square step in IRLS algo, here replaced by (zero or nonzero 
    # centered ridge regression to
    # add automatic regularisation and with possible box constraints
    # coefficient update is based on quadratic approximation of log likelihood
    res <- lmridge(X = X[good, nonzero, drop=F], 
                   y = z,
                   weights = W, # fixed bug: had W[good] but already subsetted
                   lambdas = lambdas[nonzero],
                   nonnegative = nonnegative,
                   scaling = scaling[nonzero],
                   start = beta_old[nonzero], # as.vector(crossprod(Xwaug, ywaug) / colSums(Xwaug^2)),
                   iteration = i, # IRLS iteration passed on as extra info
                   lower = lower[nonzero], 
                   upper = upper[nonzero],
                   prior.mean = center[nonzero],
                   solver = solver,
                   control.solver = control.solver,
                   method = method)
  
    if (sum(abs(res)) != 0) {
      beta[nonzero] <- res
    } 

    beta[nonzero][abs(beta[nonzero]) < thresh] = 0
    nonzero = as.vector(beta!=0) # nonzero coefficients
    if(sum(nonzero) == 0){
      conv <- TRUE
      if(verbose) warning("All coefficients penalized to zero")
      break
    }
    if (any(!is.finite(beta))) {
      conv <- FALSE
      if(verbose) warning(gettextf("Non-finite coefficients at iteration %d", i), domain = NA)
      break
    }
    
    ## calculate updated values of predictions on link & response scale η and µ with new coefs:
    # predictions on link scale η
    eta = as.matrix(X[,nonzero,drop=F] %*% beta[nonzero] + offset) 
    # predictions on response scale µ
    mu = linkinv(eta) # expected value / predictions on response scale

    # update estimated dispersion for all distributions except binomial & poisson
    # and if relative change in estimated dispersion drops below threshold
    # update dispersion to use it in calculated ridge penalty
    # we only do this update once
    # if a dispersion different than 1 is passed to L0glm.fit we keep it at that value
    if (is.null(dispersion_init)) {
      
      residuals_old = residuals # old unweighted working residuals (on adjusted z scale)
      
      if (family$link == "identity"){ # simplified formula and faster than formula in else statement
        residuals <- y-mu # updated unweighted residuals
      } else {
        residuals <- (y-mu)/gprime # updated unweighted residuals
      }
      
      n.ok <- length(good) # good observations
      rank <- sum(nonzero) # estimated rank
      df.residual <- n.ok # residual degrees of freedom
      if (rank<n.ok) { df.residual <- df.residual - rank } # to also allow for high dimensional case
      # print(df.residual)
      dispersion_old <- sum((W.ret[good] * residuals_old[good]^2))/df.residual # previous iteration's estimated dispersion # TODO sometimes causes Error longer object length is not a multiple of shorter object length
      dispersion_new <- sum((W.ret[good] * residuals[good]^2))/df.residual # updated estimate of dispersion (variance for Gaussian case)
      disp_rel_change <- abs(dispersion_new - dispersion_old) / dispersion_old
      # print(i)
      # print(disp_rel_change)
      # start updating dispersion used in calculation of lambda when estimated dispersion 
      # converged (relative change in dispersion between iterations < 10%)
      if (disp_rel_change < 0.1) {
        dispersion <- dispersion_new 
        # dispersion_init <- 1 # to stop updating dispersion in subsequent iterations, unmark if you would like to keep on updating dispersion
      }
      # print(dispersion)
    }
      
    dev_old = dev
    dev = sum(dev.resids(y, mu, weights)) # deviance
    
    dev_rel_change = sqrt(sum((beta - beta_old)^2)) / sqrt(sum(beta^2)) # relative change in L2 norm of coefficient vector
    if (dev_rel_change < tol) break
    dev_old = dev
    beta_old = beta 
  }

  # browser()
  
  return(list(coefficients = as.vector(beta), 
              weights = as.vector(W.ret), # working weights
              # residuals = as.vector(residuals.ret), # working residuals
              dispersion = dispersion, # estimated dispersion
              # pen.weights = as.vector(pen.weights), 
              lambdas = as.vector(lambdas), 
              z = as.vector(z.ret)))
}

