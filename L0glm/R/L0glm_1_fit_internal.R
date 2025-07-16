# Internal helper function with the core modified IRLS L0glm algorithm
#
# The weighted least square step of the standard IRLS Fisher-scoring based GLM
# algorithm is replaced with a weighted least square regression with adaptive
# ridge penalty (as in Frommlet & Nuel 2016), but also accounting for the
# estimated dispersion in the ridge penalty, plus adding provisions for
# nonnegativity or box constraints, and allowing for any exponential family
# error destribution plus associated link functions (canonical or non-canonical).
#
# Fits a generalized linear model with L0 regularization and option for
# nonnegativity and box constraints. It takes as input parameters the design
# matrix X (here including any intercept), the response variable y, regularization
# parameters, solver options, and control parameters. The result is a fitted
# model with coefficients estimated using the specified regularization method
# and solver. The function is called by the higher level \code{L0glm.fit}.
#
# An optional intercept term is already added by the L0glm.fit function. The
# intercept term is always left unpenalized & fit unconstrained.
# 
L0glm.fit.internal = function(X, 
                              y, 
                              weights = rep(1, nrow(X)),
                              family = poisson(identity),
                              start = NULL,
                              dispersion = NULL,
                              offset = NULL,
                              lambda = log(ncol(X))*log(log(nrow(X))), # = GIC criterion 
                              no.pen = c(),
                              normalize = TRUE,
                              solver = "eigen", 
                              method = 7L, 
                              lower = rep(-.Machine$double.xmax, ncol(X)),
                              upper = rep(.Machine$double.xmax, ncol(X)),
                              maxit = 25,
                              warm.start = 1,
                              prior.mean = NULL,
                              gamma = 2,
                              thresh = 1e-6,
                              tol = 1E-10,
                              control.solver = do.call("control.solver.gen", list(solver = solver, method = method)),
                              verbose = FALSE) { 

  # NOTE currently I haven't implemented any SCREENING RULES like
  # http://proceedings.mlr.press/v119/atamturk20a/atamturk20a.pdf
  # or sure independent screening prescreening, as implemented in abess
  # https://cran.r-project.org/web/packages/abess/abess.pdf
  # Implementing these might further speed up the algo.
  
  ## 1. sanity checks ---------------------------------------------------------
  n <- nobs <- nrow(X); p <- nvars <- ncol(X)
  
  if (length(weights) != n) stop("length(weights) must equal nrow(X)")
  if (length(offset) != n) stop("length(offset) must equal nrow(X)")
  if (length(lower) != p) stop("length(lower) must equal ncol(X)")
  if (length(upper) != p) stop("length(upper) must equal ncol(X)")
  
  if (is.null(offset)) {
    offset <- Matrix(0, n, 1, sparse = TRUE)
  }
    
  ##  2. dispersion initialisation --------------------------------------------
  #  * binomial / poisson : fixed at 1
  #  * user supplied      : keep unchanged
  #  * otherwise          : start at 1 and iteratively update
  if (is.list(dispersion)) dispersion <- unlist(dispersion)[1L]
  dispersion_init = dispersion 
  # if NULL dispersion will be estimated (=variance for gaussian), 
  # except for binomial & poisson where it is 1
  if (family$family %in% c("binomial", "poisson")) { dispersion_init <- 1 }
  if (is.null(dispersion)) { dispersion <- 1 }
  residuals = rep(0, n)                   # unweighted working residuals, here initialised at 0 for 1st iteration

  ##  3. prior mean (vector of length p) --------------------------------------  
  if (is.null(prior.mean)) { prior.mean <- rep(0, p) } 
     else if (length(prior.mean)==1) { prior.mean <- rep(prior.mean, p) } 
     else if (length(prior.mean)!=p) { stop("prior.mean must have length p") }
  
  ## 4. optional normalisation ------------------------------------------------
  scaling <- if (normalize) scaleX(X) else rep(1, p)
  if (!is.null(colnames(X))) { if (colnames(X)[1]=="Intercept") scaling[1] <- 1 }

  ## 5. family helpers --------------------------------------------------------
  if (!is(family, "family")) family = family()

  variance  <- family$variance      # variance in function of the mean µ
  linkinv   <- family$linkinv       # inverse link function, to get predictions on response scale µ from η
  linkfun   <- family$linkfun       # link function g(µ), to get predictions on link scale η from µ
  mu_eta    <- family$mu.eta        # function of derivative of the predicted values on response scale µ with respect to predicted values on link scale η = d(g-1)/dη=dμ/dη
  # dev_res   <- family$dev.resids  # deviance residuals in function of y, mu and prior weights wt 
  initexpr  <- family$initialize    # initialisation expression to initialise mustart (predictions on response scale µ) and n (nr of observations) and m (total nr of trials) for binomial data
  
  ## set mustart (and possibly n & m), requires nobs and nvars to be preset
  etastart <- NULL
  mustart <- NULL
  eval(initexpr)  
  
  eta <- linkfun(mustart)                 # g(μ)  on link scale
  mu  <- mustart                          # on response scale
  
  # dev <- sum(dev_res(y, mu, weights))
  # W.ret     <- rep(0, n)                # working weights
  # z.ret     <- rep(0, n)                # working response
  
  ## 6. set starting coefficients ---------------------------------------------
  if (!is.null(start) && all(start == 0)) start <- NULL
  if (is.null(start)) {
    start <- init_coefs(X, y, mu, weights, offset,
                        prior.mean, family, lambda) # component-wise initialisation
    start <- clamp(start, lower, upper)
    start <- warm.start * start +
      (1 - warm.start) * rep(sd(start)^(1/gamma), p)
  }
  beta <- beta_old <- start
  # names(beta) = colnames(X)
  
  # dev = sum(dev.resids(y, mu, weights)) # deviance = sum of deviance residuals

  # TODO: if warm.start=1 & start!=NULL we should kick out variables with 
  # start=0 from the very beginning, as they would get infinite penalties

  ## 7. main IRLS loop with adaptive ridge step -------------------------------
  for (iter in seq_len(maxit)) {
    ## ---------- working response & weights ---------------------------
    var = variance(mu) # variance in function of the mean µ
    
    if (family$link == "identity"){ # then use simplified formula
        W.ret <- as.numeric(weights / var) # updated 1/variance (* weights) weights with identity link, where variance = 1 with Gaussian noise, variance = mu, mu^2 and mu^3 with Poisson, Gamma & inverse Gaussian noise
        good <- which(W.ret > 0 & is.finite(W.ret))
        z.ret <- y - offset # working (adjusted) response
      } else {
        gprime <- mu_eta(eta) # derivative of the expected values µ with respect to predictions on link scale η = d(g-1)/dη=dμ/dη
        good <- which(weights != 0 & gprime != 0 & is.finite(weights))
        grad <- y - mu # gradient of log-likelihood with respect to η = ∂ℓ/∂η
        z.ret <- eta - offset + grad / gprime # working (adjusted) response = linearised quadratic approximation of log-likelihood function ℓ around η
        W.ret <- weights * as.numeric(gprime^2 / var) # updated working weights
      }
    z <- z.ret[good] # subset used for fit
    W <- W.ret[good] # subset of weights used for fit
    
    ## ---------- adaptive ridge penalties ----------------------------

    # penalty weights & penalties for the L0 norm penalty ####
    pen.w <- AR.weights(beta = beta, gamma = gamma, delta = 0, q = 0) # adaptive ridge pentalty weights, cf. Frommlet & Nuel (2016)

    lambda_tilde <- lambda*dispersion/4
    lambda_tilde[no.pen] <- 0 # terms that should be left unpenalized (e.g. intercept)
    lambdas <- lambda_tilde * pen.w
    
    # lambda_tilde ridge penalties on individual coefficients, we scale these penalties 
    # by (1/4)*dispersion (cf Frommlet & Nuel 2016) so that we 
    # can directly use the lambda penalty values implied by the AIC & BIC criteria 
    # (e.g. 2 & log(n) for AIC & BIC) (stricly speaking this was derived for 
    # the Gaussian case with dispersion=sigma^1=1 only, but this is the more
    # generic relation that is correct for all families & any estimated dispersion
    # (for binomial & poisson & negative binomial dispersion is fixed at 1)
    
    ## ----- update coefs using adaptive ridge regression  (zero or nonzero centered) ---
    nonzero <- beta != 0L          # logical mask
    if (!any(nonzero)) {           
      if (verbose) warning("All coefficients zero – nothing to fit")
      break
    }

    if (any(!is.finite(beta))) { 
      conv <- FALSE
      if (verbose) warning(paste("Non-finite coefficiente at iteration", iter))
      break
      }
    
    if (verbose) print(paste("IRLS iteration", iter))
    
    beta_sub <- lmridge(X              = X[good, which(nonzero), drop = FALSE],
                        y              = z, # already subsetted to [good]
                        weights        = W, # already subsetted to [good]
                        lambdas        = lambdas[nonzero],
                        scaling        = scaling[nonzero],
                        start          = beta[nonzero],
                        iteration      = iter,
                        lower          = lower[nonzero],
                        upper          = upper[nonzero],
                        prior.mean     = prior.mean[nonzero],
                        solver         = solver,
                        control.solver = control.solver,
                        method         = method)
    
    beta[nonzero]            <- beta_sub
    beta[(abs(beta) < max(abs(beta))*thresh)|abs(beta) < 1E-20] <- 0
    nonzero <- beta != 0L
    
    ## ----- update mu and eta --------------------------------------
    eta <- as.matrix(offset + X[, which(nonzero), drop = FALSE] %*% beta[nonzero]) # predictions on link scale η
    mu  <- linkinv(eta)    # predictions on response scale µ
    if (family$family %in% c("poisson","quasipoisson","Gamma","inverse.gaussian")) mu <- pmax(mu, .Machine$double.eps) # to avoid log(0)
    
    ## ----- dispersion update (once, if free) ----------------------
    if (is.null(dispersion_init)) {
      residuals_old <- if (family$link == "identity") y - mu else (y - mu) / gprime
      res2          <- residuals_old^2
      n.ok          <- length(good)
      rank          <- sum(nonzero) # cheap approximation of the trace of the hat matrix
      df.res        <- n.ok - min(rank, n.ok)
      disp_new      <- sum(W.ret[good] * res2[good]) / df.res
      
      if (abs(disp_new - dispersion) / dispersion < 0.1) {
        dispersion       <- disp_new
        dispersion_init  <- 1        # lock further updates
      }
    }
    
    # update estimated dispersion for all distributions except binomial & poisson
    # and if relative change in estimated dispersion drops below threshold
    # update dispersion to use it in calculated ridge penalty
    # we only do this update once
    # if a dispersion different than 1 is passed to L0glm.fit we keep it at that value
    if (is.null(dispersion_init)) {
      
      residuals_old = residuals # old unweighted working residuals (on adjusted z scale)
      
      if (family$link == "identity"){ 
        residuals <- y-mu # updated unweighted residuals
      } else {
        residuals <- (y-mu)/gprime # updated unweighted residuals
      }
      
      n.ok <- length(good) # good observations
      rank <- sum(nonzero) # estimated rank - cheap approx of trace of the hat matrix
      df.residual <- n.ok - min(rank, n.ok) + 1 # residual degrees of freedom
      dispersion_old <- sum((W.ret[good] * residuals_old[good]^2))/df.residual # previous iteration's estimated dispersion # TODO sometimes causes Error longer object length is not a multiple of shorter object length
      dispersion_new <- sum((W.ret[good] * residuals[good]^2))/df.residual # updated estimate of dispersion (variance for Gaussian case)
      disp_rel_change <- abs(dispersion_new - dispersion_old) / (dispersion_old+1E-10)
      # start updating dispersion used in calculation of lambda when estimated dispersion 
      # converged (relative change in dispersion between iterations < 10%)
      if (disp_rel_change < 0.1) {
        dispersion <- dispersion_new 
        # dispersion_init <- 1 # if we would like to stop updating dispersion in subsequent iterations, unmark if you would like to keep on updating dispersion
      }
    }
      
    
    # dev_old = dev
    # dev = sum(dev.resids(y, mu, weights)) # deviance
    
    ## ----- convergence test ---------------------------------------
    if ( sqrt(sum((beta - beta_old)^2)) /
         max(1e-12, sqrt(sum(beta^2))) < tol ) break
    beta_old <- beta
    
    # dev_rel_change = sqrt(sum((beta - beta_old)^2)) / sqrt(sum(beta^2)) # relative change in L2 norm of coefficient vector
    # if (dev_rel_change < tol) break
    # dev_old = dev
  }

  # browser()
  
  ## 8. return output ------------------------------------------------
  list(coefficients = as.numeric(beta),
       weights      = as.numeric(W.ret),
       dispersion   = dispersion,
       lambdas      = as.numeric(lambdas),
       z            = as.numeric(z.ret))
}