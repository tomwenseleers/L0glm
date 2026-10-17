# Internal helper function with the core modified IRLS L0glm algorithm
#
# The weighted least square step of the standard IRLS Fisher-scoring based GLM
# algorithm is replaced with a weighted least square regression with adaptive
# ridge penalty (as in Frommlet & Nuel 2016), but also accounting for the
# estimated dispersion in the ridge penalty, plus adding provisions for
# optional nonnegativity constraints, and allowing for any exponential family
# error destribution plus associated link functions (canonical or non-canonical).
#
# Fits a generalized linear model with L0 regularization and option for
# nonnegativity constraints. It takes as input parameters the design
# matrix X (here including any intercept), the response variable y, regularization
# parameters, solver options, and control parameters. The result is a fitted
# model with coefficients estimated using the specified regularization method
# and solver. The function is called by the higher level \code{L0glm.fit}.
#
# An optional intercept term is already added by the L0glm.fit function. The
# intercept follows the same nonnegativity constraint as the other coefficients,
# except for binomial and quasibinomial intercepts. Bounds are derived by L0glm.fit.
# 
L0glm.fit.internal = function(X, 
                              y, 
                              weights = rep(1, nrow(X)),
                              family = poisson(identity),
                              start = NULL,
                              dispersion = NULL,
                              offset = NULL,
                              lambda = log(ncol(X))*log(log(nrow(X))), # = GIC criterion 
                              always.include = NULL,
                              normalize = FALSE,
                              solver = "eigen", 
                              method = 1L, 
                              lower = rep(-.Machine$double.xmax, ncol(X)),
                              upper = rep(.Machine$double.xmax, ncol(X)),
                              maxit = 200,
                              prior.mean = NULL,
                              gamma = 2,
                              thresh = 1e-5,
                              tol = 1E-8,
                              control.solver = do.call("control.solver.gen", list(solver = solver, method = method)),
                              verbose = FALSE,
                              hard_threshold = TRUE,
                              preprocessing = NULL) {

  ## 1. sanity checks ---------------------------------------------------------
  n <- nobs <- nrow(X); p <- nvars <- ncol(X)
  
  if (is.null(offset)) {
    offset <- Matrix(0, n, 1, sparse = TRUE)
  }
  
  if (length(weights) != n) stop("length(weights) must equal nrow(X)")
  if (length(offset) != n) stop("length(offset) must equal nrow(X)")
  if (length(lower) != p) stop("length(lower) must equal ncol(X)")
  if (length(upper) != p) stop("length(upper) must equal ncol(X)")
    
  ##  2. dispersion initialisation --------------------------------------------
  if (!inherits(family, "family")) family <- family()
  dispersion <- L0glm_dispersion(dispersion, family)
  estimate_dispersion <- is.null(dispersion) || dispersion < 0
  if (!is.null(dispersion)) dispersion <- abs(dispersion)
  scale_stable <- 0L
  scale_frozen <- !estimate_dispersion
  scale_updates <- 0L
  scale_restarts <- 0L
  converged <- FALSE
  iterations <- 0L

  ##  3. prior mean (vector of length p) --------------------------------------  
  if (is.null(prior.mean)) { prior.mean <- rep(0, p) } else 
    if (length(prior.mean)==1) { prior.mean <- rep(prior.mean, p) } else 
      if (length(prior.mean)!=p) { stop("prior.mean must have length p") }
  always.include <- intersect(as.integer(always.include), seq_len(p))
  if (is.null(preprocessing))
    preprocessing <- L0glm_preprocessing(X, normalize, lower, always.include)
  if (length(start) == 1L) start <- rep(as.numeric(start), p)
  if (!is.null(start)) {
    if (length(start) != p || any(!is.finite(start)))
      stop("start must contain finite coefficients, one per design column.", call. = FALSE)
    allowed <- as.numeric(start) != 0
    allowed[always.include] <- TRUE
    if (any(!allowed)) {
      # Screen before the dispersion pilot as well as the repeated ridge solves.
      args <- mget(names(formals(L0glm.fit.internal)), envir = environment())
      args$X <- X[, allowed, drop = FALSE]
      args$start <- as.numeric(start)[allowed]
      args$prior.mean <- as.numeric(prior.mean)[allowed]
      args$lower <- lower[allowed]
      args$upper <- upper[allowed]
      args$dispersion <- if (estimate_dispersion && !is.null(dispersion)) -dispersion else dispersion
      args$always.include <- match(always.include, which(allowed))
      args$preprocessing <- L0glm_subset_preprocessing(preprocessing, which(allowed))
      fit <- do.call(L0glm.fit.internal, args)
      coefficients <- rep(0, p)
      penalties <- rep(Inf, p)
      coefficients[allowed] <- fit$coefficients
      penalties[allowed] <- fit$lambdas
      fit$coefficients <- coefficients
      fit$lambdas <- penalties
      return(fit)
    }
  }
  
  ## 4. optional normalisation ------------------------------------------------
  scaling <- preprocessing$scale
  # if (!is.null(colnames(X))) { if (colnames(X)[1]=="Intercept") scaling[1] <- 1 }

  ## 5. family helpers --------------------------------------------------------
  if (!is(family, "family")) family = family()

  variance  <- family$variance      # variance in function of the mean u
  linkinv   <- family$linkinv       # inverse link function, to get predictions on response scale u from eta
  linkfun   <- family$linkfun       # link function g(u), to get predictions on link scale eta from u
  mu_eta    <- family$mu.eta        # function of derivative of the predicted values on response scale u with respect to predicted values on link scale eta = d(g-1)/deta=dmu/deta
  # dev_res   <- family$dev.resids  # deviance residuals in function of y, mu and prior weights wt 
  initexpr  <- family$initialize    # initialisation expression to initialise mustart (predictions on response scale u) and n (nr of observations) and m (total nr of trials) for binomial data
  
  ## set mustart (and possibly n & m), requires nobs and nvars to be preset
  etastart <- NULL
  mustart <- NULL
  eval(initexpr)                          # initialise mustart and n and m (total counts) for binomial
  n_positive <- sum(is.finite(weights) & weights > 0)
  # n <- nobs
  
  eta <- linkfun(mustart)                 # g(mu)  on link scale
  mu  <- mustart                          # on response scale
  
  # dev <- sum(dev_res(y, mu, weights))
  # W.ret     <- rep(0, n)                # working weights
  # z.ret     <- rep(0, n)                # working response
  
  ## 6. set starting coefficients ---------------------------------------------
  if (is.null(start)) {
    start <- L0glm_initial_coefficients(X, y, mu, weights, offset, prior.mean,
      family, lambda, always.include, lower, preprocessing)
    start <- clamp(start, lower, upper)
  }
  start <- clamp(start, lower, upper)
  beta <- beta_old <- start
  prediction_old <- if (preprocessing$mode != "none" && length(always.include))
    as.numeric(X %*% beta) else numeric()
  if (is.null(dispersion)) {
    dispersion <- if (p) L0glm_initial_dispersion(X, y, mu, weights, offset, family,
      always.include = always.include, preprocessing = preprocessing) else
      L0glm_pearson_scale(y, linkinv(as.numeric(offset)), weights, family, 0L)
    if (!is.finite(dispersion) || dispersion <= 0) dispersion <- 1
  }
  # names(beta) = colnames(X)
  
  # dev = sum(dev.resids(y, mu, weights)) # deviance = sum of deviance residuals

  # Exact-zero starts are excluded by the active-set mask below;
  # their unsmoothed adaptive penalty weights are infinite.

  ## 7. main IRLS loop with adaptive ridge step -------------------------------
  for (iter in seq_len(maxit)) {
    iterations <- iter
    ## ---------- working response & weights ---------------------------
    var = variance(mu) # variance in function of the mean u
    
    if (family$link == "identity"){ # then use simplified formula
        W.ret <- as.numeric(weights / var) # updated 1/variance (* weights) weights with identity link, where variance = 1 with Gaussian noise, variance = mu, mu^2 and mu^3 with Poisson, Gamma & inverse Gaussian noise
        good <- which(W.ret > 0 & is.finite(W.ret))
        z.ret <- y - offset # working (adjusted) response
      } else {
        gprime <- mu_eta(eta) # derivative of the expected values u with respect to predictions on link scale eta = d(g-1)/deta=dmu/deta
        good <- which(weights != 0 & gprime != 0 & is.finite(weights))
        grad <- y - mu # gradient of log-likelihood with respect to eta
        z.ret <- eta - offset + grad / gprime # working response from the local quadratic approximation around eta
        W.ret <- weights * as.numeric(gprime^2 / var) # updated working weights
      }
    z <- z.ret[good] # subset used for fit
    W <- W.ret[good] # subset of weights used for fit
    
    ## ---------- adaptive ridge penalties ----------------------------

    # penalty weights & penalties for the L0 norm penalty ####
    pen.w <- AR.weights(beta = beta, gamma = gamma, delta = 0, q = 0)

    lambda_tilde <- lambda*dispersion/4
    lambdas <- lambda_tilde * pen.w
    lambdas[beta == 0] <- Inf # excluded coefficients, including when lambda=0
    lambdas[always.include] <- 0

    # lambda_tilde ridge penalties on individual coefficients, we scale these penalties 
    # by (1/4)*dispersion (cf Frommlet & Nuel 2016) so that we 
    # can directly use the lambda penalty values implied by the AIC & BIC criteria 
    # (e.g. 2 & log(n) for AIC & BIC). This matches the orthogonal Gaussian
    # selection threshold for known dispersion; its use for general GLMs is
    # an approximation, not a guarantee of global IC minimization.
    # For binomial, Poisson and negative binomial, dispersion is fixed at 1.
    
    ## ----- update coefs using adaptive ridge regression  (zero or nonzero centered) ---
    nonzero <- beta != 0L          # logical mask
    nonzero[always.include] <- TRUE
    if (verbose) { print(paste("IRLS iteration", iter))
                   print(paste(sum(nonzero), "coefficients nonzero; range:", min(beta), "-", max(beta))) }
    if (!any(nonzero)) {           
      if (verbose) warning("All coefficients zero - nothing to fit")
      mu <- linkinv(as.numeric(offset))
      converged <- TRUE
      break
    }

    if (any(!is.finite(beta))) { 
      conv <- FALSE
      if (verbose) warning(paste("Non-finite coefficiente at iteration", iter))
      break
      }
    
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
    # Relative selection cutoff plus an absolute numerical floor.
    # The final relaxed refit applies neither cutoff.
    if (hard_threshold) {
      magnitude <- abs(beta) * scaling
      cutoff <- max(c(0, magnitude[setdiff(seq_len(p), always.include)])) * thresh
      remove <- (magnitude < cutoff) | magnitude < 1E-20
      remove[always.include] <- FALSE
      beta[remove] <- 0
    }
    nonzero <- beta != 0L
    nonzero[always.include] <- TRUE
    previous_active <- beta_old != 0L
    previous_active[always.include] <- TRUE
    support_changed <- any(nonzero != previous_active)
    if (estimate_dispersion && hard_threshold && support_changed) {
      if (scale_frozen) scale_restarts <- scale_restarts + 1L
      scale_frozen <- FALSE
      scale_stable <- 0L
    }
    
    ## ----- update mu and eta --------------------------------------
    eta <- as.matrix(offset + X[, which(nonzero), drop = FALSE] %*% beta[nonzero]) # predictions on link scale eta
    mu  <- linkinv(eta)    # predictions on response scale u
    if (family$family %in% c("poisson","quasipoisson","Gamma","inverse.gaussian")) mu <- pmax(mu, .Machine$double.eps) # to avoid log(0)
    
    ## ----- bounded Pearson scale update, without another fit ------
    scale_estimable <- sum(nonzero) < 0.7 * n_positive
    if (!scale_frozen && hard_threshold && iter >= 3L && scale_estimable) {
      target <- L0glm_pearson_scale(y, mu, weights, family, sum(nonzero))
      if (is.finite(target) && target > 0) {
        stable <- !support_changed && abs(target / dispersion - 1) < 0.01
        scale_stable <- if (stable) scale_stable + 1L else 0L
        dispersion <- min(2 * dispersion, max(dispersion / 2, target))
        scale_updates <- scale_updates + 1L
        scale_frozen <- scale_stable >= 3L
      }
    }
      
    
    # dev_old = dev
    # dev = sum(dev.resids(y, mu, weights)) # deviance
    
    ## ----- convergence test ---------------------------------------
    change <- if (preprocessing$mode == "none") sqrt(sum((beta - beta_old)^2)) /
      max(1e-12, sqrt(sum(beta^2))) else
        L0glm_coefficient_change(beta, beta_old, scaling, always.include,
          as.numeric(eta) - as.numeric(offset), prediction_old, weights)
    if (change < tol &&
         (scale_frozen || !hard_threshold || lambda == 0 ||
          (iter >= 3L && !scale_estimable)) ) {
      converged <- TRUE
      break
    }
    beta_old <- beta
    if (length(prediction_old)) prediction_old <- as.numeric(eta) - as.numeric(offset)
    
    # dev_rel_change = sqrt(sum((beta - beta_old)^2)) / sqrt(sum(beta^2)) # relative change in L2 norm of coefficient vector
    # if (dev_rel_change < tol) break
    # dev_old = dev
  }

  # browser()
  
  ## 8. return output ------------------------------------------------
  # A relaxed fit holds scale fixed during the solves, then reports its own
  # residual scale. Supplied positive dispersions remain unchanged throughout.
  if (estimate_dispersion && (!hard_threshold || lambda == 0 || !any(beta != 0))) {
    target <- L0glm_pearson_scale(y, mu, weights, family,
                                 sum(beta != 0) + sum(beta[always.include] == 0))
    if (is.finite(target) && target > 0) dispersion <- target
  }
  list(coefficients = as.numeric(beta),
       weights      = as.numeric(W.ret),
       dispersion   = dispersion,
       convergence  = data.frame(iterations = iterations, converged = converged,
         dispersion_updates = scale_updates, dispersion_restarts = scale_restarts),
       lambdas      = as.numeric(lambdas),
       z            = as.numeric(z.ret))
}
