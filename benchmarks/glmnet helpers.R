# SOME GLMNET HELPER FUNCTIONS TO BE ABLE TO TUNE lambda USING IC ####

# Compute -2 log-likelihood from a GLM family given y, mu, weights, and deviance
.min2LL_from_family <- function(y, mu, wt, fam, dev) {
  # Clamp mu into the family domain to avoid NaNs in dev/aic
  eps <- .Machine$double.eps
  fam_name <- fam$family
  if (fam_name %in% c("poisson", "quasipoisson", "Gamma", "inverse.gaussian")) {
    mu <- pmax(mu, eps)
  } else if (fam_name == "binomial") {
    mu <- pmin(pmax(mu, eps), 1 - eps)
  }
  # For binomial, the 'n' argument to aic is the number of trials; for others it's unused.
  n_arg <- if (fam_name == "binomial") wt else rep.int(1, length(wt))
  fam$aic(y = y, n = n_arg, mu = mu, wt = wt, dev = dev)  # returns -2 log L
}

# Deviance (-2 log L up to constant) after optional unpenalized GLM refit
# Intercept handling: set intercept=TRUE to include an intercept in the refit/scoring.
.dev_and_k_for_set <- function(X, y, wt, offset, fam, idx, intercept = FALSE, refit = TRUE) {
  n <- NROW(y)
  # Null model (no predictors)
  if (length(idx) == 0L) {
    # Fit/interpolate the null mean on the link scale:
    # With intercept=TRUE, the null GLM has eta = constant; otherwise eta=0.
    if (refit && intercept) {
      fit0 <- glm(y ~ 1, weights = wt, offset = offset, family = fam)
      mu0  <- fit0$fitted.values
      dev0 <- sum(fam$dev.resids(y = y, mu = mu0, wt = wt))
      m2LL <- .min2LL_from_family(y, mu0, wt, fam, dev0)
    } else {
      eta0 <- if (intercept) fam$linkfun(weighted.mean(y, w = if (is.null(wt)) rep(1, n) else wt))
      else 0
      mu0  <- fam$linkinv(eta0)
      dev0 <- sum(fam$dev.resids(y = y, mu = rep_len(mu0, n), wt = wt))
      m2LL <- .min2LL_from_family(y, rep_len(mu0, n), wt, fam, dev0)
    }
    return(list(min2LL = m2LL, k = 0L))
  }
  
  # With predictors
  Xs <- X[, idx, drop = FALSE]
  if (refit) {
    # Unpenalized GLM on active set
    dat <- data.frame(y = y)
    fit <- try(glm(y ~ Xs, weights = wt, offset = offset,
                   family = fam, x = FALSE, y = FALSE, model = FALSE),
               silent = TRUE)
    if (inherits(fit, "try-error")) {
      # Fallback to one Newton step on penalized eta
      eta <- drop(Xs %*% numeric(NCOL(Xs))) + if (is.null(offset)) 0 else offset
      mu  <- fam$linkinv(eta)
      dev <- sum(fam$dev.resids(y = y, mu = mu, wt = wt))
      m2L <- .min2LL_from_family(y, mu, wt, fam, dev)
      return(list(min2LL = m2L, k = NCOL(Xs)))
    }
    mu  <- fit$fitted.values
    dev <- sum(fam$dev.resids(y = y, mu = mu, wt = wt))
    m2L <- .min2LL_from_family(y, mu, wt, fam, dev)
    # Effective df for IC: number of slope coefficients (exclude intercept)
    k   <- length(coef(fit)) - 1L
    return(list(min2LL = m2L, k = k))
  } else {
    # Use penalized linear predictor from current coefficients
    # (Caller must supply 'beta' if they want the exact penalized eta.)
    stop("refit = FALSE requires supplying the penalized eta/mu; use refit = TRUE for robust scoring.")
  }
}

# Information criteria on the (-2 log L) scale
.IC_values <- function(min2LL, k, n, p_total, 
                       gamma_ebic = max(0, 1 - log(n)/(2*log(p_total))), 
                       Cn_hbic = log(log(n))) {
  c(
    aic  = min2LL + 2 * k,
    bic  = min2LL + log(n) * k,
    ebic = min2LL + (log(n) + 2 * gamma_ebic * log(p_total)) * k,
    hbic = min2LL + (Cn_hbic * log(p_total)) * k,          # deviance-scale HBIC
    gic  = min2LL + (log(p_total) * log(log(n))) * k,
    mbic = min2LL + log(n * p_total^2 / 16) * k
  )
}

# ---------------------------------------------------------------------------
# Main: score a glmnet path using family$dev.resids and family$aic
# Supports gaussian(), binomial(), poisson(); (multinomial/mgaussian/cox are excluded here)
# ---------------------------------------------------------------------------
score_ic_glmnet <- function(fit, X, y,
                                   family = gaussian(),
                                   weights = NULL,
                                   offset  = NULL,
                                   intercept = FALSE,
                                   refit = TRUE, # refit selected features using unpenalised GLM to debias coefs?
                                   p_total = ncol(X),
                                   gamma_ebic = max(0, 1 - log(NROW(y))/(2*log(ncol(X)))),
                                   Cn_hbic = log(log(NROW(y)))) {
  n <- NROW(y); p <- p_total
  if (is.null(weights)) weights <- rep(1, n)
  if (is.null(offset))  offset  <- rep(0, n)
  
  # Choose EBIC gamma (default = Chen–Chen lower bound, rounded up)
  if (is.null(gamma_ebic)) {
    gamma_ebic <- max(0, 1 - log(n) / (2 * log(p)))
  }
  
  # Active sets along the lambda path
  B    <- fit$beta                     # p x L sparse matrix (lasso/elastic-net path)
  L    <- if (inherits(B, "dgCMatrix")) ncol(B) else NCOL(B)
  out  <- matrix(NA_real_, nrow = L, ncol = 6,
                 dimnames = list(NULL, c("aic", "bic", "ebic", "hbic", "gic", "mbic")))
  kvec <- integer(L)
  
  for (j in seq_len(L)) {
    Sj <- if (inherits(B, "dgCMatrix")) which(B[, j] != 0) else which(B[, j] != 0)
    # Compute -2 log L and k for this active set (unpenalized refit)
    tmp <- .dev_and_k_for_set(X, y, wt = weights, offset = offset,
                              fam = family, idx = Sj,
                              intercept = intercept, refit = refit)
    kvec[j] <- tmp$k
    out[j, ] <- .IC_values(tmp$min2LL, tmp$k, n, p, gamma_ebic, Cn_hbic)
  }
  
  data.frame(lambda = fit$lambda, k = kvec, out, check.names = FALSE)
}
