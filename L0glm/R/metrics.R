#' @importFrom stats coef cor is.empty.model lm lm.fit logLik model.matrix model.response rnorm sd
#' @importFrom utils data
#' @importFrom qlcMatrix rowMax

# function to calculate fit metrics to tune lambda ####
get.metrics <- function(beta, X, y, weights=matrix(1, nrow = nrow(y), ncol = ncol(y)), 
                        offset = Matrix(0, nrow=nrow(X), ncol=ncol(y), sparse = T), 
                        family = gaussian(identity), 
                        beta_true = NULL, 
                        p = NULL, # original design dimension if a subset was fitted
                        dispersion = NULL,
                        mu_floor = .Machine$double.eps
                        # edf=sum(beta!=0)
) {
  
  if (!is(family, "family")) family = family()
  
  # alwayse coerce beta & beta_true to a sparse dgCMatrix if either is a dense vector
  if (is.vector(beta_true)) beta_true <- Matrix(beta_true, ncol=1, sparse=T)
  if (is.vector(beta)) beta <- Matrix(beta, ncol=1, sparse=T)
  if (!is.null(beta_true) && !issparse(beta_true)) beta_true <- Matrix(beta_true, sparse=T)
  if (!issparse(beta)) beta <- Matrix(beta, sparse=T)
  
  y <- as.matrix(y)
  weights <- matrix(as.numeric(weights), nrow(y), ncol(y))
  offset <- matrix(as.numeric(offset), nrow(y), ncol(y))
  
  # nr of tasks/outcomes, variables & observations
  m = ncol(y) # nr of tasks/outcomes
  # nr of variables, cases & tasks & effective degrees of freedom ####
  if (is.null(p)) p = ncol(X) # nr of variables
  
  # A scale parameter adds to model dimension, not to the predictor pool p.
  dispersion <- L0glm_dispersion(dispersion, family, m)
  free_scale <- if (is.null(dispersion)) rep(TRUE, m) else dispersion < 0
  scale_parameters <- if (family$family %in% c("gaussian", "Gamma", "inverse.gaussian",
                                               "quasibinomial", "quasipoisson")) sum(free_scale) else 0L
  
  n = nrow(X) # nr of cases
  nobs = sum(weights != 0)
  # in glm source code calculated as
  # nobs.glm <- function(object, ...)
  #             if(!is.null(w <- object$prior.weights)) sum(w != 0) else length(object$residuals) # with residuals <- (y - mu)/mu.eta(eta)
  
  # effective degrees of freedom
  edf = sum(beta!=0) # as in abess package, note we assume beta includes intercept if present
  # allow for estimated dispersion, cf https://github.com/SurajGupta/r-source/blob/master/src/library/stats/R/glm.R
  edf <- edf + scale_parameters
  
  # beta[beta<lower] = lower
  # beta[beta>upper] = upper
  
  # metrics related to support recovery ####
  # true support & estimated support & all variable sets
  if (!is.null(beta_true)) vars_true = which(rowMax(abs(beta_true))!=0) else vars_true = NA
  vars_estimated = which(rowMax(abs(beta))!=0)
  
  # support size & ratio of estimated support size over true support size
  if (!is.null(beta_true)) support_size_true = length(vars_true) else support_size_true = NA
  support_size_estimated = length(vars_estimated) 
  if (!is.null(beta_true)) estimated_over_true_support_size = support_size_estimated / support_size_true else estimated_over_true_support_size = NA
  
  # Calculate TP, TN, FP, FN
  if (!is.null(beta_true)) {
    TP = length(intersect(vars_true, vars_estimated)) # Intersection of true and estimated sets
    FP = length(setdiff(vars_estimated, vars_true)) # In estimated but not in true
    TN = ncol(X)-length(union(vars_true, vars_estimated)) # In neither
    FN = length(setdiff(vars_true, vars_estimated)) # In true but not in estimated
    
    # calculate precision, recall, sensitivity, accuracy & F1 score to recover true support from these
    precision = TP / (TP+FP) # = positive predictive value
    recall = TP / (TP+FN) # = sensitivity = TPR
    specificity = TN / (TN+FP) # = TNR
    accuracy = (TP+TN) / (TP+TN+FP+FN)
    F1_score = if (2 * TP + FP + FN > 0) 2 * TP / (2 * TP + FP + FN) else NA_real_
    MCC = support_mcc(TP, FP, TN, FN)
    
    # metrics related to accuracy of coefficients (among all coefficients) ####
    MD_coefs = mean(beta-beta_true) # mean deviation
    MAD_coefs = mean(abs(beta-beta_true)) # mean absolute deviation
    sMAPE_coefs =  mean(2 * abs(beta - beta_true) / (abs(beta) + abs(beta_true) + 1E-10)) # symmetric mean absolute percentage error, 0 = perfect agreement
    bias_coefs = as.vector(coef(lm(as.vector(beta_true) ~ 0 + as.vector(beta)))[[1]]) # regression beta true on beta (closer to 1 is better, >1 means you underestimate true coefs, <1 means you overestimate true coefs)
    R2_coefs = summary(lm(as.vector(beta_true) ~ 0+as.vector(beta)))$r.squared # I use R2 from zero intercept regression of beta_true on beta
    
    # metrics related to accuracy of coefficients (among selected coefficients) ####
    sel = beta!=0
    MD_selcoefs = mean(beta[sel]-beta_true[sel]) # mean deviation
    MAD_selcoefs = mean(abs(beta[sel]-beta_true[sel])) # mean absolute deviation
    sMAPE_selcoefs =  mean(2 * abs(beta[sel] - beta_true[sel]) / (abs(beta[sel]) + abs(beta_true[sel]) + 1E-10)) # symmetric mean absolute percentage error, 0 = perfect agreement
    bias_selcoefs = ifelse(sum(sel)!=0, as.vector(coef(lm(as.vector(beta_true[sel]) ~ 0 + as.vector(beta[sel])))[[1]]), 0) # regression beta true on beta (closer to 1 is better, >1 means you underestimate true coefs, <1 means you overestimate true coefs)
    R2_selcoefs = ifelse(sum(sel)!=0, summary(lm(as.vector(beta_true[sel]) ~ 0+as.vector(beta[sel])))$r.squared, 0) # I use R2 from zero intercept regression of beta_true on beta
    
  } else {
    TP <- FP <- FN <- TN <- precision <- recall <- sensitivity <- accuracy <- specificity <- F1_score <- MCC <- bias_coefs <- MD_coefs <- MAD_coefs <- sMAPE_coefs <- R2_coefs <- NA
    MD_selcoefs <- MAD_selcoefs <- sMAPE_selcoefs <- bias_selcoefs <- R2_selcoefs <- NA_real_
  }
  
  # metrics related to prediction accuracy ####
  # calculate deviance
  if (is.null(beta_true)) { # calculated weights based on expected value
    y = as.matrix(y)
    # Compute the fitted mean values using the fitted coefficients
    nonzero = which(rowMax(abs(beta))!=0) # nonzero coefficients
    eta = as.matrix(X[,nonzero,drop=F] %*% beta[nonzero,,drop=F] + offset) # fitted eta = fitted values on link scale
    mu = family$linkinv(eta=eta) # fitted values on response scale
    if (family$family %in% c("poisson","quasipoisson","Gamma","inverse.gaussian")) mu <- pmax(mu, mu_floor)
    # Compute the weights at the fitted values
    gprime = family$mu.eta(eta=eta) # derivative of the expected values
    varg = family$variance(mu=mu) # variance
    W = as.vector(weights)*(gprime^2)/varg # weights
    good = !is.na(W) & W != 0 & !is.infinite(W) & gprime!=0 & as.vector(weights)!=0 # valid observations
    W[!good] = 0
    dim(W) = dim(eta)
    mu_true = rep(NA, length(mu))
    deviance_true = NA
    
  } else { # calculate weights based on the "true" responses if beta_true is given
    # calculate true response used in simulation given true coefficients beta_true
    nonzero = which(rowMax(abs(beta_true))!=0) # nonzero coefficients
    eta_true = as.matrix(X[,nonzero,drop=F] %*% beta_true[nonzero,,drop=F] + offset) # true eta = true values on link scale
    mu_true <- as.matrix(family$linkinv(eta=eta_true)) # true mu = true values on response scale
    
    # Compute the weights based on the "true" responses
    gprime_true = family$mu.eta(eta=eta_true) # derivative of the true values
    varg_true = family$variance(mu=mu_true)
    W = as.vector(weights)*(gprime_true^2)/varg_true
    good = is.finite(W) & W != 0 & as.vector(weights)!=0 # !is.na(W) & W != 0 & !is.infinite(W) & gprime_true!=0 & as.vector(weights)!=0
    W[!good] = 0
    dim(W) = dim(eta_true)
    
    # Compute the fitted mean values using the actual fitted coefficients
    nonzero = which(rowMax(abs(beta))!=0) # nonzero coefficients
    eta = as.matrix(X[,nonzero,drop=F] %*% beta[nonzero,,drop=F] + offset) # fitted eta = fitted values on link scale
    mu = family$linkinv(eta=eta) # fitted values on response scale
    if (family$family %in% c("poisson","quasipoisson","Gamma","inverse.gaussian")) mu <- pmax(mu, mu_floor)

    # calculated true theoretical deviance compared to ground truth
    dev_res_true = family$dev.resids(y=as.vector(mu_true), mu=as.vector(mu), wt=as.vector(weights))
    keep = which(as.vector(weights) != 0)
    deviance_true = sum(dev_res_true[keep]) # deviance compared to ground truth = true deviance = how well the model fits the latent noise-free mean, cf McCullagh, P. & Nelder, J. (1989) Generalized Linear Models, section 2.3, section 12.2.
  }
  
  
  # Some inverse links (notably logit) return a vector even for matrix input.
  dim(mu) <- c(n, m)
  # calculated observed residual deviance
  dev_res = family$dev.resids(y=as.vector(y), mu=as.vector(mu), wt=as.vector(weights))
  keep = which(as.vector(weights) != 0)
  deviance = sum(dev_res[keep]) # observed residual deviance = how good model fits the realised noisy data
  
  # calculate some information criteria, cf. Zhang & Chen 2010 ####
  # see also https://github.com/abess-team/abess/blob/master/R-package/R/deviance.abess.R
  # get -2 * log-likelihood & log likelihood from family$aic slot & observed deviance
  # Score each response channel at its own scale. Positive supplied dispersions
  # are known; otherwise Gaussian likelihood profiles phi = weighted RSS / N.
  # family$aic includes an extra 2 for estimated scale in the continuous families.
  min2LL <- sum(vapply(seq_len(m), function(j) {
    keep_j <- weights[, j] > 0 & is.finite(weights[, j])
    yy <- y[keep_j, j]; mm <- mu[keep_j, j]; ww <- weights[keep_j, j]
    if (!length(yy)) return(0)
    phi <- if (is.null(dispersion) || free_scale[j]) NULL else dispersion[j]
    if (family$family == "gaussian") {
      rss <- sum(ww * (yy - mm)^2)
      if (is.null(phi)) phi <- rss / length(yy)
      if (phi == 0) return(-Inf) # unbounded likelihood at an exact interpolant
      return(length(yy) * log(2 * pi * phi) - sum(log(ww)) + rss / phi)
    }
    if (family$family == "poisson") yy <- round(yy)
    if (family$family == "Gamma" && !is.null(phi))
      return(-2 * sum(ww * stats::dgamma(yy, shape = 1 / phi, scale = mm * phi, log = TRUE)))
    if (family$family == "inverse.gaussian" && !is.null(phi))
      return(sum(ww * (log(2 * pi * phi * yy^3) + (yy - mm)^2 / (phi * yy * mm^2))))
    dev <- sum(family$dev.resids(yy, mm, ww))
    value <- family$aic(y = yy, n = ww, mu = mm, wt = ww, dev = dev)
    if (family$family %in% c("Gamma", "inverse.gaussian")) value <- value - 2
    value
  }, numeric(1)))
  loglik = -0.5 * min2LL # plain log-likelihood
  IC = get.IC(loglik = loglik, n = n, p = p, edf = edf) 
  
  # weighted_rmse = sqrt(sum((predicted-actual)^2*weight)/sum(weight))
  # z = (eta - offset) + (y - mu) / mu.prime # working (adjusted) response
  # z.res <- (y - mu)/fam$mu.eta(eta) # the residuals on the adjusted scale 
  # = z.res = z - z.fit = eta + (y - mu)/g'(eta) - eta = (y - mu)/g'(eta)
  
  # TODO: still add explained variance or adjusted explained variance on adjusted z scale?
  # # define functions to calculate weighted R2 and adjusted weighted R2 :
  # get.R2 = function (z_obs, # observed data on adjusted response scale (or true data without noise on z scale if beta_true given)
  #                    z_pred, # fitted data on adjusted response scale
  #                    weights, # weights from last IRLS GLM iteration, e.g. 1/variance for identity link GLMs (or true weights if beta_true given)
  #                    intercept=FALSE) { 
  #   good = as.vector(weights<1) # to remove infinite & unreliable weights
  #   ss_residual = sum(weights[good] * (z_obs[good] - z_pred[good]) ^ 2)
  #   if (intercept) { wmean = weighted.mean(z_obs[good], weights[good]) } else { wmean = 0 }
  #   ss_total = sum(weights[good] * (z_obs[good] - wmean) ^ 2)
  #   r.squared = 1 - ss_residual / ss_total 
  #   return(r.squared) }
  # 
  # get.adjR2 = function (z_obs, # observed data on adjusted response scale
  #                       z_pred, # fitted data on adjusted response scale
  #                       weights, # weights from last IRLS GLM iteration, e.g. 1/variance for identity link GLMs
  #                       intercept=FALSE,
  #                       p # nr of variables in the model (nr of nonzero coefficients) (incl intercept)
  # ) { 
  #   good = as.vector(weights<1) # to remove infinite & unreliable weights
  #   n = sum(good) # nr of observations
  #   ss_residual = sum(weights[good] * (z_obs[good] - z_pred[good]) ^ 2)
  #   if (intercept) { wmean = weighted.mean(z_obs[good], weights[good]) } else { wmean = 0 }
  #   ss_total = sum(weights[good] * (z_obs[good] - weighted.mean(z_obs[good], weights[good])) ^ 2)
  #   df_residual = n-p-1
  #   df_total = n-1
  #   r.squared = 1 - (ss_residual/df_residual) / (ss_total/df_total) 
  #   return(r.squared) }
  # 
  # # adjusted weighted R2 between predicted signal & observed signal on adjusted z scale
  # get.adjR2(z_obs = y,
  #           z_pred = y_pred,
  #           weights = weights_est,
  #           p = sum(coefs!=0),
  #           intercept = FALSE) 
  # 
  # # Ideally I think what we would like to have is an unbiased estimator of the
  # # weighted R2 between the predicted signal & the true underlying signal without noise 
  # # on the adjusted z scale (this one you could return if beta_true is given) :
  # get.R2(z_obs = y_true, 
  #        z_pred = y_pred,
  #        weights = weights_true,
  #        intercept = FALSE) 
  
  # TODO: also add analytically calculated leave one out CV error, cf
  # earlier code Chris
  # z.res <- (y.val - y.hat)/family$mu.eta(eta.hat) # the residuals on the adjusted scale = z.res = z - z.fit = eta + (y - mu)/g'(eta) - eta = (y - mu)/g'(eta)
  # w <- w.val * as.vector(family$mu.eta(eta.hat)^2 / family$variance(y.hat)) # IWLS weights from last iteration
  # H.diag <- compute.diagH(fit = fit, X = x.val, w = w)
  # loocv <- 1/length(y.val) * sum(w * (z.res/(1 - H.diag))^2) 
  # see Syed 2011, p.7 in https://scholarworks.gsu.edu/cgi/viewcontent.cgi?referer=https://scholar.google.com/&httpsredir=1&article=1100&context=math_theses
  
  
  return(list(p=p, n=n, m=m, 
              deviance=deviance, deviance_true=deviance_true, loglik=loglik, 
              MD_coefs=MD_coefs, MAD_coefs=MAD_coefs, sMAPE_coefs=sMAPE_coefs, 
              bias_coefs=bias_coefs, R2_coefs=R2_coefs,
              MD_selcoefs=MD_selcoefs, MAD_selcoefs=MAD_selcoefs, sMAPE_selcoefs=sMAPE_selcoefs, 
              bias_selcoefs=bias_selcoefs, R2_selcoefs=R2_selcoefs,
              support_size_true=support_size_true, 
              support_size_estimated=support_size_estimated, 
              estimated_over_true_support_size=estimated_over_true_support_size,
              TP=TP, FP=FP, TN=TN, FN=FN, 
              precision=precision, recall=recall, specificity=specificity, 
              accuracy=accuracy, F1_score=F1_score, MCC=MCC,
              gic=IC$gic, mebic=IC$mebic, ebic=IC$ebic,
              bic=IC$bic, aic=IC$aic, mbic=IC$mbic) )
}

# Matthews correlation coefficient (Pearson's phi) of binary row support.
# Cast before multiplication to avoid integer overflow on large designs.
# A constant truth or prediction has undefined correlation: return NA.
support_mcc <- function(TP, FP, TN, FN) {
  counts <- as.double(c(TP, FP, TN, FN))
  TP <- counts[1L]; FP <- counts[2L]; TN <- counts[3L]; FN <- counts[4L]
  margins <- c(TP + FP, TP + FN, TN + FP, TN + FN)
  if (any(!is.finite(margins)) || any(margins <= 0)) return(NA_real_)
  (TP * TN - FP * FN) / prod(sqrt(margins))
}


# return information criteria based on log likelihood, n, p & edf
# for a good discussion on information criteria see
# https://www.sciencedirect.com/science/article/pii/S016794731100171X
get.IC = function(loglik, n, p, edf, 
                  gamma_ebic = max(0, 1 - log(n)/(2*log(p))), # reduces to BIC when p <= sqrt(n)
                  Cn_hbic = log(log(n))) {
  min2LL = -2 * loglik
  
  # p=nr of vars, n=nr of cases, edf=effective degrees of freedom
  aic = min2LL + 2 * edf # Akaike information criterion
  # aicc = ifelse(n > edf, min2LL + 2*edf*n/(n-edf), NA) # small sample AIC
  bic = min2LL + log(n) * edf # Bayesian information criterion
  gic = min2LL + log(p) * log(log(n)) * edf # generalized information criterion GIC = SIC in https://www.pnas.org/doi/10.1073/pnas.2014241117
  mebic = min2LL + 2 * log(p) * log(log(n)) * edf # Zhang (2024), doi:10.1016/j.csda.2023.107860
   
  ebic = min2LL + (log(n) + 2 * gamma_ebic * log(p)) * edf # linearized EBIC; Chen and Chen (2008), Biometrika 95, 759-771
  hbic = min2LL + Cn_hbic * log(p) * edf # HBIC, recommended for ultra-high dimensional settings, cf https://arxiv.org/pdf/1311.4981 & https://www.annualreviews.org/docserver/fulltext/statistics/7/1/annurev-statistics-030718-105038.pdf
  mbic = min2LL + log(n * (p ^ 2) / 16) * edf # multiplicity-adjusted mBIC; Bogdan et al. (2008), doi:10.1002/qre.936; distinct from mBIC2
  # TODO maybe still add mbic2, cf Frommlet paper
  return(data.frame(aic = aic, bic = bic, ebic = ebic, hbic = hbic, gic = gic, mebic = mebic, mbic = mbic))
}


# Shared, case-insensitive names for fitting, scoring and tuning interfaces.
L0glm_ic_name <- function(IC, allow_cv = FALSE) {
  choices <- c("aic", "bic", "gic", "mebic", "ebic", "mbic")
  if (allow_cv) choices <- c(choices, "cv")
  if (!is.character(IC) || length(IC) != 1L || is.na(IC) ||
      !tolower(IC) %in% choices) {
    stop(if (allow_cv) "tune must be one of " else "Information criterion must be one of ",
         paste(shQuote(choices), collapse = ", "),
         " (case-insensitive).", call. = FALSE)
  }
  tolower(IC)
}

L0glm_resolve_lambda <- function(lambda, p, n) {
  if (is.character(lambda)) lambda <- preset.lambda(lambda, p = p, n = n)
  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda)) {
    stop("lambda must be a finite numeric scalar or an information-criterion name; ",
         "see help('information_criteria', package = 'L0glm').", call. = FALSE)
  }
  lambda
}


# preset.lambda ####
# Per-coefficient penalty on the -2*logLik scale; see information_criteria.
preset.lambda <- function(IC="gic", p, n,
                          gamma_ebic = max(0, 1 - log(n)/(2*log(p)))) { # reduces to BIC when p <= sqrt(n)
  IC <- L0glm_ic_name(IC)
  lambda <- switch(IC,
                   gic = log(p)*log(log(n)),
                   mebic = 2*log(p)*log(log(n)),
                   aic = 2,
                   bic = log(n),
                   ebic = log(n) + 2*gamma_ebic*log(p),
                   mbic = log(n*(p^2)/16))
  return(lambda)
}
