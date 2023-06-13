# function to calculate fit metrics to tune lambda ####
get.metrics <- function(beta, X, y, weights=matrix(1, nrow=nrow(y), ncol=ncol(y)), 
                        offset=Matrix(0, nrow=nrow(X), ncol=ncol(y), sparse=T), 
                        family=gaussian(identity), 
                        beta_true=NULL, 
                        model_name="L0glm",
                        edf=sum(beta!=0), # alt option: sum(rowMin(abs(beta))!=0)*ncol(beta); sum(rowMin(abs(beta))!=0) didn't work
                        lower=rep(0, ncol(X)),
                        upper=rep(Inf,ncol(X)) ) {
  
  # nr of variables, cases & tasks & effective degrees of freedom ####
  p = ncol(X) # nr of variables
  n = nrow(X) # nr of cases
  m = ncol(y) # nr of tasks
  # edf = sum(rowMin(abs(beta))!=0)*m # we use this as effective degrees of freedom, note: abess just used sum(beta!=0) instead (with nonnegativity constraints this is smaller)
  
  beta[beta<lower] = lower
  beta[beta>upper] = upper
  
  # metrics related to support recovery ####
  # true support & estimated support & all variable sets
  if (!is.null(beta_true)) vars_true = which(rowMin(abs(beta_true))!=0) else vars_true = NA
  vars_estimated = which(rowMin(abs(beta))!=0)
  
  # support size & ratio of estimated support size over true support size
  if (!is.null(beta_true)) support_size_true = length(vars_true) else support_size_true = NA
  support_size_estimated = length(vars_estimated) 
  if (!is.null(beta_true)) estimated_over_true_support_size = support_size_estimated / support_size_true else estimated_over_true_support_size = NA
  
  # Calculate TP, TN, FP, FN
  if (!is.null(beta_true)) {
    TP = length(intersect(vars_true, vars_estimated)) # Intersection of true and estimated sets
    FP = length(setdiff(vars_estimated, vars_true)) # In estimated but not in true
    TN = p-length(union(vars_true, vars_estimated)) # In neither
    FN = length(setdiff(vars_true, vars_estimated)) # In true but not in estimated
    
    # calculate precision, recall, sensitivity, accuracy & F1 score to recover true support from these
    precision = TP / (TP+FP) # = positive predictive value
    recall = TP / (TP+FN) # = sensitivity = TPR
    specificity = TN / (TN+FP) # = TNR
    accuracy = (TP+TN) / (TP+TN+FP+FN)
    F1_score = 2*precision*recall / (precision + recall) # = harmonic mean of precision and recall
    # metrics related to accuracy of coefficients ####
    MD_coefs = mean(beta-beta_true) # mean deviation
    MAD_coefs = mean(abs(beta-beta_true)) # mean absolute deviation
    sMAPE_coefs =  mean(2 * abs(beta - beta_true) / (abs(beta) + abs(beta_true) + 1E-10)) # symmetric mean absolute percentage error, 0 = perfect agreement
    bias_coefs = as.vector(coef(lm(as.vector(beta_true) ~ 0 + as.vector(beta)))[[1]]) # regression beta true on beta (closer to 1 is better, >1 means you underestimate true coefs, <1 means you overestimate true coefs)
    R2_coefs = summary(lm(as.vector(beta_true) ~ 0+as.vector(beta)))$r.squared # I use R2 from zero intercept regression of beta_true on beta
  } else {
    TP <- FP <- FN <- TN <- precision <- recall <- sensitivity <- accuracy <- specificity <- F1_score <- bias_coefs <- MD_coefs <- MAD_coefs <- sMAPE_coefs <- R2_coefs <- NA
  }
  
  # metrics related to prediction accuracy ####
  # calculate deviance
  if (is.null(beta_true)) { # calculated weights based on expected value
    y = as.matrix(y)
    # Compute the fitted mean values using the fitted coefficients
    nonzero = as.vector(rowMax(abs(beta))!=0) # nonzero coefficients
    eta = X[,nonzero,drop=F] %*% beta[nonzero,,drop=F] + offset # fitted eta = fitted values on link scale
    mu = family$linkinv(eta) # fitted values on response scale
    # Compute the weights at the fitted values
    gprime = family$mu.eta(eta) # derivative of the expected values
    varg = family$variance(mu) # variance
    W = as.vector(weights)*(gprime^2)/varg # weights
    good = !is.na(W) & W != 0 & !is.infinite(W) & gprime!=0 & as.vector(weights)!=0 # valid observations
    W[!good] = 0
    dim(W) = dim(eta)
    
    # calculated deviance or true theoretical deviance compared to ground truth
    dev_res = family$dev.resids(as.vector(y), as.vector(mu), as.vector(W))
    dev_res[is.infinite(dev_res)|is.na(dev_res)] = 0
    deviance = sum(dev_res) # deviance 
    
  } else { # calculate weights based on the "true" responses if beta_true is given
    # calculate true response used in simulation given true coefficients beta_true
    nonzero = as.vector(rowMax(abs(beta_true))!=0) # nonzero coefficients
    eta_true = X[,nonzero,drop=F] %*% beta_true[nonzero,,drop=F] + offset # true eta = true values on link scale
    mu_true <- as.matrix(family$linkinv(eta_true)) # true mu = true values on response scale
    
    # Compute the weights based on the "true" responses
    gprime_true = family$mu.eta(eta_true) # derivative of the true values
    varg_true = family$variance(mu_true)
    W = as.vector(weights)*(gprime_true^2)/varg_true
    good = !is.na(W) & W != 0 & !is.infinite(W) & gprime_true!=0 & as.vector(weights)!=0
    W[!good] = 0
    dim(W) = dim(eta_true)
    
    # Compute the fitted mean values using the actual fitted coefficients
    nonzero = as.vector(rowMax(abs(beta))!=0) # nonzero coefficients
    eta = X[,nonzero,drop=F] %*% beta[nonzero,,drop=F] + offset # fitted eta = fitted values on link scale
    mu = family$linkinv(eta) # fitted values on response scale
    
    # calculated true theoretical deviance compared to ground truth
    dev_res = family$dev.resids(as.vector(mu_true), as.vector(mu), as.vector(W))
    dev_res[is.infinite(dev_res)|is.na(dev_res)] = 0
    deviance = sum(dev_res) # deviance 
  }
  
  # calculate some information criteria See Zhang & Chen 2010 ####
  # see also https://github.com/abess-team/abess/blob/master/R-package/R/deviance.abess.R
  # Get -2 * log-likelihood & log likelihood
  y = as.vector(y)
  # y = as.vector(mu_true)
  min2LL = family$aic(y = as.vector(y), mu = as.vector(mu), wt = as.vector(W), # TODO for family poisson gives Error in dpois(y, mu, log = TRUE) : Non-numeric argument to mathematical function
                      n = NULL, # TODO check what n should be for family binomial(logit)
                      dev = as.vector(deviance))
  loglik = -1/2 * min2LL
  
  # p=nr of vars, n=nr of cases, k=effective degrees of freedom
  gic = min2LL + log(p)*log(log(n))*edf # generalized information criterion GIC = SIC in https://www.pnas.org/doi/10.1073/pnas.2014241117 
  aic = min2LL + 2*edf # Akaike information criterion
  # aicc = ifelse(n > edf, min2LL + 2*edf*n/(n-edf), NA) # small sample AIC
  bic = min2LL + log(n)*edf # Bayesian information criterion
  ebic = min2LL + (log(n) + 2*(1-log(n)/(2*log(p)))*log(p))*edf # extended BIC, Chen, J. and Chen, Z. (2008). Extended Bayesian information criterion for model selection with large model space. Biometrika, 94, 759-771., https://arxiv.org/abs/1107.2502
  mbic = min2LL + log(n*(p^2)/16)*edf # see Frommlet & Nuel 2016 & Bogdan M, Ghosh JK, Zak-Szatkowska M (2008) Selecting explanatory variables with the modified version of bayesian information criterion. Quality and Reliability Engineering International 24: 627–641. doi: 10.1002/qre.936
  # c = 2 # TODO default, add as a parameter
  # hq = min2LL + c*log(log(n))*edf # Hannan and Quinnn information criterion - maybe worth to keep?
  # ric = min2LL + 2 * log(p) # risk inflation criterion
  # mric = min2LL + 2 * sum(log(p/(1:edf))) # modified risk inflation criterion
  # cic = min2LL + 4 * sum(log(p/(1:edf))) # covariance inflation criterion
  # See Xu thesis 2010, Model Selection with Information Criteria
  # g = 1 # TODO Default, add as argument
  # bicg =  min2LL + log(n)*edf + 2*g*lchoose(p,round(edf))
  # q = 0.25 # TODO Default, add as argument
  # bicq =  min2LL + log(n)*edf - 2*edf*log(q/(1-q))
  
  # weighted_rmse = sqrt(sum((predicted-actual)^2*weight)/sum(weight))
  # z = (eta - offset) + (y - mu) / mu.prime # working (adjusted) response
  # z.res <- (y - mu)/fam$mu.eta(eta) # the residuals on the adjusted scale = z.res = z - z.fit = eta + (y - mu)/g'(eta) - eta = (y - mu)/g'(eta)
  # still add explained variance or adjusted explained variance on adjusted z scale?
  
  return(list(model_name=model_name, p=p, n=n, m=m, 
              deviance=deviance, loglik=loglik,
              MD_coefs=MD_coefs, MAD_coefs=MAD_coefs, sMAPE_coefs=sMAPE_coefs, 
              bias_coefs=bias_coefs, 
              R2_coefs=R2_coefs,
              support_size_true=support_size_true, 
              support_size_estimated=support_size_estimated, 
              estimated_over_true_support_size=estimated_over_true_support_size,
              TP=TP, FP=FP, TN=TN, FN=FN, 
              precision=precision, recall=recall, specificity=specificity, 
              accuracy=accuracy, F1_score=F1_score,
              gic=gic, ebic=ebic, bic=bic, aic=aic, mbic=mbic) )
  # hq=hq, aicc=aicc, ric=ric, mric=mric, cic=cic, bicg=bicg, bicq=bicq))
}
