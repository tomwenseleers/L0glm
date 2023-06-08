# library(remotes)
# remotes::install_github("tomwenseleers/L0glm/L0glm")
# library(L0glm)
devtools::load_all()

# calculate fit metrics to tune lambda
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
    TP <- FP <- FN <- TN <- precision <- recall <- sensitivity <- accuracy <- F1_score <- bias_coefs <- MAD_coefs <- sMAPE_coefs <- R2_coefs <- NA
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
    good = !is.na(W) & W != 0 & !is.infinite(W) & gprime_true!=0 & as.vector(weights)!=0 # valid observations
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
  gic = min2LL + log(p)*log(log(n))*edf # generalized information criterion GIC 
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
  
  return( list(model_name=model_name, p=p, n=n, m=m, 
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



# EXAMPLE MULTICHANNEL DECONVOLUTION ####
# simulate blurred multichannel spike train ####
set.seed(1)
s <- 0.1 # sparsity (% of timepoints where there is a peak)
p <- 200 # 500 variables
# simulate multichannel blurred spike train with Gaussian noise
sd_noise <- 1
sim <- simulate_spike_train(n=p, p=p, k=round(s*p), 
                            mean_beta = 10000,
                            sd_logbeta = 1,
                            family="gaussian",
                            sd_noise = sd_noise,
                            multichannel=TRUE, 
                            sparse=TRUE)
X <- sim$X # n x p matrix
Y <- sim$y # n x m matrix  
colnames(X) = paste0("x", 1:ncol(X)) # needed for abess
colnames(Y) = paste0("y", 1:ncol(Y))
true_coefs <- sim$beta_true
m <- ncol(Y) # nr of tasks
n <- nrow(X) # nr of observations
p <- ncol(X) # nr of independent variables (shifted copies of point spread functions)
W <- 1/(Y+0.1) # approx 1/variance Poisson observation weights with family="poisson", n x m matrix


# relaxed group L0 performance metrics over lambda sequence  ####

# lambda values that in principle would be expected to maximize aic, bic, gic, ebic or mbic 
# can be calculated using
# (NOTE: with high dispersion/variance these come out way too high though, except if sd_noise=1 - check this)
library(glmnet)
# if dispersion/variance is known
lambda_preset_aic_knowndisp = preset.lambda(IC="aic", dispersion=sd_noise^2)  
lambda_preset_aic_knowndisp # 0.5
lambda_preset_bic_knowndisp = preset.lambda(IC="bic", dispersion=sd_noise^2)
lambda_preset_bic_knowndisp # 1.32
lambda_preset_gic_knowndisp = preset.lambda(IC="gic", dispersion=sd_noise^2)  
lambda_preset_gic_knowndisp # 2.21
lambda_preset_ebic_knowndisp = preset.lambda(IC="ebic", dispersion=sd_noise^2) 
lambda_preset_ebic_knowndisp # 2.65
lambda_preset_ebic_knowndisp = preset.lambda(IC="mbic", dispersion=sd_noise^2) # 3.28
lambda_preset_ebic_knowndisp # 3.28
# if dispersion/variance is not known it can first be estimated
lambda_preset_aic = preset.lambda(IC="aic", y=sim$y[,which.max(colMeans(abs(sim$y)))], 
              X=sim$X, weights=rep(1, nrow(sim$y)), 
              family=gaussian(identity))  
lambda_preset_aic # 0.7955769
lambda_preset_bic = preset.lambda(IC="bic", y=sim$y[,which.max(colMeans(abs(sim$y)))], 
              X=sim$X, weights=rep(1, nrow(sim$y)), 
              family=gaussian(identity))
lambda_preset_bic # 2.11
lambda_preset_gic = preset.lambda(IC="gic", y=sim$y[,which.max(colMeans(abs(sim$y)))], 
              X=sim$X, weights=rep(1, nrow(sim$y)), 
              family=gaussian(identity))
lambda_preset_gic # 3.51
lambda_preset_ebic = preset.lambda(IC="ebic", y=sim$y[,which.max(colMeans(abs(sim$y)))], 
              X=sim$X, weights=rep(1, nrow(sim$y)), 
              family=gaussian(identity))
lambda_preset_ebic # 4.21
lambda_preset_mbic = preset.lambda(IC="mbic", y=sim$y[,which.max(colMeans(abs(sim$y)))], 
              X=sim$X, weights=rep(1, nrow(sim$y)), 
              family=gaussian(identity))  # 5.219924
lambda_preset_mbic # 5.22


# below we fit over sequence of lambda values to see if these IC are indeed optimal at these values
# NOTE: we could also use golden section search
# using https://search.r-project.org/CRAN/refmans/cmna/html/goldsect.html
# as in abess

# for the lambda sequence we will below use lambda_preset_aic and the theoretical
# lambda that would penalize all variables to zero with an L0 penalty as maximum
# for this max lambda
# see https://stats.stackexchange.com/questions/416144/minimum-and-maximum-regularization-in-l0-pseudonorm-penalized-regression
# max lambda = max(diag(1/crossprod(X, X)) * (crossprod(X, y)^2))
# = max((1/colSums(X^2)) * (crossprod(X, y)^2))
maxlambda = max((1/colSums(sim$X^2)) * (crossprod(sim$X, sim$y)^2)) # lambda that would penalize all variables to 0 in L0 regression (but maybe factor of 4 of Frommlet still needs to be taken into account - check with orthogonal variable problem)

lambdas = 10^seq((round(log10(lambda_preset_aic))-3), round(log10(maxlambda))+1, by=1)

library(parallel)
library(qlcMatrix)
system.time(warm_start <- L0glm.fit_new(X = sim$X, y = sim$y,
                               family = gaussian(identity), # poisson(identity),
                               lambda = min(lambdas), #  # or 0.1 or 1
                               lower = rep(0, ncol(sim$X)),
                               upper = rep(Inf, ncol(sim$X)),
                               nonnegative = TRUE,
                               thresh = 1,
                               delta = 1E-2,
                               group = T, # use group L0 penalty
                               group.maxit = 5,
                               relaxedGroupL0 = T, # return unpenalized fit for selected variables at the end, similar to concept of relaxed LASSO
                               warmupiters = 5,
                               solver = "eigen",
                               method = 7L # using Eigen's sparse Cholesky solver
                               )) # 5s
sum(rowMin(abs(warm_start))!=0) # 21 variables selected

no_warm_start = matrix(1, nrow=nrow(warm_start), ncol=ncol(warm_start))
w = 0.95 # 0 = no warm starts, 1 = full warm starts
system.time(coefficients <- lapply(lambdas, function (lambda) {
                                new_coefs <- L0glm.fit_new(X = sim$X, y = sim$y, 
                                                                            family = gaussian(identity), # poisson(identity),
                                                                            start = w*warm_start+(1-w)*no_warm_start, # use warm starts
                                                                            lambda = lambda, # needs to be larger with larger sd_noise for family=gaussian, e.g. 0.5 with sd_noise=1 & 5 with sd_noise=5 
                                                                            lower = rep(0, ncol(sim$X)),
                                                                            upper = rep(Inf, ncol(sim$X)),
                                                                            nonnegative = TRUE,
                                                                            thresh = 1, # coefficients below this threshold are set to zero
                                                                            delta = 0, # TODO: if we have the thresh argument we can get rid of delta & set delta to zero, no need to have both
                                                                            group = T, # use group L0 penalty
                                                                            group.maxit = 6,
                                                                            relaxedGroupL0 = T, # return unpenalized fit for selected variables at the end, similar to concept of relaxed LASSO
                                                                            warmupiters = 5,
                                                                            solver = "eigen",
                                                                            method = 7L, # using Eigen's sparse Cholesky solver
                                                                            parallellize = T 
                                                                            )
                                warm_start <<- new_coefs
                                return(new_coefs)
}
)) # 66s, TODO check how warm starts would affect solution quality & timings

# metrics relaxed group L0 over sequence of lambda values
library(qlcMatrix)
metrics <- data.frame(do.call(rbind, lapply(coefficients, function (coefs) get.metrics(beta=coefs, X=sim$X, y=sim$y, 
                                                                                       weights=matrix(1, nrow=nrow(sim$y), ncol=ncol(sim$y)), 
                                                                                       offset=Matrix(0, nrow=nrow(sim$X), ncol=ncol(sim$y), sparse=T), 
                                                                                       family=gaussian(identity), 
                                                                                       beta_true=sim$beta_true,
                                                                                       edf=sum(sim$beta!=0) # sum(rowMin(abs(sim$beta))!=0)*ncol(sim$beta) # sum(rowMin(abs(sim$beta))!=0)*ncol(sim$beta) # sum(rowMin(abs(sim$beta))!=0) # sum(sim$beta!=0) # sum(rowMin(abs(sim$beta))!=0) # sum(rowMin(abs(sim$beta))!=0)*ncol(sim$beta) # sum(sim$beta!=0) # alt option: sum(rowMin(abs(beta))!=0)*ncol(beta); sum(rowMin(abs(beta))!=0) didn't work
) 
)))

# best empirical lambda based on calculated GIC over lambda sequence:
metric = "gic"
bestlambda = lambdas[which.min(metrics[,metric])]
bestlambda # 100

bestlambda = lambdas[which.min(abs(lambdas-bestlambda))]
beta_relaxedgroupL0 = as.matrix(coefficients[[which(lambdas==bestlambda)]])
# beta_relaxedgroupL0[beta_relaxedgroupL0<10] = 0

image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=as.matrix(sim$beta_true)^0.1, 
      col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="True coefficients")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim$X))[as.vector(rowMin(beta_relaxedgroupL0)!=0)], col="cyan")

image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=beta_relaxedgroupL0^0.1, 
      col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="nonnegative relaxed group L0")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim$X))[as.vector(rowMin(beta_relaxedgroupL0)!=0)], col="cyan")
sum(rowMax(sim$beta_true)!=0) # 20 true peaks
sum(rowMin(beta_relaxedgroupL0)!=0) # 20 peaks detected
length(intersect(which(rowMax(sim$beta_true)!=0), 
                 which(rowMin(beta_relaxedgroupL0)!=0))) # 20 out of 20 peaks correctly detected
length(which(sim$beta_true!=0)) # 1606 coefs nonzero
length(as.vector(sim$beta_true)) # 29800 coefs in total
length(intersect(which(sim$beta_true!=0), 
                 which(beta_relaxedgroupL0!=0))) # 1424 true pos coefs
length(intersect(which(sim$beta_true==0), 
                 which(beta_relaxedgroupL0!=0))) # 126 false pos coefs
plot(as.vector(sim$beta_true), as.vector(beta_relaxedgroupL0), pch=16)
cor(as.vector(sim$beta_true), as.vector(beta_relaxedgroupL0)) # 0.9999866


plot(x=log10(lambdas), y=metrics$support_size_estimated, pch=16, type="o")
abline(v=log10(bestlambda), col="red")
abline(h=metrics$support_size_true, col="blue")
plot(x=log10(lambdas), y=metrics$deviance, pch=16, type="o", 
     ylim=c(0,max(unlist(metrics$deviance)))) # theoretical deviance relative to true underlying signal without noise
lambda_bestdeviance = lambdas[which.min(metrics$deviance)]
lambda_bestdeviance
abline(v=log10(bestlambda), col="red")
abline(v=log10(lambda_bestdeviance), col="blue")
plot(x=log10(lambdas), y=metrics$loglik, pch=16, type="o") # theoretical real log likelihood
lambda_bestloglik = lambdas[which.max(metrics$loglik)]
lambda_bestloglik
abline(v=log10(bestlambda), col="red")
abline(v=log10(lambda_bestloglik), col="blue")
plot(x=log10(lambdas), y=metrics$MD_coefs, pch=16, type="o")
abline(v=log10(bestlambda), col="red")
plot(x=log10(lambdas), y=metrics$MAD_coefs, pch=16, type="o") # optimal at lambda=0.1
lambdas[which.min(metrics$MAD_coefs)] # 0.1
abline(v=log10(bestlambda), col="red")
plot(x=log10(lambdas), y=metrics$sMAPE_coefs, pch=16, type="o") # optimal at lambda=2.828427 - ok metric for bias on coefs?
lambdas[which.min(metrics$sMAPE_coefs)] # 2.828427 - seems too high?
abline(v=log10(bestlambda), col="red")
plot(x=log10(lambdas), y=metrics$bias_coefs, pch=16, type="o", ylim=c(0.9,1.1)) # closest to 1 best
lambda_minimal_bias = min(lambdas[which(abs(unlist(metrics$bias_coefs)-1)<=min(abs(unlist(metrics$bias_coefs)-1), na.rm=T)*1.3)]) # 1
lambda_minimal_bias # 0.71
abline(v=log10(bestlambda), col="red")
abline(h=1, col="blue")
plot(x=log10(lambdas), y=metrics$R2_coefs, pch=16, type="o", ylim=c(0, 1)) # highest lambda before sudden drop in R2 optimal
abline(v=log10(bestlambda), col="red") # best would be lambda of 1
plot(x=log10(lambdas), y=metrics$estimated_over_true_support_size, pch=16, type="o") # lambda 1 close to optimal
lambda_correct_support_size = min(lambdas[which.min(abs(unlist(metrics$estimated_over_true_support_size)-1))]) #
lambda_correct_support_size # 2
abline(v=log10(bestlambda), col="red") # lambda of 1 close to optimum
abline(v=log10(lambda_correct_support_size), col="blue")
plot(x=log10(lambdas), y=metrics$precision, pch=16, type="o")
abline(v=log10(bestlambda), col="red")
plot(x=log10(lambdas), y=metrics$recall, pch=16, type="o")
abline(v=log10(bestlambda), col="red")
plot(x=log10(lambdas), y=metrics$specificity, pch=16, type="o")
abline(v=log10(bestlambda), col="red")
plot(x=log10(lambdas), y=metrics$accuracy, pch=16, type="o") # GOOD FOR THEORETICAL BEST SUPPORT RECOVERY, has max at lambda=1
lambda_best_accuracy = lambdas[which.max(metrics$accuracy)] # 
lambda_best_accuracy
abline(v=log10(bestlambda), col="red")
abline(v=log10(lambda_best_accuracy), col="blue")
plot(x=log10(lambdas), y=metrics$F1_score, pch=16, type="o") # BEST FOR THEORETICAL BEST SUPPORT RECOVERY, has max at lambda=1
lambda_best_F1score = lambdas[which.max(metrics$F1_score)] # 
lambda_best_F1score # 2.8
abline(v=log10(bestlambda), col="red")
abline(v=log10(lambdas[which.max(metrics$F1_score)]), col="blue")
plot(x=log10(lambdas), y=metrics$ebic, pch=16, type="o") # keep, cf https://www.researchgate.net/publication/46468889_Consistency_of_the_Empirical_Bayes_Information_Criterion_for_Selecting_Variables_in_Linear_Mixed_Models#:~:text=Compared%20to%20mAIC%2C%20cAIC%20and%20BIC%2C%20the%20EBIC,of%20covariates%20and%20the%20number%20of%20clusters%20increase.
lambda_best_EBIC = lambdas[which.min(metrics$ebic)] #
lambda_best_EBIC # 1
abline(v=log10(bestlambda), col="red")
abline(v=log10(lambda_best_EBIC), col="blue") # corresponds to lambda giving optimal support recovery
plot(x=log10(lambdas), y=metrics$gic, pch=16, type="o") # keep
lambda_best_gic = lambdas[which.min(metrics$gic)] # 1
abline(v=log10(bestlambda), col="red")
abline(v=log10(lambda_best_gic), col="blue") # corresponds to lambda giving optimal support recovery
plot(x=log10(lambdas), y=metrics$aic, pch=16, type="o") # keep - no good though
lambda_best_aic = lambdas[which.min(metrics$aic)] # 0.088 - no good
abline(v=log10(bestlambda), col="red")
abline(v=log10(lambda_best_aic), col="blue") # 
plot(x=log10(lambdas), y=metrics$bic, pch=16, type="o") # keep
lambda_best_bic = lambdas[which.min(metrics$bic)] # 1
abline(v=log10(bestlambda), col="red")
abline(v=log10(lambda_best_aic), col="blue") # keep
plot(x=log10(lambdas), y=metrics$mbic, pch=16, type="o") # keep
lambda_best_mbic = lambdas[which.min(metrics$mbic)] # 1
abline(v=log10(bestlambda), col="red")
abline(v=log10(lambda_best_mbic), col="blue")





# abess multitask mgaussian fit over lambda sequence
library(abess)
allvars = 1:ncol(X)
preselvars = which(rowMin(beta_mgaussian_glmnet)!=0) # pre-select variables based on nonnegative mgaussian group LASSO
system.time(abess_fit <- abess(x = X[,preselvars,drop=F], 
                               y = Y, 
                               # weights = 1/(Y+0.1), # works in abess
                               family = "mgaussian",
                               # normalize = 3, # or 1?
                               tune.path = "sequence", # gsection", # or sequence
                               support.size = c(1:150),
                               # gs.range = c(1,200),
                               lambda=0, 
                               # init.active.set=,  # vector of integers
                               # splicing.type = 2,
                               max.splicing.iter = 20,
                               warm.start = TRUE,
                               tune.type = "ebic") # or cv or bic or aic or gic
) # 1.78s, 27.68s with p=10000
# abess_fit
extract(abess_fit)$support.size # 19
sum(rowMin(abs(sim$beta_true))!=0) # true support 20
plot(abess_fit)
plot(abess_fit, type="tune")

# metrics multitask mgaussian abess fit
support.size <- abess_fit$support.size
abess_fit$edf
coefficients <- abess_fit$beta # lapply(abess_fit$beta, function (coefs) { out = Matrix(0, nrow=p, ncol=m, sparse=T) } )
library(qlcMatrix)
metrics <- data.frame(do.call(rbind, lapply(1:length(coefficients), function (i) { coefs=coefficients[[i]]
                                                                             # edf=sum(rowMin(abs(coefs))!=0) # *ncol(coefs)
                                                                             coefs[coefs<0] = 0;
                                                                             get.metrics(beta=coefs, X=sim$X, y=sim$y, 
                                                                                       weights=matrix(1, nrow=nrow(sim$y), ncol=ncol(sim$y)), 
                                                                                       offset=Matrix(0, nrow=nrow(sim$X), ncol=ncol(sim$y), sparse=T), 
                                                                                       family=gaussian(identity), 
                                                                                       beta_true=sim$beta_true,
                                                                                       edf=sum(rowMin(abs(coefs))!=0)*ncol(coefs), # sum(coefs!=0), # abess_fit$edf[[i]], # sum(coefs!=0),
                                                                                       lower=rep(0, ncol(sim$X)),
                                                                                       upper=rep(Inf,ncol(sim$X)) ) } )))
sel=1:100

plot(x=support.size[sel], y=metrics$deviance[sel], pch=16, type="o", ylim=c(0,max(unlist(metrics$deviance[sel])))) # highest lambda before sudden jump in deviance corresponds with optimal lambda
abline(v=metrics$support_size_true, col="red")
plot(x=support.size[sel], y=metrics$loglik[sel], pch=16, type="o")
abline(v=metrics$support_size_true, col="red")
plot(x=support.size[sel], y=metrics$MD_coefs[sel], pch=16, type="o")
abline(v=metrics$support_size_true, col="red")
plot(x=support.size[sel], y=metrics$MAD_coefs[sel], pch=16, type="o") # optimal at support=11
support.size[which.min(metrics$MAD_coefs)] # 11
abline(v=metrics$support_size_true, col="red")
plot(x=support.size[sel], y=metrics$sMAPE_coefs[sel], pch=16, type="o") # optimal at support=11
support.size[which.min(metrics$sMAPE_coefs)] # 11
abline(v=metrics$support_size_true, col="red")
plot(x=support.size[sel], y=metrics$bias_coefs[sel], pch=16, type="o", ylim=c(0.1,2)) #
support.size[which.min(metrics$bias_coefs)] # 499
abline(v=metrics$support_size_true, col="red")
abline(h=1, col="blue")
plot(x=support.size[sel], y=metrics$R2_coefs[sel], pch=16, type="o", ylim=c(0, 1)) #
abline(v=metrics$support_size_true, col="red")
plot(x=support.size[sel], y=metrics$precision[sel], pch=16, type="o")
abline(v=metrics$support_size_true, col="red")
plot(x=support.size[sel], y=metrics$recall[sel], pch=16, type="o")
abline(v=metrics$support_size_true, col="red")
plot(x=support.size[sel], y=metrics$specificity[sel], pch=16, type="o")
abline(v=metrics$support_size_true, col="red")
plot(x=support.size[sel], y=metrics$accuracy[sel], pch=16, type="o") # BEST FOR SUPPORT RECOVERY, has max at support=33
support.size[which.max(metrics$accuracy)] # 33
abline(v=metrics$support_size_true, col="red")
abline(v=support.size[sel][which.max(as.vector(metrics$accuracy[sel]))], col="blue")
plot(x=support.size[sel], y=metrics$F1_score[sel], pch=16, type="o") # BEST FOR SUPPORT RECOVERY, has max at support=38
support.size[which.max(metrics$F1_score)] # 38
abline(v=metrics$support_size_true, col="red")
abline(v=support.size[sel][which.max(as.vector(metrics$F1_score[sel]))], col="blue")
plot(x=support.size[sel], y=metrics$gic[sel], pch=16, type="o") # keep
support.size[which.min(metrics$gic)] # 11
abline(v=metrics$support_size_true, col="red")
plot(x=support.size[sel], y=metrics$aic[sel], pch=16, type="o") # keep
support.size[which.min(metrics$aic)] # 11
abline(v=metrics$support_size_true, col="red")
plot(x=support.size[sel], y=metrics$aicc[sel], pch=16, type="o") # doesn't work
abline(v=metrics$support_size_true, col="red")
plot(x=support.size[sel], y=metrics$bic[sel], pch=16, type="o") # keep
support.size[which.min(metrics$bic)] # 11
abline(v=metrics$support_size_true, col="red")
plot(x=support.size[sel], y=metrics$ebic[sel], pch=16, type="o") # keep
support.size[which.min(metrics$ebic)] # 11
abline(v=metrics$support_size_true, col="red")
plot(x=support.size[sel], y=metrics$hq[sel], pch=16, type="o") # keep?
support.size[which.min(metrics$hq)] # 11
abline(v=metrics$support_size_true, col="red")
plot(x=support.size[sel], y=metrics$ric[sel], pch=16, type="o") # doesn't work
abline(v=metrics$support_size_true, col="red")
plot(x=support.size[sel], y=metrics$mric[sel], pch=16, type="o") # doesn't work
abline(v=metrics$support_size_true, col="red")
plot(x=support.size[sel], y=metrics$cic[sel], pch=16, type="o") # doesn't work
abline(v=metrics$support_size_true, col="red")
plot(x=support.size[sel], y=metrics$bicg[sel], pch=16, type="o") # doesn't work
abline(v=metrics$support_size_true, col="red")
plot(x=support.size[sel], y=metrics$bicq[sel], pch=16, type="o") # 
support.size[which.min(metrics$bicq)] # 11
abline(v=metrics$support_size_true, col="red")





