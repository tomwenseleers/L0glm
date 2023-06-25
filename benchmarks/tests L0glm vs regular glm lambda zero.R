# Check of metrics (AIC, BIC, loglik) L0glm package & comparison regular GLM ####

library(remotes)
# remotes::install_github("tomwenseleers/L0glm/L0glm")
library(L0glm)
# devtools::load_all(".")
library(export)
source("..//benchmarks/generate data.R")

# check consistency of L0glm output with glm for gaussian identity link, poisson log link, Gamma log link & binomial logit link
n = 10000
p = 100
n_nonzero = 50
dispersion = 1
crit = "gic"

# 1. GAUSSIAN GLM ####

set.seed(1)
data_gaus  = simulate_glm(n, p, n_nonzero, family = "gaussian") 
X = data_gaus$X
y = data_gaus$y
family = gaussian(identity)

# glm fit
system.time(glm_gaus <- glm(data_gaus$y ~ 0+data_gaus$X, family=gaussian(identity))) # 0.16s
glm_gaus$deviance 
glm_gaus$aic
AIC(glm_gaus)
BIC(glm_gaus)
logLik(glm_gaus)

# L0 glm fit with lambda close to zero
system.time(L0glm_gaus_lambda0 <- L0glm.fit(X = Matrix(data_gaus$X, sparse=TRUE), 
                                            y = data_gaus$y, 
                                            family = gaussian(identity), 
                                            lambda = 1E-30, 
                                            nonzero.centered = TRUE, # use nonzero centered ridge regression?
                                            nonnegative = FALSE,
                                            warmupiters = 1,
                                            # start = rep(1, ncol(data_gaus$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                                            solver = "eigen", # 
                                            method = 3L, #
                                            thresh = 1E-10
))  # 0.08s
# TODO: FIX BUG: passing dense X = data_gaus$X gives error Error: No such slot: i. 
# TODO: now returns error Error in X %*% start :
# Cholmod error 'X and/or Y have wrong dimensions' at file ../MatrixOps/cholmod_sdmult.c, line 88

max(abs(L0glm_gaus_lambda0$coefficients-coef(glm_gaus))) # 1.138658e-09, correct
# check if metrics (loglik & AIC & BIC) are correct
metrics_gaus = data.frame(do.call(cbind, get.metrics(beta = L0glm_gaus_lambda0$coefficients,  # some metrics
                                                     X = Matrix(data_gaus$X, sparse=TRUE), 
                                                     y = data_gaus$y, 
                                                     family = gaussian(identity), 
                                                     beta_true = as.matrix(data_gaus$beta_true))))
logLik(glm_gaus)-metrics_gaus$loglik # 0, correct
AIC(glm_gaus)-metrics_gaus$aic # 0, correct
BIC(glm_gaus)-metrics_gaus$bic # 0, correct


# L0 glm fit with lambda preset to optimise given IC
lambda_preset_gaus = preset.lambda(IC = crit, # options aic, bic, gic, ebic, mbic 
                              y = data_gaus$y,
                              X = data_gaus$X, 
                              family = gaussian(identity),
                              dispersion = dispersion # dispersion here was 1, here we pass this value to not have to estimate it 
                              ) 
lambda_preset_gaus # 2.556246

system.time(L0glm_gaus_presetlambda <- L0glm.fit(X = Matrix(data_gaus$X, sparse=TRUE), 
                                            y = data_gaus$y, 
                                            family = gaussian(identity), 
                                            lambda = lambda_preset_gaus,
                                            nonzero.centered = TRUE, # use nonzero centered ridge regression?
                                            nonnegative = FALSE,
                                            warmupiters = 1,
                                            # start = rep(1, ncol(data_gaus$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                                            solver = "eigen", # 
                                            method = 3L, #
                                            thresh = 1E-10
)) # 0.74s
sum(L0glm_gaus_presetlambda$coefficients!=0) # 47


# coefficients & metrics over a range of lambda values
# calculate max lambda where all coefs would be penalized to 0
weights <- L0glm_gaus_presetlambda$weights
Xw <- Diagonal(x=sqrt(weights)) %*% data_gaus$X
sel_task <- which.max(colMeans(abs(as.matrix(data_gaus$y)))) # for multitask select column with max mean abs value?
yw <- as.matrix(data_gaus$y)[, sel_task,drop=F] * sqrt(weights) 
# TODO note that yw should really be zw, i.e. for non-identity link you should 
# be working on the adjusted response scale - maybe you should also make 
# L0glm.fit return z? Or make a convenience function that returns z and w 
# given a particular given family and coefficient estimates?
maxlambda <- max((1/colSums(Xw^2)) * (crossprod(Xw, yw)^2), na.rm=T)
if (family$family == "gaussian") maxlambda <- maxlambda/4
# could be that this factor of 1/4 applies to all families & that this is that 1/4 factor mentioned in Frommlet & Nuel
mincoef <- min(abs(L0glm_gaus_presetlambda$coefficients))

lambdas <- 10^seq(round(log10(lambda_preset_gaus))-3, log10(maxlambda)+1, by=0.2)
coefmat <- Matrix(0, nrow=length(lambdas), ncol=ncol(data_gaus$X))
warmstart <- 0.1 # % by which L0glm_gaus_presetlambda$coefficients will be used as start, lowest value that led to stable coefficient paths
system.time(metrics <- data.frame(do.call(rbind, lapply(1:length(lambdas), function (i) { 
  lambda <- lambdas[[i]]
  start <- (L0glm_gaus_presetlambda$coefficients)*warmstart+rep(mincoef/2,ncol(data_gaus$X))*(1-warmstart)
  # start[L0glm_gaus_presetlambda$coefficients==0] <- 0
  L0glm.out <- L0glm.fit(Matrix(data_gaus$X, sparse=TRUE), 
                         as.matrix(data_gaus$y), 
                         family = gaussian(identity), 
                         lambda = lambda, 
                         nonzero.centered = TRUE, # use nonzero centered ridge regression? - TODO: change that if you specify starting values start that it does not keep on shrinking coefs to those values but instead shrink to coef values from previous iteration
                         nonnegative = FALSE,
                         warmupiters = 1,
                         # start = start,
                         solver = "eigen", # eigen, glmnet or osqp, # note: constrained solvers glmnet & osqp better than eigen for very small lambda values
                         method = 3L, 
                         # gamma = 2, 
                         thresh = 1E-10,
                         relaxed = FALSE # TODO check : setting this argument to TRUE doesn't affect solution, while it should (it should do an L0glm fit & then refit an unpenalized LS model for the selected coefficients)
  )$coefficients 
  coefmat[i,] <<- as.vector(L0glm.out)
  return( unlist(c(lambda=lambda, get.metrics(beta = L0glm.out, 
                                              X = X, y = y, 
                                              family = gaussian(identity), 
                                              beta_true = data_gaus$beta_true))) )
} )))) # 20s
matplot(x=log10(metrics$lambda), 
        y=as.matrix(coefmat), type="l", lty=1,
        xlab="Log10(lambda)", ylab="Coefficient estimates", main="Regularisation path L0 Gaussian GLM")
abline(v=log10(maxlambda), col="red") # max lambda where all coefs would be penalized to 0

plot(log10(metrics$lambda), metrics$deviance, type="l")
metrics$lambda[which.min(metrics$deviance)] # 0.001 # NOTE this is deviance relative to observed signal, not deviance relative to ground truth though
plot(log10(metrics$lambda), metrics$loglik, type="l")
metrics$lambda[which.max(metrics$loglik)] # 0.3981072
plot(log10(metrics$lambda), metrics$MD_coefs, type="l")
plot(log10(metrics$lambda), metrics$MAD_coefs, type="l")
metrics$lambda[which.min(metrics$MAD_coefs)] # 2.511886
plot(log10(metrics$lambda), metrics$sMAPE_coefs, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(metrics$sMAPE_coefs)] # 2.511886
plot(log10(metrics$lambda), metrics$bias_coefs, type="l")
metrics$lambda[which.min(abs(metrics$bias_coefs-1))] # 0.001
plot(log10(metrics$lambda), metrics$R2_coefs, type="l")
plot(log10(metrics$lambda), metrics$estimated_over_true_support_size, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(abs(metrics$estimated_over_true_support_size-1))] # 1.584893
plot(log10(metrics$lambda), metrics$TP, type="l")
plot(log10(metrics$lambda), metrics$FP, type="l")
plot(log10(metrics$lambda), metrics$TN, type="l")
plot(log10(metrics$lambda), metrics$FN, type="l")
plot(log10(metrics$lambda), metrics$precision, type="l")
plot(log10(metrics$lambda), metrics$recall, type="l")
plot(log10(metrics$lambda), metrics$specificity, type="l")
plot(log10(metrics$lambda), metrics$accuracy, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.max(metrics$accuracy)] # 2.511886
plot(log10(metrics$lambda), metrics$F1_score, type="l")
metrics$lambda[which.max(metrics$F1_score)] # 2.511886
plot(log10(metrics$lambda), metrics$aic, type="l")
metrics$lambda[which.min(metrics$aic)] # 1
plot(log10(metrics$lambda), metrics$bic, type="l")
metrics$lambda[which.min(metrics$bic)] # 2.51
plot(log10(metrics$lambda), metrics$gic, type="l")
metrics$lambda[which.min(metrics$gic)] # 2.51
plot(log10(metrics$lambda), metrics$ebic, type="l")
metrics$lambda[which.min(metrics$ebic)] # 2.51
plot(log10(metrics$lambda), metrics$mbic, type="l")
metrics$lambda[which.min(metrics$mbic)] # 3.98

# reasonable match with lambda values that were expected to optimise each IC :
# note: using warm starts in the L0glm.fit function would shift the
# empirically observed optimal values

preset.lambda(IC = "aic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_gaus$y), 
              X = Matrix(data_gaus$X, sparse=TRUE), 
              family = gaussian(identity),
              dispersion = dispersion) # 0.5
preset.lambda(IC = "bic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_gaus$y), 
              X = Matrix(data_gaus$X, sparse=TRUE), 
              family = gaussian(identity),
              dispersion = dispersion) # 2.30
preset.lambda(IC = "gic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_gaus$y), 
              X = Matrix(data_gaus$X, sparse=TRUE), 
              family = gaussian(identity),
              dispersion = dispersion) # 2.56
preset.lambda(IC = "ebic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_gaus$y), 
              X = Matrix(data_gaus$X, sparse=TRUE), 
              family = gaussian(identity),
              dispersion = dispersion) # 2.30
preset.lambda(IC = "mbic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_gaus$y), 
              X = Matrix(data_gaus$X, sparse=TRUE), 
              family = gaussian(identity),
              dispersion = dispersion) # 3.91



# 2. POISSON LOG LINK GLM ####

set.seed(1)
data_pois  = simulate_glm(n, p, n_nonzero, family = "poisson") 
X = data_pois$X
y = data_pois$y

# glm fit
system.time(glm_pois <- glm(data_pois$y ~ 0+data_pois$X, family=poisson(log))) # 0.37s
glm_pois$deviance 
glm_pois$aic
AIC(glm_pois)
BIC(glm_pois)
logLik(glm_pois)

# L0 glm fit with lambda close to zero
system.time(L0glm_pois_lambda0 <- L0glm.fit(X = Matrix(data_pois$X, sparse=TRUE), 
                                            y = data_pois$y, 
                                            family = poisson(log), 
                                            lambda = 1, 
                                            nonzero.centered = TRUE, # use nonzero centered ridge regression?
                                            nonnegative = FALSE,
                                            warmupiters = 4,
                                            # start = rep(1, ncol(data_pois$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                                            solver = "eigen", # 
                                            method = 3L, #
                                            thresh = 1E-10
))  #
# TODO: fix bugs if start given: Error in weights[good] : invalid subscript type 'S4'
# if start not given : Error in X[good, nonzero, drop = F] :
# invalid or not-yet-implemented 'Matrix' subsetting

max(abs(L0glm_pois_lambda0$coefficients-coef(glm_pois))) # 
# check if metrics (loglik & AIC & BIC) are correct
metrics_pois = data.frame(do.call(cbind, get.metrics(beta = L0glm_pois_lambda0$coefficients,  # some metrics
                                                     X = Matrix(data_pois$X, sparse=TRUE), 
                                                     y = as.matrix(data_pois$y), 
                                                     family = gaussian(identity), 
                                                     beta_true = as.matrix(data_pois$beta_true))))
logLik(glm_pois)-metrics_pois$loglik # 
AIC(glm_pois)-metrics_pois$aic # 
BIC(glm_pois)-metrics_pois$bic # 


# L0 glm fit with lambda preset to optimise given IC
lambda_preset_pois = preset.lambda(IC = crit, # options aic, bic, gic, ebic, mbic 
                                   y = data_pois$y, 
                                   X = Matrix(data_pois$X, sparse=TRUE), 
                                   family = gaussian(identity),
                                   dispersion = dispersion) # dispersion here was 1, here we pass this value to not have to estimate it 
lambda_preset_pois # 2.556246

system.time(L0glm_pois_presetlambda <- L0glm.fit(X = Matrix(data_pois$X, sparse=TRUE), 
                                                 y = data_pois$y, 
                                                 family = poisson(log), 
                                                 lambda = lambda_preset_pois,
                                                 nonzero.centered = TRUE, # use nonzero centered ridge regression?
                                                 nonnegative = FALSE,
                                                 warmupiters = 4,
                                                 # start = rep(1, ncol(data_pois$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                                                 solver = "eigen", # 
                                                 method = 3L, #
                                                 thresh = 1E-10
)) # 
# TODO: fix bugs if start given: Error in weights[good] : invalid subscript type 'S4'
# if start not given : Error in X[good, nonzero, drop = F] :
# invalid or not-yet-implemented 'Matrix' subsetting

# coefficients & metrics over a range of lambda values
# calculate max lambda where all coefs would be penalized to 0
weights <- L0glm_pois_presetlambda$weights
Xw <- Diagonal(x=sqrt(weights)) %*% data_pois$X
sel_task <- which.max(colMeans(abs(as.matrix(data_pois$y)))) # for multitask select column with max mean abs value?
yw <- as.matrix(data_pois$y)[, sel_task,drop=F] * sqrt(weights) 
# TODO note that yw should really be zw, i.e. for non-identity link you should 
# be working on the adjusted response scale - maybe you should also make L0glm.fit 
# return z? Or make a convenience function that returns z and w given a particular 
# given family and coefficient estimates?
maxlambda <- max((1/colSums(Xw^2)) * (crossprod(Xw, yw)^2), na.rm=T)
if (family$family == "gaussian") maxlambda <- maxlambda/4 # could be that this factor of 1/4 applies to all families & that this is that 1/4 factor mentioned in Frommlet & Nuel
mincoef <- min(abs(L0glm_pois_presetlambda$coefficients)) # *sign(L0glm_pois_presetlambda$coefficients)

lambdas <- 10^seq(round(log10(lambda_preset_pois))-3, log10(maxlambda)+1, by=0.2)
coefmat <- Matrix(0, nrow=length(lambdas), ncol=ncol(data_pois$X))
warmstart <- 0.1 # % by which L0glm_pois_presetlambda$coefficients will be used as start, lowest value that led to stable coefficient paths
system.time(metrics <- data.frame(do.call(rbind, lapply(1:length(lambdas), function (i) { 
  lambda <- lambdas[[i]]
  start <- L0glm_pois_presetlambda$coefficients*warmstart+rep(mincoef/2,ncol(data_pois$X))*(1-warmstart)
  # start[L0glm_pois_presetlambda$coefficients==0] <- 0
  L0glm.out <- L0glm.fit(X = Matrix(data_pois$X, sparse=TRUE), 
                         y = data_pois$y, 
                         family = poisson(log), 
                         lambda = lambda, 
                         nonzero.centered = TRUE, # use nonzero centered ridge regression? - TODO: change that if you specify starting values start that it does not keep on shrinking coefs to those values but instead shrink to coef values from previous iteration
                         nonnegative = FALSE,
                         warmupiters = 1,
                         # start = start,
                         solver = "eigen", # eigen, glmnet or osqp, # note: constrained solvers glmnet & osqp better than eigen for very small lambda values
                         method = 3L, 
                         # gamma = 2, 
                         thresh = 1E-10,
                         relaxed = FALSE # TODO check : setting this argument to TRUE doesn't affect solution, while it should (it should do an L0glm fit & then refit an unpenalized LS model for the selected coefficients)
  )$coefficients 
  coefmat[i,] <<- as.vector(L0glm.out)
  return( unlist(c(lambda=lambda, get.metrics(beta = L0glm.out, 
                                              X = X, y = y, 
                                              family = gaussian(identity), 
                                              beta_true = data_pois$beta_true))) )
} )))) # 
matplot(x=log10(metrics$lambda), 
        y=as.matrix(coefmat), type="l", lty=1,
        xlab="Log10(lambda)", ylab="Coefficient estimates", main="Regularisation path L0 Gaussian GLM")
abline(v=log10(maxlambda), col="red") # max lambda where all coefs would be penalized to 0

plot(log10(metrics$lambda), metrics$deviance, type="l")
metrics$lambda[which.min(metrics$deviance)] # XXX # NOTE this is deviance relative to observed signal, not deviance relative to ground truth though
plot(log10(metrics$lambda), metrics$loglik, type="l")
metrics$lambda[which.max(metrics$loglik)] # XXX
plot(log10(metrics$lambda), metrics$MD_coefs, type="l")
plot(log10(metrics$lambda), metrics$MAD_coefs, type="l")
metrics$lambda[which.min(metrics$MAD_coefs)] # XXX
plot(log10(metrics$lambda), metrics$sMAPE_coefs, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(metrics$sMAPE_coefs)] # XXX
plot(log10(metrics$lambda), metrics$bias_coefs, type="l")
metrics$lambda[which.min(abs(metrics$bias_coefs-1))] # XXX
plot(log10(metrics$lambda), metrics$R2_coefs, type="l")
plot(log10(metrics$lambda), metrics$estimated_over_true_support_size, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(abs(metrics$estimated_over_true_support_size-1))] # XXX
plot(log10(metrics$lambda), metrics$TP, type="l")
plot(log10(metrics$lambda), metrics$FP, type="l")
plot(log10(metrics$lambda), metrics$TN, type="l")
plot(log10(metrics$lambda), metrics$FN, type="l")
plot(log10(metrics$lambda), metrics$precision, type="l")
plot(log10(metrics$lambda), metrics$recall, type="l")
plot(log10(metrics$lambda), metrics$specificity, type="l")
plot(log10(metrics$lambda), metrics$accuracy, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.max(metrics$accuracy)] # XXX
plot(log10(metrics$lambda), metrics$F1_score, type="l")
metrics$lambda[which.max(metrics$F1_score)] # XXX
plot(log10(metrics$lambda), metrics$aic, type="l")
metrics$lambda[which.min(metrics$aic)] # XXX
plot(log10(metrics$lambda), metrics$bic, type="l")
metrics$lambda[which.min(metrics$bic)] # XXX
plot(log10(metrics$lambda), metrics$gic, type="l")
metrics$lambda[which.min(metrics$gic)] # XXX
plot(log10(metrics$lambda), metrics$ebic, type="l")
metrics$lambda[which.min(metrics$ebic)] # XXX
plot(log10(metrics$lambda), metrics$mbic, type="l")
metrics$lambda[which.min(metrics$mbic)] # XXX

# reasonable match with lambda values that were expected to optimise each IC :
# note: using warm starts in the L0glm.fit function would shift the
# empirically observed optimal values

preset.lambda(IC = "aic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_gaus$y), 
              X = Matrix(data_gaus$X, sparse=TRUE), 
              family = gaussian(identity),
              dispersion = dispersion) # XXX
preset.lambda(IC = "bic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_gaus$y), 
              X = Matrix(data_gaus$X, sparse=TRUE), 
              family = gaussian(identity),
              dispersion = dispersion) # XXX
preset.lambda(IC = "gic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_gaus$y), 
              X = Matrix(data_gaus$X, sparse=TRUE), 
              family = gaussian(identity),
              dispersion = dispersion) # XXX
preset.lambda(IC = "ebic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_gaus$y), 
              X = Matrix(data_gaus$X, sparse=TRUE), 
              family = gaussian(identity),
              dispersion = dispersion) # XXX
preset.lambda(IC = "mbic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_gaus$y), 
              X = Matrix(data_gaus$X, sparse=TRUE), 
              family = gaussian(identity),
              dispersion = dispersion) # XXX










