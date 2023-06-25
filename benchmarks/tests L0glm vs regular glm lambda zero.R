# Check of metrics (AIC, BIC, loglik) L0glm package & comparison regular GLM ####

library(remotes)
# remotes::install_github("tomwenseleers/L0glm/L0glm")
library(L0glm)
library(export)
# devtools::load_all(".")
source("..//benchmarks/generate data.R")

# check consistency of L0glm output with glm for gaussian identity link, poisson log link, Gamma log link & binomial logit link
n = 10000
p = 100
n_nonzero = 100
dispersion = 1
crit = "gic"

# 1. GAUSSIAN GLM ####

set.seed(1)
data_gaus  = simulate_glm(n, p, n_nonzero, family = "gaussian") 
X = data_gaus$X
y = data_gaus$y

# glm fit
system.time(glm_gaus <- glm(data_gaus$y ~ 0+data_gaus$X, family=gaussian(identity))) # 1.83s
glm_gaus$deviance 
glm_gaus$aic
AIC(glm_gaus)
BIC(glm_gaus)
logLik(glm_gaus)

# L0 glm fit with lambda close to zero
system.time(L0glm_gaus_lambda0 <- L0glm.fit(Matrix(data_gaus$X, sparse=TRUE), as.matrix(data_gaus$y), 
                                            family = gaussian(identity), 
                                            lambda = 1E-30, 
                                            nonzero.centered = TRUE, # use nonzero centered ridge regression?
                                            nonnegative = FALSE,
                                            warmupiters = 1,
                                            start = rep(1, ncol(data_gaus$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                                            solver = "eigen", # 
                                            method = 3L, #
                                            thresh = 1E-10
))  #
max(abs(L0glm_gaus_lambda0$coefficients-coef(glm_gaus))) # 2.412845e-10, correct
# check if metrics (loglik & AIC & BIC) are correct
metrics_gaus = data.frame(do.call(cbind, get.metrics(beta = L0glm_gaus_lambda0$coefficients,  # some metrics
                                                     X = Matrix(data_gaus$X, sparse=TRUE), 
                                                     y = as.matrix(data_gaus$y), 
                                                     family = gaussian(identity), 
                                                     beta_true = as.matrix(data_gaus$beta_true))))
logLik(glm_gaus)-metrics_gaus$loglik # 0, correct
AIC(glm_gaus)-metrics_gaus$aic # 0, correct
BIC(glm_gaus)-metrics_gaus$bic # 0, correct


# L0 glm fit with lambda preset to optimise given IC
lambda_preset_gaus = preset.lambda(IC = crit, # options aic, bic, gic, ebic, mbic 
                              y = as.matrix(data_gaus$y), 
                              X = Matrix(data_gaus$X, sparse=TRUE), 
                              family = gaussian(identity),
                              dispersion = dispersion) # dispersion here was 1, here we pass this value to not have to estimate it 
lambda_preset_gaus # 2.556246

system.time(L0glm_gaus_presetlambda <- L0glm.fit(Matrix(data_gaus$X, sparse=TRUE), 
                                            as.matrix(data_gaus$y), 
                                            family = gaussian(identity), 
                                            lambda = lambda_preset_gaus,
                                            nonzero.centered = TRUE, # use nonzero centered ridge regression?
                                            nonnegative = FALSE,
                                            warmupiters = 1,
                                            start = rep(1, ncol(data_gaus$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                                            solver = "eigen", # 
                                            method = 3L, #
                                            thresh = 1E-10
)) # 1s

# coefficients & metrics over a range of lambda values
# calculate max lambda where all coefs would be penalized to 0
weights <- L0glm_gaus_presetlambda$weights
Xw <- Diagonal(x=sqrt(weights)) %*% data_gaus$X
sel_task <- which.max(colMeans(abs(as.matrix(data_gaus$y)))) # for multitask select column with max mean abs value?
yw <- as.matrix(data_gaus$y)[, sel_task,drop=F] * sqrt(weights) # TODO note that yw should really be zw, i.e. for non-identity link you should be working on the adjusted response scale - maybe you should also make L0glm.fit return z? Or make a convenience function that returns z and w given a particular given family and coefficient estimates?
maxlambda <- max((1/colSums(Xw^2)) * (crossprod(Xw, yw)^2), na.rm=T)
if (family$family == "gaussian") maxlambda <- maxlambda/4
mincoef <- min(abs(L0glm_gaus_presetlambda$coefficients)) # *sign(L0glm_gaus_presetlambda$coefficients)

lambdas <- 10^seq(round(log10(lambda_preset_gaus))-3, log10(maxlambda)+1, by=0.2)
coefmat <- Matrix(0, nrow=length(lambdas), ncol=ncol(data_gaus$X))
warmstart <- 0.1 # % by which L0glm_gaus_presetlambda$coefficients will be used as start, lowest value that led to stable coefficient paths
system.time(metrics <- data.frame(do.call(rbind, lapply(1:length(lambdas), function (i) { 
  lambda <- lambdas[[i]]
  start <- L0glm_gaus_presetlambda$coefficients*warmstart+rep(mincoef/2,ncol(data_gaus$X))*(1-warmstart)
  # start[L0glm_gaus_presetlambda$coefficients==0] <- 0
  L0glm.out <- L0glm.fit(Matrix(data_gaus$X, sparse=TRUE), 
                         as.matrix(data_gaus$y), 
                         family = gaussian(identity), 
                         lambda = lambda, 
                         nonzero.centered = FALSE, # use nonzero centered ridge regression? - TODO: change that if you specify starting values start that it does not keep on shrinking coefs to those values but instead shrink to coef values from previous iteration
                         nonnegative = FALSE,
                         warmupiters = 1,
                         start = start,
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
} )))) # 
matplot(x=log10(metrics$lambda), 
        y=as.matrix(coefmat), type="l", lty=1,
        xlab="Log10(lambda)", ylab="Coefficient estimates", main="Regularisation path L0 Gaussian GLM")
abline(v=log10(maxlambda), col="red") # max lambda where all coefs would be penalized to 0

plot(log10(metrics$lambda), metrics$deviance, type="l")
plot(log10(metrics$lambda), metrics$loglik, type="l")
plot(log10(metrics$lambda), metrics$MD_coefs, type="l")
plot(log10(metrics$lambda), metrics$MAD_coefs, type="l")
plot(log10(metrics$lambda), metrics$sMAPE_coefs, type="l") # good one, optimal for lambda between 1 & 10
plot(log10(metrics$lambda), metrics$bias_coefs, type="l")
plot(log10(metrics$lambda), metrics$R2_coefs, type="l")
plot(log10(metrics$lambda), metrics$estimated_over_true_support_size, type="l") # good one, optimal for lambda between 1 & 10
plot(log10(metrics$lambda), metrics$TP, type="l")
plot(log10(metrics$lambda), metrics$FP, type="l")
plot(log10(metrics$lambda), metrics$TN, type="l")
plot(log10(metrics$lambda), metrics$FN, type="l")
plot(log10(metrics$lambda), metrics$precision, type="l")
plot(log10(metrics$lambda), metrics$recall, type="l")
plot(log10(metrics$lambda), metrics$specificity, type="l")
plot(log10(metrics$lambda), metrics$accuracy, type="l") # good one, optimal for lambda between 1 & 10
plot(log10(metrics$lambda), metrics$F1_score, type="l")
plot(log10(metrics$lambda), metrics$aic, type="l")
metrics$lambda[which.min(metrics$aic)] # 1
plot(log10(metrics$lambda), metrics$bic, type="l")
metrics$lambda[which.min(metrics$bic)] # 1.584893
plot(log10(metrics$lambda), metrics$gic, type="l")
metrics$lambda[which.min(metrics$gic)] # 1.584893
plot(log10(metrics$lambda), metrics$ebic, type="l")
metrics$lambda[which.min(metrics$ebic)] # 1.584893
plot(log10(metrics$lambda), metrics$mbic, type="l")
metrics$lambda[which.min(metrics$mbic)] # 3.981072

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









