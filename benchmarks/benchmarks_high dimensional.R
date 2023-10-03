# Test high dimensional variable selection ####

# fix to set working package repository when using Microsoft R Open with Intel MKL
options (repos = c(CRAN="https://cloud.r-project.org/"))
library(remotes)
# remotes::install_github("tomwenseleers/L0glm/L0glm")
# library(L0glm)
devtools::load_all(".")
library(export)
source("..//benchmarks/generate data.R")


# 1. HIGH DIMENSIONAL EXAMPLE IN ORDINIS PACKAGE ####

# In the example below MCP seems to give a better solution than L0 glm :-(

# 1.1 nonnegative MCP penalized regression calculated using ordinis ####

# devtools::install_github("jaredhuling/ordinis")
library(ordinis)

# compute the full solution path, n > p
set.seed(123)
n <- 500
p <- 1000000
n_nonzero <- 50 # nr of nonzero coefs
sd_noise <- 3
beta_true <- matrix(c(runif(n_nonzero), rep(0, p - n_nonzero))) # 50 000 coefs
X <- matrix(rnorm(n * p, sd = sd_noise), n, p)
y <- drop(X %*% beta_true) + rnorm(n)

system.time(mod <- ordinis(X, y, 
                           penalty = "mcp", # using MCP penalty
                           intercept = TRUE,
                           lower.limits = rep(0, p), # force all coefficients to be positive
                           penalty.factor = c(0, 0, rep(1, p-2)), # don't penalize first two coefficients
                           gamma = 1.1, # ifelse(penalty == "scad", 3.7, 1.4)
                           alpha = 0.95))  # use elastic net with alpha = 0.95
# 3.28s for n=500 & p=50000, 67s for n=500 & p=1000000
plot(mod)
logLik(mod) # 
max(logLik(mod)) # -731.8202 for SCAD with gamma=3.7, -704.5205 for MCP with gamma=1.4, -692.5572 for MCP with gamma=1.1
# min(logLik(mod)) # -1960.879
mod$lambda[which.max(logLik(mod))] # best lambda = 0.1869138
min(AIC(mod)) # 1485.114
min(BIC(mod)) # 1695.845
coefs_ordinis = mod$beta[,which.min(BIC(mod))]
coefs_ordinis[coefs_ordinis<0] = 0
sum(coefs_ordinis!=0) # 49, closely matches expected value with n_nonzero=50, but not e.g. with n_nonzero=200
library(scales)
abline(h=as.vector(beta_true[beta_true!=0]), col=alpha(1:6,0.2)) # note colours don't match the ones of ordinis

plot(x=coefs_ordinis, y=c(0,beta_true), pch=16, col="blue")



# 1.2 L0glm solution ####

Xscaled = normalizeX(X) # we include intercept & L2 norm normalize the columns of X
scaling = attr(Xscaled, "scaling")
beta_true_scaled = beta_true * scaling
family = gaussian(identity)
dispersion = sd_noise^2
crit = "bic" 
# note: BIC only works well with good warm start, e.g. based on
# MPC prescreening, which makes problem low dimensional again,
# otherwise use mBIC, appropriate for high dimensional settings 
# (that one is supposed to control the FWER)

# using lambda expected to optimise given IC
lambda_preset = preset.lambda(IC = crit, # options aic, bic, gic, ebic, mbic 
                              y = y, 
                              X = X, 
                              family = family, 
                              dispersion = dispersion) 
# with dispersion = NULL: estimating dispersion here doesn't seem to work well & 
# takes a very long time
# TODO : fix : dispersion=NULL doesn't work
lambda_preset # 13.98287, lambda that would be expected to optimise the expected value of the given IC
# TODO extremely small lambda_preset returned with dispersion = NULL

warm_start = 0.9 # 1 would essentially use MCP as pre-screen
# 0.9-1 all work well, but 1 should be faster as it should kick out 0 coefficients from the start
# TODO check in code that it is doing this in an optimal way, from the very start

# TODO: add warm_start as option in L0glm.fit & use this like below if start!=NULL
gamma = 2
start = coefs_ordinis[-1]*scaling # TODO incorporate scaling if normalize=T in L0glm.fit.internal
start = warm_start*start + (1-warm_start)*rep(sd(start)^(1/gamma), ncol(X))

# A warm_start value of 1 would indicate using the provided starting coefficients 
# directly (full warm start), while a warm_start between 0 and 1 would represent 
# a "cooled" or "tempered" start, blending the provided coefficients with  
# their standard deviations raised to the power of 1/gamma. In this way, one can
# specify one's prior belief in the starting coefficients being correct vs.
# implying identical priors on al coefficients (but still allowing 
# the observed variability in the starting coefficients to affect the chosen
# constant SD on the implied Gaussian prior of all coefficients).
system.time(L0glm.out_preset <- L0glm.fit(Xscaled, y, 
                                          family = family, 
                                          lambda = lambda_preset, 
                                          # nonzero.centered = FALSE, # use nonzero centered ridge regression?
                                          nonnegative = TRUE,
                                          lower = rep(0, ncol(X)),
                                          warmupiters = 1,
                                          no.pen = NULL,
                                          # start = rep(sqrt(0.05), ncol(X)),
                                          start = start, # rep(sqrt(50), ncol(X)), # NULL, # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                                          solver = "eigen", # note: solver osqp better than eigen for very small lambda values))
                                          method = 7L,
                                          verbose = T,
                                          thresh = 1E-20, # had to set this to low value here # TODO thresh should also be multiplied by scaling if normalize=T
                                          maxit = 20,
                                          tol = 1E-5,
                                          gamma = gamma 
)) 
# glmnet: 6s with n=500 and p=50000, but runs out of memory with n=500 & p=1000000; biglasso or picasso would be good alternative
# eigen: 4.4s, 4.0s, 3.3s, 3.1, 4.4s, 4.5s, 3.4s, 3.8s for solvers 1 to 8  with n=500 & p=500000
# eigen: 78s for solvers 3 & 4 (LSCG) with n=500 & p=1000000
# NOTE: eigen solvers 5, 6, 7 & 8 & glmnet all give an excellent solution
# eigen solver 3 & 4 seems to give back the initial values
# eigen solver 1 & 2 slightly lower quality solution
# fails for osqp 
# TODO: check in code that with warmstart=1 the coefficients with abs(start)<thresh
# are kicked out from the very beginning
L0glm.out_preset_coefs = L0glm.out_preset$coefficients/scaling # TODO do this scaling automatically if normalize=T
sum(L0glm.out_preset_coefs!=0) # 47
plot(x=L0glm.out_preset_coefs, y=beta_true, pch=16, col="blue")
# NOTE: a BIC pre-picked lambda here worked well since the MCP
# pre-screening already changed this from a high dimensional to a low dimensional
# problem
# Using MCP, SCAD and LASSO coefficients as warm start all work, but MCP here 
# gives a slightly better solution than SCAD or LASSO
# (using MCP, SCAD & LASSO as warm start returns 47, 46 & 42 nonzero coefficients)
# LASSO is quite good in that there is widely available fast solvers available,
# e.g. glmnet and for the unconstrained case biglasso
# ALASSO as warm start did not work well - only 37 nonzero coefs then & solution
# much worse
# Also note that with warm_start=0 a constrained solver like glmnet here performs
# better than the unconstrained Eigen solvers. With warm_start=0.9 performance of
# both is near-identical though.

# coefficients & metrics over a range of lambda values
# calculate max lambda where all coefs would be penalized to 0
maxlam <- maxlambda(X=Xscaled, weights=L0glm.out_preset$weights, 
                    z=y, no.pen=NULL, coefs=L0glm.out_preset$coefficients,
                    family=family)
maxlam #

# lambdas <- 10^seq(round(log10(lambda_preset))-1, log10(maxlam), by=0.2)
lambdas <- 10^seq(round(log10(lambda_preset))-2, log10(lambda_preset)+1, by=0.1)
coefmat <- Matrix(0, nrow=length(lambdas), ncol=ncol(X))
warm_start <- 0.9 # % by which L0glm.out_preset will be used as start, e.g. 0.9,
# or as a pre-screen warm_start could be set 1
system.time(metrics <- data.frame(do.call(rbind, lapply(1:length(lambdas), function (i) { 
  lambda <- lambdas[[i]]
  gamma = 2
  start = coefs_ordinis[-1]*scaling
  # start = L0glm.out_preset_coefs*scaling # TODO incorporate scaling if normalize=T in L0glm.fit.internal
  start = warm_start*start + (1-warm_start)*rep(sd(start)^(1/gamma), ncol(X))
  # start <- L0glm.out_preset_coefs*warmstart+start_coefs(X, family)*(1-warmstart)
  # start[L0glm.out_preset_coefs==0] <- 0
  L0glm.out <- L0glm.fit(Xscaled, y, 
                         family = gaussian(identity), 
                         lambda = lambda, 
                         # nonzero.centered = FALSE, # use nonzero centered ridge regression?
                         # ridge.center = NULL,
                         nonnegative = TRUE,
                         lower = rep(0, ncol(X)),
                         warmupiters = 1,
                         no.pen = NULL,
                         start = start, # NULL, # rep(0.1, ncol(sim$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                         solver = "eigen", # note: solver osqp better than eigen for very small lambda values))
                         method = 7L,
                         relaxed = FALSE,
                         thresh = 1E-20, # had to set this to low value here # TODO thresh should also be multiplied by scaling if normalize=T
                         maxit = 20,
                         tol = 1E-5,
                         gamma = gamma
  )$coefficients/scaling
  coefmat[i,] <<- as.vector(L0glm.out)
  return( unlist(c(lambda=lambda, get.metrics(beta = L0glm.out, 
                                              X = X, y = y, 
                                              family = gaussian(identity), 
                                              beta_true = beta_true))) )
} )))) # 126s
matplot(x=log10(metrics$lambda), y=as.matrix(coefmat), type="l", lty=1,
        xlab="Log10(lambda)", ylab="Coefficient estimates", main="Regularisation path L0 Gaussian GLM")
abline(v=log10(maxlam), col="red") # max lambda where all coefs would be penalized to 0
library(scales)
abline(h=as.vector(beta_true[beta_true!=0]), col=alpha(1:6,0.2))

# coefficients at best lambda, here the one optimising the chosen IC
# best = which.min(abs(metrics$estimated_over_true_support_size-1))
best = which.min(metrics[,"mbic"]) # here using mBIC worked best
L0glm.out_best_coefs = coefmat[best,]
plot(L0glm.out_best_coefs, y=beta_true, pch=16, col="blue")

# NOTE: for high dimensional problems with an extremely sparse solution
# pre-picking variables based on an MPC, SCAD or LASSO regression works well
# and gives a better solution than using a Gaussian non-adaptive prior to
# start the algorithm. This is because of the heavier tail in the implied
# prior distribution, e.g. Laplace instead of Gaussian in the case of LASSO.
# For high dimensional problems with a less sparse solution, a Gaussian non-adaptive prior
# might be expected to perform better than MCP, SCAD or LASSO though (e.g. with higher n_nonzero value).



plot(log10(metrics$lambda), metrics$deviance, type="l")
metrics$lambda[which.min(metrics$deviance)] # 1 # NOTE this is deviance relative to observed signal, not deviance relative to ground truth though
plot(log10(metrics$lambda), metrics$deviance_true, type="l")
metrics$lambda[which(metrics$deviance_true==min(metrics$deviance_true))] # 1 # deviance relative to ground truth
plot(log10(metrics$lambda), metrics$loglik, type="l")
abline(v=log10(maxlam), col="red") # max lambda where all coefs would be penalized to 0
metrics$lambda[which.max(metrics$loglik)] # 1
plot(log10(metrics$lambda), metrics$MD_coefs, type="l")
plot(log10(metrics$lambda), metrics$MAD_coefs, type="l")
metrics$lambda[which.min(metrics$MAD_coefs)] # 50.11872
plot(log10(metrics$lambda), metrics$sMAPE_coefs, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(metrics$sMAPE_coefs)] # 50.11872
plot(log10(metrics$lambda), metrics$bias_coefs, type="l")
metrics$lambda[which.min(abs(metrics$bias_coefs-1))] # 125.8925
plot(log10(metrics$lambda), metrics$R2_coefs, type="l")
plot(log10(metrics$lambda), metrics$estimated_over_true_support_size, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(abs(metrics$estimated_over_true_support_size-1))] # 1
plot(log10(metrics$lambda), metrics$TP, type="l")
plot(log10(metrics$lambda), metrics$FP, type="l")
plot(log10(metrics$lambda), metrics$TN, type="l")
plot(log10(metrics$lambda), metrics$FN, type="l")
plot(log10(metrics$lambda), metrics$precision, type="l")
plot(log10(metrics$lambda), metrics$recall, type="l")
plot(log10(metrics$lambda), metrics$specificity, type="l")
plot(log10(metrics$lambda), metrics$accuracy, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.max(metrics$accuracy)] # 50.11872
plot(log10(metrics$lambda), metrics$F1_score, type="l")
metrics$lambda[which.max(metrics$F1_score)] # 50.11872
plot(log10(metrics$lambda), metrics$aic, type="l")
metrics$lambda[which.min(metrics$aic)] # 1
plot(log10(metrics$lambda), metrics$bic, type="l")
metrics$lambda[which.min(metrics$bic)] # 1 # ca 10x higher than expected value
plot(log10(metrics$lambda), metrics$ebic, type="l")
metrics$lambda[which.min(metrics$ebic)] # 50.11872
plot(log10(metrics$lambda), metrics$gic, type="l")
metrics$lambda[which.min(metrics$gic)] # 50.11872
abline(v=log10(metrics$lambda[which.min(metrics$gic)]), col="blue") 
plot(log10(metrics$lambda), metrics$mbic, type="l")
metrics$lambda[which(metrics$mbic==min(metrics$mbic))] # 63.09573

# lambda values that a priori are expected to optimise each IC (with start = NULL)
# note: using warm starts in the L0glm.fit function would shift the
# empirically observed optimal values
preset.lambda(IC = "aic", # options aic, bic, gic, ebic, mbic 
              y = y, 
              X = X, 
              family = gaussian(identity),
              dispersion = dispersion) # 0.5
preset.lambda(IC = "bic", # options aic, bic, gic, ebic, mbic 
              y = y, 
              X = X, 
              family = gaussian(identity),
              dispersion = dispersion) # 1.151293
preset.lambda(IC = "ebic", # options aic, bic, gic, ebic, mbic 
              y = y, 
              X = X, 
              family = gaussian(identity),
              dispersion = dispersion) # 2.302585
preset.lambda(IC = "gic", # options aic, bic, gic, ebic, mbic 
              y = y, 
              X = X, 
              family = gaussian(identity),
              dispersion = dispersion) # 1.758231
preset.lambda(IC = "mbic", # options aic, bic, gic, ebic, mbic 
              y = y, 
              X = X, 
              family = gaussian(identity),
              dispersion = dispersion) # 2.76073




# comparison of metrics

# metrics for L0 glm fit with preset lambda
do.call(cbind, get.metrics(beta = L0glm.out_preset_coefs, 
                           X = Xscaled, y = y, 
                           family = family, 
                           beta_true = as.vector(beta_true_scaled)))
#          p   n m deviance deviance_true   loglik     MD_coefs  MAD_coefs sMAPE_coefs bias_coefs R2_coefs
# [1,] 50001 500 1 11794.78      11676.16 48453.33 -0.009864239 0.02036309  0.00206614   1.103877 0.737787
#      support_size_true support_size_estimated estimated_over_true_support_size TP FP    TN FN precision recall
# [1,]                50                     47                             0.94 24 23 49927 26 0.5106383   0.48
#      specificity accuracy  F1_score       gic      ebic       bic       aic      mbic
# [1,]   0.9995395  0.99902 0.4948454 -95957.86 -95867.96 -96608.36 -96810.66 -95702.74

# metrics for L0 glm fit with best fit criterion (e.g. best BIC, support size closest to ground truth, etc)
do.call(cbind, get.metrics(beta = L0glm.out_best_coefs, 
                           X = Xscaled, y = y, 
                           family = family, 
                           beta_true = as.vector(beta_true_scaled)))
#          p   n m deviance deviance_true  loglik     MD_coefs   MAD_coefs sMAPE_coefs bias_coefs R2_coefs
# [1,] 50001 500 1 3502.147      3157.237 48754.9 -0.004455869 0.008630098 0.001376975   1.046205  0.93792
#      support_size_true support_size_estimated estimated_over_true_support_size TP FP    TN FN precision recall
# [1,]                50                     49                             0.98 33 16 49934 17 0.6734694   0.66
#      specificity accuracy  F1_score       gic      ebic       bic      aic      mbic
# [1,]   0.9996797  0.99934 0.6666667 -96521.46 -96427.82 -97199.07 -97409.8 -96255.72

# metrics for ordinis nonnegative MCP fit
do.call(cbind, get.metrics(beta = coefs_ordinis, 
                           X = cbind(1, X), y = y, 
                           family = family, 
                           beta_true = c(0, beta_true)) )
#          p   n m deviance deviance_true   loglik     MD_coefs    MAD_coefs  sMAPE_coefs bias_coefs  R2_coefs
# [1,] 50002 500 1 490.1997      77.88129 49246.48 -9.47124e-06 1.690956e-05 0.0002676603   1.014516 0.9989278
#      support_size_true support_size_estimated estimated_over_true_support_size TP FP    TN FN precision recall
# [1,]                50                     50                                1 47  3 49948  3      0.94   0.94
#      specificity accuracy F1_score       gic      ebic       bic       aic      mbic
# [1,]   0.9999399  0.99988     0.94 -97484.85 -97389.34 -98176.01 -98390.96 -97213.79



# 2. TEST HIGH DIMENSIONAL SUPERRESOLUTION DECONVOLUTION ####

## 2.1 Single channel deconvolution assuming Gaussian noise ####
n = 1000
p = 10000
s <- 0.1 # sparsity, nr of timepoints with a peak 
k <- round(n*s)
mean_beta <- 2000
sd_logbeta <- 1
sd_noise <- 5 # SD of noise for Gaussian noise case
crit <- "gic" # choose lambda based on what would be expected to optimise this IC
# options in order of increasing regularisation: aic, bic, gic, ebic, mbic
family = gaussian(identity)


set.seed(1)
sim <- simulate_spike_train(n = n, 
                            p = p,
                            k = k,
                            family = "gaussian", 
                            mean_beta = mean_beta,
                            sd_logbeta = sd_logbeta,
                            sd_noise = sd_noise)
X <- sim$X
y <- sim$y
family <- gaussian(identity)

# using lambda expected to optimise given IC
lambda_preset = preset.lambda(IC = crit, # options aic, bic, gic, ebic, mbic 
                              y = y, 
                              X = X, 
                              family = family) 
lambda_preset # lambda that would be expected to optimise given IC

L0glm.out_preset <- L0glm.fit(X, y, 
                              family = family, 
                              lambda = lambda_preset, 
                              nonzero.centered = FALSE, # use nonzero centered ridge regression?
                              nonnegative = TRUE,
                              warmupiters = 1,
                              start = NULL, # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                              solver = "eigen", # note: solver osqp better than eigen for very small lambda values))
                              method = 7L # Cholesky
)
L0glm.out_preset_coefs = L0glm.out_preset$coefficients
points((1:length(L0glm.out_preset_coefs))*(length(sim$x)/length(L0glm.out_preset_coefs)), 
       L0glm.out_preset_coefs, pch = 16, col= 'blue')

# some metrics
do.call(cbind, get.metrics(beta = L0glm.out_preset_coefs, 
                           X = X, y = y, 
                           family = family, 
                           beta_true = sim$beta_true))

# coefficients & metrics over a range of lambda values
# calculate max lambda where all coefs would be penalized to 0
maxlam <- maxlambda(X=X, weights=L0glm.out_preset$weights, 
                    z=y, no.pen=NULL, coefs=L0glm.out_preset$coefficients,
                    family=family)
maxlam #

lambdas <- 10^seq(round(log10(lambda_preset))-3, log10(maxlam)+1, by=0.05)
coefmat <- Matrix(0, nrow=length(lambdas), ncol=ncol(sim$X))
warmstart <- 0.1 # % by which L0glm.out_preset will be used as start, lowest value that led to stable coefficient paths
metrics <- data.frame(do.call(rbind, lapply(1:length(lambdas), function (i) { 
  lambda <- lambdas[[i]]
  start <- L0glm.out_preset_coefs*warmstart+start_coefs(X, family)*(1-warmstart)
  # start[L0glm.out_preset_coefs==0] <- 0
  L0glm.out <- L0glm.fit(X, y, 
                         family = gaussian(identity), 
                         lambda = lambda, 
                         nonzero.centered = FALSE, # use nonzero centered ridge regression?
                         # ridge.center = NULL,
                         nonnegative = TRUE,
                         warmupiters = 1,
                         warmupiters.nzc = 5,
                         start = NULL, # start, # rep(0.1, ncol(sim$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                         solver = "eigen", # note: solver osqp better than eigen for very small lambda values))
                         method = 7L,
                         relaxed = FALSE,
                         thresh = 1E-5,
  )$coefficients
  coefmat[i,] <<- as.vector(L0glm.out)
  return( unlist(c(lambda=lambda, get.metrics(beta = L0glm.out, 
                                              X = X, y = y, 
                                              family = gaussian(identity), 
                                              beta_true = sim$beta_true))) )
} )))
matplot(x=log10(metrics$lambda), y=as.matrix(coefmat), type="l", lty=1,
        xlab="Log10(lambda)", ylab="Coefficient estimates", main="Regularisation path L0 Gaussian GLM")
abline(v=log10(maxlam), col="red") # max lambda where all coefs would be penalized to 0
library(scales)
abline(h=as.vector(sim$beta_true), col=alpha(1:6,0.2))

plot(log10(metrics$lambda), metrics$deviance, type="l")
metrics$lambda[which.min(metrics$deviance)] # 1.412538 # NOTE this is deviance relative to observed signal, not deviance relative to ground truth though
plot(log10(metrics$lambda), metrics$deviance_true, type="l")
metrics$lambda[which(metrics$deviance_true==min(metrics$deviance_true))] # 1.412538 # deviance relative to ground truth
plot(log10(metrics$lambda), metrics$loglik, type="l")
abline(v=log10(maxlam), col="red") # max lambda where all coefs would be penalized to 0
metrics$lambda[which.max(metrics$loglik)] # 1.412538
plot(log10(metrics$lambda), metrics$MD_coefs, type="l")
plot(log10(metrics$lambda), metrics$MAD_coefs, type="l")
metrics$lambda[which.min(metrics$MAD_coefs)] # 79.43282
plot(log10(metrics$lambda), metrics$sMAPE_coefs, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(metrics$sMAPE_coefs)] # 79.43282
plot(log10(metrics$lambda), metrics$bias_coefs, type="l")
metrics$lambda[which.min(abs(metrics$bias_coefs-1))] # 1.995262
plot(log10(metrics$lambda), metrics$R2_coefs, type="l")
plot(log10(metrics$lambda), metrics$estimated_over_true_support_size, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(abs(metrics$estimated_over_true_support_size-1))] # 251.1886
plot(log10(metrics$lambda), metrics$TP, type="l")
plot(log10(metrics$lambda), metrics$FP, type="l")
plot(log10(metrics$lambda), metrics$TN, type="l")
plot(log10(metrics$lambda), metrics$FN, type="l")
plot(log10(metrics$lambda), metrics$precision, type="l")
plot(log10(metrics$lambda), metrics$recall, type="l")
plot(log10(metrics$lambda), metrics$specificity, type="l")
plot(log10(metrics$lambda), metrics$accuracy, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.max(metrics$accuracy)] # 79.43282
plot(log10(metrics$lambda), metrics$F1_score, type="l")
metrics$lambda[which.max(metrics$F1_score)] # 79.43282
plot(log10(metrics$lambda), metrics$aic, type="l")
metrics$lambda[which.min(metrics$aic)] # 1.412538
plot(log10(metrics$lambda), metrics$bic, type="l")
metrics$lambda[which.min(metrics$bic)] # 56.23413 # ca 10x higher than expected value
plot(log10(metrics$lambda), metrics$ebic, type="l")
metrics$lambda[which.min(metrics$ebic)] # 56.23413
plot(log10(metrics$lambda), metrics$gic, type="l")
metrics$lambda[which.min(metrics$gic)] # 56.23413
abline(v=log10(metrics$lambda[which.min(metrics$gic)]), col="blue") 
plot(log10(metrics$lambda), metrics$mbic, type="l")
metrics$lambda[which(metrics$mbic==min(metrics$mbic))] # 223.8721

# lambda values that a priori are expected to optimise each IC (with start = NULL)
# note: using warm starts in the L0glm.fit function would shift the
# empirically observed optimal values
preset.lambda(IC = "aic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(sim$y), 
              X = Matrix(sim$X, sparse=TRUE), 
              family = gaussian(identity),
              dispersion = dispersion) # 0.5
preset.lambda(IC = "bic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(sim$y), 
              X = Matrix(sim$X, sparse=TRUE), 
              family = gaussian(identity),
              dispersion = dispersion) # 1.151293
preset.lambda(IC = "ebic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(sim$y), 
              X = Matrix(sim$X, sparse=TRUE), 
              family = gaussian(identity),
              dispersion = dispersion) # 2.302585
preset.lambda(IC = "gic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(sim$y), 
              X = Matrix(sim$X, sparse=TRUE), 
              family = gaussian(identity),
              dispersion = dispersion) # 1.758231
preset.lambda(IC = "mbic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(sim$y), 
              X = Matrix(sim$X, sparse=TRUE), 
              family = gaussian(identity),
              dispersion = dispersion) # 2.76073






# 3. TEST HIGH DIMENSIONAL ORTHOGONAL PROBLEM ####

n = 100
p = 200
n_nonzero = 50
beta_sd = 1000
error_sd = 1
dispersion = error_sd^2
crit = "mbic"

data_gaus  = simulate_standardnormal(n=n, p=p, n_nonzero=n_nonzero, beta_sd=beta_sd, 
                                     error_sd=error_sd,
                                     seed.problem=1,
                                     seed.noise=1) 
X = normalizeX(data_gaus$X) # we first L2 norm normalize the columns of X
y = data_gaus$y
scaling = attr(X, "scaling")
beta_true = data_gaus$beta_true * scaling
family = gaussian(identity)


# glm fit
# system.time(glm_gaus <- glm.fit(y=y, x=X, intercept=FALSE, family=gaussian(identity))) 
# Error: cannot allocate vector of size XX Gb

# L0 glm fit with pre-chosen lambda chosen to optimise given IC
lambda_preset_gaus = preset.lambda(IC = crit, # options aic, bic, gic, ebic, mbic 
                                   y = y,
                                   X = X, 
                                   family = gaussian(identity),
                                   dispersion = dispersion # dispersion here was 1, here we pass this value to not have to estimate it 
) 
lambda_preset_gaus 

system.time(L0glm_gaus_presetlambda <- L0glm.fit(X = Matrix(X, sparse=T), 
                                                 y = y, 
                                                 family = gaussian(identity), 
                                                 lambda = lambda_preset_gaus,
                                                 # nonzero.centered = TRUE, # use nonzero centered ridge regression?
                                                 # ridge.center = rep(sqrt(100), ncol(X))/scaleX(data_gaus$X),
                                                 no.pen = 1,
                                                 nonnegative = FALSE,
                                                 normalize = FALSE,
                                                 warmupiters = 1,
                                                 start = rep(sqrt(100), ncol(X)), # rep(sqrt(100), ncol(X)), # /scaleX(data_gaus$X), NULL, # =rep(sqrt(5), ncol(data_gaus$X))/scaling, # NULL, # rep(sqrt(5), ncol(data_gaus$X))*scaling, # NULL, # start_coefs(X, family=gaussian(identity)), # = rep(sqrt(5), ncol(data_gaus$X)) / scaling, 
                                                 solver = "eigen", # 
                                                 method = 7L, #
                                                 relaxed = FALSE,
                                                 thresh = min(abs(beta_true[!beta_true==0])) # 1E-5
)) # XXs with method 1L
sum(L0glm_gaus_presetlambda$coefficients!=0)  # 51
plot(L0glm_gaus_presetlambda$coefficients, beta_true, pch=16, col="blue")
cor(L0glm_gaus_presetlambda$coefficients, beta_true) # 0.9940267
plot(x=L0glm_gaus_presetlambda$coefficients, y=data_gaus$beta_true, pch=16, col="blue")
points(x=L0glm_gaus_presetlambda$coefficients[data_gaus$beta_true==0], y=data_gaus$beta_true[data_gaus$beta_true==0], 
       pch=16, col="red")

# check metrics & quality of solution with preset lambda
system.time(metrics_gaus <- data.frame(do.call(cbind, get.metrics(beta = L0glm_gaus_presetlambda$coefficients,  # some metrics
                                                                  X = X, 
                                                                  y = y,
                                                                  family = gaussian(identity), 
                                                                  beta_true = beta_true)))) 
metrics_gaus
#     p     n m deviance deviance_true    loglik    MD_coefs MAD_coefs sMAPE_coefs bias_coefs  R2_coefs
# 1 102 10000 1 10005.95      80.24103 -14138.36 0.009655372 0.9534635   0.1049711   1.005573 0.9982282
#   support_size_true support_size_estimated estimated_over_true_support_size TP FP TN FN precision    recall
# 1                51                     47                        0.9215686 47  0 50  4         1 0.9215686
#   specificity accuracy  F1_score      gic     ebic      bic      aic     mbic
# 1           1 0.960396 0.9591837 28769.63 28720.71 28718.81 28372.72 29029.73


# orthogonal coefficient initialisation for possible use as start ####
beta_init_orthogonal <- init_coefs(X = X,
                                   mu = y,
                                   family = gaussian(identity),
                                   lambda = lambda_preset_gaus)
cor(beta_init_orthogonal, data_gaus$beta_true) # 0.981208
plot(beta_init_orthogonal, data_gaus$beta_true, pch=16, col="blue")
mean(abs(beta_init_orthogonal-data_gaus$beta_true))

# coefficients & metrics over a range of lambda values
# calculate max lambda where all coefs would be penalized to 0
maxlam <- maxlambda(X=X, 
                    weights=L0glm_gaus_presetlambda$weights, 
                    z=L0glm_gaus_presetlambda$z,
                    no.pen=1,
                    coefs=L0glm_gaus_presetlambda$coefficients,
                    family=family)
maxlam # 2346.189
mincoef <- min(abs(L0glm_gaus_presetlambda$coefficients))

lambdas <- 10^seq(round(log10(lambda_preset_gaus))-1, log10(maxlam), by=0.1)
coefmat <- Matrix(0, nrow=length(lambdas), ncol=ncol(X))
warmstart <- 0.2 # % by which L0glm_gaus_presetlambda$coefficients will be used as start
system.time(metrics <- data.frame(do.call(rbind, lapply(1:length(lambdas), function (i) { 
  lambda <- lambdas[[i]]
  start <- (L0glm_gaus_presetlambda$coefficients)*warmstart+start_coefs(X, family)*(1-warmstart)
  # start[L0glm_gaus_presetlambda$coefficients==0] <- 0
  
  L0glm.out <- L0glm.fit(X = X, 
                         y = y, 
                         family = gaussian(identity), 
                         lambda = lambda, 
                         no.pen = 1,
                         nonzero.centered = FALSE, # use nonzero centered ridge regression?
                         ridge.center = NULL, # L0glm_gaus_presetlambda$coefficients, # rep(0, ncol(data_gaus$X)),
                         nonnegative = FALSE,
                         warmupiters = 1,
                         warmupiters.nzc = 5,
                         start = start, # 
                         solver = "eigen", # eigen, glmnet or osqp, # note: constrained solvers glmnet & osqp better than eigen for very small lambda values
                         method = 3L, # 3L here 2x as fast than 7L
                         # gamma = 2, 
                         thresh = 1E-5,
                         relaxed = FALSE 
  )$coefficients 
  coefmat[i,] <<- as.vector(L0glm.out)
  return( unlist(c(lambda=lambda, get.metrics(beta = L0glm.out, 
                                              X = X, y = y, 
                                              family = gaussian(identity), 
                                              beta_true = beta_true))) )
} )))) # 19s
matplot(x=log10(metrics$lambda), 
        y=as.matrix(coefmat), type="l", lty=1,
        xlab="Log10(lambda)", ylab="Coefficient estimates", main="Regularisation path L0 Gaussian GLM")
library(scales)
abline(h=beta_true, col=alpha(1:6,0.2))
abline(v=log10(maxlam), col="red") # max lambda where all coefs would be penalized to 0

plot(log10(metrics$lambda), metrics$deviance, type="l")
metrics$lambda[which.min(metrics$deviance)] # 0.1 # NOTE this is deviance relative to observed signal, not deviance relative to ground truth though
plot(log10(metrics$lambda), metrics$deviance_true, type="l")
metrics$lambda[which.min(metrics$deviance_true)] # 0.7943282 # deviance relative to ground truth
plot(log10(metrics$lambda), metrics$loglik, type="l")
metrics$lambda[which.max(metrics$loglik)] # 0.3981072
plot(log10(metrics$lambda), metrics$MD_coefs, type="l")
plot(log10(metrics$lambda), metrics$MAD_coefs, type="l")
metrics$lambda[which.min(metrics$MAD_coefs)] # 1
plot(log10(metrics$lambda), metrics$sMAPE_coefs, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(metrics$sMAPE_coefs)] # 1.995262
plot(log10(metrics$lambda), metrics$bias_coefs, type="l")
metrics$lambda[which(abs(metrics$bias_coefs-1)==min(abs(metrics$bias_coefs-1)))] # 794.3282
plot(log10(metrics$lambda), metrics$R2_coefs, type="l")
metrics$lambda[which.max(abs(metrics$R2_coefs))] # 0.7943282
plot(log10(metrics$lambda), metrics$estimated_over_true_support_size, type="l") # good one, optimal for lambda between 1 & 10
range(metrics$lambda[which(abs(metrics$estimated_over_true_support_size-1)==min(abs(metrics$estimated_over_true_support_size-1)))]) # 1.258925 1.258925
plot(log10(metrics$lambda), metrics$TP, type="l")
plot(log10(metrics$lambda), metrics$FP, type="l")
plot(log10(metrics$lambda), metrics$TN, type="l")
plot(log10(metrics$lambda), metrics$FN, type="l")
plot(log10(metrics$lambda), metrics$precision, type="l")
plot(log10(metrics$lambda), metrics$recall, type="l")
plot(log10(metrics$lambda), metrics$specificity, type="l")
plot(log10(metrics$lambda), metrics$accuracy, type="l") # good one, optimal for lambda between 1 & 10
range(metrics$lambda[which(metrics$accuracy==max(metrics$accuracy))]) # 1.995262 1.995262
plot(log10(metrics$lambda), metrics$F1_score, type="l")
range(metrics$lambda[which(metrics$F1_score==max(metrics$F1_score))]) # 1.995262 1.995262
plot(log10(metrics$lambda), metrics$aic, type="l")
metrics$lambda[which(metrics$aic==min(metrics$aic))] # 1
plot(log10(metrics$lambda), metrics$bic, type="l")
metrics$lambda[which.min(metrics$bic)] # 1
plot(log10(metrics$lambda), metrics$ebic, type="l")
metrics$lambda[which.min(metrics$ebic)] # 1
plot(log10(metrics$lambda), metrics$gic, type="l")
metrics$lambda[which.min(metrics$gic)] # 1
plot(log10(metrics$lambda), metrics$mbic, type="l")
metrics$lambda[which(metrics$mbic==min(metrics$mbic))] # 1.995262

# with warmstart=0.2 ca 2x lower than the lambda tilde values that 
# were expected to optimise each IC; close match with expected values
# for warmstart=0
# note: as expected using warm starts in the L0glm.fit function shift the
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
preset.lambda(IC = "ebic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_gaus$y), 
              X = Matrix(data_gaus$X, sparse=TRUE), 
              family = gaussian(identity),
              dispersion = dispersion) # 2.30
preset.lambda(IC = "gic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_gaus$y), 
              X = Matrix(data_gaus$X, sparse=TRUE), 
              family = gaussian(identity),
              dispersion = dispersion) # 2.56
preset.lambda(IC = "mbic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_gaus$y), 
              X = Matrix(data_gaus$X, sparse=TRUE), 
              family = gaussian(identity),
              dispersion = dispersion) # 3.92
