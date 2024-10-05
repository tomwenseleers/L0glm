# Test high dimensional variable selection - here I use MCP as a pre-screen on a large problem and then fit an L0glm on the pre-selected subset to remove some false positives & also test an adaptive L0glm algorithm ####

# fix to set working package repository when using Microsoft R Open with Intel MKL
options (repos = c(CRAN="https://cloud.r-project.org/"))
library(remotes)
# remotes::install_github("tomwenseleers/L0glm/L0glm")
# library(L0glm)
devtools::load_all(".")
library(export)
# source("..//benchmarks/generate data.R")
library(Matrix)


# 1. TEST HIGH DIMENSIONAL ORTHOGONAL PROBLEM ####

n = 500 # 500 # 1000
p = 1000000 # 100000 # 10000 # 1000-4000 works, between 5000 and 10000 solution breaks down - maybe pre-selecting variables with MCP might still work then
s = 0.00005 # 0.05, sparsity of true solution - better solution the sparser the true solution is
n_nonzero = round(p*s)
n_nonzero # 50
beta_sd = 1 # this would directly correspond to Cohen's d - 0.2, 0.5 or 0.8 would be considered small, medium or large
error_sd = 1 
dispersion = error_sd^2

set.seed(1)
data_gaus  = simulate_orthogonal(n=n, p=p, n_nonzero=n_nonzero, 
                                 beta_sd=beta_sd, 
                                 error_sd=error_sd) 
X = normalizeX(data_gaus$X) # we first L2 norm normalize the columns of X
Xsp = Matrix(X, sparse=T)
y = data_gaus$y
# X = data_gaus$X

# X = scale(data_gaus$X, center=T, scale=F)
# y = scale(data_gaus$y, center=T, scale=F)

scaling = attr(X, "scaling")
# beta_true = data_gaus$beta_true * scaling
beta_true = data_gaus$beta_true
L1norm_betatrue = sum(abs(beta_true))
L0norm_betatrue = sum(beta_true!=0)
L0overL1norm_betatrue = L0norm_betatrue/L1norm_betatrue
family = gaussian(identity)



# MCP fit using ordinis ####

library(ordinis)
# gamma_mcp = 2 # default for MCP=1.4
system.time(ordinis_fit <- ordinis(X, y, intercept=F,
               penalty = "mcp",
               # gamma = gamma_mcp, # using default gamma for MCP of 1.4
               # lower.limits = rep(0, p), # force all coefficients to be positive
               # penalty.factor = c(0, 0, rep(1, p-2)), # don't penalize first two coefficients
               alpha = 1 #,
               # maxit = 1, # could be set to 1 iteration when used for pre-screening
               # tol = 0.1
               )) # 34s with maxit 1, 64s otherwise
plot(ordinis_fit)
ordinis_fit$lambda[which.max(logLik(ordinis_fit))] # best lambda = 0.1302248
min(AIC(ordinis_fit)) #
min(BIC(ordinis_fit)) #
# coefs_ordinis = ordinis_fit$beta[,which.min(BIC(ordinis_fit))][-1]
# coefs_ordinis = ordinis_fit$beta[,which.min(AIC(ordinis_fit))][-1]
coefs_ordinis = ordinis_fit$beta[,which.max(logLik(ordinis_fit))][-1] # gives same result as AIC & BIC here
sum(coefs_ordinis!=0) # 57
sum(coefs_ordinis!=0&beta_true!=0) # 34 true pos
sum(coefs_ordinis!=0&beta_true==0) # 23 false positives still

beta_est0 = coefs_ordinis/scaling

table(beta_est0!=0, beta_true!=0, 
      dnn=c("estimated beta nonzero", "true beta nonzero"))
#                       true beta nonzero
# estimated beta nonzero  FALSE   TRUE
#                   FALSE 999927     16
#                   TRUE      23     34

cor(beta_est0, beta_true)  # 0.9822504

plot(x=beta_est0[beta_true!=0&beta_est0!=0], 
     y=beta_true[beta_true!=0&beta_est0!=0], pch=16, col="blue",
     xlab="Estimated coefficients", ylab="True coefficients",
     main="MCP penalized regression (ordinis)") # true positives
abline(a=0, b=1, col="lightgrey")
points(x=beta_est0[beta_true!=0&beta_est0==0], 
       y=beta_true[beta_true!=0&beta_est0==0], 
       pch=16, col="purple") # false negatives
points(x=beta_est0[beta_true==0&beta_est0!=0], 
       y=beta_true[beta_true==0&beta_est0!=0], 
       pch=16, col="red") # false positives
points(x=beta_est0[beta_true==0&beta_est0==0][[1]], 
       y=beta_true[beta_true==0&beta_est0==0][[1]], 
       pch=16, col="green") # true negatives

# check metrics & quality of MCP solution
system.time(metrics_MCP <- data.frame(do.call(cbind, get.metrics(beta = beta_est0,  # some metrics
                                                                  X = X, 
                                                                  y = y,
                                                                  family = gaussian(identity), 
                                                                  beta_true = beta_true)))) 
metrics_MCP
#         p   n m deviance deviance_true   loglik      MD_coefs    MAD_coefs  sMAPE_coefs bias_coefs
# 1 1000001 500 1 21401.21     0.7318924 998294.4 -5.635672e-07 8.210249e-06 8.102734e-05  0.9685639
#    R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true
# 1 0.9648158 -0.01079157    0.1083137      0.8601288     0.9685639    0.973887                50
#   support_size_estimated estimated_over_true_support_size TP FP     TN FN precision recall
# 1                     57                             1.14 34 23 999927 16 0.5964912   0.68
#   specificity accuracy F1_score      gic     ebic      bic      aic     mbic
# 1    0.999977 0.999961 0.635514 -1995125 -1994986 -1996228 -1996473 -1994787


# adaptive MCP fit ####

beta_est0B = beta_est0*scaling
penfact = 1/(abs(beta_est0B)+1E-5) 
# adaptive penalty based on initial MCP solution
# NOTE:  1/(abs(beta_est0B)^(1/4)+1E-5) returns 33 nonzero coefs and no false positives, 
# so slightly better, but that would introduce additional hyperparameter
penfact = length(penfact)*penfact/sum(penfact)

system.time(ordinis_fit2 <- ordinis(X, y, intercept=F,
                                   penalty = "mcp",
                                   gamma = 2, # default for MCP=1.4
                                   # lower.limits = rep(0, p), # force all coefficients to be positive
                                   penalty.factor = penfact, 
                                   alpha = 1 #,
                                   # maxit = 1, # could be set to 1 iteration when used for pre-screening
                                   #tol = 0.1
                                   )) # 37s with maxit 1, otherwise 60s
plot(ordinis_fit2, main="adaptive MCP")
ordinis_fit2$lambda[which.max(logLik(ordinis_fit2))] 
min(AIC(ordinis_fit2)) #
min(BIC(ordinis_fit2)) #
# coefs_ordinis2 = ordinis_fit2$beta[,which.min(BIC(ordinis_fit2))][-1]
# coefs_ordinis2 = ordinis_fit2$beta[,which.min(AIC(ordinis_fit2))][-1]
coefs_ordinis2 = ordinis_fit2$beta[,which.max(logLik(ordinis_fit2))][-1]
sum(coefs_ordinis2!=0) # 24
sum(coefs_ordinis2!=0&beta_true!=0) # 24 true pos
sum(coefs_ordinis2!=0&beta_true==0) # 0 false positives
# TODO This could also be calculated using an ordinis fit 
# using subsetted X, taking only cols with nonzero coefs, provided
# tuning criteria would be calculated as if we had been working with
# original problem size?

beta_est0B = coefs_ordinis2/scaling

table(beta_est0B!=0, beta_true!=0, 
      dnn=c("estimated beta nonzero", "true beta nonzero"))
# true beta               nonzero
# estimated beta nonzero  FALSE   TRUE
#                  FALSE 999950     26
#                  TRUE       0     24

cor(beta_est0B, beta_true) # 0.99194 

plot(x=beta_est0B[beta_true!=0&beta_est0B!=0], 
     y=beta_true[beta_true!=0&beta_est0B!=0], pch=16, col="blue",
     xlab="Estimated coefficients", ylab="True coefficients",
     main="adaptive MCP penalized regression (ordinis)") # true positives
abline(a=0, b=1, col="lightgrey")
points(x=beta_est0B[beta_true!=0&beta_est0B==0], 
       y=beta_true[beta_true!=0&beta_est0B==0], 
       pch=16, col="purple") # false negatives
points(x=beta_est0B[beta_true==0&beta_est0B!=0], 
       y=beta_true[beta_true==0&beta_est0B!=0], 
       pch=16, col="red") # false positives
points(x=beta_est0B[beta_true==0&beta_est0B==0][[1]], 
       y=beta_true[beta_true==0&beta_est0B==0][[1]], 
       pch=16, col="green") # true negatives

# check metrics & quality of adaptive MCP solution
system.time(metrics_adMCP <- data.frame(do.call(cbind, get.metrics(beta = beta_est0B,  # some metrics
                                                                 X = X, 
                                                                 y = y,
                                                                 family = gaussian(identity), 
                                                                 beta_true = beta_true)))) 
metrics_adMCP
#         p   n m deviance deviance_true   loglik      MD_coefs    MAD_coefs  sMAPE_coefs bias_coefs
# 1 1000001 500 1  21548.3      2.485176 998325.7 -4.155731e-07 9.378498e-06 5.652281e-05  0.9961065
#    R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true
# 1 0.9281978 -0.01760073    0.1311031      0.1884504     0.9961065    0.980418                50
#   support_size_estimated estimated_over_true_support_size TP FP     TN FN precision recall
# 1                     24                             0.48 24  0 999950 26         1   0.48
#   specificity accuracy  F1_score      gic     ebic      bic      aic     mbic
# 1           1 0.999974 0.6486486 -1996020 -1995961 -1996496 -1996601 -1995875


# L0glm fit using MCP estimates as starting values ####

# TODO test if orthogonal initialisation would work here too
# here L0glm is shown to be a very effective way to remove false positives from
# the MCP solution

# TODO check if MCP like penalties could be used within L0glm itself
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3212875/

startcoefs = beta_est0*scaling # MCP starting coefficients work best for high dimensional problems with very sparse solution
priormean = startcoefs
cutoff = quantile(abs(priormean)[priormean!=0], 0.5)
priormean[abs(priormean)<=cutoff] = 0
nonzero = startcoefs!=0
lambda_preset = preset.lambda(IC = "mbic", # options aic, bic, gic, ebic, mbic 
                              p = ncol(X), # note that we use original problem size
                              n = nrow(X))

system.time(L0glm_gaus_presetlambda <- L0glm.fit(
  X = X[, nonzero, drop=F], 
  y = y, 
  family = family, 
  lambda = lambda_preset, # mbic because high dimensional 
  nonnegative = FALSE,
  prior.mean = priormean[nonzero],
  start = startcoefs[nonzero], 
  warm.start = 1, 
  # TODO make sure that with warm.start=1 vars with zero coefs are efficiently removed from the start
  gamma = 2, 
  normalize = FALSE, # normalize = T doesn't work yet - too sparse then
  # dispersion = dispersion,
  solver = "eigen", # note: solver osqp better than eigen for very small lambda values))
  method = 3L, # faster than 7 here
  relaxed = FALSE)) # 1.8s with warm.start=1 on full matrix - 0.06s when subsetted - should ideally be the same!!!

L0glm_gaus_presetlambda$dispersion # 1.366601

L0glm_gaus_presetlambda_coefs = rep(0, ncol(X))
L0glm_gaus_presetlambda_coefs[nonzero] = L0glm_gaus_presetlambda$coefficients/scaling[nonzero]

beta_est1 = L0glm_gaus_presetlambda_coefs
sum(beta_est1!=0&beta_true!=0)  # 33
table(beta_est1!=0, beta_true!=0, 
      dnn=c("estimated beta nonzero", "true beta nonzero"))
#                          true beta nonzero
# estimated beta nonzero   FALSE   TRUE
#                  FALSE  999950     17
#                   TRUE       0     33

cor(beta_est1, beta_true) # 0.9919251

plot(x=beta_est1[beta_true!=0&beta_est1!=0], 
     y=beta_true[beta_true!=0&beta_est1!=0], pch=16, col="blue",
     xlab="Estimated coefficients", ylab="True coefficients",
     main="L0glm after MCP") # true positives
abline(a=0, b=1, col="lightgrey")
points(x=beta_est1[beta_true!=0&beta_est1==0], 
       y=beta_true[beta_true!=0&beta_est1==0], 
       pch=16, col="purple") # false negatives
points(x=beta_est1[beta_true==0&beta_est1!=0], 
       y=beta_true[beta_true==0&beta_est1!=0], 
       pch=16, col="red") # false positives !!!!!
points(x=beta_est1[beta_true==0&beta_est1==0][[1]], 
       y=beta_true[beta_true==0&beta_est1==0][[1]], 
       pch=16, col="green") # true negatives

# check metrics & quality of solution with preset lambda
system.time(metrics_L0 <- data.frame(do.call(cbind, get.metrics(beta = beta_est1,  # some metrics
                                                                  X = X, 
                                                                  y = y,
                                                                  family = gaussian(identity), 
                                                                  beta_true = beta_true)))) 
metrics_L0
#         p   n m deviance deviance_true   loglik      MD_coefs    MAD_coefs  sMAPE_coefs bias_coefs
# 1 1000001 500 1 21441.86     0.5321476 998317.9 -2.467697e-08 4.243811e-06 3.689456e-05  0.9917895
#    R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true
# 1 0.9839154 0.007209908   0.05737225     0.08771407     0.9917895   0.9955724                50
#   support_size_estimated estimated_over_true_support_size TP FP     TN FN precision recall
# 1                     33                             0.66 33  0 999950 17         1   0.66
#   specificity accuracy  F1_score      gic     ebic      bic      aic     mbic
# 1           1 0.999983 0.7951807 -1995778 -1995696 -1996425 -1996568 -1995579



# ADAPTIVE L0GLM FIT ####

# For the high dimensional case a better solution is obtained by using the 1st iteration
# solution and then using the L0glm output coefficients as starting values in 
# subsequent iterations. This would be sort of like an "adaptive L0 glm" then, where
# the original starting values that are used would themselves also be optimised.
# Given that we start from sensible starting values we can specify a prior mean here.
startcoefs = beta_est1*scaling

for (iter in (1:5)) {
  priormean = startcoefs
  # cutoff = quantile(abs(priormean)[priormean!=0], 0.1) # find a better smooth functional form?
  # priormean[abs(priormean)<=cutoff] = 0

sel = startcoefs!=0

# we use mbic because original problem is high dimensional
lambda_preset = preset.lambda(IC = "mbic", # options aic, bic, gic, ebic, mbic 
                              p = ncol(X), # note that we use original problem size
                              n = nrow(X))

system.time(L0glm.out_updated <- L0glm.fit(X = X[,sel,drop=F], 
                                    y = y, 
                                    family = gaussian(identity), 
                                    lambda = lambda_preset,  
                                    # no.pen = 1,
                                    nonnegative = FALSE,
                                    prior.mean = priormean[sel], # we can now specify our prior centers
                                    start = startcoefs[sel], # 
                                    warm.start = 0.95, # note 0.5-0.8 also work well
                                    # dispersion = L0glm_gaus_presetlambda$dispersion, # L0glm.out_2nditeration$dispersion,
                                    solver = "eigen", # eigen, glmnet or osqp, # note: constrained solvers glmnet & osqp better than eigen for very small lambda values
                                    method = 3L, # for dense case 3L 2x as fast than 7L
                                    gamma = 2, #
                                    thresh = 1E-5,
                                    relaxed = FALSE
))  # 
L0glm.out_updated_coefs = rep(0, ncol(X))
L0glm.out_updated_coefs[sel] = L0glm.out_updated$coefficients/scaling[sel]
print(sum(L0glm.out_updated_coefs!=0)) # 43
L0glm.out_updated$dispersion # 0.87
# TODO NOTE if support size is slightly misestimated dispersion can be massively different -
# this should not be the case - find bug??

# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# TODO for high dimensional case always do it like this - maybe even for low dimensional
# doing one extra 3d iteration still gets a better solution
# 
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

beta_est1B = L0glm.out_updated_coefs
sum(beta_est1B!=0&beta_true!=0)  # 33
table(beta_est1B!=0, beta_true!=0, 
      dnn=c("estimated beta nonzero", "true beta nonzero"))

cor(beta_est1B, beta_true) 

plot(x=beta_est1B[beta_true!=0&beta_est1B!=0], 
     y=beta_est1B[beta_true!=0&beta_est1B!=0], pch=16, col="blue",
     xlab="Estimated coefficients", ylab="True coefficients",
     main=paste0("adaptive L0glm after MCP iteration ", iter)) # true positives
abline(a=0, b=1, col="lightgrey")
points(x=beta_est1B[beta_true!=0&beta_est1B==0], 
       y=beta_true[beta_true!=0&beta_est1B==0], 
       pch=16, col="purple") # false negatives
points(x=beta_est1B[beta_true==0&beta_est1B!=0], 
       y=beta_true[beta_true==0&beta_est1B!=0], 
       pch=16, col="red") # false positives
points(x=beta_est1B[beta_true==0&beta_est1B==0][[1]], 
       y=beta_true[beta_true==0&beta_est1B==0][[1]], 
       pch=16, col="green") # true negatives

startcoefs = beta_est1B*scaling
}

sum(beta_est1B!=0) # 33
sum(beta_est1B!=0&beta_true!=0)  # 33
L0glm.out_updated$dispersion # 1.36228
table(beta_est1B!=0, beta_true!=0, 
      dnn=c("estimated beta nonzero", "true beta nonzero"))
#                         true beta nonzero
# estimated beta nonzero  FALSE   TRUE
# FALSE                  999950     17
# TRUE                        0     33

system.time(metrics_adL0 <- data.frame(do.call(cbind, get.metrics(beta = beta_est1B,  # some metrics
                                                                   X = X, 
                                                                   y = y,
                                                                   family = gaussian(identity), 
                                                                   beta_true = beta_true)))) 
metrics_adL0
#         p   n m deviance deviance_true loglik      MD_coefs    MAD_coefs  sMAPE_coefs bias_coefs
# 1 1000001 500 1 21435.76     0.4853219 998318 -1.443494e-07 4.003427e-06 3.616955e-05  0.9888823
#    R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true
# 1 0.9852959 0.003583471   0.05008787     0.06574387     0.9888823   0.9969693                50
#   support_size_estimated estimated_over_true_support_size TP FP     TN FN precision recall
# 1                     33                             0.66 33  0 999950 17         1   0.66
#   specificity accuracy  F1_score      gic     ebic      bic      aic     mbic
# 1           1 0.999983 0.7951807 -1995778 -1995697 -1996425 -1996568 -1995579


# COMPARISON OF METRICS ####

rbind(metrics_MCP,
      metrics_adMCP,
      metrics_L0,
      metrics_adL0)
#         p   n m deviance deviance_true   loglik      MD_coefs    MAD_coefs  sMAPE_coefs bias_coefs
# 1 1000001 500 1 21401.21     0.7318924 998294.4 -5.635672e-07 8.210249e-06 8.102734e-05  0.9685639
# 2 1000001 500 1 21548.30     2.4851756 998325.7 -4.155731e-07 9.378498e-06 5.652281e-05  0.9961065
# 3 1000001 500 1 21441.86     0.5321476 998317.9 -2.467697e-08 4.243811e-06 3.689456e-05  0.9917895
# 4 1000001 500 1 21435.76     0.4853219 998318.0 -1.443494e-07 4.003427e-06 3.616955e-05  0.9888823
#    R2_coefs  MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true
# 1 0.9648158 -0.010791569   0.10831365     0.86012876     0.9685639   0.9738870                50
# 2 0.9281978 -0.017600728   0.13110306     0.18845038     0.9961065   0.9804180                50
# 3 0.9839154  0.007209908   0.05737225     0.08771407     0.9917895   0.9955724                50
# 4 0.9852959  0.003583471   0.05008787     0.06574387     0.9888823   0.9969693                50
#   support_size_estimated estimated_over_true_support_size TP FP     TN FN precision recall
# 1                     57                             1.14 34 23 999927 16 0.5964912   0.68
# 2                     24                             0.48 24  0 999950 26 1.0000000   0.48
# 3                     33                             0.66 33  0 999950 17 1.0000000   0.66
# 4                     33                             0.66 33  0 999950 17 1.0000000   0.66
#   specificity accuracy  F1_score      gic     ebic      bic      aic     mbic
# 1    0.999977 0.999961 0.6355140 -1995125 -1994986 -1996228 -1996473 -1994787
# 2    1.000000 0.999974 0.6486486 -1996020 -1995961 -1996496 -1996601 -1995875
# 3    1.000000 0.999983 0.7951807 -1995778 -1995696 -1996425 -1996568 -1995579
# 4    1.000000 0.999983 0.7951807 -1995778 -1995697 -1996425 -1996568 -1995579

# CONCLUSION: deviance relative to zero noise ground truth as well as
# accuracy, F1 score, R2 between estimated & true coefficients plus MD,
# MAD & sMAPE for selected coefs all best for adaptive L0glm solution

# TODO : remove benchmarks related to bias on all coefficients, as measures of
# bias on selected variables are then confounded with the efficiency of variable
# selection / support recovery?

# orthogonal coefficient initialisation maybe also possible to use as start? ####
# this could be other option that could already be implemented in 1st iteration?



# TODO test initialisation with orthogonal solution
# TODO this is currently not optimised for p > n case...
# beta_init_orthogonal <- init_coefs(X = X,
#                                    y = y,
#                                    mu = y,
#                                    family = gaussian(identity),
#                                    lambda = lambda_preset)
# cor(beta_init_orthogonal, data_gaus$beta_true) # 0.981208 # 
# plot(beta_init_orthogonal, data_gaus$beta_true, pch=16, col="blue")
# mean(abs(beta_init_orthogonal-data_gaus$beta_true))


# coefficients & metrics over a range of lambda values
# calculate max lambda where all coefs would be penalized to 0
maxlam <- maxlambda(X=X, 
                    weights=L0glm.out_updated$weights, 
                    z=L0glm.out_updated$z,
                    # no.pen=1,
                    coefs=L0glm.out_updated$coefficients,
                    family=family)
maxlam # 
mincoef <- min(abs(L0glm.out_updated$coefficients))

lambda_preset = preset.lambda(IC = "bic", # options aic, bic, gic, ebic, mbic 
                              p = ncol(X),
                              n = nrow(X))
# lambdas <- 10^seq(log10(lambda_preset/10), log10(maxlam), length.out=10)
lambdas <- 10^seq(log10(lambda_preset)-1, log10(lambda_preset)+1, length.out=10)
coefmat <- Matrix(0, nrow=length(lambdas), ncol=ncol(X))
warmstart <- 0.1 # % by which L0glm_gaus_presetlambda$coefficients will be used as start
system.time(metrics <- data.frame(do.call(rbind, lapply(1:length(lambdas), function (i) { 
  lambda <- lambdas[[i]]
  # start <- (L0glm_gaus_presetlambda$coefficients)*warmstart+start_coefs(X, family)*(1-warmstart)
  # start[L0glm_gaus_presetlambda$coefficients==0] <- 0
  
  L0glm.out <- L0glm.fit(X = X, 
                         y = y, 
                         family = gaussian(identity), 
                         lambda = lambda, 
                         # no.pen = 1,
                         nonnegative = FALSE,
                         start = L0glm.out_updated$coefficients, #, start, # 
                         warm.start = warmstart,
                         # dispersion = dispersion,
                         solver = "eigen", # eigen, glmnet or osqp, # note: constrained solvers glmnet & osqp better than eigen for very small lambda values
                         method = 7L, # 3L here 2x as fast than 7L
                         # gamma = 2, 
                         thresh = 1E-5,
                         relaxed = FALSE 
  ) 
  coefmat[i,] <<- as.vector(L0glm.out$coefficients)/scaling
  return( unlist(c(lambda=lambda, get.metrics(beta = coefmat[i,], 
                                              X = X, y = y, 
                                              family = gaussian(identity), 
                                              beta_true = beta_true))) )
} )))) # 19s
matplot(x=log10(metrics$lambda), 
        y=as.matrix(coefmat), type="l", lty=1,
        xlab="Log10(lambda)", ylab="Coefficient estimates", main="Regularisation path L0 Gaussian GLM")
# maxlambda seems a factor of 2 wrong here on log10 scale

library(scales)
abline(h=beta_true, col=alpha(1:6,0.2))
abline(v=log10(maxlam), col="red") # max lambda where all coefs would be penalized to 0


# coefficient estimates at best lambda
best = which.min(metrics[,crit]) # gic, ebic or mbic all work here; bic gives 6 false positives, aic would produce 20 false positives
L0glm.out_best_coefs = coefmat[best,]
sum(L0glm.out_best_coefs!=0) # 200
bestlambda = metrics$lambda[best] # 
bestlambda

# CHECK what's wrong here - IC are probably wrong??


# metrics

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





# 2. HIGH DIMENSIONAL EXAMPLE IN ORDINIS PACKAGE ####

# In the example below MCP seems to give a better solution than L0 glm :-(

# 2.1 nonnegative MCP penalized regression calculated using ordinis ####

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
sum(coefs_ordinis!=0) # 50, closely matches expected value with n_nonzero=50, but not e.g. with n_nonzero=200
library(scales)
abline(h=as.vector(beta_true[beta_true!=0]), col=alpha(1:6,0.2)) # note colours don't match the ones of ordinis

plot(x=coefs_ordinis[!(coefs_ordinis==0&c(0,beta_true)==0)], 
     y=c(0,beta_true)[!(coefs_ordinis==0&c(0,beta_true)==0)], pch=16, col="blue")

TP = sum(coefs_ordinis!=0&c(0, beta_true)!=0) 
TP # 47 true positives
FP = sum(coefs_ordinis!=0&c(0, beta_true)==0) 
FP # 3 false positives
FN = sum(coefs_ordinis==0&c(0, beta_true)!=0) 
FN # 3 false negatives
FDR = FP / (FP + TP) 
FDR # 0.06 = 6%


# 2.2 L0glm fit ####

X = normalizeX(X) # we first L2 norm normalize the columns of X

scaling = attr(X, "scaling")
beta_true = beta_true * scaling
family = gaussian(identity)
crit = "mbic" # works best for high dimensional case

# using lambda expected to optimise given IC
lambda_preset = preset.lambda(IC = crit, # options aic, bic, gic, ebic, mbic 
                              p = ncol(X),
                              n = nrow(X))
lambda_preset # 31.07304


# glm fit
# system.time(glm_gaus <- glm.fit(y=y, x=X, intercept=FALSE, family=gaussian(identity))) 
# Error: cannot allocate vector of size XX Gb for large p


system.time(L0glm_gaus_presetlambda <- L0glm.fit(
  X = X, # Matrix(X, sparse=T), 
  y = y, 
  family = family, 
  lambda = "aic", 
  nonnegative = TRUE,
  lower = rep(0, ncol(X)),
  upper = rep(Inf, ncol(X)),
  start = NULL, 
  gamma = 1.1,
  # gamma = 1.1, # NOTE: gamma of 1.1 here works better than 2
  normalize = FALSE, # normalize = T doesn't work yet - too sparse then
  # dispersion = dispersion,
  solver = "eigen", # note: solver osqp better than eigen for very small lambda values))
  method = 7L,
  relaxed = FALSE,
  thresh = 1E-5))

sum(L0glm_gaus_presetlambda$coefficients!=0)  # 197
L0glm_gaus_presetlambda$dispersion # 5.3
plot(L0glm_gaus_presetlambda$coefficients[!(L0glm_gaus_presetlambda$coefficients==0|beta_true==0)], 
     beta_true[!(L0glm_gaus_presetlambda$coefficients==0|beta_true==0)], pch=16, col="blue")
cor(L0glm_gaus_presetlambda$coefficients, beta_true) # 0.9999994
# plot(x=L0glm_gaus_presetlambda$coefficients, y=beta_true, pch=16, col="blue")
# points(x=L0glm_gaus_presetlambda$coefficients[data_gaus$beta_true==0], y=data_gaus$beta_true[data_gaus$beta_true==0], 
#       pch=16, col="red")

# check metrics & quality of solution with preset lambda
system.time(metrics_gaus <- data.frame(do.call(cbind, get.metrics(beta = L0glm_gaus_presetlambda$coefficients,  # some metrics
                                                                  X = X, 
                                                                  y = y,
                                                                  family = gaussian(identity), 
                                                                  beta_true = beta_true)))) 
metrics_gaus


# for high dimensional case a better solution is obtained by using the 1st iteration
# solution as starting values
L0glm.out_2nditeration <- L0glm.fit(X = X, 
                                    y = y, 
                                    family = gaussian(identity), 
                                    lambda = "mbic", 
                                    # no.pen = 1,
                                    gamma = 1.1,
                                    nonnegative = TRUE,
                                    lower = rep(0, ncol(X)),
                                    start = L0glm_gaus_presetlambda$coefficients, #, start, # 
                                    warm.start = 0.1,
                                    dispersion = L0glm_gaus_presetlambda$dispersion,
                                    solver = "eigen", # eigen, glmnet or osqp, # note: constrained solvers glmnet & osqp better than eigen for very small lambda values
                                    method = 7L, # 3L here 2x as fast than 7L
                                    # gamma = 2, 
                                    thresh = 1E-5,
                                    relaxed = FALSE 
) 
L0glm.out_2nditeration_coefs = L0glm.out_2nditeration$coefficients
sum(L0glm.out_2nditeration_coefs!=0) # 226
L0glm.out_2nditeration$dispersion # 12.2
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# DIDN't GET SUCH A GOOD SOLUTION FOR THIS PROBLEM - HERE MCP SEEMS SUPERIOR
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

plot(x=L0glm.out_2nditeration_coefs[!(L0glm.out_2nditeration_coefs==0&beta_true==0)], 
     y=beta_true[!(L0glm.out_2nditeration_coefs==0&beta_true==0)], pch=16, col="blue",
     xlab="estimated coefficients", ylab="true coefficients")

system.time(metrics_gaus2 <- data.frame(do.call(cbind, get.metrics(beta = L0glm.out_2nditeration_coefs,  # some metrics
                                                                   X = X, 
                                                                   y = y,
                                                                   family = gaussian(identity), 
                                                                   beta_true = beta_true)))) 
metrics_gaus2




# 2.3 L0glm fit on MCP pre-selected variables ####

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

selected = which(coefs_ordinis[-1]!=0)

# using lambda expected to optimise given IC
lambda_preset = preset.lambda(IC = crit, # options aic, bic, gic, ebic, mbic 
                              p = ncol(X),
                              n = nrow(X))
lambda_preset # 6.214608

warm_start =  0.9 # 0.9 # 1 would essentially use MCP as pre-screen
# 0.9-1 all work well, but 1 should be faster as it should kick out 0 coefficients from the start
# TODO check in code that it is doing this in an optimal way, from the very start

# TODO: add warm_start as option in L0glm.fit & use this like below if start!=NULL
gamma = 2
start = coefs_ordinis[-1]*scaling # TODO incorporate scaling if normalize=T in L0glm.fit.internal
# start = warm_start*start + (1-warm_start)*rep(sd(start)^(1/gamma), ncol(X))

# A warm_start value of 1 would indicate using the provided starting coefficients 
# directly (full warm start), while a warm_start between 0 and 1 would represent 
# a "cooled" or "tempered" start, blending the provided coefficients with  
# their standard deviations raised to the power of 1/gamma. In this way, one can
# specify one's prior belief in the starting coefficients being correct vs.
# implying identical priors on al coefficients (but still allowing 
# the observed variability in the starting coefficients to affect the chosen
# constant SD on the implied Gaussian prior of all coefficients).
system.time(L0glm.out_preset <- L0glm.fit(Xscaled[,selected,drop=F], 
                                          y, 
                                          family = family, 
                                          lambda = 1E-10, # *lambda_preset, 
                                          # prior.mean = start[selected], # use nonzero centered ridge regression?
                                          nonnegative = TRUE,
                                          lower = rep(0, ncol(Xscaled[,selected,drop=F])),
                                          start = start[selected], # rep(sqrt(50), ncol(X)), # NULL, # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                                          warm.start = 0.1,
                                          gamma = 2,
                                          solver = "glmnet", # note: solver osqp better than eigen for very small lambda values))
                                          method = 7L,
                                          verbose = T,
                                          thresh = 1E-30, # had to set this to low value here # TODO thresh should also be multiplied by scaling if normalize=T
                                          maxit = 100,
                                          tol = 1E-30
                                          # control.solver = modifyList(do.call("control.solver.gen", 
                                          #                                     list(solver = "eigen", 
                                          #                                          method = 4L)), 
                                          #                             list(rowaug=TRUE))
)) 
# using only preselected variables: 0.11s

L0glm.out_preset_coefs = rep(0, ncol(X))
L0glm.out_preset_coefs[selected] = L0glm.out_preset$coefficients/scaling[selected] # TODO do this scaling automatically if normalize=T
sum(L0glm.out_preset_coefs!=0) # 50
plot(x=L0glm.out_preset_coefs[!(L0glm.out_preset_coefs==0&beta_true==0)], 
     y=beta_true[!(L0glm.out_preset_coefs==0&beta_true==0)], pch=16, col="blue")
# NOTE: a BIC pre-picked lambda here worked well since the MCP
# pre-screening already changed this from a high dimensional to a low dimensional
# problem
# Working on MCP or SCAD preselected coefficients works well here, combined with warm_start>=0.9
# LASSO did not work here, as that gave a bad solution in ordinis
# ALASSO also did not work here
# Also note that with warm_start=0 a constrained solver like glmnet here performs
# better than the unconstrained Eigen solvers. With warm_start=0.9 performance of
# both is near-identical though.

# coefficients & metrics over a range of lambda values
# calculate max lambda where all coefs would be penalized to 0
maxlam <- maxlambda(X=Xscaled, weights=L0glm.out_preset$weights, 
                    z=y, no.pen=NULL, coefs=L0glm.out_preset$coefficients,
                    family=family)
maxlam #

lambdas <- 10^seq(round(log10(lambda_preset))-1, log10(maxlam), by=0.2)
# lambdas <- 10^seq(round(log10(lambda_preset))-2, log10(lambda_preset)+1, by=0.1)
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
  L0glm.out <- L0glm.fit(Xscaled[,selected,drop=F], y, 
                         family = gaussian(identity), 
                         lambda = lambda, 
                         # prior.mean = FALSE, # use nonzero centered ridge regression?
                         # ridge.center = NULL,
                         nonnegative = TRUE,
                         lower = rep(0, ncol(Xscaled[,selected,drop=F])),
                         no.pen = NULL,
                         start = start[selected], # NULL, # rep(0.1, ncol(sim$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                         solver = "eigen", # note: solver osqp better than eigen for very small lambda values))
                         method = 5L,
                         relaxed = FALSE,
                         thresh = 1E-20, # had to set this to low value here # TODO thresh should also be multiplied by scaling if normalize=T
                         maxit = 20,
                         tol = 1E-5,
                         gamma = gamma
  )$coefficients/scaling[selected]
  coefmat[i,selected] <<- as.vector(L0glm.out)
  return( unlist(c(lambda=lambda, get.metrics(beta = coefmat[i,selected], 
                                              X = X[,selected,drop=F], y = y, 
                                              family = gaussian(identity), 
                                              beta_true = beta_true[selected],
                                              p = ncol(X)))) )
} )))) 
# 3.6s
matplot(x=log10(metrics$lambda), y=as.matrix(coefmat[,colSums(coefmat)!=0,drop=F]), 
        type="l", lty=1,
        xlab="Log10(lambda)", ylab="Coefficient estimates", main="Regularisation path L0 Gaussian GLM")
abline(v=log10(maxlam), col="red") # max lambda where all coefs would be penalized to 0
library(scales)
abline(h=as.vector(beta_true[beta_true!=0]), col=alpha(1:6,0.2))


# coefficients at best lambda, here the one optimising the chosen IC
# best = which.min(abs(metrics$estimated_over_true_support_size-1))
best = which.min(metrics[,"mbic"]) # here using mBIC worked best
L0glm.out_best_coefs = coefmat[best,]
sum(L0glm.out_best_coefs!=0) # 50
TP = sum(L0glm.out_best_coefs!=0&beta_true!=0) 
TP # 47 true positives
FP = sum(L0glm.out_best_coefs!=0&beta_true==0) 
FP # 3 false positive
FN = sum(L0glm.out_best_coefs==0&beta_true!=0) 
FN # 3 false negatives
FDR = FP / (FP + TP) 
FDR # 0.06 = approx 6%

# the familywise error rate I presume would be X/1000000
# NOTE: the mBIC criterion is supposed to control the FWER rate
# while the mBIC2 criterion is supposed to control the FDR rate
# https://www.sciencedirect.com/science/article/pii/S016794731100171X

plot(x=L0glm.out_best_coefs[!(L0glm.out_best_coefs==0&beta_true==0)], 
     y=beta_true[!(L0glm.out_best_coefs==0&beta_true==0)], pch=16, col="blue")

# NOTE: for high dimensional problems with an extremely sparse solution
# pre-picking variables based on an MPC & SCAD penalized regression work well
# in this example. This is because of the more aggressive shrinking towards zero
# from the very beginning, which works better for really sparse problem.



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
plot(x=L0glm.out_preset_coefs[!(L0glm.out_preset_coefs==0|beta_true_scaled==0)], 
     y=as.vector(beta_true_scaled[!(L0glm.out_preset_coefs==0|beta_true_scaled==0)]),
     pch=16, col="blue")

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

plot(x=L0glm.out_best_coefs[!(L0glm.out_best_coefs==0|beta_true_scaled==0)], 
     y=as.vector(beta_true_scaled[!(L0glm.out_best_coefs==0|beta_true_scaled==0)]),
     pch=16, col="blue")


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

plot(x=coefs_ordinis[-1][!(coefs_ordinis[-1]==0|beta_true_scaled==0)], 
     y=as.vector(beta_true_scaled[!(coefs_ordinis[-1]==0|beta_true_scaled==0)]),
     pch=16, col="blue")



# 3. TEST HIGH DIMENSIONAL SUPERRESOLUTION DECONVOLUTION ####

## 3.1 Single channel deconvolution assuming Gaussian noise ####
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
                              p = ncol(X), 
                              n = nrow(X)) 
lambda_preset # lambda that would be expected to optimise given IC

L0glm.out_preset <- L0glm.fit(X, y, 
                              family = family, 
                              lambda = lambda_preset, 
                              nonnegative = TRUE,
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
                         nonnegative = TRUE,
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








