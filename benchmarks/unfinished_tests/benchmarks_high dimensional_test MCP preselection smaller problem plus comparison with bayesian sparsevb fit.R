# Test high dimensional variable selection - here I use MCP as a pre-screen on a large problem and then fit an L0glm on the pre-selected subset to remove some false positives & also test an adaptive L0glm algorithm ####

# also includes comparison with sparsevb results (mean-field posterior approximations 
# for both linear and logistic regression models, including variable selection 
# via sparsity-inducing spike and slab priors

# fix to set working package repository when using Microsoft R Open with Intel MKL
options (repos = c(CRAN="https://cloud.r-project.org/"))
library(remotes)
# remotes::install_github("tomwenseleers/L0glm/L0glm")
# library(L0glm)
devtools::load_all(".")
library(export)
# source("..//benchmarks/generate data.R")
library(Matrix)


# 1. TEST HIGH DIMENSIONAL PROBLEM WITH STANDARD NORMAL COVARIATES ####

n = 500 # 500 # 1000
p = 100000 # 100000 # 10000 # 1000-4000 works, between 5000 and 10000 solution breaks down - maybe pre-selecting variables with MCP might still work then
s = 0.001 # 0.05, sparsity of true solution - better solution the sparser the true solution is
n_nonzero = round(p*s)
n_nonzero # 100
beta_sd = 1 # this would directly correspond to Cohen's d - 0.2, 0.5 or 0.8 would be considered small, medium or large
error_sd = 1 
dispersion = error_sd^2

data_gaus  = simulate_standardnormal(n=n, p=p, n_nonzero=n_nonzero, 
                                 beta_sd=beta_sd, 
                                 error_sd=error_sd,
                                 seed.problem=1,
                                 seed.noise=1) 
X = normalizeX(data_gaus$X) # we first L2 norm normalize the columns of X
Xsp = Matrix(X, sparse=T)
y = data_gaus$y
# X = data_gaus$X

# X = scale(data_gaus$X, center=T, scale=F)
# y = scale(data_gaus$y, center=T, scale=F)

scaling = attr(X, "scaling")
# beta_true = data_gaus$beta_true * scaling
beta_true = data_gaus$beta_true
beta_true_set = which(beta_true!=0)
L1norm_betatrue = sum(abs(beta_true))
L0norm_betatrue = sum(beta_true!=0)
L0overL1norm_betatrue = L0norm_betatrue/L1norm_betatrue
L0overL1norm_betatrue # 1.44
family = gaussian(identity)


# 1.1 LASSO fit using glmnet or biglasso TO PRE-SCREEN VARIABLES - WORSE THAN MCP SO NOT USED ####

library(glmnet)
system.time(glmnet_cvfit <- cv.glmnet(x = X, y = y, nfolds = 5,
                                      intercept=F, alpha=0.9)) # 0.2s
plot(glmnet_cvfit)
bestlambda1 <- glmnet_cvfit$lambda.min #
bestlambda1 # 0.07
bestlambda2 <- glmnet_cvfit$lambda.1se #
bestlambda2 # 0.13
system.time(glmnet_fit <- glmnet(x = X, y = y, intercept = F, alpha = 0.9)) # 28s
plot(glmnet_fit)
coefs_glmnet1 <- as.vector(coef(glmnet_fit, s=bestlambda1)[-1])
sum(coefs_glmnet1!=0) # 198
glmnet_set1 <- which(coefs_glmnet1!=0)
length(intersect(glmnet_set1, beta_true_set)) # 32
table(coefs_glmnet1!=0, beta_true!=0,  
      dnn=c("estimated beta nonzero", "true beta nonzero"))
#                        true beta nonzero
# estimated beta nonzero FALSE TRUE
#                  FALSE 99734    68
#                   TRUE   166    32

beta_LASSO = coefs_glmnet1/scaling
plot(x=beta_LASSO[beta_true!=0&beta_LASSO!=0], 
     y=beta_true[beta_true!=0&beta_LASSO!=0], pch=16, col="blue",
     xlab="Estimated coefficients", ylab="True coefficients",
     main=paste0("LASSO fit (glmnet)")) # true positives
abline(a=0, b=1, col="lightgrey")
points(x=beta_LASSO[beta_true!=0&beta_LASSO==0], 
       y=beta_true[beta_true!=0&beta_LASSO==0], 
       pch=16, col="purple") # false negatives
points(x=beta_LASSO[beta_true==0&beta_LASSO!=0], 
       y=beta_true[beta_true==0&beta_LASSO!=0], 
       pch=16, col="red") # false positives
points(x=beta_LASSO[beta_true==0&beta_LASSO==0][[1]], 
       y=beta_true[beta_true==0&beta_LASSO==0][[1]], 
       pch=16, col="green") # true negatives
# graph2png(file="..//benchmarks/plots/adaptive_L0glm/0_LASSO_p_1_million_n_500_50 nonzero.png", 
#          width=7, height=5)

# check metrics & quality of LASSO solution
system.time(metrics_LASSO <- data.frame(do.call(cbind, get.metrics(beta = coefs_glmnet1/scaling,  # some metrics
                                                                   X = X, 
                                                                   y = y,
                                                                   family = gaussian(identity), 
                                                                   beta_true = beta_true)))) 
metrics_LASSO
#        p   n m deviance deviance_true   loglik     MD_coefs   MAD_coefs sMAPE_coefs bias_coefs
# 1 100001 500 1  54735.1      22.93101 97918.62 7.143014e-05 0.000768463 0.005035101    1.57039
#    R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true
# 1 0.5634683  0.01086158    0.2213963       1.856112       1.57039   0.7255618               100
#   support_size_estimated estimated_over_true_support_size TP  FP    TN FN precision recall
# 1                    198                             1.98 32 166 99734 68 0.1616162   0.32
#   specificity accuracy  F1_score       gic      ebic       bic       aic      mbic
# 1   0.9983383  0.99766 0.2147651 -191651.7 -191255.1 -194600.5 -195439.2 -190570.1


# note: the cross validation I might not really need, as I can just select a model
# with max n variables
# 
# # note: adaptive LASSO fit also doesn't help to select more variables
# 
# # coefs_glmnet2 <- as.vector(coef(glmnet_fit, s=bestlambda2)[-1])
# # sum(coefs_glmnet2!=0) # 55
# # glmnet_set2 <- which(coefs_glmnet2!=0)
# # length(intersect(glmnet_set2, beta_true_set)) # 20
# 
# # library(biglasso)
# # system.time(Xbig <- as.big.matrix(X)) # 16s
# # system.time(biglasso_cvfit <- cv.biglasso(X = Xbig, y = y, nfolds = 5,
# #                                           alpha = 0.9)) # 44s
# # plot(biglasso_cvfit)
# # bestlambda1 <- biglasso_cvfit$lambda.min # 
# # bestlambda1 # 0.3278291
# # bestlambda2 <- biglasso_cvfit$lambda.1se # 
# # bestlambda2 # 0.6379187 # can maybe be roughly estimated from ratio of L1 and L0 norm & estimated dispersion/variance??
# # system.time(biglasso_fit <- biglasso(X = Xbig, y = y, alpha = 0.9)) # 14s
# # plot(biglasso_fit)
# # coefs_biglasso1 <- as.vector(coef(biglasso_fit, lambda=bestlambda1)[-1])
# # sum(coefs_biglasso1!=0) # 180
# # biglasso_set1 <- which(coefs_biglasso1!=0)
# # length(intersect(biglasso_set1, beta_true_set)) # 27
# # 
# # # coefs_biglasso2 <- as.vector(coef(biglasso_fit, lambda=bestlambda2)[-1])
# # # sum(coefs_biglasso2!=0) # 48
# # # biglasso_set2 <- which(coefs_biglasso2!=0)
# # # length(intersect(biglasso_set2, beta_true_set)) # 20


# 1.2 MCP fit using ordinis ####

library(ordinis)
gamma_mcp = 1.1 # default for MCP is 1.4, I am using value close to 1 to match best with my later L0 glm objective
system.time(ordinis_fit <- ordinis(X, y, intercept=F,
                                   penalty = "mcp",
                                   gamma = gamma_mcp, 
                                   # lower.limits = rep(0, p), # force all coefficients to be positive
                                   # penalty.factor = c(0, 0, rep(1, p-2)), # don't penalize first two coefficients
                                   nlambda = 100L, # for pre-screening I just need 1 solution, so nlambda=1 OK (and works best in fact)
                                   # lambda.min.ratio = 1E-10,
                                   dfmax = n/log(n),
                                   alpha = 1, #,
                                   maxit = 5 #, # NOTE: can be set to just 1 iteration when used for pre-screening
                                   # tol = 0.1
)) # 0.01s with maxit 1 (with nlambda 5), 64s with defaults (with nlambda 100) 
# note: later L0glm solution initialised with MCP solution is about as good
# with maxit=1 and tol = 0.1 & that's twice as fast
# plot(ordinis_fit)
ordinis_fit$lambda[which.max(logLik(ordinis_fit))] # best lambda = 0.07
min(AIC(ordinis_fit)) #
min(BIC(ordinis_fit)) #
which.min(BIC(ordinis_fit)) # NOTE: this is the smallest lambda value when nlambda=100
# coefs_ordinis = ordinis_fit$beta[,which.min(BIC(ordinis_fit))][-1]
# coefs_ordinis = ordinis_fit$beta[,which.min(AIC(ordinis_fit))][-1]
coefs_ordinis = ordinis_fit$beta[,which.max(logLik(ordinis_fit))][-1] # gives same result as AIC & BIC here
sum(coefs_ordinis!=0) # 79
sum(coefs_ordinis!=0&beta_true!=0) # 68 true pos
sum(coefs_ordinis!=0&beta_true==0) # 11 false positives still
L1norm_coefsordinis = sum(abs(coefs_ordinis))
L0norm_coefsordinis = sum(coefs_ordinis!=0)
L0overL1norm_coefsordinis = L0norm_betatrue/L1norm_betatrue
L0overL1norm_coefsordinis # 1.23

beta_est0 = coefs_ordinis/scaling
range(beta_est0) # 
range(beta_true) # 

table(beta_est0!=0, beta_true!=0, 
      dnn=c("estimated beta nonzero", "true beta nonzero"))
#                            true beta nonzero
# estimated beta nonzero         FALSE   TRUE
#                         FALSE  99889    32
#                         TRUE      11    68

cor(beta_est0, beta_true)  # 0.966

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
# graph2png(file="..//benchmarks/plots/adaptive_L0glm/1_MCP_regression_p_1_million_n_500_50 nonzero.png", 
#          width=7, height=5)
# (NOTE: graph looks very similar if higher nlambda & higher nr of maxit & lower tol is used)

# check metrics & quality of MCP solution
system.time(metrics_MCP <- data.frame(do.call(cbind, get.metrics(beta = beta_est0,  # some metrics
                                                                 X = X, 
                                                                 y = y,
                                                                 family = gaussian(identity), 
                                                                 beta_true = beta_true)))) 
metrics_MCP
#        p   n m deviance deviance_true   loglik      MD_coefs    MAD_coefs  sMAPE_coefs bias_coefs
# 1 100001 500 1 52793.99       1.29848 98046.64 -6.172947e-06 0.0001512331 0.0009241663  0.9901606
#    R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true
# 1 0.9666396 0.006914065    0.1072736      0.3597041     0.9901606   0.9847443               100
#   support_size_estimated estimated_over_true_support_size TP FP    TN FN precision recall
# 1                     79                             0.79 68 11 99889 32 0.8607595   0.68
#   specificity accuracy  F1_score       gic      ebic       bic       aic      mbic
# 1   0.9998899  0.99957 0.7597765 -194410.6 -194251.2 -195596.1 -195933.3 -193975.9




# 1.3 iterative adaptive MCP fit ####

# attempt to use iterative adaptive MCP regression to approximate
# L0 penalized regression & penalize the small coefficients to zero

startcoefs = beta_est0*scaling
maxcoef = max(abs(beta_est0[beta_est0!=0]))
delta = maxcoef/10 # use in penfact, decreased in later iterations
percentile = 0.5

plotit = TRUE

system.time(
  {
    for (iter in (1:30)) {
      penfact = 1/(abs(startcoefs)^2+delta)  # adaptive penalty based on initial MCP estimate
      # adaptive penalty based on initial MCP solution
      # normally I had expected 1/(abs(startcoefs)+delta) to work to approximate the L0 norm,
      # as that is the penalty one would use with LASSO to approximate the L0 penalty
      # (with ridge it would be 1/(abs(coefs)^2+delta^2)
      # but I had to use the ^2 to make convergence faster
      # I also let delta decrease in each iteration
      
      penfact = length(penfact)*penfact/sum(penfact)
      
      # given that we have good prior estimates for the coefficients we can 
      # shrink towards a given prior mean here instead of shrinking to zero
      # except for the small coefficients (10% lowest) - for those ones we will
      # still shrink to zero
      priormean = startcoefs
      cutoff = quantile(abs(priormean)[priormean!=0], percentile) # find a better smooth functional form instead of a hard cutoff?
      priormean[abs(priormean)<=cutoff] = 0 # still using zero prior center for small coefficients
      
      y_star = y - (X %*% priormean) 
      # cf https://stats.stackexchange.com/questions/76925/ridge-regression-in-r-where-coefficients-are-penalized-toward-numbers-other-than
      nonzero = startcoefs!=0
      
      system.time(ordinis_adMCP_fit <- ordinis(X[,nonzero,drop=F], 
                                               y_star, intercept=F,
                                               penalty = "mcp",
                                               gamma = gamma_mcp, # default for MCP=1.4
                                               # lower.limits = rep(0, p), # force all coefficients to be positive
                                               penalty.factor = penfact[nonzero],
                                               lambda = preset.lambda(IC="mbic", p=p, n=n)*L0overL1norm_coefsordinis, # my attempt to convert the LASSO like MCP lambda scale to the lambda scale of L0 penalized regression & then apply a penalty equivalent to that implied by MBIC
                                               nlambda = 100L,
                                               # lambda.min.ratio = 1E-100,
                                               dfmax = n/log(n),
                                               alpha = 1# , 
                                               # maxit = 5 #, # could be set to 1 iteration when used for pre-screening
                                               # tol = 0.1
      )) # 
      beta_star = rep(0, ncol(X))
      beta_star[nonzero] = ordinis_adMCP_fit$beta[,which.max(logLik(ordinis_adMCP_fit))][-1]
      # beta_star = ordinis_adMCP_fit$beta[,which.max(logLik(ordinis_adMCP_fit))][-1]
      sum(ordinis_adMCP_fit$beta[,which.max(logLik(ordinis_adMCP_fit))][-1]!=0) # 
      coefs_ordinis_adMCP = beta_star+priormean
      print(sum(coefs_ordinis_adMCP!=0)) # 
      print(sum(coefs_ordinis_adMCP!=0&beta_true!=0)) # true pos
      sum(coefs_ordinis_adMCP!=0&beta_true==0) # false positives
      
      beta_est0B = coefs_ordinis_adMCP/scaling
      
      table(beta_est0B!=0, beta_true!=0, 
            dnn=c("estimated beta nonzero", "true beta nonzero"))
      
      if (plotit) {
        plot(x=beta_est0B[beta_true!=0&beta_est0B!=0], 
             y=beta_true[beta_true!=0&beta_est0B!=0], pch=16, col="blue",
             xlab="Estimated coefficients", ylab="True coefficients",
             main=paste0("iterative adaptive MCP regression iteration ", iter)) # true positives
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
      }
      
      startcoefs = coefs_ordinis_adMCP
      delta = delta*0.9
      
    }
  }) # 13s without plotting

cor(beta_est0B, beta_true) # 0.9987002

plot(x=beta_est0B[beta_true!=0&beta_est0B!=0], 
     y=beta_true[beta_true!=0&beta_est0B!=0], pch=16, col="blue",
     xlab="Estimated coefficients", ylab="True coefficients",
     main="iterative adaptive MCP regression (ordinis)") # true positives
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
# graph2png(file="..//benchmarks/plots/adaptive_L0glm/2_adaptive_MCP_regression_p_1_million_n_500_50 nonzero.png", 
#          width=7, height=5)

# metr2 = do.call(rbind, lapply(2:ncol(ordinis_fit2$beta), function (col) data.frame(col=col, do.call(cbind, get.metrics(beta = ordinis_fit2$beta[,col][-1], 
#                                                       X = X, 
#                                                       y = y,
#                                                       family = gaussian(identity), 
#                                                       beta_true = beta_true))) ))
# plot(x=log(ordinis_fit2$lambda[metr2$col]), y=metr2$deviance_true,
#      xlab="log(lambda)", ylab="deviance relative to ground truth", type="l")
# plot(x=log(ordinis_fit2$lambda[metr2$col]), y=metr2$deviance,
#      xlab="log(lambda)", ylab="deviance relative to ground truth", type="l")
# plot(x=log(ordinis_fit2$lambda[metr2$col]), y=metr2$recall,
#      xlab="log(lambda)", ylab="precision", type="l")
# plot(x=log(ordinis_fit2$lambda[metr2$col]), y=metr2$estimated_over_true_support_size,
#      xlab="log(lambda)", ylab="deviance relative to ground truth", type="l")
# plot(x=log(ordinis_fit2$lambda[metr2$col]), y=metr2$support_size_estimated,
#      xlab="log(lambda)", ylab="estimated support size", type="l")
# the lambda just before the flattish bit is the right lambda
# around exp(11) = 60000
# (reaches max sensitivity there, at lower lambdas only false positives are added)


# check metrics & quality of adaptive MCP solution
system.time(metrics_adMCP <- data.frame(do.call(cbind, get.metrics(beta = beta_est0B,  # some metrics
                                                                   X = X, 
                                                                   y = y,
                                                                   family = gaussian(identity), 
                                                                   beta_true = beta_true)))) 
metrics_adMCP
#        p   n m deviance deviance_true   loglik     MD_coefs    MAD_coefs  sMAPE_coefs bias_coefs
# 1 100001 500 1 52933.67      3.101343 98058.98 2.389182e-05 0.0001896103 0.0008807537   1.007464
#    R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true
# 1 0.9514736  0.03668904    0.1574621      0.2435662      1.007464   0.9778968               100
#   support_size_estimated estimated_over_true_support_size TP FP    TN FN precision recall
# 1                     66                             0.66 64  2 99898 36  0.969697   0.64
#   specificity accuracy  F1_score       gic      ebic       bic     aic      mbic
# 1     0.99998  0.99962 0.7710843 -194708.8 -194575.2 -195701.6 -195984 -194344.6



# 1.4 L0glm fit on MCP pre-selected variables ####

# TODO check if MCP-like penalties could be used within L0glm itself
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3212875/

startcoefs = beta_est0*scaling # MCP starting coefficients work best
# startcoefs = coefs_glmnet1*scaling 
# NOTE: LASSO glmnet coefficient starting values also work & 
# results in smaller nr of true coefficients being picked in adaptive L0glm step (27 instead of 33), but requires more iterations than with MCP starting values) 

priormean = startcoefs
cutoff = quantile(abs(priormean)[priormean!=0], 0.5) 
priormean[abs(priormean)<=cutoff] = 0 # we shrink small coefficients to zero & keep larger ones at their original value

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
  relaxed = FALSE)) # 0.14s

L0glm_gaus_presetlambda$dispersion # 0.79

L0glm_gaus_presetlambda_coefs = rep(0, ncol(X))
L0glm_gaus_presetlambda_coefs[nonzero] = L0glm_gaus_presetlambda$coefficients/scaling[nonzero]

beta_est1 = L0glm_gaus_presetlambda_coefs
sum(beta_est1!=0&beta_true!=0)  # 10
table(beta_est1!=0, beta_true!=0, 
      dnn=c("estimated beta nonzero", "true beta nonzero"))
#                    true beta nonzero
# estimated beta nonzero   FALSE   TRUE
#                  FALSE   99900    37
#                   TRUE       0    63

cor(beta_est1, beta_true) # 0.9835172

plot(x=beta_est1[beta_true!=0&beta_est1!=0], 
     y=beta_true[beta_true!=0&beta_est1!=0], pch=16, col="blue",
     xlab="Estimated coefficients", ylab="True coefficients",
     main="L0glm on MCP preselected variables") # true positives
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
# graph2png(file="..//benchmarks/plots/adaptive_L0glm/3_L0glm_on_MCP_presel_vars_p_1_million_n_500_50 nonzero.png", 
#            width=7, height=5)

# check metrics & quality of solution with preset lambda
system.time(metrics_L0 <- data.frame(do.call(cbind, get.metrics(beta = beta_est1,  # some metrics
                                                                X = X, 
                                                                y = y,
                                                                family = gaussian(identity), 
                                                                beta_true = beta_true))))  # 0.5s
metrics_L0
#        p   n m deviance deviance_true   loglik     MD_coefs    MAD_coefs  sMAPE_coefs bias_coefs
# 1 100001 500 1 52905.13      2.281651 98062.12 1.257498e-05 0.0001399101 0.0008176071   1.004675
#    R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true
# 1 0.9685761  0.03296894   0.08887472      0.1231859      1.004675   0.9923362               100
#   support_size_estimated estimated_over_true_support_size TP FP    TN FN precision recall
# 1                     63                             0.63 63  0 99900 37         1   0.63
#   specificity accuracy  F1_score       gic      ebic       bic       aic      mbic
# 1           1  0.99963 0.7730061 -194778.1 -194650.6 -195726.5 -195996.2 -194430.3


# ATTEMPT AT FULL L0glm FIT ON ORIGINAL PROBLEM SIZE USING MCP SELECTED VARIABLES AS STARTING VALUES #### 

startcoefs = beta_est0*scaling # MCP starting coefficients work best
# startcoefs = coefs_glmnet1*scaling 
# NOTE: LASSO glmnet coefficient starting values also work & 
# results in smaller nr of true coefficients being picked in adaptive L0glm step (27 instead of 33), but requires more iterations than with MCP starting values) 

priormean = startcoefs
cutoff = quantile(abs(priormean)[priormean!=0], 0.5) 
priormean[abs(priormean)<=cutoff] = 0 # we shrink small coefficients to zero & keep larger ones at their original value

nonzero = startcoefs!=0
lambda_preset = preset.lambda(IC = "mbic", # options aic, bic, gic, ebic, mbic 
                              p = ncol(X), # note that we use original problem size
                              n = nrow(X))

system.time(L0glm_gaus_presetlambda_full <- L0glm.fit(
  X = X, # [, nonzero, drop=F], 
  y = y, 
  family = family, 
  lambda = "mbic", # mbic because high dimensional 
  nonnegative = FALSE,
  dispersion = dispersion,
  # prior.mean = priormean, # [nonzero],
  start = startcoefs, # [nonzero], 
  warm.start = 0.25, # 1-sum(startcoefs!=0)/p, # 0.1 would return too small a solution, e.g. with just 17 true pos, values >0.7 give 34 true pos 
  # TODO make sure that with warm.start=1 vars with zero coefs are efficiently removed from the start
  gamma = 2, 
  normalize = FALSE, # normalize = T doesn't work yet - too sparse then
  # dispersion = dispersion,
  solver = "eigen", # note: solver osqp better than eigen for very small lambda values))
  method = 3L, # faster than 7 here
  relaxed = FALSE)) # 22s

L0glm_gaus_presetlambda_full$dispersion # 1

L0glm_gaus_presetlambda_full_coefs = L0glm_gaus_presetlambda_full$coefficients/scaling

beta_est1_full = L0glm_gaus_presetlambda_full_coefs
sum(beta_est1_full!=0&beta_true!=0)  # 68
table(beta_est1_full!=0, beta_true!=0, 
      dnn=c("estimated beta nonzero", "true beta nonzero"))
#                    true beta nonzero
# estimated beta nonzero  FALSE   TRUE
#                  FALSE  99900    44
#                   TRUE      0    56

cor(beta_est1_full, beta_true) # 0.9988259

plot(x=beta_est1_full[beta_true!=0&beta_est1_full!=0], 
     y=beta_true[beta_true!=0&beta_est1_full!=0], pch=16, col="blue",
     xlab="Estimated coefficients", ylab="True coefficients",
     main="L0glm with MCP starting values and warm start=0.7") # true positives
abline(a=0, b=1, col="lightgrey")
points(x=beta_est1_full[beta_true!=0&beta_est1_full==0], 
       y=beta_true[beta_true!=0&beta_est1_full==0], 
       pch=16, col="purple") # false negatives
points(x=beta_est1_full[beta_true==0&beta_est1_full!=0], 
       y=beta_true[beta_true==0&beta_est1_full!=0], 
       pch=16, col="red") # false positives !!!!!
points(x=beta_est1_full[beta_true==0&beta_est1_full==0][[1]], 
       y=beta_true[beta_true==0&beta_est1_full==0][[1]], 
       pch=16, col="green") # true negatives
# graph2png(file="..//benchmarks/plots/adaptive_L0glm/3_fullL0glm_fit_MCP starting vals_p_1_million_n_500_50 nonzero.png", 
#          width=7, height=5)

# check metrics & quality of solution with preset lambda
system.time(metrics_L0 <- data.frame(do.call(cbind, get.metrics(beta = beta_est1_full,  # some metrics
                                                                X = X, 
                                                                y = y,
                                                                family = gaussian(identity), 
                                                                beta_true = beta_true))))  # 0.5s
metrics_L0
#        p   n m deviance deviance_true   loglik     MD_coefs    MAD_coefs  sMAPE_coefs bias_coefs
# 1 100001 500 1 52963.13      4.109432 98068.84 3.939592e-05 0.0001755479 0.0009360564    1.00724
#    R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true
# 1 0.9496736  0.03153921   0.09479256      0.1001008       1.00724   0.9929054               100
#   support_size_estimated estimated_over_true_support_size TP FP    TN FN precision recall
# 1                     56                             0.56 56  0 99900 44         1   0.56
#   specificity accuracy  F1_score       gic      ebic       bic       aic    mbic
# 1           1  0.99956 0.7179487 -194938.8 -194825.2 -195783.5 -196023.7 -194629

pos_TR <- as.numeric(beta_true != 0) # true positives
pos <- as.numeric(beta_est1_full != 0) # estimates positives
TPR <- sum(pos[which(pos_TR == 1)])/sum(pos_TR) # True positive rate
TPR # 0.56
FDR <- sum(pos[which(pos_TR != 1)])/max(sum(pos), 1) # False discovery rate
FDR # 0
MSPE <- sqrt(sum((X %*% beta_est1_full - y)^2)/n) #Mean squared prediction error
MSPE # 10.29205




# 1.5 ITERATIVE ADAPTIVE L0GLM FIT ####

# For the high dimensional case a solution with fewer false positives
# can be obtained by using the 1st iteration
# solution and then using the L0glm output coefficients as starting values in 
# subsequent iterations. This would be sort of like an "adaptive L0 glm" then, where
# the original starting values that are used would themselves also be optimised.
# Given that we start from sensible starting values we can specify a prior mean here.

startcoefs = beta_est1*scaling # we start from initial L0glm estimate
# NOTE: starting from the MCP solution also works, but takes more iterations
# to converge
# starting from the iterative adaptive MCP solution converges very quickly
# but I need to find out how to do those fits on X subsetting to nonzero coefs

system.time({
  for (iter in (1:20)) {
    priormean = startcoefs
    cutoff = quantile(abs(priormean)[priormean!=0], 0.05) # find a better smooth functional form instead of a hard cutoff?
    priormean[abs(priormean)<=cutoff] = 0 # still using zero prior center for small coefficients
    
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
                                               warm.start = 0.7, # note 0.5-0.8 also work well
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
         y=beta_true[beta_true!=0&beta_est1B!=0], pch=16, col="blue",
         xlab="Estimated coefficients", ylab="True coefficients",
         main=paste0("iterative adaptive L0glm on MCP preselected variables iter ", iter)) # true positives
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
}) # 3.6s

sum(beta_est1B!=0) # 61
sum(beta_est1B!=0&beta_true!=0)  # 61
L0glm.out_updated$dispersion # 3.486745
# these effect sizes in terms of Cohen's d end up as false negatives and are not detected
range(abs(beta_true)[beta_est1B==0&beta_true!=0]) 
# 0.004638211 0.359392207
# while these effect sizes in terms of Cohen's d are detected
range(abs(beta_true)[beta_est1B!=0&beta_true!=0]) 
# 0.3093196 2.7166332

table(beta_est1B!=0, beta_true!=0, 
      dnn=c("estimated beta nonzero", "true beta nonzero"))
#                         true beta nonzero
# estimated beta nonzero  FALSE   TRUE
# FALSE                   99900    39
# TRUE                        0    61

plot(x=beta_est1B[beta_true!=0&beta_est1B!=0], 
     y=beta_true[beta_true!=0&beta_est1B!=0], pch=16, col="blue",
     xlab="Estimated coefficients", ylab="True coefficients",
     main=paste0("iterative adaptive L0glm on MCP preselected variables")) # true positives
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
# graph2png(file="..//benchmarks/plots/adaptive_L0glm/4_L0glm_on_MCP_presel_vars_p_1_million_n_500_50 nonzero.png", 
#           width=7, height=5)

system.time(metrics_adL0 <- data.frame(do.call(cbind, get.metrics(beta = beta_est1B,  # some metrics
                                                                  X = X, 
                                                                  y = y,
                                                                  family = gaussian(identity), 
                                                                  beta_true = beta_true)))) 
metrics_adL0
#        p   n m deviance deviance_true   loglik     MD_coefs    MAD_coefs  sMAPE_coefs bias_coefs
# 1 100001 500 1 52877.26      2.523302 98064.25 9.720449e-06 0.0001467047 0.0008469979   0.993062
#   R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true
# 1  0.96492  0.02779224   0.08722527      0.1098326      0.993062   0.9928395               100
#   support_size_estimated estimated_over_true_support_size TP FP    TN FN precision recall
# 1                     61                             0.61 61  0 99900 39         1   0.61
#   specificity accuracy F1_score       gic      ebic       bic       aic      mbic
# 1           1  0.99961 0.757764 -194824.5 -194700.9 -195743.2 -196004.5 -194487.5

pos_TR <- as.numeric(beta_true != 0) # true positives
pos <- as.numeric(beta_est1B != 0) # estimates positives
TPR <- sum(pos[which(pos_TR == 1)])/sum(pos_TR) # True positive rate
TPR # 0.61
FDR <- sum(pos[which(pos_TR != 1)])/max(sum(pos), 1) # False discovery rate
FDR # 0
MSPE <- sqrt(sum((X %*% beta_est1B - y)^2)/n) #Mean squared prediction error
MSPE # 10.2837


# 1.6 abess fit on MCP pre-selected variables ####
library(abess)
startcoefs = beta_est0*scaling # MCP starting coefficients work best

selcoefs = which(startcoefs!=0)
p_selcoefs = length(selcoefs)
system.time(fit_abessafterMCP_selcoefs <- abess(x = X[,selcoefs,drop=F], 
                                                y = y, 
                                                family = "gaussian",
                                                tune.path = "sequence", 
                                                support.size = c(1:(p_selcoefs-1)),
                                                lambda = 0,
                                                ic.scale = (1/log(n))*preset.lambda(IC="mbic", p=p, n=n), # to convert BIC criterion to MBIC
                                                warm.start = TRUE, 
                                                tune.type = "bic"))  # or bic, ebic, aic or cv
# 2.78s
plot(fit_abessafterMCP_selcoefs, type="tune")  # optimal support size 
beta_abessafterMCP_selcoefs = as.matrix(extract(fit_abessafterMCP_selcoefs)$beta) # coefficient matrix for best subset
beta_abessafterMCP = Matrix(rep(0, p), sparse=T)
beta_abessafterMCP[selcoefs] = beta_abessafterMCP_selcoefs
beta_abessafterMCP = as.vector(beta_abessafterMCP)/scaling
# beta_abessafterMCP[beta_abessafterMCP<0] = 0
sum(beta_abessafterMCP!=0) # 64
plot(x=beta_abessafterMCP[beta_true!=0&beta_abessafterMCP!=0], 
     y=beta_true[beta_true!=0&beta_abessafterMCP!=0], pch=16, col="blue",
     xlab="Estimated coefficients", ylab="True coefficients",
     main=paste0("abess best subset on MCP preselected variables")) # true positives
abline(a=0, b=1, col="lightgrey")
points(x=beta_abessafterMCP[beta_true!=0&beta_abessafterMCP==0], 
       y=beta_true[beta_true!=0&beta_abessafterMCP==0], 
       pch=16, col="purple") # false negatives
points(x=beta_abessafterMCP[beta_true==0&beta_abessafterMCP!=0], 
       y=beta_true[beta_true==0&beta_abessafterMCP!=0], 
       pch=16, col="red") # false positives
points(x=beta_abessafterMCP[beta_true==0&beta_abessafterMCP==0][[1]], 
       y=beta_true[beta_true==0&beta_abessafterMCP==0][[1]], 
       pch=16, col="green") # true negatives
# graph2png(file="..//benchmarks/plots/adaptive_L0glm/5_abess_on_MCP_presel_vars_p_1_million_n_500_50 nonzero.png", 
#           width=7, height=5)
table(beta_abessafterMCP!=0, beta_true!=0, 
      dnn=c("estimated beta nonzero", "true beta nonzero"))
#                         true beta nonzero
# estimated beta nonzero  FALSE   TRUE
#                  FALSE  99900    36
#                   TRUE      0    64


# check metrics & quality of solution 
system.time(metrics_abess_afterMCP <- data.frame(do.call(cbind, get.metrics(beta = beta_abessafterMCP,  # some metrics
                                                                            X = X, 
                                                                            y = y,
                                                                            family = gaussian(identity), 
                                                                            beta_true = beta_true))))  # 0.5s
metrics_abess_afterMCP
#        p   n m deviance deviance_true   loglik     MD_coefs    MAD_coefs sMAPE_coefs bias_coefs
# 1 100001 500 1 52842.79      1.874567 98061.41 5.384027e-06 0.0001291263  0.00077672  0.9935329
#    R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true
# 1 0.9715504  0.02578059   0.07519904     0.08862496     0.9935329   0.9946001               100
#   support_size_estimated estimated_over_true_support_size TP FP    TN FN precision recall
# 1                     64                             0.64 64  0 99900 36         1   0.64
#   specificity accuracy  F1_score       gic      ebic       bic       aic      mbic
# 1           1  0.99964 0.7804878 -194755.7 -194626.1 -195718.9 -195992.8 -194402.4

pos_TR <- as.numeric(beta_true != 0) # true positives
pos <- as.numeric(beta_abessafterMCP != 0) # estimates positives
TPR <- sum(pos[which(pos_TR == 1)])/sum(pos_TR) # True positive rate
TPR # 0.64
FDR <- sum(pos[which(pos_TR != 1)])/max(sum(pos), 1) # False discovery rate
FDR # 0
MSPE <- sqrt(sum((X %*% beta_abessafterMCP - y)^2)/n) #Mean squared prediction error
MSPE # 10.28035


# 1.7 FULL abess FIT using MCP preselected set as initial active set ####
system.time(fit_abess <- abess(x = X, 
                               y = y, 
                               family = "gaussian",
                               tune.path = "sequence", 
                               support.size = c(1:(p_selcoefs-1)),
                               lambda = 0,
                               init.active.set = as.integer(which(startcoefs!=0)),
                               warm.start = TRUE, 
                               tune.type = "gic")) # 16s # gic, ebic or mbic all work
plot(fit_abess, type="tune")  # optimal support size 
beta_abess = as.vector(extract(fit_abess)$beta)/scaling 
sum(beta_abess!=0) # 35
plot(x=beta_abess[beta_true!=0&beta_abess!=0], 
     y=beta_true[beta_true!=0&beta_abess!=0], pch=16, col="blue",
     xlab="Estimated coefficients", ylab="True coefficients",
     main=paste0("abess best subset fit")) # true positives
abline(a=0, b=1, col="lightgrey")
points(x=beta_abess[beta_true!=0&beta_abess==0], 
       y=beta_true[beta_true!=0&beta_abess==0], 
       pch=16, col="purple") # false negatives
points(x=beta_abess[beta_true==0&beta_abess!=0], 
       y=beta_true[beta_true==0&beta_abess!=0], 
       pch=16, col="red") # false positives
points(x=beta_abess[beta_true==0&beta_abess==0][[1]], 
       y=beta_true[beta_true==0&beta_abess==0][[1]], 
       pch=16, col="green") # true negatives
# graph2png(file="..//benchmarks/plots/adaptive_L0glm/7_abess_fit_p_1_million_n_500_50 nonzero.png", 
#          width=7, height=5)
table(beta_abess!=0, beta_true!=0, 
      dnn=c("estimated beta nonzero", "true beta nonzero"))
# true beta nonzero
# estimated beta nonzero  FALSE   TRUE
#                  FALSE  99894    71
#                  TRUE       6    29

# check metrics & quality of solution 
system.time(metrics_abess <- data.frame(do.call(cbind, get.metrics(beta = beta_abess,  # some metrics
                                                                   X = X, 
                                                                   y = y,
                                                                   family = gaussian(identity), 
                                                                   beta_true = beta_true))))  # 0.5s
metrics_abess
#        p   n m deviance deviance_true  loglik     MD_coefs    MAD_coefs sMAPE_coefs bias_coefs
# 1 100001 500 1 53440.21      15.19178 98087.6 5.677336e-05 0.0004657301 0.001583043  0.9009427
#    R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true
# 1 0.7135312 -0.05182279     0.328558      0.4658378     0.9009427   0.9309704               100
#   support_size_estimated estimated_over_true_support_size TP FP    TN FN precision recall
# 1                     35                             0.35 29  6 99894 71 0.8285714   0.29
#   specificity accuracy  F1_score     gic      ebic       bic       aic      mbic
# 1   0.9999399  0.99923 0.4296296 -195418 -195346.3 -195951.5 -196103.2 -195222.4




# 1.8 FULL abess fit not specifying the initial active MCP selected set ####

p_max <- n/log(n) # just to make it faster
system.time(fit_abess2 <- abess(x = X, 
                                y = y, 
                                family = "gaussian",
                                tune.path = "sequence", 
                                support.size = c(1:p_max),
                                lambda = 0,
                                # init.active.set = as.integer(which(startcoefs!=0)),
                                warm.start = TRUE, 
                                tune.type = "gic")) # 4.5s
plot(fit_abess2, type="tune")  # optimal support size 
beta_abess2 = as.vector(extract(fit_abess2)$beta)/scaling 
sum(beta_abess2!=0) # 35
plot(x=beta_abess2[beta_true!=0&beta_abess2!=0], 
     y=beta_true[beta_true!=0&beta_abess2!=0], pch=16, col="blue",
     xlab="Estimated coefficients", ylab="True coefficients",
     main=paste0("abess best subset fit")) # true positives
abline(a=0, b=1, col="lightgrey")
points(x=beta_abess2[beta_true!=0&beta_abess2==0], 
       y=beta_true[beta_true!=0&beta_abess2==0], 
       pch=16, col="purple") # false negatives
points(x=beta_abess2[beta_true==0&beta_abess2!=0], 
       y=beta_true[beta_true==0&beta_abess2!=0], 
       pch=16, col="red") # false positives
points(x=beta_abess2[beta_true==0&beta_abess2==0][[1]], 
       y=beta_true[beta_true==0&beta_abess2==0][[1]], 
       pch=16, col="green") # true negatives
# graph2png(file="..//benchmarks/plots/adaptive_L0glm/7_abess_fit2_p_1_million_n_500_50 nonzero.png", 
#          width=7, height=5)
table(beta_abess2!=0, beta_true!=0, 
      dnn=c("estimated beta nonzero", "true beta nonzero"))
# result with ebic
# true beta nonzero
# estimated beta nonzero  FALSE   TRUE
#                  FALSE  99894    71
#                  TRUE       6    29


# check metrics & quality of solution 
system.time(metrics_abess2 <- data.frame(do.call(cbind, get.metrics(beta = beta_abess2,  # some metrics
                                                                    X = X, 
                                                                    y = y,
                                                                    family = gaussian(identity), 
                                                                    beta_true = beta_true))))  # 0.5s
metrics_abess2
#        p   n m deviance deviance_true  loglik     MD_coefs    MAD_coefs sMAPE_coefs bias_coefs
# 1 100001 500 1 53440.21      15.19178 98087.6 5.677336e-05 0.0004657301 0.001583043  0.9009427
#    R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true
# 1 0.7135312 -0.05182279     0.328558      0.4658378     0.9009427   0.9309704               100
#   support_size_estimated estimated_over_true_support_size TP FP    TN FN precision recall
# 1                     35                             0.35 29  6 99894 71 0.8285714   0.29
#   specificity accuracy  F1_score     gic      ebic       bic       aic      mbic
# 1   0.9999399  0.99923 0.4296296 -195418 -195346.3 -195951.5 -196103.2 -195222.4



# 1.9 L0Learn fit on MCP pre-selected variables ####
library(L0Learn)
startcoefs = beta_est0*scaling # MCP starting coefficients

selcoefs = which(startcoefs!=0)
p_selcoefs = length(selcoefs)
W = rep(1, nrow(X))

system.time(fit_L0Learn_selcoefs <- L0Learn.fit(x=as.matrix(Diagonal(x=sqrt(W)) %*% X[,selcoefs,drop=F]),  # note: no support for sparse matrices
                                                y=as.vector(y*sqrt(W)),
                                                loss = "SquaredError",
                                                penalty="L0", algorithm="CD",
                                                maxSuppSize=(p_selcoefs-1), 
                                                intercept=FALSE, 
                                                nLambda=100, # NOTE: 100 gave much better results than 20
                                                # activeSetNum = 3, maxSwaps = 100,
                                                autoLambda=TRUE
                                                # lambdaGrid=list(10^seq(-5,5,length.out=nlam)/ # DECREASING sequence recommended, but INCREASING works better under nonnegativity, see https://cran.r-project.org/web/packages/L0Learn/vignettes/L0Learn-vignette.html
                                                #                   (norm(sim$y*sqrt(1/(sim$y+1)),"2")^2)) 
))  # lambda of 1 is scaled by 1/(norm(weighted y,"2")) because L0Learn scales y by L2 norm of y
# 0.16s
plot(fit_L0Learn_selcoefs)

# log-likelihood method for gaussian L0Learn fit,
# cf logLik.ordinis in https://github.com/jaredhuling/ordinis/blob/master/R/ordinis_methods.R
logLik.L0Learn = function(object, weights, X, y) {
  n  = object$n
  p = object$p
  df = object$suppSize[[1]]
  # df = nr of nonzero coefficients (incl intercept)
  
  # Compute y_pred for all columns of object$beta at once
  y_pred = X %*% object$beta[[1]]
  
  # Compute squared residuals, multiply by weights, and sum over rows to get resid.ss
  resid.ss = colSums(weights * (matrix(rep(y, ncol(y_pred)), ncol = ncol(y_pred) ) -
                                  y_pred) ^ 2)
  
  rdf = n - df
  logL = -0.5 * n * (log(2 * pi) - log(rdf) + log(resid.ss)) - 0.5 * rdf
  attr(logL, "df") = df  # nr of nonzero coefficients
  attr(logL, "nobs") = n
  attr(logL, "p") = object$dim[[2]]
  return(logL)
}

LLs = logLik(fit_L0Learn_selcoefs, weights=W, X[,selcoefs,drop=F], y)
edfs = attr(LLs, "df")
LLs = as.vector(LLs)
crit = "mbic"
IC = get.IC(loglik = LLs, n = n, p = p_selcoefs, edf = edfs)[,crit]
beta_L0Learn_selcoefs = coef(fit_L0Learn_selcoefs,
                             lambda = fit_L0Learn_selcoefs$lambda[[1]][which.min(IC)],
                             gamma = 0)
beta_L0Learn_after_MCP = matrix(0, p)
beta_L0Learn_after_MCP[selcoefs] = beta_L0Learn_selcoefs
beta_L0Learn_after_MCP = as.vector(beta_L0Learn_after_MCP)/scaling
sum(beta_L0Learn_after_MCP!=0) # 77

plot(x=beta_L0Learn_after_MCP[beta_true!=0&beta_L0Learn_after_MCP!=0], 
     y=beta_true[beta_true!=0&beta_L0Learn_after_MCP!=0], pch=16, col="blue",
     xlab="Estimated coefficients", ylab="True coefficients",
     main=paste0("L0Learn on MCP preselected variables")) # true positives
abline(a=0, b=1, col="lightgrey")
points(x=beta_L0Learn_after_MCP[beta_true!=0&beta_L0Learn_after_MCP==0], 
       y=beta_true[beta_true!=0&beta_L0Learn_after_MCP==0], 
       pch=16, col="purple") # false negatives
points(x=beta_L0Learn_after_MCP[beta_true==0&beta_L0Learn_after_MCP!=0], 
       y=beta_true[beta_true==0&beta_L0Learn_after_MCP!=0], 
       pch=16, col="red") # false positives
points(x=beta_L0Learn_after_MCP[beta_true==0&beta_L0Learn_after_MCP==0][[1]], 
       y=beta_true[beta_true==0&beta_L0Learn_after_MCP==0][[1]], 
       pch=16, col="green") # true negatives
# graph2png(file="..//benchmarks/plots/adaptive_L0glm/6_L0Learn_on_MCP_presel_vars_p_1_million_n_500_50 nonzero.png", 
#           width=7, height=5)
table(beta_L0Learn_after_MCP!=0, beta_true!=0, 
      dnn=c("estimated beta nonzero", "true beta nonzero"))
#                         true beta nonzero
# estimated beta nonzero  FALSE   TRUE
#                  FALSE  99891    32
#                   TRUE      9    68

# check metrics & quality of solution 
system.time(metrics_L0Learn_after_MCP <- data.frame(do.call(cbind, get.metrics(beta = beta_L0Learn_after_MCP,  # some metrics
                                                                    X = X, 
                                                                    y = y,
                                                                    family = gaussian(identity), 
                                                                    beta_true = beta_true))))  # 0.5s
metrics_L0Learn_after_MCP
#        p   n m deviance deviance_true   loglik      MD_coefs    MAD_coefs  sMAPE_coefs bias_coefs
# 1 100001 500 1 52798.01      1.315228 98048.63 -6.045575e-06 0.0001467335 0.0008837386  0.9946475
#    R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true
# 1 0.9681505  0.00725907    0.1042162      0.3165437     0.9946475   0.9862835               100
#   support_size_estimated estimated_over_true_support_size TP FP    TN FN precision recall
# 1                     77                             0.77 68  9 99891 32 0.8831169   0.68
#   specificity accuracy  F1_score       gic      ebic       bic       aic      mbic
# 1   0.9999099  0.99959 0.7683616 -194456.7 -194301.2 -195612.5 -195941.3 -194032.8




# 1.10 sparsevb Spike-and-Slab Variational Bayes fit ####
library(sparsevb)
system.time(svbfit <- svb.fit(X, y, family = "linear")) # 168s
# NOTE noise is estimated using estimateSigma from selectiveInference package
beta_svb_posterior_mean <- svbfit$mu * svbfit$gamma # approximate posterior mean = variational mean mu x inclusion probability gamma

pos <- as.numeric(svbfit$gamma > 0.5) # significant coefficients (included >50% of the time)
sum(pos) # 6

### Assess the quality of the posterior estimates ###

pos_TR <- as.numeric(beta_true != 0) # true positives
TPR <- sum(pos[which(pos_TR == 1)])/sum(pos_TR) # True positive rate
TPR # 0.06
FDR <- sum(pos[which(pos_TR != 1)])/max(sum(pos), 1) # False discovery rate
FDR # 0
L2 <- sqrt(sum((beta_svb_posterior_mean - beta_true)^2)) #L_2-error 
MSPE <- sqrt(sum((X %*% beta_svb_posterior_mean - y)^2)/n) #Mean squared prediction error
MSPE # 9.015753

# check metrics & quality of solution with preset lambda
system.time(metrics_svbfit <- data.frame(do.call(cbind, get.metrics(beta = beta_svb_posterior_mean,  # some metrics
                                                                     X = X, 
                                                                     y = y,
                                                                     family = gaussian(identity), 
                                                                     beta_true = beta_true))))  # 0.5s
metrics_svbfit
#        p   n m deviance deviance_true    loglik     MD_coefs   MAD_coefs sMAPE_coefs bias_coefs
# 1 100001 500 1  40641.9      3684.849 -1808.956 -0.000770067 0.002272403    1.999945 0.09399519
#   R2_coefs  MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true
# 1 0.360535 -0.000770067  0.002272403       1.999945    0.09399519    0.360535               100
#   support_size_estimated estimated_over_true_support_size  TP    FP TN FN precision recall
# 1                  1e+05                             1000 100 99900  0  0     0.001      1
#   specificity accuracy    F1_score     gic    ebic      bic      aic    mbic
# 1           0    0.001 0.001998002 2106940 2306228 625084.9 203619.9 2650433



# 1.11 ebreg empirical Bayes method for high-dimensional linear models ####
library(ebreg)
system.time(ebregfit <- ebreg(y = y, X = X, XX = NULL, pred=FALSE)) # XXs


# 1.12 SSLASSO Spike-and-Slab LASSO fit ####
library(SSLASSO)
system.time(SSLASSOfit <- SSLASSO(X, y, penalty = "adaptive",
                                  variance = "unknown")) # 45s
length(SSLASSOfit$model) # 6
# plot(SSLASSOfit)
# beta_SSLASSO <- colMeans(SSLASSOfit$beta) # mean posterior

SSLASSOfit$beta
metrics_SLASSO <- do.call(rbind, apply(SSLASSOfit$beta, 2, function (column) data.frame(do.call(cbind, get.metrics(beta = column,  # some metrics
                                                                    X = X, 
                                                                    y = y,
                                                                    family = gaussian(identity), 
                                                                    beta_true = beta_true)))))
plot(metrics_SLASSO$deviance_true, type="l")
best = which.min(metrics_SLASSO$deviance_true)

beta_SLASSO = SSLASSOfit$beta[,best]
# check metrics & quality of solution with preset lambda
system.time(metrics_SLASSO <- data.frame(do.call(cbind, get.metrics(beta = beta_SLASSO,  # some metrics
                                                                    X = X, 
                                                                    y = y,
                                                                    family = gaussian(identity), 
                                                                    beta_true = beta_true))))  # 0.5s
metrics_SLASSO
#        p   n m deviance deviance_true   loglik     MD_coefs   MAD_coefs sMAPE_coefs bias_coefs
# 1 100001 500 1 32510.93      21108.74 98240.85 -0.002145881 0.004246386 0.001990509 0.04108737
#    R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true
# 1 0.3590199   -35.63765     59.79792       1.841813    0.04108737   0.9832617               100
#   support_size_estimated estimated_over_true_support_size TP FP    TN FN precision recall
# 1                      6                             0.06  6  0 99900 94         1   0.06
#   specificity accuracy  F1_score       gic      ebic       bic       aic      mbic
# 1           1  0.99906 0.1132075 -196334.5 -196320.5 -196438.2 -196467.7 -196296.4


# ALSO CHECK
# https://search.r-project.org/CRAN/refmans/BBSSL/html/SSLASSO_2.html
# https://search.r-project.org/CRAN/refmans/tidyfit/html/dot-fit.spikeslab.html
# https://search.r-project.org/CRAN/refmans/ScaleSpikeSlab/html/spike_slab_linear.html
# ?




# 1.10 COMPARISON OF METRICS ####

benchmarks = rbind(data.frame(method="MCP", metrics_MCP),
                   data.frame(method="iterative_adaptive_MCP", metrics_adMCP),
                   data.frame(method="L0glm_after_MCP", metrics_L0),
                   data.frame(method="iterative adaptive L0glm_after_MCP", metrics_adL0),
                   data.frame(method="abess", metrics_abess),
                   data.frame(method="abess_after_MCP", metrics_abess_afterMCP),
                   data.frame(method="L0Learn_after_MCP", metrics_L0Learn),
                   data.frame(method="LASSO", metrics_LASSO))
benchmarks = benchmarks[order(benchmarks$deviance_true),]
benchmarks
#                               method       p   n m deviance deviance_true   loglik      MD_coefs
# 3                    L0glm_after_MCP 1000001 500 1 21427.94     0.3106053 998313.1  2.650453e-07
# 6                    abess_after_MCP 1000001 500 1 21417.30     0.3581132 998308.2  9.717205e-07
# 2             iterative_adaptive_MCP 1000001 500 1 21438.59     0.4390093 998316.9  6.466617e-08
# 5                              abess 1000001 500 1 21435.51     0.4855943 998318.0 -1.466363e-07
# 4 iterative adaptive L0glm_after_MCP 1000001 500 1 21438.35     0.5088160 998317.9  4.126374e-10
# 1                                MCP 1000001 500 1 21415.61     0.6292156 998275.2 -9.841375e-07
# 7                  L0Learn_after_MCP 1000001 500 1 21397.74     0.6421191 998277.4 -7.387322e-07
# 8                              LASSO 1000001 500 1 21607.55     1.2141603 998015.0 -3.219889e-06
#      MAD_coefs  sMAPE_coefs bias_coefs  R2_coefs   MD_selcoefs MAD_selcoefs sMAPE_selcoefs
# 3 3.189592e-06 3.100965e-05  0.9915152 0.9916190 -0.0069207573   0.04859435     0.13183282
# 6 3.952252e-06 3.769237e-05  0.9949028 0.9884211  0.0156488416   0.06601067     0.31842729
# 2 3.813945e-06 3.469761e-05  0.9937225 0.9869970  0.0003857053   0.05228158     0.07934134
# 5 4.001710e-06 3.616757e-05  0.9887556 0.9852913  0.0035141707   0.05003584     0.06568395
# 4 4.139005e-06 3.659235e-05  0.9904504 0.9845558  0.0079701990   0.05419631     0.07855603
# 1 8.962904e-06 9.444037e-05  1.0200814 0.9687125 -0.0188729902   0.10724618     1.00579431
# 7 8.621173e-06 8.974424e-05  1.0093180 0.9712289 -0.0160667829   0.10552674     0.96951675
# 8 2.706061e-05 6.745924e-04  1.3427367 0.8656587 -0.0072952627   0.06775434     1.88800121
#   bias_selcoefs R2_selcoefs support_size_true support_size_estimated
# 3     0.9915152   0.9962630                50                     38
# 6     0.9949028   0.9917800                50                     43
# 2     0.9937225   0.9962767                50                     34
# 5     0.9887556   0.9969646                50                     33
# 4     0.9904504   0.9962204                50                     33
# 1     1.0200814   0.9711268                50                     76
# 7     1.0093180   0.9736495                50                     74
# 8     1.3427367   0.8954008                50                    334
#   estimated_over_true_support_size TP  FP     TN FN  precision recall specificity accuracy
# 3                             0.76 37   1 999949 13 0.97368421   0.74    0.999999 0.999986
# 6                             0.86 38   5 999945 12 0.88372093   0.76    0.999995 0.999983
# 2                             0.68 34   0 999950 16 1.00000000   0.68    1.000000 0.999984
# 5                             0.66 33   0 999950 17 1.00000000   0.66    1.000000 0.999983
# 4                             0.66 33   0 999950 17 1.00000000   0.66    1.000000 0.999983
# 1                             1.52 41  35 999915  9 0.53947368   0.82    0.999965 0.999956
# 7                             1.48 41  33 999917  9 0.55405405   0.82    0.999967 0.999958
# 8                             6.68 28 306 999644 22 0.08383234   0.56    0.999694 0.999672
#    F1_score      gic     ebic      bic      aic     mbic
# 3 0.8409091 -1995642 -1995549 -1996384 -1996548 -1995414
# 6 0.8172043 -1995506 -1995401 -1996343 -1996528 -1995249
# 2 0.8095238 -1995751 -1995667 -1996416 -1996564 -1995546
# 5 0.7951807 -1995778 -1995697 -1996425 -1996568 -1995579
# 4 0.7951807 -1995778 -1995696 -1996425 -1996568 -1995579
# 1 0.6507937 -1994607 -1994423 -1996072 -1996396 -1994158
# 7 0.6612903 -1994662 -1994483 -1996089 -1996405 -1994224
# 8 0.1458333 -1987575 -1986774 -1993948 -1995360 -1985620

# CONCLUSION: L0glm on MCP preselected variables has quite promising performance,
# iterative adaptive L0glm updates probably not worth it, though those
# iterations can successfully remove some false positives

# deviance relative to zero noise ground truth as well as
# accuracy, F1 score, R2 between estimated & true coefficients plus MD,
# MAD & sMAPE for selected coefs all best for L0glm solution
# abess and abess_after_MCP solutions are also very good

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
