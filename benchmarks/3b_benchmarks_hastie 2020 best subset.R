# TEST L0GLM FOR HASTIE ET AL 2020 PROBLEMS IN BEST SUBSET PACKAGE ####

library(remotes)
# remotes::install_github("jaredhuling/ordinis")
# remotes::install_github("tomwenseleers/L0glm/L0glm")
library(L0glm)
# devtools::load_all(".") # if you forked & cloned L0glm first to your laptop
library(export)
library(Matrix)
library(qlcMatrix)
library(matrixStats)
library(dplyr)
library(ncvreg)
library(SIS)
library(abess)
library(bestsubset)

# 1. TEST PROBLEM FROM BEST SUBSET PACKAGE & PAPER ####

n = 50; p = 1000 # Size of training set, and number of predictors
nval = n # Size of validation set: 50:50 split
# nrep = 10 # Number of repetitions for a given setting
seed = 0 # Random number generator seed
s <- n_nonzero <- 5 # Number of nonzero coefficients
type.vec = c(1:3,5) # Simulation settings to consider
beta.type <- type.vec[2] # simulation problem 2
rho.vec = c(0,0.35,0.7) # Pairwise predictor correlations
rho = rho.vec[1]
snr.vec = exp(seq(log(0.05),log(6),length=10)) # Signal-to-noise ratios 
snr <- snr.vec[7] # = 1.21644040
set.seed(seed)
data = bestsubset::sim.xy(n, p, nval, rho, s, beta.type, snr)
range(data$beta) # 0-1, 0-1, 0-10, -10 - 10, 0-1
range(data$x) # pos or neg for beta.type=1 2 3 4 5
apply(data$x, 2, function (col) norm(col, "2")) 
hist(data$beta)
# data$beta[data$beta<=0.5] = 0
beta_true = data$beta
beta_true_set = which(beta_true!=0)
length(beta_true_set) # 5
X = data$x
y = data$y
y_true = as.numeric(X %*% beta_true)
X_val = data$xval
y_val = data$yval
y_true_val = as.numeric(X_val %*% beta_true)
Sigma <- data$Sigma

## Population R2 ceiling (exact for these sims)
R2_max_pop <- snr / (1 + snr)
R2_max_pop # 0.5488261
R2_max_hat <- cor(y_val, y_true_val)^2 # empirical estimate
R2_max_hat #  0.5542138

w = rep(1, n) # observation weights
family = gaussian(identity)
nonnegative = all(beta_true>=0)

# error_sd = 1 
# dispersion = error_sd^2
suffix = paste0("n_",as.integer(n),"_p_",as.integer(p),"_s",s)

hist(beta_true[beta_true!=0])
range(beta_true)

# NOTE: best possible "oracle" model that would be achievable with L0glm & the gic criterion is
IC = "gic"
lambda_preset = preset.lambda(IC = IC, # options aic, bic, gic, mebic, ebic, mbic
                              p = ncol(X), # using original problem size
                              n = nrow(X))
system.time(L0glm_bestfit <- L0glm.fit(X = X, # X[,beta_true_set,drop=F], # X
                          y = y, 
                          family = gaussian(identity), 
                          lambda = lambda_preset, # "bic", # lambda_preset,  
                          nonnegative = nonnegative,
                          # prior.mean = beta_true, # [beta_true_set], # we can now specify our prior centers
                          start = beta_true, # rep(1, length(beta_true_set)), # beta_true, 
                          normalize = FALSE, # TODO fix normalize=TRUE
                          # dispersion = dispersion, 
                          solver = "glmnet", # eigen, glmnet or osqp, # note: constrained solvers glmnet & osqp better than eigen for very small lambda values
                          method = 1L, #
                          gamma = 2, # .1, # set it slightly above 2 for p>>n case
                          thresh = 1E-5,
                          relaxed = FALSE
))  # 10s
# beta_bestest = rep(0, ncol(X))
# beta_bestest[beta_true_set] = L0glm_bestfit$coefficients
beta_bestest = as.numeric(L0glm_bestfit$coefficients)
sum(beta_bestest!=0) # 3
# 40 with gic, hbic or ebic with p=1E6 lambda preset, 56, 53, 39 with aic, bic, mbic
min(abs(beta_bestest[beta_bestest!=0])) # 0.1891622
table(beta_bestest!=0, beta_true!=0, 
      dnn=c("estimated beta nonzero", "true beta nonzero"))

summary(lm(y_true ~ 0 + X %*% beta_bestest))$r.squared # true signal recovery R2=0.594359
cor(y_val, X_val %*% beta_bestest)^2 # validation set R2=0.33

beta_est0 = beta_bestest
plot(x=beta_est0[beta_true!=0&beta_est0!=0], 
     y=beta_true[beta_true!=0&beta_est0!=0], pch=16, col="blue",
     xlab="Estimated coefficients", ylab="True coefficients",
     main="oracle L0glm estimate (GIC tuned)") # true positives
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
method="oracle_L0glm"
# graph2png(file=paste0("..//benchmarks/plots/high_dimensional/",method,"_",suffix,".png"), 
#          width=7, height=5)





# 1.1 HYBRID ITERATIVE SURE INDEPENDENCE SCREENING (ISIS) - L0glm METHOD TO GET STABLE VARIABLE SELECTION IN p>>n REGIME, BASED ON MARGINAL LIKELIHOOD GAIN OF ADDING VARIABLE TO NULL MODEL ####

# WORKS & RECOVERS BETTER SOLUTION THAN ABESS & MCP BUT FASTER!

# IC = "gic"
# lambda_preset = preset.lambda(IC = IC, # options aic, bic, gic, mebic, ebic, mbic
#                               p = ncol(X), # note that we use original problem size
#                               n = nrow(X))  
# # iterating SIS & L0glm:
# beta_start <- rep(0, ncol(X))
# S <- NULL
# timing_SISL0glm <- system.time(for (iter in 1:1) {
#   sis <- sis_glm(X, y, family, 
#                  weights = NULL, 
#                  offset = X %*% beta_start,
#                  nsis = min(n, round(n/(log(log(n))*log(p)))), # max support size in abess = default
#                  return_stats = FALSE)
#   Snew <- sis$indices
#   S   <- union(S, Snew)
#   print(paste("Iteration", iter))
#   print(paste("Size of new SIS set:", length(Snew)))
#   print(paste("Total size of current SIS set:", length(S)))
#   for (j in 1:1) { # iterated twice to remove some false positives
#   L0glmfit_updated <- L0glm.fit(X = X[,S,drop=F], 
#                                  y = y, 
#                                  family = family, 
#                                  lambda = lambda_preset,  
#                                  nonnegative = nonnegative,
#                                  start = rep(1, length(S)), # if (iter==1) rep(1, length(S)) else beta_start[S], # better solution than NULL 
#                                  solver = "eigen", # , # eigen, glmnet or osqp, # note: constrained solvers glmnet & osqp better than eigen for very small lambda values
#                                  method = 1L, # Cholesky
#                                  gamma = 2.1 # taking gamma slightly higher than 2 for p>>n case
#   )  # 
#   beta_start = rep(0, ncol(X))
#   beta_start[S] = L0glmfit_updated$coefficients
#   S <- which(beta_start!=0)
#   }
#   print(paste("Number of nonzero coefficients:", sum(beta_start!=0)))
#   print(paste("Number of true positives:", length(intersect(S, beta_true_set))))
#   
# })[[3]] 
# timing_SISL0glm
# # 28s for 10 iter
# 
# system.time(L0glmfit_final <- L0glm.fit(X = X[,S,drop=F], 
#                               y = y, 
#                               family = family, 
#                               lambda = lambda_preset,  
#                               nonnegative = nonnegative,
#                               start = beta_start[S],  
#                               solver = "eigen", # , # eigen, glmnet or osqp, # note: constrained solvers glmnet & osqp better than eigen for very small lambda values
#                               method = 1L, # Cholesky
#                               gamma = 2
# ))  # 0.09s
# beta_final = rep(0, ncol(X))
# beta_final[S] = L0glmfit_final$coefficients
# # beta_final = L0glmfit_final$coefficients
# sum(beta_final!=0) # 38


IC <- "gic"
lambda_preset = preset.lambda(IC = IC, # options aic, bic, gic, mebic, ebic, mbic
                              p = ncol(X), # note that we use original problem size
                              n = nrow(X))  
system.time(sis <- sis_glm(X, y, family, 
                           weights = NULL, 
                           offset = NULL, 
                           nsis = n, # min(n, round(n/(log(log(n))*log(p)))), # max support size in abess = default
                           return_stats = TRUE)) # 2s
Tj=sis$stats$stat
Tj_signed=Tj*sis$stats$sign
pvals <- stats::pchisq(Tj, df = 1, lower.tail = FALSE)
pBH   <- stats::p.adjust(pvals, method = "BH")
pBonf <- stats::p.adjust(pvals, method = "bonferroni")
sel_fwer = sis_select_fwer(Tj, alpha=0.1)
length(sel_fwer$indices) # 111522
sel_fdr = sis_select_fdr(Tj, q = 0.01, cap = Inf)
length(sel_fdr$indices) # 388824
sel_knee = sis_select_knee(Tj, cap = n, floor = min(n, round(n/(log(log(n))*log(p)))) )
length(sel_knee$indices) # 35 -> works best??
length(intersect(sel_knee$indices, beta_true_set)) # TPs: 9
pBHcutoff <- 0.01
sel_knee_withFDR <- sel_knee$indices[pBH[sel_knee$indices]<pBHcutoff]
length(sel_knee_withFDR) # 35 -> does not help but also does not hurt
length(intersect(sel_knee_withFDR, beta_true_set)) # TPs: 9
pBonfcutoff <- 0.1
sel_knee_withBonf <- sel_knee$indices[pBonf[sel_knee$indices]<pBonfcutoff]
length(sel_knee_withBonf) # 35 -> does not help but also does not hurt
length(intersect(sel_knee_withBonf, beta_true_set)) # TPs: 9

# automatic way to decide on nsis?
length(sel_knee$indices) # 15
floor = min(n, round(n/(log(log(n))*log(p)))) 
floor # 20
ceiling = n
ceiling # 500
best = which.min(abs(length(sel_knee$indices)-c(floor, ceiling)))
# SIS = (best==1) # apply SIS or not? -> here TRUE
# SIS
nsis = c(floor, ceiling)[best]
nsis # 20

# iterating SIS & L0glm:
beta_start <- rep(0, ncol(X))
S <- Snew <- Sknee <- NULL
# S <- which(beta_start!=0)
# priormean <- rep(0, p)
disp <- NULL
timing_SISL0glm <- system.time(for (iter in 1:1) {
  sis <- sis_glm(X, y, family, 
                 weights = NULL, 
                 offset = X %*% beta_start,
                 nsis = nsis, # n, # min(n, round(n/(log(log(n))*log(p)))), # max support size in abess = default
                 return_stats = TRUE)
  Tj=sis$stats$stat
  Tj_signed=Tj*sis$stats$sign
  pvals <- stats::pchisq(Tj, df = 1, lower.tail = FALSE)
  pBH   <- stats::p.adjust(pvals, method = "BH")
  pBonf <- stats::p.adjust(pvals, method = "bonferroni")
  
  sel_knee = sis_select_knee(Tj, 
                             cap = n, 
                             floor = min(n, round(n/(log(log(n))*log(p)))) )
  Sknee <- sel_knee$indices
  print(paste("Iteration", iter))
  print(paste("Size of knee selected SIS set:", length(Sknee)))
  
  # pBonfcutoff <- 0.1
  # Snew <- Sknee[pBonf[Sknee]<pBonfcutoff]
  pBHcutoff <- 0.01
  Snew <- Sknee[pBH[Sknee]<pBHcutoff]
  
  print(paste("Size of knee selected SIS set with pBH<pBHcutoff:", length(Sknee[pBH[Sknee]<0.001])))
  print(paste("Size of knee selected SIS set with pBonf<pBonfcutoff:", length(Sknee[pBonf[Sknee]<0.001])))
  
  # Snew <- sis$indices
  # TODO fix Warning message:
  # In Matrix::Matrix(offset, n, 1, sparse = TRUE) :
  #  'nrow', 'ncol', 'byrow' disregarded for [mM]atrix 'data'
  S   <- union(S, Snew)
  print(paste("Size of new SIS set:", length(Snew)))
  print(paste("Number of true positives in new SIS set:", length(intersect(Snew, beta_true_set))))
  print(paste("Total size of current SIS set:", length(S)))
  print(paste("Number of true positives in current total SIS set:", length(intersect(S, beta_true_set))))
  #  for (j in 1:1) { # iterated twice to remove some false positives
  L0glmfit_updated <- L0glm.fit(X = X[,S,drop=F], 
                                y = y, 
                                family = family, 
                                lambda = lambda_preset,  
                                nonnegative = nonnegative,
                                start = Tj_signed[S]^2, # NULL, # rep(1, length(S)), # if (iter==1) rep(1, length(S)) else beta_start[S], # better solution than NULL 
                                # prior.mean = priormean[S],
                                # dispersion = disp,
                                solver = "eigen", # , # eigen, glmnet or osqp, # note: constrained solvers glmnet & osqp better than eigen for very small lambda values
                                method = 1L, # Cholesky
                                thresh = 1E-5,
                                gamma = 2.1 # taking gamma slightly higher than 2 for p>>n case
  )  # 
  sum(L0glmfit_updated$coefficients!=0) # 8
  beta_start = rep(0, ncol(X))
  beta_start[S] = L0glmfit_updated$coefficients
  print(paste("Number of nonzero coefficients:", sum(beta_start!=0)))
  print(paste("Number of true positives among nonzero L0glm coefficients:", length(intersect(which(beta_start!=0), beta_true_set))))
  S <- which(beta_start!=0)
  # largestcoefset <- which(abs(beta_start)>median(abs(beta_start[S]))) 
  # priormean[S] <- beta_start[largestcoefset]
  # disp <- L0glmfit_updated$dispersion
  #  }

})[[3]] 
timing_SISL0glm
# 33s - we should get 31 TPs
# gic & hbic return 31 TPs & 1 FP, ebic returns 29 TPs & 1 FP, while mbic returns only 16 TPs
# TODO check if it helps to set non-top n coefs to zero in L0glm.fit.internal

# system.time(L0glmfit_updated <- L0glm.fit(X = X[,S,drop=F], 
#                                         y = y, 
#                                         family = family, 
#                                         lambda = lambda_preset,  
#                                         nonnegative = nonnegative,
#                                         start = beta_start[S],  
#                                         # prior.mean = beta_start, # [S],
#                                         solver = "eigen", # , # eigen, glmnet or osqp, # note: constrained solvers glmnet & osqp better than eigen for very small lambda values
#                                         method = 1L, # Cholesky
#                                         thresh = 1E-5,
#                                         gamma = 2.1
# ))  # 0.09s

beta_final = beta_start # rep(0, ncol(X))
# beta_final[S] = L0glmfit_updated$coefficients

# beta_final[S] = L0glmfit_updated$coefficients[L0glmfit_updated$coefficients!=0]

sum(beta_final!=0) # 32


beta_SISL0glm = as.numeric(beta_final) # beta_start
table(beta_SISL0glm!=0, beta_true!=0, 
      dnn=c("estimated beta nonzero", "true beta nonzero"))

plot(x=beta_SISL0glm[beta_true!=0&beta_SISL0glm!=0], 
     y=beta_true[beta_true!=0&beta_SISL0glm!=0], pch=16, col="blue",
     xlab="Estimated coefficients", ylab="True coefficients",
     main=paste0("best subset (ISIS-L0glm, GIC tuned)")) # true positives
abline(a=0, b=1, col="lightgrey")
points(x=beta_SISL0glm[beta_true!=0&beta_SISL0glm==0], 
       y=beta_true[beta_true!=0&beta_SISL0glm==0], 
       pch=16, col="purple") # false negatives
points(x=beta_SISL0glm[beta_true==0&beta_SISL0glm!=0], 
       y=beta_true[beta_true==0&beta_SISL0glm!=0], 
       pch=16, col="red") # false positives
points(x=beta_SISL0glm[beta_true==0&beta_SISL0glm==0][[1]], 
       y=beta_true[beta_true==0&beta_SISL0glm==0][[1]], 
       pch=16, col="green") # true negatives
method="ISIS_L0glm"
graph2png(file=paste0("..//benchmarks/plots/high_dimensional/",method,"_",suffix,".png"), 
          width=7, height=5)

cor(beta_SISL0glm, beta_true)     # 0.9679615
sum(abs(beta_SISL0glm-beta_true)) # 5.539483

cor(data$y, data$x %*% beta_SISL0glm)^2 # squared correlation with noisy y=training set R2=0.3055049
cor(data$x %*% data$beta, data$x %*% beta_SISL0glm)^2 # measure of true signal recovery R2 in training set=0.1698546
cor(data$yval, data$xval %*% beta_SISL0glm)^2 # squared correlation with noisy y=validation set R2=0.01450341
cor(data$xval %*% data$beta, data$xval %*% beta_SISL0glm)^2 # measure of true signal recovery in validation set R2=0.03916114

beta_hat = beta_SISL0glm
rho_pop <- drop(t(beta_true) %*% Sigma %*% beta_hat) /
  sqrt( drop(t(beta_true) %*% Sigma %*% beta_true) *
          drop(t(beta_hat) %*% Sigma %*% beta_hat) )
rho_pop 
R2_sig_pop <- rho_pop^2
R2_sig_pop # 0.09750533 = squared population random-design signal correlation

cor(y_val, X_val %*% beta_SISL0glm)^2 # validation R2=0.01450341

fitted = X %*% beta_SISL0glm
summary(lm(y_true ~ 0 + fitted))$r.squared # signal recovery R2=0.1656509
plot(x = fitted, y = y_true, pch = 16, col = "steelblue")

fitted_val = X_val %*% beta_SISL0glm
summary(lm(y_true_val ~ 0 + fitted_val))$r.squared # validation signal recovery R2=0.02722646


# check metrics & quality of ISIS-L0glm solution
method = "ISIS_L0glm"
system.time(metrics_ISISL0glm <- data.frame(method=method,
                                            timing_s=timing_SISL0glm,
                                            do.call(cbind, get.metrics(beta = beta_SISL0glm,  # some metrics
                                                                       X = X, 
                                                                       y = y,
                                                                       family = gaussian(identity), 
                                                                       beta_true = beta_true)))) 
metrics_ISISL0glm



# NOTE: Sure Independence Screening via SIS::SIS is very slow
# and only supports canonical link & a couple of families
# and only dense design matrices & no support for weights or offset
# see Example in ?sis_glm

# TODO compare performance with simulation scenarios considered in
# https://www.jstatsoft.org/article/view/v083i02/0
# & empirical datasets



# 1.2A MCP fit using ordinis (box constraints support) ####

library(ordinis)
gamma_ord = 1.01 # default for MCP is 1.4 & 3.7 for SCAD, I am using value close to 1 to match best with my later L0 glm objective
timing_ordinis = system.time(ordinis_fit <- ordinis(X, y, 
               intercept=F,
               family="gaussian",
               penalty = "mcp",
               gamma = gamma_ord, 
               nlambda = 10L, 
               # lower.limits = rep(0, p), # force all coefficients to be positive
               # maxit = 10,
               # maxit.irls = 1,
               dfmax = n-1, # min(n-1, round(2*n/(log(log(n))*log(p)))),
               alpha = 1 #, # worse solution with 0.9 
               # tol = 0.1 # good enough to use as initialisation, default 1E-4 for gaussian, 1E-3 for rest
               ))[[3]] 
timing_ordinis 
# 42s with default tol (1E-4 for gaussian, 1E-3 for the rest), 21s with tol=0.1
plot(ordinis_fit)
max(ordinis_fit$niter) 
median(ordinis_fit$niter)
ordinis_fit$niter

# calculate ebic, mbic & gic/sic (cf. abess) criterion at each lambda in grid
LLs = logLik2(ordinis_fit, weights=w, X, y)
edfs = attr(LLs, "df")
LLs = as.vector(LLs)
IC.val = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,IC]
# ebic or gic would identity 5th lambda as best, mbic the 4th but too strict
best = which.min(IC.val) 
best # best lambda
secondbest = min(length(IC.val), best+1)
secondbest # 2nd best lambda (slightly less sparse) solution

# # NOTE: in ordinis LLs weights seem to be ignored (not important here):
# loglik = as.numeric(logLik(ordinis_fit))
# edf = colSums(ordinis_fit$beta!=0)
# IC = "ebic"
# IC.val = get.IC(loglik = LLs, n = n, p = p, edf = edf)[,IC]
# best = which.min(IC.val[-1])+1 
# best # 5th lambda best

ordinis_fit$lambda[best] # best lambda

colSums(ordinis_fit$beta!=0) # support sizes
colSums(ordinis_fit$beta!=0)[best] # selected support size
apply(ordinis_fit$beta, 2, function (col) sum(col[-1]!=0&beta_true!=0)) # TPs

coefs_ordinis = ordinis_fit$beta[,best] 
coefs_ordinis_set = which(coefs_ordinis!=0)-1
coefs_ordinis_secondbest = ordinis_fit$beta[,secondbest] 
# coefs_ordinis["V1"] # first non-intercept coef sometimes comes out as 1 or at ridiculously high value - is this a bug??
# prob a bug in ordinis? now setting it to zero, as it's usually a false pos
# coefs_ordinis[which(names(coefs_ordinis)=="V1")] = 0 
coefs_ordinis = coefs_ordinis[-which(names(coefs_ordinis)=="(Intercept)")] # dropping intercept
coefs_ordinis_secondbest = coefs_ordinis_secondbest[-which(names(coefs_ordinis_secondbest)=="(Intercept)")] # dropping intercept

# the max coefficients observed for any of the lambda values works best
# to use as starting coefficients in L0glm
coefs_ordinis_max = as.numeric(rowMax(ordinis_fit$beta))[-1] # drop intercept
# coefs_ordinis_max[2] = 0 # to drop weird first coefficient if needed

min(abs(coefs_ordinis)[coefs_ordinis!=0&beta_true!=0]) # det threshold
min(abs(coefs_ordinis_max)[coefs_ordinis_max!=0&beta_true!=0]) # det threshold

sum(coefs_ordinis!=0) # 7
sum(coefs_ordinis_max!=0) # 121
L1norm_coefsordinis = sum(abs(coefs_ordinis))
L0norm_coefsordinis = sum(coefs_ordinis!=0)
L0overL1norm_coefsordinis = L0norm_betatrue/L1norm_betatrue
L0overL1norm_coefsordinis 

range(coefs_ordinis) # -0.8558353  0.7161125
range(beta_true) # -0.7965632  0.6771529

table(coefs_ordinis!=0, beta_true!=0, 
      dnn=c("estimated beta nonzero", "true beta nonzero"))

table(coefs_ordinis_secondbest!=0, beta_true!=0, 
      dnn=c("estimated beta nonzero", "true beta nonzero"))

table(coefs_ordinis_max!=0, beta_true!=0, 
      dnn=c("estimated beta nonzero", "true beta nonzero"))

# TODO check how performance would be affected by choice of gamma
# (gamma close to 1 maybe results in objective function closest to that of L0glm?)

# TODO compare with SCAD (penalty="scad")

cor(coefs_ordinis, beta_true)              # 0.302131
cor(coefs_ordinis_secondbest, beta_true)   # 0.2087442
cor(coefs_ordinis_max, beta_true)          # 0.148738

fitted = X %*% coefs_ordinis
cor(y_true, fitted)^2 # population R2=0.3515819

fitted = X %*% coefs_ordinis_secondbest
cor(y_true, fitted)^2 # population R2=0.5877326

fitted_max = X %*% coefs_ordinis_max
cor(y_true, fitted_max)^2 # population R2=0.8956718

beta_est0 = coefs_ordinis
plot(x=beta_est0[beta_true!=0&beta_est0!=0], 
     y=beta_true[beta_true!=0&beta_est0!=0], pch=16, col="blue",
     xlab="Estimated coefficients", ylab="True coefficients",
     main="MCP penalized regression (ordinis, GIC tuned)") # true positives
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
method="ordinis_MCP"
graph2png(file=paste0("..//benchmarks/plots/high_dimensional/",method,"_",suffix,".png"), 
          width=7, height=5)

# check metrics & quality of MCP solution
system.time(metrics_ordinis_MCP <- data.frame(method=method,
                                      timing_s=timing_ordinis,
                                      do.call(cbind, get.metrics(beta = beta_est0,  # some metrics
                                                                  X = X, 
                                                                  y = y,
                                                                  family = gaussian(identity), 
                                                                  beta_true = beta_true)))) 
metrics_ordinis_MCP



# 1.2B MCP fit using ncvreg (no box constraints support & only canonical links for gaussian, binomial & poisson) ####

library(ncvreg)
gamma_ncvreg = 1.1 # default for MCP is 3 & 3.7 for SCAD, I am using value close to 1 to match best with my later L0 glm objective
timing_ncvreg = system.time(ncvreg_fit <- ncvreg(X, y, 
                                                 # note: intercept is always added
                                                 penalty = "MCP",
                                                 family = "gaussian",
                                                 gamma = gamma_ncvreg, 
                                                 nlambda = 10L, 
                                                 # lower.limits = rep(0, p), # force all coefficients to be positive
                                                 # maxit = 10,
                                                 # maxit.irls = 1,
                                                 dfmax = n-1, # min(n-1, round(2*n/(log(log(n))*log(p)))), # n-1, # 
                                                 alpha = 1 # , # worse solution with 0.9 
                                                 # eps = 1E-4 # good enough to use as initialisation, default 1E-4
))[[3]] 
timing_ncvreg # 21s with default eps, 23.5s with eps=0.1
plot(ncvreg_fit)
max(ncvreg_fit$iter) # 74
median(ncvreg_fit$iter) # 24.5
ncvreg_fit$iter

# calculate ebic, mbic & gic/sic (cf. abess) criterion at each lambda in grid
LLs = logLik(ncvreg_fit) # logLik2(ncvreg_fit, weights=w, X, y)
edfs = attr(LLs, "df")
LLs = as.vector(LLs)
# ebic or gic would identity 5th lambda as best, mbic the 4th but too strict
IC.val = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,IC]
best = which.min(IC.val) 
best # 5th lambda best
secondbest = min(length(IC.val), best+1)
secondbest # 2nd best (slightly less sparse) solution

ncvreg_fit$lambda[best] # best lambda=0.2526485

colSums(ncvreg_fit$beta!=0) # support sizes
colSums(ncvreg_fit$beta!=0)[best] # selected support size: 33
apply(ncvreg_fit$beta, 2, function (col) sum(col[-1]!=0&beta_true!=0)) # TPs

coefs_ncvreg = ncvreg_fit$beta[,best] 
coefs_ncvreg_set = (which(coefs_ncvreg!=0)-1)[-1]
coefs_ncvreg_secondbest = ncvreg_fit$beta[,secondbest] 
coefs_ncvreg = coefs_ncvreg[-which(names(coefs_ncvreg)=="(Intercept)")] # dropping intercept
coefs_ncvreg_secondbest = coefs_ncvreg_secondbest[-which(names(coefs_ncvreg_secondbest)=="(Intercept)")] # dropping intercept

# the max coefficients observed for any of the lambda values works best
# to use as starting coefficients in L0glm
coefs_ncvreg_max = as.numeric(rowMax(ncvreg_fit$beta))[-1] # drop intercept
# coefs_ncvreg_max[2] = 0 # to drop weird first coefficient if needed

min(abs(coefs_ncvreg)[coefs_ncvreg!=0&beta_true!=0]) # 0.2847526 - det threshold
min(abs(coefs_ncvreg_max)[coefs_ncvreg_max!=0&beta_true!=0]) # 0.2404359 - det threshold

sum(coefs_ncvreg!=0) # 32
sum(coefs_ncvreg_max!=0) # 101
L1norm_coefsncvreg = sum(abs(coefs_ncvreg))
L0norm_coefsncvreg = sum(coefs_ncvreg!=0)
L0overL1norm_coefsncvreg = L0norm_betatrue/L1norm_betatrue
L0overL1norm_coefsncvreg # 2.525916

range(coefs_ncvreg) # -0.8277528  0.7738667
range(beta_true) # -0.7984196  0.6614568

table(coefs_ncvreg!=0, beta_true!=0, 
      dnn=c("estimated beta nonzero", "true beta nonzero"))

table(coefs_ncvreg_secondbest!=0, beta_true!=0, 
      dnn=c("estimated beta nonzero", "true beta nonzero"))

table(coefs_ncvreg_max!=0, beta_true!=0, 
      dnn=c("estimated beta nonzero", "true beta nonzero"))

# TODO check how performance would be affected by choice of gamma
# (gamma close to 1 maybe results in objective function closest to that of L0glm?)

cor(coefs_ncvreg, beta_true)              # 0.4486028
cor(coefs_ncvreg_secondbest, beta_true)   # 0.3973878
cor(coefs_ncvreg_max, beta_true)          # 0.319071

fitted = X %*% coefs_ncvreg
cor(y_true, fitted)^2 # population R2=0.5812931

fitted = X %*% coefs_ncvreg_secondbest
cor(y_true, fitted)^2 # population R2=0.7750266

fitted_max = X %*% coefs_ncvreg_max
cor(y_true, fitted_max)^2 # population R2=0.880006

beta_est0 = coefs_ncvreg
plot(x=beta_est0[beta_true!=0&beta_est0!=0], 
     y=beta_true[beta_true!=0&beta_est0!=0], pch=16, col="blue",
     xlab="Estimated coefficients", ylab="True coefficients",
     main="MCP penalized regression (ncvreg, GIC tuned)") # true positives
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
method="ncvreg_MCP"
graph2png(file=paste0("..//benchmarks/plots/high_dimensional/",method,"_",suffix,".png"), 
          width=7, height=5)

# check metrics & quality of MCP solution
system.time(metrics_ncvreg_MCP <- data.frame(method=method,
                                      timing_s=timing_ncvreg,
                                      do.call(cbind, get.metrics(beta = beta_est0,  # some metrics
                                                                 X = X, 
                                                                 y = y,
                                                                 family = gaussian(identity), 
                                                                 beta_true = beta_true)))) 
metrics_ncvreg_MCP



# 1.2C SCAD fit using ncvreg (no box constraints support & only canonical links for gaussian, binomial & poisson) ####

library(ncvreg)
# gamma_ncvreg = 3.7 # default for MCP is 3 & 3.7 for SCAD, I am using value close to 1 to match best with my later L0 glm objective
timing_ncvreg_scad = system.time(ncvreg_scad_fit <- ncvreg(X, y, 
                                                 # note: intercept is always added
                                                 penalty = "SCAD",
                                                 family = "gaussian",
                                                 # gamma = gamma_ncvreg, 
                                                 nlambda = 10L, 
                                                 # lower.limits = rep(0, p), # force all coefficients to be positive
                                                 # maxit = 10,
                                                 # maxit.irls = 1,
                                                 dfmax = n-1, # min(n-1, round(2*n/(log(log(n))*log(p)))), # n-1, # 
                                                 alpha = 1 # , # worse solution with 0.9 
                                                 # eps = 1E-4 # good enough to use as initialisation, default 1E-4
))[[3]] 
timing_ncvreg_scad # 30s with default eps, 23.5s with eps=0.1
plot(ncvreg_scad_fit)
max(ncvreg_scad_fit$iter) # 
median(ncvreg_scad_fit$iter) # 
ncvreg_scad_fit$iter

# calculate ebic, mbic & gic/sic (cf. abess) criterion at each lambda in grid
LLs = logLik(ncvreg_scad_fit) # logLik2(ncvreg_fit, weights=w, X, y)
edfs = attr(LLs, "df")
LLs = as.vector(LLs)
# ebic or gic would identity 5th lambda as best, mbic the 4th but too strict
IC.val = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,IC]
best = which.min(IC.val) 
best # 5th lambda best
secondbest = min(length(IC.val), best+1)
secondbest # 2nd best (slightly less sparse) solution

ncvreg_scad_fit$lambda[best] # best lambda=0.2526485

colSums(ncvreg_scad_fit$beta!=0) # support sizes
colSums(ncvreg_scad_fit$beta!=0)[best] # selected support size: 33
apply(ncvreg_scad_fit$beta, 2, function (col) sum(col[-1]!=0&beta_true!=0)) # TPs

coefs_ncvreg_scad = ncvreg_scad_fit$beta[,best] 
coefs_ncvreg_scad_set = (which(coefs_ncvreg_scad!=0)-1)[-1]
coefs_ncvreg_scad_secondbest = ncvreg_scad_fit$beta[,secondbest] 
coefs_ncvreg_scad = coefs_ncvreg_scad[-which(names(coefs_ncvreg_scad)=="(Intercept)")] # dropping intercept
coefs_ncvreg_scad_secondbest = coefs_ncvreg_scad_secondbest[-which(names(coefs_ncvreg_scad_secondbest)=="(Intercept)")] # dropping intercept

# the max coefficients observed for any of the lambda values works best
# to use as starting coefficients in L0glm
coefs_ncvreg_scad_max = as.numeric(rowMax(ncvreg_scad_fit$beta))[-1] # drop intercept
# coefs_ncvreg_scad_max[2] = 0 # to drop weird first coefficient if needed

# table(coefs_ncvreg_scad!=0, beta_true!=0, 
#       dnn=c("estimated beta nonzero", "true beta nonzero"))

table(coefs_ncvreg_scad_secondbest!=0, beta_true!=0,  # BEST?
      dnn=c("estimated beta nonzero", "true beta nonzero"))

table(coefs_ncvreg_scad_max!=0, beta_true!=0, 
      dnn=c("estimated beta nonzero", "true beta nonzero"))

# TODO check how performance would be affected by choice of gamma
# (gamma close to 1 maybe results in objective function closest to that of L0glm?)

# cor(coefs_ncvreg_scad, beta_true)              # NA
cor(coefs_ncvreg_scad_secondbest, beta_true)   # 0.42
cor(coefs_ncvreg_scad_max, beta_true)          # 0.61

# fitted = X %*% coefs_ncvreg_scad
# cor(y_true, fitted)^2 # population R2=NA

fitted = X %*% coefs_ncvreg_scad_secondbest
cor(y_true, fitted)^2 # population R2=0.36

fitted_max = X %*% coefs_ncvreg_scad_max
cor(y_true, fitted_max)^2 # population R2=0.56

beta_est0 = coefs_ncvreg_scad_secondbest
plot(x=beta_est0[beta_true!=0&beta_est0!=0], 
     y=beta_true[beta_true!=0&beta_est0!=0], pch=16, col="blue",
     xlab="Estimated coefficients", ylab="True coefficients",
     main="MCP penalized regression (ncvreg, GIC tuned)") # true positives
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
method="ncvreg_SCAD"
graph2png(file=paste0("..//benchmarks/plots/high_dimensional/",method,"_",suffix,".png"), 
          width=7, height=5)

# check metrics & quality of MCP solution
system.time(metrics_ncvreg_SCAD <- data.frame(method=method,
                                             timing_s=timing_ncvreg_scad,
                                             do.call(cbind, get.metrics(beta = beta_est0,  # some metrics
                                                                        X = X, 
                                                                        y = y,
                                                                        family = gaussian(identity), 
                                                                        beta_true = beta_true)))) 
metrics_ncvreg_SCAD


# # 1.3 L0glm fit without SIS ####
#
# # without prior MCP variable preselection this doesn't work well
# (returns an empty set)

# # we use mbic because original problem is high dimensional
# lambda_preset = preset.lambda(IC = IC, # options aic, bic, gic, mebic, ebic, mbic
#                               p = ncol(X),
#                               n = nrow(X))
# timing_fullL0glm = system.time(fullL0glm_gaus_presetlambda <- L0glm.fit(
#   X = X,
#   y = y,
#   family = family,
#   lambda = lambda_preset,
#   nonnegative = FALSE,
#   start = rep(1, p),
#   gamma = 2.1,
#   normalize = FALSE,
#   # dispersion = dispersion,
#   solver = "eigen",
#   method = 1L, # Cholesky, TODO CHECK LSCG=4 should be fastest for dense but gives out of memory, either in line Matrix::tcrossprod(D, Xw) in lmridge.R#631 or rbind2(as(Matrix(Xw, sparse = T), "generalMatrix"), Diagonal(x = sqrt(lambdas))) at lmridge.R#722
#   relaxed = FALSE))[[3]]
# timing_fullL0glm
# # 34s
# 
# fullL0glm_gaus_presetlambda$dispersion # 11.91098
# 
# beta_est1 = fullL0glm_gaus_presetlambda$coefficients
# sum(beta_est1!=0&beta_true!=0)  # 0
# table(beta_est1!=0, beta_true!=0,
#       dnn=c("estimated beta nonzero", "true beta nonzero"))
# 
# cor(beta_est1, beta_true) #
# 
# plot(x=beta_est1[beta_true!=0&beta_est1!=0],
#      y=beta_true[beta_true!=0&beta_est1!=0], pch=16, col="blue",
#      xlab="Estimated coefficients", ylab="True coefficients",
#      main="L0glm") # true positives
# min(abs(beta_est1)[beta_true!=0&beta_est1!=0]) #
# abline(a=0, b=1, col="lightgrey")
# points(x=beta_est1[beta_true!=0&beta_est1==0],
#        y=beta_true[beta_true!=0&beta_est1==0],
#        pch=16, col="purple") # false negatives
# points(x=beta_est1[beta_true==0&beta_est1!=0],
#        y=beta_true[beta_true==0&beta_est1!=0],
#        pch=16, col="red") # false positives !!!!!
# points(x=beta_est1[beta_true==0&beta_est1==0][[1]],
#        y=beta_true[beta_true==0&beta_est1==0][[1]],
#        pch=16, col="green") # true negatives
# method="L0glm"
# # graph2png(file=paste0("..//benchmarks/plots/high_dimensional/",method,"_",suffix,".png"),
# #          width=7, height=5)
# 
# # check metrics & quality of solution with preset lambda
# system.time(metrics_fullL0 <- data.frame(method=method,
#                                      timing_s=timing_fullL0glm,
#                                      do.call(cbind, get.metrics(beta = beta_est1,  # some metrics
#                                                                 X = X,
#                                                                 y = y,
#                                                                 family = gaussian(identity),
#                                                                 beta_true = beta_true))))  # 0.5s
# metrics_fullL0
# 
# # TODO check if L0glm can be fixed to give a decent solution out of
# # the box without doing ISIS prescreening? 





# 1.4 abess BEST SUBSET FIT ####

library(abess)
# optionally to speed up things can use SIS-L0glm or MCP preselected set as initial active set
# startcoefs <- beta_SISL0glm
timing_abess = system.time(fit_abess_full <- abess(x = X, 
                                                y = y, 
                                                family = "gaussian",
                                                tune.path = "sequence", 
                                                support.size = 0:min(n, round(4*n/(log(log(n))*log(p))) ), # sum(startcoefs!=0), # 
                                                lambda = 0,
                                                ic.scale = (1/log(n))*preset.lambda(IC=IC, p=p, n=n), # to convert BIC criterion to other IC
                                                # init.active.set = as.integer(which(startcoefs!=0)),
                                                warm.start = TRUE, 
                                                tune.type = "bic"))[[3]]
timing_abess # 72s
plot(fit_abess_full, type="tune")  # optimal support size 
beta_abess = as.vector(extract(fit_abess_full)$beta) 
sum(beta_abess!=0) # 26

plot(x=beta_abess[beta_true!=0&beta_abess!=0], 
     y=beta_true[beta_true!=0&beta_abess!=0], pch=16, col="blue",
     xlab="Estimated coefficients", ylab="True coefficients",
     main=paste0("best subset (abess, GIC tuned)")) # true positives
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
method="abess"
graph2png(file=paste0("..//benchmarks/plots/high_dimensional/",method,"_",suffix,".png"), 
          width=7, height=5)

table(beta_abess!=0, beta_true!=0, 
      dnn=c("estimated beta nonzero", "true beta nonzero"))

# check metrics & quality of solution 
system.time(metrics_abess <- data.frame(method="abess",
                                             timing_s=timing_abess,
                                             do.call(cbind, get.metrics(beta = beta_abess,  # some metrics
                                                                   X = X, 
                                                                   y = y,
                                                                   family = gaussian(identity), 
                                                                   beta_true = beta_true))))  # 0.5s
metrics_abess



# 1.5 L0Learn FIT ####
library(L0Learn)

# startcoefs = coefs_ncvreg # MCP starting coefficients
# selcoefs = which(startcoefs!=0)
# p_selcoefs = length(selcoefs)
W = rep(1, nrow(X))

timing_L0Learn = system.time(fit_L0Learn <- L0Learn.fit(x=as.matrix(Diagonal(x=sqrt(W)) %*% X), # [,selcoefs,drop=F]),  # note: no support for sparse matrices
                                                  y=as.vector(y*sqrt(W)),
                                                  loss = "SquaredError",
                                                  penalty="L0", algorithm="CD",
                                                  maxSuppSize=(n-1), 
                                                  intercept=FALSE, 
                                                  nLambda=100, # NOTE: 100 gave much better results than 20
                                                  # activeSetNum = 3, maxSwaps = 100,
                                                  autoLambda=TRUE
                                                  # lambdaGrid=list(10^seq(-5,5,length.out=nlam)/ # DECREASING sequence recommended, but INCREASING works better under nonnegativity, see https://cran.r-project.org/web/packages/L0Learn/vignettes/L0Learn-vignette.html
                                                  #                   (norm(sim$y*sqrt(1/(sim$y+1)),"2")^2)) 
))[[3]]  # lambda of 1 is scaled by 1/(norm(weighted y,"2")) because L0Learn scales y by L2 norm of y
timing_L0Learn # 28s
# timing_L0Learn = timing_L0Learn + timing_ordinis
plot(fit_L0Learn)

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

LLs = logLik(fit_L0Learn, weights=W, X, y) # [,selcoefs,drop=F]
edfs = attr(LLs, "df")
LLs = as.vector(LLs)
ICs = get.IC(loglik = LLs, n = n, p = fit_L0Learn$p, edf = edfs)[,IC]
beta_L0Learn = as.numeric(coef(fit_L0Learn,
                                lambda = fit_L0Learn$lambda[[1]][which.min(ICs)],
                                gamma = 0))
# beta_L0Learn_after_MCP = matrix(0, p)
# beta_L0Learn_after_MCP[selcoefs] = beta_L0Learn_selcoefs
# beta_L0Learn_after_MCP = as.vector(beta_L0Learn_after_MCP)
sum(beta_L0Learn!=0) # 14

plot(x=beta_L0Learn[beta_true!=0&beta_L0Learn!=0], 
     y=beta_true[beta_true!=0&beta_L0Learn!=0], pch=16, col="blue",
     xlab="Estimated coefficients", ylab="True coefficients",
     main=paste0("L0Learn (GIC tuned)")) # true positives
abline(a=0, b=1, col="lightgrey")
points(x=beta_L0Learn[beta_true!=0&beta_L0Learn==0], 
       y=beta_true[beta_true!=0&beta_L0Learn==0], 
       pch=16, col="purple") # false negatives
points(x=beta_L0Learn[beta_true==0&beta_L0Learn!=0], 
       y=beta_true[beta_true==0&beta_L0Learn!=0], 
       pch=16, col="red") # false positives
points(x=beta_L0Learn[beta_true==0&beta_L0Learn==0][[1]], 
       y=beta_true[beta_true==0&beta_L0Learn==0][[1]], 
       pch=16, col="green") # true negatives
method="L0Learn"
graph2png(file=paste0("..//benchmarks/plots/high_dimensional/",method,"_",suffix,".png"), 
          width=7, height=5)

table(beta_L0Learn!=0, beta_true!=0, 
      dnn=c("estimated beta nonzero", "true beta nonzero"))


# check metrics & quality of solution with preset lambda
system.time(metrics_L0Learn <- data.frame(method="L0Learn",
                                          timing_s=timing_L0Learn,
                                          do.call(cbind, get.metrics(beta = beta_L0Learn,  # some metrics
                                                                X = X, 
                                                                y = y,
                                                                family = gaussian(identity), 
                                                                beta_true = beta_true))))  # 0.5s
metrics_L0Learn



# 1.6 LASSO FIT using glmnet ####

library(glmnet)
timing_cvglmnet = system.time(glmnet_cvfit <- cv.glmnet(x = X, y = y, 
                                                        nfolds = 5,
                                                        intercept=F, 
                                                        nlambda=10,
                                                        relaxed=F,
                                                        alpha=1))[[3]] 
timing_cvglmnet # 95s
plot(glmnet_cvfit)
bestlambdaCV1 <- glmnet_cvfit$lambda.min #
bestlambdaCV1 # 0.3441614
bestlambdaCV2 <- glmnet_cvfit$lambda.1se #
bestlambdaCV2 # 0.3441614
timing_glmnetfull = system.time(glmnet_fit <- glmnet(x = X, y = y, 
                                                     intercept = F, 
                                                     nlambda=10,
                                                     relaxed=F,
                                                     alpha = 1))[[3]] 
timing_glmnetfull # 8.6s
timing_glmnet = timing_cvglmnet+timing_glmnetfull
timing_glmnet # 105s
plot(glmnet_fit)

# TPs for LASSO coefficients obtained with different lambda
# 6th best?
apply(coef(glmnet_fit)[-1,], 2, function (coefs) length(intersect(which(coefs!=0), beta_true_set)))
apply(coef(glmnet_fit)[-1,], 2, function (coefs) sum(coefs!=0)) # support size

# # tune based on IC calculated from unpenalised GLM refitted LASSO selected coefs
# tab  <- score_ic_glmnet(glmnet_fit, X, y,
#                         family = gaussian(), 
#                         family_name = "gaussian",
#                         weights = rep(1, nrow(X)), offset = NULL,
#                         intercept = FALSE,
#                         refit_method = "bigGlm",
#                         bigGlm_path = TRUE,                 # warm-start if k ~ n
#                         lower.limits = -Inf, upper.limits = Inf,
#                         p_total = ncol(X),
#                         gamma_ebic = max(0, 1 - log(nrow(X))/(2*log(ncol(X)))),                    # or omit to use CC lower bound
#                         Cn_hbic = log(log(nrow(X))), 
#                         min_resid_df = 3L # prevent saturation
#                         )                 
# tab
# # IC-selected lambda and refitted coefficients
# bestIC    <- which.min(tab$bic)
# bestIC # 6 # BIC gives set 6, EBIC would return 1 (empty set)
# lam  <- tab$lambda[bestIC]
# S    <- tab$active[[bestIC]]
# length(intersect(S, beta_true_set)) # 23
# beta_refit <- as.vector(tab$coef_full[[bestIC]])  # refitted unpenalised coefs 
# # coefficients on S (slope-only)
# coefs_glmnetIC = as.vector(coef(glmnet_fit, 
#                                 s = lam)[-1])
# sum(coefs_glmnetIC!=0) # 432
# length(intersect(which(coefs_glmnetIC!=0), beta_true_set)) # 23
# table(coefs_glmnetIC!=0, beta_true!=0,  
#       dnn=c("estimated beta nonzero", "true beta nonzero"))

coefs_glmnetCV1 <- as.vector(coef(glmnet_fit, s=bestlambdaCV1)[-1])
sum(coefs_glmnetCV1!=0) # 135
length(intersect(which(coefs_glmnetCV1!=0), beta_true_set)) # 20
table(coefs_glmnetCV1!=0, beta_true!=0,  
      dnn=c("estimated beta nonzero", "true beta nonzero"))


beta_LASSO = coefs_glmnetCV1
plot(x=beta_LASSO[beta_true!=0&beta_LASSO!=0], 
     y=beta_true[beta_true!=0&beta_LASSO!=0], pch=16, col="blue",
     xlab="Estimated coefficients", ylab="True coefficients",
     main=paste0("LASSO (glmnet, 5-fold CV tuned)")) # true positives
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
method="glmnet_LASSO_CV"
graph2png(file=paste0("..//benchmarks/plots/high_dimensional/",method,"_",suffix,".png"), 
          width=7, height=5)

# check metrics & quality of LASSO solution tuned via CV
system.time(metrics_LASSO_CV <- data.frame(method=method,
                                        timing_s=timing_glmnet,
                                        do.call(cbind, get.metrics(beta = coefs_glmnetCV1,  # some metrics
                                                                   X = X, 
                                                                   y = y,
                                                                   family = gaussian(identity), 
                                                                   beta_true = beta_true)))) 
metrics_LASSO_CV

# note: the cross validation I might not really need, as I can just select a model
# with max n variables, but performance worse than MCP
 
# note: adaptive LASSO fit also doesn't help to select more variables



# 1.7 (fully) relaxed LASSO fit using glmnet ####
# see https://cran.r-project.org/web/packages/glmnet/vignettes/relax.pdf

library(glmnet)
timing_cvglmnetrelaxed = system.time(glmnet_cvfitrelaxed <- cv.glmnet(x = X, y = y, 
                                                                      nfolds = 5,
                                                        intercept=F, 
                                                        alpha = 1,
                                                        nlambda = 10,
                                                        gamma = 0, # fully relaxed fit
                                                        relax = T # relaxed LASSO
                                                        ))[[3]] 
timing_cvglmnetrelaxed # 247s
plot(glmnet_cvfitrelaxed)
bestlambdaCV1 <- glmnet_cvfitrelaxed$relaxed$lambda.min #
bestlambdaCV1 # 0.5740958
bestlambdaCV2 <- glmnet_cvfitrelaxed$relaxed$lambda.1se #
bestlambdaCV2 # 0.5740958
# bestgamma <- glmnet_cvfitrelaxed$relaxed$gamma.min
# bestgamma 

timing_glmnetfullrelaxed = system.time(glmnet_fitrelaxed <- glmnet(x = X, y = y, 
                                                            intercept = F, 
                                                            alpha = 1,
                                                            gamma = 0,
                                                            nlambda = 10,
                                                            relax = T,
                                                            maxp = n-1))[[3]] # 
timing_glmnetfullrelaxed # 26s
timing_glmnetrelaxed_CV = timing_cvglmnetrelaxed+timing_glmnetfullrelaxed
timing_glmnetrelaxed_BIC = timing_glmnetfullrelaxed

plot(glmnet_fitrelaxed$relaxed)
abline(v=log(bestlambdaCV1), col="red")

# coefficients obtained by 5-fold CV tuning
coefs_glmnetrelaxed_CV <- as.vector(coef(glmnet_fitrelaxed$relaxed, 
                                         s=bestlambdaCV1, gamma=0)[-1])
sum(coefs_glmnetrelaxed_CV!=0) # 26
glmnet_set1relaxed <- which(coefs_glmnetrelaxed_CV!=0)
length(intersect(glmnet_set1relaxed, beta_true_set)) # 13
table(coefs_glmnetrelaxed_CV!=0, beta_true!=0,  
      dnn=c("estimated beta nonzero", "true beta nonzero"))

# coefficients obtained by BIC tuning
LLs <- L0glm::logLik.glmnet(glmnet_fitrelaxed$relaxed, weights=w, X=X, y=y)
edfs = attr(LLs, "df")
# aic gives same result here, all other IC return intercept only model
IC.val = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,"bic"]
IC.val
bestBIC = which.min(IC.val)
bestBIC # 4
bestlambda_BIC = glmnet_fitrelaxed$relaxed$lambda[bestBIC] 
bestlambda_BIC # 0.2063193
coefs_glmnetrelaxed_BIC = as.vector(coef(glmnet_fitrelaxed$relaxed, 
                                         s = bestlambda_BIC, gamma=0)[-1])
table(coefs_glmnetrelaxed_BIC!=0, beta_true!=0,  
      dnn=c("estimated beta nonzero", "true beta nonzero"))


# CONCLUSION: here BIC yield more TPs but also very large nr of FPs

# plots & metrics for coefficients tuned via CV
if (sum(coefs_glmnetrelaxed_CV!=0)!=0) {
beta_LASSOrelaxedCV = coefs_glmnetrelaxed_CV
plot(x=beta_LASSOrelaxedCV[beta_true!=0&beta_LASSOrelaxedCV!=0], 
     y=beta_true[beta_true!=0&beta_LASSOrelaxedCV!=0], pch=16, col="blue",
     xlab="Estimated coefficients", ylab="True coefficients",
     main=paste0("relaxed LASSO (glmnet, 5-fold CV tuned)")) # true positives
abline(a=0, b=1, col="lightgrey")
points(x=beta_LASSOrelaxedCV[beta_true!=0&beta_LASSOrelaxedCV==0], 
       y=beta_true[beta_true!=0&beta_LASSOrelaxedCV==0], 
       pch=16, col="purple") # false negatives
points(x=beta_LASSOrelaxedCV[beta_true==0&beta_LASSOrelaxedCV!=0], 
       y=beta_true[beta_true==0&beta_LASSOrelaxedCV!=0], 
       pch=16, col="red") # false positives
points(x=beta_LASSOrelaxedCV[beta_true==0&beta_LASSOrelaxedCV==0][[1]], 
       y=beta_true[beta_true==0&beta_LASSOrelaxedCV==0][[1]], 
       pch=16, col="green") # true negatives
method="glmnet_relaxed_LASSO_CV"
graph2png(file=paste0("..//benchmarks/plots/high_dimensional/",method,"_",suffix,".png"), 
          width=7, height=5)

# check metrics & quality of LASSO solution
system.time(metrics_relaxedLASSO_CV <- data.frame(method=method,
                                        timing_s=timing_glmnetrelaxed_CV,
                                        do.call(cbind, get.metrics(beta = beta_LASSOrelaxedCV,  # some metrics
                                                                   X = X, 
                                                                   y = y,
                                                                   family = gaussian(identity), 
                                                                   beta_true = beta_true)))) 
metrics_relaxedLASSO_CV
}

# plots & metrics for coefficients tuned based on BIC
beta_LASSOrelaxedBIC = coefs_glmnetrelaxed_BIC
plot(x=beta_LASSOrelaxedBIC[beta_true!=0&beta_LASSOrelaxedBIC!=0], 
     y=beta_true[beta_true!=0&beta_LASSOrelaxedBIC!=0], pch=16, col="blue",
     xlab="Estimated coefficients", ylab="True coefficients",
     main=paste0("relaxed LASSO (glmnet, BIC tuned)")) # true positives
abline(a=0, b=1, col="lightgrey")
points(x=beta_LASSOrelaxedBIC[beta_true!=0&beta_LASSOrelaxedBIC==0], 
       y=beta_true[beta_true!=0&beta_LASSOrelaxedBIC==0], 
       pch=16, col="purple") # false negatives
points(x=beta_LASSOrelaxedBIC[beta_true==0&beta_LASSOrelaxedBIC!=0], 
       y=beta_true[beta_true==0&beta_LASSOrelaxedBIC!=0], 
       pch=16, col="red") # false positives
points(x=beta_LASSOrelaxedBIC[beta_true==0&beta_LASSOrelaxedBIC==0][[1]], 
       y=beta_true[beta_true==0&beta_LASSOrelaxedBIC==0][[1]], 
       pch=16, col="green") # true negatives
method="glmnet_relaxed_LASSO_BIC"
graph2png(file=paste0("..//benchmarks/plots/high_dimensional/",method,"_",suffix,".png"), 
          width=7, height=5)

# check metrics & quality of LASSO solution
system.time(metrics_relaxedLASSO_BIC <- data.frame(method=method,
                                                  timing_s=timing_glmnetrelaxed_BIC,
                                                  do.call(cbind, get.metrics(beta = beta_LASSOrelaxedBIC,  # some metrics
                                                                             X = X, 
                                                                             y = y,
                                                                             family = gaussian(identity), 
                                                                             beta_true = beta_true)))) 
metrics_relaxedLASSO_BIC



# 1.8 COMPARISON OF METRICS ####

benchmarks = bind_rows(metrics_LASSO_CV, 
                       # metrics_relaxedLASSO_CV,
                       metrics_relaxedLASSO_BIC,
                       metrics_ncvreg_MCP,
                       metrics_ncvreg_SCAD,
                       metrics_ordinis_MCP,
                       metrics_ISISL0glm,
                       metrics_abess,
                       metrics_L0Learn)
benchmarks = benchmarks[order(benchmarks$deviance_true, decreasing=F),]
benchmarks
#                    method timing_s       p   n m  deviance deviance_true     loglik
# 6               ISIS_L0glm    28.44 1000001 500 1  783.5221      365.5327  -822.7671
# 2 glmnet_relaxed_LASSO_BIC    23.33 1000001 500 1  267.1526      673.5117  -553.7723
# 7                    abess    60.59 1000001 500 1 1919.2327     1714.5316 -1046.7374
# 1          glmnet_LASSO_CV    61.13 1000001 500 1 2558.7304     2193.3494 -1118.6339
# 8                  L0Learn    28.28 1000001 500 1 3505.9349     2941.0004 -1197.3704
# 3               ncvreg_MCP    19.84 1000001 500 1 3380.5188     2987.8866 -1188.2634
# 5              ordinis_MCP    42.69 1000001 500 1 5339.3327     4574.5681 -1302.5312
# 4              ncvreg_SCAD    30.07 1000001 500 1 6904.1217     6051.6176 -1366.7857
#       MD_coefs    MAD_coefs  sMAPE_coefs bias_coefs   R2_coefs  MD_selcoefs MAD_selcoefs
# 6 1.470970e-06 5.539483e-06 4.851911e-05  1.0063953 0.93694967  0.002829568   0.05567598
# 2 5.417878e-06 3.789864e-05 4.387157e-04  0.8599231 0.19630980  0.008861523   0.13481324
# 7 5.898213e-06 1.828169e-05 9.898516e-05  0.7589768 0.49834459  0.132191790   0.21878243
# 1 6.316097e-06 2.810534e-05 4.436671e-04  1.7913683 0.28414221  0.013639284   0.08272126
# 8 3.523403e-06 1.877065e-05 9.852560e-05  0.7904977 0.42087555 -0.017563061   0.18551307
# 3 3.809486e-06 2.452741e-05 1.173056e-04  0.5624954 0.20124643 -0.091567537   0.34896751
# 5 5.155387e-06 2.461274e-05 1.183236e-04  0.5663658 0.09128508 -0.120785775   0.34934859
# 4 7.846708e-06 2.436942e-05 1.275653e-04  3.6792121 0.18717472  0.293897069   0.32130534
#   sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true support_size_estimated
# 6      0.1189240     1.0063953   0.9872315                60                     38
# 2      1.8974241     0.8599231   0.3215343                60                    188
# 7      0.6532755     0.7589768   0.8209510                60                     26
# 1      1.9237609     1.7913683   0.4653944                60                    188
# 8      0.3232571     0.7904977   0.9381563                60                     14
# 3      0.9565987     0.5624954   0.6508893                60                     16
# 5      0.9033776     0.5663658   0.5778110                60                      7
# 4      1.6588675     3.6792121   0.7404923                60                     13
#   estimated_over_true_support_size TP  FP     TN FN precision     recall specificity
# 6                        0.6333333 38   0 999940 22 1.0000000 0.63333333    1.000000
# 2                        3.1333333 19 169 999771 41 0.1010638 0.31666667    0.999831
# 7                        0.4333333 19   7 999933 41 0.7307692 0.31666667    0.999993
# 1                        3.1333333 19 169 999771 41 0.1010638 0.31666667    0.999831
# 8                        0.2333333 13   1 999939 47 0.9285714 0.21666667    0.999999
# 3                        0.2666667  9   7 999933 51 0.5625000 0.15000000    0.999993
# 5                        0.1166667  4   3 999937 56 0.5714286 0.06666667    0.999997
# 4                        0.2166667  7   6 999934 53 0.5384615 0.11666667    0.999994
#   accuracy  F1_score      gic     ebic      bic      aic     mbic
# 6 0.999978 0.7755102 2629.878 2723.144 1887.904 1723.534 2857.383
# 2 0.999790 0.1532258 5877.828 6329.808 2282.105 1485.545 6980.350
# 7 0.999952 0.4418605 2774.944 2839.512 2261.269 2147.475 2932.447
# 1 0.999790 0.1532258 7007.551 7459.531 3411.829 2615.268 8110.073
# 8 0.999952 0.3513514 2773.335 2809.206 2487.960 2424.741 2860.836
# 3 0.999942 0.2368421 2805.600 2846.254 2482.175 2410.527 2904.768
# 5 0.999941 0.1194030 2806.979 2826.111 2654.779 2621.062 2853.647
# 4 0.999941 0.1917808 3086.926 3120.406 2820.576 2761.571 3168.594
# CONCLUSION: ISIS_L0glm clear winner

benchmarks[benchmarks$method %in% c("abess", "ncvreg_MCP", "ncvreg_SCAD",
                                    "ISIS_L0glm", "glmnet_relaxed_LASSO_BIC", "glmnet_LASSO_CV"),
           c("method","timing_s","deviance_true", "TP", "FP", "TN", "FN", "F1_score", "gic")]
#                     method timing_s deviance_true TP  FP     TN FN  F1_score      gic
# 6               ISIS_L0glm    28.44      365.5327 38   0 999940 22 0.7755102 2629.878
# 2 glmnet_relaxed_LASSO_BIC    23.33      673.5117 19 169 999771 41 0.1532258 5877.828
# 7                    abess    60.59     1714.5316 19   7 999933 41 0.4418605 2774.944
# 1          glmnet_LASSO_CV    61.13     2193.3494 19 169 999771 41 0.1532258 7007.551
# 3               ncvreg_MCP    19.84     2987.8866  9   7 999933 51 0.2368421 2805.600
# 4              ncvreg_SCAD    30.07     6051.6176  7   6 999934 53 0.1917808 3086.926

# TODO: add most promising methods mentioned in https://onlinelibrary.wiley.com/doi/10.1111/insr.12469
# e.g. from https://cran.r-project.org/web/packages/flare/index.html &
# https://cran.r-project.org/web/packages/scalreg/index.html & 
# https://cran.r-project.org/web/packages/fda.usc/index.html
# e.g. AdapL.1se (adaptive lasso with 1SE rule), Dant (flare::slim with
# method = "lq"), DC.VS (fda.usc::fregre.glm.vs)
# see https://github.com/LauraFreiG/Covariates_selection/blob/main/Linear%20Regression/Real%20data%20sets/script_real_data.R
# implement ISIS method in L0glm or test if cutoff on coefficients in
# L0glm.fit, e.g. after 3d or 4th iteration, could already lead to good
# performance out of the box

# TODO : remove benchmarks related to bias on all coefficients, as measures of
# bias on selected variables are then confounded with the efficiency of variable
# selection / support recovery?



#### TEST CODE NOT CURRENTLY USED ####

# XX TEST OF SOME POSSIBLE UNIVARIATE REGRESSION STARTING VALUES ####

# these were all pretty bad - so not worth to do this

# univariate zero intercept least square regressions as starting values ####
# no good...
system.time(coefs_univarLS <- apply(X, 2, function (colX) sum(w*colX*y)/sum(w*colX^2))) # 12s
# NOTE: should be done in Rcpp & multithreaded using OpenMP
beta_univarLS = coefs_univarLS
indices_n_largest <- sort(abs(beta_univarLS), decreasing=T, index.return=T)$ix[1:n]
beta_univarLS2 <- rep(0, p)
beta_univarLS2[indices_n_largest] <- beta_univarLS[indices_n_largest]
# beta_univarLS <- beta_univarLS2
cor(beta_univarLS2, beta_true) # 0.22
plot(beta_univarLS2, beta_true)
sum(beta_univarLS2!=0) # =n
sum(beta_univarLS2!=0&beta_true!=0) # 17
sum(beta_univarLS2!=0&beta_true==0) # 483
table(beta_univarLS2!=0, beta_true!=0, 
      dnn=c("estimated beta nonzero", "true beta nonzero"))

# univariate zero intercept weighted quantile (median) regressions as starting values ####
# even worse
library(robsurvey)
system.time(coefs_univarwquant <- apply(X, 2, function (colX) weighted_quantile(x=y/colX, 
                                                                                w=abs(w * colX), 
                                                                                probs=0.5))) # 37s
# should be implemented in C++ & multithreaded using OpenMP etc
beta_univarwquant = coefs_univarwquant
indices_n_largest <- sort(abs(beta_univarwquant), decreasing=T, index.return=T)$ix[1:n]
beta_univarwquant2 <- rep(0, p)
beta_univarwquant2[indices_n_largest] <- beta_univarwquant[indices_n_largest]
# beta_univarwquant <- beta_univarwquant2
cor(beta_univarwquant2, beta_true) # 0.19
plot(beta_univarwquant2, beta_true)
sum(beta_univarwquant2!=0)
sum(beta_univarwquant2!=0&beta_true!=0) # 15
sum(beta_univarwquant2!=0&beta_true==0) # 485
table(beta_univarwquant2!=0, beta_true!=0, 
      dnn=c("estimated beta nonzero", "true beta nonzero"))



# XX get bagged MCP coefs ####
# currently not used so marked out

library(ordinis)
gamma_mcp = 1.1 # default for MCP is 1.4, I am using value close to 1 to match best with my later L0 glm objective
nboot = 5
coefs_ordinis_boot = Matrix(0, ncol=ncol(X), nrow=nboot, sparse=T)
for (i in 1:nboot) {
sel = sample(1:nrow(X), replace=T)
ordinis_fit = ordinis(X[sel,,drop=F], y[sel],
                                     intercept=F,
                                   penalty = "mcp",
                                   gamma = gamma_mcp,
                                   # lower.limits = rep(0, p), # force all coefficients to be positive
                                   nlambda = 100L, # 100L,
                                   # lambda.min.ratio = 1E-10,
                                   dfmax = n,
                                   alpha = 1 #,
                                   # maxit = 1,
                                   # tol = 0.1
)
# 7s with maxit 1 (with nlambda 1), 60s with defaults (with nlambda 100)
coefs_ordinis_b = ordinis_fit$beta[,which.max(logLik(ordinis_fit))][-1] # gives same result as AIC & BIC here
sum(coefs_ordinis_b!=0) #
sum(coefs_ordinis_b!=0&beta_true!=0) # true pos
sum(coefs_ordinis_b!=0&beta_true==0) # false pos

beta_est0_b = coefs_ordinis_b

coefs_ordinis_boot[i,] = Matrix(beta_est0_b, sparse=T)
# table(beta_est0!=0, beta_true!=0,
#      dnn=c("estimated beta nonzero", "true beta nonzero"))
}

# avg coefficients over all bootstrap replicates
coefs_ordinis_boot_avg = colMeans(rbind(coefs_ordinis, coefs_ordinis_boot))
sum(coefs_ordinis_boot_avg!=0) # 300 # best to have this < n?
sum(coefs_ordinis_boot_avg!=0&beta_true!=0) # 44 true pos - nr of true pos did not increase by using bagged MCP estimates...
sum(coefs_ordinis_boot_avg!=0&beta_true==0) # 256 false positives still

beta_est0_boot = coefs_ordinis_boot_avg

table(beta_est0_boot!=0, beta_true!=0,
      dnn=c("estimated beta nonzero", "true beta nonzero"))

# so nr of true pos did not increase by using bagged MCP estimates instead of original MCP estimate

# check metrics & quality of MCP solution
system.time(metrics_MCP_boot <- data.frame(method="ordinis_MCP",
                                      timing_s=60,
                                      do.call(cbind, get.metrics(beta = beta_est0_boot,  # some metrics
                                                                 X = X,
                                                                 y = y,
                                                                 family = gaussian(identity),
                                                                 beta_true = beta_true))))
metrics_MCP_boot


# XX iterative adaptive MCP fit ####

# attempt to use iterative adaptive MCP regression to approximate
# L0 penalized regression & penalize the small coefficients to zero

startcoefs = beta_est0_boot # beta_est0
maxcoef = max(abs(beta_est0[beta_est0!=0]))
delta = maxcoef/10 # use in penfact, decreased in later iterations
percentile = 0.1

plotit = FALSE

system.time(
  {
    for (iter in (1:40)) {
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
                                               nlambda = 1L,
                                               # lambda.min.ratio = 1E-100,
                                               dfmax = n,
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
      
      beta_est0B = coefs_ordinis_adMCP
      
      table(beta_est0B!=0, beta_true!=0, 
            dnn=c("estimated beta nonzero", "true beta nonzero"))
      #                          true beta nonzero
      # estimated beta nonzero   FALSE   TRUE
      #                   FALSE  999949     10
      #                    TRUE       1     40
      
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

cor(beta_est0B, beta_true) # 0.9863264

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
       pch=16, col="red") # false positives # TODO check where this false pos is coming from
points(x=beta_est0B[beta_true==0&beta_est0B==0][[1]], 
       y=beta_true[beta_true==0&beta_est0B==0][[1]], 
       pch=16, col="green") # true negatives
graph2png(file="..//benchmarks/plots/high_dimensional/2_adaptive_MCP_regression_p_1_million_n_500_50 nonzero.png", 
          width=7, height=5)

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
system.time(metrics_adMCP <- data.frame(method="iterative_adaptive_MCP",
                                        timing_s=10+14,
                                        do.call(cbind, get.metrics(beta = beta_est0B,  # some metrics
                                                                   X = X, 
                                                                   y = y,
                                                                   family = gaussian(identity), 
                                                                   beta_true = beta_true)))) 
metrics_adMCP
# method timing_s       p   n m deviance deviance_true   loglik      MD_coefs
# 1 iterative_adaptive_MCP       24 1000001 500 1 21418.35     0.2038977 998310.2 -1.750374e-07
# MAD_coefs  sMAPE_coefs bias_coefs  R2_coefs  MD_selcoefs MAD_selcoefs sMAPE_selcoefs
# 1 2.968202e-06 2.618622e-05  0.9880304 0.9932962 -0.003334165   0.05345283      0.1508834
# bias_selcoefs R2_selcoefs support_size_true support_size_estimated
# 1     0.9880304   0.9951838                50                     41
# estimated_over_true_support_size TP FP     TN FN precision recall specificity accuracy  F1_score
# 1                             0.82 40  1 999949 10 0.9756098    0.8    0.999999 0.999989 0.8791209
# gic     ebic      bic      aic     mbic
# 1 -1995560 -1995460 -1996359 -1996536 -1995315


# XX FULL L0Learn fit on all variables ####
# library(L0Learn)
# startcoefs = beta_est0 # MCP starting coefficients
# 
# selcoefs = 1:p # which(startcoefs!=0)
# p_selcoefs = length(selcoefs)
# W = rep(1, nrow(X))
# 
# system.time(fit_L0Learn_full <- L0Learn.fit(x=as.matrix(Diagonal(x=sqrt(W)) %*% X[,selcoefs,drop=F]),  # note: no support for sparse matrices
#                                                 y=as.vector(y*sqrt(W)),
#                                                 loss = "SquaredError",
#                                                 penalty="L0", algorithm="CD",
#                                                 maxSuppSize=(p_selcoefs-1), 
#                                                 intercept=FALSE, 
#                                                 nLambda=100, # NOTE: 100 gave much better results than 20
#                                                 # activeSetNum = 3, maxSwaps = 100,
#                                                 autoLambda=TRUE
#                                                 # lambdaGrid=list(10^seq(-5,5,length.out=nlam)/ # DECREASING sequence recommended, but INCREASING works better under nonnegativity, see https://cran.r-project.org/web/packages/L0Learn/vignettes/L0Learn-vignette.html
#                                                 #                   (norm(sim$y*sqrt(1/(sim$y+1)),"2")^2)) 
# ))  # lambda of 1 is scaled by 1/(norm(weighted y,"2")) because L0Learn scales y by L2 norm of y
# # 54s
# # plot(fit_L0Learn_full)
# 
# # log-likelihood method for gaussian L0Learn fit,
# # cf logLik.ordinis in https://github.com/jaredhuling/ordinis/blob/master/R/ordinis_methods.R
# logLik.L0Learn = function(object, weights, X, y) {
#   n  = object$n
#   p = object$p
#   df = object$suppSize[[1]]
#   # df = nr of nonzero coefficients (incl intercept)
#   
#   # Compute y_pred for all columns of object$beta at once
#   y_pred = X %*% object$beta[[1]]
#   
#   # Compute squared residuals, multiply by weights, and sum over rows to get resid.ss
#   resid.ss = colSums(weights * (matrix(rep(y, ncol(y_pred)), ncol = ncol(y_pred) ) -
#                                   y_pred) ^ 2)
#   
#   rdf = n - df
#   logL = -0.5 * n * (log(2 * pi) - log(rdf) + log(resid.ss)) - 0.5 * rdf
#   attr(logL, "df") = df  # nr of nonzero coefficients
#   attr(logL, "nobs") = n
#   attr(logL, "p") = object$dim[[2]]
#   return(logL)
# }
# 
# LLs = logLik(fit_L0Learn_full, weights=W, X[,selcoefs,drop=F], y)
# edfs = attr(LLs, "df")
# LLs = as.vector(LLs)
# crit = "mbic"
# IC = get.IC(loglik = LLs, n = n, p = p_selcoefs, edf = edfs)[,crit]
# beta_L0Learn_full = coef(fit_L0Learn_full,
#                              lambda = fit_L0Learn_full$lambda[[1]][which.min(IC)],
#                              gamma = 0)
# beta_L0Learn_full = as.vector(beta_L0Learn_full)
# sum(beta_L0Learn_full!=0) # 33
# 
# plot(x=beta_L0Learn_full[beta_true!=0&beta_L0Learn_full!=0], 
#      y=beta_true[beta_true!=0&beta_L0Learn_full!=0], pch=16, col="blue",
#      xlab="Estimated coefficients", ylab="True coefficients",
#      main=paste0("L0Learn")) # true positives
# abline(a=0, b=1, col="lightgrey")
# points(x=beta_L0Learn_full[beta_true!=0&beta_L0Learn_full==0], 
#        y=beta_true[beta_true!=0&beta_L0Learn_full==0], 
#        pch=16, col="purple") # false negatives
# points(x=beta_L0Learn_full[beta_true==0&beta_L0Learn_full!=0], 
#        y=beta_true[beta_true==0&beta_L0Learn_full!=0], 
#        pch=16, col="red") # false positives
# points(x=beta_L0Learn_full[beta_true==0&beta_L0Learn_full==0][[1]], 
#        y=beta_true[beta_true==0&beta_L0Learn_full==0][[1]], 
#        pch=16, col="green") # true negatives
# graph2png(file="..//benchmarks/plots/high_dimensional/6b_L0Learn_full_vars_p_1_million_n_500_50 nonzero.png", 
#           width=7, height=5)
# table(beta_L0Learn_full!=0, beta_true!=0, 
#       dnn=c("estimated beta nonzero", "true beta nonzero"))
# #                         true beta nonzero
# # estimated beta nonzero  FALSE   TRUE
# #                  FALSE 999950     17
# #                   TRUE      0     33
# 
# 
# # check metrics & quality of solution with preset lambda
# system.time(metrics_L0Learn_full <- data.frame(method="L0Learn",
#                                           timing_s=54,
#                                           do.call(cbind, get.metrics(beta = beta_L0Learn_full,  # some metrics
#                                                                      X = X, 
#                                                                      y = y,
#                                                                      family = gaussian(identity), 
#                                                                      beta_true = beta_true))))  # 0.5s
# metrics_L0Learn_full


