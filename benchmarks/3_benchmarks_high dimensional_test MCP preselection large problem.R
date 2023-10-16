# TEST HIGH DIMENSIONAL VARIABLE SELECTION USING L0GLM ####
# here I use MCP as a pre-screen on a large problem and then fit an L0glm on the 
# pre-selected subset to remove some false positives & also test an adaptive L0glm algorithm

# fix to set working package repository when using Microsoft R Open with Intel MKL
options(repos = c(CRAN="https://cloud.r-project.org/"))
library(remotes)
# remotes::install_github("jaredhuling/ordinis")
# remotes::install_github("tomwenseleers/L0glm/L0glm")
# library(L0glm)
devtools::load_all(".") # if you forked & cloned L0glm first to your laptop
library(export)
# source("..//benchmarks/generate data.R")
library(Matrix)


# 1. TEST HIGH DIMENSIONAL PROBLEM ####

n = 500 # 500 obs
p = 1000000 # 1 million variables
s = 0.00008 # prop of nonnegative coefficients
n_nonzero = round(p*s)
n_nonzero # 80
suffix = paste0("n_",as.integer(n),"_p_",as.integer(p),"_n_nonzero",n_nonzero)
beta_sd = 1 # this would directly correspond to Cohen's d - 0.2, 0.5 or 0.8 would be considered small, medium or large
error_sd = 1 
dispersion = error_sd^2

data_gaus  = simulate_standardnormal(n=n, p=p, n_nonzero=n_nonzero, 
                                 beta_sd=beta_sd, 
                                 error_sd=error_sd,
                                 seed.problem=1,
                                 seed.noise=1) 
# X = normalizeX(data_gaus$X) # we first L2 norm normalize the columns of X
# TODO: this scaling should be done within L0glm
X = data_gaus$X
Xsp = Matrix(X, sparse=T)
y = data_gaus$y
w = rep(1, n) # observation weights

beta_true = data_gaus$beta_true
beta_true_set = which(beta_true!=0)
L1norm_betatrue = sum(abs(beta_true))
L0norm_betatrue = sum(beta_true!=0)
L0overL1norm_betatrue = L0norm_betatrue/L1norm_betatrue
L0overL1norm_betatrue # 1.29
family = gaussian(identity)


# 1.1 MCP fit using ordinis ####

library(ordinis)
gamma_mcp = 1.1 # default for MCP is 1.4, I am using value close to 1 to match best with my later L0 glm objective
timing_ordinis = system.time(ordinis_fit <- ordinis(X, y, 
               intercept=T,
               penalty = "mcp",
               gamma = gamma_mcp, 
               # lower.limits = rep(0, p), # force all coefficients to be positive
               # penalty.factor = c(0, 0, rep(1, p-2)), # don't penalize first two coefficients
               nlambda = 100L, # 100L, # for pre-screening I just need 1 solution, so nlambda=1 OK (and works best in fact)
               lambda.min.ratio = 0.01, # NULL, # =0.05 when p>n
               dfmax = n,
               alpha = 1,
               # maxit = 10,
               # tol = 0.1
               ))[[3]] # 7s with maxit 1 (with nlambda 1), 63s with defaults (with nlambda 100) 
# note: to initialise L0glm solution with MCP solution you could use fewer
# nlmabda & higher tol & smaller nr of maxit
# plot(ordinis_fit)
max(ordinis_fit$niter) # 115
which.max(logLik(ordinis_fit)) # 100
which.min(BIC(ordinis_fit)) # 100

ordinis_fit$lambda[which.max(logLik(ordinis_fit))] # best lambda

min(AIC(ordinis_fit)) #
min(BIC(ordinis_fit)) #
# coefs_ordinis = ordinis_fit$beta[,which.min(BIC(ordinis_fit))][-1]
# coefs_ordinis = ordinis_fit$beta[,which.min(AIC(ordinis_fit))][-1]
coefs_ordinis = ordinis_fit$beta[-1,which.max(logLik(ordinis_fit))] # gives same result as AIC & BIC here
coefs_ordinis[1] = 0 # first coef sometimes comes out as 1 for some reason - prob a bug in ordinis? now setting it to zero, as it's usually a false pos

sum(coefs_ordinis!=0) # 149
sum(coefs_ordinis!=0&beta_true!=0) # 24 true pos
sum(coefs_ordinis!=0&beta_true==0) # 125 false pos
L1norm_coefsordinis = sum(abs(coefs_ordinis))
L0norm_coefsordinis = sum(coefs_ordinis!=0)
L0overL1norm_coefsordinis = L0norm_betatrue/L1norm_betatrue
L0overL1norm_coefsordinis # 1.35051

beta_est0 = coefs_ordinis
range(beta_est0) # -1.990271  2.492420
range(beta_true) # -2.267722  2.365529

table(beta_est0!=0, beta_true!=0, 
      dnn=c("estimated beta nonzero", "true beta nonzero"))
#                            true beta nonzero
# estimated beta nonzero         FALSE   TRUE
#                         FALSE 999795     56
#                         TRUE     125     24

# TODO check what tol, maxit & nlambda would be sufficient to use as starting vals in our L0glm algo
# TODO check how performance would be affected by choice of gamma
# (gamma close to 1 maybe results in objective function closest to that of L0glm?)
# TODO compare with SCAD (penalty="scad")

cor(beta_est0, beta_true)  # 0.7609432

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
method="MCP_ordinis"
graph2png(file=paste0("..//benchmarks/plots/high_dimensional/",method,"_",suffix,".png"), 
          width=7, height=5)

# check metrics & quality of MCP solution
system.time(metrics_MCP <- data.frame(method=method,
                                      timing_s=timing_ordinis,
                                      do.call(cbind, get.metrics(beta = beta_est0,  # some metrics
                                                                  X = X, 
                                                                  y = y,
                                                                  family = gaussian(identity), 
                                                                  beta_true = beta_true)))) 
metrics_MCP


# 1.2 L0glm fit on MCP pre-selected variables ####

# NOTE: fitting an L0glm directly on the whole dataset results in much
# reduced sensitivity, so for ultra high dimensional situation clearly 
# best to start from MCP coefficients

# TODO check if MCP-like penalties could be used within L0glm itself
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3212875/

startcoefs = beta_est0 # MCP starting coefficients work best
# startcoefs = beta_est0_boot # bagged MCP starting coefficients (avg across 5 bootstrap replicates)
# startcoefs = coefs_glmnet1 
# NOTE: LASSO glmnet coefficient starting values also work & 
# results in smaller nr of true coefficients being picked in adaptive L0glm step (27 instead of 33), but requires more iterations than with MCP starting values) 

# priormean = startcoefs
# cutoff = quantile(abs(priormean)[priormean!=0], 0.8) 
# priormean[abs(priormean)<=cutoff] = 0 # we shrink small coefficients to zero & keep larger ones at their original value

nonzero = startcoefs!=0
# we use mbic because original problem is high dimensional
lambda_preset = preset.lambda(IC = "mbic", # options aic, bic, gic, ebic, mbic 
                              p = ncol(X), # note that we use original problem size
                              n = nrow(X))
timing_L0glmafterMCP = system.time(L0glm_gaus_presetlambda <- L0glm.fit(
  X = X[, nonzero, drop=F], 
  y = y, 
  family = family, 
  lambda = lambda_preset, 
  nonnegative = FALSE,
  # prior.mean = priormean[nonzero],
  start = startcoefs[nonzero], 
  warm.start = 0.1, 
  # TODO make sure that with warm.start=1 vars with zero coefs are efficiently removed from the start
  gamma = 2, 
  normalize = FALSE, # normalize = T doesn't work yet - too sparse then
  # dispersion = dispersion,
  solver = "eigen", # note: solver osqp better than eigen for very small lambda values))
  method = 7L, # use 7 or maybe 5?
  relaxed = FALSE))[[3]] # 0.15s
timing_L0glmafterMCP = timing_L0glmafterMCP+timing_ordinis

L0glm_gaus_presetlambda$dispersion # 
L0glm_gaus_presetlambda_coefs = rep(0, ncol(X))
L0glm_gaus_presetlambda_coefs[nonzero] = L0glm_gaus_presetlambda$coefficients

beta_est1 = L0glm_gaus_presetlambda_coefs
sum(beta_est1!=0&beta_true!=0)  # 24
table(beta_est1!=0, beta_true!=0, 
      dnn=c("estimated beta nonzero", "true beta nonzero"))
#                    true beta nonzero
# estimated beta nonzero  FALSE   TRUE
#                  FALSE 999911     56
#                   TRUE      9     24

cor(beta_est1, beta_true) # 0.8563481

plot(x=beta_est1[beta_true!=0&beta_est1!=0], 
     y=beta_true[beta_true!=0&beta_est1!=0], pch=16, col="blue",
     xlab="Estimated coefficients", ylab="True coefficients",
     main="L0glm on MCP preselected variables") # true positives
min(abs(beta_est1)[beta_true!=0&beta_est1!=0]) # 
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
method="L0glm_after_MCP"
graph2png(file=paste0("..//benchmarks/plots/high_dimensional/",method,"_",suffix,".png"), 
          width=7, height=5)

# check metrics & quality of solution with preset lambda
system.time(metrics_L0 <- data.frame(method=method,
                                     timing_s=timing_L0glmafterMCP,
                                     do.call(cbind, get.metrics(beta = beta_est1,  # some metrics
                                                                  X = X, 
                                                                  y = y,
                                                                  family = gaussian(identity), 
                                                                  beta_true = beta_true))))  # 0.5s
metrics_L0



# 1.3 EXPERIMENTAL ITERATIVE ADAPTIVE L0GLM ALGO TO FURTHER REDUCE FALSE POS RATE ####

# For the high dimensional case a solution with fewer false positives
# can be obtained by using the 1st iteration
# solution and then using the L0glm output coefficients as starting values in 
# subsequent iterations. This would be sort of like an "adaptive L0 glm" then, where
# the original starting values that are used would themselves also be optimised.
# Given that we start from sensible starting values we can specify a prior mean here.

# startcoefs = beta_est0_boot # we start from bagged MCP starting values
startcoefs = beta_est0 # we start from MCP starting values

# NOTE: starting from the MCP solution also works, but takes more iterations
# to converge
# starting from the iterative adaptive MCP solution converges very quickly
# but I need to find out how to do those fits on X subsetting to nonzero coefs

timing_adL0glmafterMCP = system.time({
 for (iter in (1:15)) {
  priormean = startcoefs # means of prior towards we will shrink
  cutoff = quantile(abs(priormean)[priormean!=0], 0.2) # find a better smooth functional form instead of a hard cutoff?
  priormean[abs(priormean)<=cutoff] = 0 # still using zero prior center for small coefficients

sel = startcoefs!=0

# we use mbic because original problem is high dimensional
lambda_preset = preset.lambda(IC = "mbic", # options aic, bic, gic, ebic, mbic 
                              p = ncol(X), # note that we use original problem size
                              n = nrow(X))

L0glm.out_updated <- L0glm.fit(X = X[,sel,drop=F], 
                                    y = y, 
                                    family = gaussian(identity), 
                                    lambda = lambda_preset,  
                                    # no.pen = 1,
                                    nonnegative = FALSE,
                                    prior.mean = priormean[sel], # we can now specify our prior centers
                                    start = startcoefs[sel], # 
                                    warm.start = 0.8, # 0.5-0.8 work well
                                    # dispersion = L0glm_gaus_presetlambda$dispersion, # L0glm.out_2nditeration$dispersion,
                                    solver = "eigen", # eigen, glmnet or osqp, # note: constrained solvers glmnet & osqp better than eigen for very small lambda values
                                    method = 7L, #
                                    gamma = 2, #
                                    thresh = 1E-5,
                                    relaxed = FALSE
)  # 

L0glm.out_updated_coefs = rep(0, ncol(X))
L0glm.out_updated_coefs[sel] = L0glm.out_updated$coefficients
print(sum(L0glm.out_updated_coefs!=0)) # 23
L0glm.out_updated$dispersion # 15.99876

beta_est1B = L0glm.out_updated_coefs

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

startcoefs = beta_est1B
}
})[[3]] 

timing_adL0glmafterMCP = timing_adL0glmafterMCP + timing_ordinis # to add cost of MCP fit

sum(beta_est1B!=0) # 23
sum(beta_est1B!=0&beta_true!=0)  # 22
L0glm.out_updated$dispersion # 
# these effect sizes in terms of Cohen's d end up as false negatives and are not detected
range(abs(beta_true)[beta_est1B==0&beta_true!=0]) 
# 0.03105975 1.03903399
# while these effect sizes in terms of Cohen's d are detected
range(abs(beta_true)[beta_est1B!=0&beta_true!=0]) 
# 0.7615661 2.3655293

table(beta_est1B!=0, beta_true!=0, 
      dnn=c("estimated beta nonzero", "true beta nonzero"))
#               true beta nonzero
# estimated beta nonzero  FALSE   TRUE
#                  FALSE 999919     58
#                   TRUE      1     22

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
# which(beta_true==0&beta_est1B!=0)
points(x=beta_est1B[beta_true==0&beta_est1B==0][[1]], 
       y=beta_true[beta_true==0&beta_est1B==0][[1]], 
       pch=16, col="green") # true negatives
method="iterative_adaptive_L0glm_after_MCP"
graph2png(file=paste0("..//benchmarks/plots/high_dimensional/",method,"_",suffix,".png"), 
          width=7, height=5)

system.time(metrics_adL0 <- data.frame(method=method,
                                       timing_s=timing_adL0glmafterMCP,
                                       do.call(cbind, get.metrics(beta = beta_est1B,  # some metrics
                                                                   X = X, 
                                                                   y = y,
                                                                   family = gaussian(identity), 
                                                                   beta_true = beta_true)))) 
metrics_adL0



# 1.4 abess fit on MCP pre-selected variables ####

library(abess)
startcoefs = beta_est0 # MCP starting coefficients work best

selcoefs = which(startcoefs!=0)
p_selcoefs = length(selcoefs)
timing_abessafterMCP = system.time(fit_abessafterMCP_selcoefs <- abess(x = X[,selcoefs,drop=F], 
                                                  y = y, 
                                                  family = "gaussian",
                                                  tune.path = "sequence", 
                                                  support.size = c(1:(p_selcoefs-1)),
                                                  lambda = 0,
                                                  ic.scale = (1/log(n))*preset.lambda(IC="mbic", p=p, n=n), # in this way supporting all IC that L0glm supports
                                                  warm.start = TRUE, 
                                                  tune.type = "bic"))[[3]]  
# 0.36s
timing_abessafterMCP = timing_abessafterMCP + timing_ordinis
plot(fit_abessafterMCP_selcoefs, type="tune")  # optimal support size 
beta_abessafterMCP_selcoefs = as.matrix(extract(fit_abessafterMCP_selcoefs)$beta) # coefficient matrix for best subset
beta_abessafterMCP = Matrix(rep(0, p), sparse=T)
beta_abessafterMCP[selcoefs] = beta_abessafterMCP_selcoefs
beta_abessafterMCP = as.vector(beta_abessafterMCP)
# beta_abessafterMCP[beta_abessafterMCP<0] = 0
sum(beta_abessafterMCP!=0) # 22
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
method="abess_after_MCP"
graph2png(file=paste0("..//benchmarks/plots/high_dimensional/",method,"_",suffix,".png"), 
          width=7, height=5)

table(beta_abessafterMCP!=0, beta_true!=0, 
      dnn=c("estimated beta nonzero", "true beta nonzero"))
#                         true beta nonzero
# estimated beta nonzero  FALSE   TRUE
#                  FALSE 999920     58
#                   TRUE      0     22


# check metrics & quality of solution 
system.time(metrics_abess_afterMCP <- data.frame(method="abess_after_MCP",
                                                 timing_s=timing_abessafterMCP,
                                                 do.call(
                                                         cbind, get.metrics(beta = beta_abessafterMCP,  # some metrics
                                                                     X = X, 
                                                                     y = y,
                                                                     family = gaussian(identity), 
                                                                     beta_true = beta_true))))  # 0.5s
metrics_abess_afterMCP


# 1.5 FULL abess FIT using MCP preselected set as initial active set ####
timing_abessMCPinitialset = system.time(fit_abess_full <- abess(x = X, 
                                                y = y, 
                                                family = "gaussian",
                                                tune.path = "sequence", 
                                                support.size = c(1:(p_selcoefs-1)),
                                                lambda = 0,
                                                ic.scale = (1/log(n))*preset.lambda(IC="mbic", p=p, n=n), # to convert BIC criterion to MBIC
                                                init.active.set = as.integer(which(startcoefs!=0)),
                                                warm.start = TRUE, 
                                                tune.type = "bic"))[[3]] # 55s
timing_abessMCPinitialset = timing_abessMCPinitialset + timing_ordinis
# NOTE: full abess fit not specifying this starting set gives nearly identical
# result but is ca 4x slower
plot(fit_abess_full, type="tune")  # optimal support size 
beta_abess = as.vector(extract(fit_abess_full)$beta) 
sum(beta_abess!=0) # 22
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
method="abess_MCP_initial_set"
graph2png(file=paste0("..//benchmarks/plots/high_dimensional/",method,"_",suffix,".png"), 
          width=7, height=5)

table(beta_abess!=0, beta_true!=0, 
      dnn=c("estimated beta nonzero", "true beta nonzero"))
# true beta nonzero
# estimated beta nonzero  FALSE   TRUE
#                  FALSE 999919     59
#                  TRUE       1     21

# check metrics & quality of solution 
system.time(metrics_abess_full <- data.frame(method="abess_MCP_initial_set",
                                             timing_s=timing_abessMCPinitialset,
                                             do.call(cbind, get.metrics(beta = beta_abess,  # some metrics
                                                                   X = X, 
                                                                   y = y,
                                                                   family = gaussian(identity), 
                                                                   beta_true = beta_true))))  # 0.5s
metrics_abess_full



# 1.6 L0Learn fit on MCP pre-selected variables ####
library(L0Learn)
startcoefs = beta_est0 # MCP starting coefficients

selcoefs = which(startcoefs!=0)
p_selcoefs = length(selcoefs)
W = rep(1, nrow(X))

timing_L0LearnafterMCP = system.time(fit_L0Learn_selcoefs <- L0Learn.fit(x=as.matrix(Diagonal(x=sqrt(W)) %*% X[,selcoefs,drop=F]),  # note: no support for sparse matrices
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
))[[3]]  # lambda of 1 is scaled by 1/(norm(weighted y,"2")) because L0Learn scales y by L2 norm of y
timing_L0LearnafterMCP = timing_L0LearnafterMCP + timing_ordinis
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
beta_L0Learn_after_MCP = as.vector(beta_L0Learn_after_MCP)
sum(beta_L0Learn_after_MCP!=0) # 69

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
method="L0Learn_after_MCP"
graph2png(file=paste0("..//benchmarks/plots/high_dimensional/",method,"_",suffix,".png"), 
          width=7, height=5)

table(beta_L0Learn_after_MCP!=0, beta_true!=0, 
      dnn=c("estimated beta nonzero", "true beta nonzero"))
#                         true beta nonzero
# estimated beta nonzero  FALSE   TRUE
#                  FALSE 999875     56
#                   TRUE     45     24


# check metrics & quality of solution with preset lambda
system.time(metrics_L0Learn <- data.frame(method="timing_L0LearnafterMCP",
                                          timing_s=timing_L0LearnafterMCP,
                                          do.call(cbind, get.metrics(beta = beta_L0Learn_after_MCP,  # some metrics
                                                                X = X, 
                                                                y = y,
                                                                family = gaussian(identity), 
                                                                beta_true = beta_true))))  # 0.5s
metrics_L0Learn
#          method timing_s       p   n m deviance deviance_true   loglik     MD_coefs
# 1 timing_L0LearnafterMCP    64.59 1000001 500 1 1197.029      922.3419 999003.3 6.641015e-06
#     MAD_coefs  sMAPE_coefs bias_coefs  R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs
# 1 5.36834e-05 0.0002062479  0.8013139 0.5918034  0.03001058    0.4203605       1.365911
#   bias_selcoefs R2_selcoefs support_size_true support_size_estimated
# 1     0.8013139   0.7542867                80                     69
#   estimated_over_true_support_size TP FP     TN FN precision recall specificity accuracy
# 1                           0.8625 24 45 999875 56 0.3478261    0.3    0.999955 0.999899
#    F1_score      gic     ebic      bic      aic     mbic
# 1 0.3221477 -1996240 -1996072 -1997572 -1997867 -1995831



# 1.7 LASSO fit using glmnet ####

library(glmnet)
timing_cvglmnet = system.time(glmnet_cvfit <- cv.glmnet(x = X, y = y, nfolds = 5,
                                      intercept=F, alpha=1))[[3]] # 136s
plot(glmnet_cvfit)
bestlambda1 <- glmnet_cvfit$lambda.min #
bestlambda1 # 0.6624504
bestlambda2 <- glmnet_cvfit$lambda.1se #
bestlambda2 # 1.394394
timing_glmnetfull = system.time(glmnet_fit <- glmnet(x = X, y = y, intercept = F, 
                                                     alpha = 1))[[3]] # 25s
timing_glmnet = timing_cvglmnet+timing_glmnetfull
plot(glmnet_fit)
coefs_glmnet1 <- as.vector(coef(glmnet_fit, s=bestlambda1)[-1])
sum(coefs_glmnet1!=0) # 205
glmnet_set1 <- which(coefs_glmnet1!=0)
length(intersect(glmnet_set1, beta_true_set)) # 23
table(coefs_glmnet1!=0, beta_true!=0,  
      dnn=c("estimated beta nonzero", "true beta nonzero"))
#                    true beta nonzero
# estimated beta nonzero  FALSE   TRUE
# FALSE                  999738     57
# TRUE                      182     23

beta_LASSO = coefs_glmnet1
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
method="glmnet_LASSO"
graph2png(file=paste0("..//benchmarks/plots/high_dimensional/",method,"_",suffix,".png"), 
          width=7, height=5)

# check metrics & quality of LASSO solution
system.time(metrics_LASSO <- data.frame(method="glmnet_LASSO",
                                        timing_s=timing_glmnet,
                                        do.call(cbind, get.metrics(beta = coefs_glmnet1,  # some metrics
                                                                   X = X, 
                                                                   y = y,
                                                                   family = gaussian(identity), 
                                                                   beta_true = beta_true)))) 
metrics_LASSO

# note: the cross validation I might not really need, as I can just select a model
# with max n variables
# 
# # note: adaptive LASSO fit also doesn't help to select more variables


# 1.8 relaxed LASSO fit using glmnet ####
# see https://cran.r-project.org/web/packages/glmnet/vignettes/relax.pdf

library(glmnet)
timing_cvglmnetrelaxed = system.time(glmnet_cvfitrelaxed <- cv.glmnet(x = X, y = y, nfolds = 5,
                                                        intercept=F, alpha = 1,
                                                        gamma = c(0, 0.25, 0.5, 0.75, 1), # fully relaxed fit for 5 gamma values (0=fully relaxed)
                                                        relax = T, # relaxed LASSO
                                                        maxp=round(n/2)))[[3]] # 629s
plot(glmnet_cvfitrelaxed)
bestlambda1 <- glmnet_cvfitrelaxed$relaxed$lambda.min #
bestlambda1 # 
bestlambda2 <- glmnet_cvfitrelaxed$relaxed$lambda.1se #
bestlambda2 # 
bestgamma <- glmnet_cvfitrelaxed$relaxed$gamma.min
bestgamma # 0.25 - close to fully relaxed (0)
timing_glmnetfullrelaxed = system.time(glmnet_fitrelaxed <- glmnet(x = X, y = y, 
                                                            intercept = F, 
                                                            alpha = 1,
                                                            gamma = c(0, 0.25, 0.5, 0.75, 1),
                                                            relax = T,
                                                            maxp=round(n/2)))[[3]] # 122s
timing_glmnetrelaxed = timing_cvglmnetrelaxed+timing_glmnetfullrelaxed
plot(glmnet_fitrelaxed)
abline(v=log(bestlambda1), col="red")
coefs_glmnet1relaxed <- as.vector(coef(glmnet_fitrelaxed$relaxed, 
                                       s=bestlambda1, gamma=bestgamma)[-1])
sum(coefs_glmnet1relaxed!=0) # 24
glmnet_set1relaxed <- which(coefs_glmnet1relaxed!=0)
length(intersect(glmnet_set1relaxed, beta_true_set)) # 12
table(coefs_glmnet1relaxed!=0, beta_true!=0,  
      dnn=c("estimated beta nonzero", "true beta nonzero"))
#                    true beta nonzero
# estimated beta nonzero  FALSE   TRUE
# FALSE                  999908     68
# TRUE                       12     12

beta_LASSOrelaxed = coefs_glmnet1relaxed
plot(x=beta_LASSOrelaxed[beta_true!=0&beta_LASSOrelaxed!=0], 
     y=beta_true[beta_true!=0&beta_LASSOrelaxed!=0], pch=16, col="blue",
     xlab="Estimated coefficients", ylab="True coefficients",
     main=paste0("relaxed LASSO fit (glmnet)")) # true positives
abline(a=0, b=1, col="lightgrey")
points(x=beta_LASSOrelaxed[beta_true!=0&beta_LASSOrelaxed==0], 
       y=beta_true[beta_true!=0&beta_LASSOrelaxed==0], 
       pch=16, col="purple") # false negatives
points(x=beta_LASSOrelaxed[beta_true==0&beta_LASSOrelaxed!=0], 
       y=beta_true[beta_true==0&beta_LASSOrelaxed!=0], 
       pch=16, col="red") # false positives
points(x=beta_LASSOrelaxed[beta_true==0&beta_LASSOrelaxed==0][[1]], 
       y=beta_true[beta_true==0&beta_LASSOrelaxed==0][[1]], 
       pch=16, col="green") # true negatives
method="glmnet_relaxed_LASSO"
graph2png(file=paste0("..//benchmarks/plots/high_dimensional/",method,"_",suffix,".png"), 
          width=7, height=5)

# check metrics & quality of LASSO solution
system.time(metrics_relaxedLASSO <- data.frame(method=method,
                                        timing_s=161,
                                        do.call(cbind, get.metrics(beta = coefs_glmnet1relaxed,  # some metrics
                                                                   X = X, 
                                                                   y = y,
                                                                   family = gaussian(identity), 
                                                                   beta_true = beta_true)))) 
metrics_relaxedLASSO



# 1.9 COMPARISON OF METRICS ####

library(dplyr)
benchmarks = bind_rows(metrics_LASSO, 
                       metrics_relaxedLASSO,
                       metrics_MCP,
                       metrics_L0,
                       metrics_adL0,
                       metrics_abess_afterMCP,
                       metrics_abess_full,
                       metrics_L0Learn)
benchmarks = benchmarks[order(benchmarks$F1_score, decreasing=T),]
benchmarks
#                               method timing_s       p   n m   deviance deviance_true   loglik
# 6                    abess_after_MCP    64.90 1000001 500 1  7684.0557     7135.9699 998585.5
# 5 iterative_adaptive_L0glm_after_MCP    68.15 1000001 500 1  7631.4103     7167.0018 998586.2
# 4                    L0glm_after_MCP    64.82 1000001 500 1  6403.8240     5923.6036 998620.0
# 7              abess_MCP_initial_set   169.41 1000001 500 1  7440.8172     7942.4573 998593.5
# 8             timing_L0LearnafterMCP    64.59 1000001 500 1  1197.0288      922.3419 999003.3
# 2               glmnet_relaxed_LASSO   161.00 1000001 500 1 10933.4488    10534.1682 998495.3
# 3                        MCP_ordinis    64.54 1000001 500 1   714.5266      214.2935 999052.3
# 1                       glmnet_LASSO   164.55 1000001 500 1  7853.6100     7695.3407 998397.0
#       MD_coefs    MAD_coefs  sMAPE_coefs bias_coefs  R2_coefs  MD_selcoefs MAD_selcoefs
# 6 5.350454e-06 3.053394e-05 0.0001186759  0.9727985 0.7418590 0.0427150500    0.1772086
# 5 7.515924e-06 3.180458e-05 0.0001218726  0.9867476 0.7276358 0.0483509400    0.2301083
# 4 8.176878e-06 3.322297e-05 0.0001364081  1.0715341 0.7333310 0.1092906943    0.2589226
# 7 5.555376e-06 3.201309e-05 0.0001222228  0.9500401 0.7162749 0.0866463333    0.2098259
# 8 6.641015e-06 5.368340e-05 0.0002062479  0.8013139 0.5918034 0.0300105804    0.4203605
# 2 6.989105e-06 5.201936e-05 0.0001626092  0.9034791 0.3866253 0.0159113220    0.5419753
# 3 8.375465e-06 6.366255e-05 0.0003657603  0.7699911 0.5790339 0.0255381181    0.2616378
# 1 6.152113e-06 6.348239e-05 0.0005061737  1.8737438 0.4359878 0.0009656156    0.1825365
#   sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true support_size_estimated
# 6      0.1216335     0.9727985   0.9808786                80                     22
# 5      0.2553290     0.9867476   0.9575115                80                     23
# 4      0.7396381     1.0715341   0.9346716                80                     33
# 7      0.1919461     0.9500401   0.9578239                80                     22
# 8      1.3659110     0.8013139   0.7542867                80                     69
# 2      1.1087186     0.9034791   0.7446841                80                     24
# 3      1.7030892     0.7699911   0.7380113                80                    149
# 1      1.9130426     1.8737438   0.5757616                80                    205
#   estimated_over_true_support_size TP  FP     TN FN precision recall specificity accuracy
# 6                           0.2750 22   0 999920 58 1.0000000 0.2750    1.000000 0.999942
# 5                           0.2875 22   1 999919 58 0.9565217 0.2750    0.999999 0.999941
# 4                           0.4125 24   9 999911 56 0.7272727 0.3000    0.999991 0.999935
# 7                           0.2750 21   1 999919 59 0.9545455 0.2625    0.999999 0.999940
# 8                           0.8625 24  45 999875 56 0.3478261 0.3000    0.999955 0.999899
# 2                           0.3000 12  12 999908 68 0.5000000 0.1500    0.999988 0.999920
# 3                           1.8625 24 125 999795 56 0.1610738 0.3000    0.999875 0.999819
# 1                           2.5625 23 182 999738 57 0.1121951 0.2875    0.999818 0.999761
#    F1_score      gic     ebic      bic      aic     mbic
# 6 0.4313725 -1996590 -1996535 -1997028 -1997125 -1996456
# 5 0.4271845 -1996567 -1996509 -1997023 -1997124 -1996427
# 4 0.4247788 -1996382 -1996301 -1997029 -1997172 -1996184
# 7 0.4117647 -1996606 -1996551 -1997044 -1997141 -1996472
# 8 0.3221477 -1996240 -1996072 -1997572 -1997867 -1995831
# 2 0.2307692 -1996360 -1996300 -1996835 -1996941 -1996214
# 3 0.2096070 -1994319 -1993960 -1997172 -1997805 -1993444
# 1 0.1614035 -1991595 -1991102 -1995514 -1996382 -1990393

# CONCLUSION: in terms of F score to recover the true support,
# abess after MCP & adaptive L0glm are both very good
# LASSO & MCP far too many false positives
# relaxed LASSO has poor 

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
system.time(metrics_MCP_boot <- data.frame(method="MCP_ordinis",
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


