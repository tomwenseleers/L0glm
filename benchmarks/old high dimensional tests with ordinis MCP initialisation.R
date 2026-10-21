## HIGH DIMENSIONAL TESTS ####
# cf https://github.com/jaredhuling/ordinis

library(matrixStats)
library(qlcMatrix)

set.seed(123)
n <- 500
p <- 1E6
sd <- 1
x <- matrix(rnorm(n*p, sd = sd), n, p)  # orthonormal design matrix
k <- 50 # 50 nonzero coefs
beta <- rep(0, p)
beta[sample(1:p, k, replace=F)] = runif(k, 0, 1)
disp <- 1 # dispersion = variance
noise <- rnorm(n, sd = sqrt(disp))
y <- x%*%beta+noise
disp_est_crude <- var(y-mean(y)) # crude estimate of error variance if all model coefs were zero
disp_est_crude # 19.48596
# this is a factor 
sqrt(disp_est_crude)/sqrt(disp) # 4.41429
# too high

# L0 penalty used in L0glm & expected detection threshold
lambdaGIC <- log(p)*log(log(n)) # for GIC/SIC as in abess
# expected detection threshold with lambdaGIC L0 penalty
lambdaGICdetthresh <- sqrt(disp*lambdaGIC/n)
lambdaGICdetthresh # 0.2246757
lambdaGICdetthresh_crude <- sqrt(disp_est_crude*lambdaGIC/n)
lambdaGICdetthresh_crude # 0.9917838

# universal penalty level used in MCP pen regression by Zhang 2010
# that would enjoy variable selection consistency (Donoho & Johnston 1994)
lambda_universal = sqrt(disp*log(p)*2/n)
lambda_universal # 0.2350788
# NOTE: this lambda_universal is less strict than 
# lambda_GIC=sqrt(disp*lambdaGIC/n)=sqrt(disp*log(p)*log(log(n))/n) for
# n < exp(exp(2)) = 1618 but 
# slightly more strict for larger n (e.g. for n=1E6 log(log(n)) it is 2.6 instead of 2)
# given that we don't know the dispersion a priori the best we could do initially is use
lambda_universal_crude = sqrt(disp_est_crude*log(p)*2/n)
lambda_universal_crude # 1.037706


# library(selectiveInference)
# system.time(disp_est <- estimateSigma(x, y, intercept=F, standardize=T)$sigmahat) # 1.49s
# disp_est # 4.981075

# MCP fit over lambda grid
# to use as starting values in L0glm fit
library(ordinis)
system.time(mcpmod <- ordinis(x, y, 
                              intercept = F, 
                              penalty = "mcp",
                              gamma = 1.1, # default param in MCP objective = 1.4, gamma->1 converges to hard thresholding, i.e. an L0 penalty
                              # lambda = c(bestlambdapreset*4, bestlambdapreset*2, bestlambdapreset),
                              # lambda = c(bestlambdapreset_crude*4, bestlambdapreset_crude*2, bestlambdapreset_crude, bestlambdapreset_crude/2, bestlambdapreset_crude/4, bestlambdapreset_crude/6, bestlambdapreset_crude/8, bestlambdapreset_crude/10, bestlambdapreset_crude/20, bestlambdapreset_crude/100, bestlambdapreset_crude/1000),
                              nlambda = 10,
                              lower.limits = rep(0, p), # force all coefficients to be positive
                              maxit = 5,
                              maxit.irls = 1,
                              dfmax = n-1, # min(n-1, round(2*n/(log(log(n))*log(p)))),
                              alpha = 1,
                              tol = 1E-3))  # 21s
plot(mcpmod)
mcpmod$lambda
lambda_universal # 0.2350788 - this would be lambda 5, but we do not know disp in advance...
max(which(mcpmod$lambda>lambda_universal)) # 5
max(which(mcpmod$lambda>lambdaGICdetthresh)) # 6
colSums(mcpmod$beta!=0)
# the max coefficients observed for any of the lambda values works best
# to use as starting coefficients in L0glm
mcpmaxcoefs = as.numeric(rowMax(mcpmod$beta)[-1])

# loglik = as.numeric(logLik(mcpmod))
# edf = colSums(mcpmod$beta!=0)
# GIC = 2*edf - log(p)*log(log(n)) * loglik
# which.min(GIC) # 10 - not quite right (overfitted) - would a data splitting approach work?
# EBIC = 2*edf - (log(n) + 2*(1-log(n)/(2*log(p)))*log(p)) * loglik
# which.min(EBIC)
# MBIC = 2*edf - log(n*(p^2)/16) * loglik
# which.min(MBIC)
# bestlambdaid = which.min(abs(mcpmod$lambda-bestlambda))

bestlambdaid = max(which(colSums(mcpmod$beta!=0)<n)) 
# we could also take estimates with smallest lambda value that has <n coefficients retained 
# ALT: use solution with best IC, e.g. GIC 
ncol(mcpmod$beta) 
mcpcoefs <- as.numeric(mcpmod$beta[,bestlambdaid])[-1]
cor(mcpcoefs, beta) # 0.961909
summary(lm(beta~0+mcpcoefs))$r.squared # 0.925271
table(mcpcoefs!=0, beta!=0)
#         FALSE  TRUE
# FALSE  999949     11
# TRUE       56     39
table(mcpmaxcoefs!=0, beta!=0)
#        FALSE   TRUE
# FALSE 999892     11
# TRUE      58     39
plot(x=mcpcoefs[mcpcoefs!=0|beta!=0], y=beta[mcpcoefs!=0|beta!=0], 
     pch=16, col="steelblue")
min(abs(mcpcoefs[beta!=0&mcpcoefs>0])) # 0.2585947


library(glmnet)
# system.time(cfit <- cv.glmnet(x, y, relax = TRUE, family="gaussian", # standardize = F,
#                               lower.limits = rep(0, p))) # 2812s=46 mins - VERY SLOW
# plot(cfit)
# bestlambda <- cfit$relaxed$lambda.min
# bestlambda  # 1.841157
# bestgamma <- cfit$relaxed$gamma.min
# bestgamma # 0 - so fully relaxed
system.time(relaxlassofit <- glmnet(x, y, 
                                    relax = TRUE, 
                                    gamma = 0,
                                    intercept = FALSE,
                                    family="gaussian", 
                                    standardize = FALSE,
                                    standardize.response = FALSE,    
                                    dfmax = n-1,
                                    lower.limits = rep(0, p), 
                                    nlambda=10)) # 28s
plot(relaxlassofit, gamma=0)
relaxlassofit$relaxed$lambda
supsize <- colSums(relaxlassofit$relaxed$beta!=0)
dfmax = n-1 # min(n-1, 6*n/(log(log(n))*log(p)))
best = 7 # max(which(supsize<=dfmax))
best # 7
relaxlassocoefs <- as.numeric(relaxlassofit$relaxed$beta[,best])
cor(relaxlassocoefs, beta) # 0.7475808
# the max coefficients observed for any of the lambda values to use as
# starting coefficients in L0glm
relaxlassomaxcoefs = as.numeric(rowMax(relaxlassofit$relaxed$beta))
cor(relaxlassomaxcoefs, beta) # 0.7710445
table(relaxlassocoefs!=0, beta!=0)
#        FALSE   TRUE
# FALSE 999538     24
# TRUE     412     26
table(relaxlassomaxcoefs!=0, beta!=0)
#        FALSE   TRUE
# FALSE 999503     24
# TRUE     447     26
plot(x=relaxlassocoefs[relaxlassocoefs!=0|beta!=0], 
     y=beta[relaxlassocoefs!=0|beta!=0], pch=16, col="steelblue")
min(abs(relaxlassocoefs[beta!=0&relaxlassocoefs!=0])) # 0.05134241

library(abess)
system.time(res <- abess(x, as.vector(y),
                         support.size=1:min(n, round(3*n/(log(log(n))*log(p)))), # 
                         # init.active.set=which(mcpcoefs!=0),
                         fit.intercept=F,
                         beta.low=0,
                         family="gaussian", tune.type="gic")) 
# 56s
abesscoefs = as.numeric(coef(res, support.size=res$best.size))[-1]
cor(abesscoefs, beta) # 0.9856229
plot(x=abesscoefs[abesscoefs!=0|beta!=0], 
     y=beta[abesscoefs!=0|beta!=0], pch=16, col="steelblue")
table(abesscoefs!=0, beta!=0)
#        FALSE   TRUE
# FALSE 999948     11
# TRUE       2     39

library(L0glm)
library(matrixStats)
# start = relaxlassomaxcoefs  # not a great solution
start = mcpmaxcoefs  # WORKS BEST
# start[mcpcoefs!=0] = 1 # to give selected variables equal weight - this also works 
# but can result in small nr of false pos which would need selecting out via iterating 
# the algo via an "adaptive L0glm fit"
nonzero = which(mcpmaxcoefs!=0)
mcpfitted = as.numeric(x %*% mcpmaxcoefs)
disp_est = sum((mcpfitted - y)^2)/(n - sum(mcpmaxcoefs!=0))
disp_est # 0.5469727 -  too small as it's overfitted
system.time(res_nonzero <- L0glm.fit(x[,nonzero,drop=F], y, 
                                     family=gaussian(identity),
                                     start=start[nonzero],
                                     lambda=lambdaGIC,
                                     nonnegative=TRUE,
                                     # always.include=1:p,
                                     # normalize=TRUE,
                                     # dispersion=disp_est,
                                     solver="eigen", 
                                     method=4L # LSCG, fastest for large dense
)) # 0.06s - TODO kick out zero coefs from model from the start in L0glm
L0glmcoefs = rep(0, p)
L0glmcoefs[nonzero] = as.numeric(res_nonzero$coefficients)
cor(L0glmcoefs, beta) # 0.9881105
table(L0glmcoefs!=0, beta!=0)
#        FALSE   TRUE
# FALSE 999948     11
# TRUE       2     39 
# PS sometimes some false positives, caused by dispersion being estimated 
# too low
# TODO check if this can be improved
# could be solved via iterative adaptive L0glm or setting gamma slightly above 2
# NOTE: we could also threshold coefs < sqrt(res_nonzero$dispersion*lambdaGIC/n)
# (the expected detection threshold when using lambdaGIC for orthogonal design)
# to zero??
res_nonzero$dispersion # 0.9036417 # 23.24771
plot(x=L0glmcoefs[L0glmcoefs!=0|beta!=0], 
     y=beta[L0glmcoefs!=0|beta!=0], pch=16, col="steelblue")
min(abs(L0glmcoefs[beta!=0&L0glmcoefs>0])) # 0.2667079
sqrt(disp*lambdaGIC/n) # 0.2246757 - approx matches MCP detection threshold lambda_universal
lambda_universal # 0.2350788
summary(lm(beta~0+L0glmcoefs))$r.squared # 0.976363
# with start = mcpmaxcoefs
sqrt(res_nonzero$dispersion*lambdaGIC/n) # 0.2135769 - lambdaGIC expected detection threshold
sqrt(res_nonzero$dispersion*log(p)*2/n) # 0.2234661 - universa
