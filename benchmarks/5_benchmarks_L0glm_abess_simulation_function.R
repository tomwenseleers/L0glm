# SOME L0GLM BENCHMARKS BASED ON THE ABESS BENCHMARKS ####
# see https://github.com/abess-team/abess/blob/master/docs/simulation/R/run_benchmark.R
# and https://github.com/abess-team/abess

library(abess)
devtools::load_all(".")

# 1. MULTICHANNEL MULTIVARIATE GAUSSIAN EXAMPLE ####

n = 500 # obs # or 1000 for logistic in https://github.com/abess-team/abess/blob/master/docs/simulation/R/run_benchmark.R
p = 1000 # features
y.dim = 10 # nr of tasks/channels
k = 100 # features with nonzero coefs
rho = 0.7 # 0.1 or 0.7 in https://github.com/abess-team/abess/blob/master/docs/simulation/R/run_benchmark.R
cortype = ifelse(rho <= 0.5, 2, 3) # 1 (independent), 2 (exponential structure) or 3 (constant structure)
# sigma = 0.1 # sigma does not work with generate.data for multichannel data
snr = 0.1 # signal to noise ratio
family = "mgaussian"
crit = "mbic" # using mbic as p > n
seed = 1

data.mc.gaus = generate.data(n, 
                             p, 
                             support.size = k, 
                     family = "mgaussian", # here simulating multitask case with Gaussian noise, but also interested in multitask case with Poisson noise (& identity link)
                     y.dim = y.dim,
                     cortype = cortype, # 1 (independent), 2 (exponential structure) or 3 (constant structure)
                     rho = rho, # pairwise correlation in predictors
                     snr = snr, # signal to noise ratio
                     # sigma = sigma, # does not appear to work with multichannel data
                     seed = seed)
coef_real_mc = data.mc.gaus$beta
range(coef_real_mc) # 0.00000 83.10339

image(x=1:p, y=1:y.dim, z=as.matrix(data.mc.gaus$beta)^1, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Feature", ylab="Task", main="true coefficients")
abline(v=(1:p)[as.vector(rowMin(data.mc.gaus$beta)!=0)], col="red")


## 1.1 Multichannel abess fit ####

system.time(abess_fit_mc <- abess(
  x=data.mc.gaus$x, 
  y=data.mc.gaus$y,
  # weight=1/(data.mc.gaus$y+0.1), # I'd like to use this to approximate identity link Poisson case - now gives error "Error in weight.glm(para) : is.vector(para$weight) is not TRUE"
  family="mgaussian",
  beta.low=0, # nonnegativity constraints
  support.size=1:round(n*0.2), # 1:(n-1), # 
  fit.intercept=FALSE,
  tune.path="sequence",
  tune.type = "bic", # trick to support any IC that L0glm supports
  ic.scale = preset.lambda(IC = crit, # options aic, bic, gic, mebic, ebic, mbic
                           p = p, 
                           n = n)/log(n) # trick to support any IC that L0glm supports
)) # 20s with support.size=1:(n-1), 0.31s with support.size=1:round(n*0.2) 

coef_abess_mc = coef(abess_fit_mc, 
                  support.size = abess_fit_mc[["best.size"]], 
                  sparse = FALSE)[[1]][-1,]
# NOTE: extract(abess_fit_mc)
# gives error with fit.intercept=FALSE
# Error in h(simpleError(msg, call)) : 
# error in evaluating the argument 'x' in selecting a method for function 'as.matrix': incorrect number of dimensions
# this looks like a bug

R2 = round(cor(as.vector(rowMeans(coef_real_mc)), as.vector(rowMeans(coef_abess_mc))), digits = 4)^2 # Calculate R2
title = paste0("True vs estimated coefficients\nnonnegative mgaussian abess fit\nR2=", R2)
plot(as.vector(rowMeans(coef_abess_mc)), 
     as.vector(rowMeans(coef_real_mc)), pch=16, 
     xlab="estimated coefficients",
     ylab="true coefficients", col="steelblue",
     main=title)

image(x=1:p, y=1:y.dim, z=as.matrix(coef_abess_mc)^1, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Feature", ylab="Task", main="nonnegative mgaussian abess fit\n(red = true nonzero features, cyan = estimated nonzero features)")
abline(v=(1:p)[as.vector(rowMin(coef_real_mc)!=0)], col="red")
abline(v=(1:p)[as.vector(rowMin(coef_abess_mc)!=0)], col="cyan")

cor(as.vector(as.vector(coef_abess_mc)), as.vector(coef_real_mc)) 
# 0.8500504

table(rowMeans(coef_real_mc)!=0,
      rowMeans(coef_abess_mc)!=0,
      dnn=c("estimated mean beta nonzero", "true mean beta nonzero"))
#                             true mean beta nonzero
# estimated mean beta nonzero FALSE TRUE
# FALSE                         900    0
# TRUE                           27   73




## 1.2 Multichannel Gaussian L0glm fit ####

# system.time(tune.chol.mgaus <- L0glm.tune(X=data.mc.gaus$x, y=data.mc.gaus$y, 
#                                        IC = crit, 
#                                        plot = FALSE, coefs.true = coef_real_mc, 
#                                        grouped = TRUE,
#                                        relaxed = TRUE, 
#                                        family = gaussian(identity), 
#                                        nonnegative = TRUE, 
#                                        parallelize = TRUE, 
#                                        solver = "eigen", method = 1L))
# # XXs
# # TODO fix bug
# # Error in lm.fit(x, y, offset = offset, singular.ok = singular.ok, ...) : 
# # 0 (non-NA) cases
# # In addition: There were 50 or more warnings (use warnings() to see the first 50)

# just using a preset lambda, preset to optimise the given IC
system.time(L0glm.chol.mgaus <- L0glm.fit(data.mc.gaus$x, data.mc.gaus$y, 
                                         family = gaussian(identity), 
                                         lambda = crit,
                                         start = NULL, # TODO try different possible initialisations
                                         nonnegative = TRUE,
                                         grouped = TRUE,
                                         relaxed = TRUE, # gives best result for multichannel
                                         parallelize = TRUE,
                                         solver = "eigen", # osqp gives very similar solution
                                         method = 1L))  # use eigen Cholesky LLT solver
# 12s

coef_L0glm.chol.mgaus = L0glm.chol.mgaus$coefficients

R2 = round(cor(as.vector(rowMeans(coef_real_mc)), as.vector(rowMeans(coef_L0glm.chol.mgaus))), digits = 4)^2 # Calculate R2
title = paste0("True vs estimated coefficients\nrelaxed nonnegative L0 grouped GLM\nR2=", R2)
plot(as.vector(rowMeans(coef_L0glm.chol.mgaus)), 
     as.vector(rowMeans(coef_real_mc)), pch=16, 
     xlab="estimated coefficients",
     ylab="true coefficients", col="steelblue",
     main=title)

image(x=1:p, y=1:y.dim, z=as.matrix(coef_L0glm.chol.mgaus)^1, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Feature", ylab="Task", main="relaxed nonnegative group L0glm\n(red = true nonzero features, cyan = estimated nonzero features)")
abline(v=(1:p)[as.vector(rowMin(coef_real_mc)!=0)], col="red")
abline(v=(1:p)[as.vector(rowMin(coef_L0glm.chol.mgaus)!=0)], col="cyan")

cor(as.vector(as.vector(coef_L0glm.chol.mgaus)), as.vector(coef_real_mc)) 
# 0.9728172

table(rowMeans(coef_real_mc)!=0,
      rowMeans(L0glm.chol.mgaus$coefficients)!=0,
      dnn=c("estimated mean beta nonzero", "true mean beta nonzero"))
#                             true mean beta nonzero
# estimated mean beta nonzero FALSE TRUE
# FALSE                         883   17
# TRUE                            0  100


# TODO also add glmnet mgaussian fit


## 1.3 BENCHMARK TABLE

# XXX

# performance for multichannel case clearly much better than abess
# faster & higher TPR & lower FPR & better predictive performance



# 2. SINGLE CHANNEL GAUSSIAN EXAMPLE ####

n = 500 # obs
p = 1000 # features
k = 10 # features with nonzero coefs
rho = 0.7 # 0.1 or 0.7 in https://github.com/abess-team/abess/blob/master/docs/simulation/R/run_benchmark.R
cortype = ifelse(rho <= 0.5, 2, 3) # 1 (independent), 2 (exponential structure) or 3 (constant structure)
# sigma = 1
snr = 100 # signal to noise ratio (10 is default)
family = "gaussian"
crit = "mbic" # using mbic as p > n
seed = 1

data.sc.gaus = generate.data(n, 
                             p, 
                             support.size = k, 
                             family = "gaussian", # here simulating case with Gaussian noise, but also interested in case with Poisson noise (& identity link)
                             cortype = cortype, # 1 (independent), 2 (exponential structure) or 3 (constant structure)
                             rho = rho, # pairwise correlation in predictors
                             # sigma = sigma,
                             snr = snr, # signal to noise ratio
                             seed = seed)
coef_real_sc = data.sc.gaus$beta
range(coef_real_sc) # 0.00000 75.63361


## 2.1 Single channel abess fit ####

system.time(abess_fit_sc <- abess(
  x=data.sc.gaus$x, 
  y=data.sc.gaus$y,
  # weight=as.vector(1/(pmax(data.sc.gaus$y,0)+0.1)), 
  # family="gaussian",
  # beta.low=0, # nonnegativity constraints
  # fit.intercept=FALSE,
  support.size=1:round(n*0.2), # NOTE: 1:(n-1) would be general case & is much slower: 21s
  tune.path = "sequence",
  # tune.type = "cv", num.threads = 5 # slower this tuning based on IC
  tune.type = "bic", # trick to support any IC that L0glm supports
  ic.scale = preset.lambda(IC = crit, # options aic, bic, gic, mebic, ebic, mbic
                           p = p, 
                           n = n)/log(n) # trick to support any IC that L0glm supports
)) # 21s with support.size=1:(n-1), 0.2s with support.size=1:round(n*0.2)

coef_abess_sc = coef(abess_fit_sc, 
                  support.size = abess_fit_sc[["best.size"]], 
                  sparse = FALSE)[-1,]
# NOTE: extract(coef_abess_sc)$beta also works

R2 = round(cor(as.vector(coef_real_sc), as.vector(coef_abess_sc)), digits = 4)^2 # Calculate R2
title = paste0("True vs estimated coefficients\nabess\nR2=", R2)
plot(as.vector(coef_abess_sc), 
     as.vector(coef_real_sc), pch=16, 
     xlab="estimated coefficients",
     ylab="true coefficients", col="steelblue",
     main=title)

table(coef_abess_sc!=0, 
      coef_real_sc!=0,
      dnn=c("estimated coefs nonzero", "true coefs nonzero"))
#                         true coefs nonzero
# estimated coefs nonzero FALSE TRUE
#                   FALSE   988    1
#                    TRUE     2    9



## 2.2 Single channel L0glm fit ####

# system.time(tune.chol.gaus <- L0glm.tune(X=data.sc.gaus$x, y=data.sc.gaus$y, 
#                                        IC = crit, 
#                                        plot = FALSE, coefs.true = coef_real_sc, 
#                                        grouped = FALSE,
#                                        relaxed = FALSE, 
#                                        family = gaussian(identity), 
#                                        nonnegative = TRUE, 
#                                        parallelize = TRUE, 
#                                        solver = "eigen", method = 1L))
# # XXs
# # TODO fix bug
# # Error in lm.fit(x, y, offset = offset, singular.ok = singular.ok, ...) : 
# # 0 (non-NA) cases
# # In addition: There were 50 or more warnings (use warnings() to see the first 50)

# just using a preset lambda, preset to optimise the given IC
system.time(L0glm.chol.gaus <- L0glm.fit(data.sc.gaus$x, data.sc.gaus$y, 
                                          family = gaussian(identity), 
                                          lambda = crit,
                                          start = NULL, # rep(1, p), # NULL, # TODO try different possible initialisations
                                          nonnegative = TRUE,
                                          grouped = FALSE,
                                          relaxed = FALSE, # gives best result for multichannel
                                          parallelize = FALSE,
                                          solver = "eigen", # osqp gives very similar solution
                                          method = 1L))  # use Cholesky LLT solver
# 0.26s

coef_L0glm.chol.gaus = L0glm.chol.gaus$coefficients

R2 = round(cor(as.vector(coef_real_sc), as.vector(coef_L0glm.chol.gaus)), digits = 4)^2 # Calculate R2
title = paste0("True vs estimated coefficients\nnonnegative L0 GLM\nR2=", R2)
plot(as.vector(coef_L0glm.chol.gaus), 
     as.vector(coef_real_sc), pch=16, 
     xlab="estimated coefficients",
     ylab="true coefficients", col="steelblue",
     main=title)

table(coef_L0glm.chol.gaus!=0, 
      coef_real_sc!=0,
      dnn=c("estimated coefs nonzero", "true coefs nonzero"))
# true coefs nonzero
# estimated coefs nonzero FALSE TRUE
#                   FALSE   989    1
#                    TRUE     1    9


## 2.3 BENCHMARK TABLE

# XXX

# performance for single channel similar to abess in terms of solution quality
# but faster for support size 1:(n-1)
