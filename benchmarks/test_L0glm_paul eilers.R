# Minimal tests of L0glm package ####

library(L0glm)
library(export)

n <- p <- 100 # to check performance use e.g. n=p=10000, otherwise use e.g. n=p=100-300
s <- 0.1 # sparsity, nr of timepoints with a peak 
k <- round(n*s)
mean_beta <- 2000
sd_logbeta <- 1
sd_noise <- 5 # SD of noise for Gaussian noise case
crit <- "bic" # choose lambda based on what would be expected to optimise this IC
# options in order of increasing regularisation: aic, bic, gic, ebic, mbic

# 1. SINGLE CHANNEL DECONVOLUTION ####
## 1.1 Single channel deconvolution assuming Poisson noise ####
set.seed(1)
sim <- simulate_spike_train(n = n, 
                            p = p,
                            k = k,
                            family = "poisson", 
                            mean_beta = mean_beta,
                            sd_logbeta = sd_logbeta)
X <- sim$X
y <- sim$y

# using lambda 1
L0glm.out <- L0glm.fit_new(X, y, 
                           family = poisson(identity), 
                           lambda = 1, 
                           nonzero.centered = FALSE, # use nonzero centered ridge regression?
                           nonnegative = TRUE)
# points(sim$x, L0glm.out, pch = 16, col= 'blue')

# using lambda expected to optimise given IC
lambda_preset = preset.lambda(IC = crit, # options aic, bic, gic, ebic, mbic 
                                  y = y, 
                                  X = X, 
                                  family = poisson(identity)) 
lambda_preset # lambda = 1.15 would be expected to optimise BIC with n=p=100, lambda = 1.76 would optimise GIC

L0glm.out <- L0glm.fit_new(X, y, 
                               family = poisson(identity), 
                               lambda = lambda_preset, 
                               nonzero.centered = FALSE, # use nonzero centered ridge regression?
                               nonnegative = TRUE)
points(sim$x, L0glm.out, pch = 16, col= 'blue')
graph2pdf(file="..//benchmarks/plots/1_1 single channel deconvolution Poisson.pdf", width=7, height=5)

# some metrics
do.call(cbind, get.metrics(beta = L0glm.out, 
                           X = X, y = y, 
                           family = poisson(identity), 
                           beta_true = sim$beta_true))


## 1.2 Single channel deconvolution assuming Gaussian noise ####
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

# using lambda 1 here would only work if sd_noise would be close to 1
# L0glm.out <- L0glm.fit_new(X, y, 
#                           family = gaussian(identity), 
#                           lambda = 1, 
#                           nonzero.centered = FALSE, # use nonzero centered ridge regression?
#                           nonnegative = TRUE)
# points(sim$x, L0glm.out, pch = 16, col= 'blue')

# using lambda expected to optimise given IC
lambda_preset = preset.lambda(IC = crit, # options aic, bic, gic, ebic, mbic 
                                  y = y, 
                                  X = X, 
                                  family = gaussian(identity)) 
lambda_preset # lambda that would be expected to optimise given IC

L0glm.out <- L0glm.fit_new(X, y, 
                               family = gaussian(identity), 
                               lambda = lambda_preset, 
                               nonzero.centered = FALSE, # use nonzero centered ridge regression?
                               nonnegative = TRUE)
points(sim$x, L0glm.out, pch = 16, col= 'blue')

# some metrics
do.call(cbind, get.metrics(beta = L0glm.out, 
                           X = X, y = y, 
                           family = poisson(identity), 
                           beta_true = sim$beta_true))
graph2pdf(file="..//benchmarks/plots/1_2 single channel deconvolution Gaussian.pdf", width=7, height=5)



# 2. MULTICHANNEL DECONVOLUTION ####
## 2.1 Multichannel deconvolution assuming Poisson noise ####
set.seed(1)
sim <- simulate_spike_train(n = n, 
                            p = p,
                            k = ifelse(round(n/10)<=100, round(n/10), 100),
                            family="poisson", 
                            mean_beta = mean_beta,
                            sd_logbeta = sd_logbeta,
                            multichannel=TRUE)
X <- sim$X
y <- sim$y

# using lambda expected to optimise given IC (lambda = 1.15)
lambda_preset = preset.lambda(IC = crit, # options aic, bic, gic, ebic, mbic 
                                  y = y, 
                                  X = X, 
                                  family = poisson(identity)) 
lambda_preset # lambda = 1.15 would be expected to optimise BIC, lambda = 1.76 would optimise GIC

L0glm.out <- L0glm.fit_new(X, y, 
                               family = poisson(identity), 
                               lambda = lambda_preset, 
                               nonzero.centered = FALSE, # use nonzero centered ridge regression?
                               nonnegative = TRUE)
abline(v=which(rowMin(L0glm.out)!=0), col = 'cyan')
graph2pdf(file="..//benchmarks/plots/2_1 multichannel deconvolution Poisson.pdf", width=7, height=5)

# some metrics
do.call(cbind, get.metrics(beta = L0glm.out, 
                           X = X, y = y, 
                           family = poisson(identity), 
                           beta_true = sim$beta_true))



## 2.2 Multichannel deconvolution assuming Gaussian noise ####
set.seed(1)
sim <- simulate_spike_train(n = n, 
                            p = p,
                            k = ifelse(round(n/10)<=100, round(n/10), 100),
                            family = "gaussian", 
                            mean_beta = mean_beta,
                            sd_logbeta = sd_logbeta,
                            sd_noise = sd_noise, 
                            multichannel = TRUE)
X <- sim$X
y <- sim$y

# using lambda expected to optimise given IC
lambda_preset = preset.lambda(IC = crit, # options aic, bic, gic, ebic, mbic 
                                  y = y, 
                                  X = X,  
                                  family = gaussian(identity)) 
lambda_preset # lambda that would be expected to optimise given IC

L0glm.out <- L0glm.fit_new(X, y, 
                               family = gaussian(identity), 
                               lambda = lambda_preset, 
                               nonzero.centered = FALSE, # use nonzero centered ridge regression?
                               nonnegative = TRUE)
abline(v=which(rowMin(L0glm.out)!=0), col = 'cyan')
graph2pdf(file="..//benchmarks/plots/2_2 multichannel deconvolution Gaussian.pdf", width=7, height=5)

# some metrics
do.call(cbind, get.metrics(beta = L0glm.out, 
                           X = X, y = y, 
                           family = gaussian(identity), 
                           beta_true = sim$beta_true))

