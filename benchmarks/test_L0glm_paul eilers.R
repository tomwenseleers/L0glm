# Minimal tests of L0glm package ####

library(L0glm)
library(export)

n <- p <- 100 # to check performance use e.g. n=p=10000, otherwise use e.g. n=p=100
s <- 0.1 # sparsity, nr of timepoints with a peak
k <- round(n*s)
set.seed(1)

# 1. SINGLE CHANNEL DECONVOLUTION ####
## 1.1 Single channel deconvolution assuming Poisson noise ####
sim <- simulate_spike_train(n = n, 
                            p = p,
                            k = k,
                            family = "poisson", 
                            mean_beta = 10000,
                            sd_logbeta = 1)
X <- sim$X
y <- sim$y

# using lambda 1
L0glm.out <- L0glm.fit_new(X, y, 
                           family = poisson(identity), 
                           lambda = 1, 
                           nonnegative = TRUE)
# points(1:100, L0glm.out, pch = 16, col= 'blue')

# using lambda expected to optimise BIC (lambda = 1.15)
lambda_preset_bic = preset.lambda(IC = "bic", # options aic, bic, gic, ebic, mbic 
                                  y = y, 
                                  X = X, 
                                  family = poisson(identity)) 
lambda_preset_bic # lambda = 1.15 would be expected to optimise BIC

L0glm.out_bic <- L0glm.fit_new(X, y, 
                               family = poisson(identity), 
                               lambda = lambda_preset_bic, 
                               nonnegative = TRUE)
points(sim$x, L0glm.out_bic, pch = 16, col= 'blue')
graph2pdf(file="..//benchmarks/plots/1_1 single channel deconvolution Poisson.pdf", width=7, height=5)

# some metrics
do.call(cbind, get.metrics(beta = L0glm.out_bic, 
                           X = X, y = y, 
                           family = poisson(identity), 
                           beta_true = sim$beta_true))


## 1.2 Single channel deconvolution assuming Gaussian noise ####
sim <- simulate_spike_train(n = n, 
                            p = p,
                            k = k,
                            family = "gaussian", 
                            mean_beta = 10000,
                            sd_logbeta = 1,
                            sd_noise = 10)
X <- sim$X
y <- sim$y

# using lambda 1 here would only work if sd_noise would be close to 1
# L0glm.out <- L0glm.fit_new(X, y, 
#                           family = gaussian(identity), 
#                           lambda = 1, 
#                           nonnegative = TRUE)
# points(1:100, L0glm.out, pch = 16, col= 'blue')

# using lambda expected to optimise BIC
lambda_preset_bic = preset.lambda(IC = "bic", # options aic, bic, gic, ebic, mbic 
                                  y = y, 
                                  X = X, 
                                  family = gaussian(identity)) 
lambda_preset_bic # lambda = 80.81514 would be expected to optimise BIC

L0glm.out_bic <- L0glm.fit_new(X, y, 
                               family = gaussian(identity), 
                               lambda = lambda_preset_bic, 
                               nonnegative = TRUE)
points(sim$x, L0glm.out_bic, pch = 16, col= 'blue')

# some metrics
do.call(cbind, get.metrics(beta = L0glm.out_bic, 
                           X = X, y = y, 
                           family = poisson(identity), 
                           beta_true = sim$beta_true))
graph2pdf(file="..//benchmarks/plots/1_2 single channel deconvolution Gaussian.pdf", width=7, height=5)



# 2. MULTICHANNEL DECONVOLUTION ####
## 2.1 Single channel deconvolution assuming Poisson noise ####
sim <- simulate_spike_train(n = n, 
                            p = p,
                            k = ifelse(round(n/10)<=100, round(n/10), 100),
                            family="poisson", 
                            mean_beta = 10000,
                            sd_logbeta = 1,
                            multichannel=TRUE)
X <- sim$X
y <- sim$y

# using lambda expected to optimise BIC (lambda = 1.15)
lambda_preset_bic = preset.lambda(IC = "bic", # options aic, bic, gic, ebic, mbic 
                                  y = y, 
                                  X = X, 
                                  family = poisson(identity)) 
lambda_preset_bic # lambda = 1.15 would be expected to optimise BIC

L0glm.out_bic <- L0glm.fit_new(X, y, 
                               family = poisson(identity), 
                               lambda = lambda_preset_bic, 
                               nonnegative = TRUE)
abline(v=which(rowMin(L0glm.out_bic)!=0), col = 'cyan')
graph2pdf(file="..//benchmarks/plots/2_1 multichannel deconvolution Poisson.pdf", width=7, height=5)


# some metrics
do.call(cbind, get.metrics(beta = L0glm.out_bic, 
                           X = X, y = y, 
                           family = poisson(identity), 
                           beta_true = sim$beta_true))



## 2.2 Multichannel deconvolution assuming Gaussian noise ####
sim <- simulate_spike_train(n = n, 
                            p = p,
                            k = ifelse(round(n/10)<=100, round(n/10), 100),
                            family = "gaussian", 
                            mean_beta = 10000,
                            sd_logbeta = 1,
                            sd_noise = 10, 
                            multichannel = TRUE)
X <- sim$X
y <- sim$y

# using lambda expected to optimise BIC
lambda_preset_bic = preset.lambda(IC = "bic", # options aic, bic, gic, ebic, mbic 
                                  y = y, 
                                  X = X,  
                                  family = gaussian(identity)) 
lambda_preset_bic # lambda = 70.14565 would be expected to optimise BIC

L0glm.out_bic <- L0glm.fit_new(X, y, 
                               family = gaussian(identity), 
                               lambda = lambda_preset_bic, 
                               nonnegative = TRUE)
abline(v=which(rowMin(L0glm.out_bic)!=0), col = 'cyan')
graph2pdf(file="..//benchmarks/plots/2_2 multichannel deconvolution Gaussian.pdf", width=7, height=5)

# some metrics
do.call(cbind, get.metrics(beta = L0glm.out_bic, 
                           X = X, y = y, 
                           family = gaussian(identity), 
                           beta_true = sim$beta_true))
