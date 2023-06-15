# Simple tests of L0glm package ####

library(L0glm)
library(glmnet)
devtools::load_all(".")

# 1. SINGLE CHANNEL DECONVOLUTION ####
## 1.1 Single channel deconvolution assuming Poisson noise ####

sim <- simulate_spike_train(family="poisson")
X <- sim$X
y <- sim$y

# using lambda 1
L0glm.out <- L0glm.fit_new(X, y, family = poisson(identity), 
                 lambda = 1, nonnegative = TRUE)
points(1:100, L0glm.out, pch = 16, col= 'blue')

# using lambda expected to optimise BIC (lambda = 1.15)
lambda_preset_bic = preset.lambda(IC="bic", # options aic, bic, gic, ebic, mbic 
                                  y=sim$y[,which.max(colMeans(abs(sim$y)))], 
                                  X=sim$X, weights=rep(1, nrow(sim$y)), 
                                  family=poisson(identity)) 
lambda_preset_bic # lambda = 1.15 would be expected to optimise BIC

L0glm.out_bic <- L0glm.fit_new(X, y, family = poisson(identity), 
                           lambda = lambda_preset_bic, nonnegative = TRUE,
                           control.solve = ctrl.fit)
points(1:100, L0glm.out_bic, pch = 16, col= 'red')

# some metrics
get.metrics(beta=L0glm.out_bic, X=sim$X, y=sim$y, 
            family=poisson(identity), 
            beta_true=sim$beta_true)


## 1.2 Single channel deconvolution assuming Gaussian noise ####

sim <- simulate_spike_train(family="gaussian", sd_noise=50)
X <- sim$X
y <- sim$y

# using lambda 1 here would only work is sd_noise would be close to 1
L0glm.out <- L0glm.fit_new(X, y, family = gaussian(identity), 
                           lambda = 1, nonnegative = TRUE)
# points(1:100, L0glm.out, pch = 16, col= 'blue')

# using lambda expected to optimise BIC
lambda_preset_bic = preset.lambda(IC="bic", # options aic, bic, gic, ebic, mbic 
                                  y=sim$y[,which.max(colMeans(abs(sim$y)))], 
                                  X=sim$X, weights=rep(1, nrow(sim$y)), 
                                  family=gaussian(identity)) 
lambda_preset_bic # lambda = 1508.727 would be expected to optimise BIC

L0glm.out_bic <- L0glm.fit_new(X, y, family = gaussian(identity), 
                               lambda = lambda_preset_bic, nonnegative = TRUE,
                               control.solve = ctrl.fit)
points(1:100, L0glm.out_bic, pch = 16, col= 'red')

# some metrics
get.metrics(beta=L0glm.out_bic, X=sim$X, y=sim$y, 
            family=poisson(identity), 
            beta_true=sim$beta_true)



# 2. MULTICHANNEL DECONVOLUTION ####
## 2.1 Single channel deconvolution assuming Poisson noise ####

sim <- simulate_spike_train(family="poisson", multichannel=TRUE)
X <- sim$X
y <- sim$y

# using lambda 1
L0glm.out <- L0glm.fit_new(X, y, family = poisson(identity), 
                           lambda = 1, nonnegative = TRUE)
points(1:100, L0glm.out, pch = 16, col= 'magenta')

# using lambda expected to optimise BIC (lambda = 1.15)
lambda_preset_bic = preset.lambda(IC="bic", # options aic, bic, gic, ebic, mbic 
                                  y=sim$y[,which.max(colMeans(abs(sim$y)))], 
                                  X=sim$X, weights=rep(1, nrow(sim$y)), 
                                  family=poisson(identity)) 
lambda_preset_bic # lambda = 1.15 would be expected to optimise BIC

L0glm.out_bic <- L0glm.fit_new(X, y, family = poisson(identity), 
                               lambda = lambda_preset_bic, nonnegative = TRUE,
                               control.solve = ctrl.fit)
points(1:100, L0glm.out_bic, pch = 16, col= 'cyan')

# some metrics
get.metrics(beta=L0glm.out_bic, X=sim$X, y=sim$y, 
            family=poisson(identity), 
            beta_true=sim$beta_true)


## 2.2 Multichannel deconvolution assuming Gaussian noise ####

sim <- simulate_spike_train(family="gaussian", sd_noise=50)
X <- sim$X
y <- sim$y

# using lambda 1 here would only work is sd_noise would be close to 1
L0glm.out <- L0glm.fit_new(X, y, family = gaussian(identity), 
                           lambda = 1, nonnegative = TRUE)
# points(1:100, L0glm.out, pch = 16, col= 'blue')

# using lambda expected to optimise BIC
lambda_preset_bic = preset.lambda(IC="bic", # options aic, bic, gic, ebic, mbic 
                                  y=sim$y[,which.max(colMeans(abs(sim$y)))], 
                                  X=sim$X, weights=rep(1, nrow(sim$y)), 
                                  family=gaussian(identity)) 
lambda_preset_bic # lambda = 1508.727 would be expected to optimise BIC

L0glm.out_bic <- L0glm.fit_new(X, y, family = gaussian(identity), 
                               lambda = lambda_preset_bic, nonnegative = TRUE,
                               control.solve = ctrl.fit)
points(1:100, L0glm.out_bic, pch = 16, col= 'red')

# some metrics
get.metrics(beta=L0glm.out_bic, X=sim$X, y=sim$y, 
            family=poisson(identity), 
            beta_true=sim$beta_true)









