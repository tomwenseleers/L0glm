# Simple test of L0glm package

library(L0glm)
library(glmnet)
devtools::load_all(".")

# single channel for Poisson case

sim <- simulate_spike_train(family="poisson")
X <- sim$X
y <- sim$y

ctrl.fit <- control.fit.gen() # default

L0glm.out <- L0glm.fit_new(X, y, family = gaussian(identity), 
                 lambda = 1, nonnegative = TRUE,
                 control.solve = ctrl.fit)
points(1:100, L0glm.out, pch = 16, col= 'blue')

lambda_preset_gic = preset.lambda(IC="gic", y=sim$y[,which.max(colMeans(abs(sim$y)))], 
                                  X=sim$X, weights=rep(1, nrow(sim$y)), 
                                  family=poisson(identity))



