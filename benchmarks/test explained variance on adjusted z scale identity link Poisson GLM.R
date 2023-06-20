# Minimal tests of L0glm package ####

library(remotes)
# remotes::install_github("tomwenseleers/L0glm/L0glm")
library(L0glm)
library(export)

n <- p <- 100 # nr of variables p = nr of observations n in this problem
s <- 0.1 # sparsity, nr of timepoints with a peak 
k <- round(n*s)
mean_beta <- 200
sd_logbeta <- 1
sd_noise <- 5 # SD of noise for Gaussian noise case
crit <- "bic" # choose lambda based on what would be expected to optimise this IC, options in order of increasing regularisation: aic, bic, gic, ebic, mbic

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
beta_true <- sim$beta_true
y_true <- sim$y_true
weights_true <- 1/y_true # true 1/variance weights for identity link Poisson case

# coefficients estimated by L0 glm (approximating GLM that optimises BIC here, i.e. best subset)
coefs <- L0glm.fit_new(X, y, 
                       family = poisson(identity), 
                       lambda = preset.lambda(IC = crit, # options aic, bic, gic, ebic, mbic 
                                              y = y, 
                                              X = X, 
                                              family = poisson(identity)), 
                        nonzero.centered = TRUE, # use nonzero centered ridge regression for iterative adaptive ridge regression updates
                        nonnegative = TRUE) # enforce coefficients to be nonnegative
points(sim$x, coefs, pch = 16, col= 'blue') # estimated coefficients
# graph2pdf(file="..//benchmarks/plots/1_1 single channel deconvolution Poisson.pdf", width=7, height=5)

get.R2 = function (z_obs, # observed data on adjusted response scale
                   z_pred, # fitted data on adjusted response scale
                   weights) { # weights from last IRLS GLM iteration, e.g. 1/variance for identity link GLMs
  good = (weights!=0)&(!is.infinite(weights))
  ss_residual = sum(weights[good] * (z_obs[good] - z_pred[good]) ^ 2)
  ss_total = sum(weights[good] * (z_obs[good] - weighted.mean(z_obs[good], weights[good])) ^ 2)
  r.squared = 1 - ss_residual / ss_total 
  return(r.squared) }

get.R2(z_obs = sim$y,
       z_pred = coefs %*% sim$X,
       weights = weights_true) # this is the true R2 I would like to estimate: 0.9971443

get.adjR2 = function (z_obs, # observed data on adjusted response scale
                      z_pred, # fitted data on adjusted response scale
                      weights, # weights from last IRLS GLM iteration, e.g. 1/variance for identity link GLMs
                      p # nr of variables in the model (nr of nonzero coefficients)
                      ) { 
  good = (weights!=0)&(!is.infinite(weights))
  n = sum(good) # nr of observations
  ss_residual = sum(weights[good] * (z_obs[good] - z_pred[good]) ^ 2)
  ss_total = sum(weights[good] * (z_obs[good] - weighted.mean(z_obs[good], weights[good])) ^ 2)
  df_residual = n-p
  df_total = n-1
  r.squared = 1 - (ss_residual/df_residual) / (ss_total/df_total) 
  return(r.squared) }

get.adjR2(z_obs = sim$y,
          z_pred = as.vector(coefs %*% sim$X),
          weights = 1/as.vector(coefs %*% sim$X),
          p = sum(coefs!=0)) # 0.996848


# some metrics
do.call(cbind, get.metrics(beta = L0glm.out, 
                           X = X, y = y, 
                           family = poisson(identity), 
                           beta_true = sim$beta_true))
