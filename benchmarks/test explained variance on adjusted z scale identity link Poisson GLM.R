# Tests R2 for GLMs with best subset selection ####

library(remotes)
# remotes::install_github("tomwenseleers/L0glm/L0glm")
library(L0glm)
library(export)

## Single channel deconvolution of a blurred spike train with Poisson noise ####

# simulated dataset: blurred superimposed spike train (spike train convolution with given point spread function)
# with Poisson noise on it
set.seed(1)
sim <- simulate_spike_train(n = 100, # nr of observations 
                            p = p, # nr of variables
                            k = 10, # true nr of nonzero coefficients
                            family = "poisson", 
                            mean_beta = 30, # true nonzero coefficients have lognormal distribution
                            sd_logbeta = 1)
X <- sim$X # design matrix = temporally shifted copies of known point spread function
y <- sim$y # observed signal (with Poisson noise)
beta_true <- sim$beta_true # true coefficients used in simulation
y_true <- sim$y_true # true underlying signal without Poisson noise
weights_true <- 1/y_true # true 1/variance weights for identity link Poisson case

# coefficients estimated by L0 glm (approximating GLM that optimises BIC, i.e. best subset)
coefs <- L0glm.fit_new(X, y, 
                       family = poisson(identity), 
                       lambda = preset.lambda(IC = crit, # use lambda that optimises given IC, options aic, bic, gic, ebic, mbic
                                              y = y,
                                              X = X,
                                              family = poisson(identity)),
                        # lambda = 100,
                        nonzero.centered = TRUE, # use nonzero centered ridge regression for iterative adaptive ridge regression updates to approximate best subset
                        nonnegative = TRUE)$coefficients # enforce coefficients to be nonnegative
# fitted values
y_pred <- as.vector(coefs %*% sim$X)
# estimated 1/variance weights
weights_est <- 1/y_pred

points(sim$x, coefs, pch = 16, col= 'blue') # estimated coefficients
graph2png(file="..//benchmarks/plots/single channel deconvolution Poisson.png", width=7, height=5, dpi=150)

# define functions to calculate weighted R2 and adjusted weighted R2 :
get.R2 = function (z_obs, # observed data on adjusted response scale
                   z_pred, # fitted data on adjusted response scale
                   weights, # weights from last IRLS GLM iteration, e.g. 1/variance for identity link GLMs
                   intercept=FALSE) { 
  good = as.vector(weights<1) # to remove infinite & unreliable weights
  ss_residual = sum(weights[good] * (z_obs[good] - z_pred[good]) ^ 2)
  if (intercept) { wmean = weighted.mean(z_obs[good], weights[good]) } else { wmean = 0 }
  ss_total = sum(weights[good] * (z_obs[good] - wmean) ^ 2)
  r.squared = 1 - ss_residual / ss_total 
  return(r.squared) }

get.adjR2 = function (z_obs, # observed data on adjusted response scale
                      z_pred, # fitted data on adjusted response scale
                      weights, # weights from last IRLS GLM iteration, e.g. 1/variance for identity link GLMs
                      intercept=FALSE,
                      p # nr of variables in the model (nr of nonzero coefficients) (incl intercept)
                      ) { 
  good = as.vector(weights<1) # to remove infinite & unreliable weights
  n = sum(good) # nr of observations
  ss_residual = sum(weights[good] * (z_obs[good] - z_pred[good]) ^ 2)
  if (intercept) { wmean = weighted.mean(z_obs[good], weights[good]) } else { wmean = 0 }
  ss_total = sum(weights[good] * (z_obs[good] - weighted.mean(z_obs[good], weights[good])) ^ 2)
  df_residual = n-p-1
  df_total = n-1
  r.squared = 1 - (ss_residual/df_residual) / (ss_total/df_total) 
  return(r.squared) }

# adjusted weighted R2 between predicted signal & observed signal
get.adjR2(z_obs = y,
          z_pred = y_pred,
          weights = weights_est,
          p = sum(coefs!=0),
          intercept = FALSE) 
# 0.96794

# Ideally I think what I would like to have is an unbiased estimator of the
# weighted R2 between the predicted signal & the true underlying signal without noise :

# weighted R2 between predicted signal & true underlying signal without noise, 
# using true 1/variance weights from simulation
# Is the adjusted weighted R2 between predicted signal & observed signal an unbiased
# estimator of this? (for a wide range of lambda values it seems it is in the right
# ballpark, but for this particular lambda & when it gets really close to 1 it is 
# higher than the estimated value above, I presume because for sim$y we only measure 1
# realisation of the Poisson noise & the discrete nature of the noise?)
get.R2(z_obs = y_true, 
       z_pred = y_pred,
       weights = weights_true,
       intercept = FALSE
       ) 
# 0.9967505


# some metrics
do.call(cbind, get.metrics(beta = coefs, 
                           X = X, y = y, 
                           family = poisson(identity), 
                           beta_true = beta_true))
