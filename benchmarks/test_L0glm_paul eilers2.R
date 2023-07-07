# Minimal tests of L0glm package ####

library(remotes)
# remotes::install_github("tomwenseleers/L0glm/L0glm")
library(L0glm)
library(export)
# devtools::load_all(".")

# 1. SINGLE CHANNEL DECONVOLUTION ####

n <- p <- 100 # 10000 # to check performance use e.g. n=p=10000, otherwise use e.g. n=p=100-300
s <- 0.1 # sparsity, nr of timepoints with a peak 
k <- round(n*s)
mean_beta <- 2000
sd_logbeta <- 1
sd_noise <- 5 # SD of noise for Gaussian noise case
crit <- "gic" # choose lambda based on what would be expected to optimise this IC
# options in order of increasing regularisation: aic, bic, gic, ebic, mbic

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
dispersion <- 1
family <- poisson(identity)

# using lambda 1
system.time(L0glm.out <- L0glm.fit(X, y, 
                           family = poisson(identity), 
                           lambda = 1, 
                           nonzero.centered = FALSE, # use nonzero centered ridge regression?
                           # ridge.center = start, # ridge.center
                           nonnegative = TRUE,
                           warmupiters = 1,
                           # warmupiters.nzc = 1,
                           # start = rep(10, ncol(sim$X)), # coefficient starting values
                           solver = "eigen", # note: solver osqp better than eigen for very small lambda values)
                           method = 7L # 1 & 2 = QR, 3 & 4 = LSCG, 5 & 6 = CG, 7 & 8 = Cholesky
)$coefficients) # 0.25s for n=p=10000 with eigen solver methods 7 & 8 (Cholesky), 0.75s with methods 5 & 6 (CG), 1.9s with methods 3, 4 (LSCG), 219s with methods 1 & 2 (QR), 5s for n=p=10000 with osqp solver, 2.27s for n=p=10000 with glmnet solver
# points(sim$x, L0glm.out, pch = 16, col= 'blue')

# using lambda expected to optimise given IC
lambda_preset = preset.lambda(IC = crit, # options aic, bic, gic, ebic, mbic 
                                  y = y, 
                                  X = X, 
                                  family = family) 
lambda_preset # lambda = 1.15 would be expected to optimise BIC with n=p=100, lambda = 1.76 would optimise GIC

system.time(L0glm.out_preset <- L0glm.fit(X, y, 
                               family = family, 
                               lambda = lambda_preset, 
                               nonzero.centered = FALSE, # use nonzero centered ridge regression?
                               # ridge.center = start, # ridge.center
                               nonnegative = TRUE,
                               warmupiters = 1,
                               # warmupiters.nzc = 1,
                               # start = rep(1, ncol(sim$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                               solver = "eigen", # note: solver osqp & glmnet better than eigen for very small lambda values)
                               method = 7L, # 0.27s for n=p=10000 with eigen solver methods 7 & 8 (Cholesky), 0.61s with methods 5 & 6 (CG), 1.27s with methods 3, 4 (LSCG), 219s with methods 1 & 2 (QR), 4.8s for n=p=10000 with osqp solver, 2s for n=p=10000 with glmnet solver
                               thresh = 1E-5
)) # 0.02s for n=p=100, 0.25s for n=p=10000
L0glm.out_preset_coefs <- L0glm.out_preset$coefficients
sum(L0glm.out_preset_coefs!=0)  # 1004
points(sim$x, L0glm.out_preset_coefs, pch = 16, col= 'blue')
graph2pdf(file="..//benchmarks/plots/1_1 single channel deconvolution Poisson.pdf", width=7, height=5)

mincoef <- min(abs(L0glm.out_preset_coefs)[L0glm.out_preset_coefs!=0]) # nonzero coefficient with smallest absolute value
do.call(cbind, get.metrics(beta = L0glm.out_preset_coefs,  # some metrics
                           X = X, y = y, 
                           family = family, 
                           beta_true = sim$beta_true))

plot(x=L0glm.out_preset_coefs, y=sim$beta_true, pch=16, col="blue")


# coefficients & metrics over a range of lambda values
# calculate max lambda where all coefs would be penalized to 0
maxlam <- maxlambda(X=X, weights=L0glm.out_preset$weights, 
                    z=y, no.pen=NULL, coefs=L0glm.out_preset$coefficients,
                    family=family)
maxlam # 

lambdas <- 10^seq(round(log10(lambda_preset))-1, log10(maxlam)+1, by=0.05) #
coefmat <- Matrix(0, nrow=length(lambdas), ncol=ncol(sim$X))
warmstart <- 0.001 # % by which L0glm.out_preset will be used as start, here used lowest value that still gave stable coef paths
system.time(metrics <- data.frame(do.call(rbind, lapply(1:length(lambdas), function (i) { 
  lambda <- lambdas[[i]]
  start <- L0glm.out_preset_coefs*warmstart+start_coefs(X, family)*(1-warmstart)
  # start[L0glm.out_preset_coefs==0] <- 0 # slightly faster & stable with this
  L0glm.out <- L0glm.fit(X, y, 
                             family = family, 
                             # family = gaussian(identity),
                             # weights = 1/(sim$y+0.1),
                             lambda = lambda, 
                             nonzero.centered = FALSE, # use nonzero centered ridge regression?
                             ridge.center = NULL, 
                             nonnegative = TRUE,
                             warmupiters = 1,
                             warmupiters.nzc = 5,
                             start = start, # 
                             solver = "eigen", # eigen, glmnet or osqp, # note: constrained solvers glmnet & osqp better than eigen for very small lambda values
                             method = 7L, 
                             # gamma = 2, 
                             thresh = 1E-5,
                             relaxed = FALSE
                             ) 
  L0glm.out_coefs = L0glm.out$coefficients
  coefmat[i,] <<- as.vector(L0glm.out_coefs)
  return( unlist(c(lambda=lambda, get.metrics(beta = L0glm.out_coefs,
                                                       X = X, y = y,
                                                       family = family, 
                                                       beta_true = sim$beta_true) ))

          )
  } )))) # 4s for n=p=100 & 29s for n=p=10000 with eigen solver methods 7 & 8 (Cholesky), 165s for n=p=10000 with glmnet solver, 455s for n=p=10000 with osqp solver
matplot(x=log10(metrics$lambda), y=as.matrix(coefmat), type="l", lty=1,
        xlab="Log10(lambda)", ylab="Coefficient estimates", main="Regularisation path L0 Poisson GLM")
abline(v=log10(maxlam), col="red") # max lambda where all coefs would be penalized to 0
library(scales)
abline(h=as.vector(sim$beta_true), col=alpha(1:6,0.2))

plot(log10(metrics$lambda), metrics$deviance, type="l")
metrics$lambda[which.min(metrics$deviance)] # 1 for n=p=10000, 0.1 for n=p=100 # NOTE this is deviance relative to observed signal, not deviance relative to ground truth though
plot(log10(metrics$lambda), metrics$deviance_true, type="l")
metrics$lambda[which(metrics$deviance_true==min(metrics$deviance_true))] # XX for n=p=10000, 10 for n=p=100 # NOTE this is deviance relative to observed signal, not deviance relative to ground truth though
plot(log10(metrics$lambda), metrics$loglik, type="l")
metrics$lambda[which(metrics$loglik==max(metrics$loglik))] # 1 for n=p=10000, 0.1 for n=p=100
plot(log10(metrics$lambda), metrics$MD_coefs, type="l")
plot(log10(metrics$lambda), metrics$MAD_coefs, type="l")
min(metrics$lambda[which(round(metrics$MAD_coefs,0)==min(round(metrics$MAD_coefs,0)))]) # 2.238721 for n=p=10000, 0.5623413 for n=p=100 
plot(log10(metrics$lambda), metrics$sMAPE_coefs, type="l") # good one, optimal for lambda between 1 & 10
min(metrics$lambda[which(round(metrics$sMAPE_coefs,2)==min(round(metrics$sMAPE_coefs,2)))]) # 11.22018 for n=p=10000, 0.5623413 for n=p=100
plot(log10(metrics$lambda), metrics$bias_coefs, type="l")
metrics$lambda[which(abs(metrics$bias_coefs-1)==min(abs(metrics$bias_coefs-1)))] # 
plot(log10(metrics$lambda), metrics$R2_coefs, type="l")
plot(log10(metrics$lambda), metrics$estimated_over_true_support_size, type="l") # good one, optimal for lambda between 1 & 10
range(metrics$lambda[which(abs(metrics$estimated_over_true_support_size-1)==min(abs(metrics$estimated_over_true_support_size-1)))]) # 2.238721 for n=p=10000, 0.5623413-25.1188643 for n=p=100
plot(log10(metrics$lambda), metrics$TP, type="l")
plot(log10(metrics$lambda), metrics$FP, type="l")
plot(log10(metrics$lambda), metrics$TN, type="l")
plot(log10(metrics$lambda), metrics$FN, type="l")
plot(log10(metrics$lambda), metrics$precision, type="l")
plot(log10(metrics$lambda), metrics$recall, type="l")
plot(log10(metrics$lambda), metrics$specificity, type="l")
plot(log10(metrics$lambda), metrics$accuracy, type="l") # good one, optimal for lambda between 1 & 10
min(metrics$lambda[which(metrics$accuracy==max(metrics$accuracy))]) # 12.58925 for n=p=10000, 0.5623413 for n=p=100
range(metrics$lambda[which(metrics$accuracy==max(metrics$accuracy))]) # 12.58925 12.58925 for n=p=10000, 0.5623413-25.1188643 for n=p=100
plot(log10(metrics$lambda), metrics$F1_score, type="l")
metrics$lambda[which.max(metrics$F1_score)] # 12.58925 for n=p=10000, 0.5623413 for n=p=100
plot(log10(metrics$lambda), metrics$aic, type="l")
min(metrics$lambda[which(metrics$aic==min(metrics$aic))]) # 1 for n=p=10000, 0.5623413 for n=p=10000 (optimal at same value as where loglik is optimal)
plot(log10(metrics$lambda), metrics$bic, type="l")
max(metrics$lambda[which(metrics$bic==min(metrics$bic))]) # 1.584893 for n=p=10000, 0.5623413 for n=p=100
plot(log10(metrics$lambda), metrics$ebic, type="l")
max(metrics$lambda[which(metrics$ebic==min(metrics$ebic))]) # 1.584893 for n=p=10000, 0.5623413 for n=p=100
plot(log10(metrics$lambda), metrics$gic, type="l")
max(metrics$lambda[which(metrics$gic==min(metrics$gic))]) # 1.584893 for n=p=10000, 0.5623413 for n=p=100
abline(v=log10(max(metrics$lambda[which(metrics$gic==min(metrics$gic))])), col="red") 
plot(log10(metrics$lambda), metrics$mbic, type="l")
min(metrics$lambda[which(metrics$mbic==min(metrics$mbic))]) # 1.995262 for n=p=10000, 0.5623413 for n=p=100


# lambda values that are expected to optimise each IC (with start = rep(1, ncol(X)))
# note: using warm starts in the L0glm.fit function would shift the
# empirically observed optimal values
preset.lambda(IC = "aic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(sim$y), 
              X = Matrix(sim$X, sparse=TRUE), 
              family = poisson(identity),
              dispersion = dispersion) # 0.5 for n=p=100 or n=p=10000
preset.lambda(IC = "bic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(sim$y), 
              X = Matrix(sim$X, sparse=TRUE), 
              family = poisson(identity),
              dispersion = dispersion) # 1.151293 for n=p=100, 2.302585 for n=p=10000
preset.lambda(IC = "ebic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(sim$y), 
              X = Matrix(sim$X, sparse=TRUE), 
              family = poisson(identity),
              dispersion = dispersion) # 2.302585 for n=p=100, 4.60517 for n=p=10000
preset.lambda(IC = "gic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(sim$y), 
              X = Matrix(sim$X, sparse=TRUE), 
              family = poisson(identity),
              dispersion = dispersion) # 1.758231 for n=p=100, 5.112491 for n=p=10000
preset.lambda(IC = "mbic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(sim$y), 
              X = Matrix(sim$X, sparse=TRUE), 
              family = poisson(identity),
              dispersion = dispersion) # 2.76073 for n=p=100, 6.214608 for n=p=10000




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
family <- gaussian(identity)

# using lambda expected to optimise given IC
lambda_preset = preset.lambda(IC = crit, # options aic, bic, gic, ebic, mbic 
                                  y = y, 
                                  X = X, 
                                  family = family) 
lambda_preset # lambda that would be expected to optimise given IC

L0glm.out_preset <- L0glm.fit(X, y, 
                               family = family, 
                               lambda = lambda_preset, 
                               nonzero.centered = FALSE, # use nonzero centered ridge regression?
                               nonnegative = TRUE,
                               warmupiters = 1,
                               start = NULL, # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                               solver = "eigen", # note: solver osqp better than eigen for very small lambda values))
                               method = 7L # Cholesky
                              )
L0glm.out_preset_coefs = L0glm.out_preset$coefficients
points(sim$x, L0glm.out_preset_coefs, pch = 16, col= 'blue')

# some metrics
do.call(cbind, get.metrics(beta = L0glm.out_preset_coefs, 
                           X = X, y = y, 
                           family = family, 
                           beta_true = sim$beta_true))
graph2pdf(file="..//benchmarks/plots/1_2 single channel deconvolution Gaussian.pdf", width=7, height=5)



# coefficients & metrics over a range of lambda values
# calculate max lambda where all coefs would be penalized to 0
maxlam <- maxlambda(X=X, weights=L0glm.out_preset$weights, 
                    z=y, no.pen=NULL, coefs=L0glm.out_preset$coefficients,
                    family=family)
maxlam #

lambdas <- 10^seq(round(log10(lambda_preset))-1, log10(maxlam)+1, by=0.05)
coefmat <- Matrix(0, nrow=length(lambdas), ncol=ncol(sim$X))
warmstart <- 0.1 # % by which L0glm.out_preset will be used as start, lowest value that led to stable coefficient paths
metrics <- data.frame(do.call(rbind, lapply(1:length(lambdas), function (i) { 
  lambda <- lambdas[[i]]
  start <- L0glm.out_preset_coefs*warmstart+start_coefs(X, family)*(1-warmstart)
  # start[L0glm.out_preset_coefs==0] <- 0
  L0glm.out <- L0glm.fit(X, y, 
                             family = gaussian(identity), 
                             lambda = lambda, 
                             nonzero.centered = FALSE, # use nonzero centered ridge regression?
                             # ridge.center = NULL,
                             nonnegative = TRUE,
                             warmupiters = 1,
                             warmupiters.nzc = 5,
                             start = start, # rep(0.1, ncol(sim$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                             solver = "eigen", # note: solver osqp better than eigen for very small lambda values))
                             method = 7L,
                             relaxed = FALSE,
                             thresh = 1E-5,
)$coefficients
  coefmat[i,] <<- as.vector(L0glm.out)
  return( unlist(c(lambda=lambda, get.metrics(beta = L0glm.out, 
                                              X = X, y = y, 
                                              family = gaussian(identity), 
                                              beta_true = sim$beta_true))) )
} )))
matplot(x=log10(metrics$lambda), y=as.matrix(coefmat), type="l", lty=1,
        xlab="Log10(lambda)", ylab="Coefficient estimates", main="Regularisation path L0 Gaussian GLM")
abline(v=log10(maxlam), col="red") # max lambda where all coefs would be penalized to 0
library(scales)
abline(h=as.vector(sim$beta_true), col=alpha(1:6,0.2))

plot(log10(metrics$lambda), metrics$deviance, type="l")
metrics$lambda[which.min(metrics$deviance)] # 1 # NOTE this is deviance relative to observed signal, not deviance relative to ground truth though
plot(log10(metrics$lambda), metrics$deviance_true, type="l")
metrics$lambda[which(metrics$deviance_true==min(metrics$deviance_true))] # 100 # deviance relative to ground truth
plot(log10(metrics$lambda), metrics$loglik, type="l")
abline(v=log10(maxlam), col="red") # max lambda where all coefs would be penalized to 0
metrics$lambda[which.max(metrics$loglik)] # 4.466836
plot(log10(metrics$lambda), metrics$MD_coefs, type="l")
plot(log10(metrics$lambda), metrics$MAD_coefs, type="l")
metrics$lambda[which.min(metrics$MAD_coefs)] # 4.466836
plot(log10(metrics$lambda), metrics$sMAPE_coefs, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(metrics$sMAPE_coefs)] # 4.466836
plot(log10(metrics$lambda), metrics$bias_coefs, type="l")
metrics$lambda[which.min(abs(metrics$bias_coefs-1))] # 1.778279
plot(log10(metrics$lambda), metrics$R2_coefs, type="l")
plot(log10(metrics$lambda), metrics$estimated_over_true_support_size, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(abs(metrics$estimated_over_true_support_size-1))] # 4.466836
plot(log10(metrics$lambda), metrics$TP, type="l")
plot(log10(metrics$lambda), metrics$FP, type="l")
plot(log10(metrics$lambda), metrics$TN, type="l")
plot(log10(metrics$lambda), metrics$FN, type="l")
plot(log10(metrics$lambda), metrics$precision, type="l")
plot(log10(metrics$lambda), metrics$recall, type="l")
plot(log10(metrics$lambda), metrics$specificity, type="l")
plot(log10(metrics$lambda), metrics$accuracy, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.max(metrics$accuracy)] # 4.466836
plot(log10(metrics$lambda), metrics$F1_score, type="l")
metrics$lambda[which.max(metrics$F1_score)] # 
plot(log10(metrics$lambda), metrics$aic, type="l")
metrics$lambda[which.min(metrics$aic)] # 4.466836
plot(log10(metrics$lambda), metrics$bic, type="l")
metrics$lambda[which.min(metrics$bic)] # 4.466836
plot(log10(metrics$lambda), metrics$ebic, type="l")
metrics$lambda[which.min(metrics$ebic)] # 4.466836
plot(log10(metrics$lambda), metrics$gic, type="l")
metrics$lambda[which.min(metrics$gic)] # 4.466836
abline(v=log10(metrics$lambda[which.min(metrics$gic)]), col="blue") 
plot(log10(metrics$lambda), metrics$mbic, type="l")
metrics$lambda[which(metrics$mbic==min(metrics$mbic))] # 4.466836

# lambda values that a priori are expected to optimise each IC (with start = NULL)
# note: using warm starts in the L0glm.fit function would shift the
# empirically observed optimal values
preset.lambda(IC = "aic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(sim$y), 
              X = Matrix(sim$X, sparse=TRUE), 
              family = gaussian(identity),
              dispersion = dispersion) # 0.5
preset.lambda(IC = "bic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(sim$y), 
              X = Matrix(sim$X, sparse=TRUE), 
              family = gaussian(identity),
              dispersion = dispersion) # 1.151293
preset.lambda(IC = "ebic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(sim$y), 
              X = Matrix(sim$X, sparse=TRUE), 
              family = gaussian(identity),
              dispersion = dispersion) # 2.302585
preset.lambda(IC = "gic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(sim$y), 
              X = Matrix(sim$X, sparse=TRUE), 
              family = gaussian(identity),
              dispersion = dispersion) # 1.758231
preset.lambda(IC = "mbic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(sim$y), 
              X = Matrix(sim$X, sparse=TRUE), 
              family = gaussian(identity),
              dispersion = dispersion) # 2.76073




# 2. MULTICHANNEL DECONVOLUTION ####

n <- p <- 500 # 10000 # to check performance use e.g. n=p=10000, otherwise use e.g. n=p=500
s <- 0.1 # sparsity, nr of timepoints with a peak 
k <- round(n*s)
mean_beta <- 2000
sd_logbeta <- 1
sd_noise <- 5 # SD of noise for Gaussian noise case
crit <- "gic" # choose lambda based on what would be expected to optimise this IC
# options in order of increasing regularisation: aic, bic, gic, ebic, mbic

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
family <- poisson(identity)
m <- ncol(y) # nr of tasks/channels

# using lambda expected to optimise given IC
lambda_preset = preset.lambda(IC = crit, # options aic, bic, gic, ebic, mbic 
                                  y = y, 
                                  X = X, 
                                  family = family) 
lambda_preset 

system.time(L0glm.out_preset <- L0glm.fit(X, y, 
                               family = family, 
                               lambda = lambda_preset, 
                               nonzero.centered = FALSE, # use nonzero centered ridge regression?
                               # ridge.center = NULL,
                               grouped = TRUE,
                               relaxed = TRUE, # NOTE: CHECK setting relaxed to TRUE does not seem to improve solution quality
                               nonnegative = TRUE,
                               warmupiters = 1,
                               parallellize = TRUE,
                               start = NULL, # rep(10, ncol(sim$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                               solver = "eigen", # note: solver osqp better than eigen for very small lambda values))
                               method = 7L,
                               thresh = 1E-5
)) # 10.46s for n=p=100, 12s for n=p=500, 30s for n=p=10000
L0glm.out_preset_coefs = L0glm.out_preset$coefficients
abline(v=which(rowMin(L0glm.out_preset_coefs)!=0), col = 'cyan')
graph2pdf(file="..//benchmarks/plots/2_1 multichannel deconvolution Poisson.pdf", width=7, height=5)

# some metrics
metrics_L0glm = do.call(cbind, get.metrics(beta = L0glm.out_preset_coefs, 
                           X = X, y = y, 
                           family = poisson(identity), 
                           beta_true = sim$beta_true))
metrics_L0glm
#        p   n   m deviance deviance_true loglik     MD_coefs MAD_coefs sMAPE_coefs bias_coefs  R2_coefs
# [1,] 500 500 157    1e+11      2.71e+13   -Inf -0.001275467 0.1924362  0.03985753   0.999548 0.9950662
#      support_size_true support_size_estimated estimated_over_true_support_size TP FP  TN FN precision recall
# [1,]                50                     50                                1 48  2 448  2      0.96   0.96
#      specificity accuracy F1_score gic ebic bic aic mbic
# [1,]   0.9955556    0.992     0.96 Inf  Inf Inf Inf  Inf


# coefficients & metrics over a range of lambda values
# calculate max lambda where all coefs would be penalized to 0
maxlam <- maxlambda(X=X, weights=L0glm.out_preset$weights, 
                    z=y, no.pen=NULL, coefs=L0glm.out_preset$coefficients,
                    family=family)
maxlam # 

lambdas <- 10^seq(round(log10(lambda_preset))-1, log10(maxlam), by=0.2) #
coefmat <- Matrix(0, nrow=length(lambdas), ncol=ncol(sim$X)) # matrix with avg coefficients over tasks/channels
system.time(metrics <- data.frame(do.call(rbind, lapply(1:length(lambdas), function (i) { 
  lambda <- lambdas[[i]] # for i=8 gives error at L0glm.R#742 number of rows of result is not a multiple of vector length (arg 8)
  L0glm.out <- L0glm.fit(X, y, 
                             family = poisson(identity), 
                             lambda = lambda, 
                             nonzero.centered = FALSE, # use nonzero centered ridge regression?
                             grouped = TRUE,
                             relaxed = TRUE, 
                             nonnegative = TRUE,
                             warmupiters = 1,                               
                             start = NULL, # rep(1, ncol(sim$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                             solver = "eigen", # note: solver osqp better than eigen for very small lambda values))
                             method = 7L,
                             parallellize = TRUE,
                             thresh = 1E-5
  )$coefficients
  coefmat[i,] <<- rowMeans(L0glm.out) # avg coefficients across tasks/channels
  return( unlist(c(lambda=lambda, get.metrics(beta = L0glm.out, 
                                              X = X, y = y, 
                                              family = poisson(identity), 
                                              beta_true = sim$beta_true) )) )
} ) ))) # 263s for n=p=100, 

matplot(x=log10(metrics$lambda), y=as.matrix(coefmat), type="l", lty=1,
        xlab="Log10(lambda)", ylab="Coefficient estimates", main="Regularisation path Group L0 Poisson GLM")

plot(log10(metrics$lambda), metrics$deviance, type="l")
metrics$lambda[which.min(metrics$deviance)] # 0.1584893 # NOTE this is deviance relative to observed signal, not deviance relative to ground truth though
plot(log10(metrics$lambda), metrics$deviance_true, type="l")
metrics$lambda[which(metrics$deviance_true==min(metrics$deviance_true))] # 0.1 # deviance relative to ground truth
plot(log10(metrics$lambda), metrics$loglik, type="l")
metrics$lambda[which.max(metrics$loglik)] # 0.1584893
plot(log10(metrics$lambda), metrics$MD_coefs, type="l")
plot(log10(metrics$lambda), metrics$MAD_coefs, type="l")
metrics$lambda[which.min(metrics$MAD_coefs)] # 0.2511886
plot(log10(metrics$lambda), metrics$sMAPE_coefs, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(metrics$sMAPE_coefs)] # 0.2511886
plot(log10(metrics$lambda), metrics$bias_coefs, type="l")
metrics$lambda[which.min(abs(metrics$bias_coefs-1))] # 0.2511886
plot(log10(metrics$lambda), metrics$R2_coefs, type="l")
plot(log10(metrics$lambda), metrics$estimated_over_true_support_size, type="l") # good one, optimal for lambda between 1 & 10
range(metrics$lambda[which(abs(metrics$estimated_over_true_support_size-1)==min(abs(metrics$estimated_over_true_support_size-1)))]) # 0.2511886 1.0000000
plot(log10(metrics$lambda), metrics$TP, type="l")
plot(log10(metrics$lambda), metrics$FP, type="l")
plot(log10(metrics$lambda), metrics$TN, type="l")
plot(log10(metrics$lambda), metrics$FN, type="l")
plot(log10(metrics$lambda), metrics$precision, type="l")
plot(log10(metrics$lambda), metrics$recall, type="l")
plot(log10(metrics$lambda), metrics$specificity, type="l")
plot(log10(metrics$lambda), metrics$accuracy, type="l") # good one
metrics$lambda[which.max(metrics$accuracy)] # 0.2511886
plot(log10(metrics$lambda), metrics$F1_score, type="l")
metrics$lambda[which.max(metrics$F1_score)] # 0.2511886
plot(log10(metrics$lambda), metrics$aic, type="l")
metrics$lambda[which(metrics$aic==min(metrics$aic))] # 0.1584893
plot(log10(metrics$lambda), metrics$bic, type="l")
metrics$lambda[which.min(metrics$bic)] # 0.1584893
plot(log10(metrics$lambda), metrics$ebic, type="l")
metrics$lambda[which.min(metrics$ebic)] # 0.1584893
plot(log10(metrics$lambda), metrics$gic, type="l")
metrics$lambda[which.min(metrics$gic)] # 0.1584893
abline(v=log10(metrics$lambda[which.min(metrics$gic)]), col="blue") 
plot(log10(metrics$lambda), metrics$mbic, type="l")
metrics$lambda[which(metrics$mbic==min(metrics$mbic))] # 0.1584893

# these lambda values are ca 3-17x lower than the ones that might be expected to optimise these IC
preset.lambda(IC = "aic", # options aic, bic, gic, ebic, mbic 
              y = y, 
              X = X, 
              family = family) # 0.5
preset.lambda(IC = "gic", # options aic, bic, gic, ebic, mbic 
              y = y, 
              X = X, 
              family = family) # 1.758231
preset.lambda(IC = "mbic", # options aic, bic, gic, ebic, mbic 
              y = y, 
              X = X, 
              family = family) # 2.76073


### 2.2.1 Comparison with abess relaxed group L0 weighted mgaussian best subset glm #### 

# This was our L0 glm metrics:
metrics_L0glm = do.call(cbind, get.metrics(beta = L0glm.out_preset_coefs, 
                                           X = X, y = y, 
                                           family = family, 
                                           beta_true = sim$beta_true))
metrics_L0glm
#        p   n   m    deviance deviance_true loglik     MD_coefs MAD_coefs sMAPE_coefs bias_coefs  R2_coefs
# [1,] 100 100 124 40000006978      6.28e+12   -Inf -0.005259792 0.1925821  0.04782878   0.995435 0.9960168
#      support_size_true support_size_estimated estimated_over_true_support_size TP FP TN FN precision recall
# [1,]                10                      9                              0.9  9  0 90  1         1    0.9
#      specificity accuracy  F1_score gic ebic bic aic mbic
# [1,]           1     0.99 0.9473684 Inf  Inf Inf Inf  Inf
# TODO fix Inf values in metrics

# abess mgaussian best subset fit plus metrics:
# best subset multitask learning using family="mgaussian" using abess
library(abess)
system.time(abess_fit <- abess(x = X, 
                               y = y, 
                               weights = 1/(y+0.1),  # approx 1/variance weights to approximate family="poisson", which is not supported
                               family = "mgaussian",
                               tune.path = "sequence", 
                               support.size = c(1:(p-1)),
                               lambda=0, 
                               warm.start = TRUE,
                               tune.type = "gic")) # 0.28s, or cv or ebic or bic or aic
# note: currently no support for nonnegativity constraints
plot(abess_fit, type="tune")  # optimal support size 
beta_abess = as.matrix(extract(abess_fit)$beta) # coefficient matrix for best subset
library(qlcMatrix)
image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=beta_abess^0.1, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="abess multitask weighted mgaussian (red=true peaks)")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim$X))[as.vector(rowMin(beta_abess)!=0)], col="cyan")
sum(rowMax(sim$beta_true)!=0) # 10 true peaks
sum(rowMin(beta_abess)!=0) # 8 peaks detected
# note: for larger p, e.g. 500, solution quality quickly becomes a lot worse
# see https://github.com/abess-team/abess/issues/510

metrics_abess = do.call(cbind, get.metrics(beta = beta_abess, 
                                           X = X, y = y, 
                                           family = poisson(identity), 
                                           beta_true = sim$beta_true))
metrics_abess
#        p   n   m deviance deviance_true loglik    MD_coefs MAD_coefs sMAPE_coefs bias_coefs  R2_coefs
# [1,] 100 100 124 2.05e+12     1.034e+13    NaN -0.04160188 0.4694801   0.1366435  0.9838769 0.9650847
#      support_size_true support_size_estimated estimated_over_true_support_size TP FP TN FN precision recall
# [1,]                10                      8                              0.8  7  1 89  3     0.875    0.7
#      specificity accuracy  F1_score gic ebic bic aic mbic
# [1,]   0.9888889     0.96 0.7777778 NaN  NaN NaN NaN  NaN
# TODO check cause of NaNs in IC in metrics



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
family <- gaussian(identity)
m <- ncol(y) # nr of tasks/channels

# using lambda expected to optimise given IC
lambda_preset = preset.lambda(IC = crit, # options aic, bic, gic, ebic, mbic 
                                  y = y, 
                                  X = X,  
                                  family = family) 
lambda_preset # lambda that would be expected to optimise given IC

system.time(L0glm.out_preset <- L0glm.fit(X, y, 
                               family = family, 
                               lambda = lambda_preset, 
                               nonzero.centered = FALSE, # use nonzero centered ridge regression?
                               grouped = TRUE,
                               relaxed = TRUE,
                               nonnegative = TRUE,
                               warmupiters = 1,                               
                               start = NULL, # rep(1, ncol(sim$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                               solver = "eigen", # note: solver osqp better than eigen for very small lambda values))
                               method = 7L,
                               thresh = 1E-5
)) # 9.79s for n=p=100, 33s for n=p=10000
L0glm.out_preset_coefs = L0glm.out_preset$coefficients
abline(v=which(rowMin(L0glm.out_preset_coefs)!=0), col = 'cyan')
graph2pdf(file="..//benchmarks/plots/2_2 multichannel deconvolution Gaussian.pdf", width=7, height=5)

# some metrics
do.call(cbind, get.metrics(beta = L0glm.out_preset_coefs, 
                           X = X, y = y, 
                           family = family, 
                           beta_true = sim$beta_true))

# coefficients & metrics over a range of lambda values
# calculate max lambda where all coefs would be penalized to 0
maxlam <- maxlambda(X=X, weights=L0glm.out_preset$weights, 
                    z=y, no.pen=NULL, coefs=L0glm.out_preset$coefficients,
                    family=family)
maxlam # 

lambdas <- 10^seq(round(log10(lambda_preset))-1, log10(maxlam), by=0.2) #
coefmat <- Matrix(0, nrow=length(lambdas), ncol=ncol(sim$X)) # matrix with avg coefficients over tasks/channels
system.time(metrics <- data.frame(do.call(rbind, lapply(1:length(lambdas), function (i) { 
  lambda <- lambdas[[i]]
  L0glm.out <- L0glm.fit(X, y, 
                             family = gaussian(identity), 
                             lambda = lambda, 
                             nonzero.centered = FALSE, # use nonzero centered ridge regression?
                             grouped = TRUE,
                             relaxed = TRUE, 
                             nonnegative = TRUE,
                             warmupiters = 1,                               
                             start = NULL, # rep(10, ncol(sim$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                             solver = "eigen", # note: solver osqp better than eigen for very small lambda values))
                             method = 7L,
                             parallellize = TRUE, 
                             thresh = 1E-5
  )$coefficients
  coefmat[i,] <<- rowMeans(L0glm.out) # avg coefficients across tasks/channels
  return( unlist(c(lambda=lambda, get.metrics(beta = L0glm.out, 
                                              X = X, y = y, 
                                              family = family, 
                                              beta_true = sim$beta_true))) )
} )))) # 300s

matplot(x=log10(metrics$lambda), y=as.matrix(coefmat), type="l", lty=1,
        xlab="Log10(lambda)", ylab="Coefficient estimates", main="Regularisation path Group L0 Gaussian GLM")

plot(log10(metrics$lambda), metrics$deviance, type="l")
metrics$lambda[which.min(metrics$deviance)] # 1.584893 # NOTE this is deviance relative to observed signal, not deviance relative to ground truth though
plot(log10(metrics$lambda), metrics$deviance_true, type="l")
metrics$lambda[which(metrics$deviance_true==min(metrics$deviance_true))] # 10.00000 15.84893 # deviance relative to ground truth
plot(log10(metrics$lambda), metrics$loglik, type="l")
metrics$lambda[which.max(metrics$loglik)] # 1.584893
plot(log10(metrics$lambda), metrics$MD_coefs, type="l")
plot(log10(metrics$lambda), metrics$MAD_coefs, type="l")
metrics$lambda[which.min(metrics$MAD_coefs)] # 10
plot(log10(metrics$lambda), metrics$sMAPE_coefs, type="l") # good one
metrics$lambda[which.min(metrics$sMAPE_coefs)] # 10
plot(log10(metrics$lambda), metrics$bias_coefs, type="l")
metrics$lambda[which.min(abs(metrics$bias_coefs-1))] # 2511.886
plot(log10(metrics$lambda), metrics$R2_coefs, type="l")
plot(log10(metrics$lambda), metrics$estimated_over_true_support_size, type="l") # good one
range(metrics$lambda[which(abs(metrics$estimated_over_true_support_size-1)==min(abs(metrics$estimated_over_true_support_size-1)))]) # 10.00000 15.84893
plot(log10(metrics$lambda), metrics$TP, type="l")
plot(log10(metrics$lambda), metrics$FP, type="l")
plot(log10(metrics$lambda), metrics$TN, type="l")
plot(log10(metrics$lambda), metrics$FN, type="l")
plot(log10(metrics$lambda), metrics$precision, type="l")
plot(log10(metrics$lambda), metrics$recall, type="l")
plot(log10(metrics$lambda), metrics$specificity, type="l")
plot(log10(metrics$lambda), metrics$accuracy, type="l") # good one
metrics$lambda[which.max(metrics$accuracy)] # 10
plot(log10(metrics$lambda), metrics$F1_score, type="l")
metrics$lambda[which.max(metrics$F1_score)] # 10
plot(log10(metrics$lambda), metrics$aic, type="l")
metrics$lambda[which(metrics$aic==min(metrics$aic))] # 10.00000 15.84893
plot(log10(metrics$lambda), metrics$bic, type="l")
metrics$lambda[which.min(metrics$bic)] # 10
plot(log10(metrics$lambda), metrics$ebic, type="l")
metrics$lambda[which.min(metrics$ebic)] # 10
plot(log10(metrics$lambda), metrics$gic, type="l")
metrics$lambda[which.min(metrics$gic)] # 10
abline(v=log10(metrics$lambda[which.min(metrics$gic)]), col="blue") 
plot(log10(metrics$lambda), metrics$mbic, type="l")
metrics$lambda[which(metrics$mbic==min(metrics$mbic))] # 63.09573 100.00000 158.48932 251.18864 398.10717 630.95734

# these lambda values are roughly in the expected to optimise these IC
preset.lambda(IC = "aic", # options aic, bic, gic, ebic, mbic 
              y = y, 
              X = X, 
              family = family) # 7.106486
preset.lambda(IC = "bic", # options aic, bic, gic, ebic, mbic 
              y = y, 
              X = X, 
              family = family) # 16.12445
preset.lambda(IC = "gic", # options aic, bic, gic, ebic, mbic 
              y = y, 
              X = X, 
              family = family) # 25.48548
preset.lambda(IC = "mbic", # options aic, bic, gic, ebic, mbic 
              y = y, 
              X = X, 
              family = family) # 40.01668


### 2.2.1 Comparison with abess relaxed group L0 mgaussian best subset glm #### 

# This was our L0 glm metrics:
metrics_L0glm = do.call(cbind, get.metrics(beta = L0glm.out_preset_coefs, 
                             X = X, y = y, 
                             family = family, 
                             beta_true = sim$beta_true))
metrics_L0glm
#        p   n   m    deviance deviance_true loglik     MD_coefs MAD_coefs sMAPE_coefs bias_coefs  R2_coefs
# [1,] 100 100 124 40000006978      6.28e+12   -Inf -0.005259792 0.1925821  0.04782878   0.995435 0.9960168
#      support_size_true support_size_estimated estimated_over_true_support_size TP FP TN FN precision recall
# [1,]                10                      9                              0.9  9  0 90  1         1    0.9
#      specificity accuracy  F1_score gic ebic bic aic mbic
# [1,]           1     0.99 0.9473684 Inf  Inf Inf Inf  Inf
# TODO fix Inf values in metrics

# abess mgaussian best subset fit plus metrics:
# best subset multitask learning using family="mgaussian" using abess
library(abess)
system.time(abess_fit <- abess(x = X, 
                   y = y, 
                   # weights = 1/(Y+0.1),  # NOTE: can be used as approx 1/variance weights instead of using family="poisson"
                   family = "mgaussian",
                   tune.path = "sequence", 
                   support.size = c(1:(p-1)),
                   lambda=0, 
                   warm.start = TRUE,
                   tune.type = "gic")) # 0.31s, or cv or ebic or bic or aic
# note: currently no support for nonnegativity constraints
plot(abess_fit, type="tune")  # optimal support size 
beta_abess = as.matrix(extract(abess_fit)$beta) # coefficient matrix for best subset
library(qlcMatrix)
image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=beta_abess^0.1, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="abess multitask mgaussian (red=true peaks)")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim$X))[as.vector(rowMin(beta_abess)!=0)], col="cyan")
sum(rowMax(sim$beta_true)!=0) # 10 true peaks
sum(rowMin(beta_abess)!=0) # 8 peaks detected
# note: for larger p, e.g. 500, solution quality quickly becomes a lot worse
# see https://github.com/abess-team/abess/issues/510

metrics_abess = do.call(cbind, get.metrics(beta = beta_abess, 
                                           X = X, y = y, 
                                           family = family, 
                                           beta_true = sim$beta_true))
metrics_abess
#        p   n   m deviance deviance_true loglik    MD_coefs MAD_coefs sMAPE_coefs bias_coefs  R2_coefs
# [1,] 100 100 124 2.05e+12     1.034e+13    NaN -0.04160188 0.4694801   0.1366435  0.9838769 0.9650847
#      support_size_true support_size_estimated estimated_over_true_support_size TP FP TN FN precision recall
# [1,]                10                      8                              0.8  7  1 89  3     0.875    0.7
#      specificity accuracy  F1_score gic ebic bic aic mbic
# [1,]   0.9888889     0.96 0.7777778 NaN  NaN NaN NaN  NaN
# TODO check cause of NaNs in IC in metrics


### 2.2.2 Comparison with glmnet nonnegative mgaussian multivariate group LASSO fit #### 

library(glmnet)
system.time(cvfit <- cv.glmnet(X, 
                               y, 
                               alpha = 1,  # lasso
                               family = "mgaussian", 
                               nlambda = 100,
                               nfolds = 5, 
                               standardize = FALSE,
                               standardize.response = FALSE,
                               intercept = FALSE, 
                               relax = FALSE,
                               lower.limits = rep(0, ncol(X)+1))) # for nonnegativity constraints
# 6.32s

plot(cvfit)

# Fit the nonnegative group LASSO glmnet model over a sequence of lambda values
system.time(fit_glmnet <- glmnet(X, 
                                 y,
                                 alpha = 1, # lasso
                                 family = "mgaussian", 
                                 nlambda = 100,
                                 standardize = FALSE,
                                 standardize.response = FALSE,
                                 intercept = FALSE, 
                                 relax = FALSE, 
                                 lower.limits = rep(0, ncol(X)+1))) # for nonnegativity constraints
# 3.11s

# Select the best lambda value (I am using cvfit$lambda.1se to get slightly sparser solution)
best_lambda <- cvfit$lambda.1se

plot(fit_glmnet, xvar="lambda")

# get the coefficients for each task for the best lambda value
coefs <- coef(fit_glmnet, s = best_lambda)
beta_mgaussian_glmnet <- do.call(cbind, lapply(seq_len(m), function (channel) as.matrix(coefs[[channel]][-1,,drop=F])))
colnames(beta_mgaussian_glmnet) = 1:ncol(beta_mgaussian_glmnet)
beta_mgaussian_glmnet[abs(beta_mgaussian_glmnet)<1E-5] = 0 # slight amount of thresholding of small coefficients using same thresh as used in L0glm

image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=beta_mgaussian_glmnet^0.1, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="nonnegative group Lasso glmnet (red=true support)")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim$X))[as.vector(rowMin(beta_mgaussian_glmnet)!=0)], col="cyan")
sum(rowMax(sim$beta_true)!=0) # 10 true peaks
sum(rowMin(beta_mgaussian_glmnet)!=0) # 13 peaks detected

metrics_glmnet = do.call(cbind, get.metrics(beta = beta_mgaussian_glmnet, 
                                           X = X, y = y, 
                                           family = family, 
                                           beta_true = sim$beta_true))
metrics_glmnet
#        p   n   m deviance deviance_true    loglik    MD_coefs MAD_coefs sMAPE_coefs bias_coefs  R2_coefs
# [1,] 100 100 124 8403.412      4.28e+12 -15396.65 -0.07910432 0.4553128  0.09181394   1.013359 0.9695373
#      support_size_true support_size_estimated estimated_over_true_support_size TP FP TN FN precision recall
# [1,]                10                     13                              1.3 10  3 87  0 0.7692308      1
#      specificity accuracy  F1_score      gic     ebic      bic      aic     mbic
# [1,]   0.9666667     0.97 0.8695652 38205.99 40500.99 35647.14 32901.29 42432.53

# note: for larger p, e.g. 500, solution quality quickly becomes a lot worse
# see https://github.com/abess-team/abess/issues/510