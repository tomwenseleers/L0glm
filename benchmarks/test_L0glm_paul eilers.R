# Minimal tests of L0glm package ####

library(remotes)
# remotes::install_github("tomwenseleers/L0glm/L0glm")
library(L0glm)
library(export)
# devtools::load_all(".")

n <- p <- 100 # to check performance use e.g. n=p=10000, otherwise use e.g. n=p=100-300
s <- 0.1 # sparsity, nr of timepoints with a peak 
k <- round(n*s)
mean_beta <- 2000
sd_logbeta <- 1
sd_noise <- 5 # SD of noise for Gaussian noise case
crit <- "gic" # choose lambda based on what would be expected to optimise this IC
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
system.time(L0glm.out <- L0glm.fit(X, y, 
                           family = poisson(identity), 
                           lambda = 1, 
                           nonzero.centered = TRUE, # use nonzero centered ridge regression?
                           nonnegative = TRUE,
                           warmupiters = 4,
                           # start = rep(10, ncol(sim$X)), # coefficient starting values
                           solver = "eigen", # note: solver osqp better than eigen for very small lambda values)
                           method = 7L # 1 & 2 = QR, 3 & 4 = LSCG, 5 & 6 = CG, 7 & 8 = Cholesky
)$coefficients) # 0.29s for n=p=10000 with eigen solver methods 7 & 8 (Cholesky), 0.75s with methods 5 & 6 (CG), 1.9s with methods 3, 4 (LSCG), 219s with methods 1 & 2 (QR), 5s for n=p=10000 with osqp solver, 2.27s for n=p=10000 with glmnet solver
# points(sim$x, L0glm.out, pch = 16, col= 'blue')

# using lambda expected to optimise given IC
lambda_preset = preset.lambda(IC = crit, # options aic, bic, gic, ebic, mbic 
                                  y = y, 
                                  X = X, 
                                  family = poisson(identity)) 
lambda_preset # lambda = 1.15 would be expected to optimise BIC with n=p=100, lambda = 1.76 would optimise GIC

system.time(L0glm.out_preset <- L0glm.fit(X, y, 
                               family = poisson(identity), 
                               lambda = lambda_preset, 
                               nonzero.centered = TRUE, # use nonzero centered ridge regression?
                               nonnegative = TRUE,
                               warmupiters = 4,
                               # start = rep(10, ncol(sim$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                               solver = "eigen", # note: solver osqp & glmnet better than eigen for very small lambda values)
                               method = 7L # 0.27s for n=p=10000 with eigen solver methods 7 & 8 (Cholesky), 0.61s with methods 5 & 6 (CG), 1.27s with methods 3, 4 (LSCG), 219s with methods 1 & 2 (QR), 4.8s for n=p=10000 with osqp solver, 2s for n=p=10000 with glmnet solver
)) 
L0glm.out_preset_coefs <- L0glm.out_preset$coefficients
points(sim$x, L0glm.out_preset_coefs, pch = 16, col= 'blue')

mincoef <- min(abs(L0glm.out_preset_coefs)[L0glm.out_preset_coefs!=0]) # nonzero coefficient with smallest absolute value
do.call(cbind, get.metrics(beta = L0glm.out_preset_coefs,  # some metrics
                           X = X, y = y, 
                           family = poisson(identity), 
                           beta_true = sim$beta_true))
graph2pdf(file="..//benchmarks/plots/1_1 single channel deconvolution Poisson.pdf", width=7, height=5)

plot(x=L0glm.out_preset_coefs, y=sim$beta_true, pch=16, col="blue")


# coefficients & metrics over a range of lambda values
# calculate max lambda where all coefs would be penalized to 0
weights <- L0glm.out_preset$weights
Xw <- Diagonal(x=sqrt(weights)) %*% X
sel_task <- which.max(colMeans(abs(sim$y))) # for multitask select column with max mean abs value?
yw <- y[, sel_task] * sqrt(weights) # TODO note that yw should really be zw, i.e. for non-identity link you should be working on the adjusted response scale - maybe you should also make L0glm.fit return z? Or make a convenience function that returns z and w given a particular given family and coefficient estimates?
maxlambda <- max((1/colSums(Xw^2)) * (crossprod(Xw, yw)^2), na.rm=T)
maxlambda <- maxlambda*4
  
lambdas <- 10^seq(round(log10(lambda_preset))-3, log10(maxlambda)+1, by=0.1)
coefmat <- Matrix(0, nrow=length(lambdas), ncol=ncol(sim$X))
warmstart <- 0.1 # % by which L0glm.out_preset will be used as start, lowest value that led to stable coefficient paths
system.time(metrics <- data.frame(do.call(rbind, lapply(1:length(lambdas), function (i) { 
  lambda <- lambdas[[i]]
  start <- L0glm.out_preset_coefs*warmstart+rep(mincoef/2,ncol(sim$X))*(1-warmstart)
  start[L0glm.out_preset_coefs==0] <- 0
  L0glm.out <- L0glm.fit(X, y, 
                             family = poisson(identity), 
                             lambda = lambda, 
                             nonzero.centered = FALSE, # use nonzero centered ridge regression?
                             nonnegative = TRUE,
                             warmupiters = 4,
                             start = rep(1, ncol(X)), # start,
                             solver = "eigen", # eigen, glmnet or osqp, # note: constrained solvers glmnet & osqp better than eigen for very small lambda values
                             method = 7L, 
                             # gamma = 2, 
                             # thresh = 10
                             relaxed = FALSE
                             )$coefficients 
  coefmat[i,] <<- as.vector(L0glm.out)
  return( unlist(c(lambda=lambda, get.metrics(beta = L0glm.out, 
                                          X = X, y = y, 
                                          family = poisson(identity), 
                                          beta_true = sim$beta_true))) )
  } )))) # 32s for n=p=10000 with eigen solver methods 7 & 8 (Cholesky), 165s for n=p=10000 with glmnet solver, 455s for n=p=10000 with osqp solver
matplot(x=log10(metrics$lambda), y=as.matrix(coefmat), type="l", lty=1,
        xlab="Log10(lambda)", ylab="Coefficient estimates", main="Regularisation path L0 Poisson GLM")
abline(v=log10(maxlambda), col="red") # max lambda where all coefs would be penalized to 0
# TODO: looks like a bug - before we had all coefficients penalized to zero at this value, now we don't

plot(log10(metrics$lambda), metrics$deviance, type="l")
plot(log10(metrics$lambda), metrics$loglik, type="l")
plot(log10(metrics$lambda), metrics$MD_coefs, type="l")
plot(log10(metrics$lambda), metrics$MAD_coefs, type="l")
plot(log10(metrics$lambda), metrics$sMAPE_coefs, type="l") # good one, optimal for lambda between 1 & 10
plot(log10(metrics$lambda), metrics$bias_coefs, type="l")
plot(log10(metrics$lambda), metrics$R2_coefs, type="l")
plot(log10(metrics$lambda), metrics$estimated_over_true_support_size, type="l") # good one, optimal for lambda between 1 & 10
plot(log10(metrics$lambda), metrics$TP, type="l")
plot(log10(metrics$lambda), metrics$FP, type="l")
plot(log10(metrics$lambda), metrics$TN, type="l")
plot(log10(metrics$lambda), metrics$FN, type="l")
plot(log10(metrics$lambda), metrics$precision, type="l")
plot(log10(metrics$lambda), metrics$recall, type="l")
plot(log10(metrics$lambda), metrics$specificity, type="l")
plot(log10(metrics$lambda), metrics$accuracy, type="l") # good one, optimal for lambda between 1 & 10
plot(log10(metrics$lambda), metrics$F1_score, type="l")
plot(log10(metrics$lambda), metrics$gic, type="l")
plot(log10(metrics$lambda), metrics$ebic, type="l")
plot(log10(metrics$lambda), metrics$bic, type="l")
plot(log10(metrics$lambda), metrics$aic, type="l")
plot(log10(metrics$lambda), metrics$mbic, type="l")


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
# L0glm.out <- L0glm.fit(X, y, 
#                           family = gaussian(identity), 
#                           lambda = 1, 
#                           nonzero.centered = TRUE, # use nonzero centered ridge regression?
#                           nonnegative = TRUE,
#                           warmupiters = 1,
#                           start = rep(10, ncol(sim$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
#                           solver = "eigen" # note: solver osqp better than eigen for very small lambda values))
# )
# points(sim$x, L0glm.out, pch = 16, col= 'blue')

# using lambda expected to optimise given IC
lambda_preset = preset.lambda(IC = crit, # options aic, bic, gic, ebic, mbic 
                                  y = y, 
                                  X = X, 
                                  family = gaussian(identity)) 
lambda_preset # lambda that would be expected to optimise given IC

L0glm.out <- L0glm.fit(X, y, 
                               family = gaussian(identity), 
                               lambda = lambda_preset, 
                               nonzero.centered = TRUE, # use nonzero centered ridge regression?
                               nonnegative = TRUE,
                               warmupiters = 1,
                               # start = rep(10, ncol(sim$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                               solver = "eigen" # note: solver osqp better than eigen for very small lambda values))
)$coefficients
points(sim$x, L0glm.out, pch = 16, col= 'blue')

# some metrics
do.call(cbind, get.metrics(beta = L0glm.out, 
                           X = X, y = y, 
                           family = gaussian(identity), 
                           beta_true = sim$beta_true))
graph2pdf(file="..//benchmarks/plots/1_2 single channel deconvolution Gaussian.pdf", width=7, height=5)

# coefficients & metrics over a range of lambda values
lambdas <- 10^seq(-3, 7, by=0.1)
coefmat <- Matrix(0, nrow=length(lambdas), ncol=ncol(sim$X))
metrics <- data.frame(do.call(rbind, lapply(1:length(lambdas), function (i) { 
  lambda <- lambdas[[i]]
  L0glm.out <- L0glm.fit(X, y, 
                             family = gaussian(identity), 
                             lambda = lambda, 
                             nonzero.centered = TRUE, # use nonzero centered ridge regression?
                             nonnegative = TRUE,
                             warmupiters = 1,
                             # start = rep(10, ncol(sim$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                             solver = "eigen" # note: solver osqp better than eigen for very small lambda values))
)$coefficients
  coefmat[i,] <<- as.vector(L0glm.out)
  return( unlist(c(lambda=lambda, get.metrics(beta = L0glm.out, 
                                              X = X, y = y, 
                                              family = gaussian(identity), 
                                              beta_true = sim$beta_true))) )
} )))
matplot(x=log10(metrics$lambda), y=as.matrix(coefmat), type="l", lty=1,
        xlab="Log10(lambda)", ylab="Coefficient estimates", main="Regularisation path L0 Gaussian GLM")
plot(log10(metrics$lambda), metrics$deviance, type="l")
plot(log10(metrics$lambda), metrics$loglik, type="l")
plot(log10(metrics$lambda), metrics$MD_coefs, type="l")
plot(log10(metrics$lambda), metrics$MAD_coefs, type="l")
plot(log10(metrics$lambda), metrics$sMAPE_coefs, type="l") # good one, optimal for lambda around 10-100
plot(log10(metrics$lambda), metrics$bias_coefs, type="l")
plot(log10(metrics$lambda), metrics$R2_coefs, type="l")
plot(log10(metrics$lambda), metrics$estimated_over_true_support_size, type="l") # good one, optimal for lambda between 10-100
plot(log10(metrics$lambda), metrics$TP, type="l")
plot(log10(metrics$lambda), metrics$FP, type="l")
plot(log10(metrics$lambda), metrics$TN, type="l")
plot(log10(metrics$lambda), metrics$FN, type="l")
plot(log10(metrics$lambda), metrics$precision, type="l")
plot(log10(metrics$lambda), metrics$recall, type="l")
plot(log10(metrics$lambda), metrics$specificity, type="l")
plot(log10(metrics$lambda), metrics$accuracy, type="l") # good one, optimal for lambda between 10-100
plot(log10(metrics$lambda), metrics$F1_score, type="l")
plot(log10(metrics$lambda), metrics$gic, type="l") # optimum around 10
plot(log10(metrics$lambda), metrics$ebic, type="l")
plot(log10(metrics$lambda), metrics$bic, type="l")
plot(log10(metrics$lambda), metrics$aic, type="l")
plot(log10(metrics$lambda), metrics$mbic, type="l")


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

L0glm.out <- L0glm.fit(X, y, 
                               family = poisson(identity), 
                               lambda = lambda_preset, 
                               nonzero.centered = FALSE, # use nonzero centered ridge regression?
                               grouped = TRUE,
                               relaxed = TRUE, 
                               nonnegative = TRUE,
                               warmupiters = 4,                               
                             # start = rep(10, ncol(sim$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                               solver = "eigen" # note: solver osqp better than eigen for very small lambda values))
)$coefficients
abline(v=which(rowMin(L0glm.out)!=0), col = 'cyan')
graph2pdf(file="..//benchmarks/plots/2_1 multichannel deconvolution Poisson.pdf", width=7, height=5)

# some metrics
do.call(cbind, get.metrics(beta = L0glm.out, 
                           X = X, y = y, 
                           family = poisson(identity), 
                           beta_true = sim$beta_true))

# coefficients & metrics over a range of lambda values
lambdas <- 10^seq(-3, 7, by=0.5)
coefmat <- Matrix(0, nrow=length(lambdas), ncol=ncol(sim$X)) # matrix with avg coefficients over tasks/channels
metrics <- data.frame(do.call(rbind, lapply(1:length(lambdas), function (i) { 
  lambda <- lambdas[[i]]
  L0glm.out <- L0glm.fit(X, y, 
                             family = poisson(identity), 
                             lambda = lambda, 
                             nonzero.centered = FALSE, # use nonzero centered ridge regression?
                             grouped = TRUE,
                             relaxed = TRUE, 
                             nonnegative = TRUE,
                             warmupiters = 4,                               
                             # start = rep(10, ncol(sim$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                             solver = "eigen" # note: solver osqp better than eigen for very small lambda values))
  )$coefficients
  coefmat[i,] <<- rowMeans(L0glm.out) # avg coefficients across tasks/channels
  return( unlist(c(lambda=lambda, get.metrics(beta = L0glm.out, 
                                              X = X, y = y, 
                                              family = poisson(identity), 
                                              beta_true = sim$beta_true))) )
} )))
# TODO : returns error Error in replCmat4(x, i1 = if (iMi) seq.int(from = 0L, length.out = di[1L]) else .ind.prep2(i,  : 
# number of items to replace is not a multiple of replacement length

matplot(x=log10(metrics$lambda), y=as.matrix(coefmat), type="l", lty=1,
        xlab="Log10(lambda)", ylab="Coefficient estimates", main="Regularisation path Group L0 Poisson GLM")
plot(log10(metrics$lambda), metrics$deviance, type="l")
plot(log10(metrics$lambda), metrics$loglik, type="l")
plot(log10(metrics$lambda), metrics$MD_coefs, type="l")
plot(log10(metrics$lambda), metrics$MAD_coefs, type="l") # good one, optimal for lambda of 10
plot(log10(metrics$lambda), metrics$sMAPE_coefs, type="l") # good one, optimal for lambda of 10
plot(log10(metrics$lambda), metrics$bias_coefs, type="l")
plot(log10(metrics$lambda), metrics$R2_coefs, type="l")
plot(log10(metrics$lambda), metrics$estimated_over_true_support_size, type="l") # good one, optimal for lambda of 1-10
plot(log10(metrics$lambda), metrics$TP, type="l")
plot(log10(metrics$lambda), metrics$FP, type="l")
plot(log10(metrics$lambda), metrics$TN, type="l")
plot(log10(metrics$lambda), metrics$FN, type="l")
plot(log10(metrics$lambda), metrics$precision, type="l")
plot(log10(metrics$lambda), metrics$recall, type="l")
plot(log10(metrics$lambda), metrics$specificity, type="l")
plot(log10(metrics$lambda), metrics$accuracy, type="l") # good one, optimal for lambda of 10
plot(log10(metrics$lambda), metrics$F1_score, type="l")
plot(log10(metrics$lambda), metrics$gic, type="l")
plot(log10(metrics$lambda), metrics$ebic, type="l")
plot(log10(metrics$lambda), metrics$bic, type="l")
plot(log10(metrics$lambda), metrics$aic, type="l")
plot(log10(metrics$lambda), metrics$mbic, type="l")




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

L0glm.out <- L0glm.fit(X, y, 
                               family = gaussian(identity), 
                               lambda = lambda_preset, 
                               nonzero.centered = FALSE, # use nonzero centered ridge regression?
                               grouped = TRUE,
                               relaxed = TRUE,
                               nonnegative = TRUE,
                               warmupiters = 1,                               
                             # start = rep(10, ncol(sim$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                               solver = "eigen" # note: solver osqp better than eigen for very small lambda values))
)$coefficients
abline(v=which(rowMin(L0glm.out)!=0), col = 'cyan')
graph2pdf(file="..//benchmarks/plots/2_2 multichannel deconvolution Gaussian.pdf", width=7, height=5)

# some metrics
do.call(cbind, get.metrics(beta = L0glm.out, 
                           X = X, y = y, 
                           family = gaussian(identity), 
                           beta_true = sim$beta_true))

# coefficients & metrics over a range of lambda values
lambdas <- 10^seq(-3, 7, by=0.5)
coefmat <- Matrix(0, nrow=length(lambdas), ncol=ncol(sim$X)) # matrix with avg coefficients over tasks/channels
metrics <- data.frame(do.call(rbind, lapply(1:length(lambdas), function (i) { 
  lambda <- lambdas[[i]]
  L0glm.out <- L0glm.fit(X, y, 
                             family = gaussian(identity), 
                             lambda = lambda, 
                             nonzero.centered = FALSE, # use nonzero centered ridge regression?
                             grouped = TRUE,
                             relaxed = TRUE, 
                             nonnegative = TRUE,
                             warmupiters = 1,                               
                           # start = rep(10, ncol(sim$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                             solver = "eigen" # note: solver osqp better than eigen for very small lambda values))
  )$coefficients
  coefmat[i,] <<- rowMeans(L0glm.out) # avg coefficients across tasks/channels
  return( unlist(c(lambda=lambda, get.metrics(beta = L0glm.out, 
                                              X = X, y = y, 
                                              family = gaussian(identity), 
                                              beta_true = sim$beta_true))) )
} )))
# TODO produces error
# Error in .local(x, i, j, ..., value) :
# RHS 'value' (class list) matches 'ANY', but must match matrix class dgCMatrix
matplot(x=log10(metrics$lambda), y=as.matrix(coefmat), type="l", lty=1,
        xlab="Log10(lambda)", ylab="Coefficient estimates", main="Regularisation path Group L0 Gaussian GLM")
plot(log10(metrics$lambda), metrics$deviance, type="l")
plot(log10(metrics$lambda), metrics$loglik, type="l")
plot(log10(metrics$lambda), metrics$MD_coefs, type="l")
plot(log10(metrics$lambda), metrics$MAD_coefs, type="l") # good one, optimal for lambda of ca 10-100
plot(log10(metrics$lambda), metrics$sMAPE_coefs, type="l") # good one, optimal for lambda of ca 10
plot(log10(metrics$lambda), metrics$bias_coefs, type="l")
plot(log10(metrics$lambda), metrics$R2_coefs, type="l")
plot(log10(metrics$lambda), metrics$estimated_over_true_support_size, type="l") # good one, optimal for lambda of ca 10-100
plot(log10(metrics$lambda), metrics$TP, type="l")
plot(log10(metrics$lambda), metrics$FP, type="l")
plot(log10(metrics$lambda), metrics$TN, type="l")
plot(log10(metrics$lambda), metrics$FN, type="l")
plot(log10(metrics$lambda), metrics$precision, type="l")
plot(log10(metrics$lambda), metrics$recall, type="l")
plot(log10(metrics$lambda), metrics$specificity, type="l")
plot(log10(metrics$lambda), metrics$accuracy, type="l") # good one, optimal for lambda of ca 10-100
plot(log10(metrics$lambda), metrics$F1_score, type="l") # good one, optimal for lambda of ca 10-100
plot(log10(metrics$lambda), metrics$gic, type="l")
plot(log10(metrics$lambda), metrics$ebic, type="l")
plot(log10(metrics$lambda), metrics$bic, type="l")
plot(log10(metrics$lambda), metrics$aic, type="l")
plot(log10(metrics$lambda), metrics$mbic, type="l")
