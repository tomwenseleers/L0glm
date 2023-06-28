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
dispersion <- 1

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
                               nonzero.centered = FALSE, # use nonzero centered ridge regression?
                               # ridge.center = start, # ridge.center
                               nonnegative = TRUE,
                               warmupiters = 1,
                               # warmupiters.nzc = 1,
                               # start = rep(1, ncol(sim$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                               solver = "eigen", # note: solver osqp & glmnet better than eigen for very small lambda values)
                               method = 7L, # 0.27s for n=p=10000 with eigen solver methods 7 & 8 (Cholesky), 0.61s with methods 5 & 6 (CG), 1.27s with methods 3, 4 (LSCG), 219s with methods 1 & 2 (QR), 4.8s for n=p=10000 with osqp solver, 2s for n=p=10000 with glmnet solver
                               thresh = 1E-10
)) 
L0glm.out_preset_coefs <- L0glm.out_preset$coefficients
sum(L0glm.out_preset_coefs!=0)  # 999
points(sim$x, L0glm.out_preset_coefs, pch = 16, col= 'blue')
graph2pdf(file="..//benchmarks/plots/1_1 single channel deconvolution Poisson.pdf", width=7, height=5)

mincoef <- min(abs(L0glm.out_preset_coefs)[L0glm.out_preset_coefs!=0]) # nonzero coefficient with smallest absolute value
do.call(cbind, get.metrics(beta = L0glm.out_preset_coefs,  # some metrics
                           X = X, y = y, 
                           family = poisson(identity), 
                           beta_true = sim$beta_true))

plot(x=L0glm.out_preset_coefs, y=sim$beta_true, pch=16, col="blue")


# # couple of possible coefficient initialisations to use as starting values in L0glm.fit which could be used for start ####
# # weights = L0glm.out_preset$weights
# weights = 1/(sim$y+0.1)
# Xw = Matrix(Diagonal(x=sqrt(weights)) %*% sim$X, sparse=F) # = X * sqrt(weights)
# Xwaug = Matrix(rbind(Xw, Diagonal(x=sqrt(rep(1,ncol(sim$X))))), sparse=F) # row augmented covariate matrix
# ywaug = as.vector(c(sim$y, rep(0, ncol(sim$X)))) # 
# 
# # orthogonal initialisation (assuming variables are orthogonal) ####
# beta_init_orthogonal <- as.vector(crossprod(Xwaug, ywaug) / colSums(Xwaug^2)) 
# cor(beta_init_orthogonal, as.vector(sim$beta_true)) # 0.3363417
# plot(beta_init_orthogonal, as.vector(sim$beta_true), pch=16, col="blue")
# mean(abs(beta_init_orthogonal-sim$beta_true)) # 485.0457
# # colSums(Xwaug^2) = diag(crossprod(Xwaug)))
# # PS original variables in columns of X can be very far from orthogonal, 
# # but the augmentation with a diagonal matrix with sqrt(lambdas) along the diagonal
# # makes the assumption hold better
# 
# # diagonal initialisation (assuming system is dominated by diagonal of X) ####
# # for signal deconvolution problem with n=p this gives the best initial solution
# # only works if X is square though
# beta_init_diagonal <- y*sqrt(weights) / diag(Xw)
# beta_init_diagonal[is.na(beta_init_diagonal)] = 0
# cor(beta_init_diagonal, as.vector(sim$beta_true)) # 0.4712314
# plot(beta_init_diagonal, as.vector(sim$beta_true), pch=16, col="blue")
# mean(abs(beta_init_diagonal-sim$beta_true)) # 1834.576

# coefficients & metrics over a range of lambda values
# calculate max lambda where all coefs would be penalized to 0
weights = L0glm.out_preset$weights
Xw <- Diagonal(x=sqrt(weights)) %*% X[,,drop=F]
sel_task <- which.max(colMeans(abs(sim$y))) # for multitask select column with max mean abs value?
yw <- y[, sel_task] * sqrt(weights) 
# TODO note that yw should really be zw, i.e. for non-identity link you should be 
# working on the adjusted response scale - maybe you should also make L0glm.fit return z? 
# Or make a convenience function that returns z and w given a particular given family and 
# coefficient estimates?
maxlambda <- max((1/colSums(Xw^2)) * (crossprod(Xw, yw)^2), na.rm=T)
# maxlambda <- maxlambda/4
maxlambda
  
lambdas <- 10^seq(round(log10(lambda_preset))-3, log10(maxlambda), by=0.05)
coefmat <- Matrix(0, nrow=length(lambdas), ncol=ncol(sim$X))
warmstart <- 0.001 # % by which L0glm.out_preset will be used as start, lowest value that led to stable coefficient paths
system.time(metrics <- data.frame(do.call(rbind, lapply(1:length(lambdas), function (i) { 
  lambda <- lambdas[[i]]
  start <- L0glm.out_preset_coefs*warmstart+rep(1,ncol(sim$X))*(1-warmstart)
  start[L0glm.out_preset_coefs==0] <- 0 # slightly faster & stable with this
  L0glm.out <- L0glm.fit(X, y, 
                             family = poisson(identity), 
                             lambda = lambda, 
                             nonzero.centered = FALSE, # use nonzero centered ridge regression?
                             # ridge.center = L0glm.out_preset_coefs*0.001, 
                             nonnegative = TRUE,
                             warmupiters = 1,
                             warmupiters.nzc = 1,
                             start = start, # rep(1, ncol(X)), # start, # rep(1, ncol(X)),  # start, # 
                             solver = "eigen", # eigen, glmnet or osqp, # note: constrained solvers glmnet & osqp better than eigen for very small lambda values
                             method = 7L, 
                             # gamma = 2, 
                             thresh = 1E-5,
                             relaxed = FALSE
                             )$coefficients 
  coefmat[i,] <<- as.vector(L0glm.out)
  return( unlist(c(lambda=lambda, get.metrics(beta = L0glm.out, 
                                          X = X, y = y, 
                                          family = poisson(identity), 
                                          beta_true = sim$beta_true))) )
  } )))) # 29s for n=p=10000 with eigen solver methods 7 & 8 (Cholesky), 165s for n=p=10000 with glmnet solver, 455s for n=p=10000 with osqp solver
matplot(x=log10(metrics$lambda), y=as.matrix(coefmat), type="l", lty=1,
        xlab="Log10(lambda)", ylab="Coefficient estimates", main="Regularisation path L0 Poisson GLM")
abline(v=log10(maxlambda), col="red") # max lambda where all coefs would be penalized to 0
library(scales)
abline(h=as.vector(sim$beta_true), col=alpha(1:6,0.2))

plot(log10(metrics$lambda), metrics$deviance, type="l")
metrics$lambda[which.min(metrics$deviance)] # # NOTE this is deviance relative to observed signal, not deviance relative to ground truth though
plot(log10(metrics$lambda), metrics$loglik, type="l")
metrics$lambda[which.max(metrics$loglik)] # 0
plot(log10(metrics$lambda), metrics$MD_coefs, type="l")
plot(log10(metrics$lambda), metrics$MAD_coefs, type="l")
metrics$lambda[which.min(metrics$MAD_coefs)] # 
plot(log10(metrics$lambda), metrics$sMAPE_coefs, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(metrics$sMAPE_coefs)] # 
plot(log10(metrics$lambda), metrics$bias_coefs, type="l")
metrics$lambda[which.min(abs(metrics$bias_coefs-1))] # 
plot(log10(metrics$lambda), metrics$R2_coefs, type="l")
plot(log10(metrics$lambda), metrics$estimated_over_true_support_size, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(abs(metrics$estimated_over_true_support_size-1))] #
plot(log10(metrics$lambda), metrics$TP, type="l")
plot(log10(metrics$lambda), metrics$FP, type="l")
plot(log10(metrics$lambda), metrics$TN, type="l")
plot(log10(metrics$lambda), metrics$FN, type="l")
plot(log10(metrics$lambda), metrics$precision, type="l")
plot(log10(metrics$lambda), metrics$recall, type="l")
plot(log10(metrics$lambda), metrics$specificity, type="l")
plot(log10(metrics$lambda), metrics$accuracy, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.max(metrics$accuracy)] # 
plot(log10(metrics$lambda), metrics$F1_score, type="l")
metrics$lambda[which.max(metrics$F1_score)] # 
plot(log10(metrics$lambda), metrics$aic, type="l")
metrics$lambda[which.min(metrics$aic)] #
plot(log10(metrics$lambda), metrics$bic, type="l")
metrics$lambda[which.min(metrics$bic)] # 
plot(log10(metrics$lambda), metrics$ebic, type="l")
metrics$lambda[which.min(metrics$ebic)] # 
plot(log10(metrics$lambda), metrics$gic, type="l")
metrics$lambda[which.min(metrics$gic)] # 
plot(log10(metrics$lambda), metrics$mbic, type="l")
metrics$lambda[which.min(metrics$mbic)] #

# lambda values that are expected to optimise each IC (with start = rep(1, ncol(X)))
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

# using lambda expected to optimise given IC
lambda_preset = preset.lambda(IC = crit, # options aic, bic, gic, ebic, mbic 
                                  y = y, 
                                  X = X, 
                                  family = gaussian(identity)) 
lambda_preset # lambda that would be expected to optimise given IC

L0glm.out_preset <- L0glm.fit(X, y, 
                               family = gaussian(identity), 
                               lambda = lambda_preset, 
                               nonzero.centered = FALSE, # use nonzero centered ridge regression?
                               nonnegative = TRUE,
                               warmupiters = 1,
                               # start = rep(10, ncol(sim$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                               solver = "eigen" # note: solver osqp better than eigen for very small lambda values))
)
L0glm.out_preset_coefs = L0glm.out_preset$coefficients
points(sim$x, L0glm.out_preset_coefs, pch = 16, col= 'blue')

# some metrics
do.call(cbind, get.metrics(beta = L0glm.out_preset_coefs, 
                           X = X, y = y, 
                           family = gaussian(identity), 
                           beta_true = sim$beta_true))
graph2pdf(file="..//benchmarks/plots/1_2 single channel deconvolution Gaussian.pdf", width=7, height=5)



# coefficients & metrics over a range of lambda values
# calculate max lambda where all coefs would be penalized to 0
weights <- L0glm.out_preset$weights
Xw <- Diagonal(x=sqrt(weights)) %*% X
sel_task <- which.max(colMeans(abs(sim$y))) # for multitask select column with max mean abs value?
yw <- y[, sel_task] * sqrt(weights) 
# TODO note that yw should really be zw, i.e. for non-identity link you should be 
# working on the adjusted response scale - maybe you should also make L0glm.fit return z? 
# Or make a convenience function that returns z and w given a particular given family and 
# coefficient estimates?
maxlambda <- max((1/colSums(Xw^2)) * (crossprod(Xw, yw)^2), na.rm=T)
maxlambda <- maxlambda*(1/4) # factor of 1/4 of Frommlet & Nuel to convert lambda to lambda tilde for Gaussian case

lambdas <- 10^seq(round(log10(lambda_preset))-3, log10(maxlambda), by=0.05)
coefmat <- Matrix(0, nrow=length(lambdas), ncol=ncol(sim$X))
warmstart <- 0.2 # % by which L0glm.out_preset will be used as start, lowest value that led to stable coefficient paths
metrics <- data.frame(do.call(rbind, lapply(1:length(lambdas), function (i) { 
  lambda <- lambdas[[i]]
  start <- L0glm.out_preset_coefs*warmstart+rep(1,ncol(sim$X))*(1-warmstart)
  start[L0glm.out_preset_coefs==0] <- 0
  L0glm.out <- L0glm.fit(X, y, 
                             family = gaussian(identity), 
                             lambda = lambda, 
                             nonzero.centered = FALSE, # use nonzero centered ridge regression?
                             # ridge.center = start;
                             nonnegative = TRUE,
                             warmupiters = 1,
                             start = start, # rep(1, ncol(sim$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                             solver = "eigen", # note: solver osqp better than eigen for very small lambda values))
                             method = 7L,
                             relaxed = FALSE,
                             thresh = 1E-10,
)$coefficients
  coefmat[i,] <<- as.vector(L0glm.out)
  return( unlist(c(lambda=lambda, get.metrics(beta = L0glm.out, 
                                              X = X, y = y, 
                                              family = gaussian(identity), 
                                              beta_true = sim$beta_true))) )
} )))
matplot(x=log10(metrics$lambda), y=as.matrix(coefmat), type="l", lty=1,
        xlab="Log10(lambda)", ylab="Coefficient estimates", main="Regularisation path L0 Gaussian GLM")
abline(v=log10(maxlambda), col="red") # max lambda where all coefs would be penalized to 0
library(scales)
abline(h=as.vector(sim$beta_true), col=alpha(1:6,0.2))

plot(log10(metrics$lambda), metrics$deviance, type="l")
metrics$lambda[which.min(metrics$deviance)] # # NOTE this is deviance relative to observed signal, not deviance relative to ground truth though
plot(log10(metrics$lambda), metrics$loglik, type="l")
metrics$lambda[which.max(metrics$loglik)] # 0
# TODO: check why loglikelihood comes out so high for the very high lambda values - does not look right
plot(log10(metrics$lambda), metrics$MD_coefs, type="l")
plot(log10(metrics$lambda), metrics$MAD_coefs, type="l")
metrics$lambda[which.min(metrics$MAD_coefs)] # 
plot(log10(metrics$lambda), metrics$sMAPE_coefs, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(metrics$sMAPE_coefs)] # 
plot(log10(metrics$lambda), metrics$bias_coefs, type="l")
metrics$lambda[which.min(abs(metrics$bias_coefs-1))] # 
plot(log10(metrics$lambda), metrics$R2_coefs, type="l")
plot(log10(metrics$lambda), metrics$estimated_over_true_support_size, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(abs(metrics$estimated_over_true_support_size-1))] #
plot(log10(metrics$lambda), metrics$TP, type="l")
plot(log10(metrics$lambda), metrics$FP, type="l")
plot(log10(metrics$lambda), metrics$TN, type="l")
plot(log10(metrics$lambda), metrics$FN, type="l")
plot(log10(metrics$lambda), metrics$precision, type="l")
plot(log10(metrics$lambda), metrics$recall, type="l")
plot(log10(metrics$lambda), metrics$specificity, type="l")
plot(log10(metrics$lambda), metrics$accuracy, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.max(metrics$accuracy)] # 
plot(log10(metrics$lambda), metrics$F1_score, type="l")
metrics$lambda[which.max(metrics$F1_score)] # 
plot(log10(metrics$lambda), metrics$aic, type="l")
metrics$lambda[which.min(metrics$aic)] #
plot(log10(metrics$lambda), metrics$bic, type="l")
metrics$lambda[which.min(metrics$bic)] # 
plot(log10(metrics$lambda), metrics$ebic, type="l")
metrics$lambda[which.min(metrics$ebic)] # 
plot(log10(metrics$lambda), metrics$gic, type="l")
metrics$lambda[which.min(metrics$gic)] # 
plot(log10(metrics$lambda), metrics$mbic, type="l")
metrics$lambda[which.min(metrics$mbic)] #

# lambda values that a priori are expected to optimise each IC (with start = rep(1, ncol(X)))
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

system.time(L0glm.out_preset <- L0glm.fit(X, y, 
                               family = poisson(identity), 
                               lambda = lambda_preset, 
                               nonzero.centered = FALSE, # use nonzero centered ridge regression?
                               # ridge.center = XX,
                               grouped = TRUE,
                               relaxed = TRUE, # NOTE: CHECK setting relaxed to TRUE does not seem to improve solution quality
                               nonnegative = TRUE,
                               warmupiters = 1,                               
                             # start = rep(10, ncol(sim$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                               solver = "eigen", # note: solver osqp better than eigen for very small lambda values))
                               method = 7L,
                               thresh = 1E-10
)) # 14.97s
L0glm.out_preset_coefs = L0glm.out_preset$coefficients
abline(v=which(rowMin(L0glm.out_preset_coefs)!=0), col = 'cyan')
graph2pdf(file="..//benchmarks/plots/2_1 multichannel deconvolution Poisson.pdf", width=7, height=5)

# some metrics
do.call(cbind, get.metrics(beta = L0glm.out_preset_coefs, 
                           X = X, y = y, 
                           family = poisson(identity), 
                           beta_true = sim$beta_true))

# coefficients & metrics over a range of lambda values
lambdas <- 10^seq(-3, 7, by=1)
coefmat <- Matrix(0, nrow=length(lambdas), ncol=ncol(sim$X)) # matrix with avg coefficients over tasks/channels
metrics <- data.frame(do.call(rbind, lapply(1:length(lambdas), function (i) { 
  print(i)
  lambda <- lambdas[[i]] # for i=8 gives error at L0glm.R#742 number of rows of result is not a multiple of vector length (arg 8)
  L0glm.out <- L0glm.fit(X, y, 
                             family = poisson(identity), 
                             lambda = lambda, 
                             nonzero.centered = FALSE, # use nonzero centered ridge regression?
                             grouped = TRUE,
                             relaxed = TRUE, 
                             nonnegative = TRUE,
                             warmupiters = 1,                               
                             # start = rep(1, ncol(sim$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                             solver = "eigen", # note: solver osqp better than eigen for very small lambda values))
                             method = 7L,
                             parallelize = FALSE, # TODO produces error with parallellize=TRUE
                             thresh = 1E-10
  )$coefficients
  coefmat[i,] <<- rowMeans(L0glm.out) # avg coefficients across tasks/channels
  return( unlist(c(lambda=lambda, get.metrics(beta = L0glm.out, 
                                              X = X, y = y, 
                                              family = poisson(identity), 
                                              beta_true = sim$beta_true) )) )
} ) ))
# TODO fix bug with parallelize = TRUE returns Error in replCmat4(x, i1 = if (iMi) seq.int(from = 0L, length.out = di[1L]) else .ind.prep2(i,  : 
# number of items to replace is not a multiple of replacement length at L0glm.R#742
# In addition: Warning message:
# In (function (..., deparse.level = 1)  :
# number of rows of result is not a multiple of vector length (arg 8)
  
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

L0glm.out_preset <- L0glm.fit(X, y, 
                               family = gaussian(identity), 
                               lambda = lambda_preset, 
                               nonzero.centered = FALSE, # use nonzero centered ridge regression?
                               grouped = TRUE,
                               relaxed = TRUE,
                               nonnegative = TRUE,
                               warmupiters = 1,                               
                             # start = rep(1, ncol(sim$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                               solver = "eigen", # note: solver osqp better than eigen for very small lambda values))
                               method = 7L,
                               thresh = 1E-10
)
L0glm.out_preset_coefs = L0glm.out_preset$coefficients
abline(v=which(rowMin(L0glm.out_preset_coefs)!=0), col = 'cyan')
graph2pdf(file="..//benchmarks/plots/2_2 multichannel deconvolution Gaussian.pdf", width=7, height=5)

# some metrics
do.call(cbind, get.metrics(beta = L0glm.out_preset_coefs, 
                           X = X, y = y, 
                           family = gaussian(identity), 
                           beta_true = sim$beta_true))

# coefficients & metrics over a range of lambda values
lambdas <- 10^seq(-3, 7, by=1)
coefmat <- Matrix(0, nrow=length(lambdas), ncol=ncol(sim$X)) # matrix with avg coefficients over tasks/channels
metrics <- data.frame(do.call(rbind, lapply(1:length(lambdas), function (i) { 
  print(i)
  # TODO fix bug for i=
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
                             solver = "eigen", # note: solver osqp better than eigen for very small lambda values))
                             method = 7L,
                             parallellize = FALSE, # TODO produces error with parallellize=TRUE
                             thresh = 1E-10
  )$coefficients
  coefmat[i,] <<- rowMeans(L0glm.out) # avg coefficients across tasks/channels
  return( unlist(c(lambda=lambda, get.metrics(beta = L0glm.out, 
                                              X = X, y = y, 
                                              family = gaussian(identity), 
                                              beta_true = sim$beta_true))) )
} )))
# TODO fix bug : for i=8 and with parallellize = TRUE produces error at L0glm.R#742
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
