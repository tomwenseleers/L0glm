# Check of metrics (AIC, BIC, loglik) L0glm package & comparison regular GLM ####

library(remotes)
# remotes::install_github("tomwenseleers/L0glm/L0glm")
library(L0glm)
# devtools::load_all(".")
library(export)
source("..//benchmarks/generate data.R")

# check consistency of L0glm output with glm for gaussian identity link, poisson log link, Gamma log link & binomial logit link
n = 10000
p = 100
n_nonzero = 50
intercept = 1
dispersion = 1
ntot = 10 # nr of trials per row in dataset for binomial case
crit = "gic"

# 1. GAUSSIAN GLM ####

set.seed(1)
data_gaus  = simulate_glm(n = n, p = p, n_nonzero = n_nonzero, intercept = intercept, family = "gaussian") 
X = data_gaus$X
y = data_gaus$y
family = gaussian(identity)

# glm fit
system.time(glm_gaus <- glm(y ~ 0+X, family=gaussian(identity))) # 0.16s
system.time(glm_gaus2 <- glm.fit(y=y, x=X, intercept=FALSE, family=gaussian(identity))) # 0.14s
glm_gaus$deviance 
glm_gaus$aic
AIC(glm_gaus)
BIC(glm_gaus)
logLik(glm_gaus)
cor(coef(glm_gaus), data_gaus$beta_true) # 0.9986453
plot(x=coef(glm_gaus), y=data_gaus$beta_true, pch=16, col="blue")

# L0 glm fit with lambda close to zero
system.time(L0glm_gaus_lambda0 <- L0glm.fit(X = X, 
                                            y = y, 
                                            family = gaussian(identity), 
                                            lambda = 1E-30, 
                                            nonzero.centered = FALSE, # use nonzero centered ridge regression?
                                            nonnegative = FALSE,
                                            warmupiters = 1,
                                            start = rep(0.1, ncol(data_gaus$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                                            solver = "eigen", # 
                                            method = 3L, 
                                            relaxed = FALSE,
                                            thresh = 1E-10
))  # 0.06s with method=3L, 0.14s with method=7L

# NOTE: note that L0glm can fit GLMs with perfectly collinear variables, e.g.
# with X = cbind(data_gaus$X, data_gaus$X)
# in that case, it would divide the total effect over the different collinear versions
# of each variable

max(abs(L0glm_gaus_lambda0$coefficients-coef(glm_gaus))) # 8.246142e-11, correct
cor(L0glm_gaus_lambda0$coefficients, data_gaus$beta_true) # 0.9986713
plot(x=L0glm_gaus_lambda0$coefficients, y=data_gaus$beta_true, pch=16, col="blue")


# check if metrics (loglik & AIC & BIC) are correct
system.time(metrics_gaus <- data.frame(do.call(cbind, get.metrics(beta = L0glm_gaus_lambda0$coefficients,  # some metrics
                                                     X = data_gaus$X, 
                                                     y = data_gaus$y,
                                                     family = gaussian(identity), 
                                                     beta_true = as.matrix(data_gaus$beta_true))))) # 0.02s
glm_gaus$deviance-metrics_gaus$deviance # 0, correct
logLik(glm_gaus)-metrics_gaus$loglik # 0, correct
AIC(glm_gaus)-metrics_gaus$aic # 0, correct
BIC(glm_gaus)-metrics_gaus$bic # 0, correct


# L0 glm fit with lambda preset to optimise given IC
lambda_preset_gaus = preset.lambda(IC = crit, # options aic, bic, gic, ebic, mbic 
                              y = data_gaus$y,
                              X = data_gaus$X, 
                              family = gaussian(identity),
                              dispersion = dispersion # dispersion here was 1, here we pass this value to not have to estimate it 
                              ) 
lambda_preset_gaus # 2.561769

system.time(L0glm_gaus_presetlambda <- L0glm.fit(X = data_gaus$X, 
                                            y = data_gaus$y, 
                                            family = gaussian(identity), 
                                            lambda = lambda_preset_gaus,
                                            nonzero.centered = FALSE, # use nonzero centered ridge regression?
                                            # ridge.center = XX,
                                            nonnegative = FALSE,
                                            warmupiters = 1,
                                            start = rep(0.1, ncol(data_gaus$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                                            solver = "eigen", # 
                                            method = 3L, #
                                            relaxed = TRUE,
                                            thresh = 1E-10
)) # 0.61s
sum(L0glm_gaus_presetlambda$coefficients!=0) 
# with nonzero.centered=FALSE 48, so OK
# TODO: check : not OK with nonzero.centered=TRUE - then gives 100, both if 
# start values are or are not given, ie a non-sparse solution
# when before solution was sparse also with nonzero.centered=TRUE
cor(L0glm_gaus_presetlambda$coefficients, data_gaus$beta_true) # 0.9992506
plot(x=L0glm_gaus_presetlambda$coefficients, y=data_gaus$beta_true, pch=16, col="blue")

# couple of possible coefficient initialisations to use as starting values in L0glm.fit which could be used for start ####
weights = L0glm_gaus_presetlambda$weights
Xw = Matrix(Diagonal(x=sqrt(weights)) %*% data_gaus$X, sparse=F) # = X * sqrt(weights)
Xwaug = Matrix(rbind(Xw, Diagonal(x=sqrt(rep(lambda_preset_gaus,ncol(data_gaus$X))))), sparse=F) # row augmented covariate matrix
ywaug = as.vector(c(data_gaus$y, rep(0, ncol(data_gaus$X)))) # 

# orthogonal initialisation (assuming variables are orthogonal) ####
beta_init_orthogonal <- as.vector(crossprod(Xwaug, ywaug) / colSums(Xwaug^2)) 
cor(beta_init_orthogonal, data_gaus$beta_true) # 0.9862637
plot(beta_init_orthogonal, data_gaus$beta_true, pch=16, col="blue")
mean(abs(beta_init_orthogonal-data_gaus$beta_true)) # 0.3878541
# colSums(Xwaug^2) = diag(crossprod(Xwaug)))
# PS original variables in columns of X can be very far from orthogonal, 
# but the augmentation with a diagonal matrix with sqrt(lambdas) along the diagonal
# makes the assumption hold better

# diagonal initialisation (assuming system is dominated by diagonal of X) ####
# for signal deconvolution problem with n=p this gives the best initial solution
# only works if X is square though
# beta_init_diagonal <- y*sqrt(weights) / diag(Xw)


# coefficients & metrics over a range of lambda values
# calculate max lambda where all coefs would be penalized to 0
weights <- L0glm_gaus_presetlambda$weights
Xw <- Diagonal(x=sqrt(weights)) %*% data_gaus$X
sel_task <- which.max(colMeans(abs(as.matrix(data_gaus$y)))) # for multitask select column with max mean abs value?
yw <- as.matrix(data_gaus$y)[, sel_task,drop=F] * sqrt(weights) 
# TODO note that yw should really be zw, i.e. for non-identity link you should 
# be working on the adjusted response scale - maybe you should also make 
# L0glm.fit return z? Or make a convenience function that returns z and w 
# given a particular given family and coefficient estimates?
maxlambda <- max((1/colSums(Xw^2)) * (crossprod(Xw, yw)^2), na.rm=T)
if (family$family == "gaussian") maxlambda <- maxlambda/4
# NOTE: Frommlet & Nuel factor of 1/4 to convert lambda to lambda tilde, 
# NOTE: could be that this factor of 1/4 applies to all families - CHECK 
mincoef <- min(abs(L0glm_gaus_presetlambda$coefficients))

lambdas <- 10^seq(round(log10(lambda_preset_gaus))-3, log10(maxlambda)+1, by=0.1)
# lambdas <- 10^seq(log10(lambda_preset_gaus*0.9), log10(lambda_preset_gaus*1.1), by=0.005)
coefmat <- Matrix(0, nrow=length(lambdas), ncol=ncol(data_gaus$X))
warmstart <- 0.1 # % by which L0glm_gaus_presetlambda$coefficients will be used as start, lowest value that led to stable coefficient paths
system.time(metrics <- data.frame(do.call(rbind, lapply(1:length(lambdas), function (i) { 
  lambda <- lambdas[[i]]
  start <- (L0glm_gaus_presetlambda$coefficients)*warmstart+rep(mean(L0glm_gaus_presetlambda$coefficients), ncol(data_gaus$X))*(1-warmstart)
  # start[L0glm_gaus_presetlambda$coefficients==0] <- 0
  L0glm.out <- L0glm.fit(X = data_gaus$X, 
                         y = as.matrix(data_gaus$y), 
                         family = gaussian(identity), 
                         lambda = lambda, 
                         nonzero.centered = FALSE, # use nonzero centered ridge regression?
                         ridge.center = NULL, # L0glm_gaus_presetlambda$coefficients, # rep(0, ncol(data_gaus$X)),
                         nonnegative = FALSE,
                         warmupiters = 1,
                         warmupiters.nzc = 5,
                         start = rep(1, ncol(data_gaus$X)), # beta_init_orthogonal, # rep(1, ncol(data_gaus$X)), # beta_init_orthogonal, # start # rep(1, ncol(data_gaus$X)), # 
                         solver = "eigen", # eigen, glmnet or osqp, # note: constrained solvers glmnet & osqp better than eigen for very small lambda values
                         method = 3L, # 3L here 2x as fast than 7L
                         # gamma = 2, 
                         thresh = 1E-10,
                         relaxed = FALSE 
  )$coefficients 
  coefmat[i,] <<- as.vector(L0glm.out)
  return( unlist(c(lambda=lambda, get.metrics(beta = L0glm.out, 
                                              X = X, y = y, 
                                              family = gaussian(identity), 
                                              beta_true = data_gaus$beta_true))) )
} )))) # 42s
matplot(x=log10(metrics$lambda), 
        y=as.matrix(coefmat), type="l", lty=1,
        xlab="Log10(lambda)", ylab="Coefficient estimates", main="Regularisation path L0 Gaussian GLM")
library(scales)
abline(h=data_gaus$beta_true, col=alpha(1:6,0.2))
abline(v=log10(maxlambda*4), col="red") # max lambda where all coefs would be penalized to 0

plot(log10(metrics$lambda), metrics$deviance, type="l")
metrics$lambda[which.min(metrics$deviance)] # 0.001 # NOTE this is deviance relative to observed signal, not deviance relative to ground truth though
plot(log10(metrics$lambda), metrics$deviance_true, type="l")
metrics$lambda[which.min(metrics$deviance_true)] # 2.511886 # deviance relative to ground truth
plot(log10(metrics$lambda), metrics$loglik, type="l")
metrics$lambda[which.max(metrics$loglik)] # 0.3981072
plot(log10(metrics$lambda), metrics$MD_coefs, type="l")
plot(log10(metrics$lambda), metrics$MAD_coefs, type="l")
metrics$lambda[which.min(metrics$MAD_coefs)] # 1.995262
plot(log10(metrics$lambda), metrics$sMAPE_coefs, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(metrics$sMAPE_coefs)] # 3.162278
plot(log10(metrics$lambda), metrics$bias_coefs, type="l")
metrics$lambda[which.min(abs(metrics$bias_coefs-1))] # 0.001258925
plot(log10(metrics$lambda), metrics$R2_coefs, type="l")
metrics$lambda[which.max(abs(metrics$R2_coefs))] # 1.995262
plot(log10(metrics$lambda), metrics$estimated_over_true_support_size, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(abs(metrics$estimated_over_true_support_size-1))] # 1.258925
plot(log10(metrics$lambda), metrics$TP, type="l")
plot(log10(metrics$lambda), metrics$FP, type="l")
plot(log10(metrics$lambda), metrics$TN, type="l")
plot(log10(metrics$lambda), metrics$FN, type="l")
plot(log10(metrics$lambda), metrics$precision, type="l")
plot(log10(metrics$lambda), metrics$recall, type="l")
plot(log10(metrics$lambda), metrics$specificity, type="l")
plot(log10(metrics$lambda), metrics$accuracy, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.max(metrics$accuracy)] # 3.162278
plot(log10(metrics$lambda), metrics$F1_score, type="l")
metrics$lambda[which.max(metrics$F1_score)] # 3.162278
plot(log10(metrics$lambda), metrics$aic, type="l")
metrics$lambda[which.min(metrics$aic)] # 1
plot(log10(metrics$lambda), metrics$bic, type="l")
metrics$lambda[which.min(metrics$bic)] # 1.995262
plot(log10(metrics$lambda), metrics$ebic, type="l")
metrics$lambda[which.min(metrics$ebic)] # 1.995262
plot(log10(metrics$lambda), metrics$gic, type="l")
metrics$lambda[which.min(metrics$gic)] # 3.162278
plot(log10(metrics$lambda), metrics$mbic, type="l")
metrics$lambda[which.min(metrics$mbic)] # 3.162278

# reasonable match with lambda values that were expected to optimise each IC :
# note: using warm starts in the L0glm.fit function would shift the
# empirically observed optimal values

preset.lambda(IC = "aic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_gaus$y), 
              X = Matrix(data_gaus$X, sparse=TRUE), 
              family = gaussian(identity),
              dispersion = dispersion) # 0.5
preset.lambda(IC = "bic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_gaus$y), 
              X = Matrix(data_gaus$X, sparse=TRUE), 
              family = gaussian(identity),
              dispersion = dispersion) # 2.30
preset.lambda(IC = "ebic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_gaus$y), 
              X = Matrix(data_gaus$X, sparse=TRUE), 
              family = gaussian(identity),
              dispersion = dispersion) # 2.30
preset.lambda(IC = "gic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_gaus$y), 
              X = Matrix(data_gaus$X, sparse=TRUE), 
              family = gaussian(identity),
              dispersion = dispersion) # 2.56
preset.lambda(IC = "mbic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_gaus$y), 
              X = Matrix(data_gaus$X, sparse=TRUE), 
              family = gaussian(identity),
              dispersion = dispersion) # 3.92



# 2. POISSON LOG LINK GLM ####

set.seed(1)
data_pois  = simulate_glm(n = n, p = p, n_nonzero = n_nonzero, intercept = intercept, family = "poisson") 
X = data_pois$X
y = data_pois$y

# glm fit
system.time(glm_pois <- glm(data_pois$y ~ 0+data_pois$X, family=poisson(log))) # 0.41s
glm_pois$deviance 
glm_pois$aic
AIC(glm_pois)
BIC(glm_pois)
logLik(glm_pois)
cor(coef(glm_pois), data_pois$beta_true) # 0.9998249
plot(x=coef(glm_pois), y=data_pois$beta_true, pch=16, col="blue")

# Calculate L2 norm of columns
X.n <- sqrt(Matrix::colSums(X^2))
# X_normalized <- X %*% Matrix::Diagonal(x = 1 / X.n)
X_normalized <- sweep(X, MARGIN=2, X.n, "/")


# L0 glm fit with lambda close to zero
system.time(L0glm_pois_lambda0 <- L0glm.fit(X = X_normalized, 
                                            y = y, 
                                            family = poisson(log), 
                                            lambda = 1E-10, 
                                            nonzero.centered = FALSE, # use nonzero centered ridge regression?
                                            ridge.center = rep(0, ncol(X)),
                                            nonnegative = FALSE,
                                            warmupiters = 1,
                                            start = rep(1, ncol(data_pois$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                                            solver = "eigen", # 
                                            method = 7L, #
                                            thresh = 1E-10
))  
# TODO: start = rep(1, ncol(data_pois$X)) gives bad solution, rep(0.05, ncol(data_pois$X)) works OK
# CHECK how to fix this in an automatic way - 
# maybe scale columns of X or normalize penalty weights to sum to nrows(data)?

L0glm_pois_lambda0_coefs = L0glm_pois_lambda0$coefficients / X.n
max(abs(L0glm_pois_lambda0_coefs-coef(glm_pois))) # 8.484768e-12
# TODO gives bad solution with start=rep(1E-2, ncol(data_pois$X)), but 
# good solution with rep(1, ncol(data_pois$X)) - fix
plot(L0glm_pois_lambda0$coefficients, coef(glm_pois), pch=16, col="blue")

cor(L0glm_pois_lambda0$coefficients, data_pois$beta_true) # 0.9998249
plot(x=L0glm_pois_lambda0$coefficients, y=data_pois$beta_true, pch=16, col="blue")

# check if metrics (loglik & AIC & BIC) are correct
metrics_pois = data.frame(do.call(cbind, get.metrics(beta = L0glm_pois_lambda0$coefficients,  # some metrics
                                                     X = Matrix(data_pois$X, sparse=TRUE), 
                                                     y = as.matrix(data_pois$y), 
                                                     family = gaussian(identity), 
                                                     beta_true = as.matrix(data_pois$beta_true))))
metrics_pois
glm_pois$deviance-metrics_pois$deviance # does not match - TODO check source of discrepancy
logLik(glm_pois)-metrics_pois$loglik # does not match because R's glm function kicks out infinites & I don't
AIC(glm_pois)-metrics_pois$aic # 
BIC(glm_pois)-metrics_pois$bic # 


# L0 glm fit with lambda preset to optimise given IC
lambda_preset_pois = preset.lambda(IC = crit, # options aic, bic, gic, ebic, mbic 
                                   y = data_pois$y, 
                                   X = Matrix(data_pois$X, sparse=TRUE), 
                                   family = poisson(log),
                                   dispersion = dispersion) # dispersion here was 1, here we pass this value to not have to estimate it 
lambda_preset_pois # 2.561769

system.time(L0glm_pois_presetlambda <- L0glm.fit(X = Matrix(data_pois$X, sparse=TRUE), 
                                                 y = data_pois$y, 
                                                 family = poisson(log), 
                                                 lambda = lambda_preset_pois,
                                                 nonzero.centered = FALSE, # use nonzero centered ridge regression?
                                                 nonnegative = FALSE,
                                                 warmupiters = 1,
                                                 start = rep(0.05, ncol(data_pois$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                                                 solver = "eigen", # 
                                                 method = 1L, #
                                                 thresh = 1E-10
)) # 1.92s
mean(abs(L0glm_pois_presetlambda$coefficients-coef(glm_pois))) # 0.008970691
plot(L0glm_pois_presetlambda$coefficients, coef(glm_pois), pch=16, col="blue")

cor(L0glm_pois_presetlambda$coefficients, data_pois$beta_true) # 0.9996334
plot(x=L0glm_pois_presetlambda$coefficients, y=data_pois$beta_true, pch=16, col="blue")


# coefficients & metrics over a range of lambda values
# calculate max lambda where all coefs would be penalized to 0
weights <- L0glm_pois_presetlambda$weights
Xw <- Diagonal(x=sqrt(weights)) %*% data_pois$X
sel_task <- which.max(colMeans(abs(as.matrix(data_pois$y)))) # for multitask select column with max mean abs value?
yw <- as.matrix(data_pois$y)[, sel_task,drop=F] * sqrt(weights) 
# TODO note that yw should really be zw, i.e. for non-identity link you should 
# be working on the adjusted response scale - maybe you should also make L0glm.fit 
# return z? Or make a convenience function that returns z and w given a particular 
# given family and coefficient estimates?
maxlambda <- max((1/colSums(Xw^2)) * (crossprod(Xw, yw)^2), na.rm=T)
if (family$family == "gaussian") maxlambda <- maxlambda/4 # could be that this factor of 1/4 applies to all families & that this is that 1/4 factor mentioned in Frommlet & Nuel
mincoef <- min(abs(L0glm_pois_presetlambda$coefficients)) # *sign(L0glm_pois_presetlambda$coefficients)

lambdas <- 10^seq(round(log10(lambda_preset_pois))-3, log10(maxlambda)+1, by=0.2)
coefmat <- Matrix(0, nrow=length(lambdas), ncol=ncol(data_pois$X))
warmstart <- 0.1 # % by which L0glm_pois_presetlambda$coefficients will be used as start, lowest value that led to stable coefficient paths
system.time(metrics <- data.frame(do.call(rbind, lapply(1:length(lambdas), function (i) { 
  lambda <- lambdas[[i]]
  start <- L0glm_pois_presetlambda$coefficients*warmstart+rep(mincoef/2,ncol(data_pois$X))*(1-warmstart)
  # start[L0glm_pois_presetlambda$coefficients==0] <- 0
  L0glm.out <- L0glm.fit(X = Matrix(data_pois$X, sparse=TRUE), 
                         y = data_pois$y, 
                         family = poisson(log), 
                         lambda = lambda, 
                         nonzero.centered = FALSE, # use nonzero centered ridge regression? 
                         nonnegative = FALSE,
                         warmupiters = 1,
                         start = rep(0.05, ncol(data_pois$X)),
                         solver = "eigen", # eigen, glmnet or osqp, # note: constrained solvers glmnet & osqp better than eigen for very small lambda values
                         method = 7L, 
                         # gamma = 2, 
                         thresh = 1E-10,
                         relaxed = FALSE 
  )$coefficients 
  coefmat[i,] <<- as.vector(L0glm.out)
  return( unlist(c(lambda=lambda, get.metrics(beta = L0glm.out, 
                                              X = X, y = y, 
                                              family = poisson(log), 
                                              beta_true = data_pois$beta_true))) )
} )))) # 136s
matplot(x=log10(metrics$lambda), 
        y=as.matrix(coefmat), type="l", lty=1,
        xlab="Log10(lambda)", ylab="Coefficient estimates", main="Regularisation path L0 Gaussian GLM")
abline(v=log10(maxlambda), col="red") # max lambda where all coefs would be penalized to 0
# TODO calculate correctly using adjusted zw score instead of yw

plot(log10(metrics$lambda), metrics$deviance, type="l")
metrics$lambda[which.min(metrics$deviance)] # 0.001 # NOTE this is deviance relative to observed signal, not deviance relative to ground truth though
plot(log10(metrics$lambda), metrics$deviance_true, type="l")
metrics$lambda[which.min(metrics$deviance_true)] # 0.6309573 # deviance relative to ground truth
plot(log10(metrics$lambda), metrics$loglik, type="l")
metrics$lambda[which.max(metrics$loglik)] # 0.6309573
plot(log10(metrics$lambda), metrics$MD_coefs, type="l")
plot(log10(metrics$lambda), metrics$MAD_coefs, type="l")
metrics$lambda[which.min(metrics$MAD_coefs)] # 0.6309573
plot(log10(metrics$lambda), metrics$sMAPE_coefs, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(metrics$sMAPE_coefs)] # 0.6309573
plot(log10(metrics$lambda), metrics$bias_coefs, type="l")
metrics$lambda[which.min(abs(metrics$bias_coefs-1))] # 0.1
plot(log10(metrics$lambda), metrics$R2_coefs, type="l")
plot(log10(metrics$lambda), metrics$estimated_over_true_support_size, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(abs(metrics$estimated_over_true_support_size-1))] # 0.6309573
plot(log10(metrics$lambda), metrics$TP, type="l")
plot(log10(metrics$lambda), metrics$FP, type="l")
plot(log10(metrics$lambda), metrics$TN, type="l")
plot(log10(metrics$lambda), metrics$FN, type="l")
plot(log10(metrics$lambda), metrics$precision, type="l")
plot(log10(metrics$lambda), metrics$recall, type="l")
plot(log10(metrics$lambda), metrics$specificity, type="l")
plot(log10(metrics$lambda), metrics$accuracy, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.max(metrics$accuracy)] # 0.6309573
plot(log10(metrics$lambda), metrics$F1_score, type="l")
metrics$lambda[which.max(metrics$F1_score)] # 0.6309573
plot(log10(metrics$lambda), metrics$aic, type="l")
metrics$lambda[which.min(metrics$aic)] # 0.6309573
plot(log10(metrics$lambda), metrics$bic, type="l")
metrics$lambda[which.min(metrics$bic)] # 0.6309573
plot(log10(metrics$lambda), metrics$gic, type="l")
metrics$lambda[which.min(metrics$gic)] # 0.6309573
plot(log10(metrics$lambda), metrics$ebic, type="l")
metrics$lambda[which.min(metrics$ebic)] # 0.6309573
plot(log10(metrics$lambda), metrics$mbic, type="l")
metrics$lambda[which(metrics$mbic==min(metrics$mbic))] # 0.6309573

# reasonable match with lambda values that were expected to optimise each IC :
# note: using warm starts in the L0glm.fit function would shift the
# empirically observed optimal values

preset.lambda(IC = "aic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_pois$y), 
              X = Matrix(data_pois$X, sparse=TRUE), 
              family = poisson(log),
              dispersion = dispersion) # 0.5
preset.lambda(IC = "bic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_pois$y), 
              X = Matrix(data_pois$X, sparse=TRUE), 
              family = poisson(log),
              dispersion = dispersion) # 2.30
preset.lambda(IC = "gic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_pois$y), 
              X = Matrix(data_pois$X, sparse=TRUE), 
              family = poisson(log),
              dispersion = dispersion) # 2.56
preset.lambda(IC = "ebic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_pois$y), 
              X = Matrix(data_pois$X, sparse=TRUE), 
              family = poisson(log),
              dispersion = dispersion) # 2.30
preset.lambda(IC = "mbic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_pois$y), 
              X = Matrix(data_pois$X, sparse=TRUE), 
              family = poisson(log),
              dispersion = dispersion) # 3.91


# 3. BINOMIAL LOGIT LINK GLM ####

set.seed(1)
data_binom  = simulate_glm(n = n, p = p, n_nonzero = n_nonzero, intercept = intercept, family = "binomial", ntot = ntot)
X = data_binom$X
y = data_binom$y # proportions
weights = data_binom$weights # total sample size per row in dataset

# glm fit
system.time(glm_binom <- glm(data_binom$y ~ 0+data_binom$X, weights=data_binom$weights, family=binomial(logit))) # 0.37s
glm_binom$deviance 
glm_binom$aic
AIC(glm_binom)
BIC(glm_binom)
logLik(glm_binom)
cor(coef(glm_binom), data_binom$beta_true) # 0.9983967
plot(x=coef(glm_binom), y=data_binom$beta_true, pch=16, col="blue")

# L0 glm fit with lambda close to zero
system.time(L0glm_binom_lambda0 <- L0glm.fit.internal(X = X, 
                                            y = y, 
                                            weights = weights,
                                            family = binomial(logit), 
                                            lambda = 1E-10, 
                                            nonzero.centered = FALSE, # use nonzero centered ridge regression?
                                            ridge.center = rep(0, ncol(X)),
                                            nonnegative = FALSE,
                                            warmupiters = 1,
                                            start = rep(0.01, ncol(data_binom$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                                            no.pen = 1,
                                            solver = "eigen", # 
                                            method = 3L, #
                                            thresh = 1E-10,
                                            lower = rep(-Inf, ncol(X)),
                                            upper = rep(Inf, ncol(X))
))  # 0.22s
# TODO fix L0glm.fit gives error Error in linkinv(eta) : LENGTH or similar applied to S4 object
mean(abs(coef(glm_binom)-L0glm_binom_lambda0$coefficients)) # 2.523878e-11
plot(x=L0glm_binom_lambda0$coefficients, y=data_binom$beta_true, pch=16, col="blue")

cor(L0glm_binom_lambda0$coefficients, data_binom$beta_true) # 0.9983967
plot(x=L0glm_binom_lambda0$coefficients, y=data_binom$beta_true, pch=16, col="blue")

# check if metrics (loglik & AIC & BIC) are correct
metrics_binom = data.frame(do.call(cbind, get.metrics(beta = as.matrix(L0glm_binom_lambda0$coefficients),  # some metrics
                                                     X = Matrix(data_binom$X, sparse=TRUE), 
                                                     y = as.matrix(data_binom$y), 
                                                     weights = as.matrix(weights),
                                                     family = binomial(logit), 
                                                     beta_true = as.matrix(data_binom$beta_true) ) ) )
metrics_binom
glm_binom$deviance-metrics_binom$deviance # 0, correct
logLik(glm_binom)-metrics_binom$loglik # 0, correct
AIC(glm_binom)-metrics_binom$aic # 0, correct
BIC(glm_binom)-metrics_binom$bic # 0, correct


# L0 glm fit with lambda preset to optimise given IC
lambda_preset_binom = preset.lambda(IC = crit, # options aic, bic, gic, ebic, mbic 
                                   y = data_binom$y, 
                                   X = Matrix(data_binom$X, sparse=TRUE), 
                                   family = binomial(logit),
                                   dispersion = dispersion) # dispersion here was 1, here we pass this value to not have to estimate it 
lambda_preset_binom # 2.561769

system.time(L0glm_binom_presetlambda <- L0glm.fit.internal(X = X, 
                                                      y = y, 
                                                      weights = weights,
                                                      family = binomial(logit), 
                                                      lambda = lambda_preset_binom, 
                                                      nonzero.centered = FALSE, # use nonzero centered ridge regression?
                                                      ridge.center = rep(0, ncol(X)),
                                                      nonnegative = FALSE,
                                                      warmupiters = 1,
                                                      start = rep(0.01, ncol(data_binom$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                                                      no.pen = 1,
                                                      solver = "eigen", # 
                                                      method = 3L, #
                                                      thresh = 1E-10,
                                                      lower = rep(-Inf, ncol(X)),
                                                      upper = rep(Inf, ncol(X))
)) # 0.36s
# TODO fix L0glm.fit gives error Error in linkinv(eta) : LENGTH or similar applied to S4 object
mean(abs(L0glm_binom_presetlambda$coefficients-coef(glm_binom))) # 0.03297451 
plot(L0glm_binom_presetlambda$coefficients, coef(glm_binom), pch=16, col="blue")

cor(L0glm_binom_presetlambda$coefficients, data_binom$beta_true) # 0.9967551
plot(x=L0glm_binom_presetlambda$coefficients, y=data_binom$beta_true, pch=16, col="blue")


# coefficients & metrics over a range of lambda values
# calculate max lambda where all coefs would be penalized to 0
weights <- L0glm_binom_presetlambda$weights
Xw <- Diagonal(x=sqrt(weights)) %*% data_binom$X
sel_task <- which.max(colMeans(abs(as.matrix(data_binom$y)))) # for multitask select column with max mean abs value?
yw <- as.matrix(data_binom$y)[, sel_task,drop=F] * sqrt(weights) 
# TODO note that yw should really be zw, i.e. for non-identity link you should 
# be working on the adjusted response scale - maybe you should also make L0glm.fit 
# return z? Or make a convenience function that returns z and w given a particular 
# given family and coefficient estimates?
maxlambda <- max((1/colSums(Xw^2)) * (crossprod(Xw, yw)^2), na.rm=T)
if (family$family == "gaussian") maxlambda <- maxlambda/4 # could be that this factor of 1/4 applies to all families & that this is that 1/4 factor mentioned in Frommlet & Nuel
mincoef <- min(abs(L0glm_binom_presetlambda$coefficients)) # *sign(L0glm_binom_presetlambda$coefficients)

lambdas <- 10^seq(round(log10(lambda_preset_binom))-3, log10(maxlambda)+1, by=0.2)
coefmat <- Matrix(0, nrow=length(lambdas), ncol=ncol(data_binom$X))
warmstart <- 0.1 # % by which L0glm_binom_presetlambda$coefficients will be used as start, lowest value that led to stable coefficient paths
system.time(metrics <- data.frame(do.call(rbind, lapply(1:length(lambdas), function (i) { 
  lambda <- lambdas[[i]]
  start <- L0glm_binom_presetlambda$coefficients*warmstart+rep(mincoef/2,ncol(data_binom$X))*(1-warmstart)
  # start[L0glm_binom_presetlambda$coefficients==0] <- 0
  L0glm.out <- L0glm.fit(X = Matrix(data_binom$X, sparse=TRUE), 
                         y = data_binom$y, 
                         family = binomial(logit), 
                         lambda = lambda, 
                         nonzero.centered = FALSE, # use nonzero centered ridge regression? 
                         nonnegative = FALSE,
                         warmupiters = 1,
                         start = rep(0.05, ncol(data_binom$X)),
                         solver = "eigen", # eigen, glmnet or osqp, # note: constrained solvers glmnet & osqp better than eigen for very small lambda values
                         method = 7L, 
                         # gamma = 2, 
                         thresh = 1E-10,
                         relaxed = FALSE 
  )$coefficients 
  coefmat[i,] <<- as.vector(L0glm.out)
  return( unlist(c(lambda=lambda, get.metrics(beta = L0glm.out, 
                                              X = X, y = y, 
                                              family = binomial(logit), 
                                              beta_true = data_binom$beta_true))) )
} )))) # 91s
matplot(x=log10(metrics$lambda), 
        y=as.matrix(coefmat), type="l", lty=1,
        xlab="Log10(lambda)", ylab="Coefficient estimates", main="Regularisation path L0 Gaussian GLM")
abline(v=log10(maxlambda), col="red") # max lambda where all coefs would be penalized to 0
# TODO calculate correctly using adjusted zw score instead of yw

plot(log10(metrics$lambda), metrics$deviance, type="l")
metrics$lambda[which.min(metrics$deviance)] # 0.03981072 # NOTE this is deviance relative to observed signal, not deviance relative to ground truth though
plot(log10(metrics$lambda), metrics$loglik, type="l")
metrics$lambda[which.max(metrics$loglik)] # 0.3981072
plot(log10(metrics$lambda), metrics$MD_coefs, type="l")
plot(log10(metrics$lambda), metrics$MAD_coefs, type="l")
metrics$lambda[which.min(metrics$MAD_coefs)] # 1
plot(log10(metrics$lambda), metrics$sMAPE_coefs, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(metrics$sMAPE_coefs)] # 1
plot(log10(metrics$lambda), metrics$bias_coefs, type="l")
metrics$lambda[which.min(abs(metrics$bias_coefs-1))] # 0.6309573
plot(log10(metrics$lambda), metrics$R2_coefs, type="l")
plot(log10(metrics$lambda), metrics$estimated_over_true_support_size, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(abs(metrics$estimated_over_true_support_size-1))] # 1
plot(log10(metrics$lambda), metrics$TP, type="l")
plot(log10(metrics$lambda), metrics$FP, type="l")
plot(log10(metrics$lambda), metrics$TN, type="l")
plot(log10(metrics$lambda), metrics$FN, type="l")
plot(log10(metrics$lambda), metrics$precision, type="l")
plot(log10(metrics$lambda), metrics$recall, type="l")
plot(log10(metrics$lambda), metrics$specificity, type="l")
plot(log10(metrics$lambda), metrics$accuracy, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.max(metrics$accuracy)] # 1
plot(log10(metrics$lambda), metrics$F1_score, type="l")
metrics$lambda[which.max(metrics$F1_score)] # 1
plot(log10(metrics$lambda), metrics$aic, type="l")
metrics$lambda[which.min(metrics$aic)] # 1
plot(log10(metrics$lambda), metrics$bic, type="l")
metrics$lambda[which.min(metrics$bic)] # 1
plot(log10(metrics$lambda), metrics$gic, type="l")
metrics$lambda[which.min(metrics$gic)] # 1
plot(log10(metrics$lambda), metrics$ebic, type="l")
metrics$lambda[which.min(metrics$ebic)] # 1
plot(log10(metrics$lambda), metrics$mbic, type="l")
metrics$lambda[which(metrics$mbic==min(metrics$mbic))] # 1

# reasonable match with lambda values that were expected to optimise each IC :
# note: using warm starts in the L0glm.fit function would shift the
# empirically observed optimal values

preset.lambda(IC = "aic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_binom$y), 
              X = Matrix(data_binom$X, sparse=TRUE), 
              family = binomial(logit),
              dispersion = dispersion) # 0.5
preset.lambda(IC = "bic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_binom$y), 
              X = Matrix(data_binom$X, sparse=TRUE), 
              family = binomial(logit),
              dispersion = dispersion) # 2.30
preset.lambda(IC = "gic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_binom$y), 
              X = Matrix(data_binom$X, sparse=TRUE), 
              family = binomial(logit),
              dispersion = dispersion) # 2.56
preset.lambda(IC = "ebic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_binom$y), 
              X = Matrix(data_binom$X, sparse=TRUE), 
              family = binomial(logit),
              dispersion = dispersion) # 2.30
preset.lambda(IC = "mbic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_binom$y), 
              X = Matrix(data_binom$X, sparse=TRUE), 
              family = binomial(logit),
              dispersion = dispersion) # 3.91




# 4. GAMMA INVERSE LINK GLM ####

set.seed(1)
data_gamma  = simulate_glm(n, p, n_nonzero, family = "gamma") 
X = data_gamma$X
y = data_gamma$y

# glm fit
system.time(glm_gamma <- glm(data_gamma$y ~ 0+data_gamma$X, family=Gamma(inverse))) # 0.45s
glm_gamma$deviance 
glm_gamma$aic
AIC(glm_gamma)
BIC(glm_gamma)
logLik(glm_gamma)
cor(coef(glm_gamma), data_gamma$beta_true)
plot(x=coef(glm_gamma), y=data_gamma$beta_true, pch=16, col="blue")


# L0 glm fit with lambda close to zero
system.time(L0glm_gamma_lambda0 <- L0glm.fit(X = X, 
                                            y = y, 
                                            family = Gamma(inverse), 
                                            lambda = 1E-10, 
                                            nonzero.centered = FALSE, # use nonzero centered ridge regression?
                                            ridge.center = rep(0, ncol(X)),
                                            nonnegative = FALSE,
                                            warmupiters = 1,
                                            start = rep(0.05, ncol(data_gamma$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                                            solver = "eigen", # 
                                            method = 7L, #
                                            thresh = 1E-10
))  # 1.55s
# TODO: start = rep(1, ncol(data_gamma$X)) gives bad solution, rep(1E-2, ncol(data_gamma$X)) OK

max(abs(L0glm_gamma_lambda0$coefficients-coef(glm_gamma))) # 2.121917e-08, TODO gives bad solution with start=rep(1E-2, ncol(data_gamma$X)), but good solution with rep(1, ncol(data_gamma$X)) - fix
plot(L0glm_gamma_lambda0$coefficients, coef(glm_gamma), pch=16, col="blue")
# check if metrics (loglik & AIC & BIC) are correct
metrics_gamma = data.frame(do.call(cbind, get.metrics(beta = L0glm_gamma_lambda0$coefficients,  # some metrics
                                                     X = Matrix(data_gamma$X, sparse=TRUE), 
                                                     y = as.matrix(data_gamma$y), 
                                                     family = Gamma(inverse), 
                                                     beta_true = as.matrix(data_gamma$beta_true))))
metrics_gamma
logLik(glm_gamma)-metrics_gamma$loglik # 0, correct
AIC(glm_gamma)-metrics_gamma$aic # 0
BIC(glm_gamma)-metrics_gamma$bic # 0


# L0 glm fit with lambda preset to optimise given IC
lambda_preset_gamma = preset.lambda(IC = crit, # options aic, bic, gic, ebic, mbic 
                                   y = data_gamma$y, 
                                   X = Matrix(data_gamma$X, sparse=TRUE), 
                                   family = Gamma(inverse),
                                   dispersion = dispersion) # dispersion here was 1, here we pass this value to not have to estimate it 
lambda_preset_gamma # 2.556246

system.time(L0glm_gamma_presetlambda <- L0glm.fit(X = Matrix(data_gamma$X, sparse=TRUE), 
                                                 y = data_gamma$y, 
                                                 family = Gamma(inverse), 
                                                 lambda = lambda_preset_gamma,
                                                 nonzero.centered = FALSE, # use nonzero centered ridge regression?
                                                 nonnegative = FALSE,
                                                 warmupiters = 1,
                                                 start = rep(0.05, ncol(data_gamma$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                                                 solver = "eigen", # 
                                                 method = 1L, #
                                                 thresh = 1E-10
)) # 1.62s
mean(abs(L0glm_gamma_presetlambda$coefficients-coef(glm_gamma))) # 0.5051875
plot(L0glm_gamma_presetlambda$coefficients, coef(glm_gamma), pch=16, col="blue")
# 


# coefficients & metrics over a range of lambda values
# calculate max lambda where all coefs would be penalized to 0
weights <- L0glm_gamma_presetlambda$weights
Xw <- Diagonal(x=sqrt(weights)) %*% data_gamma$X
sel_task <- which.max(colMeans(abs(as.matrix(data_gamma$y)))) # for multitask select column with max mean abs value?
yw <- as.matrix(data_gamma$y)[, sel_task,drop=F] * sqrt(weights) 
# TODO note that yw should really be zw, i.e. for non-identity link you should 
# be working on the adjusted response scale - maybe you should also make L0glm.fit 
# return z? Or make a convenience function that returns z and w given a particular 
# given family and coefficient estimates?
maxlambda <- max((1/colSums(Xw^2)) * (crossprod(Xw, yw)^2), na.rm=T)
if (family$family == "gaussian") maxlambda <- maxlambda/4 # could be that this factor of 1/4 applies to all families & that this is that 1/4 factor mentioned in Frommlet & Nuel
mincoef <- min(abs(L0glm_gamma_presetlambda$coefficients)) # *sign(L0glm_gamma_presetlambda$coefficients)

lambdas <- 10^seq(round(log10(lambda_preset_gamma))-3, log10(maxlambda)+3, by=0.2)
coefmat <- Matrix(0, nrow=length(lambdas), ncol=ncol(data_gamma$X))
warmstart <- 0.1 # % by which L0glm_gamma_presetlambda$coefficients will be used as start, lowest value that led to stable coefficient paths
system.time(metrics <- data.frame(do.call(rbind, lapply(1:length(lambdas), function (i) { 
  lambda <- lambdas[[i]]
  start <- L0glm_gamma_presetlambda$coefficients*warmstart+rep(mincoef/2,ncol(data_gamma$X))*(1-warmstart)
  # start[L0glm_gamma_presetlambda$coefficients==0] <- 0
  L0glm.out <- L0glm.fit(X = Matrix(data_gamma$X, sparse=TRUE), 
                         y = data_gamma$y, 
                         family = Gamma(inverse), 
                         lambda = lambda, 
                         nonzero.centered = FALSE, # use nonzero centered ridge regression? 
                         nonnegative = FALSE,
                         warmupiters = 1,
                         start = rep(0.05, ncol(data_gamma$X)),
                         solver = "eigen", # eigen, glmnet or osqp, # note: constrained solvers glmnet & osqp better than eigen for very small lambda values
                         method = 7L, 
                         # gamma = 2, 
                         thresh = 1E-10,
                         relaxed = FALSE 
  )$coefficients 
  coefmat[i,] <<- as.vector(L0glm.out)
  return( unlist(c(lambda=lambda, get.metrics(beta = L0glm.out, 
                                              X = X, y = y, 
                                              family = Gamma(inverse), 
                                              beta_true = data_gamma$beta_true))) )
} )))) # 54.06s
matplot(x=log10(metrics$lambda), 
        y=as.matrix(coefmat), type="l", lty=1,
        xlab="Log10(lambda)", ylab="Coefficient estimates", main="Regularisation path L0 Gaussian GLM")
abline(v=log10(maxlambda), col="red") # max lambda where all coefs would be penalized to 0
# TODO calculate correctly using adjusted zw score instead of yw

plot(log10(metrics$lambda), metrics$deviance, type="l")
metrics$lambda[which.min(metrics$deviance)] # 0.001 # NOTE this is deviance relative to observed signal, not deviance relative to ground truth though
plot(log10(metrics$lambda), metrics$deviance_true, type="l")
metrics$lambda[which.min(metrics$deviance_true)] # 0.1584893 # deviance compared to true signal
plot(log10(metrics$lambda), metrics$loglik, type="l")
metrics$lambda[which.max(metrics$loglik)] # 0.6309573
plot(log10(metrics$lambda), metrics$loglik_true, type="l")
metrics$lambda[which.max(metrics$loglik_true)] # 0.1584893
plot(log10(metrics$lambda), metrics$MD_coefs, type="l")
plot(log10(metrics$lambda), metrics$MAD_coefs, type="l")
metrics$lambda[which.min(metrics$MAD_coefs)] # 0.2511886
plot(log10(metrics$lambda), metrics$sMAPE_coefs, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(metrics$sMAPE_coefs)] # 0.2511886
plot(log10(metrics$lambda), metrics$bias_coefs, type="l")
metrics$lambda[which.min(abs(metrics$bias_coefs-1))] # 0.1584893
plot(log10(metrics$lambda), metrics$R2_coefs, type="l")
plot(log10(metrics$lambda), metrics$estimated_over_true_support_size, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(abs(metrics$estimated_over_true_support_size-1))] # 0.06309573
plot(log10(metrics$lambda), metrics$TP, type="l")
plot(log10(metrics$lambda), metrics$FP, type="l")
plot(log10(metrics$lambda), metrics$TN, type="l")
plot(log10(metrics$lambda), metrics$FN, type="l")
plot(log10(metrics$lambda), metrics$precision, type="l")
plot(log10(metrics$lambda), metrics$recall, type="l")
plot(log10(metrics$lambda), metrics$specificity, type="l")
plot(log10(metrics$lambda), metrics$accuracy, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.max(metrics$accuracy)] # 0.2511886
plot(log10(metrics$lambda), metrics$F1_score, type="l")
metrics$lambda[which.max(metrics$F1_score)] # 0.1
plot(log10(metrics$lambda), metrics$aic, type="l")
metrics$lambda[which.min(metrics$aic)] # 1.584893
plot(log10(metrics$lambda), metrics$bic, type="l")
metrics$lambda[which.min(metrics$bic)] # 3.981072
plot(log10(metrics$lambda), metrics$gic, type="l")
metrics$lambda[which.min(metrics$gic)] # 3.981072
plot(log10(metrics$lambda), metrics$ebic, type="l")
metrics$lambda[which.min(metrics$ebic)] # 3.981072
plot(log10(metrics$lambda), metrics$mbic, type="l")
metrics$lambda[which(metrics$mbic==min(metrics$mbic))] # 3.981072

# reasonable match with lambda values that were expected to optimise each IC :
# note: using warm starts in the L0glm.fit function would shift the
# empirically observed optimal values

preset.lambda(IC = "aic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_gamma$y), 
              X = Matrix(data_gamma$X, sparse=TRUE), 
              family = poisson(log),
              dispersion = dispersion) # 0.5
preset.lambda(IC = "bic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_gamma$y), 
              X = Matrix(data_gamma$X, sparse=TRUE), 
              family = poisson(log),
              dispersion = dispersion) # 2.30
preset.lambda(IC = "gic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_gamma$y), 
              X = Matrix(data_gamma$X, sparse=TRUE), 
              family = poisson(log),
              dispersion = dispersion) # 2.56
preset.lambda(IC = "ebic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_gamma$y), 
              X = Matrix(data_gamma$X, sparse=TRUE), 
              family = poisson(log),
              dispersion = dispersion) # 2.30
preset.lambda(IC = "mbic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_gamma$y), 
              X = Matrix(data_gamma$X, sparse=TRUE), 
              family = poisson(log),
              dispersion = dispersion) # 3.91


# 5. INVERSE GAUSSIAN 1/mu^2 GLM ####

set.seed(1)
data_invgaus  = simulate_glm(n, p, n_nonzero, family = "inverse.gaussian") 
X = data_invgaus$X
y = data_invgaus$y

# glm fit
system.time(glm_invgaus <- glm(data_invgaus$y ~ 0+data_invgaus$X, 
                               family=inverse.gaussian(link="1/mu^2"))) # 0.45s
# TODO: check logic of simulation for this family
# now returns error Error: no valid set of coefficients has been found: please supply starting values
glm_invgaus$deviance 
glm_invgaus$aic
AIC(glm_invgaus)
BIC(glm_invgaus)
logLik(glm_invgaus)
cor(coef(glm_invgaus), data_invgaus$beta_true)
plot(x=coef(glm_invgaus), y=data_invgaus$beta_true, pch=16, col="blue")


# L0 glm fit with lambda close to zero
system.time(L0glm_invgaus_lambda0 <- L0glm.fit(X = X, 
                                             y = y, 
                                             family = Gamma(inverse), 
                                             lambda = 1E-10, 
                                             nonzero.centered = FALSE, # use nonzero centered ridge regression?
                                             ridge.center = rep(0, ncol(X)),
                                             nonnegative = FALSE,
                                             warmupiters = 1,
                                             start = rep(0.05, ncol(data_invgaus$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                                             solver = "eigen", # 
                                             method = 7L, #
                                             thresh = 1E-10
))  # 1.55s
# TODO: start = rep(1, ncol(data_invgaus$X)) gives bad solution, rep(1E-2, ncol(data_invgaus$X)) OK

max(abs(L0glm_invgaus_lambda0$coefficients-coef(glm_invgaus))) # 2.121917e-08, TODO gives bad solution with start=rep(1E-2, ncol(data_invgaus$X)), but good solution with rep(1, ncol(data_invgaus$X)) - fix
plot(L0glm_invgaus_lambda0$coefficients, coef(glm_invgaus), pch=16, col="blue")
# check if metrics (loglik & AIC & BIC) are correct
metrics_invgaus = data.frame(do.call(cbind, get.metrics(beta = L0glm_invgaus_lambda0$coefficients,  # some metrics
                                                      X = Matrix(data_invgaus$X, sparse=TRUE), 
                                                      y = as.matrix(data_invgaus$y), 
                                                      family = Gamma(inverse), 
                                                      beta_true = as.matrix(data_invgaus$beta_true))))
metrics_invgaus
logLik(glm_invgaus)-metrics_invgaus$loglik # 0, correct
AIC(glm_invgaus)-metrics_invgaus$aic # 0
BIC(glm_invgaus)-metrics_invgaus$bic # 0


# L0 glm fit with lambda preset to optimise given IC
lambda_preset_invgaus = preset.lambda(IC = crit, # options aic, bic, gic, ebic, mbic 
                                    y = data_invgaus$y, 
                                    X = Matrix(data_invgaus$X, sparse=TRUE), 
                                    family = Gamma(inverse),
                                    dispersion = dispersion) # dispersion here was 1, here we pass this value to not have to estimate it 
lambda_preset_invgaus # 2.556246

system.time(L0glm_invgaus_presetlambda <- L0glm.fit(X = Matrix(data_invgaus$X, sparse=TRUE), 
                                                  y = data_invgaus$y, 
                                                  family = Gamma(inverse), 
                                                  lambda = lambda_preset_invgaus,
                                                  nonzero.centered = FALSE, # use nonzero centered ridge regression?
                                                  nonnegative = FALSE,
                                                  warmupiters = 1,
                                                  start = rep(0.05, ncol(data_invgaus$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                                                  solver = "eigen", # 
                                                  method = 1L, #
                                                  thresh = 1E-10
)) # 1.62s
mean(abs(L0glm_invgaus_presetlambda$coefficients-coef(glm_invgaus))) # 0.5051875
plot(L0glm_invgaus_presetlambda$coefficients, coef(glm_invgaus), pch=16, col="blue")
# 


# coefficients & metrics over a range of lambda values
# calculate max lambda where all coefs would be penalized to 0
weights <- L0glm_invgaus_presetlambda$weights
Xw <- Diagonal(x=sqrt(weights)) %*% data_invgaus$X
sel_task <- which.max(colMeans(abs(as.matrix(data_invgaus$y)))) # for multitask select column with max mean abs value?
yw <- as.matrix(data_invgaus$y)[, sel_task,drop=F] * sqrt(weights) 
# TODO note that yw should really be zw, i.e. for non-identity link you should 
# be working on the adjusted response scale - maybe you should also make L0glm.fit 
# return z? Or make a convenience function that returns z and w given a particular 
# given family and coefficient estimates?
maxlambda <- max((1/colSums(Xw^2)) * (crossprod(Xw, yw)^2), na.rm=T)
if (family$family == "gaussian") maxlambda <- maxlambda/4 # could be that this factor of 1/4 applies to all families & that this is that 1/4 factor mentioned in Frommlet & Nuel
mincoef <- min(abs(L0glm_invgaus_presetlambda$coefficients)) # *sign(L0glm_invgaus_presetlambda$coefficients)

lambdas <- 10^seq(round(log10(lambda_preset_invgaus))-3, log10(maxlambda)+3, by=0.2)
coefmat <- Matrix(0, nrow=length(lambdas), ncol=ncol(data_invgaus$X))
warmstart <- 0.1 # % by which L0glm_invgaus_presetlambda$coefficients will be used as start, lowest value that led to stable coefficient paths
system.time(metrics <- data.frame(do.call(rbind, lapply(1:length(lambdas), function (i) { 
  lambda <- lambdas[[i]]
  start <- L0glm_invgaus_presetlambda$coefficients*warmstart+rep(mincoef/2,ncol(data_invgaus$X))*(1-warmstart)
  # start[L0glm_invgaus_presetlambda$coefficients==0] <- 0
  L0glm.out <- L0glm.fit(X = Matrix(data_invgaus$X, sparse=TRUE), 
                         y = data_invgaus$y, 
                         family = Gamma(inverse), 
                         lambda = lambda, 
                         nonzero.centered = FALSE, # use nonzero centered ridge regression? 
                         nonnegative = FALSE,
                         warmupiters = 1,
                         start = rep(0.05, ncol(data_invgaus$X)),
                         solver = "eigen", # eigen, glmnet or osqp, # note: constrained solvers glmnet & osqp better than eigen for very small lambda values
                         method = 7L, 
                         # gamma = 2, 
                         thresh = 1E-10,
                         relaxed = FALSE 
  )$coefficients 
  coefmat[i,] <<- as.vector(L0glm.out)
  return( unlist(c(lambda=lambda, get.metrics(beta = L0glm.out, 
                                              X = X, y = y, 
                                              family = Gamma(inverse), 
                                              beta_true = data_invgaus$beta_true))) )
} )))) # 54.06s
matplot(x=log10(metrics$lambda), 
        y=as.matrix(coefmat), type="l", lty=1,
        xlab="Log10(lambda)", ylab="Coefficient estimates", main="Regularisation path L0 Gaussian GLM")
abline(v=log10(maxlambda), col="red") # max lambda where all coefs would be penalized to 0
# TODO calculate correctly using adjusted zw score instead of yw

plot(log10(metrics$lambda), metrics$deviance, type="l")
metrics$lambda[which.min(metrics$deviance)] # 0.001 # NOTE this is deviance relative to observed signal, not deviance relative to ground truth though
plot(log10(metrics$lambda), metrics$deviance_true, type="l")
metrics$lambda[which.min(metrics$deviance_true)] # 0.1584893 # deviance compared to true signal
plot(log10(metrics$lambda), metrics$loglik, type="l")
metrics$lambda[which.max(metrics$loglik)] # 0.6309573
plot(log10(metrics$lambda), metrics$loglik_true, type="l")
metrics$lambda[which.max(metrics$loglik_true)] # 0.1584893
plot(log10(metrics$lambda), metrics$MD_coefs, type="l")
plot(log10(metrics$lambda), metrics$MAD_coefs, type="l")
metrics$lambda[which.min(metrics$MAD_coefs)] # 0.2511886
plot(log10(metrics$lambda), metrics$sMAPE_coefs, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(metrics$sMAPE_coefs)] # 0.2511886
plot(log10(metrics$lambda), metrics$bias_coefs, type="l")
metrics$lambda[which.min(abs(metrics$bias_coefs-1))] # 0.1584893
plot(log10(metrics$lambda), metrics$R2_coefs, type="l")
plot(log10(metrics$lambda), metrics$estimated_over_true_support_size, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(abs(metrics$estimated_over_true_support_size-1))] # 0.06309573
plot(log10(metrics$lambda), metrics$TP, type="l")
plot(log10(metrics$lambda), metrics$FP, type="l")
plot(log10(metrics$lambda), metrics$TN, type="l")
plot(log10(metrics$lambda), metrics$FN, type="l")
plot(log10(metrics$lambda), metrics$precision, type="l")
plot(log10(metrics$lambda), metrics$recall, type="l")
plot(log10(metrics$lambda), metrics$specificity, type="l")
plot(log10(metrics$lambda), metrics$accuracy, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.max(metrics$accuracy)] # 0.2511886
plot(log10(metrics$lambda), metrics$F1_score, type="l")
metrics$lambda[which.max(metrics$F1_score)] # 0.1
plot(log10(metrics$lambda), metrics$aic, type="l")
metrics$lambda[which.min(metrics$aic)] # 1.584893
plot(log10(metrics$lambda), metrics$bic, type="l")
metrics$lambda[which.min(metrics$bic)] # 3.981072
plot(log10(metrics$lambda), metrics$gic, type="l")
metrics$lambda[which.min(metrics$gic)] # 3.981072
plot(log10(metrics$lambda), metrics$ebic, type="l")
metrics$lambda[which.min(metrics$ebic)] # 3.981072
plot(log10(metrics$lambda), metrics$mbic, type="l")
metrics$lambda[which(metrics$mbic==min(metrics$mbic))] # 3.981072

# reasonable match with lambda values that were expected to optimise each IC :
# note: using warm starts in the L0glm.fit function would shift the
# empirically observed optimal values

preset.lambda(IC = "aic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_invgaus$y), 
              X = Matrix(data_invgaus$X, sparse=TRUE), 
              family = poisson(log),
              dispersion = dispersion) # 0.5
preset.lambda(IC = "bic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_invgaus$y), 
              X = Matrix(data_invgaus$X, sparse=TRUE), 
              family = poisson(log),
              dispersion = dispersion) # 2.30
preset.lambda(IC = "gic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_invgaus$y), 
              X = Matrix(data_invgaus$X, sparse=TRUE), 
              family = poisson(log),
              dispersion = dispersion) # 2.56
preset.lambda(IC = "ebic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_invgaus$y), 
              X = Matrix(data_invgaus$X, sparse=TRUE), 
              family = poisson(log),
              dispersion = dispersion) # 2.30
preset.lambda(IC = "mbic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_invgaus$y), 
              X = Matrix(data_invgaus$X, sparse=TRUE), 
              family = poisson(log),
              dispersion = dispersion) # 3.91