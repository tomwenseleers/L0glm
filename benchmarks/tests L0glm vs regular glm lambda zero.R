# Check of metrics (AIC, BIC, loglik) L0glm package & comparison regular GLM ####

library(remotes)
# remotes::install_github("tomwenseleers/L0glm/L0glm")
# library(L0glm)
devtools::load_all(".")
library(export)
source("..//benchmarks/generate data.R")

# check consistency of L0glm output with glm for gaussian identity link, poisson log link, Gamma log link & binomial logit link
n = 10000
p = 100
n_nonzero = 50
dispersion = 1
crit = "gic"

# 1. GAUSSIAN GLM ####

set.seed(1)
data_gaus  = simulate_glm(n, p, n_nonzero, family = "gaussian") 
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

# L0 glm fit with lambda close to zero
system.time(L0glm_gaus_lambda0 <- L0glm.fit(X = X, 
                                            y = y, 
                                            family = gaussian(identity), 
                                            lambda = 1E-30, 
                                            nonzero.centered = FALSE, # use nonzero centered ridge regression?
                                            nonnegative = FALSE,
                                            warmupiters = 1,
                                            # start = rep(1, ncol(data_gaus$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                                            solver = "eigen", # 
                                            method = 3L, 
                                            relaxed = FALSE,
                                            thresh = 1E-10
))  # 0.06s

# NOTE: note that L0glm can fit GLMs with perfectly collinear variables, e.g.
# with X = cbind(data_gaus$X, data_gaus$X)
# in that case, it would divide the total effect over the different collinear versions
# of each variable

max(abs(L0glm_gaus_lambda0$coefficients-coef(glm_gaus))) # 1.476034e-09, correct

# check if metrics (loglik & AIC & BIC) are correct
metrics_gaus = data.frame(do.call(cbind, get.metrics(beta = L0glm_gaus_lambda0$coefficients,  # some metrics
                                                     X = data_gaus$X, 
                                                     y = data_gaus$y,
                                                     family = gaussian(identity), 
                                                     beta_true = as.matrix(data_gaus$beta_true))))
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
lambda_preset_gaus # 2.556246

system.time(L0glm_gaus_presetlambda <- L0glm.fit(X = data_gaus$X, 
                                            y = data_gaus$y, 
                                            family = gaussian(identity), 
                                            lambda = lambda_preset_gaus,
                                            nonzero.centered = FALSE, # use nonzero centered ridge regression?
                                            # ridge.center = XX,
                                            nonnegative = FALSE,
                                            warmupiters = 1,
                                            # start = rep(1, ncol(data_gaus$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                                            solver = "eigen", # 
                                            method = 3L, #
                                            relaxed = TRUE,
                                            thresh = 1E-10
)) # 0.77s
sum(L0glm_gaus_presetlambda$coefficients!=0) 
# with nonzero.centered=FALSE 47, so OK
# TODO: check : not OK with nonzero.centered=TRUE - then gives 100, both if 
# start values are or are not given, ie a non-sparse solution
# when before solution was sparse also with nonzero.centered=TRUE

# couple of possible coefficient initialisations to use as starting values in L0glm.fit which could be used for start ####
weights = L0glm_gaus_presetlambda$weights
Xw = Matrix(Diagonal(x=sqrt(weights)) %*% data_gaus$X, sparse=F) # = X * sqrt(weights)
Xwaug = Matrix(rbind(Xw, Diagonal(x=sqrt(rep(1,ncol(data_gaus$X))))), sparse=F) # row augmented covariate matrix
ywaug = as.vector(c(data_gaus$y, rep(0, ncol(data_gaus$X)))) # 

# orthogonal initialisation (assuming variables are orthogonal) ####
beta_init_orthogonal <- as.vector(crossprod(Xwaug, ywaug) / colSums(Xwaug^2)) 
cor(beta_init_orthogonal, data_gaus$beta_true) # 0.9923448
plot(beta_init_orthogonal, data_gaus$beta_true, pch=16, col="blue")
mean(abs(beta_init_orthogonal-data_gaus$beta_true)) # 1.629831
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
# NOTE: could be that this factor of 1/4 applies to all families & that this is that 1/4 factor mentioned in Frommlet & Nuel
mincoef <- min(abs(L0glm_gaus_presetlambda$coefficients))

lambdas <- 10^seq(round(log10(lambda_preset_gaus))-3, log10(maxlambda)+1, by=0.2)
coefmat <- Matrix(0, nrow=length(lambdas), ncol=ncol(data_gaus$X))
warmstart <- 0.01 # % by which L0glm_gaus_presetlambda$coefficients will be used as start, lowest value that led to stable coefficient paths
system.time(metrics <- data.frame(do.call(rbind, lapply(1:length(lambdas), function (i) { 
  lambda <- lambdas[[i]]
  start <- (L0glm_gaus_presetlambda$coefficients)*warmstart+rep(1, ncol(data_gaus$X))*(1-warmstart)
  # start[L0glm_gaus_presetlambda$coefficients==0] <- 0
  L0glm.out <- L0glm.fit(X = data_gaus$X, 
                         y = as.matrix(data_gaus$y), 
                         family = gaussian(identity), 
                         lambda = lambda, 
                         nonzero.centered = FALSE, # use nonzero centered ridge regression?
                         # ridge.center = L0glm_gaus_presetlambda$coefficients, # rep(0, ncol(data_gaus$X)),
                         nonnegative = FALSE,
                         warmupiters = 1,
                         start = start, # rep(1, ncol(data_gaus$X)), # beta_init_orthogonal, # rep(1, ncol(data_gaus$X)), # beta_init_orthogonal, # start # rep(1, ncol(data_gaus$X)), # 
                         solver = "eigen", # eigen, glmnet or osqp, # note: constrained solvers glmnet & osqp better than eigen for very small lambda values
                         method = 3L, 
                         # gamma = 2, 
                         thresh = 1E-10,
                         relaxed = FALSE 
  )$coefficients 
  coefmat[i,] <<- as.vector(L0glm.out)
  return( unlist(c(lambda=lambda, get.metrics(beta = L0glm.out, 
                                              X = X, y = y, 
                                              family = gaussian(identity), 
                                              beta_true = data_gaus$beta_true))) )
} )))) # 26s
matplot(x=log10(metrics$lambda), 
        y=as.matrix(coefmat), type="l", lty=1,
        xlab="Log10(lambda)", ylab="Coefficient estimates", main="Regularisation path L0 Gaussian GLM")
library(scales)
abline(h=data_gaus$beta_true, col=alpha(1:6,0.2))
abline(v=log10(maxlambda), col="red") # max lambda where all coefs would be penalized to 0

plot(log10(metrics$lambda), metrics$deviance, type="l")
metrics$lambda[which.min(metrics$deviance)] # 0.001 # NOTE this is deviance relative to observed signal, not deviance relative to ground truth though
plot(log10(metrics$lambda), metrics$loglik, type="l")
metrics$lambda[which.max(metrics$loglik)] # 0.3981072
plot(log10(metrics$lambda), metrics$MD_coefs, type="l")
plot(log10(metrics$lambda), metrics$MAD_coefs, type="l")
metrics$lambda[which.min(metrics$MAD_coefs)] # 2.511886
plot(log10(metrics$lambda), metrics$sMAPE_coefs, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(metrics$sMAPE_coefs)] # 2.511886
plot(log10(metrics$lambda), metrics$bias_coefs, type="l")
metrics$lambda[which.min(abs(metrics$bias_coefs-1))] # 0.001
plot(log10(metrics$lambda), metrics$R2_coefs, type="l")
plot(log10(metrics$lambda), metrics$estimated_over_true_support_size, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(abs(metrics$estimated_over_true_support_size-1))] # 1.584893
plot(log10(metrics$lambda), metrics$TP, type="l")
plot(log10(metrics$lambda), metrics$FP, type="l")
plot(log10(metrics$lambda), metrics$TN, type="l")
plot(log10(metrics$lambda), metrics$FN, type="l")
plot(log10(metrics$lambda), metrics$precision, type="l")
plot(log10(metrics$lambda), metrics$recall, type="l")
plot(log10(metrics$lambda), metrics$specificity, type="l")
plot(log10(metrics$lambda), metrics$accuracy, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.max(metrics$accuracy)] # 2.511886
plot(log10(metrics$lambda), metrics$F1_score, type="l")
metrics$lambda[which.max(metrics$F1_score)] # 2.511886
plot(log10(metrics$lambda), metrics$aic, type="l")
metrics$lambda[which.min(metrics$aic)] # 1
plot(log10(metrics$lambda), metrics$bic, type="l")
metrics$lambda[which.min(metrics$bic)] # 2.511886
plot(log10(metrics$lambda), metrics$ebic, type="l")
metrics$lambda[which.min(metrics$ebic)] # 2.511886
plot(log10(metrics$lambda), metrics$gic, type="l")
metrics$lambda[which.min(metrics$gic)] # 2.511886
plot(log10(metrics$lambda), metrics$mbic, type="l")
metrics$lambda[which.min(metrics$mbic)] # 3.981072

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
              dispersion = dispersion) # 3.91



# 2. POISSON LOG LINK GLM ####

set.seed(1)
data_pois  = simulate_glm(n, p, n_nonzero, family = "poisson") 
X = data_pois$X
y = data_pois$y

# glm fit
system.time(glm_pois <- glm(data_pois$y ~ 0+data_pois$X, family=poisson(log))) # 0.37s
glm_pois$deviance 
glm_pois$aic
AIC(glm_pois)
BIC(glm_pois)
logLik(glm_pois)

# L0 glm fit with lambda close to zero
system.time(L0glm_pois_lambda0 <- L0glm.fit(X = X, 
                                            y = y, 
                                            family = poisson(log), 
                                            lambda = 1E-10, 
                                            nonzero.centered = FALSE, # use nonzero centered ridge regression?
                                            ridge.center = rep(0, ncol(X)),
                                            nonnegative = FALSE,
                                            warmupiters = 1,
                                            start = rep(1, ncol(data_pois$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                                            solver = "eigen", # 
                                            method = 3L, #
                                            thresh = 1E-10
))  
# TODO: gives bad solution - check what is going wrong...

max(abs(L0glm_pois_lambda0$coefficients-coef(glm_pois))) # TODO gives bad solution - fix
plot(L0glm_pois_lambda0$coefficients, coef(glm_pois), pch=16, col="blue")
# check if metrics (loglik & AIC & BIC) are correct
metrics_pois = data.frame(do.call(cbind, get.metrics(beta = L0glm_pois_lambda0$coefficients,  # some metrics
                                                     X = Matrix(data_pois$X, sparse=TRUE), 
                                                     y = as.matrix(data_pois$y), 
                                                     family = gaussian(identity), 
                                                     beta_true = as.matrix(data_pois$beta_true))))
metrics_pois
logLik(glm_pois)-metrics_pois$loglik # 
AIC(glm_pois)-metrics_pois$aic # 
BIC(glm_pois)-metrics_pois$bic # 


# L0 glm fit with lambda preset to optimise given IC
lambda_preset_pois = preset.lambda(IC = crit, # options aic, bic, gic, ebic, mbic 
                                   y = data_pois$y, 
                                   X = Matrix(data_pois$X, sparse=TRUE), 
                                   family = gaussian(identity),
                                   dispersion = dispersion) # dispersion here was 1, here we pass this value to not have to estimate it 
lambda_preset_pois # 2.556246

system.time(L0glm_pois_presetlambda <- L0glm.fit(X = Matrix(data_pois$X, sparse=TRUE), 
                                                 y = data_pois$y, 
                                                 family = poisson(log), 
                                                 lambda = lambda_preset_pois,
                                                 nonzero.centered = FALSE, # use nonzero centered ridge regression?
                                                 nonnegative = FALSE,
                                                 warmupiters = 1,
                                                 # start = rep(1, ncol(data_pois$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                                                 solver = "eigen", # 
                                                 method = 3L, #
                                                 thresh = 1E-10
)) # 
max(abs(L0glm_pois_presetlambda$coefficients-coef(glm_pois))) # TODO gives bad solution - fix
plot(L0glm_pois_presetlambda$coefficients, coef(glm_pois), pch=16, col="blue")



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
                         nonzero.centered = TRUE, # use nonzero centered ridge regression? 
                         nonnegative = FALSE,
                         warmupiters = 1,
                         # start = start,
                         solver = "eigen", # eigen, glmnet or osqp, # note: constrained solvers glmnet & osqp better than eigen for very small lambda values
                         method = 3L, 
                         # gamma = 2, 
                         thresh = 1E-10,
                         relaxed = FALSE 
  )$coefficients 
  coefmat[i,] <<- as.vector(L0glm.out)
  return( unlist(c(lambda=lambda, get.metrics(beta = L0glm.out, 
                                              X = X, y = y, 
                                              family = gaussian(identity), 
                                              beta_true = data_pois$beta_true))) )
} )))) # 
matplot(x=log10(metrics$lambda), 
        y=as.matrix(coefmat), type="l", lty=1,
        xlab="Log10(lambda)", ylab="Coefficient estimates", main="Regularisation path L0 Gaussian GLM")
abline(v=log10(maxlambda), col="red") # max lambda where all coefs would be penalized to 0

plot(log10(metrics$lambda), metrics$deviance, type="l")
metrics$lambda[which.min(metrics$deviance)] # XXX # NOTE this is deviance relative to observed signal, not deviance relative to ground truth though
plot(log10(metrics$lambda), metrics$loglik, type="l")
metrics$lambda[which.max(metrics$loglik)] # XXX
plot(log10(metrics$lambda), metrics$MD_coefs, type="l")
plot(log10(metrics$lambda), metrics$MAD_coefs, type="l")
metrics$lambda[which.min(metrics$MAD_coefs)] # XXX
plot(log10(metrics$lambda), metrics$sMAPE_coefs, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(metrics$sMAPE_coefs)] # XXX
plot(log10(metrics$lambda), metrics$bias_coefs, type="l")
metrics$lambda[which.min(abs(metrics$bias_coefs-1))] # XXX
plot(log10(metrics$lambda), metrics$R2_coefs, type="l")
plot(log10(metrics$lambda), metrics$estimated_over_true_support_size, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(abs(metrics$estimated_over_true_support_size-1))] # XXX
plot(log10(metrics$lambda), metrics$TP, type="l")
plot(log10(metrics$lambda), metrics$FP, type="l")
plot(log10(metrics$lambda), metrics$TN, type="l")
plot(log10(metrics$lambda), metrics$FN, type="l")
plot(log10(metrics$lambda), metrics$precision, type="l")
plot(log10(metrics$lambda), metrics$recall, type="l")
plot(log10(metrics$lambda), metrics$specificity, type="l")
plot(log10(metrics$lambda), metrics$accuracy, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.max(metrics$accuracy)] # XXX
plot(log10(metrics$lambda), metrics$F1_score, type="l")
metrics$lambda[which.max(metrics$F1_score)] # XXX
plot(log10(metrics$lambda), metrics$aic, type="l")
metrics$lambda[which.min(metrics$aic)] # XXX
plot(log10(metrics$lambda), metrics$bic, type="l")
metrics$lambda[which.min(metrics$bic)] # XXX
plot(log10(metrics$lambda), metrics$gic, type="l")
metrics$lambda[which.min(metrics$gic)] # XXX
plot(log10(metrics$lambda), metrics$ebic, type="l")
metrics$lambda[which.min(metrics$ebic)] # XXX
plot(log10(metrics$lambda), metrics$mbic, type="l")
metrics$lambda[which.min(metrics$mbic)] # XXX

# reasonable match with lambda values that were expected to optimise each IC :
# note: using warm starts in the L0glm.fit function would shift the
# empirically observed optimal values

preset.lambda(IC = "aic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_gaus$y), 
              X = Matrix(data_gaus$X, sparse=TRUE), 
              family = gaussian(identity),
              dispersion = dispersion) # XXX
preset.lambda(IC = "bic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_gaus$y), 
              X = Matrix(data_gaus$X, sparse=TRUE), 
              family = gaussian(identity),
              dispersion = dispersion) # XXX
preset.lambda(IC = "gic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_gaus$y), 
              X = Matrix(data_gaus$X, sparse=TRUE), 
              family = gaussian(identity),
              dispersion = dispersion) # XXX
preset.lambda(IC = "ebic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_gaus$y), 
              X = Matrix(data_gaus$X, sparse=TRUE), 
              family = gaussian(identity),
              dispersion = dispersion) # XXX
preset.lambda(IC = "mbic", # options aic, bic, gic, ebic, mbic 
              y = as.matrix(data_gaus$y), 
              X = Matrix(data_gaus$X, sparse=TRUE), 
              family = gaussian(identity),
              dispersion = dispersion) # XXX