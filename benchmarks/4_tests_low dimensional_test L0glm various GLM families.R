# Check of metrics (AIC, BIC, loglik) L0glm package & comparison with regular GLM ####

# this was an attempt to test L0glm for some different noise families
# but script not up to date yet
# maybe benchmarks in abess and abess:generate.data might be good for
# benchmarking different families too (maybe with some modification, 
# e.g. to also allow simulated sparse covariate matrices)

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
intercept = 1
dispersion = 1
ntot = 10 # nr of trials per row in dataset for binomial case
crit = "gic"

# 1. GAUSSIAN GLM ####

set.seed(1)
data_gaus  = simulate_glm(n = n, p = p, n_nonzero = n_nonzero, intercept = intercept, family = "gaussian") 
X = normalizeX(data_gaus$X) # we first L2 norm normalize the columns of X
y = data_gaus$y
scaling = attr(X, "scaling")
beta_true = data_gaus$beta_true * scaling
family = gaussian(identity)


# glm fit
system.time(glm_gaus <- glm(y ~ 0+X, family=gaussian(identity))) # 0.16s
system.time(glm_gaus2 <- glm.fit(y=y, x=X, intercept=FALSE, family=gaussian(identity))) # 0.14s
glm_gaus$deviance 
glm_gaus$aic
AIC(glm_gaus)
BIC(glm_gaus)
logLik(glm_gaus)
cor(coef(glm_gaus), beta_true) # 0.998713
plot(x=coef(glm_gaus), y=beta_true, pch=16, col="blue")

# L0 glm fit with lambda close to zero
system.time(L0glm_gaus_lambda0 <- L0glm.fit(X = X, 
                                            y = y, 
                                            family = gaussian(identity), 
                                            lambda = 1E-30,
                                            nonnegative = FALSE,
                                            normalize = FALSE,
                                            no.pen = 1,
                                            start = NULL, # rep(sqrt(5), ncol(data_gaus$X)), # NULL, # 
                                            solver = "eigen", # 
                                            method = 4L,
                                            relaxed = FALSE,
                                            thresh = 1E-5
))  # historical timings: 0.06s with LSCG, 0.14s with Cholesky LLT

# NOTE: note that L0glm can fit GLMs with perfectly collinear variables, e.g.
# with X = cbind(data_gaus$X, data_gaus$X)
# in that case, it would divide the total effect over the different collinear versions
# of each variable

max(abs(L0glm_gaus_lambda0$coefficients-coef(glm_gaus))) # 4.780842e-09, correct
cor(L0glm_gaus_lambda0$coefficients, beta_true) # 0.998713
mean(abs(L0glm_gaus_lambda0$coefficients-beta_true)) # 1.648505
plot(x=L0glm_gaus_lambda0$coefficients, y=beta_true, pch=16, col="blue")


# check if loglik & AIC & BIC are correct
system.time(metrics_gaus <- data.frame(do.call(cbind, get.metrics(beta = L0glm_gaus_lambda0$coefficients,  # some metrics
                                                     X = X, 
                                                     y = y,
                                                     family = gaussian(identity), 
                                                     beta_true = as.matrix(data_gaus$beta_true))))) # 0.02s
glm_gaus$deviance-metrics_gaus$deviance # 0, correct
logLik(glm_gaus)-metrics_gaus$loglik # 0, correct
AIC(glm_gaus)-metrics_gaus$aic # 0, correct
BIC(glm_gaus)-metrics_gaus$bic # 0, correct


# L0 glm fit with lambda preset to optimise given IC
lambda_preset_gaus = preset.lambda(IC = crit, # options aic, bic, gic, mebic, ebic, mbic
                              y = y,
                              X = X, 
                              family = gaussian(identity),
                              dispersion = dispersion # dispersion here was 1, here we pass this value to not have to estimate it 
                              ) 
lambda_preset_gaus # 2.561769

system.time(L0glm_gaus_presetlambda <- L0glm.fit(X = X, 
                                            y = y, 
                                            family = gaussian(identity), 
                                            lambda = lambda_preset_gaus,
                                            no.pen = 1,
                                            nonnegative = FALSE,
                                            normalize = FALSE,
                                            start = NULL, # =rep(sqrt(5), ncol(data_gaus$X))/scaling, # NULL, # rep(sqrt(5), ncol(data_gaus$X))*scaling, # NULL, # start_coefs(X, family=gaussian(identity)), # = rep(sqrt(5), ncol(data_gaus$X)) / scaling, 
                                            solver = "eigen", # 
                                            method = 1L, #
                                            relaxed = FALSE,
                                            thresh = 1E-5
)) # 0.41s with method 3L, 0.72s with method 7L
# TODO: check solver method 1 here fails
sum(L0glm_gaus_presetlambda$coefficients!=0)  # 47
cor(L0glm_gaus_presetlambda$coefficients, data_gaus$beta_true) # 0.9940267
plot(x=L0glm_gaus_presetlambda$coefficients, y=data_gaus$beta_true, pch=16, col="blue")
points(x=L0glm_gaus_presetlambda$coefficients[data_gaus$beta_true==0], y=data_gaus$beta_true[data_gaus$beta_true==0], 
     pch=16, col="red")

# check metrics & quality of solution with preset lambda
system.time(metrics_gaus <- data.frame(do.call(cbind, get.metrics(beta = L0glm_gaus_presetlambda$coefficients,  # some metrics
                                                                  X = X, 
                                                                  y = y,
                                                                  family = gaussian(identity), 
                                                                  beta_true = beta_true)))) 
metrics_gaus
#     p     n m deviance deviance_true    loglik    MD_coefs MAD_coefs sMAPE_coefs bias_coefs  R2_coefs
# 1 102 10000 1 10005.95      80.24103 -14138.36 0.009655372 0.9534635   0.1049711   1.005573 0.9982282
#   support_size_true support_size_estimated estimated_over_true_support_size TP FP TN FN precision    recall
# 1                51                     47                        0.9215686 47  0 50  4         1 0.9215686
#   specificity accuracy  F1_score      gic     ebic      bic      aic     mbic
# 1           1 0.960396 0.9591837 28769.63 28720.71 28718.81 28372.72 29029.73


# orthogonal coefficient initialisation for possible use as start ####
beta_init_orthogonal <- init_coefs(X = X,
           mu = y,
           family = gaussian(identity),
           lambda = lambda_preset_gaus)
cor(beta_init_orthogonal, data_gaus$beta_true) # 0.981208
plot(beta_init_orthogonal, data_gaus$beta_true, pch=16, col="blue")
mean(abs(beta_init_orthogonal-data_gaus$beta_true))

# coefficients & metrics over a range of lambda values
# calculate max lambda where all coefs would be penalized to 0
maxlam <- maxlambda(X=X, 
                    weights=L0glm_gaus_presetlambda$weights, 
                    z=L0glm_gaus_presetlambda$z,
                    no.pen=1,
                    coefs=L0glm_gaus_presetlambda$coefficients,
                    family=family)
maxlam # 2346.189
mincoef <- min(abs(L0glm_gaus_presetlambda$coefficients))

lambdas <- 10^seq(round(log10(lambda_preset_gaus))-1, log10(maxlam), by=0.1)
coefmat <- Matrix(0, nrow=length(lambdas), ncol=ncol(X))
warmstart <- 0.2 # % by which L0glm_gaus_presetlambda$coefficients will be used as start
system.time(metrics <- data.frame(do.call(rbind, lapply(1:length(lambdas), function (i) { 
  lambda <- lambdas[[i]]
  start <- (L0glm_gaus_presetlambda$coefficients)*warmstart+start_coefs(X, family)*(1-warmstart)
  # start[L0glm_gaus_presetlambda$coefficients==0] <- 0
  
  L0glm.out <- L0glm.fit(X = X, 
                         y = y, 
                         family = gaussian(identity), 
                         lambda = lambda, 
                         no.pen = 1,
                         nonnegative = FALSE,
                         start = start, # 
                         solver = "eigen", # eigen, glmnet or osqp, # note: constrained solvers glmnet & osqp better than eigen for very small lambda values
                         method = 4L, # historically LSCG was twice as fast as Cholesky LLT here
                         # gamma = 2, 
                         thresh = 1E-5,
                         relaxed = FALSE 
  )$coefficients 
  coefmat[i,] <<- as.vector(L0glm.out)
  return( unlist(c(lambda=lambda, get.metrics(beta = L0glm.out, 
                                              X = X, y = y, 
                                              family = gaussian(identity), 
                                              beta_true = beta_true))) )
} )))) # 19s
matplot(x=log10(metrics$lambda), 
        y=as.matrix(coefmat), type="l", lty=1,
        xlab="Log10(lambda)", ylab="Coefficient estimates", main="Regularisation path L0 Gaussian GLM")
library(scales)
abline(h=beta_true, col=alpha(1:6,0.2))
abline(v=log10(maxlam), col="red") # max lambda where all coefs would be penalized to 0

plot(log10(metrics$lambda), metrics$deviance, type="l")
metrics$lambda[which.min(metrics$deviance)] # 0.1 # NOTE this is deviance relative to observed signal, not deviance relative to ground truth though
plot(log10(metrics$lambda), metrics$deviance_true, type="l")
metrics$lambda[which.min(metrics$deviance_true)] # 0.7943282 # deviance relative to ground truth
plot(log10(metrics$lambda), metrics$loglik, type="l")
metrics$lambda[which.max(metrics$loglik)] # 0.3981072
plot(log10(metrics$lambda), metrics$MD_coefs, type="l")
plot(log10(metrics$lambda), metrics$MAD_coefs, type="l")
metrics$lambda[which.min(metrics$MAD_coefs)] # 1
plot(log10(metrics$lambda), metrics$sMAPE_coefs, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(metrics$sMAPE_coefs)] # 1.995262
plot(log10(metrics$lambda), metrics$bias_coefs, type="l")
metrics$lambda[which(abs(metrics$bias_coefs-1)==min(abs(metrics$bias_coefs-1)))] # 794.3282
plot(log10(metrics$lambda), metrics$R2_coefs, type="l")
metrics$lambda[which.max(abs(metrics$R2_coefs))] # 0.7943282
plot(log10(metrics$lambda), metrics$estimated_over_true_support_size, type="l") # good one, optimal for lambda between 1 & 10
range(metrics$lambda[which(abs(metrics$estimated_over_true_support_size-1)==min(abs(metrics$estimated_over_true_support_size-1)))]) # 1.258925 1.258925
plot(log10(metrics$lambda), metrics$TP, type="l")
plot(log10(metrics$lambda), metrics$FP, type="l")
plot(log10(metrics$lambda), metrics$TN, type="l")
plot(log10(metrics$lambda), metrics$FN, type="l")
plot(log10(metrics$lambda), metrics$precision, type="l")
plot(log10(metrics$lambda), metrics$recall, type="l")
plot(log10(metrics$lambda), metrics$specificity, type="l")
plot(log10(metrics$lambda), metrics$accuracy, type="l") # good one, optimal for lambda between 1 & 10
range(metrics$lambda[which(metrics$accuracy==max(metrics$accuracy))]) # 1.995262 1.995262
plot(log10(metrics$lambda), metrics$F1_score, type="l")
range(metrics$lambda[which(metrics$F1_score==max(metrics$F1_score))]) # 1.995262 1.995262
plot(log10(metrics$lambda), metrics$aic, type="l")
metrics$lambda[which(metrics$aic==min(metrics$aic))] # 1
plot(log10(metrics$lambda), metrics$bic, type="l")
metrics$lambda[which.min(metrics$bic)] # 1
plot(log10(metrics$lambda), metrics$ebic, type="l")
metrics$lambda[which.min(metrics$ebic)] # 1
plot(log10(metrics$lambda), metrics$gic, type="l")
metrics$lambda[which.min(metrics$gic)] # 1
plot(log10(metrics$lambda), metrics$mbic, type="l")
metrics$lambda[which(metrics$mbic==min(metrics$mbic))] # 1.995262

# with warmstart=0.2 ca 2x lower than the lambda tilde values that 
# were expected to optimise each IC; close match with expected values
# for warmstart=0
# note: as expected using warm starts in the L0glm.fit function shift the
# empirically observed optimal values

preset.lambda(IC = "aic", # options aic, bic, gic, mebic, ebic, mbic
              y = as.matrix(data_gaus$y), 
              X = Matrix(data_gaus$X, sparse=TRUE), 
              family = gaussian(identity),
              dispersion = dispersion) # 0.5
preset.lambda(IC = "bic", # options aic, bic, gic, mebic, ebic, mbic
              y = as.matrix(data_gaus$y), 
              X = Matrix(data_gaus$X, sparse=TRUE), 
              family = gaussian(identity),
              dispersion = dispersion) # 2.30
preset.lambda(IC = "ebic", # options aic, bic, gic, mebic, ebic, mbic
              y = as.matrix(data_gaus$y), 
              X = Matrix(data_gaus$X, sparse=TRUE), 
              family = gaussian(identity),
              dispersion = dispersion) # 2.30
preset.lambda(IC = "gic", # options aic, bic, gic, mebic, ebic, mbic
              y = as.matrix(data_gaus$y), 
              X = Matrix(data_gaus$X, sparse=TRUE), 
              family = gaussian(identity),
              dispersion = dispersion) # 2.56
preset.lambda(IC = "mbic", # options aic, bic, gic, mebic, ebic, mbic
              y = as.matrix(data_gaus$y), 
              X = Matrix(data_gaus$X, sparse=TRUE), 
              family = gaussian(identity),
              dispersion = dispersion) # 3.92



# 2. LOG LINK POISSON GLM ####

set.seed(1)
data_pois  = simulate_glm(n = n, p = p, n_nonzero = n_nonzero, intercept = intercept, family = "poisson") 
X = normalizeX(data_pois$X) # we first L2 norm normalize the columns of X
y = data_pois$y
scaling = attr(X, "scaling")
beta_true = data_pois$beta_true * scaling
family = poisson(log)

# glm fit
system.time(glm_pois <- glm(y ~ 0+X, family=poisson(log))) # 0.39s
glm_pois$deviance 
glm_pois$aic
AIC(glm_pois)
BIC(glm_pois)
logLik(glm_pois)
cor(coef(glm_pois), beta_true) # 0.9997832
plot(x=coef(glm_pois), y=beta_true, pch=16, col="blue")

# L0 glm fit with lambda close to zero
system.time(L0glm_pois_lambda0 <- L0glm.fit(X = X, 
                                            y = y, 
                                            family = poisson(log), 
                                            lambda = 1E-30, 
                                            no.pen = 1,
                                            nonnegative = FALSE,
                                            normalize = FALSE,
                                            start = NULL, # = start_coefs(X=X, family=family) = rep(sqrt(5), ncol(X))/scaling,
                                            solver = "eigen", # 
                                            method = 1L, #
                                            thresh = 1E-5
))   # 0.99s # note 3L is 2x faster but a bit more inaccurate
L0glm_pois_lambda0_coefs = L0glm_pois_lambda0$coefficients
max(abs(L0glm_pois_lambda0_coefs-coef(glm_pois))) # 1.109015e-10
plot(L0glm_pois_lambda0$coefficients, coef(glm_pois), pch=16, col="blue")

cor(L0glm_pois_lambda0$coefficients, beta_true) # 0.9997832
plot(x=L0glm_pois_lambda0$coefficients, y=beta_true, pch=16, col="blue")

# check if loglik & AIC & BIC are correct
metrics_pois = data.frame(do.call(cbind, get.metrics(beta = L0glm_pois_lambda0$coefficients,  # some metrics
                                                     X = X, 
                                                     y = y, 
                                                     family = poisson(log), 
                                                     beta_true = beta_true)))
metrics_pois
# TODO: resolve warnings In dpois(y, mu, log = TRUE) : non-integer x
glm_pois$deviance-metrics_pois$deviance # -1.818989e-12, correct
logLik(glm_pois)-metrics_pois$loglik # 0, correct
AIC(glm_pois)-metrics_pois$aic # 0, correct
BIC(glm_pois)-metrics_pois$bic # 0, correct


# L0 glm fit with lambda preset to optimise given IC
lambda_preset_pois = preset.lambda(IC = crit, # options aic, bic, gic, mebic, ebic, mbic
                                   y = y, 
                                   X = X, 
                                   family = poisson(log),
                                   dispersion = dispersion) # dispersion here was 1, here we pass this value to not have to estimate it 
lambda_preset_pois # 2.561769

system.time(L0glm_pois_presetlambda <- L0glm.fit(X = X, 
                                                 y = y, 
                                                 family = poisson(log), 
                                                 lambda = lambda_preset_pois,
                                                 no.pen = 1,
                                                 nonnegative = FALSE,
                                                 normalize = FALSE,
                                                 start = NULL, # =start_coefs(X=X, family=family) = rep(sqrt(5), ncol(X))/scaling
                                                 solver = "eigen", # 
                                                 method = 4L, #
                                                 thresh = 1E-5
)) # 0.48s
mean(abs(L0glm_pois_presetlambda$coefficients-beta_true)) # 0.3255121
mean(abs(coef(glm_pois)-beta_true)) # 0.5939227
plot(L0glm_pois_presetlambda$coefficients, coef(glm_pois), pch=16, col="blue")

cor(L0glm_pois_presetlambda$coefficients, beta_true) # 0.9998809
plot(x=L0glm_pois_presetlambda$coefficients, y=beta_true, pch=16, col="blue")
points(x=L0glm_pois_presetlambda$coefficients[beta_true==0], y=beta_true[beta_true==0], pch=16, col="red")

# check metrics & quality of solution with preset lambda
system.time(metrics_pois <- data.frame(do.call(cbind, get.metrics(beta = L0glm_pois_presetlambda$coefficients,  # some metrics
                                                                  X = X, 
                                                                  y = y,
                                                                  family = poisson(log), 
                                                                  beta_true = beta_true)))) 
metrics_pois
#     p     n m deviance deviance_true    loglik    MD_coefs MAD_coefs sMAPE_coefs bias_coefs  R2_coefs
# 1 101 10000 1 8794.341      51.81464 -13666.31 0.003695359 0.3255121  0.05111916   1.000102 0.9997617
# support_size_true support_size_estimated estimated_over_true_support_size TP FP TN FN precision    recall
# 1                51                     49                        0.9607843 49  0 50  2         1 0.9607843
#   specificity accuracy F1_score      gic    ebic      bic      aic     mbic
# 1           1 0.980198     0.98 27834.73 27784.9 27783.93 27430.62 28100.35
# TODO fix warnings e.g. In dpois(y, mu, log = TRUE) : non-integer x = 0.322116
# add a round() somewhere?


# coefficients & metrics over a range of lambda values
# calculate max lambda where all coefs would be penalized to 0
maxlam <- maxlambda(X=X, weights=L0glm_pois_presetlambda$weights, 
                    z=L0glm_pois_presetlambda$z,
                    no.pen=1, coefs=L0glm_pois_presetlambda$coefficients,
                    family=family)
maxlam # 57297.43
mincoef <- min(abs(L0glm_pois_presetlambda$coefficients))

lambdas <- 10^seq(round(log10(lambda_preset_pois))-1, log10(maxlam), by=0.1)
coefmat <- Matrix(0, nrow=length(lambdas), ncol=ncol(X))
warmstart <- 0.2 # here using a higher value than for Gaussian case
system.time(metrics <- data.frame(do.call(rbind, lapply(1:length(lambdas), function (i) { 
  # print(i)
  lambda <- lambdas[[i]]
  start <- L0glm_pois_presetlambda$coefficients*warmstart+start_coefs(X, family)*(1-warmstart)
  # start[L0glm_pois_presetlambda$coefficients==0] <- 0

  L0glm.out <- L0glm.fit(X = Matrix(X, sparse=T), 
                         y = y, 
                         family = poisson(log), 
                         lambda = lambda, 
                         no.pen = 1,
                         nonnegative = FALSE,
                         normalize = FALSE,
                         start = start, 
                         solver = "eigen", # eigen, glmnet or osqp, # note: constrained solvers glmnet & osqp better than eigen for very small lambda values
                         method = 4L,
                         # gamma = 2, 
                         thresh = 1E-5,
                         relaxed = FALSE 
  )$coefficients 
  coefmat[i,] <<- as.vector(L0glm.out)
  return( unlist(c(lambda=lambda, get.metrics(beta = L0glm.out, 
                                              X = X, y = y, 
                                              family = poisson(log), 
                                              beta_true = beta_true))) )
} )))) # 56s
# TODO with method = 1L and dense X = X returns Error in h(simpleError(msg, call)) :
# error in evaluating the argument 'x' in selecting a method for function 'as.matrix': non-conformable arguments
# error in solve_dense_chol
# Add some error catching in solvers?
matplot(x=log10(metrics$lambda), 
        y=as.matrix(coefmat), type="l", lty=1,
        xlab="Log10(lambda)", ylab="Coefficient estimates", main="Regularisation path L0 Poisson GLM")
abline(v=log10(maxlam), col="red") 
# max lambda where all coefs would be penalized to 0
library(scales)
abline(h=beta_true, col=alpha(1:6,0.2))

plot(log10(metrics$lambda), metrics$deviance, type="l")
metrics$lambda[which.min(metrics$deviance)] # 0.1 # NOTE this is deviance relative to observed signal, not deviance relative to ground truth though
plot(log10(metrics$lambda), metrics$deviance_true, type="l")
metrics$lambda[which.min(metrics$deviance_true)] # 1 # deviance relative to ground truth
plot(log10(metrics$lambda), metrics$loglik, type="l")
metrics$lambda[which.max(metrics$loglik)] # 0.5011872
plot(log10(metrics$lambda), metrics$MD_coefs, type="l")
plot(log10(metrics$lambda), metrics$MAD_coefs, type="l")
metrics$lambda[which.min(metrics$MAD_coefs)] # 1
plot(log10(metrics$lambda), metrics$sMAPE_coefs, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(metrics$sMAPE_coefs)] # 1
plot(log10(metrics$lambda), metrics$bias_coefs, type="l")
metrics$lambda[which.min(abs(metrics$bias_coefs-1))] # 0.2511886
plot(log10(metrics$lambda), metrics$R2_coefs, type="l")
plot(log10(metrics$lambda), metrics$estimated_over_true_support_size, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(abs(metrics$estimated_over_true_support_size-1))] # 0.5011872
plot(log10(metrics$lambda), metrics$TP, type="l")
plot(log10(metrics$lambda), metrics$FP, type="l")
plot(log10(metrics$lambda), metrics$TN, type="l")
plot(log10(metrics$lambda), metrics$FN, type="l")
plot(log10(metrics$lambda), metrics$precision, type="l")
plot(log10(metrics$lambda), metrics$recall, type="l")
plot(log10(metrics$lambda), metrics$specificity, type="l")
plot(log10(metrics$lambda), metrics$accuracy, type="l") # good one, optimal for lambda between 1 & 10
range(metrics$lambda[which(metrics$accuracy==max(metrics$accuracy))]) # 1
plot(log10(metrics$lambda), metrics$F1_score, type="l")
range(metrics$lambda[which(metrics$F1_score==max(metrics$F1_score))]) # 1
plot(log10(metrics$lambda), metrics$aic, type="l")
metrics$lambda[which(metrics$aic==min(metrics$aic))] # 1
plot(log10(metrics$lambda), metrics$bic, type="l")
metrics$lambda[which(metrics$bic==min(metrics$bic))] # 1.258925
plot(log10(metrics$lambda), metrics$gic, type="l")
metrics$lambda[which(metrics$gic==min(metrics$gic))] # 1.258925
plot(log10(metrics$lambda), metrics$ebic, type="l")
metrics$lambda[which.min(metrics$ebic)] # 1.258925
plot(log10(metrics$lambda), metrics$mbic, type="l")
metrics$lambda[which(metrics$mbic==min(metrics$mbic))] # 1.258925

# with warmstart=0 or 0.2 ca 2x lower than the lambda tilde values that 
# were expected to optimise each IC
# note: as expected using warm starts in the L0glm.fit function shift the
# empirically observed optimal values

preset.lambda(IC = "aic", # options aic, bic, gic, mebic, ebic, mbic
              y = as.matrix(data_pois$y), 
              X = Matrix(data_pois$X, sparse=TRUE), 
              family = poisson(log),
              dispersion = dispersion) # 0.5
preset.lambda(IC = "bic", # options aic, bic, gic, mebic, ebic, mbic
              y = as.matrix(data_pois$y), 
              X = Matrix(data_pois$X, sparse=TRUE), 
              family = poisson(log),
              dispersion = dispersion) # 2.30
preset.lambda(IC = "gic", # options aic, bic, gic, mebic, ebic, mbic
              y = as.matrix(data_pois$y), 
              X = Matrix(data_pois$X, sparse=TRUE), 
              family = poisson(log),
              dispersion = dispersion) # 2.56
preset.lambda(IC = "ebic", # options aic, bic, gic, mebic, ebic, mbic
              y = as.matrix(data_pois$y), 
              X = Matrix(data_pois$X, sparse=TRUE), 
              family = poisson(log),
              dispersion = dispersion) # 2.30
preset.lambda(IC = "mbic", # options aic, bic, gic, mebic, ebic, mbic
              y = as.matrix(data_pois$y), 
              X = Matrix(data_pois$X, sparse=TRUE), 
              family = poisson(log),
              dispersion = dispersion) # 3.91


# 3. LOGIT LINK BINOMIAL GLM ####

set.seed(1)
data_binom  = simulate_glm(n = n, p = p, n_nonzero = n_nonzero, intercept = intercept, family = "binomial", ntot = ntot)
X = normalizeX(data_binom$X) # we first L2 norm normalize the columns of X
y = data_binom$y
scaling = attr(X, "scaling")
beta_true = data_binom$beta_true * scaling
family = binomial(logit)
weights = data_binom$weights # total sample size per row in dataset

# glm fit
system.time(glm_binom <- glm(y ~ 0+X, weights=data_binom$weights, family=binomial(logit))) # 0.37s
glm_binom$deviance 
glm_binom$aic
AIC(glm_binom)
BIC(glm_binom)
logLik(glm_binom)
cor(coef(glm_binom), beta_true) # 0.9971224
plot(x=coef(glm_binom), y=beta_true, pch=16, col="blue")

# L0 glm fit with lambda close to zero
system.time(L0glm_binom_lambda0 <- L0glm.fit(X = X, 
                                            y = y, 
                                            weights = weights,
                                            family = binomial(logit), 
                                            lambda = 1E-30, 
                                            nonnegative = FALSE,
                                            normalize = FALSE,
                                            start = NULL, # = rep(sqrt(5), ncol(data_binom$X))/scaling,
                                            no.pen = 1,
                                            solver = "eigen", # 
                                            method = 4L, #
                                            thresh = 1E-5))
# 0.61s
mean(abs(coef(glm_binom)-L0glm_binom_lambda0$coefficients)) # 1.271946e-09
plot(x=L0glm_binom_lambda0$coefficients, y=beta_true, pch=16, col="blue")

cor(L0glm_binom_lambda0$coefficients, beta_true) # 0.9971224

# check if loglik & AIC & BIC are correct
metrics_binom = data.frame(do.call(cbind, get.metrics(beta = as.matrix(L0glm_binom_lambda0$coefficients),  # some metrics
                                                     X = X, 
                                                     y = y, 
                                                     weights = weights,
                                                     family = binomial(logit), 
                                                     beta_true = beta_true ) ) )
metrics_binom
glm_binom$deviance-metrics_binom$deviance # 0, correct
logLik(glm_binom)-metrics_binom$loglik # 0, correct
AIC(glm_binom)-metrics_binom$aic # 0, correct
BIC(glm_binom)-metrics_binom$bic # 0, correct


# L0 glm fit with lambda preset to optimise given IC
lambda_preset_binom = preset.lambda(IC = crit, # options aic, bic, gic, mebic, ebic, mbic
                                   y = y, 
                                   X = X, 
                                   family = binomial(logit),
                                   dispersion = dispersion) # dispersion here was 1, here we pass this value to not have to estimate it 
lambda_preset_binom # 2.561769

system.time(L0glm_binom_presetlambda <- L0glm.fit(X = X, 
                                                      y = y, 
                                                      weights = weights,
                                                      family = binomial(logit), 
                                                      lambda = lambda_preset_binom, 
                                                      nonnegative = FALSE,
                                                      normalize = FALSE,
                                                      start = NULL, # = rep(sqrt(5), ncol(X))/scaling
                                                      no.pen = 1,
                                                      solver = "eigen", # 
                                                      method = 1L, #
                                                      thresh = 1E-5)) 
# 0.89s
mean(abs(L0glm_binom_presetlambda$coefficients-coef(glm_binom))) # 0.9392265
plot(L0glm_binom_presetlambda$coefficients, coef(glm_binom), pch=16, col="blue")

cor(L0glm_binom_presetlambda$coefficients, beta_true) # 0.9982796
mean(abs(L0glm_binom_presetlambda$coefficients-beta_true)) # 0.9830373
mean(abs(coef(glm_binom)-beta_true)) # 1.678541
plot(x=L0glm_binom_presetlambda$coefficients, y=beta_true, pch=16, col="blue")
points(x=L0glm_binom_presetlambda$coefficients[beta_true==0], y=beta_true[beta_true==0], pch=16, col="red")


# check metrics & quality of solution with preset lambda
system.time(metrics_binom <- data.frame(do.call(cbind, get.metrics(beta = L0glm_binom_presetlambda$coefficients,  # some metrics
                                                                  X = X, 
                                                                  y = y,
                                                                  family = binomial(logit), 
                                                                  beta_true = beta_true)))) 
metrics_binom
#     p     n m deviance deviance_true    loglik   MD_coefs MAD_coefs sMAPE_coefs bias_coefs
# 1 101 10000 1 1054.091      8.248035 -2952.732 0.03663581 0.9830373  0.08576971   1.012569
#            R2_coefs support_size_true support_size_estimated estimated_over_true_support_size TP FP TN
# 1         0.9965617                51                     48                        0.9411765 48  0 50
# FN precision    recall specificity accuracy F1_score      gic     ebic     bic      aic
# 1          3         1   0.9411765                 1 0.970297 0.969697 6397.324 6348.516 6347.56 6001.464
#       mbic
# 1 6657.528


# coefficients & metrics over a range of lambda values
# calculate max lambda where all coefs would be penalized to 0
maxlam <- maxlambda(X=X, weights=L0glm_binom_presetlambda$weights, 
                    z=L0glm_binom_presetlambda$z, no.pen=1,
                    coefs=L0glm_binom_presetlambda$coefficients,
                    family=family)
maxlam # 1460.132
mincoef <- min(abs(L0glm_binom_presetlambda$coefficients)) # *sign(L0glm_binom_presetlambda$coefficients)

lambdas <- 10^seq(round(log10(lambda_preset_binom))-1, log10(maxlam), by=0.1)
coefmat <- Matrix(0, nrow=length(lambdas), ncol=ncol(X))
warmstart <- 0.2 # % by which L0glm_binom_presetlambda$coefficients will be used as start, lowest value that led to stable coefficient paths
system.time(metrics <- data.frame(do.call(rbind, lapply(1:length(lambdas), function (i) { 
  lambda <- lambdas[[i]]
  start <- L0glm_binom_presetlambda$coefficients*warmstart+start_coefs(X, family)*(1-warmstart) # (rep(sqrt(100), ncol(X))/scaling)
  # start[L0glm_binom_presetlambda$coefficients==0] <- 0
  L0glm.out <- L0glm.fit(X = X, 
                         y = y, 
                         family = binomial(logit), 
                         lambda = lambda, 
                         nonnegative = FALSE,
                         normalize = FALSE,
                         no.pen = 1,
                         start = start, # 
                         solver = "eigen", # eigen, glmnet or osqp, # note: constrained solvers glmnet & osqp better than eigen for very small lambda values
                         method = 4L,
                         # gamma = 2, 
                         thresh = 1E-5,
                         relaxed = FALSE 
  )$coefficients 
  coefmat[i,] <<- as.vector(L0glm.out)
  return( unlist(c(lambda=lambda, get.metrics(beta = L0glm.out, 
                                              X = X, y = y, 
                                              family = binomial(logit), 
                                              beta_true = beta_true))) )
} )))) # 19s
# TODO debug warnings In eval(init) : non-integer #successes in a binomial glm!
# TODO with method=1L and dense X=X returns error Error in h(simpleError(msg, call)) :
# error in evaluating the argument 'x' in selecting a method for function 'as.matrix': non-conformable arguments
# In Cholesky dense solver. Add some error catching in solvers?
matplot(x=log10(metrics$lambda), 
        y=as.matrix(coefmat), type="l", lty=1,
        xlab="Log10(lambda)", ylab="Coefficient estimates", main="Regularisation path L0 Binomial GLM")
abline(v=log10(maxlam), col="red") # max lambda where all coefs would be penalized to 0
library(scales)
abline(h=beta_true, col=alpha(1:6,0.2))

plot(log10(metrics$lambda), metrics$deviance, type="l")
metrics$lambda[which.min(metrics$deviance)] # 0.1 # NOTE this is deviance relative to observed signal, not deviance relative to ground truth though
plot(log10(metrics$lambda), metrics$loglik, type="l")
metrics$lambda[which.max(metrics$loglik)] # 0.1258925
plot(log10(metrics$lambda), metrics$MD_coefs, type="l")
plot(log10(metrics$lambda), metrics$MAD_coefs, type="l")
metrics$lambda[which.min(metrics$MAD_coefs)] # 0.1258925
plot(log10(metrics$lambda), metrics$sMAPE_coefs, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(metrics$sMAPE_coefs)] # 0.1258925
plot(log10(metrics$lambda), metrics$bias_coefs, type="l")
metrics$lambda[which.min(abs(metrics$bias_coefs-1))] # 25.11886
plot(log10(metrics$lambda), metrics$R2_coefs, type="l")
plot(log10(metrics$lambda), metrics$estimated_over_true_support_size, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which(abs(metrics$estimated_over_true_support_size-1)==min(abs(metrics$estimated_over_true_support_size-1)))] # 0.1258925
plot(log10(metrics$lambda), metrics$TP, type="l")
plot(log10(metrics$lambda), metrics$FP, type="l")
plot(log10(metrics$lambda), metrics$TN, type="l")
plot(log10(metrics$lambda), metrics$FN, type="l")
plot(log10(metrics$lambda), metrics$precision, type="l")
plot(log10(metrics$lambda), metrics$recall, type="l")
plot(log10(metrics$lambda), metrics$specificity, type="l")
plot(log10(metrics$lambda), metrics$accuracy, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.max(metrics$accuracy)] # 0.1258925
plot(log10(metrics$lambda), metrics$F1_score, type="l")
metrics$lambda[which.max(metrics$F1_score)] # 0.1258925
plot(log10(metrics$lambda), metrics$aic, type="l")
metrics$lambda[which.min(metrics$aic)] # 0.1258925
plot(log10(metrics$lambda), metrics$bic, type="l")
metrics$lambda[which.min(metrics$bic)] # 0.1584893
plot(log10(metrics$lambda), metrics$gic, type="l")
metrics$lambda[which.min(metrics$gic)] # 0.1584893
plot(log10(metrics$lambda), metrics$ebic, type="l")
metrics$lambda[which.min(metrics$ebic)] # 0.1584893
plot(log10(metrics$lambda), metrics$mbic, type="l")
metrics$lambda[which(metrics$mbic==min(metrics$mbic))] # 0.6309573

# ca 4x lower than lambda values that were expected to optimise each IC :
# note: as expected using warm starts in the L0glm.fit function shifts the
# empirically observed optimal values

preset.lambda(IC = "aic", # options aic, bic, gic, mebic, ebic, mbic
              y = as.matrix(data_binom$y), 
              X = Matrix(data_binom$X, sparse=TRUE), 
              family = binomial(logit),
              dispersion = dispersion) # 0.5
preset.lambda(IC = "bic", # options aic, bic, gic, mebic, ebic, mbic
              y = as.matrix(data_binom$y), 
              X = Matrix(data_binom$X, sparse=TRUE), 
              family = binomial(logit),
              dispersion = dispersion) # 2.30
preset.lambda(IC = "gic", # options aic, bic, gic, mebic, ebic, mbic
              y = as.matrix(data_binom$y), 
              X = Matrix(data_binom$X, sparse=TRUE), 
              family = binomial(logit),
              dispersion = dispersion) # 2.56
preset.lambda(IC = "ebic", # options aic, bic, gic, mebic, ebic, mbic
              y = as.matrix(data_binom$y), 
              X = Matrix(data_binom$X, sparse=TRUE), 
              family = binomial(logit),
              dispersion = dispersion) # 2.30
preset.lambda(IC = "mbic", # options aic, bic, gic, mebic, ebic, mbic
              y = as.matrix(data_binom$y), 
              X = Matrix(data_binom$X, sparse=TRUE), 
              family = binomial(logit),
              dispersion = dispersion) # 3.91




# 4. LOG LINK GAMMA GLM ####

set.seed(1)
data_gamma  = simulate_glm(n, p, n_nonzero, family = "gamma") 
X = normalizeX(data_gamma$X) # we first L2 norm normalize the columns of X
y = data_gamma$y
scaling = attr(X, "scaling")
beta_true = data_gamma$beta_true * scaling
family = Gamma(log)

# glm fit
system.time(glm_gamma <- glm(y ~ 0+X, family=Gamma(log))) # 0.53s
glm_gamma$deviance 
glm_gamma$aic
AIC(glm_gamma)
BIC(glm_gamma)
logLik(glm_gamma)
cor(coef(glm_gamma), beta_true) # 0.9987695
plot(x=coef(glm_gamma), y=beta_true, pch=16, col="blue")


# L0 glm fit with lambda close to zero
system.time(L0glm_gamma_lambda0 <- L0glm.fit(X = X, 
                                            y = y, 
                                            family = Gamma(log), 
                                            lambda = 1E-30, 
                                            nonnegative = FALSE,
                                            normalize = FALSE,
                                            no.pen = 1,
                                            start = NULL, # rep(sqrt(5), ncol(data_gamma$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                                            solver = "eigen", # 
                                            method = 1L, #
                                            thresh = 1E-5
))  # 2.19s

max(abs(L0glm_gamma_lambda0$coefficients-coef(glm_gamma))) # 0.001024
plot(L0glm_gamma_lambda0$coefficients, coef(glm_gamma), pch=16, col="blue")

# check if loglik & AIC & BIC are correct
metrics_gamma = data.frame(do.call(cbind, get.metrics(beta = L0glm_gamma_lambda0$coefficients,  # some metrics
                                                     X = X, 
                                                     y = y, 
                                                     family = Gamma(log), 
                                                     beta_true = beta_true)))
metrics_gamma
logLik(glm_gamma)-metrics_gamma$loglik # -4.90345e-08, correct
AIC(glm_gamma)-metrics_gamma$aic # 9.806899e-08
BIC(glm_gamma)-metrics_gamma$bic # 9.806899e-08


# L0 glm fit with lambda preset to optimise given IC
lambda_preset_gamma = preset.lambda(IC = crit, # options aic, bic, gic, mebic, ebic, mbic
                                   y = y, 
                                   X = X, 
                                   family = Gamma(log),
                                   dispersion = dispersion) # dispersion here was 1, here we pass this value to not have to estimate it 
lambda_preset_gamma # 2.561769

system.time(L0glm_gamma_presetlambda <- L0glm.fit(X = X, 
                                                 y = y, 
                                                 family = Gamma(log), 
                                                 lambda = lambda_preset_gamma,
                                                 nonnegative = FALSE,
                                                 normalize = FALSE,
                                                 no.pen = 1,
                                                 start = NULL, # start_coefs(X,family)/scaleX(X), # rep(sqrt(50), ncol(X)), # rep(sqrt(6), ncol(data_gamma$X)), # rep(sqrt(6), ncol(data_gamma$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                                                 solver = "eigen", # 
                                                 method = 1L, #
                                                 thresh = 1E-5
)) # 1.33s
mean(abs(L0glm_gamma_presetlambda$coefficients-coef(glm_gamma))) # 0.9436883
mean(abs(L0glm_gamma_presetlambda$coefficients-beta_true)) # 0.8023027
mean(abs(coef(glm_gamma)-beta_true)) # 1.616754

plot(x=L0glm_gamma_presetlambda$coefficients, y=beta_true, pch=16, col="blue")
points(x=L0glm_gamma_presetlambda$coefficients[beta_true==0], y=beta_true[beta_true==0], pch=16, col="red")


# check metrics & quality of solution with preset lambda
system.time(metrics_gamma <- data.frame(do.call(cbind, get.metrics(beta = L0glm_gamma_presetlambda$coefficients,  # some metrics
                                                                   X = X, 
                                                                   y = y,
                                                                   family = Gamma(log), 
                                                                   beta_true = beta_true)))) 
metrics_gamma
#     p     n m deviance deviance_true    loglik      MD_coefs MAD_coefs sMAPE_coefs bias_coefs
# 1 102 10000 1 11346.36      52.70173 -8960.633 -0.0002905924 0.8023027  0.06617226   1.004495
#     R2_coefs support_size_true support_size_estimated estimated_over_true_support_size TP FP TN
# 1  0.9988544                51                     49                        0.9607843 49  0 50
#   FN precision    recall specificity accuracy F1_score      gic     ebic      bic      aic
# 1  2         1 0.9607843           1 0.980198     0.98 18434.71 18383.76 18381.78 18021.27
#       mbic
# 1 18705.65



# coefficients & metrics over a range of lambda values
# calculate max lambda where all coefs would be penalized to 0
maxlam <- maxlambda(X=X, weights=L0glm_gamma_presetlambda$weights, 
                    z=L0glm_gamma_presetlambda$z, no.pen=1,
                    coefs=L0glm_gamma_presetlambda$coefficients,
                    family=family)
maxlam # 2269.353
mincoef <- min(abs(L0glm_gamma_presetlambda$coefficients)) # *sign(L0glm_gamma_presetlambda$coefficients)

lambdas <- 10^seq(round(log10(lambda_preset_gamma))-1, log10(maxlam), by=0.1)
coefmat <- Matrix(0, nrow=length(lambdas), ncol=ncol(X))
warmstart <- 0.2 # % by which L0glm_gamma_presetlambda$coefficients will be used as start
system.time(metrics <- data.frame(do.call(rbind, lapply(1:length(lambdas), function (i) { 
  lambda <- lambdas[[i]]
  start <- L0glm_gamma_presetlambda$coefficients*warmstart+start_coefs(X, family)*(1-warmstart)
  # start[L0glm_gamma_presetlambda$coefficients==0] <- 0
  L0glm.out <- L0glm.fit(X = X, 
                         y = y, 
                         family = Gamma(log), 
                         lambda = lambda, 
                         nonnegative = FALSE,
                         no.pen = 1,
                         start = start,
                         solver = "eigen", # eigen, glmnet or osqp, # note: constrained solvers glmnet & osqp better than eigen for very small lambda values
                         method = 4L,
                         # gamma = 2, 
                         thresh = 1E-5,
                         relaxed = FALSE 
  )$coefficients 
  coefmat[i,] <<- as.vector(L0glm.out)
  return( unlist(c(lambda=lambda, get.metrics(beta = L0glm.out, 
                                              X = X, y = y, 
                                              family = Gamma(log), 
                                              beta_true = beta_true))) )
} )))) # 25s
# TODO: with solver method=1L returns error Error in h(simpleError(msg, call)) :
# error in evaluating the argument 'x' in selecting a method for function 'as.matrix': non-conformable arguments
# in Cholesky dense solver. Add some error catching in solvers or lmridge?
matplot(x=log10(metrics$lambda), 
        y=as.matrix(coefmat), type="l", lty=1,
        xlab="Log10(lambda)", ylab="Coefficient estimates", main="Regularisation path L0 Gamma GLM")
abline(v=log10(maxlam), col="red") # max lambda where all coefs would be penalized to 0
library(scales)
abline(h=beta_true, col=alpha(1:6,0.2))

plot(log10(metrics$lambda), metrics$deviance, type="l")
metrics$lambda[which.min(metrics$deviance)] # 0.1 # NOTE this is deviance relative to observed signal, not deviance relative to ground truth though
plot(log10(metrics$lambda), metrics$deviance_true, type="l")
metrics$lambda[which.min(metrics$deviance_true)] # 1.995262 # deviance compared to true signal
plot(log10(metrics$lambda), metrics$loglik, type="l")
metrics$lambda[which.max(metrics$loglik)] # 0.5011872
plot(log10(metrics$lambda), metrics$MD_coefs, type="l")
plot(log10(metrics$lambda), metrics$MAD_coefs, type="l")
metrics$lambda[which.min(metrics$MAD_coefs)] # 1.995262
plot(log10(metrics$lambda), metrics$sMAPE_coefs, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(metrics$sMAPE_coefs)] # 1.995262
plot(log10(metrics$lambda), metrics$bias_coefs, type="l")
metrics$lambda[which.min(abs(metrics$bias_coefs-1))] # 1
plot(log10(metrics$lambda), metrics$R2_coefs, type="l")
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
metrics$lambda[which.max(metrics$accuracy)] # 1.995262
plot(log10(metrics$lambda), metrics$F1_score, type="l")
metrics$lambda[which.max(metrics$F1_score)] # 1.995262
plot(log10(metrics$lambda), metrics$aic, type="l")
metrics$lambda[which.min(metrics$aic)] # 0.6309573
plot(log10(metrics$lambda), metrics$bic, type="l")
metrics$lambda[which.min(metrics$bic)] # 1.995262
plot(log10(metrics$lambda), metrics$gic, type="l")
metrics$lambda[which.min(metrics$gic)] # 1.995262
plot(log10(metrics$lambda), metrics$ebic, type="l")
metrics$lambda[which.min(metrics$ebic)] # 1.995262
plot(log10(metrics$lambda), metrics$mbic, type="l")
metrics$lambda[which(metrics$mbic==min(metrics$mbic))] # 1.995262


# very close to values that were expected to optimise each IC :
# note: using warm starts in the L0glm.fit function shift the
# empirically observed optimal values

preset.lambda(IC = "aic", # options aic, bic, gic, mebic, ebic, mbic
              y = as.matrix(data_gamma$y), 
              X = Matrix(data_gamma$X, sparse=TRUE), 
              family = Gamma(log),
              dispersion = dispersion) # 0.5
preset.lambda(IC = "bic", # options aic, bic, gic, mebic, ebic, mbic
              y = as.matrix(data_gamma$y), 
              X = Matrix(data_gamma$X, sparse=TRUE), 
              family = Gamma(log),
              dispersion = dispersion) # 2.30
preset.lambda(IC = "gic", # options aic, bic, gic, mebic, ebic, mbic
              y = as.matrix(data_gamma$y), 
              X = Matrix(data_gamma$X, sparse=TRUE), 
              family = Gamma(log),
              dispersion = dispersion) # 2.56
preset.lambda(IC = "ebic", # options aic, bic, gic, mebic, ebic, mbic
              y = as.matrix(data_gamma$y), 
              X = Matrix(data_gamma$X, sparse=TRUE), 
              family = Gamma(log),
              dispersion = dispersion) # 2.30
preset.lambda(IC = "mbic", # options aic, bic, gic, mebic, ebic, mbic
              y = as.matrix(data_gamma$y), 
              X = Matrix(data_gamma$X, sparse=TRUE), 
              family = Gamma(log),
              dispersion = dispersion) # 3.91


# 5. LOG LINK INVERSE GAUSSIAN GLM ####

set.seed(1)
data_invgaus  = simulate_glm(n, p, n_nonzero, family = "inverse.gaussian") 
X = normalizeX(data_invgaus$X) # we first L2 norm normalize the columns of X
y = data_invgaus$y
scaling = attr(X, "scaling")
beta_true = data_invgaus$beta_true * scaling
family = inverse.gaussian(log)

# glm fit
system.time(glm_invgaus <- glm(y ~ 0+X, 
                               family=inverse.gaussian(log))) # 0.52s
glm_invgaus$deviance 
glm_invgaus$aic
AIC(glm_invgaus)
BIC(glm_invgaus)
logLik(glm_invgaus)
cor(coef(glm_invgaus), data_invgaus$beta_true) # 0.9954206
plot(x=coef(glm_invgaus), y=beta_true, pch=16, col="blue")

# L0 glm fit with lambda close to zero
system.time(L0glm_invgaus_lambda0 <- L0glm.fit(X = X, 
                                             y = y, 
                                             family = inverse.gaussian(log), 
                                             lambda = 1E-30, 
                                             nonnegative = FALSE,
                                             normalize = FALSE,
                                             no.pen = 1,
                                             start = NULL, # rep(sqrt(10), ncol(data_invgaus$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                                             solver = "eigen", # 
                                             method = 1L, #
                                             thresh = 1E-5
))  # 2.56s

max(abs(L0glm_invgaus_lambda0$coefficients-coef(glm_invgaus))) # 0.001652324
plot(L0glm_invgaus_lambda0$coefficients, coef(glm_invgaus), pch=16, col="blue")
# check if loglik & AIC & BIC are correct
metrics_invgaus = data.frame(do.call(cbind, get.metrics(beta = L0glm_invgaus_lambda0$coefficients,  # some metrics
                                                      X = X, 
                                                      y = y, 
                                                      family = inverse.gaussian(log), 
                                                      beta_true = beta_true)))
metrics_invgaus
logLik(glm_invgaus)-metrics_invgaus$loglik # -7.439357e-07, correct
AIC(glm_invgaus)-metrics_invgaus$aic # 1.487871e-06
BIC(glm_invgaus)-metrics_invgaus$bic # 1.487871e-06


# L0 glm fit with lambda preset to optimise given IC
lambda_preset_invgaus = preset.lambda(IC = crit, # options aic, bic, gic, mebic, ebic, mbic
                                    y = y, 
                                    X = X, 
                                    family = inverse.gaussian(log),
                                    dispersion = dispersion) # dispersion here was 1, here we pass this value to not have to estimate it 
lambda_preset_invgaus # 2.561769

system.time(L0glm_invgaus_presetlambda <- L0glm.fit(X = X, 
                                                  y = y, 
                                                  family = inverse.gaussian(log), 
                                                  lambda = lambda_preset_invgaus,
                                                  nonnegative = FALSE,
                                                  normalize = FALSE,
                                                  no.pen = 1,
                                                  start = NULL, # (abs(coef(glm_invgaus)))^(1/3), # rep(sqrt(3)/4, ncol(data_invgaus$X)), # NULL, # rep(sqrt(10), ncol(data_invgaus$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                                                  solver = "eigen", # 
                                                  method = 4L, #
                                                  thresh = 1E-5)) 
# 0.92s

mean(abs(L0glm_invgaus_presetlambda$coefficients-coef(glm_invgaus))) # 0.3814394
mean(abs(L0glm_invgaus_presetlambda$coefficients-beta_true)) # 0.3533345
mean(abs(coef(glm_invgaus)-beta_true)) # 0.6385165

plot(x=L0glm_invgaus_presetlambda$coefficients, y=beta_true, pch=16, col="blue")
points(x=L0glm_invgaus_presetlambda$coefficients[beta_true==0], y=beta_true[beta_true==0], pch=16, col="red")


# check metrics & quality of solution with preset lambda
system.time(metrics_invgaus <- data.frame(do.call(cbind, get.metrics(beta = L0glm_invgaus_presetlambda$coefficients,  # some metrics
                                                                   X = X, 
                                                                   y = y,
                                                                   family = inverse.gaussian(log), 
                                                                   beta_true = beta_true)))) 
metrics_invgaus
#     p     n m deviance deviance_true    loglik  MD_coefs MAD_coefs sMAPE_coefs bias_coefs
# 1 102 10000 1 9988.176      49.01657 -3258.538 0.0123419 0.3533345  0.03335892  0.9985487
#    R2_coefs support_size_true support_size_estimated estimated_over_true_support_size TP FP TN
# 1 0.9997829                51                     50                        0.9803922 50  0 50
# FN precision    recall specificity accuracy F1_score      gic     ebic      bic      aic
#  1         1         1 0.9803922           1 0.990099 0.990099 7040.793 6988.824 6986.804 6619.076
# mbic
# 1 7317.149



# coefficients & metrics over a range of lambda values
# calculate max lambda where all coefs would be penalized to 0
maxlam <- maxlambda(X=X, weights=L0glm_invgaus_presetlambda$weights, 
                    z=L0glm_invgaus_presetlambda$z, no.pen=1,
                    coefs=L0glm_invgaus_presetlambda$coefficients,
                    family=family)
maxlam # 232265.1
mincoef <- min(abs(L0glm_invgaus_presetlambda$coefficients)) # *sign(L0glm_invgaus_presetlambda$coefficients)

lambdas <- 10^seq(round(log10(lambda_preset_invgaus))-1, log10(maxlam), by=0.1)
coefmat <- Matrix(0, nrow=length(lambdas), ncol=ncol(X))
warmstart <- 0.2 # % by which L0glm_invgaus_presetlambda$coefficients will be used as start, lowest value that led to stable coefficient paths
system.time(metrics <- data.frame(do.call(rbind, lapply(1:length(lambdas), function (i) { 
  lambda <- lambdas[[i]]
  start <- (L0glm_invgaus_presetlambda$coefficients*warmstart+start_coefs(X, family)*(1-warmstart)) # rep(sqrt(1),ncol(X))
  # start[L0glm_invgaus_presetlambda$coefficients==0] <- 0
  L0glm.out <- L0glm.fit(X = X, 
                         y = y, 
                         family = inverse.gaussian(log), 
                         lambda = lambda, 
                         nonnegative = FALSE,
                         normalize = FALSE,
                         no.pen = 1,
                         start = start, # coef(glm_invgaus), # NULL, # rep(sqrt(10), ncol(data_invgaus$X)),
                         solver = "eigen", # eigen, glmnet or osqp, # note: constrained solvers glmnet & osqp better than eigen for very small lambda values
                         method = 4L,
                         # gamma = 2, 
                         thresh = 1E-5,
                         relaxed = FALSE 
  )$coefficients 
  coefmat[i,] <<- as.vector(L0glm.out)
  return( unlist(c(lambda=lambda, get.metrics(beta = L0glm.out, 
                                              X = X, y = y, 
                                              family = inverse.gaussian(log), 
                                              beta_true = beta_true))) )
} )))) # 32s
coefmat = as.matrix(coefmat)
coefmat[abs(coefmat)>2*max(abs(beta_true))] = 0
# TODO check where instability in coefficients at high lambda values comes from
matplot(x=log10(metrics$lambda), 
        y=as.matrix(coefmat), type="l", lty=1,
        xlab="Log10(lambda)", ylab="Coefficient estimates", main="Regularisation path L0 Inverse Gaussian GLM")
abline(v=log10(maxlam), col="red") # max lambda where all coefs would be penalized to 0
library(scales)
abline(h=beta_true, col=alpha(1:6,0.2))

# TODO check why values explode for lambda = 50 and 200

plot(log10(metrics$lambda), metrics$deviance, type="l")
metrics$lambda[which.min(metrics$deviance)] # 0.1 # NOTE this is deviance relative to observed signal, not deviance relative to ground truth though
plot(log10(metrics$lambda), metrics$deviance_true, type="l")
metrics$lambda[which.min(metrics$deviance_true)] # 1 # deviance compared to true signal
plot(log10(metrics$lambda), metrics$loglik, type="l")
metrics$lambda[which.max(metrics$loglik)] # 0.5011872
plot(log10(metrics$lambda), metrics$MD_coefs, type="l")
plot(log10(metrics$lambda), metrics$MAD_coefs, type="l")
metrics$lambda[which.min(metrics$MAD_coefs)] # 1
plot(log10(metrics$lambda), metrics$sMAPE_coefs, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(metrics$sMAPE_coefs)] # 1
plot(log10(metrics$lambda), metrics$bias_coefs, type="l")
metrics$lambda[which.min(abs(metrics$bias_coefs-1))] # 6.309573
plot(log10(metrics$lambda), metrics$R2_coefs, type="l")
plot(log10(metrics$lambda), metrics$estimated_over_true_support_size, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(abs(metrics$estimated_over_true_support_size-1))] # 0.7943282
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
metrics$lambda[which(metrics$mbic==min(metrics$mbic))] # 3.981072

# reasonable match to lambda values that were expected to optimise each IC :
# note: using warm starts in the L0glm.fit function shift the
# empirically observed optimal values

preset.lambda(IC = "aic", # options aic, bic, gic, mebic, ebic, mbic
              y = as.matrix(data_invgaus$y), 
              X = Matrix(data_invgaus$X, sparse=TRUE), 
              family = inverse.gaussian(log),
              dispersion = dispersion) # 0.5
preset.lambda(IC = "bic", # options aic, bic, gic, mebic, ebic, mbic
              y = as.matrix(data_invgaus$y), 
              X = Matrix(data_invgaus$X, sparse=TRUE), 
              family = inverse.gaussian(log),
              dispersion = dispersion) # 2.30
preset.lambda(IC = "gic", # options aic, bic, gic, mebic, ebic, mbic
              y = as.matrix(data_invgaus$y), 
              X = Matrix(data_invgaus$X, sparse=TRUE), 
              family = inverse.gaussian(log),
              dispersion = dispersion) # 2.56
preset.lambda(IC = "ebic", # options aic, bic, gic, mebic, ebic, mbic
              y = as.matrix(data_invgaus$y), 
              X = Matrix(data_invgaus$X, sparse=TRUE), 
              family = inverse.gaussian(log),
              dispersion = dispersion) # 2.30
preset.lambda(IC = "mbic", # options aic, bic, gic, mebic, ebic, mbic
              y = as.matrix(data_invgaus$y), 
              X = Matrix(data_invgaus$X, sparse=TRUE), 
              family = inverse.gaussian(log),
              dispersion = dispersion) # 3.91