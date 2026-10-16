# TEST NONNEGATIVE IDENTITY LINK POISSON GLM glmnet, L0glm, abess ####

# Load required packages
library(glmnet)
library(L0glm)
library(abess)

# simulate a blurred spike train with Poisson noise
sim = simul(method=1, n=1000, p=1000, s=0.1, sparse=T, seed=1)
X <- sim$X # banded matrix with time-shifted copies of Gaussian peak shapes (sparse)
Y <- sim$y
beta <- c(0, as.vector(sim$beta_true)) # intercept is zero


# 1. fit nonnegative identity link Poisson LASSO model ####
# Cross-validation to select the best lambda
system.time(cvfit1 <- cv.glmnet(X, Y, alpha=1, family = poisson("identity"), lower.limits = 0, intercept=T))
# 69s
plot(cvfit1)
optlambda1 <- cvfit1$lambda.min
# Get coefficients at optimal lambda
system.time(fit1 <- glmnet(X, Y, alpha=1, family = poisson(link = "identity"), lower.limits = 0, intercept=T))
# 8s 
# Plot the coefficient paths
plot(fit1, xvar = "lambda", label = TRUE)

coefs1 <- as.vector(coef(fit1, s = optlambda1))
cor(coefs1, beta) # 0.9812176
table(beta!=0, coefs1!=0)
#       FALSE TRUE
# FALSE   793  108
# TRUE      2   98



# 2. fit nonnegative weighted least square LASSO model approximating nonnegative identity link Poisson LASSO ####
# Cross-validation to select the best lambda
system.time(cvfit2 <- cv.glmnet(X, Y, alpha=1, 
                   family = "gaussian",
                   weights = 1/(Y+0.1), # approx 1/variance weights for identity link Poisson
                   lower.limits = 0, intercept=T))
# 0.2s
plot(cvfit2)
optlambda2 <- cvfit2$lambda.min
# Get coefficients at optimal lambda
system.time(fit2 <- glmnet(X, Y, alpha=1, 
              family = "gaussian", 
              weights = 1/(Y+0.1), # approx 1/variance weights for identity link Poisson
              lower.limits = 0, intercept=T))
# 0.02s
# Plot the coefficient paths
plot(fit2, xvar = "lambda", label = TRUE)

coefs2 <- as.vector(coef(fit2, s = optlambda2))
cor(coefs2, beta) # 0.928145
table(beta!=0, coefs2!=0)
#       FALSE TRUE
# FALSE   767  134
# TRUE      6   94
plot(coefs2, beta)


# 3. fit nonnegative identity link Poisson relaxed LASSO model ####
# Cross-validation to select the best lambda
system.time(cvfit3 <- cv.glmnet(X, Y, alpha=1, 
                                gamma=c(0, 0.25, 0.5, 0.75, 1), relax=T,  # using relaxed LASSO here
                                family = poisson(link="identity"),
                                lower.limits = 0, intercept=T))
# 274s
plot(cvfit3)
optlambda3 <- cvfit3$relaxed$lambda.min
optgamma3 <- cvfit3$relaxed$gamma.min 
optgamma3 # 0

# Get coefficients at optimal lambda
system.time(fit3 <- glmnet(X, Y, alpha=1, 
                           gamma=optgamma3, relax=T,  # using relaxed LASSO here
                           family = poisson(link="identity"), 
                           lower.limits = 0, intercept=T))
# 27s
# Plot the coefficient paths
plot(fit3, xvar = "lambda", gamma = optgamma3, label = TRUE)

coefs3 <- as.vector(coef(fit3, s = optlambda3, gamma = optgamma3))
cor(coefs3, beta) # 0.9884491
table(beta!=0, coefs3!=0)
#       FALSE TRUE
# FALSE   812   89
# TRUE      2   98
plot(coefs3, beta)
# good solution but slooooowwww


# 4. fit nonnegative weighted least square relaxed LASSO model approximating nonnegative identity link Poisson relaxed LASSO ####
# Cross-validation to select the best lambda
system.time(cvfit4 <- cv.glmnet(X, Y, alpha=1, 
                                gamma=c(0, 0.25, 0.5, 0.75, 1), relax=T,  # using relaxed LASSO here
                                family = "gaussian",
                                weights = 1/(Y+0.1), # approx 1/variance weights for identity link Poisson
                                lower.limits = 0, intercept=T))
# 2.78s
plot(cvfit4)
optlambda4 <- cvfit4$relaxed$lambda.min
optgamma4 <- cvfit4$relaxed$gamma.min 
optgamma4 # 0

# Get coefficients at optimal lambda
system.time(fit4 <- glmnet(X, Y, alpha=1, 
                           gamma=optgamma4, relax=T,  # using fully relaxed LASSO here
                           family = "gaussian", 
                           weights = 1/(Y+0.1), # approx 1/variance weights for identity link Poisson
                           lower.limits = 0, intercept=T))
# 0.18s
# Plot the coefficient paths
plot(fit4, xvar = "lambda", gamma = optgamma4, label = TRUE)

coefs4 <- as.vector(coef(fit4, s = optlambda4, gamma = optgamma4))
cor(coefs4, beta) # 0.9206419
table(beta!=0, coefs4!=0)
#       FALSE TRUE
# FALSE   824   77
# TRUE     16   84
plot(coefs4, beta)


# fit nonnegative L0 penalized identity link Poisson GLM, fit using an iterative adaptive ridge procedure ####
system.time(L0glmfit <- L0glm.fit(X=X, y=Y, 
                      family=poisson(identity), 
                      nonnegative=TRUE, 
                      start = coefs3[-1],
                      lambda='gic'))
# 0.03s
coefs5 <- as.vector(L0glmfit$coefficients)
cor(coefs5, beta[-1]) # 0.9926351
table(beta[-1]!=0, coefs5!=0)
#       FALSE TRUE
# FALSE   893    7
# TRUE      7   93
plot(coefs5, beta[-1])
# BEST SOLUTION, but best solution obtained if initialised with coefs1 or coefs3 (which is slow)

system.time(L0glmfit <- L0glm.fit(X=X, y=Y, 
                                  family=poisson(identity), 
                                  nonnegative=TRUE, 
                                  start = rep(1, ncol(X)), 
                                  lambda='gic'))
# 0.05s
coefs5b <- as.vector(L0glmfit$coefficients)
cor(coefs5b, beta[-1]) # 0.9304038
table(beta[-1]!=0, coefs5b!=0)
#       FALSE TRUE
# FALSE   864   36
# TRUE     16   84
plot(coefs5b, beta[-1])




# fit nonnegative L0 penalized weighted least square GLM, fit using an iterative adaptive ridge procedure ####
system.time(L0glmfit2 <- L0glm.fit(X=X, y=Y, weights=(1/(Y+0.1)), 
                                  family=gaussian(identity), 
                                  nonnegative=TRUE, 
                                  start = coefs3[-1], # NULL, # rep(1, ncol(X)), # 
                                  lambda='gic'))
# 0.01s
coefs6 <- as.vector(L0glmfit2$coefficients)
cor(coefs6, beta[-1]) # 0.9802931
table(beta[-1]!=0, coefs6!=0)
#       FALSE TRUE
# FALSE   894    8
# TRUE     11   89
plot(coefs6, beta[-1])




# fit nonnegative L0 penalized weighted least square GLM, fit using abess ####
system.time(abessfit <- abess(x=X, y=as.vector(Y), weight=(1/(as.vector(Y)+0.1)), 
                                   family="gaussian", 
                                   normalize=F,
                                   beta.low=0, 
                                   fit.intercept = F,
                                   tune.type='gic'))
# 0.06s
coefs7 <- as.vector(coef(abessfit, support.size = abessfit[["best.size"]]))
cor(coefs7, beta) # 0.6050387
table(beta!=0, coefs7!=0)
#       FALSE TRUE
# FALSE   874   27
# TRUE     86   14
plot(coefs7, beta)
# poor solution


