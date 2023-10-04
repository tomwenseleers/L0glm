# BENCHMARKS L0GLM ####

library(remotes)
# remotes::install_github("jaredhuling/ordinis")
# remotes::install_github("tomwenseleers/L0glm/L0glm")
# branch of L0Learn with nonnegativity constraints
# remotes::install_github("hazimehh/L0Learn", ref="nng") 
# library(L0glm)
devtools::load_all(".")
library(export)
library(Matrix)
library(ordinis)
library(L0Learn)
library(abess)
library(glmnet)


# 1. SINGLE CHANNEL ####
# Simulate SC n=p=500
n <- p <- 500
sim.gaus <- simul(method = 1, n = n, p = p, sparse = T)
sim.pois <- simul(method = 2, n = n, p = p, sparse = T)

## 1.1. Poisson noise ####
## Tune lambda relaxed = F with gic
system.time(tune.osqp.sc <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "gic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = FALSE, family = poisson(identity), nonnegative = TRUE, solver = "osqp"))
# bestlambda = 10, time = 2.7s
system.time(tune.glmnet.sc <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "gic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = FALSE, family = poisson(identity), nonnegative = TRUE, solver = "glmnet"))
# bestlambda = 10, time = 8.1s
system.time(tune.chol.sc <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "gic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = FALSE, family = poisson(identity), nonnegative = TRUE, solver = "eigen", method = 7L))
# bestlambda = 10, time = 1.1s
system.time(tune.QR.sc <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "gic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = FALSE, family = poisson(identity), nonnegative = TRUE, solver = "eigen", method = 1L))
# bestlambda = 10, time = 1.6s
system.time(tune.CG.sc <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "ebic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = FALSE, family = poisson(identity), nonnegative = TRUE, solver = "eigen", method = 5L))
# bestlambda = 10, time = 1.5s
system.time(tune.LSCG.sc <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "ebic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = FALSE, family = poisson(identity), nonnegative = TRUE, solver = "eigen", method = 3L))
# bestlambda = 10, time = 1.9s

## 1.1.1 Fit relaxed = FALSE ####
# OSQP ####
system.time(L0glm.osqp <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                    family = poisson(identity), 
                                    lambda = "ebic", # tune.osqp.sc$best.lambda, # choose preset or tuned
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "osqp")) 
# 0.09s
cor(as.vector(as.vector(L0glm.osqp$coefficients)), as.vector(sim.pois$beta_true)) # 0.9780096
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.osqp$coefficients) != 0)) # TP = 45
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.osqp$coefficients) != 0)) # FP = 7
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.osqp$coefficients) == 0)) # TN = 443
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.osqp$coefficients) == 0)) # FN = 5

# GLMNET ####
system.time(L0glm.glmnet <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                      family = poisson(identity), 
                                      lambda = "ebic", # tune.glmnet.sc$best.lambda,
                                      nonnegative = TRUE,
                                      relaxed = FALSE,
                                      solver = "glmnet")) 
# 0.36s
cor(as.vector(as.vector(L0glm.glmnet$coefficients)), as.vector(sim.pois$beta_true)) # 0.9777347
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.glmnet$coefficients) != 0)) # TP = 45
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.glmnet$coefficients) != 0)) # FP = 7
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.glmnet$coefficients) == 0)) # TN = 443
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.glmnet$coefficients) == 0)) # FN = 5

# QR ####
system.time(L0glm.QR <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                  family = poisson(identity), 
                                  lambda = "ebic", # tune.QR.sc$best.lambda,
                                  nonnegative = TRUE,
                                  relaxed = FALSE,
                                  solver = "eigen",
                                  method = 1L)) 
# 0.07s
cor(as.vector(as.vector(L0glm.QR$coefficients)), as.vector(sim.pois$beta_true)) # 0.9781847
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.QR$coefficients) != 0)) # TP = 45
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.QR$coefficients) != 0)) # FP = 7
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.QR$coefficients) == 0)) # TN = 443
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.osqp$coefficients) == 0)) # FN = 5

# Cholesky ####
system.time(L0glm.chol <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                    family = poisson(identity), 
                                    lambda = "ebic", # tune.chol.sc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "eigen",
                                    method = 7L)) 
# 0.03s
cor(as.vector(as.vector(L0glm.chol$coefficients)), as.vector(sim.pois$beta_true)) # 0.9781847
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.chol$coefficients) != 0)) # TP = 45
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.chol$coefficients) != 0)) # FP = 7
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.chol$coefficients) == 0)) # TN = 443
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.chol$coefficients) == 0)) # FN = 5

# LSCG ####
system.time(L0glm.LSCG <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                    family = poisson(identity), 
                                    lambda = "ebic", # tune.LSCG.sc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "eigen",
                                    method = 3L)) 
# 0.16s
cor(as.vector(as.vector(L0glm.LSCG$coefficients)), as.vector(sim.pois$beta_true)) # 0.9790543
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.LSCG$coefficients) != 0)) # TP = 45
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.LSCG$coefficients) != 0)) # FP = 8
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.LSCG$coefficients) == 0)) # TN = 442
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.LSCG$coefficients) == 0)) # FN = 5

# CG ####
system.time(L0glm.CG <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                  family = poisson(identity), 
                                  lambda = "ebic", # tune.CG.sc$best.lambda,
                                  nonnegative = TRUE,
                                  relaxed = TRUE,
                                  solver = "eigen",
                                  method = 5L)) 
# 0.09s
cor(as.vector(as.vector(L0glm.CG$coefficients)), as.vector(sim.pois$beta_true)) # 0.9796666
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.CG$coefficients) != 0)) # TP = 45
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.CG$coefficients) != 0)) # FP = 8
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.CG$coefficients) == 0)) # TN = 442
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.CG$coefficients) == 0)) # FN = 5

## Tune lambda relaxed = T with gic
system.time(tune.osqp.sc <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "gic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = TRUE, family = poisson(identity), nonnegative = TRUE, solver = "osqp"))
# bestlambda = 10, time = 2.9s
system.time(tune.glmnet.sc <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "gic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = TRUE, family = poisson(identity), nonnegative = TRUE, solver = "glmnet"))
# bestlambda = 10, time = 9.2s
system.time(tune.chol.sc <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "gic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = TRUE, family = poisson(identity), nonnegative = TRUE, solver = "eigen", method = 7L))
# bestlambda = 10, time = 1.3s
system.time(tune.QR.sc <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "gic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = TRUE, family = poisson(identity), nonnegative = TRUE, solver = "eigen", method = 1L))
# bestlambda = 10, time = 2.0s
system.time(tune.CG.sc <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "gic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = TRUE, family = poisson(identity), nonnegative = TRUE, solver = "eigen", method = 5L))
# bestlambda = 10, time = 1.6s
system.time(tune.LSCG.sc <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "gic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = TRUE, family = poisson(identity), nonnegative = TRUE, solver = "eigen", method = 3L))
# bestlambda = 10, time = 2.2s

## 1.1.2 Fit relaxed = TRUE ####
# OSQP ####
system.time(L0glm.osqp <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                    family = poisson(identity), 
                                    lambda = "ebic", # tune.osqp.sc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = TRUE,
                                    solver = "osqp")) 
# 0.09s
cor(as.vector(as.vector(L0glm.osqp$coefficients)), as.vector(sim.pois$beta_true)) # 0.9764704
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.osqp$coefficients) != 0)) # TP = 45
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.osqp$coefficients) != 0)) # FP = 7
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.osqp$coefficients) == 0)) # TN = 443
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.osqp$coefficients) == 0)) # FN = 5


# GLMNET ####
system.time(L0glm.glmnet <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                      family = poisson(identity), 
                                      lambda = "ebic", # tune.glmnet.sc$best.lambda,
                                      nonnegative = TRUE,
                                      relaxed = TRUE,
                                      solver = "glmnet")) 
# 0.4s
cor(as.vector(as.vector(L0glm.glmnet$coefficients)), as.vector(sim.pois$beta_true)) # 0.9762409
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.glmnet$coefficients) != 0)) # TP = 45
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.glmnet$coefficients) != 0)) # FP = 7
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.glmnet$coefficients) == 0)) # TN = 443
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.glmnet$coefficients) == 0)) # FN = 5

# QR ####
system.time(L0glm.QR <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                  family = poisson(identity), 
                                  lambda = "ebic", # tune.QR.sc$best.lambda,
                                  nonnegative = TRUE,
                                  relaxed = TRUE,
                                  solver = "eigen",
                                  method = 1L)) 
# 0.06s
cor(as.vector(as.vector(L0glm.QR$coefficients)), as.vector(sim.pois$beta_true)) # 0.9768456
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.QR$coefficients) != 0)) # TP = 45
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.QR$coefficients) != 0)) # FP = 7
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.QR$coefficients) == 0)) # TN = 443
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.osqp$coefficients) == 0)) # FN = 5

# Cholesky ####
system.time(L0glm.chol <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                    family = poisson(identity), 
                                    lambda = "ebic", # tune.chol.sc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = TRUE,
                                    solver = "eigen",
                                    method = 7L)) 
# 0.04s
cor(as.vector(as.vector(L0glm.chol$coefficients)), as.vector(sim.pois$beta_true)) # 0.9768456
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.chol$coefficients) != 0)) # TP = 45
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.chol$coefficients) != 0)) # FP = 7
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.chol$coefficients) == 0)) # TN = 443
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.chol$coefficients) == 0)) # FN = 5

# LSCG ####
system.time(L0glm.LSCG <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                    family = poisson(identity), 
                                    lambda = "ebic", # tune.LSCG.sc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = TRUE,
                                    solver = "eigen",
                                    method = 3L)) 
# 0.18s
cor(as.vector(as.vector(L0glm.LSCG$coefficients)), as.vector(sim.pois$beta_true)) # 0.9796581
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.LSCG$coefficients) != 0)) # TP = 45
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.LSCG$coefficients) != 0)) # FP = 8
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.LSCG$coefficients) == 0)) # TN = 442
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.LSCG$coefficients) == 0)) # FN = 5

# CG ####
system.time(L0glm.CG <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                  family = poisson(identity), 
                                  lambda = "ebic", # tune.CG.sc$best.lambda,
                                  nonnegative = TRUE,
                                  relaxed = TRUE,
                                  solver = "eigen",
                                  method = 5L)) 
# 0.08s
cor(as.vector(as.vector(L0glm.CG$coefficients)), as.vector(sim.pois$beta_true)) # 0.9796666
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.CG$coefficients) != 0)) # TP = 45
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.CG$coefficients) != 0)) # FP = 8
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.CG$coefficients) == 0)) # TN = 442
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.CG$coefficients) == 0)) # FN = 5



## 1.2. Gaussian noise ####

# TODO: benchmarks below return error
# FIX BUG
# Error in lm.fit(x, y, offset = offset, singular.ok = singular.ok, ...) : 
# 0 (non-NA) cases

## Tune lambda relaxed = F 
system.time(tune.osqp.sc <- L0glm.tune(X=sim.gaus$X, y=sim.gaus$y, IC = "gic", plot = FALSE, coefs.true = sim.gaus$beta_true, relaxed = FALSE, family = gaussian(identity), nonnegative = TRUE, solver = "osqp"))
# bestlambda = 3.162277, time = 4.614s
system.time(tune.glmnet.sc <- L0glm.tune(X=sim.gaus$X, y=sim.gaus$y, IC = "gic", plot = FALSE, coefs.true = sim.gaus$beta_true, relaxed = FALSE, family = gaussian(identity), nonnegative = TRUE, solver = "glmnet"))
# bestlambda = 316.2277, time = 5.411s
system.time(tune.chol.sc <- L0glm.tune(X=sim.gaus$X, y=sim.gaus$y, IC = "gic", plot = FALSE, coefs.true = sim.gaus$beta_true, relaxed = FALSE, family = gaussian(identity), nonnegative = TRUE, solver = "eigen", method = 7L))
# bestlambda = 3.162277, time = 3.052s
system.time(tune.QR.sc <- L0glm.tune(X=sim.gaus$X, y=sim.gaus$y, IC = "gic", plot = FALSE, coefs.true = sim.gaus$beta_true, relaxed = TRUE, family = gaussian(identity), nonnegative = TRUE, solver = "eigen", method = 1L))
# bestlambda = 3.162277, time = 15.925s
system.time(tune.CG.sc <- L0glm.tune(X=sim.gaus$X, y=sim.gaus$y, IC = "gic", plot = FALSE, coefs.true = sim.gaus$beta_true, relaxed = FALSE, family = gaussian(identity), nonnegative = TRUE, solver = "eigen", method = 5L))
# bestlambda = 3.162277, time = 7.951s
system.time(tune.LSCG.sc <- L0glm.tune(X=sim.gaus$X, y=sim.gaus$y, IC = "gic", plot = TRUE, coefs.true = sim.gaus$beta_true, relaxed = FALSE, family = gaussian(identity), nonnegative = TRUE, solver = "eigen", method = 3L))
# bestlambda = 3.162277, time = 10.929s

## 1.2.1 Fit relaxed = FALSE ####
# OSQP ####
system.time(L0glm.osqp <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                    family = gaussian(identity), 
                                    lambda = "gic", # tune.osqp.sc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "osqp")) 
# 0.2s
cor(as.vector(as.vector(L0glm.osqp$coefficients)), as.vector(sim.gaus$beta_true)) #  0.9836706
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.osqp$coefficients) != 0)) # TP = 48
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.osqp$coefficients) != 0)) # FP = 7
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.osqp$coefficients) == 0)) # TN = 443
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.osqp$coefficients) == 0)) # FN = 2

# GLMNET ####
system.time(L0glm.glmnet <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                      family = gaussian(identity), 
                                      lambda = "gic", # tune.glmnet.sc$best.lambda,
                                      nonnegative = TRUE,
                                      relaxed = FALSE,
                                      solver = "glmnet")) 
# 1.5s
cor(as.vector(as.vector(L0glm.glmnet$coefficients)), as.vector(sim.gaus$beta_true)) # 0.9637935
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.glmnet$coefficients) != 0)) # TP = 47
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.glmnet$coefficients) != 0)) # FP = 24
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.glmnet$coefficients) == 0)) # TN = 426
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.glmnet$coefficients) == 0)) # FN = 3

# QR ####
system.time(L0glm.QR <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                  family = gaussian(identity), 
                                  lambda = "gic", # tune.QR.sc$best.lambda,
                                  nonnegative = TRUE,
                                  relaxed = FALSE,
                                  solver = "eigen",
                                  method = 1L)) 
# 0.04s
cor(as.vector(as.vector(L0glm.QR$coefficients)), as.vector(sim.gaus$beta_true)) # 0.997474
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.QR$coefficients) != 0)) # TP = 49
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.QR$coefficients) != 0)) # FP = 4
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.QR$coefficients) == 0)) # TN = 446
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.osqp$coefficients) == 0)) # FN = 2

# Cholesky ####
system.time(L0glm.chol <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                    family = gaussian(identity), 
                                    lambda = "gic", # tune.chol.sc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "eigen",
                                    method = 7L)) 
# 0.03s
cor(as.vector(as.vector(L0glm.chol$coefficients)), as.vector(sim.gaus$beta_true)) # 0.997474
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.chol$coefficients) != 0)) # TP = 49
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.chol$coefficients) != 0)) # FP = 4
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.chol$coefficients) == 0)) # TN = 446
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.chol$coefficients) == 0)) # FN = 1

# LSCG ####
system.time(L0glm.LSCG <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                    family = gaussian(identity), 
                                    lambda = "gic", # tune.CG.sc$best.lambda,
                                    nonnegative = TRUE,
                                    start = rep(1, 500),
                                    relaxed = FALSE,
                                    solver = "eigen",
                                    method = 3L)) 
# 0.12s
cor(as.vector(as.vector(L0glm.LSCG$coefficients)), as.vector(sim.gaus$beta_true)) # 0.9777964
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.LSCG$coefficients) != 0)) # TP = 47
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.LSCG$coefficients) != 0)) # FP = 12
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.LSCG$coefficients) == 0)) # TN = 438
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.LSCG$coefficients) == 0)) # FN = 3

# CG ####
system.time(L0glm.CG <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                  family = gaussian(identity), 
                                  lambda = "gic", # tune.CG.sc$best.lambda,
                                  nonnegative = TRUE,
                                  relaxed = TRUE,
                                  solver = "eigen",
                                  method = 5L)) 
# 0.06s
cor(as.vector(as.vector(L0glm.CG$coefficients)), as.vector(sim.gaus$beta_true)) # 0.9871489
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.CG$coefficients) != 0)) # TP = 46
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.CG$coefficients) != 0)) # FP = 7
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.CG$coefficients) == 0)) # TN = 443
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.CG$coefficients) == 0)) # FN = 4


# 2. MULTICHANNEL ####
# Simulate MC n=p=500
n <- p <- 500
sim.gaus.mc <- simul(method = 3, n = n, p = p, seed = 1)
sim.pois.mc <- simul(method = 4, n = n, p = p, seed = 2)


## 2.1. Poisson noise ####
## 2.1.1 Relaxed = FALSE, grouped = FALSE ####
## Tune lambda relaxed = FALSE, grouped = FALSE
system.time(tune.osqp.mc <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, IC = "gic", 
                                       plot = FALSE, coefs.true = sim.pois.mc$beta_true, 
                                       relaxed = FALSE, family = poisson(identity), 
                                       nonnegative = TRUE, 
                                       parallelize = TRUE, grouped = FALSE, solver = "osqp"))
# bestlambda = 3.16227766016838, time = 127s

system.time(tune.chol.mc <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, IC = "gic", 
                                       plot = FALSE, coefs.true = sim.pois.mc$beta_true, 
                                       relaxed = FALSE, family = poisson(identity), 
                                       nonnegative = TRUE, 
                                       parallelize = TRUE, grouped = FALSE,
                                       solver = "eigen", method = 7L))
# bestlambda = 3.16227766016838, time = 113s

## Fit relaxed = FALSE, grouped = FALSE
# OSQP ####
system.time(L0glm.osqp <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                    family = poisson(identity), 
                                    lambda = "ebic", # tune.osqp.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallelize = TRUE,
                                    grouped = FALSE,
                                    solver = "osqp")) 
# 7.1s
abline(v=which(rowMin(L0glm.osqp$coefficients)!=0), col = 'cyan')

cor(as.vector(as.vector(L0glm.osqp$coefficients)), as.vector(sim.pois.mc$beta_true)) # 0.9941838
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)!=0) # TP = 50
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)!=0) # FP = 169
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)==0) # TN = 281
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)==0) # FN = 0

# Cholesky ####
system.time(L0glm.chol <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                    family = poisson(identity), 
                                    lambda = "ebic", # tune.chol.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallelize = TRUE,
                                    grouped = FALSE,
                                    solver = "eigen",
                                    method = 7L)) 
# 6.5s
abline(v=which(rowMin(L0glm.chol$coefficients)!=0), col = 'cyan')

R2 = round(cor(as.vector(rowMeans(sim.pois.mc$beta_true)), as.vector(rowMeans(L0glm.chol$coefficients))), digits = 4)^2 # Calculate R2
title = paste("True vs estimated coefficients\nnonnegative L0 independent GLM\nR2=", R2)
plot(as.vector(rowMeans(L0glm.chol$coefficients)), as.vector(rowMeans(sim.pois.mc$beta_true)), pch=16, 
     xlab="estimated coefficients nonnegative L0 group GLM",
     ylab="true coefficients", col="steelblue",
     main=title)

cor(as.vector(as.vector(L0glm.chol$coefficients)), as.vector(sim.pois.mc$beta_true)) # 0.9946084
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.chol$coefficients)!=0) # TP = 50
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.chol$coefficients)!=0) # FP = 169
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.chol$coefficients)==0) # TN = 281
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.chol$coefficients)==0) # FN = 0


## 2.1.2 Relaxed = FALSE, grouped = TRUE ####
## Tune lambda relaxed = FALSE, grouped = TRUE
system.time(tune.osqp.mc <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, IC = "gic", 
                                       plot = FALSE, coefs.true = sim.pois.mc$beta_true, 
                                       relaxed = FALSE, family = poisson(identity), 
                                       nonnegative = TRUE, 
                                       parallelize = TRUE, grouped = TRUE, solver = "osqp"))

# bestlambda = 3.16227766016838, time = 253s
system.time(tune.chol.mc <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, IC = "gic", 
                                       plot = FALSE, coefs.true = sim.pois.mc$beta_true, 
                                       relaxed = FALSE, family = poisson(identity), 
                                       nonnegative = TRUE, 
                                       parallelize = TRUE, grouped = TRUE,
                                       solver = "eigen", method = 7L))
# bestlambda = 3.16227766016838, time = 190s

## Fit relaxed = FALSE, grouped = TRUE
# OSQP ####
system.time(L0glm.osqp <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                    family = poisson(identity), 
                                    lambda = "ebic", # tune.osqp.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallelize = TRUE,
                                    grouped = TRUE,
                                    solver = "osqp"))
# 14s
abline(v=which(rowMin(L0glm.osqp$coefficients)!=0), col = 'cyan')

cor(as.vector(as.vector(L0glm.osqp$coefficients)), as.vector(sim.pois.mc$beta_true)) # 0.9992435
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)!=0) # TP = 49
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)!=0) # FP = 1
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)==0) # TN = 449
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)==0) # FN = 1

# Cholesky ####
system.time(L0glm.chol <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                    family = poisson(identity), 
                                    lambda = "ebic", # tune.chol.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallelize = TRUE,
                                    grouped = TRUE,
                                    solver = "eigen",
                                    method = 7L)) 
# 12s
abline(v=which(rowMin(L0glm.chol$coefficients)!=0), col = 'cyan')

cor(as.vector(as.vector(L0glm.chol$coefficients)), as.vector(sim.pois.mc$beta_true)) # 0.9992436
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.chol$coefficients)!=0) # TP = 49
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.chol$coefficients)!=0) # FP = 1
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.chol$coefficients)==0) # TN = 449
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.chol$coefficients)==0) # FN = 1

## 2.1.3 Relaxed = TRUE, grouped = TRUE ####
# TODO: this should fit the subsequent iterations using grouped=TRUE & relaxed=FALSE
# & only as an extra iteration then refit using grouped=TRUE & relaxed=TRUE
# check that this is how it is done now

## Tune lambda relaxed = TRUE, grouped = TRUE
system.time(tune.osqp.mc <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, IC = "gic", 
                                       plot = FALSE, coefs.true = sim.pois.mc$beta_true, 
                                       relaxed = TRUE, family = poisson(identity), 
                                       nonnegative = TRUE, 
                                       parallelize = TRUE, grouped = TRUE, solver = "osqp"))

# bestlambda = 3.16227766016838, time = 275s
system.time(tune.chol.mc <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, IC = "gic", 
                                       plot = FALSE, coefs.true = sim.pois.mc$beta_true, 
                                       relaxed = TRUE, family = poisson(identity), 
                                       nonnegative = TRUE, 
                                       parallelize = TRUE, grouped = TRUE,
                                       solver = "eigen", method = 7L))
# bestlambda = 3.16227766016838, time = 206s

## Fit relaxed = TRUE, grouped = TRUE
# OSQP ####
system.time(L0glm.osqp <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                    family = poisson(identity), 
                                    lambda = "ebic", # tune.osqp.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = TRUE,
                                    parallelize = TRUE,
                                    grouped = TRUE,
                                    solver = "osqp")) 
# 15s
abline(v=which(rowMin(L0glm.osqp$coefficients)!=0), col = 'cyan')

cor(as.vector(as.vector(L0glm.osqp$coefficients)), as.vector(sim.pois.mc$beta_true)) # 0.9995664
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)!=0) # TP = 49
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)!=0) # FP = 1
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)==0) # TN = 449
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)==0) # FN = 1

# Cholesky ####
system.time(L0glm.chol <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                    family = poisson(identity), 
                                    lambda = "ebic", # tune.chol.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = TRUE,
                                    parallelize = TRUE,
                                    grouped = TRUE,
                                    solver = "eigen",
                                    method = 7L)) 
# 14s
abline(v=which(rowMin(L0glm.chol$coefficients)!=0), col = 'cyan')

R2 = round(cor(as.vector(rowMeans(sim.pois.mc$beta_true)), as.vector(rowMeans(L0glm.chol$coefficients))), digits = 4)^2 # Calculate R2
title = paste("True vs estimated coefficients\nnonnegative L0 grouped relaxed GLM\nR2=", R2)
plot(as.vector(rowMeans(L0glm.chol$coefficients)), as.vector(rowMeans(sim.pois.mc$beta_true)), pch=16, 
     xlab="estimated coefficients nonnegative L0 group GLM",
     ylab="true coefficients", col="steelblue",
     main=title)

cor(as.vector(as.vector(L0glm.chol$coefficients)), as.vector(sim.pois.mc$beta_true)) # 0.9995638
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.chol$coefficients)!=0) # TP = 49
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.chol$coefficients)!=0) # FP = 1
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.chol$coefficients)==0) # TN = 449
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.chol$coefficients)==0) # FN = 1


## 2.2. Gaussian noise ####

# TODO fix bugs that appear below
# Error in lm.fit(x, y, offset = offset, singular.ok = singular.ok, ...) : 
# 0 (non-NA) cases
# In addition: There were 17 warnings (use warnings() to see them)
# probably something going wrong with passing the dispersion

## 2.2.1 Relaxed = FALSE, grouped = FALSE ####
## Tune lambda relaxed = FALSE, grouped = FALSE
system.time(tune.osqp.mc <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, IC = "gic", 
                                       plot = FALSE, coefs.true = sim.gaus.mc$beta_true, 
                                       relaxed = FALSE, family = gaussian(identity), 
                                       nonnegative = TRUE, 
                                       parallelize = TRUE, grouped = FALSE, solver = "osqp"))
# bestlambda = XXX, time = XXXs
# TODO fix bug
system.time(tune.chol.mc <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, IC = "gic", 
                                       plot = FALSE, coefs.true = sim.gaus.mc$beta_true, 
                                       relaxed = FALSE, family = gaussian(identity), 
                                       nonnegative = TRUE, 
                                       parallelize = TRUE, grouped = FALSE,
                                       solver = "eigen", method = 7L))
# bestlambda = XXX, time = XXXs
# TODO fix bug

## Fit relaxed = FALSE, grouped = FALSE
# OSQP ####
system.time(L0glm.osqp <- L0glm.fit(sim.gaus.mc$X, sim.gaus.mc$y, 
                                    family = gaussian(identity), 
                                    lambda = "ebic", # tune.osqp.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallelize = TRUE,
                                    grouped = FALSE,
                                    solver = "osqp")) 
# 9s
# TODO fix warning
# number of rows of result is not a multiple of vector length (arg 6)

abline(v=which(rowMin(L0glm.osqp$coefficients)!=0), col = 'cyan')

cor(as.vector(as.vector(L0glm.osqp$coefficients)), as.vector(sim.gaus.mc$beta_true)) # 0.9898752
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)!=0) # TP = 50
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)!=0) # FP = 409
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)==0) # TN = 41
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)==0) # FN = 0

# Cholesky ####
system.time(L0glm.chol <- L0glm.fit(sim.gaus.mc$X, sim.gaus.mc$y, 
                                    family = gaussian(identity), 
                                    lambda = "ebic", # tune.chol.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallelize = TRUE,
                                    grouped = FALSE,
                                    solver = "eigen",
                                    method = 7L)) 
# 6.3s
# TODO fix warning
# number of rows of result is not a multiple of vector length (arg 6)

abline(v=which(rowMin(L0glm.chol$coefficients)!=0), col = 'cyan')

R2 = round(cor(as.vector(rowMeans(sim.gaus.mc$beta_true)), as.vector(rowMeans(L0glm.chol$coefficients))), digits = 4)^2 # Calculate R2
title = paste("True vs estimated coefficients\nnonnegative L0 independent GLM\nR2=", R2)
plot(as.vector(rowMeans(L0glm.chol$coefficients)), as.vector(rowMeans(sim.gaus.mc$beta_true)), pch=16, 
     xlab="estimated coefficients nonnegative L0 group GLM",
     ylab="true coefficients", col="steelblue",
     main=title)

cor(as.vector(as.vector(L0glm.chol$coefficients)), as.vector(sim.gaus.mc$beta_true)) # 0.9936296
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.chol$coefficients)!=0) # TP = 50
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.chol$coefficients)!=0) # FP = 219
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.chol$coefficients)==0) # TN = 231
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.chol$coefficients)==0) # FN = 0


## 2.2.2 Relaxed = FALSE, grouped = TRUE ####
## Tune lambda relaxed = FALSE, grouped = TRUE
system.time(tune.osqp.mc <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, IC = "gic", 
                                       plot = TRUE, coefs.true = sim.gaus.mc$beta_true, 
                                       relaxed = FALSE, family = gaussian(identity), 
                                       nonnegative = TRUE, 
                                       parallelize = TRUE, grouped = TRUE, solver = "osqp"))
# bestlambda = XXX, time = XXXs
system.time(tune.chol.mc <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, IC = "gic", 
                                       plot = FALSE, coefs.true = sim.gaus.mc$beta_true, 
                                       relaxed = FALSE, family = gaussian(identity), 
                                       nonnegative = TRUE, 
                                       parallelize = TRUE, grouped = TRUE,
                                       solver = "eigen", method = 7L))
# bestlambda = XXX, time = XXX 

## Fit relaxed = FALSE, grouped = TRUE
# OSQP ####
system.time(L0glm.osqp <- L0glm.fit(sim.gaus.mc$X, sim.gaus.mc$y, 
                                    family = gaussian(identity), 
                                    lambda = "ebic", # tune.osqp.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallelize = TRUE,
                                    grouped = TRUE,
                                    solver = "osqp")) 
# 23s
# TODO fix bug In (function (..., deparse.level = 1)  :
# number of rows of result is not a multiple of vector length (arg 6)
abline(v=which(rowMin(L0glm.osqp$coefficients)!=0), col = 'cyan')

cor(as.vector(as.vector(L0glm.osqp$coefficients)), as.vector(sim.gaus.mc$beta_true)) # 0.9999461
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)!=0) # TP = 50
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)!=0) # FP = 0
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)==0) # TN = 450
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)==0) # FN = 0

# Cholesky ####
system.time(L0glm.chol <- L0glm.fit(sim.gaus.mc$X, sim.gaus.mc$y, 
                                    family = gaussian(identity), 
                                    lambda = "ebic", # tune.chol.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallelize = TRUE,
                                    grouped = TRUE,
                                    solver = "eigen",
                                    method = 7L)) 
# 12s
# TODO fix bug In (function (..., deparse.level = 1)  :
# number of rows of result is not a multiple of vector length (arg 6)
abline(v=which(rowMin(L0glm.chol$coefficients)!=0), col = 'cyan')

R2 = round(cor(as.vector(rowMeans(sim.gaus.mc$beta_true)), as.vector(rowMeans(L0glm.chol$coefficients))), digits = 4)^2 # Calculate R2
title = paste("True vs estimated coefficients\nnonnegative L0 grouped GLM\nR2=", R2)
plot(as.vector(rowMeans(L0glm.chol$coefficients)), as.vector(rowMeans(sim.gaus.mc$beta_true)), pch=16, 
     xlab="estimated coefficients nonnegative L0 group GLM",
     ylab="true coefficients", col="steelblue",
     main=title)

cor(as.vector(as.vector(L0glm.chol$coefficients)), as.vector(sim.gaus.mc$beta_true)) # 0.9999453
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.chol$coefficients)!=0) # TP = 50
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.chol$coefficients)!=0) # FP = 0
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.chol$coefficients)==0) # TN = 450
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.chol$coefficients)==0) # FN = 0

## 2.2.3 Relaxed = TRUE, grouped = TRUE ####
## Tune lambda relaxed = TRUE, grouped = TRUE
system.time(tune.osqp.mc <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, IC = "gic", 
                                       plot = TRUE, coefs.true = sim.gaus.mc$beta_true, 
                                       relaxed = TRUE, family = gaussian(identity), 
                                       nonnegative = TRUE, 
                                       parallelize = TRUE, grouped = TRUE, solver = "osqp"))

# bestlambda = XX, time = XX
system.time(tune.chol.mc <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, IC = "gic", 
                                       plot = FALSE, coefs.true = sim.gaus.mc$beta_true, 
                                       relaxed = TRUE, family = gaussian(identity), 
                                       nonnegative = TRUE, 
                                       parallelize = TRUE, grouped = TRUE,
                                       solver = "eigen", method = 7L))
# bestlambda = XX, time = XX

## Fit relaxed = TRUE, grouped = TRUE
# OSQP ####
system.time(L0glm.osqp <- L0glm.fit(sim.gaus.mc$X, sim.gaus.mc$y, 
                                    family = gaussian(identity), 
                                    lambda = "ebic", # tune.osqp.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = TRUE,
                                    parallelize = TRUE,
                                    grouped = TRUE,
                                    solver = "osqp")) 
# 23s
# TODO fix bug In (function (..., deparse.level = 1)  :
# number of rows of result is not a multiple of vector length (arg 6)
abline(v=which(rowMin(L0glm.osqp$coefficients)!=0), col = 'cyan')

cor(as.vector(as.vector(L0glm.osqp$coefficients)), as.vector(sim.gaus.mc$beta_true)) # 0.9999514
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)!=0) # TP = 50
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)!=0) # FP = 0
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)==0) # TN = 450
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)==0) # FN = 0

# Cholesky ####
system.time(L0glm.chol <- L0glm.fit(sim.gaus.mc$X, sim.gaus.mc$y, 
                                    family = gaussian(identity), 
                                    lambda = "ebic", # tune.chol.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = TRUE,
                                    parallelize = TRUE,
                                    grouped = TRUE,
                                    solver = "eigen",
                                    method = 7L)) 
# 30.235s
# TODO fix bug In (function (..., deparse.level = 1)  :
# number of rows of result is not a multiple of vector length (arg 6)
abline(v=which(rowMin(L0glm.chol$coefficients)!=0), col = 'cyan')

R2 = round(cor(as.vector(rowMeans(sim.gaus.mc$beta_true)), as.vector(rowMeans(L0glm.chol$coefficients))), digits = 4)^2 # Calculate R2
title = paste("True vs estimated coefficients\nnonnegative L0 grouped relaxed GLM\nR2=", R2)
plot(as.vector(rowMeans(L0glm.chol$coefficients)), as.vector(rowMeans(sim.gaus.mc$beta_true)), pch=16, 
     xlab="estimated coefficients nonnegative L0 group GLM",
     ylab="true coefficients", col="steelblue",
     main=title)

cor(as.vector(as.vector(L0glm.chol$coefficients)), as.vector(sim.gaus.mc$beta_true)) # 0.999951
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.chol$coefficients)!=0) # TP = 50
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.chol$coefficients)!=0) # FP = 0
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.chol$coefficients)==0) # TN = 450
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.chol$coefficients)==0) # FN = 0


# 3. COMPARISON WITH OTHER METHODS ####
# Take our best working method and compare it against others: glmnet, ordinis, L0learn, abess.
# 3.1 SINGLE CHANNEL POISSON ####
n <- p <- 500
sim.pois <- simul(method = 2, n = n, p = p, seed = 2)

# Tune and fit with Cholesky decomposition ####
system.time(tune.chol.sc <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "gic", 
                                       plot = FALSE, coefs.true = sim.pois$beta_true, 
                                       relaxed = FALSE, family = poisson(identity), 
                                       nonnegative = TRUE, 
                                       solver = "eigen", method = 7L))
# 1s

system.time(L0glm.chol <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                    family = poisson(identity), 
                                    lambda = "ebic", # tune.chol.sc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "eigen",
                                    method = 7L)) 
# 0.01s
cor1 <- cor(as.vector(as.vector(L0glm.chol$coefficients)), as.vector(sim.pois$beta_true)) # 0.9649188
cor1
TP1 <- sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.chol$coefficients) != 0)) # TP = 42
TP1
FP1 <- sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.chol$coefficients) != 0)) # FP = 8
FP1
TN1 <- sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.chol$coefficients) == 0)) # TN = 442
TN1
FN1 <- sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.chol$coefficients) == 0)) # FN = 8
FN1
plot(x=L0glm.chol$coefficients, y=sim.pois$beta_true, pch=16, col='blue',
     main="L0 glm identity link Poisson with lambda tuned to optimise gic",
     xlab="Estimated coefficients", ylab="True coefficients")

# Nonnegative weighted lasso regression in glmnet ####
system.time(glmnet.sc.pois <- fit.GLM(sim = sim.pois, method = 2, family = poisson(identity), plot = TRUE)) # 0.158
# 0.15s

cor2 <- cor(as.vector(glmnet.sc.pois[-1]), as.vector(sim.pois$beta_true)) # 0.7715938
cor2
TP2 <- sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(glmnet.sc.pois[-1]) != 0)) # TP = 36
TP2
FP2 <- sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(glmnet.sc.pois[-1]) != 0)) # FP = 57
FP2
TN2 <- sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(glmnet.sc.pois[-1]) == 0)) # TN = 393
TN2
FN2 <- sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(glmnet.sc.pois[-1]) == 0)) # FN = 14
FN2

# Nonnegative adaptive lasso regression in glmnet ####
system.time(glmnet.al.sc.pois <- fit.GLM(sim = sim.pois, method = 3, family = poisson(identity), plot = T)) # 0.249 
# 0.23s
cor3 <- cor(as.vector(glmnet.al.sc.pois[-1]), as.vector(sim.pois$beta_true)) # 0.823307
cor3
TP3 <- sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(glmnet.al.sc.pois[-1]) != 0)) # TP = 29
TP3
FP3 <- sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(glmnet.al.sc.pois[-1]) != 0)) # FP = 32
FP3
TN3 <- sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(glmnet.al.sc.pois[-1]) == 0)) # TN = 418
TN3
FN3 <- sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(glmnet.al.sc.pois[-1]) == 0)) # FN = 21
FN3

# nonnegative L0Learn best subset ####
system.time(L0learn.sc.pois <- fit.GLM(sim = sim.pois, method = 4, 
                                       family = poisson(identity), plot = T))  
# 0.90s

cor(as.vector(L0learn.sc.pois), as.vector(sim.pois$beta_true)) # 0.5168346
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0learn.sc.pois) != 0)) # TP = 14
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0learn.sc.pois) != 0)) # FP = 20
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0learn.sc.pois) == 0)) # TN = 430
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0learn.sc.pois) == 0)) # FN = 36

# nonnegative L0Learn best subset fit after glmnet nonnegative LASSO ####
# Same results as regular L0learn but faster
system.time(L0learn.after.lasso.sc.pois <- fit.GLM(sim = sim.pois, method = 5, family = poisson(identity), plot = T)) 
# 0.43s
cor(as.vector(L0learn.after.lasso.sc.pois), as.vector(sim.pois$beta_true)) # 0.5240371
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0learn.after.lasso.sc.pois) != 0)) # TP = 16
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0learn.after.lasso.sc.pois) != 0)) # FP = 21
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0learn.after.lasso.sc.pois) == 0)) # TN = 429
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0learn.after.lasso.sc.pois) == 0)) # FN = 34

# abess best subset fit ####
# TODO: check what is wrong - getting negative coefficients, 
# which should not be allowed with provided option beta.low = 0
# Maybe raise an issue on github with the authors if this cannot be resolved?
system.time(abess.sc.pois <- fit.GLM(sim = sim.pois, method = 6, 
                                     family = poisson(identity), 
                                     plot = T, IC = "gic"))  
# 47 s
cor(as.vector(abess.sc.pois), as.vector(sim.pois$beta_true)) # 0.09423815
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(abess.sc.pois) != 0)) # TP = 28
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(abess.sc.pois) != 0)) # FP = 159
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(abess.sc.pois) == 0)) # TN = 291
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(abess.sc.pois) == 0)) # FN = 22


# abess best subset fit after nonnegative LASSO ####
system.time(abess.after.lasso.sc.pois <- fit.GLM(sim = sim.pois, method = 7, family = poisson(identity), plot = T, IC = "gic"))  
# 0.11 s
cor(as.vector(abess.after.lasso.sc.pois), as.vector(sim.pois$beta_true))  # 0.8849495
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(abess.after.lasso.sc.pois) != 0)) # TP = 29
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(abess.after.lasso.sc.pois) != 0)) # FP = 15
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(abess.after.lasso.sc.pois) == 0)) # TN = 435
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(abess.after.lasso.sc.pois) == 0)) # FN = 21

# nonnegative weighted MCP penalized regression in ordinis ####
# TODO: check - coefficients now all come out as zero
system.time(ordinis.mcp.sc.pois <- fit.GLM(sim = sim.pois, method = 8, 
                                           family = poisson(identity), plot = T, IC = "gic"))  
# 17s
cor(as.vector(ordinis.mcp.sc.pois[-1]), as.vector(sim.pois$beta_true)) # NA
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(ordinis.mcp.sc.pois[-1]) != 0)) # TP = 0
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(ordinis.mcp.sc.pois[-1]) != 0)) # FP = 0
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(ordinis.mcp.sc.pois[-1]) == 0)) # TN = 450
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(ordinis.mcp.sc.pois[-1]) == 0)) # FN = 50

# nonnegative weighted SCAD penalized regression in ordinis ####
# TODO: check - coefficients now all come out as zero
system.time(ordinis.scad.sc.pois <- fit.GLM(sim = sim.pois, method = 9, family = poisson(identity), plot = T, IC = "gic"))  
# 15s
cor(as.vector(ordinis.scad.sc.pois[-1]), as.vector(sim.pois$beta_true)) # NA
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(ordinis.scad.sc.pois[-1]) != 0)) # TP = 0
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(ordinis.scad.sc.pois[-1]) != 0)) # FP = 0
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(ordinis.scad.sc.pois[-1]) == 0)) # TN = 450
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(ordinis.scad.sc.pois[-1]) == 0)) # FN = 50


X <- sim.pois$X
y <- sim.pois$y
W <- 1/(y+1E-10)

benchmarks_1_1 = rbind(
  data.frame(method="L0glm_tuned_lambda", do.call(cbind, get.metrics(beta = L0glm.chol$coefficients,  # some metrics
                                                                     X = X, y = y, 
                                                                     family = family, 
                                                                     beta_true = sim.pois$beta_true))),
  data.frame(method="nnalasso_glmnet", do.call(cbind, get.metrics(beta = glmnet.al.sc.pois,  # some metrics
                                                                  X = cbind(1, X), y = y, 
                                                                  weights = W,
                                                                  family = gaussian(identity), 
                                                                  beta_true = c(0, as.vector(sim.pois$beta_true))) )),
  data.frame(method="nnL0Learn_after_glmnet_nnlasso", do.call(cbind, get.metrics(beta = L0learn.after.lasso.sc.pois,  # some metrics
                                                                                 X = X, y = y, 
                                                                                 weights = W,
                                                                                 family = gaussian(identity),
                                                                                 beta_true = as.vector(sim.pois$beta_true)) )),
  data.frame(method="abess", do.call(cbind, get.metrics(beta = abess.sc.pois,  # some metrics
                                                        X = X, y = y, 
                                                        weights = W,
                                                        family = gaussian(identity),
                                                        beta_true = as.vector(sim.pois$beta_true)) )),
  data.frame(method="abess_after_glmnet_nnlasso", do.call(cbind, get.metrics(beta = abess.after.lasso.sc.pois,  # some metrics
                                                                             X = X, y = y, 
                                                                             weights = W,
                                                                             family = gaussian(identity),
                                                                             beta_true = as.vector(sim.pois$beta_true)) )),
  data.frame(method="nnL0Learn", do.call(cbind, get.metrics(beta = L0learn.sc.pois,  # some metrics
                                                            X = X, y = y, 
                                                            weights = W,
                                                            family = gaussian(identity), 
                                                            beta_true = as.vector(sim.pois$beta_true)) )),
  data.frame(method="nnlasso_glmnet", do.call(cbind, get.metrics(beta = glmnet.sc.pois,  # some metrics
                                                                 X = cbind(1, X), y = y, 
                                                                 weights = W,
                                                                 family = gaussian(identity), 
                                                                 beta_true = c(0, as.vector(sim.pois$beta_true))) )),
  data.frame(method="nnMCP_ordinis", do.call(cbind, get.metrics(beta = ordinis.mcp.sc.pois,  # some metrics
                                                                X = cbind(1, X), y = y, 
                                                                weights = W,
                                                                family = gaussian(identity),
                                                                beta_true = c(0, as.vector(sim.pois$beta_true))) )),
  data.frame(method="nnSCAD_ordinis", do.call(cbind, get.metrics(beta = ordinis.scad.sc.pois,  # some metrics
                                                                 X = cbind(1, X), y = y, 
                                                                 weights = W,
                                                                 family = gaussian(identity), 
                                                                 beta_true = c(0, as.vector(sim.pois$beta_true))) )))

benchmarks_1_1 = benchmarks_1_1[order(benchmarks_1_1[,crit]),] # sort by IC
benchmarks_1_1

# 3.2 SINGLE CHANNEL GAUSSIAN ####
n <- p <- 500
sim.gaus <- simul(method = 1, n = n, p = p, seed = 2)



# Tune and fit with Cholesky decomposition ####
system.time(tune.chol.sc <- L0glm.tune(X=sim.gaus$X, y=sim.gaus$y, IC = "gic", 
                                       plot = FALSE, coefs.true = sim.gaus$beta_true, 
                                       relaxed = FALSE, family = gaussian(identity), 
                                       nonnegative = TRUE, 
                                       solver = "eigen", method = 7L))
# TODO: fix bug Error in lm.fit(x, y, offset = offset, singular.ok = singular.ok, ...) : 
# 0 (non-NA) cases

system.time(L0glm.chol <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                    family = gaussian(identity), 
                                    lambda = "ebic", # tune.chol.sc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "eigen",
                                    method = 7L)) 
# 0.01s
cor(as.vector(as.vector(L0glm.chol$coefficients)), as.vector(sim.gaus$beta_true)) # 0.9782595
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.chol$coefficients) != 0)) # TP = 44
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.chol$coefficients) != 0)) # FP = 9
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.chol$coefficients) == 0)) # TN = 441
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.chol$coefficients) == 0)) # FN = 6

plot(x=L0glm.chol$coefficients, y=sim.gaus$beta_true, pch=16, col='blue',
     main="L0 glm identity link gaussian with lambda tuned to optimise gic",
     xlab="Estimated coefficients", ylab="True coefficients")

# Nonnegative weighted lasso regression in glmnet ####
system.time(glmnet.sc.pois <- fit.GLM(sim = sim.gaus, method = 2, 
                                      family = gaussian(identity), plot = TRUE)) 
# 0.14s 
cor(as.vector(glmnet.sc.pois[-1]), as.vector(sim.gaus$beta_true)) # 0.9611292
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(glmnet.sc.pois[-1]) != 0)) # TP = 49
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(glmnet.sc.pois[-1]) != 0)) # FP = 68
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(glmnet.sc.pois[-1]) == 0)) # TN = 382
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(glmnet.sc.pois[-1]) == 0)) # FN = 1



system.time(glmnet.al.sc.pois <- fit.GLM(sim = sim.gaus, method = 3, family = gaussian(identity), plot = T)) 
# 0.3s 
cor(as.vector(glmnet.al.sc.pois[-1]), as.vector(sim.gaus$beta_true)) # 0.9626742
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(glmnet.al.sc.pois[-1]) != 0)) # TP = 37
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(glmnet.al.sc.pois[-1]) != 0)) # FP = 12
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(glmnet.al.sc.pois[-1]) == 0)) # TN = 438
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(glmnet.al.sc.pois[-1]) == 0)) # FN = 13


# nonnegative L0Learn best subset ####
system.time(L0learn.sc.pois <- fit.GLM(sim = sim.gaus, method = 4, 
                                       family = gaussian(identity), plot = T))  
#  1.02s
cor(as.vector(L0learn.sc.pois), as.vector(sim.gaus$beta_true)) # 0.7419012
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0learn.sc.pois) != 0)) # TP = 30
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0learn.sc.pois) != 0)) # FP = 38
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0learn.sc.pois) == 0)) # TN = 412
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0learn.sc.pois) == 0)) # FN = 20

# nonnegative L0Learn best subset fit after glmnet nonnegative LASSO ####
# Same results as regular L0learn but faster
system.time(L0learn.after.lasso.sc.pois <- fit.GLM(sim = sim.gaus, method = 5, 
                                                   family = gaussian(identity), plot = T))  
# 0.55s 
cor(as.vector(L0learn.after.lasso.sc.pois), as.vector(sim.gaus$beta_true)) # 0.8748054
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0learn.after.lasso.sc.pois) != 0)) # TP = 40
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0learn.after.lasso.sc.pois) != 0)) # FP = 23
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0learn.after.lasso.sc.pois) == 0)) # TN = 427
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0learn.after.lasso.sc.pois) == 0)) # FN = 10

# abess best subset fit ####
# NOTE: here nonnegativity constraints are working OK
system.time(abess.sc.pois <- fit.GLM(sim = sim.gaus, method = 6, 
                                     family = gaussian(identity), 
                                     plot = T, IC = "gic"))  
# 1.01s 
cor(as.vector(abess.sc.pois), as.vector(sim.gaus$beta_true)) # 0.8707492
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(abess.sc.pois) != 0)) # TP = 28
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(abess.sc.pois) != 0)) # FP = 19
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(abess.sc.pois) == 0)) # TN = 431
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(abess.sc.pois) == 0)) # FN = 22


# abess best subset fit after nonnegative LASSO ####
system.time(abess.after.lasso.sc.pois <- fit.GLM(sim = sim.gaus, method = 7, 
                                                 family = gaussian(identity), 
                                                 plot = T, IC = "gic"))  
# 0.15s
cor(as.vector(abess.after.lasso.sc.pois), as.vector(sim.gaus$beta_true))  # 0.8814479
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(abess.after.lasso.sc.pois) != 0)) # TP = 29
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(abess.after.lasso.sc.pois) != 0)) # FP = 19
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(abess.after.lasso.sc.pois) == 0)) # TN = 431
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(abess.after.lasso.sc.pois) == 0)) # FN = 21


# nonnegative weighted MCP penalized regression in ordinis ####
system.time(ordinis.mcp.sc.pois <- fit.GLM(sim = sim.gaus, method = 8, 
                                           family = gaussian(identity), 
                                           plot = T, IC = "gic"))  
# 0.24s
cor(as.vector(ordinis.mcp.sc.pois[-1]), as.vector(sim.gaus$beta_true)) # 0.7612053
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(ordinis.mcp.sc.pois[-1]) != 0)) # TP = 24
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(ordinis.mcp.sc.pois[-1]) != 0)) # FP = 41
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(ordinis.mcp.sc.pois[-1]) == 0)) # TN = 409
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(ordinis.mcp.sc.pois[-1]) == 0)) # FN = 26

# nonnegative weighted SCAD penalized regression in ordinis ####
system.time(ordinis.scad.sc.pois <- fit.GLM(sim = sim.gaus, method = 9, 
                                            family = gaussian(identity), 
                                            plot = T, IC = "gic"))  
# 0.27s 
cor(as.vector(ordinis.scad.sc.pois[-1]), as.vector(sim.gaus$beta_true)) # 0.9205098
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(ordinis.scad.sc.pois[-1]) != 0)) # TP = 34
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(ordinis.scad.sc.pois[-1]) != 0)) # FP = 22
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(ordinis.scad.sc.pois[-1]) == 0)) # TN = 428
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(ordinis.scad.sc.pois[-1]) == 0)) # FN = 16


X <- sim.gaus$X
y <- sim.gaus$y
W <- as.matrix(rep(1,ncol(X)))

benchmarks_1_2 = rbind(
  data.frame(method="L0glm_tuned_lambda", do.call(cbind, get.metrics(beta = L0glm.chol$coefficients,  # some metrics
                                                                     X = X, y = y, 
                                                                     family = family, 
                                                                     beta_true = sim.gaus$beta_true))),
  data.frame(method="nnalasso_glmnet", do.call(cbind, get.metrics(beta = glmnet.al.sc.pois,  # some metrics
                                                                  X = cbind(1, X), y = y, 
                                                                  weights = W,
                                                                  family = gaussian(identity), 
                                                                  beta_true = c(0, as.vector(sim.gaus$beta_true))) )),
  data.frame(method="nnL0Learn_after_glmnet_nnlasso", do.call(cbind, get.metrics(beta = L0learn.after.lasso.sc.pois,  # some metrics
                                                                                 X = X, y = y, 
                                                                                 weights = W,
                                                                                 family = gaussian(identity),
                                                                                 beta_true = as.vector(sim.gaus$beta_true)) )),
  data.frame(method="abess_after_glmnet_nnlasso", do.call(cbind, get.metrics(beta = abess.after.lasso.sc.pois,  # some metrics
                                                                             X = X, y = y, 
                                                                             weights = W,
                                                                             family = gaussian(identity),
                                                                             beta_true = as.vector(sim.gaus$beta_true)) )),
  data.frame(method="nnL0Learn", do.call(cbind, get.metrics(beta = L0learn.sc.pois,  # some metrics
                                                            X = X, y = y, 
                                                            weights = W,
                                                            family = gaussian(identity), 
                                                            beta_true = as.vector(sim.gaus$beta_true)) )),
  data.frame(method="abess", do.call(cbind, get.metrics(beta = abess.sc.pois,  # some metrics
                                                        X = X, y = y, 
                                                        weights = W,
                                                        family = gaussian(identity), 
                                                        beta_true = as.vector(sim.gaus$beta_true)) )),
  data.frame(method="nnlasso_glmnet", do.call(cbind, get.metrics(beta = glmnet.sc.pois,  # some metrics
                                                                 X = cbind(1, X), y = y, 
                                                                 weights = W,
                                                                 family = gaussian(identity), 
                                                                 beta_true = c(0, as.vector(sim.gaus$beta_true))) )),
  data.frame(method="nnMCP_ordinis", do.call(cbind, get.metrics(beta = ordinis.mcp.sc.pois,  # some metrics
                                                                X = cbind(1, X), y = y, 
                                                                weights = W,
                                                                family = gaussian(identity),
                                                                beta_true = c(0, as.vector(sim.gaus$beta_true))) )),
  data.frame(method="nnSCAD_ordinis", do.call(cbind, get.metrics(beta = ordinis.scad.sc.pois,  # some metrics
                                                                 X = cbind(1, X), y = y, 
                                                                 weights = W,
                                                                 family = gaussian(identity), 
                                                                 beta_true = c(0, as.vector(sim.gaus$beta_true))) )))

benchmarks_1_2 = benchmarks_1_2[order(benchmarks_1_2[,crit]),] # sort by IC
benchmarks_1_2



# 3.3 MULTICHANNEL POISSON ####
n <- p <- 500
sim.pois.mc <- simul(method = 4, n = n, p = p, seed = 2)

# Tune and fit with Cholesky decomposition ####
system.time(tune.chol.mc <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, IC = "gic", 
                                       plot = FALSE, coefs.true = sim.pois.mc$beta_true, 
                                       relaxed = FALSE, family = poisson(identity), 
                                       nonnegative = TRUE, 
                                       parallelize = TRUE, grouped = TRUE,
                                       solver = "eigen", method = 7L)) 
# 197s

system.time(L0glm.chol.pois <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                         family = poisson(identity), 
                                         lambda = "ebic",
                                         nonnegative = TRUE,
                                         relaxed = TRUE,
                                         parallelize = TRUE,
                                         grouped = TRUE,
                                         solver = "eigen",
                                         method = 7L)) 
# 10.3s
abline(v=which(rowMin(L0glm.chol.pois$coefficients)!=0), col = 'cyan')

R2 = round(cor(as.vector(rowMeans(sim.pois.mc$beta_true)), as.vector(rowMeans(L0glm.chol.pois$coefficients))), digits = 4)^2 # Calculate R2
title = paste("True vs estimated coefficients\nnonnegative L0 grouped GLM\nR2=", R2)
plot(as.vector(rowMeans(L0glm.chol.pois$coefficients)), as.vector(rowMeans(sim.pois.mc$beta_true)), pch=16, 
     xlab="estimated coefficients nonnegative L0 group GLM",
     ylab="true coefficients", col="steelblue",
     main=title)


image(x=1:nrow(sim.pois.mc$y), y=1:ncol(sim.pois.mc$y), z=as.matrix(L0glm.chol.pois$coefficients)^0.1, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="L0glm\n(red = true peaks, cyan = estimated peaks)")
abline(v=(1:nrow(sim.pois.mc$X))[as.vector(rowMax(sim.pois.mc$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim.pois.mc$X))[as.vector(rowMin(L0glm.chol.pois$coefficients)!=0)], col="cyan")

cor(as.vector(as.vector(L0glm.chol.pois$coefficients)), as.vector(sim.pois.mc$beta_true)) # 0.9995638
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.chol.pois$coefficients)!=0) # TP = 49
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.chol.pois$coefficients)!=0) # FP = 1
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.chol.pois$coefficients)==0) # TN = 449
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.chol.pois$coefficients)==0) # FN = 1


# best subset multitask learning using family="mgaussian" using abess ####
system.time(bsgs.abess.poisson <- fit.GLM.mc(sim = sim.pois.mc, method = 1, 
                                             family = poisson(identity), 
                                             plot = T, IC = "gic")) 
# XXX
# TODO: find bug - returns error
# Error in h(simpleError(msg, call)) : 
# error in evaluating the argument 'x' in selecting a method for function 
# 'as.matrix': error in evaluating the argument 'x' in selecting a method for 
# function 'as.matrix': incorrect number of dimensions

cor(as.vector(as.vector(sim.pois.mc$beta_true)), as.vector(bsgs.abess.poisson)) # XXX
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(bsgs.abess.poisson)!=0) # TP = XXX
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(bsgs.abess.poisson)!=0) # FP = XXX
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(bsgs.abess.poisson)==0) # TN = XXX
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(bsgs.abess.poisson)==0) # FN = XXX

# 3.4 MULTICHANNEL GAUSSIAN ####
n <- p <- 500
sim.gaus.mc <- simul(method = 3, n = n, p = p, seed = 1)

# TODO fix bug
# Error in lm.fit(x, y, offset = offset, singular.ok = singular.ok, ...) : 
# 0 (non-NA) cases
# In addition: There were 50 or more warnings (use warnings() to see the first 50)

# Tune and fit with Cholesky decomposition ####
system.time(tune.chol.mc <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, IC = "gic", 
                                       plot = FALSE, coefs.true = sim.gaus.mc$beta_true, 
                                       relaxed = FALSE, family = gaussian(identity), 
                                       nonnegative = TRUE, 
                                       parallelize = TRUE, grouped = TRUE,
                                       solver = "eigen", method = 7L)) 
# XXs
system.time(L0glm.chol.gaus <- L0glm.fit(sim.gaus.mc$X, sim.gaus.mc$y, 
                                    family = gaussian(identity), 
                                    lambda = "ebic", 
                                    nonnegative = TRUE,
                                    relaxed = TRUE,
                                    parallelize = TRUE,
                                    grouped = TRUE,
                                    solver = "eigen",
                                    method = 7L)) 
# 11s
# TODO fix warning
# Warning messages:
# 1: In (function (..., deparse.level = 1)  :
#         number of rows of result is not a multiple of vector length (arg 6)
       
# abline(v=which(rowMin(L0glm.chol.gaus$coefficients)!=0), col = 'cyan')

R2 = round(cor(as.vector(rowMeans(sim.gaus.mc$beta_true)), as.vector(rowMeans(L0glm.chol.gaus$coefficients))), digits = 4)^2 # Calculate R2
title = paste("True vs estimated coefficients\nnonnegative L0 grouped GLM\nR2=", R2)
plot(as.vector(rowMeans(L0glm.chol.gaus$coefficients)), as.vector(rowMeans(sim.gaus.mc$beta_true)), pch=16, 
     xlab="estimated coefficients nonnegative L0 group GLM",
     ylab="true coefficients", col="steelblue",
     main=title)


image(x=1:nrow(sim.gaus.mc$y), y=1:ncol(sim.gaus.mc$y), z=as.matrix(L0glm.chol.gaus$coefficients)^0.1, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="L0glm\n(red = true peaks, cyan = estimated peaks)")
abline(v=(1:nrow(sim.gaus.mc$X))[as.vector(rowMax(sim.gaus.mc$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim.gaus.mc$X))[as.vector(rowMin(L0glm.chol.gaus$coefficients)!=0)], col="cyan")

cor(as.vector(as.vector(L0glm.chol.gaus$coefficients)), as.vector(sim.gaus.mc$beta_true)) # 0.999951
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.chol.gaus$coefficients)!=0) # TP = 50
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.chol.gaus$coefficients)!=0) # FP = 0
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.chol.gaus$coefficients)==0) # TN = 450
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.chol.gaus$coefficients)==0) # FN = 0


# best subset multitask learning using family="mgaussian" using abess ####
system.time(bsgs.abess.gaus <- fit.GLM.mc(sim = sim.gaus.mc, method = 1, 
                                          family = gaussian(identity), plot = T, 
                                          IC = "gic")) 
# XXs
# TODO: fix bug 
# Error in h(simpleError(msg, call)) : 
# error in evaluating the argument 'x' in selecting a method for 
# function 'as.matrix': error in evaluating the argument 'x' in selecting a 
# method for function 'as.matrix': incorrect number of dimensions

cor(as.vector(as.vector(sim.gaus.mc$beta_true)), as.vector(bsgs.abess.gaus)) # XX
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(bsgs.abess.gaus)!=0) # TP = XX
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(bsgs.abess.gaus)!=0) # FP = XX
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(bsgs.abess.gaus)==0) # TN = XX
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(bsgs.abess.gaus)==0) # FN = XX


# glmnet nonnegative mgaussian multivariate group LASSO fit
system.time(glmnet.mc <- fit.GLM.mc(sim = sim.gaus.mc, method = 2, 
                                    family = gaussian(identity), plot = T, IC = "gic"))
# TODO fix, outputs wrong dimensions
# bug Error in fit.GLM.mc(sim = sim.gaus.mc, method = 2, family = gaussian(identity),  : 
# object 'm' not found
# XXs

image(x=1:nrow(sim$y), y=1:ncol(sim.gaus.mc$y), z=as.matrix(glmnet.mc)^0.1, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="L0glm\n(red = true peaks, cyan = estimated peaks)")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim.gaus.mc$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim$X))[as.vector(rowMin(glmnet.mc)!=0)], col="cyan")


cor(as.vector(as.vector(sim.gaus.mc$beta_true)), as.vector(glmnet.mc)) #
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(glmnet.mc)!=0) #
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(glmnet.mc)!=0) #
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(glmnet.mc)==0) #
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(glmnet.mc)==0) #


benchmarks_1_3 = rbind(
  data.frame(method="Grouped_L0glm_poisson", do.call(cbind, get.metrics(beta = L0glm.chol.pois$coefficients,  # some metrics
                                                                        X = sim.pois.mc$X, y = sim.pois.mc$y, 
                                                                        family = poisson(identity), 
                                                                        beta_true = sim.pois.mc$beta_true))),
  data.frame(method="BSGS_abess_poisson", do.call(cbind, get.metrics(beta = bsgs.abess.poisson,  # some metrics
                                                                     X = sim.pois.mc$X, y = sim.pois.mc$y, 
                                                                     family = gaussian(identity), 
                                                                     beta_true = sim.pois.mc$beta_true))),
  data.frame(method="Grouped_L0glm_gaussian", do.call(cbind, get.metrics(beta = L0glm.chol.gaus$coefficients,  # some metrics
                                                                         X = sim.gaus.mc$X, y = sim.gaus.mc$y, 
                                                                         weights = W,
                                                                         family = gaussian(identity),
                                                                         beta_true = sim.gaus.mc$beta_true))),
  data.frame(method="BSGS_abess_gaussian", do.call(cbind, get.metrics(beta = bsgs.abess.gaus,  # some metrics
                                                                      X = sim.gaus.mc$X, y = sim.gaus.mc$y, 
                                                                      family = gaussian(identity),
                                                                      beta_true = sim.gaus.mc$beta_true))))

benchmarks_1_3
