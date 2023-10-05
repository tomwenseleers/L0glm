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
# 0.16s
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
# 0.18s
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
# 0.03s
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
# Take our best working solver & method and compare it against others: 
# glmnet (LASSO), ordinis (MCP & SCAD), L0Learn & abess

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
# 0.03s
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
system.time(glmnet.sc.pois <- fit.GLM(sim = sim.pois, method = 2, 
                                      family = poisson(identity), plot = TRUE)) # 0.158
# 0.29s

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
system.time(glmnet.al.sc.pois <- fit.GLM(sim = sim.pois, method = 3, 
                                         family = poisson(identity), plot = T)) # 0.249 
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
# Slightly better than regular L0Learn & faster
system.time(L0learn.after.lasso.sc.pois <- fit.GLM(sim = sim.pois, method = 5, family = poisson(identity), plot = T)) 
# 0.43s
cor(as.vector(L0learn.after.lasso.sc.pois), as.vector(sim.pois$beta_true)) # 0.5240371
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0learn.after.lasso.sc.pois) != 0)) # TP = 16
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0learn.after.lasso.sc.pois) != 0)) # FP = 21
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0learn.after.lasso.sc.pois) == 0)) # TN = 429
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0learn.after.lasso.sc.pois) == 0)) # FN = 34

# abess best subset fit ####
system.time(abess.sc.pois <- fit.GLM(sim = sim.pois, method = 6, 
                                     family = poisson(identity), 
                                     plot = T, IC = "ebic"))  
# 1.0s
cor(as.vector(abess.sc.pois), as.vector(sim.pois$beta_true)) # 0.8960064
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(abess.sc.pois) != 0)) # TP = 29
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(abess.sc.pois) != 0)) # FP = 16
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(abess.sc.pois) == 0)) # TN = 434
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(abess.sc.pois) == 0)) # FN = 21


# abess best subset fit after nonnegative LASSO ####
system.time(abess.after.lasso.sc.pois <- fit.GLM(sim = sim.pois, method = 7, 
                                                 family = poisson(identity), plot = T, 
                                                 IC = "ebic"))  
# 0.14s
cor(as.vector(abess.after.lasso.sc.pois), as.vector(sim.pois$beta_true))  # 0.8849495
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(abess.after.lasso.sc.pois) != 0)) # TP = 29
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(abess.after.lasso.sc.pois) != 0)) # FP = 15
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(abess.after.lasso.sc.pois) == 0)) # TN = 435
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(abess.after.lasso.sc.pois) == 0)) # FN = 21

# nonnegative weighted MCP penalized regression in ordinis ####
system.time(ordinis.mcp.sc.pois <- fit.GLM(sim = sim.pois, method = 8, 
                                           family = poisson(identity), plot = T,
                                           IC = "ebic"))  
# 17s
cor(as.vector(ordinis.mcp.sc.pois[-1]), as.vector(sim.pois$beta_true)) # 0.4934751
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(ordinis.mcp.sc.pois[-1]) != 0)) # TP = 31
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(ordinis.mcp.sc.pois[-1]) != 0)) # FP = 74
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(ordinis.mcp.sc.pois[-1]) == 0)) # TN = 376
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(ordinis.mcp.sc.pois[-1]) == 0)) # FN = 19

# nonnegative weighted SCAD penalized regression in ordinis ####
system.time(ordinis.scad.sc.pois <- fit.GLM(sim = sim.pois, method = 9, 
                                            family = poisson(identity), 
                                            plot = T, IC = "ebic"))  
# 15s
cor(as.vector(ordinis.scad.sc.pois[-1]), as.vector(sim.pois$beta_true)) # 0.3871479
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(ordinis.scad.sc.pois[-1]) != 0)) # TP = 34
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(ordinis.scad.sc.pois[-1]) != 0)) # FP = 140
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(ordinis.scad.sc.pois[-1]) == 0)) # TN = 310
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(ordinis.scad.sc.pois[-1]) == 0)) # FN = 16


X <- sim.pois$X
y <- sim.pois$y
W <- 1/(y+1E-10)
family <- poisson(identity)
IC <- "ebic"

benchmarks_1_1 = rbind(
  data.frame(method="L0glm", do.call(cbind, get.metrics(beta = L0glm.chol$coefficients,  # some metrics
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
  data.frame(method="nnabess", do.call(cbind, get.metrics(beta = abess.sc.pois,  # some metrics
                                                        X = X, y = y, 
                                                        weights = W,
                                                        family = gaussian(identity),
                                                        beta_true = as.vector(sim.pois$beta_true)) )),
  data.frame(method="nnabess_after_glmnet_nnlasso", do.call(cbind, get.metrics(beta = abess.after.lasso.sc.pois,  # some metrics
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

benchmarks_1_1 = benchmarks_1_1[order(benchmarks_1_1[,"deviance_true"]),] # sort by deviance relative to ground truth
benchmarks_1_1
#                           method   p   n m     deviance deviance_true    loglik    MD_coefs
# 1                          L0glm 500 500 1 4.438390e+02  1.000000e+10 -1691.096   0.4454946
# 5   nnabess_after_glmnet_nnlasso 501 500 1 5.574415e+09  3.386834e+10 -5291.828  -1.6789370
# 8                  nnMCP_ordinis 502 500 1 1.340793e+05  5.646967e+10 -2693.012 -13.2227345
# 7                 nnlasso_glmnet 502 500 1 3.674281e+04  5.647872e+10 -2357.389 -11.7842135
# 2                nnalasso_glmnet 502 500 1 3.687858e+04  5.647888e+10 -2326.312 -11.4425444
# 3 nnL0Learn_after_glmnet_nnlasso 501 500 1 4.899196e+04  5.647890e+10 -2374.318 -15.1283751
# 6                      nnL0Learn 501 500 1 5.950925e+04  5.647891e+10 -2419.937 -18.4846836
# 9                 nnSCAD_ordinis 502 500 1 1.000234e+05  5.647895e+10 -2688.755 -47.7944838
# 4                        nnabess 501 500 1 2.085163e+12  1.684732e+12 -6773.932   1.3860283
#   MAD_coefs sMAPE_coefs bias_coefs  R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs
# 1  64.41875  0.06789582  0.9831815 0.9343798   207.51031     441.1321      0.3589582
# 5 181.21174  0.15463477  0.8764853 0.7927635   939.74340    1100.4021      0.8026678
# 8 418.02641  0.41876542  0.8094662 0.2773025   310.20936    1621.2825      1.6362045
# 7 262.54815  0.32664096  0.8779815 0.6136579   196.19003    1154.6996      1.4585712
# 2 223.39696  0.23222499  0.8760065 0.6921514   633.52652    1107.2794      1.2187659
# 3 379.78663  0.22925527  0.5287055 0.2972245  2443.80208    2484.0122      1.2602063
# 6 396.71492  0.23377023  0.5213261 0.2894750  2774.58026    2787.6291      1.3201505
# 9 464.67008  0.69166974  1.0847489 0.1905280    75.58191    1124.7323      1.8076238
# 4 166.91188  0.15662347  0.8552278 0.8113303   928.71875     941.2580      0.8069275
#   bias_selcoefs R2_selcoefs support_size_true support_size_estimated
# 1     0.9831815   0.9494126                50                     50
# 5     0.8764853   0.8887902                50                     44
# 8     0.8094662   0.3186301                50                    105
# 7     0.8779815   0.6468717                50                     93
# 2     0.8760065   0.8090143                50                     61
# 3     0.5287055   0.5675801                50                     37
# 6     0.5213261   0.5680066                50                     34
# 9     1.0847489   0.2201641                50                    174
# 4     0.8552278   0.9070003                50                     45
#   estimated_over_true_support_size TP  FP  TN FN precision recall specificity  accuracy
# 1                             1.00 42   8 442  8 0.8400000   0.84   0.9822222 0.9680000
# 5                             0.88 29  15 435 21 0.6590909   0.58   0.9666667 0.9280000
# 8                             2.10 31  74 377 19 0.2952381   0.62   0.8359202 0.8143713
# 7                             1.86 36  57 394 14 0.3870968   0.72   0.8736142 0.8582834
# 2                             1.22 29  32 419 21 0.4754098   0.58   0.9290466 0.8942116
# 3                             0.74 16  21 429 34 0.4324324   0.32   0.9533333 0.8900000
# 6                             0.68 14  20 430 36 0.4117647   0.28   0.9555556 0.8880000
# 9                             3.48 34 140 311 16 0.1954023   0.68   0.6895787 0.6886228
# 4                             0.90 29  16 434 21 0.6444444   0.58   0.9644444 0.9260000
#    F1_score       gic      ebic       bic       aic      mbic
# 1 0.8400000  3949.866  4003.652  3692.922  3482.191  4175.753
# 5 0.6170213 11094.728 11143.151 10863.314 10673.657 11298.042
# 8 0.4000000  6590.266  6704.367  6044.772  5598.023  7069.221
# 7 0.5034965  5782.692  5883.876  5298.952  4902.779  6207.426
# 2 0.5225225  5356.991  5423.729  5037.929  4776.623  5637.135
# 3 0.3678161  5180.207  5221.098  4984.791  4824.636  5351.894
# 6 0.3333333  5237.373  5275.036  5057.385  4909.873  5395.507
# 9 0.3035714  7365.646  7554.020  6465.066  5727.510  8156.373
# 4 0.6105263 14070.292 14119.792 13833.736 13639.864 14278.125

# TODO: should also include timings (but L0glm is much faster than all the rest)


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
# 0.03s
cor(as.vector(as.vector(L0glm.chol$coefficients)), as.vector(sim.gaus$beta_true)) # 0.9782595
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.chol$coefficients) != 0)) # TP = 44
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.chol$coefficients) != 0)) # FP = 9
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.chol$coefficients) == 0)) # TN = 441
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.chol$coefficients) == 0)) # FN = 6

plot(x=L0glm.chol$coefficients, y=sim.gaus$beta_true, pch=16, col='blue',
     main="L0 glm identity link gaussian with lambda tuned to optimise gic",
     xlab="Estimated coefficients", ylab="True coefficients")

# Nonnegative weighted lasso regression in glmnet ####
system.time(glmnet.sc.gaus <- fit.GLM(sim = sim.gaus, method = 2, 
                                      family = gaussian(identity), plot = TRUE)) 
# 0.19s 
cor(as.vector(glmnet.sc.gaus[-1]), as.vector(sim.gaus$beta_true)) # 0.9611292
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(glmnet.sc.gaus[-1]) != 0)) # TP = 49
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(glmnet.sc.gaus[-1]) != 0)) # FP = 68
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(glmnet.sc.gaus[-1]) == 0)) # TN = 382
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(glmnet.sc.gaus[-1]) == 0)) # FN = 1


system.time(glmnet.al.sc.gaus <- fit.GLM(sim = sim.gaus, method = 3, 
                                         family = gaussian(identity), plot = T)) 
# 0.22s 
cor(as.vector(glmnet.al.sc.gaus[-1]), as.vector(sim.gaus$beta_true)) # 0.9626742
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(glmnet.al.sc.gaus[-1]) != 0)) # TP = 37
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(glmnet.al.sc.gaus[-1]) != 0)) # FP = 12
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(glmnet.al.sc.gaus[-1]) == 0)) # TN = 438
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(glmnet.al.sc.gaus[-1]) == 0)) # FN = 13


# nonnegative L0Learn best subset ####
system.time(L0learn.sc.gaus <- fit.GLM(sim = sim.gaus, method = 4, 
                                       family = gaussian(identity), plot = T))  
#  1.02s
cor(as.vector(L0learn.sc.gaus), as.vector(sim.gaus$beta_true)) # 0.7419012
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0learn.sc.gaus) != 0)) # TP = 30
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0learn.sc.gaus) != 0)) # FP = 38
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0learn.sc.gaus) == 0)) # TN = 412
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0learn.sc.gaus) == 0)) # FN = 20

# nonnegative L0Learn best subset fit after glmnet nonnegative LASSO ####
# Better than regular L0Learn & faster
system.time(L0learn.after.lasso.sc.gaus <- fit.GLM(sim = sim.gaus, method = 5, 
                                                   family = gaussian(identity), plot = T))  
# 0.55s 
cor(as.vector(L0learn.after.lasso.sc.pois), as.vector(sim.gaus$beta_true)) # 0.8748054
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0learn.after.lasso.sc.pois) != 0)) # TP = 40
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0learn.after.lasso.sc.pois) != 0)) # FP = 23
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0learn.after.lasso.sc.pois) == 0)) # TN = 427
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0learn.after.lasso.sc.pois) == 0)) # FN = 10

# nonnegative abess best subset fit ####
system.time(abess.sc.gaus <- fit.GLM(sim = sim.gaus, method = 6, 
                                     family = gaussian(identity), 
                                     plot = T, IC = "gic"))  
# 1.01s 
cor(as.vector(abess.sc.gaus), as.vector(sim.gaus$beta_true)) # 0.8707492
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(abess.sc.gaus) != 0)) # TP = 28
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(abess.sc.gaus) != 0)) # FP = 20
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(abess.sc.gaus) == 0)) # TN = 430
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(abess.sc.gaus) == 0)) # FN = 22


# nonnegative abess best subset fit after nonnegative LASSO ####
system.time(abess.after.lasso.sc.gaus <- fit.GLM(sim = sim.gaus, method = 7, 
                                                 family = gaussian(identity), 
                                                 plot = T, IC = "gic"))  
# 0.19s
cor(as.vector(abess.after.lasso.sc.gaus), as.vector(sim.gaus$beta_true))  # 0.8814479
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(abess.after.lasso.sc.gaus) != 0)) # TP = 29
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(abess.after.lasso.sc.gaus) != 0)) # FP = 19
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(abess.after.lasso.sc.gaus) == 0)) # TN = 431
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(abess.after.lasso.sc.gaus) == 0)) # FN = 21


# nonnegative weighted MCP penalized regression in ordinis ####
system.time(ordinis.mcp.sc.gaus <- fit.GLM(sim = sim.gaus, method = 8, 
                                           family = gaussian(identity), 
                                           plot = T, IC = "ebic"))  
# 0.30s
cor(as.vector(ordinis.mcp.sc.gaus[-1]), as.vector(sim.gaus$beta_true)) # 0.7612053
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(ordinis.mcp.sc.gaus[-1]) != 0)) # TP = 24
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(ordinis.mcp.sc.gaus[-1]) != 0)) # FP = 41
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(ordinis.mcp.sc.gaus[-1]) == 0)) # TN = 409
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(ordinis.mcp.sc.gaus[-1]) == 0)) # FN = 26

# nonnegative weighted SCAD penalized regression in ordinis ####
system.time(ordinis.scad.sc.gaus <- fit.GLM(sim = sim.gaus, method = 9, 
                                            family = gaussian(identity), 
                                            plot = T, IC = "ebic"))  
# 0.33s 
cor(as.vector(ordinis.scad.sc.gaus[-1]), as.vector(sim.gaus$beta_true)) # 0.9205098
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(ordinis.scad.sc.gaus[-1]) != 0)) # TP = 34
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(ordinis.scad.sc.gaus[-1]) != 0)) # FP = 22
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(ordinis.scad.sc.gaus[-1]) == 0)) # TN = 428
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(ordinis.scad.sc.gaus[-1]) == 0)) # FN = 16

X <- sim.gaus$X
y <- sim.gaus$y
W <- as.matrix(rep(1,ncol(X)))
family <- gaussian(identity)
IC <- "ebic"

benchmarks_1_2 = rbind(
  data.frame(method="L0glm", do.call(cbind, get.metrics(beta = L0glm.chol$coefficients,  # some metrics
                                                                     X = X, y = y, 
                                                                     family = family, 
                                                                     beta_true = sim.gaus$beta_true))),
  data.frame(method="nnalasso_glmnet", do.call(cbind, get.metrics(beta = glmnet.al.sc.gaus,  # some metrics
                                                                  X = cbind(1, X), y = y, 
                                                                  weights = W,
                                                                  family = gaussian(identity), 
                                                                  beta_true = c(0, as.vector(sim.gaus$beta_true))) )),
  data.frame(method="nnL0Learn_after_glmnet_nnlasso", do.call(cbind, get.metrics(beta = L0learn.after.lasso.sc.gaus,  # some metrics
                                                                                 X = X, y = y, 
                                                                                 weights = W,
                                                                                 family = gaussian(identity),
                                                                                 beta_true = as.vector(sim.gaus$beta_true)) )),
  data.frame(method="abess_after_glmnet_nnlasso", do.call(cbind, get.metrics(beta = abess.after.lasso.sc.gaus,  # some metrics
                                                                             X = X, y = y, 
                                                                             weights = W,
                                                                             family = gaussian(identity),
                                                                             beta_true = as.vector(sim.gaus$beta_true)) )),
  data.frame(method="nnL0Learn", do.call(cbind, get.metrics(beta = L0learn.sc.gaus,  # some metrics
                                                            X = X, y = y, 
                                                            weights = W,
                                                            family = gaussian(identity), 
                                                            beta_true = as.vector(sim.gaus$beta_true)) )),
  data.frame(method="abess", do.call(cbind, get.metrics(beta = abess.sc.gaus,  # some metrics
                                                        X = X, y = y, 
                                                        weights = W,
                                                        family = gaussian(identity), 
                                                        beta_true = as.vector(sim.gaus$beta_true)) )),
  data.frame(method="nnlasso_glmnet", do.call(cbind, get.metrics(beta = glmnet.sc.gaus,  # some metrics
                                                                 X = cbind(1, X), y = y, 
                                                                 weights = W,
                                                                 family = gaussian(identity), 
                                                                 beta_true = c(0, as.vector(sim.gaus$beta_true))) )),
  data.frame(method="nnMCP_ordinis", do.call(cbind, get.metrics(beta = ordinis.mcp.sc.gaus,  # some metrics
                                                                X = cbind(1, X), y = y, 
                                                                weights = W,
                                                                family = gaussian(identity),
                                                                beta_true = c(0, as.vector(sim.gaus$beta_true))) )),
  data.frame(method="nnSCAD_ordinis", do.call(cbind, get.metrics(beta = ordinis.scad.sc.gaus,  # some metrics
                                                                 X = cbind(1, X), y = y, 
                                                                 weights = W,
                                                                 family = gaussian(identity), 
                                                                 beta_true = c(0, as.vector(sim.gaus$beta_true))) )))

benchmarks_1_2 = benchmarks_1_2[order(benchmarks_1_2[,"deviance_true"]),] # sort by deviance relative to ground truth
benchmarks_1_2
#                           method   p   n m   deviance deviance_true    loglik     MD_coefs
# 1                          L0glm 501 500 1   60715.87      54193.76 -1462.307   0.02734418
# 3 nnL0Learn_after_glmnet_nnlasso 501 500 1  183682.11     169463.44 -1749.058   0.29231018
# 9                 nnSCAD_ordinis 502 500 1  508208.93     491714.23 -1995.479   0.60143749
# 5                      nnL0Learn 501 500 1  915515.42     892280.13 -2155.628   0.97776208
# 8                  nnMCP_ordinis 502 500 1 1273567.56    1249035.32 -2234.150   1.12336285
# 4     abess_after_glmnet_nnlasso 501 500 1 3258284.46    3269706.78 -2452.995   0.94014763
# 2                nnalasso_glmnet 502 500 1 3665597.68    3667508.25 -2482.443  -6.09600492
# 7                 nnlasso_glmnet 502 500 1 4338167.47    4350827.91 -2592.558 -11.60169096
# 6                          abess 501 500 1 7116162.10    7116773.33 -2648.287   2.07653849
#   MAD_coefs sMAPE_coefs bias_coefs  R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs
# 1  45.79231  0.06576185  0.9736408 0.9590305   122.13727     310.1236      0.3939798
# 3 159.66279  0.15079786  0.9322298 0.7763735   350.89061     918.5943      0.8793481
# 9 131.55853  0.15941428  0.9207443 0.8543947   576.38337     605.9764      0.8547599
# 5 253.69498  0.24459843  0.7897581 0.5698949   844.41919    1028.1745      1.2102826
# 8 260.28405  0.27135676  0.8101574 0.5979304   959.01220    1055.8357      1.2915344
# 4 176.14620  0.16691999  0.8598522 0.7867479   916.22988     928.4196      0.8637499
# 2  95.93336  0.11343157  0.9764912 0.9301852   223.27211     695.2690      0.6291677
# 7 125.95294  0.33233782  1.0491874 0.9274381   -39.80944     529.4674      1.4059936
# 6 186.68587  0.17576504  0.8381777 0.7684228   979.49260     986.7825      0.9142191
#   bias_selcoefs R2_selcoefs support_size_true support_size_estimated
# 1     0.9736408   0.9658211                50                     53
# 3     0.9322298   0.8213999                50                     63
# 9     0.9207443   0.9355237                50                     56
# 5     0.7897581   0.7753111                50                     68
# 8     0.8101574   0.8073822                50                     65
# 4     0.8598522   0.8957397                50                     48
# 2     0.9764912   0.9447399                50                     49
# 7     1.0491874   0.9285669                50                    117
# 6     0.8381777   0.8803839                50                     48
#   estimated_over_true_support_size TP FP  TN FN precision recall specificity  accuracy
# 1                             1.06 44  9 441  6 0.8301887   0.88   0.9800000 0.9700000
# 3                             1.26 40 23 427 10 0.6349206   0.80   0.9488889 0.9340000
# 9                             1.12 34 22 429 16 0.6071429   0.68   0.9512195 0.9241517
# 5                             1.36 30 38 412 20 0.4411765   0.60   0.9155556 0.8840000
# 8                             1.30 24 41 410 26 0.3692308   0.48   0.9090909 0.8662675
# 4                             0.96 29 19 431 21 0.6041667   0.58   0.9577778 0.9200000
# 2                             0.98 37 12 439 13 0.7551020   0.74   0.9733925 0.9500998
# 7                             2.34 49 68 383  1 0.4188034   0.98   0.8492239 0.8622754
# 6                             0.96 28 20 430 22 0.5833333   0.56   0.9555556 0.9160000
#    F1_score      gic     ebic      bic      aic     mbic
# 1 0.8543689 3537.900 3596.008 3260.204 3032.615 3781.877
# 3 0.7079646 4224.972 4293.841 3895.850 3626.115 4514.130
# 9 0.6415094 4638.523 4699.879 4345.191 4104.958 4896.074
# 5 0.5084746 5094.898 5169.147 4740.064 4449.256 5406.647
# 8 0.4173913 5218.112 5289.156 4878.465 4600.301 5516.329
# 4 0.5918367 5462.490 5515.218 5210.506 5003.990 5683.877
# 2 0.7474747 5532.924 5586.745 5275.616 5064.885 5758.846
# 7 0.5868263 6525.687 6652.705 5918.440 5421.116 7058.864
# 6 0.5714286 5853.074 5905.801 5601.090 5394.574 6074.460

# TODO: should also include timings (but L0glm is much faster than all the rest)


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
# 13.7s

R2 = round(cor(as.vector(rowMeans(sim.pois.mc$beta_true)), as.vector(rowMeans(L0glm.chol.pois$coefficients))), digits = 4)^2 # Calculate R2
title = paste("True vs estimated coefficients\nnonnegative relaxed L0 grouped GLM\nR2=", R2)
plot(as.vector(rowMeans(L0glm.chol.pois$coefficients)), as.vector(rowMeans(sim.pois.mc$beta_true)), pch=16, 
     xlab="estimated coefficients nonnegative relaxed L0 group GLM",
     ylab="true coefficients", col="steelblue",
     main=title)


image(x=1:nrow(sim.pois.mc$y), y=1:ncol(sim.pois.mc$y), z=as.matrix(L0glm.chol.pois$coefficients)^0.1, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", 
      main="nonnegative relaxed L0 group GLM\n(red = true peaks, cyan = estimated peaks)")
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
                                             plot = T, IC = "ebic")) 
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
