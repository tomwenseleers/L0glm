# 1. SINGLE CHANNEL####
# Simulate SC n=p=500
n <- p <- 500
sim.gaus <- simulate(method = 1, n = n, p = p, sparse = F)
sim.pois <- simulate(method = 2, n = n, p = p)

## 1.1. Poisson####
## Tune lambda relaxed = F with AIC
system.time(tune.osqp.sc <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "aic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = FALSE, family = poisson(identity), nonnegative = TRUE, warmupiters = 1, solver = "osqp"))
# bestlambda = 1, time = 4.614s
system.time(tune.glmnet.sc <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "aic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = FALSE, family = poisson(identity), nonnegative = TRUE, warmupiters = 1, solver = "glmnet"))
# bestlambda = 0.3162, time = 5.411s
system.time(tune.chol.sc <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "aic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = FALSE, family = poisson(identity), nonnegative = TRUE, warmupiters = 1, solver = "eigen", method = 7L))
# bestlambda = 1, time = 3.052s
system.time(tune.QR.sc <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "aic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = FALSE, family = poisson(identity), nonnegative = TRUE, warmupiters = 1, solver = "eigen", method = 1L))
# bestlambda = 1, time = 15.925s
system.time(tune.CG.sc <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "bic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = FALSE, family = poisson(identity), nonnegative = TRUE, warmupiters = 1, solver = "eigen", method = 5L))
# bestlambda = 10, time = 7.951s
system.time(tune.LSCG.sc <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "bic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = FALSE, family = poisson(identity), nonnegative = TRUE, warmupiters = 1, solver = "eigen", method = 3L))
# bestlambda = 10, time = 10.929s

## 1.1.1 Fit relaxed = FALSE####
# OSQP####
system.time(L0glm.osqp <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                    family = poisson(identity), 
                                    lambda = tune.osqp.sc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "osqp")) 
# 0.125s
cor(as.vector(as.vector(L0glm.osqp$coefficients)), as.vector(sim.pois$beta_true)) # 0.9804646
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.osqp$coefficients) != 0)) # TP = 44
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.osqp$coefficients) != 0)) # FP = 40
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.osqp$coefficients) == 0)) # TN = 410
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.osqp$coefficients) == 0)) # FN = 6

# GLMNET####
system.time(L0glm.glmnet <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                      family = poisson(identity), 
                                      lambda = tune.glmnet.sc$best.lambda,
                                      nonnegative = TRUE,
                                      relaxed = FALSE,
                                      solver = "glmnet")) 
# 0.268s
cor(as.vector(as.vector(L0glm.glmnet$coefficients)), as.vector(sim.pois$beta_true)) # 0.9866356
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.glmnet$coefficients) != 0)) # TP = 45
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.glmnet$coefficients) != 0)) # FP = 21
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.glmnet$coefficients) == 0)) # TN = 429
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.glmnetp$coefficients) == 0)) # FN = 5

# QR####
system.time(L0glm.QR <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                  family = poisson(identity), 
                                  lambda = tune.QR.sc$best.lambda,
                                  nonnegative = TRUE,
                                  relaxed = FALSE,
                                  solver = "eigen",
                                  method = 1L)) 
# 0.83s
cor(as.vector(as.vector(L0glm.QR$coefficients)), as.vector(sim.pois$beta_true)) # 0.9921868
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.QR$coefficients) != 0)) # TP = 45
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.QR$coefficients) != 0)) # FP = 11
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.QR$coefficients) == 0)) # TN = 439
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.osqp$coefficients) == 0)) # FN = 6

# Chol####
system.time(L0glm.chol <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                    family = poisson(identity), 
                                    lambda = tune.chol.sc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "eigen",
                                    method = 7L)) 
# 0.122s
cor(as.vector(as.vector(L0glm.chol$coefficients)), as.vector(sim.pois$beta_true)) # 0.9921868
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.chol$coefficients) != 0)) # TP = 45
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.chol$coefficients) != 0)) # FP = 11
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.chol$coefficients) == 0)) # TN = 439
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.chol$coefficients) == 0)) # FN = 5

# LSCG####
system.time(L0glm.LSCG <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                    family = poisson(identity), 
                                    lambda = tune.LSCG.sc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "eigen",
                                    method = 3L)) 
# 0.536s
cor(as.vector(as.vector(L0glm.LSCG$coefficients)), as.vector(sim.pois$beta_true)) # 0.9308438
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.LSCG$coefficients) != 0)) # TP = 42
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.LSCG$coefficients) != 0)) # FP = 9
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.LSCG$coefficients) == 0)) # TN = 441
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.LSCG$coefficients) == 0)) # FN = 8

# CG####
system.time(L0glm.CG <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                  family = poisson(identity), 
                                  lambda = tune.CG.sc$best.lambda,
                                  nonnegative = TRUE,
                                  relaxed = TRUE,
                                  solver = "eigen",
                                  method = 5L)) 
# 0.384s
cor(as.vector(as.vector(L0glm.CG$coefficients)), as.vector(sim.pois$beta_true)) # 0.938986
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.CG$coefficients) != 0)) # TP = 42
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.CG$coefficients) != 0)) # FP = 9
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.CG$coefficients) == 0)) # TN = 441
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.CG$coefficients) == 0)) # FN = 8

## Tune lambda relaxed = T with AIC
system.time(tune.osqp.sc <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "aic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = TRUE, family = poisson(identity), nonnegative = TRUE, warmupiters = 1, solver = "osqp"))
# bestlambda = 1, time = 4.614s
system.time(tune.glmnet.sc <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "aic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = TRUE, family = poisson(identity), nonnegative = TRUE, warmupiters = 1, solver = "glmnet"))
# bestlambda = 1, time = 5.411s
system.time(tune.chol.sc <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "aic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = TRUE, family = poisson(identity), nonnegative = TRUE, warmupiters = 1, solver = "eigen", method = 7L))
# bestlambda = 1, time = 3.052s
system.time(tune.QR.sc <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "aic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = TRUE, family = poisson(identity), nonnegative = TRUE, warmupiters = 1, solver = "eigen", method = 1L))
# bestlambda = 1, time = 15.925s
system.time(tune.CG.sc <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "bic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = TRUE, family = poisson(identity), nonnegative = TRUE, warmupiters = 1, solver = "eigen", method = 5L))
# bestlambda = 10, time = 7.951s
system.time(tune.LSCG.sc <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "bic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = TRUE, family = poisson(identity), nonnegative = TRUE, warmupiters = 1, solver = "eigen", method = 3L))
# bestlambda = 10, time = 10.929s

## 1.1.2 Fit relaxed = T####
# OSQP####
system.time(L0glm.osqp <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                    family = poisson(identity), 
                                    lambda = tune.osqp.sc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "osqp")) 
# 0.125s
cor(as.vector(as.vector(L0glm.osqp$coefficients)), as.vector(sim.pois$beta_true)) # 0.9921868
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.osqp$coefficients) != 0)) # TP = 45
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.osqp$coefficients) != 0)) # FP = 12
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.osqp$coefficients) == 0)) # TN = 438
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.osqp$coefficients) == 0)) # FN = 5

# GLMNET####
system.time(L0glm.glmnet <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                      family = poisson(identity), 
                                      lambda = tune.glmnet.sc$best.lambda,
                                      nonnegative = TRUE,
                                      relaxed = FALSE,
                                      solver = "glmnet")) 
# 0.268s
cor(as.vector(as.vector(L0glm.glmnet$coefficients)), as.vector(sim.pois$beta_true)) # 0.9921868
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.glmnet$coefficients) != 0)) # TP = 45
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.glmnet$coefficients) != 0)) # FP = 12
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.glmnet$coefficients) == 0)) # TN = 438
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.glmnetp$coefficients) == 0)) # FN = 5

# QR####
system.time(L0glm.QR <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                  family = poisson(identity), 
                                  lambda = tune.QR.sc$best.lambda,
                                  nonnegative = TRUE,
                                  relaxed = FALSE,
                                  solver = "eigen",
                                  method = 1L)) 
# 0.83s
cor(as.vector(as.vector(L0glm.QR$coefficients)), as.vector(sim.pois$beta_true)) # 0.9921868
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.QR$coefficients) != 0)) # TP = 45
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.QR$coefficients) != 0)) # FP = 12
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.QR$coefficients) == 0)) # TN = 438
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.osqp$coefficients) == 0)) # FN = 5

# Chol####
system.time(L0glm.chol <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                    family = poisson(identity), 
                                    lambda = tune.chol.sc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "eigen",
                                    method = 7L)) 
# 0.122s
cor(as.vector(as.vector(L0glm.chol$coefficients)), as.vector(sim.pois$beta_true)) # 0.9921868
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.chol$coefficients) != 0)) # TP = 45
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.chol$coefficients) != 0)) # FP = 12
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.chol$coefficients) == 0)) # TN = 438
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.chol$coefficients) == 0)) # FN = 5

# LSCG####
system.time(L0glm.LSCG <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                    family = poisson(identity), 
                                    lambda = tune.LSCG.sc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "eigen",
                                    method = 3L)) 
# 0.536s
cor(as.vector(as.vector(L0glm.LSCG$coefficients)), as.vector(sim.pois$beta_true)) # 0.9308438
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.LSCG$coefficients) != 0)) # TP = 42
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.LSCG$coefficients) != 0)) # FP = 9
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.LSCG$coefficients) == 0)) # TN = 441
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.LSCG$coefficients) == 0)) # FN = 8

# CG####
system.time(L0glm.CG <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                  family = poisson(identity), 
                                  lambda = tune.CG.sc$best.lambda,
                                  nonnegative = TRUE,
                                  relaxed = TRUE,
                                  solver = "eigen",
                                  method = 5L)) 
# 0.384s
cor(as.vector(as.vector(L0glm.CG$coefficients)), as.vector(sim.pois$beta_true)) # 0.938986
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.CG$coefficients) != 0)) # TP = 42
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.CG$coefficients) != 0)) # FP = 9
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.CG$coefficients) == 0)) # TN = 441
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.CG$coefficients) == 0)) # FN = 8

## 1.2. Gaussian####
## Tune lambda relaxed = F 
system.time(tune.osqp.sc <- L0glm.tune(X=sim.gaus$X, y=sim.gaus$y, IC = "aic", plot = FALSE, coefs.true = sim.gaus$beta_true, relaxed = FALSE, family = gaussian(identity), nonnegative = TRUE, warmupiters = 1, solver = "osqp"))
# bestlambda = 3.162277, time = 4.614s
system.time(tune.glmnet.sc <- L0glm.tune(X=sim.gaus$X, y=sim.gaus$y, IC = "aic", plot = FALSE, coefs.true = sim.gaus$beta_true, relaxed = FALSE, family = gaussian(identity), nonnegative = TRUE, warmupiters = 1, solver = "glmnet"))
# bestlambda = 316.2277, time = 5.411s
system.time(tune.chol.sc <- L0glm.tune(X=sim.gaus$X, y=sim.gaus$y, IC = "aic", plot = FALSE, coefs.true = sim.gaus$beta_true, relaxed = FALSE, family = gaussian(identity), nonnegative = TRUE, warmupiters = 1, solver = "eigen", method = 7L))
# bestlambda = 3.162277, time = 3.052s
system.time(tune.QR.sc <- L0glm.tune(X=sim.gaus$X, y=sim.gaus$y, IC = "aic", plot = FALSE, coefs.true = sim.gaus$beta_true, relaxed = TRUE, family = gaussian(identity), nonnegative = TRUE, warmupiters = 1, solver = "eigen", method = 1L))
# bestlambda = 3.162277, time = 15.925s
system.time(tune.CG.sc <- L0glm.tune(X=sim.gaus$X, y=sim.gaus$y, IC = "aic", plot = FALSE, coefs.true = sim.gaus$beta_true, relaxed = FALSE, family = gaussian(identity), nonnegative = TRUE, warmupiters = 1, solver = "eigen", method = 5L))
# bestlambda = 3.162277, time = 7.951s
system.time(tune.LSCG.sc <- L0glm.tune(X=sim.gaus$X, y=sim.gaus$y, IC = "aic", plot = TRUE, coefs.true = sim.gaus$beta_true, relaxed = FALSE, family = gaussian(identity), nonnegative = TRUE, warmupiters = 1, solver = "eigen", method = 3L))
# bestlambda = 3.162277, time = 10.929s

## 1.2.1 Fit relaxed = FALSE####
# OSQP####
system.time(L0glm.osqp <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                    family = gaussian(identity), 
                                    lambda = tune.osqp.sc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "osqp")) 
# 1.085s
cor(as.vector(as.vector(L0glm.osqp$coefficients)), as.vector(sim.gaus$beta_true)) #  0.9961844
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.osqp$coefficients) != 0)) # TP = 47
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.osqp$coefficients) != 0)) # FP = 34
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.osqp$coefficients) == 0)) # TN = 416
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.osqp$coefficients) == 0)) # FN = 3

# GLMNET####
system.time(L0glm.glmnet <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                      family = gaussian(identity), 
                                      lambda = tune.glmnet.sc$best.lambda,
                                      nonnegative = TRUE,
                                      relaxed = FALSE,
                                      solver = "glmnet")) 
# 0.676s
cor(as.vector(as.vector(L0glm.glmnet$coefficients)), as.vector(sim.gaus$beta_true)) # 0.9819911
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.glmnet$coefficients) != 0)) # TP = 45
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.glmnet$coefficients) != 0)) # FP = 13
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.glmnet$coefficients) == 0)) # TN = 437
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.glmnetp$coefficients) == 0)) # FN = 5

# QR####
system.time(L0glm.QR <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                  family = gaussian(identity), 
                                  lambda = tune.QR.sc$best.lambda,
                                  nonnegative = TRUE,
                                  relaxed = FALSE,
                                  solver = "eigen",
                                  method = 1L)) 
# 0.69s
cor(as.vector(as.vector(L0glm.QR$coefficients)), as.vector(sim.gaus$beta_true)) # 0.9975929
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.QR$coefficients) != 0)) # TP = 47
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.QR$coefficients) != 0)) # FP = 26
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.QR$coefficients) == 0)) # TN = 424
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.osqp$coefficients) == 0)) # FN = 3

# Chol####
system.time(L0glm.chol <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                    family = gaussian(identity), 
                                    lambda = tune.chol.sc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "eigen",
                                    method = 7L)) 
# 0.107s
cor(as.vector(as.vector(L0glm.chol$coefficients)), as.vector(sim.gaus$beta_true)) # 0.9975929
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.chol$coefficients) != 0)) # TP = 47
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.chol$coefficients) != 0)) # FP = 26
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.chol$coefficients) == 0)) # TN = 424
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.chol$coefficients) == 0)) # FN = 3

# LSCG####
system.time(L0glm.LSCG <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                    family = gaussian(identity), 
                                    lambda = tune.CG.sc$best.lambda,
                                    nonnegative = TRUE,
                                    start = rep(1, 500),
                                    relaxed = FALSE,
                                    solver = "eigen",
                                    method = 3L)) 
# 0.536s
cor(as.vector(as.vector(L0glm.LSCG$coefficients)), as.vector(sim.gaus$beta_true)) # 0.4905367
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.LSCG$coefficients) != 0)) # TP = 50
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.LSCG$coefficients) != 0)) # FP = 301
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.LSCG$coefficients) == 0)) # TN = 149
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.LSCG$coefficients) == 0)) # FN = 0

# CG####
system.time(L0glm.CG <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                  family = gaussian(identity), 
                                  lambda = tune.CG.sc$best.lambda,
                                  nonnegative = TRUE,
                                  relaxed = TRUE,
                                  solver = "eigen",
                                  method = 5L)) 
# 0.384s
cor(as.vector(as.vector(L0glm.CG$coefficients)), as.vector(sim.gaus$beta_true)) # 0.5015072
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.CG$coefficients) != 0)) # TP = 50
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.CG$coefficients) != 0)) # FP = 288
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.CG$coefficients) == 0)) # TN = 162
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.CG$coefficients) == 0)) # FN = 0


# 2. MULTICHANNEL####
# Simulate MC n=p=500
n <- p <- 500
sim.gaus.mc <- simulate(method = 3, n = n, p = p, seed = 1)
sim.pois.mc <- simulate(method = 4, n = n, p = p, seed = 2)


## 2.1. Poisson####
## 2.1.1 Relaxed = FALSE, grouped = FALSE####
## Tune lambda relaxed = FALSE, grouped = FAKSE
system.time(tune.osqp.mc <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, IC = "bic", 
                                       plot = FALSE, coefs.true = sim.pois.mc$beta_true, 
                                       relaxed = FALSE, family = poisson(identity), 
                                       nonnegative = TRUE, warmupiters = 1, 
                                       parallellize = TRUE, grouped = FALSE, solver = "osqp"))

# bestlambda = 0.3162, time = 175.216s
system.time(tune.chol.mc <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, IC = "aic", 
                                       plot = FALSE, coefs.true = sim.pois.mc$beta_true, 
                                       relaxed = FALSE, family = poisson(identity), 
                                       nonnegative = TRUE, warmupiters = 1, 
                                       parallellize = TRUE, grouped = FALSE,
                                       solver = "eigen", method = 7L))
# bestlambda = 1, time = 135.984 s

## Fit relaxed = FALSE, grouped = FALSE
# OSQP####
system.time(L0glm.osqp <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                    family = poisson(identity), 
                                    lambda = tune.osqp.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallellize = TRUE,
                                    grouped = FALSE,
                                    solver = "osqp")) 
# 0.125s
abline(v=which(rowMin(L0glm.osqp$coefficients)!=0), col = 'cyan')

cor(as.vector(as.vector(L0glm.osqp$coefficients)), as.vector(sim.pois.mc$beta_true)) # 0.9909261
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)!=0) # TP = 50
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)!=0) # FP = 241
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)==0) # TN = 209
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)==0) # FN = 0

# Chol####
system.time(L0glm.chol <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                    family = poisson(identity), 
                                    lambda = tune.chol.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallellize = TRUE,
                                    grouped = FALSE,
                                    solver = "eigen",
                                    method = 7L)) 
# 9.932s
abline(v=which(rowMin(L0glm.chol$coefficients)!=0), col = 'cyan')

R2 = round(cor(as.vector(rowMeans(sim.pois.mc$beta_true)), as.vector(rowMeans(L0glm.chol$coefficients))), digits = 4)^2 # Calculate R2
title = paste("True vs estimated coefficients\nnonnegative L0 independent GLM\nR2=", R2)
plot(as.vector(rowMeans(L0glm.chol$coefficients)), as.vector(rowMeans(sim.pois.mc$beta_true)), pch=16, 
     xlab="estimated coefficients nonnegative L0 group GLM",
     ylab="true coefficients", col="steelblue",
     main=title)

cor(as.vector(as.vector(L0glm.chol$coefficients)), as.vector(sim.pois.mc$beta_true)) # 0.9856273
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.chol$coefficients)!=0) # TP = 50
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.chol$coefficients)!=0) # FP = 241
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.chol$coefficients)==0) # TN = 209
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.chol$coefficients)==0) # FN = 0


## 2.1.2 Relaxed = FALSE, grouped = TRUE####
## Tune lambda relaxed = FALSE, grouped = TRUE
system.time(tune.osqp.mc <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, IC = "bic", 
                                       plot = FALSE, coefs.true = sim.pois.mc$beta_true, 
                                       relaxed = FALSE, family = poisson(identity), 
                                       nonnegative = TRUE, warmupiters = 1, 
                                       parallellize = TRUE, grouped = TRUE, solver = "osqp"))

# bestlambda = 1, time = 477.367s
system.time(tune.chol.mc <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, IC = "bic", 
                                       plot = FALSE, coefs.true = sim.pois.mc$beta_true, 
                                       relaxed = FALSE, family = poisson(identity), 
                                       nonnegative = TRUE, warmupiters = 1, 
                                       parallellize = TRUE, grouped = TRUE,
                                       solver = "eigen", method = 7L))
# bestlambda = 1, time = 385.675s

## Fit relaxed = FALSE, grouped = TRUE
# OSQP####
system.time(L0glm.osqp <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                    family = poisson(identity), 
                                    lambda = tune.osqp.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallellize = TRUE,
                                    grouped = TRUE,
                                    solver = "osqp"))
# 30.511s
abline(v=which(rowMin(L0glm.osqp$coefficients)!=0), col = 'cyan')

cor(as.vector(as.vector(L0glm.osqp$coefficients)), as.vector(sim.pois.mc$beta_true)) # 0.9972466
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)!=0) # TP = 49
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)!=0) # FP = 1
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)==0) # TN = 449
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)==0) # FN = 1

# Chol####
system.time(L0glm.chol <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                    family = poisson(identity), 
                                    lambda = tune.chol.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallellize = TRUE,
                                    grouped = TRUE,
                                    solver = "eigen",
                                    method = 7L)) 
# 19.245s
abline(v=which(rowMin(L0glm.chol$coefficients)!=0), col = 'cyan')

cor(as.vector(as.vector(L0glm.chol$coefficients)), as.vector(sim.pois.mc$beta_true)) # 0.9994343
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.chol$coefficients)!=0) # TP = 49
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.chol$coefficients)!=0) # FP = 1
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.chol$coefficients)==0) # TN = 449
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.chol$coefficients)==0) # FN = 1

## 2.1.3 Relaxed = TRUE, grouped = TRUE####
## Tune lambda relaxed = TRUE, grouped = TRUE
system.time(tune.osqp.mc <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, IC = "bic", 
                                       plot = FALSE, coefs.true = sim.pois.mc$beta_true, 
                                       relaxed = TRUE, family = poisson(identity), 
                                       nonnegative = TRUE, warmupiters = 1, 
                                       parallellize = TRUE, grouped = TRUE, solver = "osqp"))

# bestlambda = 1, time = 1166.551s
system.time(tune.chol.mc <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, IC = "bic", 
                                       plot = FALSE, coefs.true = sim.pois.mc$beta_true, 
                                       relaxed = TRUE, family = poisson(identity), 
                                       nonnegative = TRUE, warmupiters = 1, 
                                       parallellize = TRUE, grouped = TRUE,
                                       solver = "eigen", method = 7L))
# bestlambda = 1, time = 570.158s

## Fit relaxed = TRUE, grouped = TRUE
# OSQP####
system.time(L0glm.osqp <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                    family = poisson(identity), 
                                    lambda = tune.osqp.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = TRUE,
                                    parallellize = TRUE,
                                    grouped = TRUE,
                                    solver = "osqp")) 
# 36.663s
abline(v=which(rowMin(L0glm.osqp$coefficients)!=0), col = 'cyan')

cor(as.vector(as.vector(L0glm.osqp$coefficients)), as.vector(sim.pois.mc$beta_true)) # 0.999498
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)!=0) # TP = 49
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)!=0) # FP = 449
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)==0) # TN = 1
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)==0) # FN = 1

# Chol####
system.time(L0glm.chol <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                    family = poisson(identity), 
                                    lambda = tune.chol.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = TRUE,
                                    parallellize = TRUE,
                                    grouped = TRUE,
                                    solver = "eigen",
                                    method = 7L)) 
# 30.235s
abline(v=which(rowMin(L0glm.chol$coefficients)!=0), col = 'cyan')

R2 = round(cor(as.vector(rowMeans(sim.pois.mc$beta_true)), as.vector(rowMeans(L0glm.chol$coefficients))), digits = 4)^2 # Calculate R2
title = paste("True vs estimated coefficients\nnonnegative L0 grouped relaxed GLM\nR2=", R2)
plot(as.vector(rowMeans(L0glm.chol$coefficients)), as.vector(rowMeans(sim.pois.mc$beta_true)), pch=16, 
     xlab="estimated coefficients nonnegative L0 group GLM",
     ylab="true coefficients", col="steelblue",
     main=title)

cor(as.vector(as.vector(L0glm.chol$coefficients)), as.vector(sim.pois.mc$beta_true)) # 0.9994725
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.chol$coefficients)!=0) # TP = 49
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.chol$coefficients)!=0) # FP = 1
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.chol$coefficients)==0) # TN = 449
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.chol$coefficients)==0) # FN = 1


## 2.2. Gaussian####
## 2.2.1 Relaxed = FALSE, grouped = FALSE####
## Tune lambda relaxed = FALSE, grouped = FALSE
system.time(tune.osqp.mc <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, IC = "bic", 
                                       plot = FALSE, coefs.true = sim.gaus.mc$beta_true, 
                                       relaxed = FALSE, family = gaussian(identity), 
                                       nonnegative = TRUE, warmupiters = 1, 
                                       parallellize = TRUE, grouped = FALSE, solver = "osqp"))

# bestlambda = 100, time = 201.820s
system.time(tune.chol.mc <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, IC = "bic", 
                                       plot = FALSE, coefs.true = sim.gaus.mc$beta_true, 
                                       relaxed = FALSE, family = gaussian(identity), 
                                       nonnegative = TRUE, warmupiters = 1, 
                                       parallellize = TRUE, grouped = FALSE,
                                       solver = "eigen", method = 7L))
# bestlambda = 100, time = 172.151s

## Fit relaxed = FALSE, grouped = FALSE
# OSQP####
system.time(L0glm.osqp <- L0glm.fit(sim.gaus.mc$X, sim.gaus.mc$y, 
                                    family = gaussian(identity), 
                                    lambda = tune.osqp.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallellize = TRUE,
                                    grouped = FALSE,
                                    solver = "osqp")) 
# 17.887s
abline(v=which(rowMin(L0glm.osqp$coefficients)!=0), col = 'cyan')

cor(as.vector(as.vector(L0glm.osqp$coefficients)), as.vector(sim.gaus.mc$beta_true)) # 0.9887807
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)!=0) # TP = 50
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)!=0) # FP = 411
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)==0) # TN = 39
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)==0) # FN = 0

# Chol####
system.time(L0glm.chol <- L0glm.fit(sim.gaus.mc$X, sim.gaus.mc$y, 
                                    family = gaussian(identity), 
                                    lambda = tune.chol.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallellize = TRUE,
                                    grouped = FALSE,
                                    solver = "eigen",
                                    method = 7L)) 
# 10.584s
abline(v=which(rowMin(L0glm.chol$coefficients)!=0), col = 'cyan')

R2 = round(cor(as.vector(rowMeans(sim.gaus.mc$beta_true)), as.vector(rowMeans(L0glm.chol$coefficients))), digits = 4)^2 # Calculate R2
title = paste("True vs estimated coefficients\nnonnegative L0 independent GLM\nR2=", R2)
plot(as.vector(rowMeans(L0glm.chol$coefficients)), as.vector(rowMeans(sim.gaus.mc$beta_true)), pch=16, 
     xlab="estimated coefficients nonnegative L0 group GLM",
     ylab="true coefficients", col="steelblue",
     main=title)

cor(as.vector(as.vector(L0glm.chol$coefficients)), as.vector(sim.gaus.mc$beta_true)) # 0.885
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.chol$coefficients)!=0) # TP = 50
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.chol$coefficients)!=0) # FP = 450
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.chol$coefficients)==0) # TN = 0
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.chol$coefficients)==0) # FN = 0


## 2.2.2 Relaxed = FALSE, grouped = TRUE####
## Tune lambda relaxed = FALSE, grouped = TRUE
system.time(tune.osqp.mc <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, IC = "bic", 
                                       plot = TRUE, coefs.true = sim.gaus.mc$beta_true, 
                                       relaxed = FALSE, family = gaussian(identity), 
                                       nonnegative = TRUE, warmupiters = 1, 
                                       parallellize = TRUE, grouped = TRUE, solver = "osqp"))

# bestlambda = 10, time = 756.785s
system.time(tune.chol.mc <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, IC = "bic", 
                                       plot = FALSE, coefs.true = sim.gaus.mc$beta_true, 
                                       relaxed = FALSE, family = gaussian(identity), 
                                       nonnegative = TRUE, warmupiters = 1, 
                                       parallellize = TRUE, grouped = TRUE,
                                       solver = "eigen", method = 7L))
# bestlambda = 10, time = 31748.790s TODO make sure this is correct?? 

## Fit relaxed = FALSE, grouped = TRUE
# OSQP####
system.time(L0glm.osqp <- L0glm.fit(sim.gaus.mc$X, sim.gaus.mc$y, 
                                    family = gaussian(identity), 
                                    lambda = tune.osqp.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallellize = TRUE,
                                    grouped = TRUE,
                                    solver = "osqp")) 
# 79.685 s
abline(v=which(rowMin(L0glm.osqp$coefficients)!=0), col = 'cyan')

cor(as.vector(as.vector(L0glm.osqp$coefficients)), as.vector(sim.gaus.mc$beta_true)) # 0.9999549
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)!=0) # TP = 50
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)!=0) # FP = 0
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)==0) # TN = 450
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)==0) # FN = 0

# Chol####
system.time(L0glm.chol <- L0glm.fit(sim.gaus.mc$X, sim.gaus.mc$y, 
                                    family = gaussian(identity), 
                                    lambda = tune.chol.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallellize = TRUE,
                                    grouped = TRUE,
                                    solver = "eigen",
                                    method = 7L)) 
# 56.365s
abline(v=which(rowMin(L0glm.chol$coefficients)!=0), col = 'cyan')

R2 = round(cor(as.vector(rowMeans(sim.gaus.mc$beta_true)), as.vector(rowMeans(L0glm.chol$coefficients))), digits = 4)^2 # Calculate R2
title = paste("True vs estimated coefficients\nnonnegative L0 grouped GLM\nR2=", R2)
plot(as.vector(rowMeans(L0glm.chol$coefficients)), as.vector(rowMeans(sim.gaus.mc$beta_true)), pch=16, 
     xlab="estimated coefficients nonnegative L0 group GLM",
     ylab="true coefficients", col="steelblue",
     main=title)

cor(as.vector(as.vector(L0glm.chol$coefficients)), as.vector(sim.gaus.mc$beta_true)) # 0.999955
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.chol$coefficients)!=0) # TP = 50
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.chol$coefficients)!=0) # FP = 0
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.chol$coefficients)==0) # TN = 450
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.chol$coefficients)==0) # FN = 0

## 2.2.3 Relaxed = TRUE, grouped = TRUE####
## Tune lambda relaxed = TRUE, grouped = TRUE
system.time(tune.osqp.mc <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, IC = "bic", 
                                       plot = TRUE, coefs.true = sim.gaus.mc$beta_true, 
                                       relaxed = TRUE, family = gaussian(identity), 
                                       nonnegative = TRUE, warmupiters = 1, 
                                       parallellize = TRUE, grouped = TRUE, solver = "osqp"))

# bestlambda = , time = 1038.3s
system.time(tune.chol.mc <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, IC = "bic", 
                                       plot = FALSE, coefs.true = sim.gaus.mc$beta_true, 
                                       relaxed = TRUE, family = gaussian(identity), 
                                       nonnegative = TRUE, warmupiters = 1, 
                                       parallellize = TRUE, grouped = TRUE,
                                       solver = "eigen", method = 7L))
# bestlambda = , time = 689.32s

## Fit relaxed = TRUE, grouped = TRUE
# OSQP####
system.time(L0glm.osqp <- L0glm.fit(sim.gaus.mc$X, sim.gaus.mc$y, 
                                    family = gaussian(identity), 
                                    lambda = tune.osqp.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = TRUE,
                                    parallellize = TRUE,
                                    grouped = TRUE,
                                    solver = "osqp")) 
# 36.663s
abline(v=which(rowMin(L0glm.osqp$coefficients)!=0), col = 'cyan')

cor(as.vector(as.vector(L0glm.osqp$coefficients)), as.vector(sim.gaus.mc$beta_true)) # 
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)!=0) # TP = 
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)!=0) # FP = 
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)==0) # TN = 
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)==0) # FN = 

# Chol####
system.time(L0glm.chol <- L0glm.fit(sim.gaus.mc$X, sim.gaus.mc$y, 
                                    family = gaussian(identity), 
                                    lambda = tune.chol.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = TRUE,
                                    parallellize = TRUE,
                                    grouped = TRUE,
                                    solver = "eigen",
                                    method = 7L)) 
# 30.235s
abline(v=which(rowMin(L0glm.chol$coefficients)!=0), col = 'cyan')

R2 = round(cor(as.vector(rowMeans(sim.gaus.mc$beta_true)), as.vector(rowMeans(L0glm.chol$coefficients))), digits = 4)^2 # Calculate R2
title = paste("True vs estimated coefficients\nnonnegative L0 grouped relaxed GLM\nR2=", R2)
plot(as.vector(rowMeans(L0glm.chol$coefficients)), as.vector(rowMeans(sim.gaus.mc$beta_true)), pch=16, 
     xlab="estimated coefficients nonnegative L0 group GLM",
     ylab="true coefficients", col="steelblue",
     main=title)

cor(as.vector(as.vector(L0glm.chol$coefficients)), as.vector(sim.gaus.mc$beta_true)) # 
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.chol$coefficients)!=0) # TP = 
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.chol$coefficients)!=0) # FP = 
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.chol$coefficients)==0) # TN = 
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.chol$coefficients)==0) # FN = 


# 3. COMPARISON WITH OTHER METHODS####
# Take our best working method and compare it against others: Glmnet, ordinis, L0learn, abess.
# 3.1 SINGLE CHANNEL POISSON####
n <- p <- 500
sim.pois <- simulate(method = 2, n = n, p = p, seed = 2)

# Tune and fit with Cholesky decomposition ####
system.time(tune.chol.sc <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "aic", 
                                       plot = FALSE, coefs.true = sim.pois$beta_true, 
                                       relaxed = FALSE, family = poisson(identity), 
                                       nonnegative = TRUE, warmupiters = 1, 
                                       solver = "eigen", method = 7L))

system.time(L0glm.chol <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                    family = poisson(identity), 
                                    lambda = tune.chol.sc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "eigen",
                                    method = 7L)) 
# 0.122s
cor1 <- cor(as.vector(as.vector(L0glm.chol$coefficients)), as.vector(sim.pois$beta_true)) # 0.9675023
TP1 <- sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.chol$coefficients) != 0)) # TP = 40
FP1 <- sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.chol$coefficients) != 0)) # FP = 16
TN1 <- sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.chol$coefficients) == 0)) # TN = 434
FN1 <- sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.chol$coefficients) == 0)) # FN = 10

plot(x=L0glm.chol$coefficients, y=sim.pois$beta_true, pch=16, col='blue',
     main="L0 glm identity link Poisson with lambda tuned to optimise AIC",
     xlab="Estimated coefficients", ylab="True coefficients")

# Nonnegative weighted lasso regression in glmnet ####
library(glmnet)
system.time(glmnet.sc.pois <- fit.GLM(sim = sim.pois, method = 2, family = poisson(identity), plot = TRUE)) # 0.158
cor2 <- cor(as.vector(glmnet.sc.pois[-1]), as.vector(sim.pois$beta_true)) # 0.8545053
TP2 <- sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(glmnet.sc.pois[-1]) != 0)) # TP = 37
FP2 <- sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(glmnet.sc.pois[-1]) != 0)) # FP = 59
TN2 <- sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(glmnet.sc.pois[-1]) == 0)) # TN = 391
FN2 <- sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(glmnet.sc.pois[-1]) == 0)) # FN = 13


# Nonnegative adaptive lasso regression in glmnet ####
system.time(glmnet.al.sc.pois <- fit.GLM(sim = sim.pois, method = 3, family = poisson(identity), plot = T)) # 0.249 
cor3 <- cor(as.vector(glmnet.al.sc.pois[-1]), as.vector(sim.pois$beta_true)) # 0.8913109
TP3 <- sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(glmnet.al.sc.pois[-1]) != 0)) # TP = 28
FP3 <- sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(glmnet.al.sc.pois[-1]) != 0)) # FP = 30
TN3 <- sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(glmnet.al.sc.pois[-1]) == 0)) # TN = 420
FN3 <- sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(glmnet.al.sc.pois[-1]) == 0)) # FN = 22


# nonnegative L0Learn best subset ####
library(L0Learn)
system.time(L0learn.sc.pois <- fit.GLM(sim = sim.pois, method = 4, family = poisson(identity), plot = T))  # 1.120 
cor(as.vector(L0learn.sc.pois), as.vector(sim.pois$beta_true)) # 0.8519745
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0learn.sc.pois) != 0)) # TP = 22
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0learn.sc.pois) != 0)) # FP = 11
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0learn.sc.pois) == 0)) # TN = 439
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0learn.sc.pois) == 0)) # FN = 28

# nonnegative L0Learn best subset fit after glmnet nonnegative LASSO ####
# Same results as regular L0learn but faster
system.time(L0learn.after.lasso.sc.pois <- fit.GLM(sim = sim.pois, method = 5, family = poisson(identity), plot = T))  # 0.242 
cor(as.vector(L0learn.after.lasso.sc.pois), as.vector(sim.pois$beta_true)) # 0.8688236
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0learn.after.lasso.sc.pois) != 0)) # TP = 26
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0learn.after.lasso.sc.pois) != 0)) # FP = 17
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0learn.after.lasso.sc.pois) == 0)) # TN = 433
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0learn.after.lasso.sc.pois) == 0)) # FN = 24

# abess best subset fit ####
# Super bad solution
library(abess)
system.time(abess.sc.pois <- fit.GLM(sim = sim.pois, method = 6, family = poisson(identity), plot = T, IC = "ebic"))  # 52.256 
cor(as.vector(abess.sc.pois), as.vector(sim.pois$beta_true)) # 0.01254078
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(abess.sc.pois) != 0)) # TP = 27
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(abess.sc.pois) != 0)) # FP = 112
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(abess.sc.pois) == 0)) # TN = 338
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(abess.sc.pois) == 0)) # FN = 23


# abess best subset fit after nonnegative LASSO ####
system.time(abess.after.lasso.sc.pois <- fit.GLM(sim = sim.pois, method = 7, family = poisson(identity), plot = T, IC = "ebic"))  # 0.272 s
cor(as.vector(abess.after.lasso.sc.pois), as.vector(sim.pois$beta_true))  # 0.8246482
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(abess.after.lasso.sc.pois) != 0)) # TP = 31
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(abess.after.lasso.sc.pois) != 0)) # FP = 24
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(abess.after.lasso.sc.pois) == 0)) # TN = 426
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(abess.after.lasso.sc.pois) == 0)) # FN = 19

# nonnegative weighted MCP penalized regression in ordinis ####
library(ordinis)
system.time(ordinis.mcp.sc.pois <- fit.GLM(sim = sim.pois, method = 8, family = poisson(identity), plot = T, IC = "ebic"))  # 33.463 s
cor(as.vector(ordinis.mcp.sc.pois[-1]), as.vector(sim.pois$beta_true)) # 0.6024787
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(ordinis.mcp.sc.pois[-1]) != 0)) # TP = 32
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(ordinis.mcp.sc.pois[-1]) != 0)) # FP = 75
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(ordinis.mcp.sc.pois[-1]) == 0)) # TN = 375
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(ordinis.mcp.sc.pois[-1]) == 0)) # FN = 18

# nonnegative weighted SCAD penalized regression in ordinis ####
system.time(ordinis.scad.sc.pois <- fit.GLM(sim = sim.pois, method = 9, family = poisson(identity), plot = T, IC = "ebic"))  # 34.296 s
cor(as.vector(ordinis.scad.sc.pois[-1]), as.vector(sim.pois$beta_true)) # 0.6079265
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(ordinis.scad.sc.pois[-1]) != 0)) # TP = 34
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(ordinis.scad.sc.pois[-1]) != 0)) # FP = 85
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(ordinis.scad.sc.pois[-1]) == 0)) # TN = 365
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(ordinis.scad.sc.pois[-1]) == 0)) # FN = 16


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
sim.gaus <- simulate(method = 1, n = n, p = p, seed = 2)

# Tune and fit with Cholesky decomposition ####
system.time(tune.chol.sc <- L0glm.tune(X=sim.gaus$X, y=sim.gaus$y, IC = "bic", 
                                       plot = FALSE, coefs.true = sim.gaus$beta_true, 
                                       relaxed = FALSE, family = gaussian(identity), 
                                       nonnegative = TRUE, warmupiters = 1, 
                                       solver = "eigen", method = 7L))

system.time(L0glm.chol <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                    family = gaussian(identity), 
                                    lambda = tune.chol.sc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "eigen",
                                    method = 7L)) 
# 0.122s
cor(as.vector(as.vector(L0glm.chol$coefficients)), as.vector(sim.gaus$beta_true)) # 0.9915894
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.chol$coefficients) != 0)) # TP = 48
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.chol$coefficients) != 0)) # FP = 3
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.chol$coefficients) == 0)) # TN = 447
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.chol$coefficients) == 0)) # FN = 2

plot(x=L0glm.chol$coefficients, y=sim.gaus$beta_true, pch=16, col='blue',
     main="L0 glm identity link gaussian with lambda tuned to optimise AIC",
     xlab="Estimated coefficients", ylab="True coefficients")

# Nonnegative weighted lasso regression in glmnet ####
library(glmnet)
system.time(glmnet.sc.pois <- fit.GLM(sim = sim.gaus, method = 2, family = gaussian(identity), plot = TRUE)) # 0.132 
cor(as.vector(glmnet.sc.pois[-1]), as.vector(sim.gaus$beta_true)) # 0.9613213
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(glmnet.sc.pois[-1]) != 0)) # TP = 49
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(glmnet.sc.pois[-1]) != 0)) # FP = 66
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(glmnet.sc.pois[-1]) == 0)) # TN = 384
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(glmnet.sc.pois[-1]) == 0)) # FN = 1



system.time(glmnet.al.sc.pois <- fit.GLM(sim = sim.gaus, method = 3, family = gaussian(identity), plot = T)) # 0.172 
cor(as.vector(glmnet.al.sc.pois[-1]), as.vector(sim.gaus$beta_true)) # 0.9557091
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(glmnet.al.sc.pois[-1]) != 0)) # TP = 37
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(glmnet.al.sc.pois[-1]) != 0)) # FP = 12
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(glmnet.al.sc.pois[-1]) == 0)) # TN = 438
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(glmnet.al.sc.pois[-1]) == 0)) # FN = 13


# nonnegative L0Learn best subset ####
library(L0Learn)
system.time(L0learn.sc.pois <- fit.GLM(sim = sim.gaus, method = 4, family = gaussian(identity), plot = T))  #  1.151 s
cor(as.vector(L0learn.sc.pois), as.vector(sim.gaus$beta_true)) # 0.887
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0learn.sc.pois) != 0)) # TP = 31
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0learn.sc.pois) != 0)) # FP = 55
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0learn.sc.pois) == 0)) # TN = 395
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0learn.sc.pois) == 0)) # FN = 19

# nonnegative L0Learn best subset fit after glmnet nonnegative LASSO ####
# Same results as regular L0learn but faster
system.time(L0learn.after.lasso.sc.pois <- fit.GLM(sim = sim.gaus, method = 5, family = gaussian(identity), plot = T))  # 1.000 
cor(as.vector(L0learn.after.lasso.sc.pois), as.vector(sim.gaus$beta_true)) # 0.9475639
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0learn.after.lasso.sc.pois) != 0)) # TP = 45
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0learn.after.lasso.sc.pois) != 0)) # FP = 44
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0learn.after.lasso.sc.pois) == 0)) # TN = 406
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0learn.after.lasso.sc.pois) == 0)) # FN = 5

# abess best subset fit ####
# Super bad solution
library(abess)
system.time(abess.sc.pois <- fit.GLM(sim = sim.gaus, method = 6, family = gaussian(identity), plot = T, IC = "ebic"))  # 1.139 
cor(as.vector(abess.sc.pois), as.vector(sim.gaus$beta_true)) # 0.584285
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(abess.sc.pois) != 0)) # TP = 32
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(abess.sc.pois) != 0)) # FP = 67
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(abess.sc.pois) == 0)) # TN = 383
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(abess.sc.pois) == 0)) # FN = 18


# abess best subset fit after nonnegative LASSO ####
system.time(abess.after.lasso.sc.pois <- fit.GLM(sim = sim.gaus, method = 7, family = gaussian(identity), plot = T, IC = "ebic"))  # 0.269 s
cor(as.vector(abess.after.lasso.sc.pois), as.vector(sim.gaus$beta_true))  # 0.9365722
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(abess.after.lasso.sc.pois) != 0)) # TP = 42
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(abess.after.lasso.sc.pois) != 0)) # FP = 34
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(abess.after.lasso.sc.pois) == 0)) # TN = 416
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(abess.after.lasso.sc.pois) == 0)) # FN = 8

# nonnegative weighted MCP penalized regression in ordinis ####
library(ordinis)
system.time(ordinis.mcp.sc.pois <- fit.GLM(sim = sim.gaus, method = 8, family = gaussian(identity), plot = T, IC = "ebic"))  # 0.529s
cor(as.vector(ordinis.mcp.sc.pois[-1]), as.vector(sim.gaus$beta_true)) # 0.7706272
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(ordinis.mcp.sc.pois[-1]) != 0)) # TP = 29
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(ordinis.mcp.sc.pois[-1]) != 0)) # FP = 34
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(ordinis.mcp.sc.pois[-1]) == 0)) # TN = 416
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(ordinis.mcp.sc.pois[-1]) == 0)) # FN = 21

# nonnegative weighted SCAD penalized regression in ordinis ####
system.time(ordinis.scad.sc.pois <- fit.GLM(sim = sim.gaus, method = 9, family = gaussian(identity), plot = T, IC = "ebic"))  # 0.889 
cor(as.vector(ordinis.scad.sc.pois[-1]), as.vector(sim.gaus$beta_true)) # 0.9208366
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
sim.pois.mc <- simulate(method = 4, n = n, p = p, seed = 2)

# Tune and fit with Cholesky decomposition ####
system.time(tune.chol.mc <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, IC = "aic", 
                                       plot = FALSE, coefs.true = sim.pois.mc$beta_true, 
                                       relaxed = FALSE, family = poisson(identity), 
                                       nonnegative = TRUE, warmupiters = 1, 
                                       parallellize = TRUE, grouped = TRUE,
                                       solver = "eigen", method = 7L)) # lambda = 1

system.time(L0glm.chol.pois <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                         family = poisson(identity), 
                                         lambda = 1,
                                         nonnegative = TRUE,
                                         relaxed = FALSE,
                                         parallellize = TRUE,
                                         grouped = TRUE,
                                         solver = "eigen",
                                         method = 7L)) 
# 34.459 s
abline(v=which(rowMin(L0glm.chol.pois$coefficients)!=0), col = 'cyan')

R2 = round(cor(as.vector(rowMeans(sim.pois.mc$beta_true)), as.vector(rowMeans(L0glm.chol.pois$coefficients))), digits = 4)^2 # Calculate R2
title = paste("True vs estimated coefficients\nnonnegative L0 grouped GLM\nR2=", R2)
plot(as.vector(rowMeans(L0glm.chol.pois$coefficients)), as.vector(rowMeans(sim.pois.mc$beta_true)), pch=16, 
     xlab="estimated coefficients nonnegative L0 group GLM",
     ylab="true coefficients", col="steelblue",
     main=title)


image(x=1:nrow(sim$y), y=1:ncol(sim.pois.mc$y), z=as.matrix(L0glm.chol.pois$coefficients)^0.1, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="L0glm\n(red = true peaks, cyan = estimated peaks)")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim.pois.mc$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim$X))[as.vector(rowMin(L0glm.chol.pois$coefficients)!=0)], col="cyan")

cor(as.vector(as.vector(L0glm.chol.pois$coefficients)), as.vector(sim.pois.mc$beta_true)) # 0.9994343
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.chol.pois$coefficients)!=0) # TP = 49
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.chol.pois$coefficients)!=0) # FP = 1
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.chol.pois$coefficients)==0) # TN = 449
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.chol.pois$coefficients)==0) # FN = 1


# best subset multitask learning using family="mgaussian" using abess ####
system.time(bsgs.abess.poisson <- fit.GLM.mc(sim = sim.pois.mc, method = 1, family = poisson(identity), plot = T, IC = "gic")) 
# 33.131 s


cor(as.vector(as.vector(sim.pois.mc$beta_true)), as.vector(bsgs.abess.poisson)) # 0.1213357
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(bsgs.abess.poisson)!=0) # TP = 34
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(bsgs.abess.poisson)!=0) # FP = 47
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(bsgs.abess.poisson)==0) # TN = 403
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(bsgs.abess.poisson)==0) # FN = 16

# 3.4 MULTICHANNEL GAUSSIAN ####
n <- p <- 500
sim.gaus.mc <- simulate(method = 3, n = n, p = p, seed = 1)

# Tune and fit with Cholesky decomposition ####
system.time(tune.chol.mc <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, IC = "aic", 
                                       plot = FALSE, coefs.true = sim.gaus.mc$beta_true, 
                                       relaxed = FALSE, family = gaussian(identity), 
                                       nonnegative = TRUE, warmupiters = 1, 
                                       parallellize = TRUE, grouped = TRUE,
                                       solver = "eigen", method = 7L)) # lambda = 10

system.time(L0glm.chol.gaus <- L0glm.fit(sim.gaus.mc$X, sim.gaus.mc$y, 
                                    family = gaussian(identity), 
                                    lambda = 10,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallellize = TRUE,
                                    grouped = TRUE,
                                    solver = "eigen",
                                    method = 7L)) 
# 35.681 s
abline(v=which(rowMin(L0glm.chol.gaus$coefficients)!=0), col = 'cyan')

R2 = round(cor(as.vector(rowMeans(sim.gaus.mc$beta_true)), as.vector(rowMeans(L0glm.chol.gaus$coefficients))), digits = 4)^2 # Calculate R2
title = paste("True vs estimated coefficients\nnonnegative L0 grouped GLM\nR2=", R2)
plot(as.vector(rowMeans(L0glm.chol.gaus$coefficients)), as.vector(rowMeans(sim.gaus.mc$beta_true)), pch=16, 
     xlab="estimated coefficients nonnegative L0 group GLM",
     ylab="true coefficients", col="steelblue",
     main=title)


image(x=1:nrow(sim$y), y=1:ncol(sim.gaus.mc$y), z=as.matrix(L0glm.chol.gaus$coefficients)^0.1, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="L0glm\n(red = true peaks, cyan = estimated peaks)")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim.gaus.mc$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim$X))[as.vector(rowMin(L0glm.chol.gaus$coefficients)!=0)], col="cyan")

cor(as.vector(as.vector(L0glm.chol.gaus$coefficients)), as.vector(sim.gaus.mc$beta_true)) # 0.9963284
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.chol.gaus$coefficients)!=0) # TP = 
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.chol.gaus$coefficients)!=0) # FP = 
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.chol.gaus$coefficients)==0) # TN = 
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.chol.gaus$coefficients)==0) # FN = 


# best subset multitask learning using family="mgaussian" using abess ####
system.time(bsgs.abess.gaus <- fit.GLM.mc(sim = sim.gaus.mc, method = 1, family = gaussian(identity), plot = T, IC = "gic")) 
# 38.391 s


cor(as.vector(as.vector(sim.gaus.mc$beta_true)), as.vector(bsgs.abess.gaus)) # 0.01352945
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(bsgs.abess.gaus)!=0) # TP = 31
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(bsgs.abess.gaus)!=0) # FP = 79
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(bsgs.abess.gaus)==0) # TN = 371
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(bsgs.abess.gaus)==0) # FN = 19

# # glmnet nonnegative mgaussian multivariate group LASSO fit
# system.time(glmnet.mc <- fit.GLM.mc(sim = sim.gaus.mc, method = 2, family = gaussian(identity), plot = T, IC = "gic")) 
# # TODO fix, outputs wrong dimensions
# # 54.257s
# 
# 
# image(x=1:nrow(sim$y), y=1:ncol(sim.gaus.mc$y), z=as.matrix(glmnet.mc)^0.1, col = topo.colors(255),
#       useRaster=TRUE,
#       xlab="Time", ylab="Channel", main="L0glm\n(red = true peaks, cyan = estimated peaks)")
# abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim.gaus.mc$beta_true)!=0)], col="red")
# abline(v=(1:nrow(sim$X))[as.vector(rowMin(glmnet.mc)!=0)], col="cyan")
# 
# 
# cor(as.vector(as.vector(sim.gaus.mc$beta_true)), as.vector(glmnet.mc)) # 
# sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(glmnet.mc)!=0) # 
# sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(glmnet.mc)!=0) # 
# sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(glmnet.mc)==0) # 
# sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(glmnet.mc)==0) # 
# 

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
