# BENCHMARKS L0GLM ####

library(remotes)
# remotes::install_github("jaredhuling/ordinis")
# remotes::install_github("tomwenseleers/L0glm/L0glm")
# library(L0glm)
devtools::load_all(".")
library(export)
library(Matrix)
library(ordinis)
library(L0Learn) # use at least version 2
library(abess)
library(glmnet)


# 1. SINGLE CHANNEL ####
# Simulate SC n=p=500
n <- p <- 500 # 10000
sim.gaus <- simul(method = 1, n = n, p = p, sparse = T)
sim.pois <- simul(method = 2, n = n, p = p, sparse = T)

## 1.1. L0GLM FIT SINGLE CHANNEL WITH POISSON NOISE ####

## Tune lambda relaxed = F with ebic
system.time(tune.osqp.sc.pois <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "ebic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = FALSE, family = poisson(identity), nonnegative = TRUE, solver = "osqp"))
# bestlambda = 10, time = 2.7s (for n=p=500), bestlambda = 10, time = 71s (for n=p=10000)
# TODO fix <sparse>[ <logic> ] : .M.sub.i.logical() maybe inefficient

system.time(tune.glmnet.sc.pois <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "ebic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = FALSE, family = poisson(identity), nonnegative = TRUE, solver = "glmnet"))
# bestlambda = 10, time = 8.1s, bestlambda = 10, time = 31s (for n=p=10000)
system.time(tune.chol.sc.pois <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "ebic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = FALSE, family = poisson(identity), nonnegative = TRUE, solver = "eigen", method = 7L))
# bestlambda = 10, time = 1.1s, bestlambda = 10, time = 5s (for n=p=10000)
system.time(tune.QR.sc.pois <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "ebic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = FALSE, family = poisson(identity), nonnegative = TRUE, solver = "eigen", method = 1L))
# bestlambda = 10, time = 1.6s, bestlambda = 10, time = 486s (for n=p=10000)
system.time(tune.CG.sc.pois <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "ebic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = FALSE, family = poisson(identity), nonnegative = TRUE, solver = "eigen", method = 5L))
# bestlambda = 10, time = 1.5s, bestlambda = 10, time = 16s (for n=p=10000)
system.time(tune.LSCG.sc.pois <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "ebic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = FALSE, family = poisson(identity), nonnegative = TRUE, solver = "eigen", method = 3L))
# bestlambda = 10, time = 1.9s, bestlambda = 10, time = 48s (for n=p=10000)

#### 1.1.1 L0GLM FIT WITH OSQP SOLVER ####
system.time(L0glm.osqp.pois <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                    family = poisson(identity), 
                                    lambda = "gic", # tune.osqp.sc$best.lambda, # choose preset or tuned
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "osqp")) 
# 0.16s with n=p=500, 5.5s with n=p=10000
cor(as.vector(as.vector(L0glm.osqp.pois$coefficients)), as.vector(sim.pois$beta_true)) # 0.9778754 with n=p=500, 0.9422652 with n=p=10000
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.osqp.pois$coefficients) != 0)) # TP = 45 with n=p=500, 839 with n=p=10000
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.osqp.pois$coefficients) != 0)) # FP = 8 with n=p=500, 173 with n=p=10000
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.osqp.pois$coefficients) == 0)) # TN = 442 with n=p=500, 8827 with n=p=10000
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.osqp.pois$coefficients) == 0)) # FN = 5 with n=p=500, 161 with n=p=10000

#### 1.1.2 L0GLM FIT WITH GLMNET/bigGlm SOLVER ####
system.time(L0glm.glmnet.pois <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                      family = poisson(identity), 
                                      lambda = "gic", # tune.glmnet.sc$best.lambda,
                                      nonnegative = TRUE,
                                      relaxed = FALSE,
                                      solver = "glmnet")) 
# 0.36s with n=p=500, 2.13s with n=p=10000
cor(as.vector(as.vector(L0glm.glmnet.pois$coefficients)), as.vector(sim.pois$beta_true)) # 0.9776087 with n=p=500, 0.938967 with n=p=10000
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.glmnet.pois$coefficients) != 0)) # TP = 45 with n=p=500, 839 with n=p=10000
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.glmnet.pois$coefficients) != 0)) # FP = 8 with n=p=500, 176 with n=p=10000
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.glmnet.pois$coefficients) == 0)) # TN = 442 with n=p=500, 8824 with n=p=10000
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.glmnet.pois$coefficients) == 0)) # FN = 5 with n=p=500, 161 with n=p=10000

#### 1.1.3 L0GLM FIT WITH QR SOLVER ####
system.time(L0glm.QR.pois <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                  family = poisson(identity), 
                                  lambda = "gic", # tune.QR.sc$best.lambda,
                                  nonnegative = TRUE,
                                  relaxed = FALSE,
                                  solver = "eigen",
                                  method = 1L)) 
# 0.07s with n=p=500, 20.0s with n=p=10000
# NOTE: TODO: check why QR solver is now slower than it used to be for n=p=10000
cor(as.vector(as.vector(L0glm.QR.pois$coefficients)), as.vector(sim.pois$beta_true)) # 0.9838164 with n=p=500, 0.9595521 with n=p=10000
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.QR.pois$coefficients) != 0)) # TP = 45 with n=p=500, 841 with n=p=10000
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.QR.pois$coefficients) != 0)) # FP = 7 with n=p=500, 173 with n=p=10000
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.QR.pois$coefficients) == 0)) # TN = 443 with n=p=500, 8827 with n=p=10000
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.QR.pois$coefficients) == 0)) # FN = 5 with n=p=500, 161 with n=p=10000

#### 1.1.4 L0GLM FIT WITH CHOLESKY SOLVER ####
system.time(L0glm.chol.pois <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                    family = poisson(identity), 
                                    lambda = "gic", # tune.chol.sc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "eigen",
                                    method = 7L)) 
# 0.03s with n=p=500, 0.3s with n=p=10000
cor(as.vector(as.vector(L0glm.chol.pois$coefficients)), as.vector(sim.pois$beta_true)) # 0.9838164 with n=p=500, 0.9595521 with n=p=10000
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.chol.pois$coefficients) != 0)) # TP = 45 with n=p=500, 841 with n=p=10000
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.chol.pois$coefficients) != 0)) # FP = 7 with n=p=500, 173 with n=p=10000
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.chol.pois$coefficients) == 0)) # TN = 443 with n=p=500, 8827 with n=p=10000
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.chol.pois$coefficients) == 0)) # FN = 5 with n=p=500, 159 with n=p=10000

#### 1.1.5 L0GLM FIT WITH LSCG SOLVER ####
system.time(L0glm.LSCG.pois <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                    family = poisson(identity), 
                                    lambda = "gic", # tune.LSCG.sc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "eigen",
                                    method = 3L)) 
# 0.16s with n=p=500, 0.94s with n=p=10000
# NOTE: LSCG solver now became faster in latest implementation
cor(as.vector(as.vector(L0glm.LSCG.pois$coefficients)), as.vector(sim.pois$beta_true)) # 0.984469 with n=p=500, 0.9595335 with n=p=10000
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.LSCG.pois$coefficients) != 0)) # TP = 45 with n=p=500, 841 with n=p=10000
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.LSCG.pois$coefficients) != 0)) # FP = 7 with n=p=500, 173 with n=p=10000
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.LSCG.pois$coefficients) == 0)) # TN = 443 with n=p=500, 8827 with n=p=10000
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.LSCG.pois$coefficients) == 0)) # FN = 5 with n=p=500, 159 with n=p=10000

#### 1.1.6 L0GLM FIT WITH CG SOLVER ####
system.time(L0glm.CG.pois <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                  family = poisson(identity), 
                                  lambda = "gic", # tune.CG.sc$best.lambda,
                                  nonnegative = TRUE,
                                  relaxed = FALSE,
                                  solver = "eigen",
                                  method = 5L)) 
# 0.09s with n=p=500, 0.51s with n=p=10000
cor(as.vector(as.vector(L0glm.CG.pois$coefficients)), as.vector(sim.pois$beta_true)) # 0.9845091 with n=p=500, 0.9593445 with n=p=10000
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.CG.pois$coefficients) != 0)) # TP = 45 with n=p=500, 841 with n=p=10000
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.CG.pois$coefficients) != 0)) # FP = 7 with n=p=500, 173 with n=p=10000
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.CG.pois$coefficients) == 0)) # TN = 443 with n=p=500, 8827 with n=p=10000
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.CG.pois$coefficients) == 0)) # FN = 5 with n=p=500, 159 with n=p=10000


#### 1.1.7 L0GLM FIT WITH NNLS SOLVER ####
system.time(L0glm.nnls.pois <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                  family = poisson(identity), 
                                  lambda = "gic", # tune.CG.sc$best.lambda,
                                  nonnegative = TRUE,
                                  relaxed = FALSE,
                                  solver = "nnls")) 
# 0.53s with n=p=500, XXs with n=p=10000
cor(as.vector(as.vector(L0glm.nnls.pois$coefficients)), as.vector(sim.pois$beta_true)) # 0.9778754 with n=p=500, XX with n=p=10000
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.nnls.pois$coefficients) != 0)) # TP = 45 with n=p=500, XX with n=p=10000
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.nnls.pois$coefficients) != 0)) # FP = 8 with n=p=500, XX with n=p=10000
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.nnls.pois$coefficients) == 0)) # TN = 442 with n=p=500, XX with n=p=10000
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.nnls.pois$coefficients) == 0)) # FN = 5 with n=p=500, XX with n=p=10000


#### 1.1.8 PLAIN NNLS FIT WITHOUT REGULARISATION ####
system.time(L0glm.nnlsnoreg.pois <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                    family = poisson(identity), 
                                    lambda = 0,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "nnls",
                                    maxit = 1)) 
# 0.08s with n=p=500, XXs with n=p=10000
cor(as.vector(as.vector(L0glm.nnlsnoreg.pois$coefficients)), as.vector(sim.pois$beta_true)) # 0.9803764 with n=p=500, XX with n=p=10000
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.nnlsnoreg.pois$coefficients) != 0)) # TP = 45 with n=p=500, XX with n=p=10000
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.nnlsnoreg.pois$coefficients) != 0)) # FP = 50 with n=p=500, XX with n=p=10000
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.nnlsnoreg.pois$coefficients) == 0)) # TN = 400 with n=p=500, XX with n=p=10000
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.nnlsnoreg.pois$coefficients) == 0)) # FN = 5 with n=p=500, XX with n=p=10000


### 1.1.9 BENCHMARKS TABLE ####

X = sim.pois$X
y = sim.pois$y
W = 1/(sim.pois$y+0.1)
benchmarks_1_1 = rbind(
  data.frame(method="L0glm_osqp", do.call(cbind, get.metrics(beta = L0glm.osqp.pois$coefficients,  # some metrics
                                                        X = X, y = y, 
                                                        weights = W,
                                                        family = gaussian(identity), 
                                                        beta_true = sim.pois$beta_true))),
  data.frame(method="L0glm_bigglm", do.call(cbind, get.metrics(beta = L0glm.glmnet.pois$coefficients,  # some metrics
                                                                  X = X, y = y, 
                                                                  weights = W,
                                                                  family = gaussian(identity), 
                                                                  beta_true = sim.pois$beta_true))),
  data.frame(method="L0glm_qr", do.call(cbind, get.metrics(beta = L0glm.QR.pois$coefficients,  # some metrics
                                                                                 X = X, y = y, 
                                                                                 weights = W,
                                                                                 family = gaussian(identity),
                                                                                 beta_true = as.vector(sim.pois$beta_true)) )),
  data.frame(method="L0glm_cholesky", do.call(cbind, get.metrics(beta = L0glm.chol.pois$coefficients,  # some metrics
                                                          X = X, y = y, 
                                                          weights = W,
                                                          family = gaussian(identity),
                                                          beta_true = sim.pois$beta_true) )),
  data.frame(method="L0glm_LSCG", do.call(cbind, get.metrics(beta = L0glm.LSCG.pois$coefficients,  # some metrics
                                                                               X = X, y = y, 
                                                                               weights = W,
                                                                               family = gaussian(identity),
                                                                               beta_true = sim.pois$beta_true) )),
  data.frame(method="L0glm_CG", do.call(cbind, get.metrics(beta = L0glm.CG.pois$coefficients,  # some metrics
                                                            X = X, y = y, 
                                                            weights = W,
                                                            family = gaussian(identity), 
                                                            beta_true = sim.pois$beta_true) )),
  data.frame(method="L0glm_nnls", do.call(cbind, get.metrics(beta = L0glm.nnls.pois$coefficients,  # some metrics
                                                                 X = X, y = y, 
                                                                 weights = W,
                                                                 family = gaussian(identity), 
                                                                 beta_true = sim.pois$beta_true) )),
  data.frame(method="nnls_noreg", do.call(cbind, get.metrics(beta = L0glm.nnlsnoreg.pois$coefficients,  # some metrics
                                                                X = X, y = y, 
                                                                weights = W,
                                                                family = gaussian(identity),
                                                                beta_true = sim.pois$beta_true) )) 
)

benchmarks_1_1 = benchmarks_1_1[order(benchmarks_1_1[,"deviance_true"]),] # sort by deviance relative to ground truth
benchmarks_1_1
#           method   p   n m deviance deviance_true    loglik  MD_coefs MAD_coefs sMAPE_coefs
# 6       L0glm_CG 501 500 1 514.2263      70.91080 -1562.391 0.2074357  37.86133  0.05302939
# 5     L0glm_LSCG 501 500 1 514.4095      71.20591 -1562.480 0.2072647  37.86417  0.05302451
# 4 L0glm_cholesky 501 500 1 518.3894      76.91713 -1564.407 0.2075084  38.14443  0.05304871
# 3       L0glm_qr 501 500 1 518.3894      76.91713 -1564.407 0.2075084  38.14443  0.05304871
# 7     L0glm_nnls 501 500 1 525.5174      84.43703 -1568.821 0.2070470  49.54378  0.05814334
# 1     L0glm_osqp 501 500 1 525.5174      84.43705 -1568.821 0.2070470  49.54378  0.05814334
# 2   L0glm_bigglm 501 500 1 525.4374      84.59868 -1568.783 0.2069567  49.79920  0.05817963
# 8     nnls_noreg 501 500 1 443.1191      95.59246 -1568.185 0.1159119  72.72359  0.23426504
#   bias_coefs  R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs
# 6   1.013969 0.9704626   112.75148     253.2943      0.3175903      1.013969   0.9788911
# 5   1.014041 0.9703866   112.74984     253.3217      0.3175433      1.014041   0.9788144
# 4   1.015115 0.9691515   112.75218     256.0164      0.3177761      1.015115   0.9775686
# 3   1.015115 0.9691515   112.75218     256.0164      0.3177761      1.015115   0.9775686
# 7   1.024967 0.9579401   110.62043     358.7269      0.3598428      1.024967   0.9662599
# 1   1.024967 0.9579401   110.62043     358.7270      0.3598428      1.024967   0.9662599
# 2   1.024829 0.9574396   110.61958     361.1366      0.3601852      1.024829   0.9657549
# 8   1.032099 0.9626282    74.87893     308.4869      1.1277107      1.032099   0.9730918
#   support_size_true support_size_estimated estimated_over_true_support_size TP FP  TN FN
# 6                50                     52                             1.04 45  7 443  5
# 5                50                     52                             1.04 45  7 443  5
# 4                50                     52                             1.04 45  7 443  5
# 3                50                     52                             1.04 45  7 443  5
# 7                50                     53                             1.06 45  8 442  5
# 1                50                     53                             1.06 45  8 442  5
# 2                50                     53                             1.06 45  8 442  5
# 8                50                     95                             1.90 45 50 400  5
#   precision recall specificity accuracy  F1_score      gic     ebic      bic      aic     mbic
# 6 0.8653846    0.9   0.9844444    0.976 0.8823529 3726.710 3783.742 3454.156 3230.782 3966.169
# 5 0.8653846    0.9   0.9844444    0.976 0.8823529 3726.888 3783.920 3454.334 3230.960 3966.347
# 4 0.8653846    0.9   0.9844444    0.976 0.8823529 3730.742 3787.774 3458.188 3234.814 3970.201
# 3 0.8653846    0.9   0.9844444    0.976 0.8823529 3730.742 3787.774 3458.188 3234.814 3970.201
# 7 0.8490566    0.9   0.9822222    0.974 0.8737864 3750.927 3809.035 3473.231 3245.642 3994.904
# 1 0.8490566    0.9   0.9822222    0.974 0.8737864 3750.927 3809.035 3473.231 3245.642 3994.904
# 2 0.8490566    0.9   0.9822222    0.974 0.8737864 3750.851 3808.959 3473.155 3245.566 3994.828
# 8 0.4736842    0.9   0.8888889    0.890 0.6206897 4226.655 4329.958 3732.972 3328.370 4660.392




## 1.2. L0GLM FIT SINGLE CHANNEL WITH GAUSSIAN NOISE ####

# TODO: tuning function below returns error
# FIX BUG
# Error in lm.fit(x, y, offset = offset, singular.ok = singular.ok, ...) : 
# 0 (non-NA) cases

## Tune lambda relaxed = F 
system.time(tune.osqp.sc.gaus <- L0glm.tune(X=sim.gaus$X, y=sim.gaus$y, IC = "ebic", plot = FALSE, coefs.true = sim.gaus$beta_true, relaxed = FALSE, family = gaussian(identity), nonnegative = TRUE, solver = "osqp"))
# bestlambda = 3.162277, time = 4.614s
system.time(tune.glmnet.sc.gaus <- L0glm.tune(X=sim.gaus$X, y=sim.gaus$y, IC = "ebic", plot = FALSE, coefs.true = sim.gaus$beta_true, relaxed = FALSE, family = gaussian(identity), nonnegative = TRUE, solver = "glmnet"))
# bestlambda = 316.2277, time = 5.411s
system.time(tune.chol.sc.gaus <- L0glm.tune(X=sim.gaus$X, y=sim.gaus$y, IC = "ebic", plot = FALSE, coefs.true = sim.gaus$beta_true, relaxed = FALSE, family = gaussian(identity), nonnegative = TRUE, solver = "eigen", method = 7L))
# bestlambda = 3.162277, time = 3.052s
system.time(tune.QR.sc.gaus <- L0glm.tune(X=sim.gaus$X, y=sim.gaus$y, IC = "ebic", plot = FALSE, coefs.true = sim.gaus$beta_true, relaxed = TRUE, family = gaussian(identity), nonnegative = TRUE, solver = "eigen", method = 1L))
# bestlambda = 3.162277, time = 15.925s
system.time(tune.CG.sc.gaus <- L0glm.tune(X=sim.gaus$X, y=sim.gaus$y, IC = "ebic", plot = FALSE, coefs.true = sim.gaus$beta_true, relaxed = FALSE, family = gaussian(identity), nonnegative = TRUE, solver = "eigen", method = 5L))
# bestlambda = 3.162277, time = 7.951s
system.time(tune.LSCG.sc.gaus <- L0glm.tune(X=sim.gaus$X, y=sim.gaus$y, IC = "ebic", plot = TRUE, coefs.true = sim.gaus$beta_true, relaxed = FALSE, family = gaussian(identity), nonnegative = TRUE, solver = "eigen", method = 3L))
# bestlambda = 3.162277, time = 10.929s

#### 1.2.1 L0GLM FIT WITH OSQP SOLVER ####
system.time(L0glm.osqp.gaus <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                    family = gaussian(identity), 
                                    lambda = "gic", # tune.osqp.sc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "osqp")) 
# 0.2s for n=p=500, 42s for n=p=10000
cor(as.vector(as.vector(L0glm.osqp.gaus$coefficients)), as.vector(sim.gaus$beta_true)) #  0.9836709
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.osqp.gaus$coefficients) != 0)) # TP = 48
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.osqp.gaus$coefficients) != 0)) # FP = 7
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.osqp.gaus$coefficients) == 0)) # TN = 443
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.osqp.gaus$coefficients) == 0)) # FN = 2

#### 1.2.2 L0GLM FIT WITH GLMNET/bigGlm SOLVER ####
system.time(L0glm.glmnet.gaus <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                      family = gaussian(identity), 
                                      lambda = "gic", # tune.glmnet.sc$best.lambda,
                                      nonnegative = TRUE,
                                      relaxed = FALSE,
                                      solver = "glmnet")) 
# 1.44s for n=p=500, 5.6s for n=p=10000 
cor(as.vector(as.vector(L0glm.glmnet.gaus$coefficients)), as.vector(sim.gaus$beta_true)) # 0.9637935
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.glmnet.gaus$coefficients) != 0)) # TP = 47
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.glmnet.gaus$coefficients) != 0)) # FP = 24
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.glmnet.gaus$coefficients) == 0)) # TN = 426
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.glmnet.gaus$coefficients) == 0)) # FN = 3

#### 1.2.3 L0GLM FIT WITH QR SOLVER ####
system.time(L0glm.QR.gaus <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                  family = gaussian(identity), 
                                  lambda = "gic", # tune.QR.sc$best.lambda,
                                  nonnegative = TRUE,
                                  relaxed = FALSE,
                                  solver = "eigen",
                                  method = 1L)) 
# 0.04s for n=p=500, 143s for n=p=10000
cor(as.vector(as.vector(L0glm.QR.gaus$coefficients)), as.vector(sim.gaus$beta_true)) # 0.997474
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.QR.gaus$coefficients) != 0)) # TP = 49
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.QR.gaus$coefficients) != 0)) # FP = 4
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.QR.gaus$coefficients) == 0)) # TN = 446
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.QR.gaus$coefficients) == 0)) # FN = 1

#### 1.2.4 L0GLM FIT WITH Cholesky SOLVER ####
system.time(L0glm.chol.gaus <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                    family = gaussian(identity), 
                                    lambda = "gic", # tune.chol.sc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "eigen",
                                    method = 7L)) 
# 0.02s for n=p=500, 0.25s for n=p=10000
cor(as.vector(as.vector(L0glm.chol.gaus$coefficients)), as.vector(sim.gaus$beta_true)) # 0.997474
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.chol.gaus$coefficients) != 0)) # TP = 49
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.chol.gaus$coefficients) != 0)) # FP = 4
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.chol.gaus$coefficients) == 0)) # TN = 446
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.chol.gaus$coefficients) == 0)) # FN = 1

#### 1.2.5 L0GLM FIT WITH LSCG SOLVER ####
system.time(L0glm.LSCG.gaus <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                    family = gaussian(identity), 
                                    lambda = "gic", # tune.CG.sc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "eigen",
                                    method = 4L)) 
# 0.13s for n=p=500, 8.5s for n=p=10000
cor(as.vector(as.vector(L0glm.LSCG.gaus$coefficients)), as.vector(sim.gaus$beta_true)) # 0.9777964
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.LSCG.gaus$coefficients) != 0)) # TP = 47
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.LSCG.gaus$coefficients) != 0)) # FP = 12
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.LSCG.gaus$coefficients) == 0)) # TN = 438
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.LSCG.gaus$coefficients) == 0)) # FN = 3

#### 1.2.6 L0GLM FIT WITH CG SOLVER ####
system.time(L0glm.CG.gaus <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                  family = gaussian(identity), 
                                  lambda = "gic", # tune.CG.sc$best.lambda,
                                  nonnegative = TRUE,
                                  relaxed = TRUE,
                                  solver = "eigen",
                                  method = 5L)) 
# 0.06s for n=p=500, 3.7s for n=p=10000
cor(as.vector(as.vector(L0glm.CG.gaus$coefficients)), as.vector(sim.gaus$beta_true)) # 0.9871489
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.CG.gaus$coefficients) != 0)) # TP = 46
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.CG.gaus$coefficients) != 0)) # FP = 7
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.CG.gaus$coefficients) == 0)) # TN = 443
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.CG.gaus$coefficients) == 0)) # FN = 4


#### 1.2.7 L0GLM FIT WITH NNLS SOLVER ####
system.time(L0glm.nnls.gaus <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                    family = gaussian(identity), 
                                    lambda = "gic", # tune.CG.sc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "nnls")) 
# 0.38s with n=p=500, XXs with n=p=10000
cor(as.vector(as.vector(L0glm.nnls.gaus$coefficients)), as.vector(sim.pois$beta_true)) # 0.9871489 with n=p=500, XX with n=p=10000
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.nnls.gaus$coefficients) != 0)) # TP = 48 with n=p=500, XX with n=p=10000
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.nnls.gaus$coefficients) != 0)) # FP = 7 with n=p=500, XX with n=p=10000
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.nnls.gaus$coefficients) == 0)) # TN = 443 with n=p=500, XX with n=p=10000
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.nnls.gaus$coefficients) == 0)) # FN = 2 with n=p=500, XX with n=p=10000


#### 1.2.8 PLAIN NNLS FIT WITHOUT REGULARISATION ####
system.time(L0glm.nnlsnoreg.gaus <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                         family = gaussian(identity), 
                                         lambda = 0,
                                         nonnegative = TRUE,
                                         relaxed = FALSE,
                                         solver = "nnls",
                                         maxit = 1)) 
# 0.1s with n=p=500, XXs with n=p=10000
cor(as.vector(as.vector(L0glm.nnlsnoreg.gaus$coefficients)), as.vector(sim.pois$beta_true)) # 0.9972479 with n=p=500, XX with n=p=10000
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.nnlsnoreg.gaus$coefficients) != 0)) # TP = 49 with n=p=500, XX with n=p=10000
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.nnlsnoreg.gaus$coefficients) != 0)) # FP = 77 with n=p=500, XX with n=p=10000
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.nnlsnoreg.gaus$coefficients) == 0)) # TN = 373 with n=p=500, XX with n=p=10000
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.nnlsnoreg.gaus$coefficients) == 0)) # FN = 1 with n=p=500, XX with n=p=10000



### 1.2.9 BENCHMARKS TABLE ####

X = sim.gaus$X
y = sim.gaus$y
W = rep(1, nrow(X))
benchmarks_1_2 = rbind(
  data.frame(method="L0glm_osqp", do.call(cbind, get.metrics(beta = L0glm.osqp.gaus$coefficients,  # some metrics
                                                             X = X, y = y, 
                                                             weights = W,
                                                             family = gaussian(identity), 
                                                             beta_true = sim.gaus$beta_true))),
  data.frame(method="L0glm_bigglm", do.call(cbind, get.metrics(beta = L0glm.glmnet.gaus$coefficients,  # some metrics
                                                               X = X, y = y, 
                                                               weights = W,
                                                               family = gaussian(identity), 
                                                               beta_true = sim.gaus$beta_true))),
  data.frame(method="L0glm_qr", do.call(cbind, get.metrics(beta = L0glm.QR.gaus$coefficients,  # some metrics
                                                           X = X, y = y, 
                                                           weights = W,
                                                           family = gaussian(identity),
                                                           beta_true = as.vector(sim.gaus$beta_true)) )),
  data.frame(method="L0glm_cholesky", do.call(cbind, get.metrics(beta = L0glm.chol.gaus$coefficients,  # some metrics
                                                                 X = X, y = y, 
                                                                 weights = W,
                                                                 family = gaussian(identity),
                                                                 beta_true = sim.gaus$beta_true) )),
  data.frame(method="L0glm_LSCG", do.call(cbind, get.metrics(beta = L0glm.LSCG.gaus$coefficients,  # some metrics
                                                             X = X, y = y, 
                                                             weights = W,
                                                             family = gaussian(identity),
                                                             beta_true = sim.gaus$beta_true) )),
  data.frame(method="L0glm_CG", do.call(cbind, get.metrics(beta = L0glm.CG.gaus$coefficients,  # some metrics
                                                           X = X, y = y, 
                                                           weights = W,
                                                           family = gaussian(identity), 
                                                           beta_true = sim.gaus$beta_true) )),
  data.frame(method="L0glm_nnls", do.call(cbind, get.metrics(beta = L0glm.nnls.gaus$coefficients,  # some metrics
                                                             X = X, y = y, 
                                                             weights = W,
                                                             family = gaussian(identity), 
                                                             beta_true = sim.gaus$beta_true) )),
  data.frame(method="nnls_noreg", do.call(cbind, get.metrics(beta = L0glm.nnlsnoreg.gaus$coefficients,  # some metrics
                                                             X = X, y = y, 
                                                             weights = W,
                                                             family = gaussian(identity),
                                                             beta_true = sim.gaus$beta_true) )) 
)

benchmarks_1_2 = benchmarks_1_2[order(benchmarks_1_2[,"deviance_true"]),] # sort by deviance relative to ground truth
benchmarks_1_2
#         method   p   n m  deviance deviance_true    loglik   MD_coefs MAD_coefs sMAPE_coefs
# 3       L0glm_qr 501 500 1 10063.744      1649.103 -1012.991 0.01229024  9.014314  0.02220367
# 4 L0glm_cholesky 501 500 1 10062.394      1650.157 -1012.957 0.01228269  9.014960  0.02220381
# 7     L0glm_nnls 501 500 1 10918.777      3177.996 -1035.377 0.02047878 26.729453  0.03747524
# 1     L0glm_osqp 501 500 1 10985.337      3409.879 -1036.896 0.02230046 32.811943  0.03810321
# 8     nnls_noreg 501 500 1  8040.643      3439.262 -1029.883 0.10402808 17.717465  0.31621381
# 5     L0glm_LSCG 501 500 1 14506.473      6307.089 -1110.405 0.03529220 41.519060  0.06255280
# 2   L0glm_bigglm 501 500 1 16220.289      8284.284 -1150.322 0.05986969 82.259515  0.12272163
# 6       L0glm_CG 501 500 1 61265.073     57983.646 -1464.559 0.10564020 31.709166  0.04516657
#   bias_coefs  R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs
# 3  0.9974662 0.9951533   16.819839     68.33680      0.1717328     0.9974662   0.9958997
# 4  0.9974666 0.9951533   16.819768     68.34290      0.1717341     0.9974666   0.9958996
# 7  1.0429526 0.9752123   54.385742    188.79545      0.2679568     1.0429526   0.9800624
# 1  1.0461334 0.9687987   54.402303    244.09082      0.2736655     1.0461334   0.9736169
# 8  1.0113971 0.9947146    7.439051     63.28116      1.2389437     1.0113971   0.9954606
# 5  1.0605777 0.9576665   54.491444    297.66408      0.4284135     1.0605777   0.9624727
# 2  1.0680164 0.9315111   40.781697    538.93228      0.7797298     1.0680164   0.9357567
# 6  1.0095223 0.9754668  104.894522    195.24516      0.2751564     1.0095223   0.9854120
#   support_size_true support_size_estimated estimated_over_true_support_size TP FP  TN FN
# 3                50                     53                             1.06 49  4 446  1
# 4                50                     53                             1.06 49  4 446  1
# 7                50                     55                             1.10 48  7 443  2
# 1                50                     55                             1.10 48  7 443  2
# 8                50                    126                             2.52 49 77 373  1
# 5                50                     59                             1.18 47 12 438  3
# 2                50                     71                             1.42 47 24 426  3
# 6                50                     53                             1.06 46  7 443  4
#   precision recall specificity accuracy  F1_score      gic     ebic      bic      aic     mbic
# 3 0.9245283   0.98   0.9911111    0.990 0.9514563 2639.267 2697.375 2361.571 2133.982 2883.244
# 4 0.9245283   0.98   0.9911111    0.990 0.9514563 2639.200 2697.308 2361.503 2133.915 2883.177
# 7 0.8727273   0.96   0.9844444    0.982 0.9142857 2706.754 2767.014 2418.772 2182.754 2959.767
# 1 0.8727273   0.96   0.9844444    0.982 0.9142857 2709.792 2770.053 2421.811 2185.793 2962.806
# 8 0.3888889   0.98   0.8288889    0.844 0.5568182 3502.123 3638.785 2849.022 2313.767 4075.921
# 5 0.7966102   0.94   0.9733333    0.970 0.8623853 2902.238 2966.802 2593.686 2340.810 3173.323
# 2 0.6619718   0.94   0.9466667    0.946 0.7768595 3118.357 3195.835 2748.095 2444.644 3443.660
# 6 0.8679245   0.92   0.9844444    0.978 0.8932039 3542.402 3600.511 3264.706 3037.117 3786.380




# 2. MULTICHANNEL ####
# Simulate MC n=p=500
n <- p <- 500
sim.gaus.mc <- simul(method = 3, n = n, p = p, seed = 1)
sim.pois.mc <- simul(method = 4, n = n, p = p, seed = 2)



## 2.1. MULTICHANNEL WITH POISSON NOISE ####

### 2.1.1 GROUPED L0GLM FIT WITH relaxed = FALSE, grouped = FALSE ####

## Tune lambda relaxed = FALSE, grouped = FALSE
system.time(tune.osqp.mc <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, IC = "ebic", 
                                       plot = FALSE, coefs.true = sim.pois.mc$beta_true, 
                                       relaxed = FALSE, family = poisson(identity), 
                                       nonnegative = TRUE, 
                                       parallelize = TRUE, grouped = FALSE, solver = "osqp"))
# bestlambda = 3.16227766016838, time = 127s

system.time(tune.chol.mc <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, IC = "ebic", 
                                       plot = FALSE, coefs.true = sim.pois.mc$beta_true, 
                                       relaxed = FALSE, family = poisson(identity), 
                                       nonnegative = TRUE, 
                                       parallelize = TRUE, grouped = FALSE,
                                       solver = "eigen", method = 7L))
# bestlambda = 3.16227766016838, time = 113s


#### OSQP ####
system.time(L0glm.osqp <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                    family = poisson(identity), 
                                    lambda = "ebic", # tune.osqp.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallelize = TRUE,
                                    grouped = FALSE,
                                    solver = "osqp")) 
# 7.1s
# abline(v=which(rowMin(L0glm.osqp$coefficients)!=0), col = 'cyan')

cor(as.vector(as.vector(L0glm.osqp$coefficients)), as.vector(sim.pois.mc$beta_true)) # 0.9941838
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)!=0) # TP = 50
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)!=0) # FP = 169
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)==0) # TN = 281
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)==0) # FN = 0

#### Cholesky ####
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
# abline(v=which(rowMin(L0glm.chol$coefficients)!=0), col = 'cyan')

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


### 2.1.2 GROUPED L0GLM FIT WITH relaxed = FALSE, grouped = TRUE ####

## Tune lambda relaxed = FALSE, grouped = TRUE
system.time(tune.osqp.mc <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, IC = "ebic", 
                                       plot = FALSE, coefs.true = sim.pois.mc$beta_true, 
                                       relaxed = FALSE, family = poisson(identity), 
                                       nonnegative = TRUE, 
                                       parallelize = TRUE, grouped = TRUE, solver = "osqp"))

# bestlambda = 3.16227766016838, time = 253s
system.time(tune.chol.mc <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, IC = "ebic", 
                                       plot = FALSE, coefs.true = sim.pois.mc$beta_true, 
                                       relaxed = FALSE, family = poisson(identity), 
                                       nonnegative = TRUE, 
                                       parallelize = TRUE, grouped = TRUE,
                                       solver = "eigen", method = 7L))
# bestlambda = 3.16227766016838, time = 190s

#### OSQP ####
system.time(L0glm.osqp <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                    family = poisson(identity), 
                                    lambda = "ebic", # tune.osqp.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallelize = TRUE,
                                    grouped = TRUE,
                                    solver = "osqp"))
# 14s
# abline(v=which(rowMin(L0glm.osqp$coefficients)!=0), col = 'cyan')

cor(as.vector(as.vector(L0glm.osqp$coefficients)), as.vector(sim.pois.mc$beta_true)) # 0.9992435
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)!=0) # TP = 49
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)!=0) # FP = 1
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)==0) # TN = 449
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)==0) # FN = 1

#### Cholesky ####
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
# abline(v=which(rowMin(L0glm.chol$coefficients)!=0), col = 'cyan')

cor(as.vector(as.vector(L0glm.chol$coefficients)), as.vector(sim.pois.mc$beta_true)) # 0.9992436
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.chol$coefficients)!=0) # TP = 49
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.chol$coefficients)!=0) # FP = 1
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.chol$coefficients)==0) # TN = 449
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.chol$coefficients)==0) # FN = 1


## 2.1.3 GROUPED L0GLM FIT WITH relaxed = TRUE, grouped = TRUE ####

# TODO: this should fit the subsequent iterations using grouped=TRUE & relaxed=FALSE
# & only as an extra iteration then refit using grouped=TRUE & relaxed=TRUE
# check that this is how it is done now

## Tune lambda relaxed = TRUE, grouped = TRUE
system.time(tune.osqp.mc <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, IC = "ebic", 
                                       plot = FALSE, coefs.true = sim.pois.mc$beta_true, 
                                       relaxed = TRUE, family = poisson(identity), 
                                       nonnegative = TRUE, 
                                       parallelize = TRUE, grouped = TRUE, solver = "osqp"))

# bestlambda = 3.16227766016838, time = 275s
system.time(tune.chol.mc <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, IC = "ebic", 
                                       plot = FALSE, coefs.true = sim.pois.mc$beta_true, 
                                       relaxed = TRUE, family = poisson(identity), 
                                       nonnegative = TRUE, 
                                       parallelize = TRUE, grouped = TRUE,
                                       solver = "eigen", method = 7L))
# bestlambda = 3.16227766016838, time = 206s

## Fit relaxed = TRUE, grouped = TRUE
#### OSQP ####
system.time(L0glm.osqp <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                    family = poisson(identity), 
                                    lambda = "ebic", # tune.osqp.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = TRUE,
                                    parallelize = TRUE,
                                    grouped = TRUE,
                                    solver = "osqp")) 
# 15s
# abline(v=which(rowMin(L0glm.osqp$coefficients)!=0), col = 'cyan')

cor(as.vector(as.vector(L0glm.osqp$coefficients)), as.vector(sim.pois.mc$beta_true)) # 0.9995664
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)!=0) # TP = 49
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)!=0) # FP = 1
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)==0) # TN = 449
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)==0) # FN = 1

#### Cholesky ####
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
# abline(v=which(rowMin(L0glm.chol$coefficients)!=0), col = 'cyan')

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


## 2.2. MULTICHANNEL WITH GAUSSIAN NOISE ####

### 2.2.1 GROUPED L0GLM FIT WITH relaxed = FALSE, grouped = FALSE ####

# TODO fix bugs that appear below
# Error in lm.fit(x, y, offset = offset, singular.ok = singular.ok, ...) : 
# 0 (non-NA) cases
# In addition: There were 17 warnings (use warnings() to see them)
# probably something going wrong with passing the dispersion

## Tune lambda relaxed = FALSE, grouped = FALSE
system.time(tune.osqp.mc <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, IC = "ebic", 
                                       plot = FALSE, coefs.true = sim.gaus.mc$beta_true, 
                                       relaxed = FALSE, family = gaussian(identity), 
                                       nonnegative = TRUE, 
                                       parallelize = TRUE, grouped = FALSE, solver = "osqp"))
# bestlambda = XXX, time = XXXs
# TODO fix bug
system.time(tune.chol.mc <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, IC = "ebic", 
                                       plot = FALSE, coefs.true = sim.gaus.mc$beta_true, 
                                       relaxed = FALSE, family = gaussian(identity), 
                                       nonnegative = TRUE, 
                                       parallelize = TRUE, grouped = FALSE,
                                       solver = "eigen", method = 7L))
# bestlambda = XXX, time = XXXs
# TODO fix bug

#### OSQP ####
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

# abline(v=which(rowMin(L0glm.osqp$coefficients)!=0), col = 'cyan')

cor(as.vector(as.vector(L0glm.osqp$coefficients)), as.vector(sim.gaus.mc$beta_true)) # 0.9898752
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)!=0) # TP = 50
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)!=0) # FP = 409
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)==0) # TN = 41
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)==0) # FN = 0

#### Cholesky ####
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

# abline(v=which(rowMin(L0glm.chol$coefficients)!=0), col = 'cyan')

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


### 2.2.2 GROUPED L0GLM FIT WITH relaxed = FALSE, grouped = TRUE ####

## Tune lambda relaxed = FALSE, grouped = TRUE
system.time(tune.osqp.mc <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, IC = "ebic", 
                                       plot = TRUE, coefs.true = sim.gaus.mc$beta_true, 
                                       relaxed = FALSE, family = gaussian(identity), 
                                       nonnegative = TRUE, 
                                       parallelize = TRUE, grouped = TRUE, solver = "osqp"))
# bestlambda = XXX, time = XXXs
system.time(tune.chol.mc <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, IC = "ebic", 
                                       plot = FALSE, coefs.true = sim.gaus.mc$beta_true, 
                                       relaxed = FALSE, family = gaussian(identity), 
                                       nonnegative = TRUE, 
                                       parallelize = TRUE, grouped = TRUE,
                                       solver = "eigen", method = 7L))
# bestlambda = XXX, time = XXX 


#### OSQP ####
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
# abline(v=which(rowMin(L0glm.osqp$coefficients)!=0), col = 'cyan')

cor(as.vector(as.vector(L0glm.osqp$coefficients)), as.vector(sim.gaus.mc$beta_true)) # 0.9999461
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)!=0) # TP = 50
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)!=0) # FP = 0
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)==0) # TN = 450
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)==0) # FN = 0

#### Cholesky ####
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
# abline(v=which(rowMin(L0glm.chol$coefficients)!=0), col = 'cyan')

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

### 2.2.3 GROUPED L0GLM FIT WITH relaxed = TRUE, grouped = TRUE ####

## Tune lambda relaxed = TRUE, grouped = TRUE
system.time(tune.osqp.mc <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, IC = "ebic", 
                                       plot = TRUE, coefs.true = sim.gaus.mc$beta_true, 
                                       relaxed = TRUE, family = gaussian(identity), 
                                       nonnegative = TRUE, 
                                       parallelize = TRUE, grouped = TRUE, solver = "osqp"))

# bestlambda = XX, time = XX
system.time(tune.chol.mc <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, IC = "ebic", 
                                       plot = FALSE, coefs.true = sim.gaus.mc$beta_true, 
                                       relaxed = TRUE, family = gaussian(identity), 
                                       nonnegative = TRUE, 
                                       parallelize = TRUE, grouped = TRUE,
                                       solver = "eigen", method = 7L))
# bestlambda = XX, time = XX


#### OSQP ####
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
# abline(v=which(rowMin(L0glm.osqp$coefficients)!=0), col = 'cyan')

cor(as.vector(as.vector(L0glm.osqp$coefficients)), as.vector(sim.gaus.mc$beta_true)) # 0.9999514
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)!=0) # TP = 50
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)!=0) # FP = 0
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)==0) # TN = 450
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)==0) # FN = 0


#### Cholesky ####
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
# abline(v=which(rowMin(L0glm.chol$coefficients)!=0), col = 'cyan')

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

## 3.1 SINGLE CHANNEL POISSON ####
n <- p <- 500
sim.pois <- simul(method = 2, n = n, p = p, seed = 2)

### 3.1.1 L0GLM VIA CHOLESKY ####
system.time(tune.chol.sc <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "ebic", 
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
# 0.02s
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
     main="L0 glm identity link Poisson with lambda tuned to optimise ebic",
     xlab="Estimated coefficients", ylab="True coefficients")

### 3.1.2 LASSO USING GLMNET ####
system.time(glmnet.sc.pois <- fit.GLM(sim = sim.pois, method = 2, 
                                      family = poisson(identity), plot = TRUE)) 
# 0.19s

cor2 <- cor(as.vector(glmnet.sc.pois[-1]), as.vector(sim.pois$beta_true)) # 0.9203877
cor2
TP2 <- sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(glmnet.sc.pois[-1]) != 0)) # TP = 48
TP2
FP2 <- sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(glmnet.sc.pois[-1]) != 0)) # FP = 81
FP2
TN2 <- sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(glmnet.sc.pois[-1]) == 0)) # TN = 369
TN2
FN2 <- sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(glmnet.sc.pois[-1]) == 0)) # FN = 2
FN2

### 3.1.3 ADAPTIVE LASSO USING GLMNET ####
system.time(glmnet.al.sc.pois <- fit.GLM(sim = sim.pois, method = 3, 
                                         family = poisson(identity), plot = T)) 
# 0.3s
cor3 <- cor(as.vector(glmnet.al.sc.pois[-1]), as.vector(sim.pois$beta_true)) # 0.9768725
cor3
TP3 <- sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(glmnet.al.sc.pois[-1]) != 0)) # TP = 42
TP3
FP3 <- sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(glmnet.al.sc.pois[-1]) != 0)) # FP = 13
FP3
TN3 <- sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(glmnet.al.sc.pois[-1]) == 0)) # TN = 437
TN3
FN3 <- sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(glmnet.al.sc.pois[-1]) == 0)) # FN = 8
FN3

### 3.1.4 L0LEARN ####
system.time(L0learn.sc.pois <- fit.GLM(sim = sim.pois, method = 4, 
                                       family = poisson(identity), plot = T))  
# 0.90s

cor(as.vector(L0learn.sc.pois), as.vector(sim.pois$beta_true)) # 0.5835151
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0learn.sc.pois) != 0)) # TP = 13
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0learn.sc.pois) != 0)) # FP = 8
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0learn.sc.pois) == 0)) # TN = 442
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0learn.sc.pois) == 0)) # FN = 37


### 3.1.5 L0LEARN AFTER LASSO ####
# Slightly better than regular L0Learn & faster
system.time(L0learn.after.lasso.sc.pois <- fit.GLM(sim = sim.pois, method = 5, family = poisson(identity), plot = T)) 
# 0.64s
cor(as.vector(L0learn.after.lasso.sc.pois), as.vector(sim.pois$beta_true)) # 0.578168
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0learn.after.lasso.sc.pois) != 0)) # TP = 10
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0learn.after.lasso.sc.pois) != 0)) # FP = 5
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0learn.after.lasso.sc.pois) == 0)) # TN = 445
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0learn.after.lasso.sc.pois) == 0)) # FN = 40


### 3.1.6 ABESS ####
system.time(abess.sc.pois <- fit.GLM(sim = sim.pois, method = 6, 
                                     family = poisson(identity), 
                                     plot = T, IC = "ebic"))  
# 1.0s
cor(as.vector(abess.sc.pois), as.vector(sim.pois$beta_true)) # 0.3533756
# TODO: NOTE: result is better if I do not use weights - then cor 0.8960064 - CHECK
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(abess.sc.pois) != 0)) # TP = 11
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(abess.sc.pois) != 0)) # FP = 15
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(abess.sc.pois) == 0)) # TN = 435
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(abess.sc.pois) == 0)) # FN = 39


### 3.1.7 ABESS AFTER LASSO ####
system.time(abess.after.lasso.sc.pois <- fit.GLM(sim = sim.pois, method = 7, 
                                                 family = poisson(identity), 
                                                 plot = T, 
                                                 IC = "ebic"))  
# 0.19s
cor(as.vector(abess.after.lasso.sc.pois), as.vector(sim.pois$beta_true))  # 0.3527641
# TODO: NOTE: result is better if I do not use weights - then cor 0.8849495 - CHECK
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(abess.after.lasso.sc.pois) != 0)) # TP = 14
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(abess.after.lasso.sc.pois) != 0)) # FP = 13
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(abess.after.lasso.sc.pois) == 0)) # TN = 437
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(abess.after.lasso.sc.pois) == 0)) # FN = 36

### 3.1.8 MCP USING ORDINIS ####
system.time(ordinis.mcp.sc.pois <- fit.GLM(sim = sim.pois, method = 8, 
                                           family = poisson(identity), plot = T,
                                           IC = "ebic"))  
# 11s
cor(as.vector(ordinis.mcp.sc.pois[-1]), as.vector(sim.pois$beta_true)) # 0.7109204
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(ordinis.mcp.sc.pois[-1]) != 0)) # TP = 47
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(ordinis.mcp.sc.pois[-1]) != 0)) # FP = 85
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(ordinis.mcp.sc.pois[-1]) == 0)) # TN = 365
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(ordinis.mcp.sc.pois[-1]) == 0)) # FN = 3

### 3.1.9 SCAD USING ORDINIS ####
system.time(ordinis.scad.sc.pois <- fit.GLM(sim = sim.pois, method = 9, 
                                            family = poisson(identity), 
                                            plot = T, IC = "ebic"))  
# 10.5s
cor(as.vector(ordinis.scad.sc.pois[-1]), as.vector(sim.pois$beta_true)) # 0.7141376
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(ordinis.scad.sc.pois[-1]) != 0)) # TP = 47
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(ordinis.scad.sc.pois[-1]) != 0)) # FP = 85
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(ordinis.scad.sc.pois[-1]) == 0)) # TN = 365
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(ordinis.scad.sc.pois[-1]) == 0)) # FN = 3


X <- sim.pois$X
y <- sim.pois$y
W <- 1/(y+0.1)
IC <- "ebic"


### 3.1.10 BENCHMARKS TABLE ####

benchmarks_3_1 = rbind(
  data.frame(method="L0glm", do.call(cbind, get.metrics(beta = L0glm.chol$coefficients,  # some metrics
                                                                     X = X, y = y, 
                                                                     weights = W,
                                                                     family = gaussian(identity), 
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

benchmarks_3_1 = benchmarks_3_1[order(benchmarks_3_1[,"deviance_true"]),] # sort by deviance relative to ground truth
benchmarks_3_1
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


## 3.2 SINGLE CHANNEL GAUSSIAN ####
n <- p <- 500
sim.gaus <- simul(method = 1, n = n, p = p, seed = 2)

### 3.2.1 L0GLM VIA CHOLESKY ####

system.time(tune.chol.sc <- L0glm.tune(X=sim.gaus$X, y=sim.gaus$y, IC = "ebic", 
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
     main="L0 glm identity link gaussian with lambda tuned to optimise ebic",
     xlab="Estimated coefficients", ylab="True coefficients")


### 3.2.2 LASSO USING GLMNET ####
system.time(glmnet.sc.gaus <- fit.GLM(sim = sim.gaus, method = 2, 
                                      family = gaussian(identity), plot = TRUE)) 
# 0.19s 
cor(as.vector(glmnet.sc.gaus[-1]), as.vector(sim.gaus$beta_true)) # 0.9611292
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(glmnet.sc.gaus[-1]) != 0)) # TP = 49
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(glmnet.sc.gaus[-1]) != 0)) # FP = 68
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(glmnet.sc.gaus[-1]) == 0)) # TN = 382
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(glmnet.sc.gaus[-1]) == 0)) # FN = 1


### 3.2.3 ADAPTIVE LASSO USING GLMNET ####
system.time(glmnet.al.sc.gaus <- fit.GLM(sim = sim.gaus, method = 3, 
                                         family = gaussian(identity), plot = T)) 
# 0.22s 
cor(as.vector(glmnet.al.sc.gaus[-1]), as.vector(sim.gaus$beta_true)) # 0.9626742
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(glmnet.al.sc.gaus[-1]) != 0)) # TP = 37
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(glmnet.al.sc.gaus[-1]) != 0)) # FP = 12
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(glmnet.al.sc.gaus[-1]) == 0)) # TN = 438
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(glmnet.al.sc.gaus[-1]) == 0)) # FN = 13


### 3.2.4 L0LEARN ####

system.time(L0learn.sc.gaus <- fit.GLM(sim = sim.gaus, method = 4, 
                                       family = gaussian(identity), plot = T))  
#  1.02s
cor(as.vector(L0learn.sc.gaus), as.vector(sim.gaus$beta_true)) # 0.7419012
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0learn.sc.gaus) != 0)) # TP = 30
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0learn.sc.gaus) != 0)) # FP = 38
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0learn.sc.gaus) == 0)) # TN = 412
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0learn.sc.gaus) == 0)) # FN = 20


### 3.2.5 L0LEARN AFTER LASSO ####
# Better than regular L0Learn & faster
system.time(L0learn.after.lasso.sc.gaus <- fit.GLM(sim = sim.gaus, method = 5, 
                                                   family = gaussian(identity), plot = T))  
# 0.55s 
cor(as.vector(L0learn.after.lasso.sc.gaus), as.vector(sim.gaus$beta_true)) # 0.8748054
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0learn.after.lasso.sc.gaus) != 0)) # TP = 40
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0learn.after.lasso.sc.gaus) != 0)) # FP = 23
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0learn.after.lasso.sc.gaus) == 0)) # TN = 427
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0learn.after.lasso.sc.gaus) == 0)) # FN = 10


### 3.2.6 ABESS ####
system.time(abess.sc.gaus <- fit.GLM(sim = sim.gaus, method = 6, 
                                     family = gaussian(identity), 
                                     plot = T, IC = "ebic"))  
# 1.01s 
cor(as.vector(abess.sc.gaus), as.vector(sim.gaus$beta_true)) # 0.8706258
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(abess.sc.gaus) != 0)) # TP = 28
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(abess.sc.gaus) != 0)) # FP = 20
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(abess.sc.gaus) == 0)) # TN = 430
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(abess.sc.gaus) == 0)) # FN = 22


### 3.2.7 ABESS AFTER LASSO ####
system.time(abess.after.lasso.sc.gaus <- fit.GLM(sim = sim.gaus, method = 7, 
                                                 family = gaussian(identity), 
                                                 plot = T, IC = "ebic"))  
# 0.19s
cor(as.vector(abess.after.lasso.sc.gaus), as.vector(sim.gaus$beta_true))  # 0.8814479
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(abess.after.lasso.sc.gaus) != 0)) # TP = 29
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(abess.after.lasso.sc.gaus) != 0)) # FP = 19
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(abess.after.lasso.sc.gaus) == 0)) # TN = 431
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(abess.after.lasso.sc.gaus) == 0)) # FN = 21


### 3.2.8 MCP USING ORDINIS ####
system.time(ordinis.mcp.sc.gaus <- fit.GLM(sim = sim.gaus, method = 8, 
                                           family = gaussian(identity), 
                                           plot = T, IC = "ebic"))  
# 0.30s
cor(as.vector(ordinis.mcp.sc.gaus[-1]), as.vector(sim.gaus$beta_true)) # 0.7612053
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(ordinis.mcp.sc.gaus[-1]) != 0)) # TP = 24
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(ordinis.mcp.sc.gaus[-1]) != 0)) # FP = 41
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(ordinis.mcp.sc.gaus[-1]) == 0)) # TN = 409
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(ordinis.mcp.sc.gaus[-1]) == 0)) # FN = 26


### 3.2.9 SCAD USING ORDINIS ####
system.time(ordinis.scad.sc.gaus <- fit.GLM(sim = sim.gaus, method = 9, 
                                            family = gaussian(identity), 
                                            plot = T, IC = "ebic"))  
# 0.33s 
cor(as.vector(ordinis.scad.sc.gaus[-1]), as.vector(sim.gaus$beta_true)) # 0.9205098
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(ordinis.scad.sc.gaus[-1]) != 0)) # TP = 34
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(ordinis.scad.sc.gaus[-1]) != 0)) # FP = 22
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(ordinis.scad.sc.gaus[-1]) == 0)) # TN = 428
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(ordinis.scad.sc.gaus[-1]) == 0)) # FN = 16


### 3.2.10 BENCHMARKS TABLE ####

X <- sim.gaus$X
y <- sim.gaus$y
W <- as.matrix(rep(1,ncol(X)))
IC <- "ebic"

benchmarks_3_2 = rbind(
  data.frame(method="L0glm", do.call(cbind, get.metrics(beta = L0glm.chol$coefficients,  # some metrics
                                                                     X = X, y = y, 
                                                                     weights = W,
                                                                     family = gaussian(identity), 
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

benchmarks_3_2 = benchmarks_3_2[order(benchmarks_3_2[,"deviance_true"]),] # sort by deviance relative to ground truth
benchmarks_3_2
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


## 3.3 MULTICHANNEL POISSON ####
n <- p <- 500
sim.pois.mc <- simul(method = 4, n = n, p = p, seed = 2)

### 3.3.1 RELAXED GROUP L0GLM VIA CHOLESKY ####
system.time(tune.chol.mc <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, IC = "ebic", 
                                       plot = FALSE, coefs.true = sim.pois.mc$beta_true, 
                                       relaxed = FALSE, family = poisson(identity), 
                                       nonnegative = TRUE, 
                                       parallelize = TRUE, grouped = TRUE,
                                       solver = "eigen", method = 7L)) 
# XX
# TODO fix bug 
# Error in h(simpleError(msg, call)) : 
# error in evaluating the argument 'x' in selecting a method for function 'mean': non-conformable matrix dimensions in .Arith.Csparse(e1, e2, .Generic, class. = "dgCMatrix")

system.time(L0glm.chol.pois <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                         family = poisson(identity), 
                                         lambda = "ebic",
                                         nonnegative = TRUE,
                                         relaxed = TRUE,
                                         parallelize = TRUE,
                                         grouped = TRUE,
                                         solver = "eigen",
                                         method = 7L)) 
# 13.0s

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


### 3.3.2 BSGS USING ABESS ####
system.time(bsgs.abess.pois <- fit.GLM.mc(sim = sim.pois.mc, method = 1, 
                                             family = poisson(identity), 
                                             plot = T, IC = "ebic")) 
# 4.9s

cor(as.vector(as.vector(sim.pois.mc$beta_true)), as.vector(bsgs.abess.pois)) # 0.9060618
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(bsgs.abess.pois)!=0) # TP = 9
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(bsgs.abess.pois)!=0) # FP = 1
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(bsgs.abess.pois)==0) # TN = 449
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(bsgs.abess.pois)==0) # FN = 41


### 3.3.3 MULTIVARIATE GROUP LASSO USING GLMNET ####
system.time(glmnet.mc.pois <- fit.GLM.mc(sim = sim.pois.mc, method = 2, 
                                    family = poisson(identity), plot = T, 
                                    IC = "ebic"))
# 22s

image(x=1:nrow(sim.pois.mc$y), y=1:ncol(sim.pois.mc$y), z=as.matrix(glmnet.mc.pois)^0.1, 
      col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", 
      main="glmnet nonnegative multivariate group LASSO\n(red = true peaks, cyan = estimated peaks)")
abline(v=(1:nrow(sim.pois.mc$X))[as.vector(rowMax(sim.pois.mc$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim.pois.mc$X))[as.vector(rowMin(glmnet.mc.pois)!=0)], col="cyan")


cor(as.vector(as.vector(sim.pois.mc$beta_true)), as.vector(glmnet.mc.pois)) # 0.987812
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(glmnet.mc.pois)!=0) # 43
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(glmnet.mc.pois)!=0) # 17
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(glmnet.mc.pois)==0) # 433
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(glmnet.mc.pois)==0) # 7


### 3.3.4 BENCHMARKS TABLE ####

W <- (1/(sim.pois.mc$y+0.1))

benchmarks_3_3 = rbind(
  data.frame(method="Grouped_L0glm_poisson", do.call(cbind, get.metrics(beta = as.matrix(L0glm.chol.pois$coefficients),  # TODO: make sure that metrics work with either dense or sparse coefs
                                                                        X = sim.pois.mc$X, y = sim.pois.mc$y,
                                                                        weights = W, 
                                                                        family = gaussian(identity), 
                                                                        beta_true = sim.pois.mc$beta_true) ) ),
  data.frame(method="BSGS_abess_poisson", do.call(cbind, get.metrics(beta = bsgs.abess.poisson, 
                                                                     X = sim.pois.mc$X, y = sim.pois.mc$y, 
                                                                     weights = W,
                                                                     family = gaussian(identity), 
                                                                     beta_true = sim.pois.mc$beta_true) ) ),
  data.frame(method="multivariate_group_LASSO_poisson", do.call(cbind, get.metrics(beta = glmnet.mc.pois,  
                                                                     X = sim.pois.mc$X, y = sim.pois.mc$y, 
                                                                     weights = W,
                                                                     family = gaussian(identity), 
                                                                     beta_true = sim.pois.mc$beta_true) ) ) 
  )
benchmarks_3_3 = benchmarks_3_3[order(benchmarks_3_3[,"deviance_true"]),] # sort by deviance relative to ground truth
benchmarks_3_3
#                             method   p   n   m   deviance deviance_true    loglik     MD_coefs
# 1            Grouped_L0glm_poisson 501 500 160   67399.34      6737.463 -102070.4 -0.005062176
# 3 multivariate_group_LASSO_poisson 501 500 160  236566.92    187659.275 -154353.3 -0.816490366
# 2               BSGS_abess_poisson 501 500 160 2415489.28   2414515.360 -242537.3 -3.180608680
#   MAD_coefs sMAPE_coefs bias_coefs  R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs
# 1 0.4429137  0.02643714  1.0012997 0.9991303    1.156908     7.993593      0.3821764
# 3 2.6313856  0.10274055  1.0381172 0.9758372   -7.819712    32.482804      1.1828171
# 2 6.2203325  0.09183572  0.9146294 0.8211056  106.532307   107.910290      1.0483753
#   bias_selcoefs R2_selcoefs support_size_true support_size_estimated
# 1     1.0012997   0.9993185                50                     50
# 3     1.0381172   0.9771360                50                     60
# 2     0.9146294   0.9430487                50                     10
#   estimated_over_true_support_size TP FP  TN FN precision recall specificity accuracy
# 1                              1.0 49  1 449  1 0.9800000   0.98   0.9977778    0.996
# 3                              1.2 43 17 433  7 0.7166667   0.86   0.9622222    0.952
# 2                              0.2  9  1 449 41 0.9000000   0.18   0.9977778    0.916
#    F1_score      gic     ebic      bic      aic     mbic
# 1 0.9800000 247627.4 251747.7 227936.6 211798.9 264927.2
# 3 0.7818182 375577.4 381913.3 345298.2 320482.6 402179.9
# 2 0.3000000 497965.0 499186.3 492128.2 487344.6 503093.0



## 3.4 MULTICHANNEL GAUSSIAN ####
n <- p <- 500
sim.gaus.mc <- simul(method = 3, n = n, p = p, seed = 1)

### 3.4.1 L0GLM VIA CHOLESKY ####

system.time(tune.chol.mc <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, IC = "ebic", 
                                       plot = FALSE, coefs.true = sim.gaus.mc$beta_true, 
                                       relaxed = FALSE, family = gaussian(identity), 
                                       nonnegative = TRUE, 
                                       parallelize = TRUE, grouped = TRUE,
                                       solver = "eigen", method = 7L)) 
# XXs
# TODO fix bug
# Error in lm.fit(x, y, offset = offset, singular.ok = singular.ok, ...) : 
# 0 (non-NA) cases
# In addition: There were 50 or more warnings (use warnings() to see the first 50)

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
       
R2 = round(cor(as.vector(rowMeans(sim.gaus.mc$beta_true)), as.vector(rowMeans(L0glm.chol.gaus$coefficients))), digits = 4)^2 # Calculate R2
title = paste("True vs estimated coefficients\nrelaxed nonnegative L0 grouped GLM\nR2=", R2)
plot(as.vector(rowMeans(L0glm.chol.gaus$coefficients)), 
     as.vector(rowMeans(sim.gaus.mc$beta_true)), pch=16, 
     xlab="estimated coefficients relaxed nonnegative L0 group GLM",
     ylab="true coefficients", col="steelblue",
     main=title)

image(x=1:nrow(sim.gaus.mc$y), y=1:ncol(sim.gaus.mc$y), z=as.matrix(L0glm.chol.gaus$coefficients)^0.1, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="relaxed nonnegative group L0glm\n(red = true peaks, cyan = estimated peaks)")
abline(v=(1:nrow(sim.gaus.mc$X))[as.vector(rowMax(sim.gaus.mc$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim.gaus.mc$X))[as.vector(rowMin(L0glm.chol.gaus$coefficients)!=0)], col="cyan")

cor(as.vector(as.vector(L0glm.chol.gaus$coefficients)), as.vector(sim.gaus.mc$beta_true)) # 0.999951
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.chol.gaus$coefficients)!=0) # TP = 50
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.chol.gaus$coefficients)!=0) # FP = 0
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.chol.gaus$coefficients)==0) # TN = 450
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.chol.gaus$coefficients)==0) # FN = 0


### 3.4.2 BSGS USING ABESS ####
system.time(bsgs.abess.gaus <- fit.GLM.mc(sim = sim.gaus.mc, method = 1, 
                                          family = gaussian(identity), plot = T, 
                                          IC = "ebic")) 
# 3.67s

cor(as.vector(as.vector(sim.gaus.mc$beta_true)), as.vector(bsgs.abess.gaus)) # 0.6662117
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(bsgs.abess.gaus)!=0) # TP = 17
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(bsgs.abess.gaus)!=0) # FP = 7
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(bsgs.abess.gaus)==0) # TN = 443
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(bsgs.abess.gaus)==0) # FN = 33


# glmnet nonnegative mgaussian multivariate group LASSO fit
system.time(glmnet.mc.gaus <- fit.GLM.mc(sim = sim.gaus.mc, method = 2, 
                                    family = gaussian(identity), plot = T, 
                                    IC = "ebic"))
# 36.69s

image(x=1:nrow(sim.gaus.mc$y), y=1:ncol(sim.gaus.mc$y), z=as.matrix(glmnet.mc.gaus)^0.1, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", 
      main="nonnegative multivariate group LASSO\n(red = true peaks, cyan = estimated peaks)")
abline(v=(1:nrow(sim.gaus.mc$X))[as.vector(rowMax(sim.gaus.mc$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim.gaus.mc$X))[as.vector(rowMin(glmnet.mc.gaus)!=0)], col="cyan")


cor(as.vector(as.vector(sim.gaus.mc$beta_true)), as.vector(glmnet.mc.gaus)) # 0.9817636
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(glmnet.mc.gaus)!=0) # 49
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(glmnet.mc.gaus)!=0) # 31
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(glmnet.mc.gaus)==0) # 419
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(glmnet.mc.gaus)==0) # 1


### 3.4.3 BENCHMARKS TABLE ####

benchmarks_3_4 = rbind(
  data.frame(method="Grouped_L0glm_gaussian", do.call(cbind, get.metrics(beta = as.matrix(L0glm.chol.gaus$coefficients),  # TODO: allow coefs to be either dense or sparse matrix in get.metrics
                                                                         X = sim.gaus.mc$X, y = sim.gaus.mc$y, 
                                                                         family = gaussian(identity),
                                                                         beta_true = sim.gaus.mc$beta_true)) ),
  data.frame(method="BSGS_abess_gaussian", do.call(cbind, get.metrics(beta = bsgs.abess.gaus,  # some metrics
                                                                      X = sim.gaus.mc$X, y = sim.gaus.mc$y, 
                                                                      family = gaussian(identity),
                                                                      beta_true = sim.gaus.mc$beta_true)) ),
  data.frame(method="multivariate_group_LASSO_gaussian", do.call(cbind, get.metrics(beta = glmnet.mc.gaus,  
                                                                                   X = sim.gaus.mc$X, y = sim.gaus.mc$y, 
                                                                                   family = gaussian(identity), 
                                                                                   beta_true = sim.gaus.mc$beta_true)) ) )

benchmarks_3_4 = benchmarks_3_4[order(benchmarks_3_4[,"deviance_true"]),] # sort by deviance relative to ground truth
benchmarks_3_4
#                              method   p   n   m  deviance deviance_true    loglik    MD_coefs
# 1            Grouped_L0glm_gaussian 501 500 165   1125238        207933 -231325.7  0.09550816
# 3 multivariate_group_LASSO_gaussian 501 500 165  10219775       9310928 -326523.8 -0.19791744
# 2               BSGS_abess_gaussian 501 500 165 312127788     311389742 -460225.2 -0.73411430
#   MAD_coefs sMAPE_coefs bias_coefs  R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs
# 1  0.254989  0.09548499  1.0000758 0.9999020    1.223229      2.92006       1.067991
# 3  2.441680  0.20675867  1.0467529 0.9640264   -1.161945     17.73998       1.495933
# 2  7.242264  0.12506039  0.6556086 0.4456978   69.614337     70.57427       1.393076
#   bias_selcoefs R2_selcoefs support_size_true support_size_estimated
# 1     1.0000758   0.9999054                50                     50
# 3     1.0467529   0.9644986                50                     80
# 2     0.6556086   0.7048314                50                     24
#   estimated_over_true_support_size TP FP  TN FN precision recall specificity accuracy
# 1                             1.00 50  0 450  0 1.0000000   1.00   1.0000000    1.000
# 3                             1.60 49 31 419  1 0.6125000   0.98   0.9311111    0.936
# 2                             0.48 17  7 443 33 0.7083333   0.34   0.9844444    0.920
#    F1_score      gic     ebic      bic      aic     mbic
# 1 1.0000000 541924.1 549435.2 506029.3 476611.3 573460.5
# 3 0.7538462 779872.7 791889.3 722446.1 675381.6 830326.3
# 2 0.4594595 963959.5 968082.0 944258.5 928112.4 981268.4