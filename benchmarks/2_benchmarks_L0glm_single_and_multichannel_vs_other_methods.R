# BENCHMARKS L0GLM: SINGLE CHANNEL VS MULTICHANNEL DECONVOLUTION OF BLURRED SPIKE TRAIN ####

library(remotes)
# remotes::install_github("jaredhuling/ordinis")
# remotes::install_github("tomwenseleers/L0glm/L0glm")
devtools::load_all(".")
# library(L0glm)
library(export)
library(Matrix)
library(ordinis)
library(L0Learn) # use at least version 2
library(abess)
library(glmnet)


# 1. SINGLE CHANNEL, SPARSE  ####
# Simulate SC n=p=500
n <- p <- 500 # 10000
sim.gaus <- simul(method = 1, n = n, p = p, sparse = T)
sim.pois <- simul(method = 2, n = n, p = p, sparse = T)

## 1.1. L0GLM FIT SINGLE CHANNEL WITH POISSON NOISE ####

## Tune lambda relaxed = F with gic
system.time(tune.osqp.sc.pois <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "gic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = FALSE, family = poisson(identity), nonnegative = TRUE, solver = "osqp"))
# bestlambda = 10, time = 2.7s (for n=p=500), bestlambda = 10, time = 71s (for n=p=10000)
system.time(tune.glmnet.sc.pois <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "gic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = FALSE, family = poisson(identity), nonnegative = TRUE, solver = "glmnet"))
# bestlambda = 10, time = 8.1s, bestlambda = 10, time = 31s (for n=p=10000)
system.time(tune.chol.sc.pois <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "gic", plot = TRUE, coefs.true = sim.pois$beta_true, relaxed = FALSE, family = poisson(identity), nonnegative = TRUE, solver = "eigen", method = 7L))
# bestlambda = 10, time = 1.1s, bestlambda = 10, time = 5s (for n=p=10000)
system.time(tune.QR.sc.pois <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "gic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = FALSE, family = poisson(identity), nonnegative = TRUE, solver = "eigen", method = 1L))
# bestlambda = 10, time = 1.6s, bestlambda = 10, time = 486s (for n=p=10000)
system.time(tune.CG.sc.pois <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "gic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = FALSE, family = poisson(identity), nonnegative = TRUE, solver = "eigen", method = 5L))
# bestlambda = 10, time = 1.5s, bestlambda = 10, time = 16s (for n=p=10000)
system.time(tune.LSCG.sc.pois <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "gic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = FALSE, family = poisson(identity), nonnegative = TRUE, solver = "eigen", method = 3L))
# bestlambda = 10, time = 1.9s, bestlambda = 10, time = 48s (for n=p=10000)

#### 1.1.1 L0GLM FIT WITH OSQP SOLVER ####
system.time(L0glm.osqp.pois <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                    family = poisson(identity), 
                                    lambda = "gic", # tune.osqp.sc$best.lambda, # choose preset or tuned
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "osqp")) 
# 0.12s with n=p=500, 5.5s with n=p=10000
cor(as.vector(as.vector(L0glm.osqp.pois$coefficients)), as.vector(sim.pois$beta_true)) # 0.9778539 with n=p=500, 0.9422652 with n=p=10000
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.osqp.pois$coefficients) != 0)) # TP = 45 with n=p=500, 839 with n=p=10000
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.osqp.pois$coefficients) != 0)) # FP = 7 with n=p=500, 173 with n=p=10000
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.osqp.pois$coefficients) == 0)) # TN = 443 with n=p=500, 8827 with n=p=10000
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.osqp.pois$coefficients) == 0)) # FN = 5 with n=p=500, 161 with n=p=10000

#### 1.1.2 L0GLM FIT WITH glmnet::bigGlm SOLVER ####
system.time(L0glm.glmnet.pois <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                      family = poisson(identity), 
                                      lambda = "gic", # tune.glmnet.sc$best.lambda,
                                      nonnegative = TRUE,
                                      relaxed = FALSE,
                                      solver = "glmnet")) 
# 0.15s with n=p=500, 2.13s with n=p=10000
cor(as.vector(as.vector(L0glm.glmnet.pois$coefficients)), as.vector(sim.pois$beta_true)) # 0.9775844 with n=p=500, 0.938967 with n=p=10000
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.glmnet.pois$coefficients) != 0)) # TP = 45 with n=p=500, 839 with n=p=10000
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.glmnet.pois$coefficients) != 0)) # FP = 7 with n=p=500, 176 with n=p=10000
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.glmnet.pois$coefficients) == 0)) # TN = 443 with n=p=500, 8824 with n=p=10000
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.glmnet.pois$coefficients) == 0)) # FN = 5 with n=p=500, 161 with n=p=10000

#### 1.1.3 L0GLM FIT WITH QR SOLVER ####
system.time(L0glm.QR.pois <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                  family = poisson(identity), 
                                  lambda = "gic", # tune.QR.sc$best.lambda,
                                  nonnegative = TRUE,
                                  relaxed = FALSE,
                                  solver = "eigen",
                                  method = 1L)) 
# 0.05s with n=p=500, 20.0s with n=p=10000
# NOTE: TODO: check why QR solver is now slower than it used to be for n=p=10000
cor(as.vector(as.vector(L0glm.QR.pois$coefficients)), as.vector(sim.pois$beta_true)) # 0.9780483 with n=p=500, 0.9595521 with n=p=10000
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.QR.pois$coefficients) != 0)) # TP = 45 with n=p=500, 841 with n=p=10000
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.QR.pois$coefficients) != 0)) # FP = 8 with n=p=500, 173 with n=p=10000
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.QR.pois$coefficients) == 0)) # TN = 442 with n=p=500, 8827 with n=p=10000
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.QR.pois$coefficients) == 0)) # FN = 5 with n=p=500, 161 with n=p=10000

#### 1.1.4 L0GLM FIT WITH CHOLESKY SOLVER ####
system.time(L0glm.chol.pois <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                    family = poisson(identity), 
                                    lambda = "gic", # tune.chol.sc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "eigen",
                                    method = 7L)) 
# 0.04s with n=p=500, 0.3s with n=p=10000
cor(as.vector(as.vector(L0glm.chol.pois$coefficients)), as.vector(sim.pois$beta_true)) # 0.9780483 with n=p=500, 0.9595521 with n=p=10000
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.chol.pois$coefficients) != 0)) # TP = 45 with n=p=500, 841 with n=p=10000
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.chol.pois$coefficients) != 0)) # FP = 8 with n=p=500, 173 with n=p=10000
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.chol.pois$coefficients) == 0)) # TN = 442 with n=p=500, 8827 with n=p=10000
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.chol.pois$coefficients) == 0)) # FN = 5 with n=p=500, 159 with n=p=10000

#### 1.1.5 L0GLM FIT WITH LSCG SOLVER ####
system.time(L0glm.LSCG.pois <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                    family = poisson(identity), 
                                    lambda = "gic", # tune.LSCG.sc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "eigen",
                                    method = 3L)) 
# 0.03s with n=p=500, 0.94s with n=p=10000
cor(as.vector(as.vector(L0glm.LSCG.pois$coefficients)), as.vector(sim.pois$beta_true)) # 0.97805 with n=p=500, 0.9595335 with n=p=10000
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.LSCG.pois$coefficients) != 0)) # TP = 45 with n=p=500, 841 with n=p=10000
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.LSCG.pois$coefficients) != 0)) # FP = 8 with n=p=500, 173 with n=p=10000
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.LSCG.pois$coefficients) == 0)) # TN = 442 with n=p=500, 8827 with n=p=10000
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.LSCG.pois$coefficients) == 0)) # FN = 5 with n=p=500, 159 with n=p=10000

#### 1.1.6 L0GLM FIT WITH CG SOLVER ####
system.time(L0glm.CG.pois <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                  family = poisson(identity), 
                                  lambda = "gic", # tune.CG.sc$best.lambda,
                                  nonnegative = TRUE,
                                  relaxed = FALSE,
                                  solver = "eigen",
                                  method = 5L)) 
# 0.05s with n=p=500, 0.51s with n=p=10000
cor(as.vector(as.vector(L0glm.CG.pois$coefficients)), as.vector(sim.pois$beta_true)) # 0.97805 with n=p=500, 0.9593445 with n=p=10000
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.CG.pois$coefficients) != 0)) # TP = 45 with n=p=500, 841 with n=p=10000
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.CG.pois$coefficients) != 0)) # FP = 8 with n=p=500, 173 with n=p=10000
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.CG.pois$coefficients) == 0)) # TN = 442 with n=p=500, 8827 with n=p=10000
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.CG.pois$coefficients) == 0)) # FN = 5 with n=p=500, 159 with n=p=10000


#### 1.1.7 L0GLM FIT WITH NNLS SOLVER (very slow for larger problem size) ####
system.time(L0glm.nnls.pois <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                  family = poisson(identity), 
                                  lambda = "gic", # tune.CG.sc$best.lambda,
                                  nonnegative = TRUE,
                                  relaxed = FALSE,
                                  solver = "nnls")) 
# 0.56s with n=p=500, XXs with n=p=10000
cor(as.vector(as.vector(L0glm.nnls.pois$coefficients)), as.vector(sim.pois$beta_true)) # 0.9778538 with n=p=500, XX with n=p=10000
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.nnls.pois$coefficients) != 0)) # TP = 45 with n=p=500, XX with n=p=10000
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.nnls.pois$coefficients) != 0)) # FP = 7 with n=p=500, XX with n=p=10000
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.nnls.pois$coefficients) == 0)) # TN = 443 with n=p=500, XX with n=p=10000
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.nnls.pois$coefficients) == 0)) # FN = 5 with n=p=500, XX with n=p=10000


#### 1.1.8 PLAIN NNLS FIT WITHOUT REGULARISATION ####
system.time(L0glm.nnlsnoreg.pois <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                    family = poisson(identity), 
                                    lambda = 0,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "nnls",
                                    maxit = 1)) 
# 0.1s with n=p=500, XXs with n=p=10000
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
#           method   p   n m deviance deviance_true    loglik  MD_coefs MAD_coefs sMAPE_coefs bias_coefs  R2_coefs MD_selcoefs
# 5     L0glm_LSCG 501 500 1 519.4664      82.81184 -1565.926 0.2072258  48.77975  0.05809695   1.022185 0.9582754   110.62212
# 6       L0glm_CG 501 500 1 519.4664      82.81185 -1565.926 0.2072265  48.77975  0.05809695   1.022185 0.9582754   110.62212
# 3       L0glm_qr 501 500 1 519.4656      82.81249 -1565.925 0.2072201  48.78212  0.05809725   1.022185 0.9582707   110.62206
# 4 L0glm_cholesky 501 500 1 519.4656      82.81249 -1565.925 0.2072201  48.78212  0.05809725   1.022185 0.9582707   110.62206
# 7     L0glm_nnls 501 500 1 525.6903      84.66144 -1567.903 0.2072207  49.55770  0.05414609   1.024999 0.9578996   112.74942
# 1     L0glm_osqp 501 500 1 525.6903      84.66152 -1567.903 0.2072200  49.55766  0.05414609   1.024999 0.9578997   112.74941
# 2   L0glm_bigglm 501 500 1 525.6329      84.85241 -1567.876 0.2071412  49.81488  0.05418274   1.024866 0.9573938   112.74865
# 8     nnls_noreg 501 500 1 443.1191      95.59246 -1568.185 0.1159119  72.72359  0.23426504   1.032099 0.9626282    74.87893
#   MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true support_size_estimated estimated_over_true_support_size TP
# 5     351.5191      0.3594052      1.022185   0.9665981                50                     53                             1.06 45
# 6     351.5192      0.3594052      1.022185   0.9665981                50                     53                             1.06 45
# 3     351.5415      0.3594081      1.022185   0.9665933                50                     53                             1.06 45
# 4     351.5415      0.3594081      1.022185   0.9665933                50                     53                             1.06 45
# 7     365.7595      0.3283278      1.024999   0.9662189                50                     52                             1.04 45
# 1     365.7591      0.3283278      1.024999   0.9662190                50                     52                             1.04 45
# 2     368.2323      0.3286802      1.024866   0.9657088                50                     52                             1.04 45
# 8     308.4869      1.1277107      1.032099   0.9730918                50                     95                             1.90 45
#   FP  TN FN precision recall specificity accuracy  F1_score      gic     gic      bic      aic     mbic
# 5  8 442  5 0.8490566    0.9   0.9822222    0.974 0.8737864 3745.137 3803.245 3467.440 3239.851 3989.114
# 6  8 442  5 0.8490566    0.9   0.9822222    0.974 0.8737864 3745.137 3803.245 3467.440 3239.851 3989.114
# 3  8 442  5 0.8490566    0.9   0.9822222    0.974 0.8737864 3745.136 3803.244 3467.439 3239.851 3989.113
# 4  8 442  5 0.8490566    0.9   0.9822222    0.974 0.8737864 3745.136 3803.244 3467.439 3239.851 3989.113
# 7  7 443  5 0.8653846    0.9   0.9844444    0.976 0.8823529 3737.735 3794.767 3465.181 3241.806 3977.194
# 1  7 443  5 0.8653846    0.9   0.9844444    0.976 0.8823529 3737.735 3794.767 3465.181 3241.806 3977.194
# 2  7 443  5 0.8653846    0.9   0.9844444    0.976 0.8823529 3737.680 3794.712 3465.126 3241.752 3977.139
# 8 50 400  5 0.4736842    0.9   0.8888889    0.890 0.6206897 4226.655 4329.958 3732.972 3328.370 4660.392




## 1.2. L0GLM FIT SINGLE CHANNEL WITH GAUSSIAN NOISE ####

## Tune lambda relaxed = F 
system.time(tune.osqp.sc.gaus <- L0glm.tune(X=sim.gaus$X, y=sim.gaus$y, IC = "gic", plot = FALSE, coefs.true = sim.gaus$beta_true, relaxed = FALSE, family = gaussian(identity), nonnegative = TRUE, solver = "osqp"))
# bestlambda = 100, time = 3s
system.time(tune.glmnet.sc.gaus <- L0glm.tune(X=sim.gaus$X, y=sim.gaus$y, IC = "gic", plot = FALSE, coefs.true = sim.gaus$beta_true, relaxed = FALSE, family = gaussian(identity), nonnegative = TRUE, solver = "glmnet"))
# bestlambda = 100, time = 8.7s
system.time(tune.chol.sc.gaus <- L0glm.tune(X=sim.gaus$X, y=sim.gaus$y, IC = "gic", plot = TRUE, coefs.true = sim.gaus$beta_true, relaxed = FALSE, family = gaussian(identity), nonnegative = TRUE, solver = "eigen", method = 7L))
# bestlambda = 316.227766016838, time = 1.1s
system.time(tune.QR.sc.gaus <- L0glm.tune(X=sim.gaus$X, y=sim.gaus$y, IC = "gic", plot = FALSE, coefs.true = sim.gaus$beta_true, relaxed = TRUE, family = gaussian(identity), nonnegative = TRUE, solver = "eigen", method = 1L))
# bestlambda = 316.227766016838, time = 2.3s
system.time(tune.CG.sc.gaus <- L0glm.tune(X=sim.gaus$X, y=sim.gaus$y, IC = "gic", plot = FALSE, coefs.true = sim.gaus$beta_true, relaxed = FALSE, family = gaussian(identity), nonnegative = TRUE, solver = "eigen", method = 5L))
# bestlambda = 316.227766016838, time = 1.17s
system.time(tune.LSCG.sc.gaus <- L0glm.tune(X=sim.gaus$X, y=sim.gaus$y, IC = "gic", plot = TRUE, coefs.true = sim.gaus$beta_true, relaxed = FALSE, family = gaussian(identity), nonnegative = TRUE, solver = "eigen", method = 3L))
# bestlambda = 316.227766016838, time = 3.1s

# TODO: in plot with regularisation path max lambda is going too high & is not correctly calculated - FIX


#### 1.2.1 L0GLM FIT WITH OSQP SOLVER ####
system.time(L0glm.osqp.gaus <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                    family = gaussian(identity), 
                                    lambda = "gic", # tune.osqp.sc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "osqp")) 
# 0.64s for n=p=500, 42s for n=p=10000
cor(as.vector(as.vector(L0glm.osqp.gaus$coefficients)), as.vector(sim.gaus$beta_true)) #  0.9900089
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.osqp.gaus$coefficients) != 0)) # TP = 48
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.osqp.gaus$coefficients) != 0)) # FP = 8
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.osqp.gaus$coefficients) == 0)) # TN = 442
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.osqp.gaus$coefficients) == 0)) # FN = 2

#### 1.2.2 L0GLM FIT WITH glmnet::bigGlm SOLVER ####
system.time(L0glm.glmnet.gaus <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                      family = gaussian(identity), 
                                      lambda = "gic", # tune.glmnet.sc$best.lambda,
                                      nonnegative = TRUE,
                                      relaxed = FALSE,
                                      solver = "glmnet")) 
# 0.56s for n=p=500, 5.6s for n=p=10000 
cor(as.vector(as.vector(L0glm.glmnet.gaus$coefficients)), as.vector(sim.gaus$beta_true)) # 0.9626836
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
# 0.05s for n=p=500, 143s for n=p=10000
cor(as.vector(as.vector(L0glm.QR.gaus$coefficients)), as.vector(sim.gaus$beta_true)) # 0.992033
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.QR.gaus$coefficients) != 0)) # TP = 47
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.QR.gaus$coefficients) != 0)) # FP = 9
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.QR.gaus$coefficients) == 0)) # TN = 441
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.QR.gaus$coefficients) == 0)) # FN = 3

#### 1.2.4 L0GLM FIT WITH Cholesky SOLVER ####
system.time(L0glm.chol.gaus <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                    family = gaussian(identity), 
                                    lambda = "gic", # tune.chol.sc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "eigen",
                                    method = 7L)) 
# 0.01s for n=p=500, 0.25s for n=p=10000
cor(as.vector(as.vector(L0glm.chol.gaus$coefficients)), as.vector(sim.gaus$beta_true)) # 0.992033
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.chol.gaus$coefficients) != 0)) # TP = 47
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.chol.gaus$coefficients) != 0)) # FP = 9
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.chol.gaus$coefficients) == 0)) # TN = 441
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.chol.gaus$coefficients) == 0)) # FN = 3

#### 1.2.5 L0GLM FIT WITH LSCG SOLVER ####
system.time(L0glm.LSCG.gaus <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                    family = gaussian(identity), 
                                    lambda = "gic", # tune.CG.sc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "eigen",
                                    method = 4L)) 
# 0.05s for n=p=500, 8.5s for n=p=10000
cor(as.vector(as.vector(L0glm.LSCG.gaus$coefficients)), as.vector(sim.gaus$beta_true)) # 0.9920338
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.LSCG.gaus$coefficients) != 0)) # TP = 47
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.LSCG.gaus$coefficients) != 0)) # FP = 9
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.LSCG.gaus$coefficients) == 0)) # TN = 441
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.LSCG.gaus$coefficients) == 0)) # FN = 3

#### 1.2.6 L0GLM FIT WITH CG SOLVER ####
system.time(L0glm.CG.gaus <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                  family = gaussian(identity), 
                                  lambda = "gic", # tune.CG.sc$best.lambda,
                                  nonnegative = TRUE,
                                  relaxed = TRUE,
                                  solver = "eigen",
                                  method = 5L)) 
# 0.10s for n=p=500, 3.7s for n=p=10000
cor(as.vector(as.vector(L0glm.CG.gaus$coefficients)), as.vector(sim.gaus$beta_true)) # 0.9920172
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.CG.gaus$coefficients) != 0)) # TP = 47
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.CG.gaus$coefficients) != 0)) # FP = 9
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.CG.gaus$coefficients) == 0)) # TN = 441
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.CG.gaus$coefficients) == 0)) # FN = 3


#### 1.2.7 L0GLM FIT WITH NNLS SOLVER (slow for large problem sizes) ####
system.time(L0glm.nnls.gaus <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                    family = gaussian(identity), 
                                    lambda = "gic", # tune.CG.sc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "nnls")) 
# 0.33s with n=p=500, XXs with n=p=10000
cor(as.vector(as.vector(L0glm.nnls.gaus$coefficients)), as.vector(sim.pois$beta_true)) # 0.9991312 with n=p=500, XX with n=p=10000
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.nnls.gaus$coefficients) != 0)) # TP = 49 with n=p=500, XX with n=p=10000
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.nnls.gaus$coefficients) != 0)) # FP = 4 with n=p=500, XX with n=p=10000
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.nnls.gaus$coefficients) == 0)) # TN = 446 with n=p=500, XX with n=p=10000
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.nnls.gaus$coefficients) == 0)) # FN = 1 with n=p=500, XX with n=p=10000


#### 1.2.8 PLAIN NNLS FIT WITHOUT REGULARISATION ####
system.time(L0glm.nnlsnoreg.gaus <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                         family = gaussian(identity), 
                                         lambda = 0,
                                         nonnegative = TRUE,
                                         relaxed = FALSE,
                                         solver = "nnls",
                                         maxit = 1)) 
# 0.09s with n=p=500, XXs with n=p=10000
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
#           method   p   n m  deviance deviance_true    loglik    MD_coefs MAD_coefs sMAPE_coefs bias_coefs  R2_coefs MD_selcoefs
# 7     L0glm_nnls 501 500 1  9813.156      1561.416 -1006.687 0.003784367  5.807666  0.02071088   1.000932 0.9983316   16.739595
# 1     L0glm_osqp 501 500 1  9654.912      2081.956 -1005.623 0.011256735 25.710028  0.04156464   1.024850 0.9808769   53.332229
# 3       L0glm_qr 501 500 1 11421.772      2947.926 -1047.636 0.025208535 22.392914  0.04937338   1.024253 0.9847315   57.320596
# 4 L0glm_cholesky 501 500 1 11420.458      2948.547 -1047.608 0.025223430 22.393325  0.04937354   1.024253 0.9847315   57.320729
# 5     L0glm_LSCG 501 500 1 11420.322      2948.759 -1047.605 0.025218801 22.392502  0.04937357   1.024250 0.9847329   57.320688
# 6       L0glm_CG 501 500 1 11271.417      3262.546 -1044.324 0.033115675 22.530908  0.04941450   1.024352 0.9847010   57.391195
# 8     nnls_noreg 501 500 1  8040.643      3439.262 -1029.883 0.104028079 17.717465  0.31621381   1.011397 0.9947146    7.439051
# 2   L0glm_bigglm 501 500 1 16993.177      9003.334 -1161.959 0.066070039 83.192194  0.12280366   1.069266 0.9294499   40.825361
# MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true support_size_estimated estimated_over_true_support_size TP
# 7     38.08541      0.1576498      1.000932   0.9990804                50                     53                             1.06 49
# 1    176.32210      0.2996843      1.024850   0.9857552                50                     56                             1.12 48
# 3    142.84121      0.3336909      1.024253   0.9896734                50                     56                             1.12 47
# 4    142.84488      0.3336924      1.024253   0.9896734                50                     56                             1.12 47
# 5    142.83754      0.3336926      1.024250   0.9896749                50                     56                             1.12 47
# 6    144.07330      0.3340581      1.024352   0.9896429                50                     56                             1.12 47
# 8     63.28116      1.2389437      1.011397   0.9954606                50                    126                             2.52 49
# 2    545.50044      0.7803075      1.069266   0.9336861                50                     71                             1.42 47
# FP  TN FN precision recall specificity accuracy  F1_score      gic     gic      bic      aic     mbic
# 7  4 446  1 0.9245283   0.98   0.9911111    0.990 0.9514563 2626.659 2684.768 2348.963 2121.374 2870.637
# 1  8 442  2 0.8571429   0.96   0.9822222    0.980 0.9056604 2658.602 2719.939 2365.478 2125.246 2916.134
# 3  9 441  3 0.8392857   0.94   0.9800000    0.976 0.8867925 2742.629 2803.966 2449.505 2209.273 3000.161
# 4  9 441  3 0.8392857   0.94   0.9800000    0.976 0.8867925 2742.572 2803.908 2449.448 2209.215 3000.103
# 5  9 441  3 0.8392857   0.94   0.9800000    0.976 0.8867925 2742.566 2803.902 2449.442 2209.209 3000.098
# 6  9 441  3 0.8392857   0.94   0.9800000    0.976 0.8867925 2736.004 2797.340 2442.880 2202.647 2993.535
# 8 77 373  1 0.3888889   0.98   0.8288889    0.844 0.5568182 3502.123 3638.785 2849.022 2313.767 4075.921
# 2 24 426  3 0.6619718   0.94   0.9466667    0.946 0.7768595 3141.632 3219.109 2771.370 2467.918 3466.935




# 2. MULTICHANNEL ####

# TODO: the multichannel algo could still be optimised by pre-clustering the outcome
# variables & summing (for poisson or binomial) or averaging (for other families)
# the outcome variables, then doing the group L0 fit on those & then refitting the
# selected variables on all the original outcome variables
# so allow for an option precluster = T or F or indices of variables to be clustered together

# Simulate MC n=p=500
n <- p <- 500
sim.gaus.mc <- simul(method = 3, n = n, p = p, seed = 1)
sim.pois.mc <- simul(method = 4, n = n, p = p, seed = 2)



## 2.1. MULTICHANNEL WITH POISSON NOISE ####

### 2.1.1 GROUPED L0GLM FIT WITH relaxed = FALSE, grouped = FALSE (not great) ####

# TODO: tuning below ends in an error - I think due to bugs in metrics.R - FIX
# TODO: fix metrics.R so that it works when beta / beta_true are either a vector or a dense or sparse matrix

## Tune lambda relaxed = FALSE, grouped = FALSE
system.time(tune.osqp.mc <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, IC = "gic", 
                                       plot = FALSE, coefs.true = sim.pois.mc$beta_true, 
                                       relaxed = FALSE, family = poisson(identity), 
                                       nonnegative = TRUE, 
                                       parallelize = TRUE, grouped = FALSE, solver = "osqp"))
# bestlambda = 3.16227766016838, time = 127s
# abline(h=tune.osqp.mc$best.lambda, col = 'cyan')
# TODO fix bug 8 nodes produced errors; first error: object 'intercept' not found

system.time(tune.chol.mc <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, IC = "gic", 
                                       plot = FALSE, coefs.true = sim.pois.mc$beta_true, 
                                       relaxed = FALSE, family = poisson(identity), 
                                       nonnegative = TRUE, 
                                       parallelize = TRUE, grouped = FALSE,
                                       solver = "eigen", method = 7L))
# bestlambda = 3.16227766016838, time = 113s
# abline(h=tune.chol.mc$best.lambda, col = 'cyan')


#### OSQP ####
system.time(L0glm.osqp <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                    family = poisson(identity), 
                                    lambda = "gic", # tune.osqp.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallelize = TRUE,
                                    grouped = FALSE,
                                    solver = "osqp")) 
# 7.1s

cor(as.vector(as.vector(L0glm.osqp$coefficients)), as.vector(sim.pois.mc$beta_true)) # 0.9935513
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)!=0) # TP = 50
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)!=0) # FP = 158  # so high because grouped=FALSE
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)==0) # TN = 292
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)==0) # FN = 0

#### Cholesky ####
system.time(L0glm.chol <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                    family = poisson(identity), 
                                    lambda = "gic", # tune.chol.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallelize = TRUE,
                                    grouped = FALSE,
                                    solver = "eigen",
                                    method = 7L)) 
# 6.5s

R2 = round(cor(as.vector(rowMeans(sim.pois.mc$beta_true)), as.vector(rowMeans(L0glm.chol$coefficients))), digits = 4)^2 # Calculate R2
title = paste("true vs estimated coefficients\nnonnegative independent L0 GLM\nR2=", R2)
plot(as.vector(rowMeans(L0glm.chol$coefficients)), as.vector(rowMeans(sim.pois.mc$beta_true)), pch=16, 
     xlab="estimated coefficients",
     ylab="true coefficients", col="steelblue",
     main=title)

cor(as.vector(as.vector(L0glm.chol$coefficients)), as.vector(sim.pois.mc$beta_true)) # 0.9947139
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.chol$coefficients)!=0) # TP = 50
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.chol$coefficients)!=0) # FP = 157  # so high because grouped=FALSE
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.chol$coefficients)==0) # TN = 293
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.chol$coefficients)==0) # FN = 0




### 2.1.2 GROUPED L0GLM FIT WITH relaxed = FALSE, grouped = TRUE (better) ####

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


#### OSQP ####
system.time(L0glm.osqp <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                    family = poisson(identity), 
                                    lambda = "gic", # tune.osqp.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallelize = TRUE,
                                    grouped = TRUE,
                                    solver = "osqp"))
# 14s

cor(as.vector(as.vector(L0glm.osqp$coefficients)), as.vector(sim.pois.mc$beta_true)) # 0.9992514
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)!=0) # TP = 49
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)!=0) # FP = 2
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)==0) # TN = 448
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)==0) # FN = 1

#### Cholesky ####
system.time(L0glm.chol <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                    family = poisson(identity), 
                                    lambda = "gic", # tune.chol.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallelize = TRUE,
                                    grouped = TRUE,
                                    solver = "eigen",
                                    method = 7L)) 
# 12s
# abline(v=which(rowMin(L0glm.chol$coefficients)!=0), col = 'cyan')

cor(as.vector(as.vector(L0glm.chol$coefficients)), as.vector(sim.pois.mc$beta_true)) # 0.9992514
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.chol$coefficients)!=0) # TP = 49
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.chol$coefficients)!=0) # FP = 1
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.chol$coefficients)==0) # TN = 449
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.chol$coefficients)==0) # FN = 1



## 2.1.3 GROUPED L0GLM FIT WITH relaxed = TRUE, grouped = TRUE (best) ####

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
#### OSQP ####
system.time(L0glm.osqp <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                    family = poisson(identity), 
                                    lambda = "gic", # tune.osqp.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = TRUE,
                                    parallelize = TRUE,
                                    grouped = TRUE,
                                    solver = "osqp")) 
# 17s
# abline(v=which(rowMin(L0glm.osqp$coefficients)!=0), col = 'cyan')

cor(as.vector(as.vector(L0glm.osqp$coefficients)), as.vector(sim.pois.mc$beta_true)) # 0.9995656
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)!=0) # TP = 49
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)!=0) # FP = 2
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)==0) # TN = 448
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)==0) # FN = 1

#### Cholesky ####
system.time(L0glm.chol <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                    family = poisson(identity), 
                                    lambda = "gic", # tune.chol.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = TRUE,
                                    parallelize = TRUE,
                                    grouped = TRUE,
                                    solver = "eigen",
                                    method = 7L)) 
# 13s
# abline(v=which(rowMin(L0glm.chol$coefficients)!=0), col = 'cyan')

R2 = round(cor(as.vector(rowMeans(sim.pois.mc$beta_true)), as.vector(rowMeans(L0glm.chol$coefficients))), digits = 4)^2 # Calculate R2
title = paste("true vs estimated coefficients\nnonnegative L0 grouped relaxed GLM\nR2=", R2)
plot(as.vector(rowMeans(L0glm.chol$coefficients)), as.vector(rowMeans(sim.pois.mc$beta_true)), pch=16, 
     xlab="estimated coefficients",
     ylab="true coefficients", col="steelblue",
     main=title)

cor(as.vector(as.vector(L0glm.chol$coefficients)), as.vector(sim.pois.mc$beta_true)) # 0.9995638
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.chol$coefficients)!=0) # TP = 49
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.chol$coefficients)!=0) # FP = 1
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(L0glm.chol$coefficients)==0) # TN = 449
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(L0glm.chol$coefficients)==0) # FN = 1


## 2.2. MULTICHANNEL WITH GAUSSIAN NOISE ####

### 2.2.1 GROUPED L0GLM FIT WITH relaxed = FALSE, grouped = FALSE (not great) ####

# TODO fix bugs that appear below
# Error in lm.fit(x, y, offset = offset, singular.ok = singular.ok, ...) : 
# 0 (non-NA) cases
# In addition: There were 17 warnings (use warnings() to see them)
# probably something going wrong with passing the dispersion

## Tune lambda relaxed = FALSE, grouped = FALSE
system.time(tune.osqp.mc <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, IC = "gic", 
                                       plot = FALSE, coefs.true = sim.gaus.mc$beta_true, 
                                       relaxed = FALSE, family = gaussian(identity), 
                                       nonnegative = TRUE, 
                                       parallelize = TRUE, grouped = FALSE, solver = "osqp"))
# bestlambda = 0.1, time = 156.45s
# TODO fix warnings - probably something in metrics.R?
# 1: In (function (..., deparse.level = 1)  ... :
# number of rows of result is not a multiple of vector length (arg 6)

system.time(tune.chol.mc <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, IC = "gic", 
                                       plot = FALSE, coefs.true = sim.gaus.mc$beta_true, 
                                       relaxed = FALSE, family = gaussian(identity), 
                                       nonnegative = TRUE, 
                                       parallelize = TRUE, grouped = FALSE,
                                       solver = "eigen", method = 7L))
# bestlambda = 0.1, time = 147.31s
# TODO fix warnings - probably something in metrics.R?

#### OSQP ####
system.time(L0glm.osqp <- L0glm.fit(sim.gaus.mc$X, sim.gaus.mc$y, 
                                    family = gaussian(identity), 
                                    lambda = "gic", # tune.osqp.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallelize = TRUE,
                                    grouped = FALSE,
                                    solver = "osqp")) 
# 9s
# TODO fix warning
# number of rows of result is not a multiple of vector length (arg 6)

# abline(v=which(rowMin(L0glm.osqp$coefficients)!=0), col = 'cyan')

cor(as.vector(as.vector(L0glm.osqp$coefficients)), as.vector(sim.gaus.mc$beta_true)) # 0.9963107
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)!=0) # TP = 50
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)!=0) # FP = 149 # so high because grouped=FALSE
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)==0) # TN = 301
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)==0) # FN = 0

#### Cholesky ####
system.time(L0glm.chol <- L0glm.fit(sim.gaus.mc$X, sim.gaus.mc$y, 
                                    family = gaussian(identity), 
                                    lambda = "gic", # tune.chol.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallelize = TRUE,
                                    grouped = FALSE,
                                    solver = "eigen",
                                    method = 7L)) 
# 7.4s
# TODO fix warning
# number of rows of result is not a multiple of vector length (arg 6)

# abline(v=which(rowMin(L0glm.chol$coefficients)!=0), col = 'cyan')

R2 = round(cor(as.vector(rowMeans(sim.gaus.mc$beta_true)), as.vector(rowMeans(L0glm.chol$coefficients))), digits = 4)^2 # Calculate R2
title = paste("true vs estimated coefficients\nnonnegative L0 independent GLM\nR2=", R2)
plot(as.vector(rowMeans(L0glm.chol$coefficients)), as.vector(rowMeans(sim.gaus.mc$beta_true)), pch=16, 
     xlab="estimated coefficients",
     ylab="true coefficients", col="steelblue",
     main=title)

cor(as.vector(as.vector(L0glm.chol$coefficients)), as.vector(sim.gaus.mc$beta_true)) # 0.9904238
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.chol$coefficients)!=0) # TP = 50
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.chol$coefficients)!=0) # FP = 248 # so high because grouped=FALSE
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.chol$coefficients)==0) # TN = 202
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.chol$coefficients)==0) # FN = 0


### 2.2.2 GROUPED L0GLM FIT WITH relaxed = FALSE, grouped = TRUE (better) ####

## Tune lambda relaxed = FALSE, grouped = TRUE
system.time(tune.osqp.mc <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, IC = "gic", 
                                       plot = TRUE, coefs.true = sim.gaus.mc$beta_true, 
                                       relaxed = FALSE, family = gaussian(identity), 
                                       nonnegative = TRUE, 
                                       parallelize = TRUE, grouped = TRUE, solver = "osqp"))
# bestlambda = 0.1, time = 282.71s
# TODO check why coefficients drop of to zero much before calculated maximum lambda
# check if I calculate maximum lambda correctly

system.time(tune.chol.mc <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, IC = "gic", 
                                       plot = FALSE, coefs.true = sim.gaus.mc$beta_true, 
                                       relaxed = FALSE, family = gaussian(identity), 
                                       nonnegative = TRUE, 
                                       parallelize = TRUE, grouped = TRUE,
                                       solver = "eigen", method = 7L))
# bestlambda = 0.1, time = 213.20s 


#### OSQP ####
system.time(L0glm.osqp <- L0glm.fit(sim.gaus.mc$X, sim.gaus.mc$y, 
                                    family = gaussian(identity), 
                                    lambda = "gic", # tune.osqp.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallelize = TRUE,
                                    grouped = TRUE,
                                    solver = "osqp")) 
# 23s for n=p=500

cor(as.vector(as.vector(L0glm.osqp$coefficients)), as.vector(sim.gaus.mc$beta_true)) # 0.9999483
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)!=0) # TP = 50
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)!=0) # FP = 0
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)==0) # TN = 450
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)==0) # FN = 0

#### Cholesky ####
system.time(L0glm.chol <- L0glm.fit(sim.gaus.mc$X, sim.gaus.mc$y, 
                                    family = gaussian(identity), 
                                    lambda = "gic", # tune.chol.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallelize = TRUE,
                                    grouped = TRUE,
                                    solver = "eigen",
                                    method = 7L)) 
# 12s for n=p=500

R2 = round(cor(as.vector(rowMeans(sim.gaus.mc$beta_true)), as.vector(rowMeans(L0glm.chol$coefficients))), digits = 4)^2 # Calculate R2
title = paste("true vs estimated coefficients\nnonnegative L0 grouped GLM\nR2=", R2)
plot(as.vector(rowMeans(L0glm.chol$coefficients)), as.vector(rowMeans(sim.gaus.mc$beta_true)), pch=16, 
     xlab="estimated coefficients",
     ylab="true coefficients", col="steelblue",
     main=title)

cor(as.vector(as.vector(L0glm.chol$coefficients)), as.vector(sim.gaus.mc$beta_true)) # 0.9999476
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.chol$coefficients)!=0) # TP = 50
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.chol$coefficients)!=0) # FP = 0
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.chol$coefficients)==0) # TN = 450
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.chol$coefficients)==0) # FN = 0


### 2.2.3 GROUPED L0GLM FIT WITH relaxed = TRUE, grouped = TRUE (best) ####

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


#### OSQP ####
system.time(L0glm.osqp <- L0glm.fit(sim.gaus.mc$X, sim.gaus.mc$y, 
                                    family = gaussian(identity), 
                                    lambda = "gic", # tune.osqp.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = TRUE,
                                    parallelize = TRUE,
                                    grouped = TRUE,
                                    solver = "osqp")) 
# 23s for n=p=500

cor(as.vector(as.vector(L0glm.osqp$coefficients)), as.vector(sim.gaus.mc$beta_true)) # 0.9999514
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)!=0) # TP = 50
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)!=0) # FP = 0
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(L0glm.osqp$coefficients)==0) # TN = 450
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(L0glm.osqp$coefficients)==0) # FN = 0


#### Cholesky ####
system.time(L0glm.chol <- L0glm.fit(sim.gaus.mc$X, sim.gaus.mc$y, 
                                    family = gaussian(identity), 
                                    lambda = "gic", # tune.chol.mc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = TRUE,
                                    parallelize = TRUE,
                                    grouped = TRUE,
                                    solver = "eigen",
                                    method = 7L)) 
# 13.71s

R2 = round(cor(as.vector(rowMeans(sim.gaus.mc$beta_true)), as.vector(rowMeans(L0glm.chol$coefficients))), digits = 4)^2 # Calculate R2
title = paste("true vs estimated coefficients\nnonnegative L0 grouped relaxed GLM\nR2=", R2)
plot(as.vector(rowMeans(L0glm.chol$coefficients)), as.vector(rowMeans(sim.gaus.mc$beta_true)), pch=16, 
     xlab="estimated coefficients",
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
system.time(tune.chol.sc <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "gic", 
                                       plot = FALSE, coefs.true = sim.pois$beta_true, 
                                       relaxed = FALSE, family = poisson(identity), 
                                       nonnegative = TRUE, 
                                       solver = "eigen", method = 7L))
# 1s

system.time(L0glm.chol <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                    family = poisson(identity), 
                                    lambda = "gic", # tune.chol.sc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "eigen",
                                    method = 7L)) 
# 0.03s
cor1 <- cor(as.vector(as.vector(L0glm.chol$coefficients)), as.vector(sim.pois$beta_true)) # 0.9672425
cor1
TP1 <- sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.chol$coefficients) != 0)) # TP = 43
TP1
FP1 <- sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.chol$coefficients) != 0)) # FP = 6
FP1
TN1 <- sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0glm.chol$coefficients) == 0)) # TN = 444
TN1
FN1 <- sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0glm.chol$coefficients) == 0)) # FN = 7
FN1
plot(x=L0glm.chol$coefficients, y=sim.pois$beta_true, pch=16, col='blue',
     main="L0 glm identity link Poisson with lambda tuned to optimise gic",
     xlab="estimatedcoefficients", ylab="true coefficients")


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
# 0.30s

cor(as.vector(L0learn.sc.pois), as.vector(sim.pois$beta_true)) # 0.8512041
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0learn.sc.pois) != 0)) # TP = 13
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0learn.sc.pois) != 0)) # FP = 1
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0learn.sc.pois) == 0)) # TN = 449
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0learn.sc.pois) == 0)) # FN = 37


### 3.1.5 L0LEARN AFTER LASSO ####
# Slightly better than regular L0Learn & faster
system.time(L0learn.after.lasso.sc.pois <- fit.GLM(sim = sim.pois, method = 5, family = poisson(identity), plot = T)) 
# 0.7s
cor(as.vector(L0learn.after.lasso.sc.pois), as.vector(sim.pois$beta_true)) # 0.8512041
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0learn.after.lasso.sc.pois) != 0)) # TP = 13
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0learn.after.lasso.sc.pois) != 0)) # FP = 1
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(L0learn.after.lasso.sc.pois) == 0)) # TN = 449
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(L0learn.after.lasso.sc.pois) == 0)) # FN = 37


### 3.1.6 ABESS ####
system.time(abess.sc.pois <- fit.GLM(sim = sim.pois, method = 6, 
                                     family = poisson(identity), 
                                     plot = T, IC = "gic"))  
# 0.8s
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
                                                 IC = "gic"))  
# 0.22s
cor(as.vector(abess.after.lasso.sc.pois), as.vector(sim.pois$beta_true))  # 0.3527641
# TODO: NOTE: result is better if I do not use weights - then cor 0.8849495 - CHECK
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(abess.after.lasso.sc.pois) != 0)) # TP = 14
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(abess.after.lasso.sc.pois) != 0)) # FP = 13
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(abess.after.lasso.sc.pois) == 0)) # TN = 437
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(abess.after.lasso.sc.pois) == 0)) # FN = 36

### 3.1.8 MCP USING ORDINIS ####
system.time(ordinis.mcp.sc.pois <- fit.GLM(sim = sim.pois, method = 8, 
                                           family = poisson(identity), plot = T,
                                           IC = "gic"))  
# 11s
cor(as.vector(ordinis.mcp.sc.pois[-1]), as.vector(sim.pois$beta_true)) # 0.7109204
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(ordinis.mcp.sc.pois[-1]) != 0)) # TP = 47
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(ordinis.mcp.sc.pois[-1]) != 0)) # FP = 85
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(ordinis.mcp.sc.pois[-1]) == 0)) # TN = 365
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(ordinis.mcp.sc.pois[-1]) == 0)) # FN = 3

### 3.1.9 SCAD USING ORDINIS ####
system.time(ordinis.scad.sc.pois <- fit.GLM(sim = sim.pois, method = 9, 
                                            family = poisson(identity), 
                                            plot = T, IC = "gic"))  
# 10.5s
cor(as.vector(ordinis.scad.sc.pois[-1]), as.vector(sim.pois$beta_true)) # 0.7141376
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(ordinis.scad.sc.pois[-1]) != 0)) # TP = 47
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(ordinis.scad.sc.pois[-1]) != 0)) # FP = 85
sum((as.vector(sim.pois$beta_true) == 0) & (as.vector(ordinis.scad.sc.pois[-1]) == 0)) # TN = 365
sum((as.vector(sim.pois$beta_true) != 0) & (as.vector(ordinis.scad.sc.pois[-1]) == 0)) # FN = 3


X <- sim.pois$X
y <- sim.pois$y
W <- 1/(y+0.1)
IC <- "gic"


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
#                           method   p   n m     deviance deviance_true    loglik      MD_coefs   MAD_coefs sMAPE_coefs  bias_coefs
# 1                          L0glm 501 500 1 4.794381e+02  5.883981e+01 -1696.678     0.4515978    54.96617  0.05471149 0.959293459
# 2                nnalasso_glmnet 502 500 1 8.799060e+02  5.003641e+02 -1853.478    -1.2191931    73.01081  0.09703150 1.010322884
# 7                 nnlasso_glmnet 502 500 1 9.412260e+02  5.779549e+02 -1944.320    -6.4230012   168.57499  0.36660963 1.156047013
# 9                 nnSCAD_ordinis 502 500 1 1.547711e+03  1.228634e+03 -2071.657    -0.2072492   322.92210  0.40784001 1.195603989
# 8                  nnMCP_ordinis 502 500 1 1.582690e+03  1.262712e+03 -2077.244    -0.2296764   325.44781  0.40850065 1.194418337
# 3 nnL0Learn_after_glmnet_nnlasso 501 500 1 2.749567e+05  2.735459e+05 -3249.616   -87.2668739   204.21553  0.16016479 0.794622032
# 6                      nnL0Learn 501 500 1 2.749567e+05  2.735459e+05 -3249.616   -87.2668739   204.21553  0.16016479 0.794622032
# 5   nnabess_after_glmnet_nnlasso 501 500 1 3.976425e+09  3.976522e+09 -5657.436 17789.8020892 18214.16146  0.25006342 0.004393986
# 4                        nnabess 501 500 1 3.987785e+09  3.987876e+09 -5657.149 17978.4522010 18430.23369  0.25890768 0.004496484
#    R2_coefs   MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true support_size_estimated
# 1 0.9385892    259.072826     306.4146      0.2725662   0.959293459   0.9719406                50                     49
# 2 0.9564822    134.245569     519.7108      0.5929597   1.010322884   0.9645622                50                     55
# 7 0.8535633     -5.307428     635.0605      1.3928018   1.156047013   0.8560704                50                    129
# 9 0.5315166     29.392924    1195.4566      1.5024837   1.195603989   0.5341256                50                    132
# 8 0.5271784     29.307803    1205.0429      1.5049911   1.194418337   0.5297660                50                    132
# 3 0.7291883   2088.368783    2088.3688      0.4344569   0.794622032   0.8924053                50                     14
# 6 0.7291883   2088.368783    2088.3688      0.4344569   0.794622032   0.8924053                50                     14
# 5 0.1385699 333370.032902  333370.0329      1.9641374   0.004393986   0.4167017                50                     27
# 4 0.1395601 350083.518163  350083.5182      1.9789938   0.004496484   0.4412216                50                     26
#   estimated_over_true_support_size TP FP  TN FN precision recall specificity  accuracy  F1_score       gic      ebic       bic
# 1                             0.98 43  6 444  7 0.8775510   0.86   0.9866667 0.9740000 0.8686869  3961.212  4015.016  3704.086
# 2                             1.10 42 13 438  8 0.7636364   0.84   0.9711752 0.9580838 0.8000000  4343.159  4403.439  4054.974
# 7                             2.58 48 81 370  2 0.3720930   0.96   0.8203991 0.8343313 0.5363128  5365.541  5505.476  4696.539
# 9                             2.64 47 85 366  3 0.3560606   0.94   0.8115299 0.8243513 0.5164835  5654.297  5797.462  4969.857
# 8                             2.64 47 85 366  3 0.3560606   0.94   0.8115299 0.8243513 0.5164835  5665.472  5808.636  4981.032
# 3                             0.28 13  1 449 37 0.9285714   0.26   0.9977778 0.9240000 0.4062500  6669.590  6685.731  6592.452
# 6                             0.28 13  1 449 37 0.9285714   0.26   0.9977778 0.9240000 0.4062500  6669.590  6685.731  6592.452
# 5                             0.54 14 13 437 36 0.5185185   0.28   0.9711111 0.9020000 0.3636364 11632.872 11663.002 11488.882
# 4                             0.52 11 15 435 39 0.4230769   0.22   0.9666667 0.8920000 0.2894737 11620.942 11649.996 11482.093
#    aic      mbic
# 1  3493.355  4187.117
# 2  3818.956  4596.192
# 7  4148.640  5952.938
# 9  4409.314  6255.250
# 8  4420.489  6266.425
# 3  6529.233  6737.361
# 6  6529.233  6737.361
# 5 11370.872 11759.379
# 4 11368.299 11742.930

# TODO: should also include timings (but L0glm is much faster than all the rest)


## 3.2 SINGLE CHANNEL GAUSSIAN ####
n <- p <- 500
sim.gaus <- simul(method = 1, n = n, p = p, seed = 2)

### 3.2.1 L0GLM VIA CHOLESKY ####

system.time(tune.chol.sc <- L0glm.tune(X=sim.gaus$X, y=sim.gaus$y, IC = "gic", 
                                       plot = FALSE, coefs.true = sim.gaus$beta_true, 
                                       relaxed = FALSE, family = gaussian(identity), 
                                       nonnegative = TRUE, 
                                       solver = "eigen", method = 7L))

system.time(L0glm.chol <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                    family = gaussian(identity), 
                                    lambda = "gic", # tune.chol.sc$best.lambda,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "eigen",
                                    method = 7L)) 
# 0.01s for n=p=500
cor(as.vector(as.vector(L0glm.chol$coefficients)), as.vector(sim.gaus$beta_true)) # 0.986257
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.chol$coefficients) != 0)) # TP = 45
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.chol$coefficients) != 0)) # FP = 8
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0glm.chol$coefficients) == 0)) # TN = 442
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0glm.chol$coefficients) == 0)) # FN = 5

plot(x=L0glm.chol$coefficients, y=sim.gaus$beta_true, pch=16, col='blue',
     main="L0 glm identity link gaussian with lambda tuned to optimise gic",
     xlab="estimatedcoefficients", ylab="true coefficients")


### 3.2.2 LASSO USING GLMNET ####
system.time(glmnet.sc.gaus <- fit.GLM(sim = sim.gaus, method = 2, 
                                      family = gaussian(identity), plot = TRUE)) 
# 0.22s for n=p=500 
cor(as.vector(glmnet.sc.gaus[-1]), as.vector(sim.gaus$beta_true)) # 0.9611292
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(glmnet.sc.gaus[-1]) != 0)) # TP = 49
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(glmnet.sc.gaus[-1]) != 0)) # FP = 68
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(glmnet.sc.gaus[-1]) == 0)) # TN = 382
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(glmnet.sc.gaus[-1]) == 0)) # FN = 1


### 3.2.3 ADAPTIVE LASSO USING GLMNET ####
system.time(glmnet.al.sc.gaus <- fit.GLM(sim = sim.gaus, method = 3, 
                                         family = gaussian(identity), plot = T)) 
# 0.20s for n=p=500 
cor(as.vector(glmnet.al.sc.gaus[-1]), as.vector(sim.gaus$beta_true)) # 0.9626742
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(glmnet.al.sc.gaus[-1]) != 0)) # TP = 37
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(glmnet.al.sc.gaus[-1]) != 0)) # FP = 12
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(glmnet.al.sc.gaus[-1]) == 0)) # TN = 438
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(glmnet.al.sc.gaus[-1]) == 0)) # FN = 13


### 3.2.4 L0LEARN ####

system.time(L0learn.sc.gaus <- fit.GLM(sim = sim.gaus, method = 4, 
                                       family = gaussian(identity), plot = T))  
#  0.24s for n=p=500
cor(as.vector(L0learn.sc.gaus), as.vector(sim.gaus$beta_true)) # 0.914644
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0learn.sc.gaus) != 0)) # TP = 34
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0learn.sc.gaus) != 0)) # FP = 33
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0learn.sc.gaus) == 0)) # TN = 417
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0learn.sc.gaus) == 0)) # FN = 16


### 3.2.5 L0LEARN AFTER LASSO ####
# Better than regular L0Learn & faster
system.time(L0learn.after.lasso.sc.gaus <- fit.GLM(sim = sim.gaus, method = 5, 
                                                   family = gaussian(identity), plot = T))  
# 0.66s for n=p=500 
cor(as.vector(L0learn.after.lasso.sc.gaus), as.vector(sim.gaus$beta_true)) # 0.9429308
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0learn.after.lasso.sc.gaus) != 0)) # TP = 44
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0learn.after.lasso.sc.gaus) != 0)) # FP = 47
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(L0learn.after.lasso.sc.gaus) == 0)) # TN = 403
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(L0learn.after.lasso.sc.gaus) == 0)) # FN = 6


### 3.2.6 ABESS ####
system.time(abess.sc.gaus <- fit.GLM(sim = sim.gaus, method = 6, 
                                     family = gaussian(identity), 
                                     plot = T, IC = "gic"))  
# 0.99s for n=p=500
cor(as.vector(abess.sc.gaus), as.vector(sim.gaus$beta_true)) # 0.8706258
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(abess.sc.gaus) != 0)) # TP = 28
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(abess.sc.gaus) != 0)) # FP = 20
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(abess.sc.gaus) == 0)) # TN = 430
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(abess.sc.gaus) == 0)) # FN = 22


### 3.2.7 ABESS AFTER LASSO ####
system.time(abess.after.lasso.sc.gaus <- fit.GLM(sim = sim.gaus, method = 7, 
                                                 family = gaussian(identity), 
                                                 plot = T, IC = "gic"))  
# 0.15s for n=p=500
cor(as.vector(abess.after.lasso.sc.gaus), as.vector(sim.gaus$beta_true))  # 0.8814479
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(abess.after.lasso.sc.gaus) != 0)) # TP = 29
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(abess.after.lasso.sc.gaus) != 0)) # FP = 19
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(abess.after.lasso.sc.gaus) == 0)) # TN = 431
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(abess.after.lasso.sc.gaus) == 0)) # FN = 21


### 3.2.8 MCP USING ORDINIS ####
system.time(ordinis.mcp.sc.gaus <- fit.GLM(sim = sim.gaus, method = 8, 
                                           family = gaussian(identity), 
                                           plot = T, IC = "gic"))  
# 0.30s for n=p=500
cor(as.vector(ordinis.mcp.sc.gaus[-1]), as.vector(sim.gaus$beta_true)) # 0.7612053
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(ordinis.mcp.sc.gaus[-1]) != 0)) # TP = 24
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(ordinis.mcp.sc.gaus[-1]) != 0)) # FP = 41
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(ordinis.mcp.sc.gaus[-1]) == 0)) # TN = 409
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(ordinis.mcp.sc.gaus[-1]) == 0)) # FN = 26


### 3.2.9 SCAD USING ORDINIS ####
system.time(ordinis.scad.sc.gaus <- fit.GLM(sim = sim.gaus, method = 9, 
                                            family = gaussian(identity), 
                                            plot = T, IC = "gic"))  
# 0.38s for n=p=500
cor(as.vector(ordinis.scad.sc.gaus[-1]), as.vector(sim.gaus$beta_true)) # 0.9205098
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(ordinis.scad.sc.gaus[-1]) != 0)) # TP = 34
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(ordinis.scad.sc.gaus[-1]) != 0)) # FP = 22
sum((as.vector(sim.gaus$beta_true) == 0) & (as.vector(ordinis.scad.sc.gaus[-1]) == 0)) # TN = 428
sum((as.vector(sim.gaus$beta_true) != 0) & (as.vector(ordinis.scad.sc.gaus[-1]) == 0)) # FN = 16


### 3.2.10 BENCHMARKS TABLE ####

X <- sim.gaus$X
y <- sim.gaus$y
W <- as.matrix(rep(1,ncol(X)))
IC <- "gic"

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
#                           method   p   n m   deviance deviance_true    loglik      MD_coefs MAD_coefs sMAPE_coefs bias_coefs
# 3 nnL0Learn_after_glmnet_nnlasso 501 500 1   13547.22      6733.678 -1125.301  -0.007830595 113.60805  0.23134485  0.9414228
# 1                          L0glm 501 500 1   30617.34     19406.542 -1291.148   0.005731796  38.64495  0.05595009  0.9882487
# 5                      nnL0Learn 501 500 1  232710.80    221496.014 -1812.205   0.411679765 146.54041  0.20577746  0.9247740
# 9                 nnSCAD_ordinis 502 500 1  508208.93    491714.228 -1995.479   0.601437492 131.55853  0.15941428  0.9207443
# 8                  nnMCP_ordinis 502 500 1 1273567.56   1249035.323 -2234.150   1.123362852 260.28405  0.27135676  0.8101574
# 4     abess_after_glmnet_nnlasso 501 500 1 3258284.46   3269706.782 -2452.995   0.940147627 176.14620  0.16691999  0.8598522
# 2                nnalasso_glmnet 502 500 1 3665597.68   3667508.250 -2482.443  -6.096004919  95.93336  0.11343157  0.9764912
# 7                 nnlasso_glmnet 502 500 1 4338167.47   4350827.911 -2592.558 -11.601690959 125.95294  0.33233782  1.0491874
# 6                          abess 501 500 1 7116162.10   7116773.335 -2648.287   2.076538493 186.68587  0.17576504  0.8381777
#    R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true support_size_estimated
# 3 0.8942910   105.25784     518.9192      1.1392574     0.9414228   0.9093948                50                     91
# 1 0.9740121    93.52144     271.1076      0.3391518     0.9882487   0.9792205                50                     53
# 5 0.8441756   473.67692     622.9805      1.0580407     0.9247740   0.9178036                50                     67
# 9 0.8543947   576.38337     605.9764      0.8547599     0.9207443   0.9355237                50                     56
# 8 0.5979304   959.01220    1055.8357      1.2915344     0.8101574   0.8073822                50                     65
# 4 0.7867479   916.22988     928.4196      0.8637499     0.8598522   0.8957397                50                     48
# 2 0.9301852   223.27211     695.2690      0.6291677     0.9764912   0.9447399                50                     49
# 7 0.9274381   -39.80944     529.4674      1.4059936     1.0491874   0.9285669                50                    117
# 6 0.7684228   979.49260     986.7825      0.9142191     0.8381777   0.8803839                50                     48
#   estimated_over_true_support_size TP FP  TN FN precision recall specificity  accuracy  F1_score      gic     ebic      bic      aic
# 3                             1.82 44 47 403  6 0.4835165   0.88   0.8955556 0.8940000 0.6241135 3295.459 3394.458 2822.347 2434.603
# 1                             1.06 45  8 442  5 0.8490566   0.90   0.9822222 0.9740000 0.8737864 3195.581 3253.689 2917.884 2690.295
# 5                             1.34 34 33 417 16 0.5074627   0.68   0.9266667 0.9020000 0.5811966 4396.695 4469.869 4047.004 3760.410
# 9                             1.12 34 22 429 16 0.6071429   0.68   0.9512195 0.9241517 0.6415094 4638.523 4699.879 4345.191 4104.958
# 8                             1.30 24 41 410 26 0.3692308   0.48   0.9090909 0.8662675 0.4173913 5218.112 5289.156 4878.465 4600.301
# 4                             0.96 29 19 431 21 0.6041667   0.58   0.9577778 0.9200000 0.5918367 5462.490 5515.218 5210.506 5003.990
# 2                             0.98 37 12 439 13 0.7551020   0.74   0.9733925 0.9500998 0.7474747 5532.924 5586.745 5275.616 5064.885
# 7                             2.34 49 68 383  1 0.4188034   0.98   0.8492239 0.8622754 0.5868263 6525.687 6652.705 5918.440 5421.116
# 6                             0.96 28 20 430 22 0.5833333   0.56   0.9555556 0.9160000 0.5714286 5853.074 5905.801 5601.090 5394.574
#       mbic
# 3 3711.124
# 1 3439.558
# 5 4703.926
# 9 4896.074
# 8 5516.329
# 4 5683.877
# 2 5758.846
# 7 7058.864
# 6 6074.460

# TODO: should also include timings (but L0glm is much faster than all the rest)
# deviance relative to ground truth best for nnL0Learn_after_glmnet_nnlasso but
# F1_score & gic better for L0glm


## 3.3 MULTICHANNEL POISSON ####
n <- p <- 500
sim.pois.mc <- simul(method = 4, n = n, p = p, seed = 2)

### 3.3.1 RELAXED GROUP L0GLM VIA CHOLESKY ####
system.time(tune.chol.mc <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, IC = "gic", 
                                       plot = FALSE, coefs.true = sim.pois.mc$beta_true, 
                                       family = poisson(identity), 
                                       nonnegative = TRUE, 
                                       grouped = TRUE,
                                       relaxed = TRUE, 
                                       parallelize = FALSE, 
                                       solver = "eigen", 
                                       method = 7L)) 
# 345s with parallelize=F, returns error with parallelize=TRUE
# probably something in L0glm_2_fit.R?

system.time(L0glm.chol.pois <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                         family = poisson(identity), 
                                         lambda = "gic",
                                         nonnegative = TRUE,
                                         grouped = TRUE,
                                         relaxed = TRUE,
                                         parallelize = FALSE,
                                         solver = "eigen",
                                         method = 7L)) 
# 21s

R2 = round(cor(as.vector(rowMeans(sim.pois.mc$beta_true)), as.vector(rowMeans(L0glm.chol.pois$coefficients))), digits = 4)^2 # Calculate R2
title = paste("true vs estimated coefficients\nnonnegative relaxed L0 grouped GLM\nR2=", R2)
plot(as.vector(rowMeans(L0glm.chol.pois$coefficients)), as.vector(rowMeans(sim.pois.mc$beta_true)), pch=16, 
     xlab="estimated coefficients",
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


### 3.3.2 BEST SUBSET OF GROUP SELECTION (BSGS) USING ABESS ####
system.time(bsgs.abess.pois <- fit.GLM.mc(sim = sim.pois.mc, method = 1, 
                                             family = poisson(identity), 
                                             plot = T, IC = "gic")) 
# 4.4s for n=p=500
# NOTE: this is now fit as if noise was Gaussian
# TODO also compare what solution & timings you get with
# https://skscope.readthedocs.io/en/latest/gallery/GeneralizedLinearModels/poisson-identity-link.html
# maybe implement this as an additional fitting method in 
# fit.GLM.mc?

cor(as.vector(as.vector(sim.pois.mc$beta_true)), as.vector(bsgs.abess.pois)) # 0.9060618
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(bsgs.abess.pois)!=0) # TP = 9
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(bsgs.abess.pois)!=0) # FP = 1
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(bsgs.abess.pois)==0) # TN = 449
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(bsgs.abess.pois)==0) # FN = 41



### 3.3.3 MULTIVARIATE GROUP LASSO USING GLMNET ####
system.time(glmnet.mc.pois <- fit.GLM.mc(sim = sim.pois.mc, method = 2, 
                                    family = poisson(identity), plot = T, 
                                    IC = "gic"))
# 19.6s for n=p=500

image(x=1:nrow(sim.pois.mc$y), y=1:ncol(sim.pois.mc$y), z=as.matrix(glmnet.mc.pois)^0.1, 
      col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", 
      main="glmnet nonnegative multivariate group LASSO\n(red = true peaks, cyan = estimated peaks)")
abline(v=(1:nrow(sim.pois.mc$X))[as.vector(rowMax(sim.pois.mc$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim.pois.mc$X))[as.vector(rowMin(glmnet.mc.pois)!=0)], col="cyan")


cor(as.vector(as.vector(sim.pois.mc$beta_true)), as.vector(glmnet.mc.pois)) # 0.9888263
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(glmnet.mc.pois)!=0) # 46
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(glmnet.mc.pois)!=0) # 22
sum(rowMin(sim.pois.mc$beta_true)==0 & rowMin(glmnet.mc.pois)==0) # 428
sum(rowMin(sim.pois.mc$beta_true)!=0 & rowMin(glmnet.mc.pois)==0) # 4


### 3.3.4 BENCHMARKS TABLE ####

W <- (1/(sim.pois.mc$y+0.1))

benchmarks_3_3 = rbind(
  data.frame(method="Relaxed_group_L0glm_poisson", do.call(cbind, get.metrics(beta = as.matrix(L0glm.chol.pois$coefficients),  # TODO: make sure that metrics work with either dense or sparse coefs
                                                                        X = sim.pois.mc$X, y = sim.pois.mc$y,
                                                                        weights = W, 
                                                                        family = gaussian(identity), 
                                                                        beta_true = sim.pois.mc$beta_true) ) ),
  data.frame(method="BSGS_abess_poisson", do.call(cbind, get.metrics(beta = bsgs.abess.pois, 
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
#                             method   p   n   m   deviance deviance_true    loglik     MD_coefs MAD_coefs sMAPE_coefs bias_coefs
# 1      Relaxed_group_L0glm_poisson 501 500 160   67399.34      6737.463 -102070.4 -0.005062176 0.4429137  0.02643714  1.0012997
# 3 multivariate_group_LASSO_poisson 501 500 160  168497.73    113862.261 -141536.9 -0.561630507 2.3569584  0.11016033  1.0299012
# 2               BSGS_abess_poisson 501 500 160 2415489.28   2414515.360 -242537.3 -3.180608680 6.2203325  0.09183572  0.9146294
#    R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true support_size_estimated
# 1 0.9991303    1.156908     7.993593      0.3821764     1.0012997   0.9993185                50                     50
# 3 0.9778393   -5.446269    27.066958      1.2092166     1.0299012   0.9783603                50                     68
# 2 0.8211056  106.532307   107.910290      1.0483753     0.9146294   0.9430487                50                     10
#   estimated_over_true_support_size TP FP  TN FN precision recall specificity accuracy F1_score      gic     ebic      bic      aic
# 1                             1.00 49  1 449  1 0.9800000   0.98   0.9977778    0.996 0.980000 247627.4 251747.7 227936.6 211798.9
# 3                             1.36 46 22 428  4 0.6764706   0.92   0.9511111    0.948 0.779661 358530.7 365680.1 324363.7 296361.9
# 2                             0.20  9  1 449 41 0.9000000   0.18   0.9977778    0.916 0.300000 497965.0 499186.3 492128.2 487344.6
#       mbic
# 1 264927.2
# 3 388548.9
# 2 503093.0




## 3.4 MULTICHANNEL GAUSSIAN ####
n <- p <- 500
sim.gaus.mc <- simul(method = 3, n = n, p = p, seed = 1)

### 3.4.1 RELAXED GROUP L0GLM VIA CHOLESKY ####

system.time(tune.chol.mc <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, IC = "gic", 
                                       plot = FALSE, coefs.true = sim.gaus.mc$beta_true, 
                                       family = gaussian(identity), 
                                       nonnegative = TRUE, 
                                       grouped = TRUE,
                                       relaxed = TRUE, 
                                       parallelize = FALSE, # TODO fix error in L0glm.tune when set to TRUE
                                       solver = "eigen", 
                                       method = 7L)) 
# 278s with parallelize=F else returns error, best lambda 3.16227766016838
# TODO fix error when option parallelize = TRUE
# Error in `[<-`(`*tmp*`, variable.indices, , value = Matrix(do.call(cbind,  : 
#                                                                      subscript out of bounds
#                                                                    In addition: There were 50 or more warnings (use warnings() to see the first 50)
#                                                                    Called from: L0glm.fit(X = X, y = y, family = family, weights = weights, lambda = lambda, 
#                                                                                           dispersion = dispersion, start = start, offset = offset, 
#                                                                                           ...)
# > traceback()
# 3: drop0(X)
# 2: as.simple_triplet_matrix(drop0(X))
# 1: rowMin(L0glm.chol.pois$coefficients) # not sure where L0glm.chol.pois$coefficients comes from...

system.time(L0glm.chol.gaus <- L0glm.fit(sim.gaus.mc$X, sim.gaus.mc$y, 
                                    family = gaussian(identity), 
                                    lambda = "gic", 
                                    nonnegative = TRUE,
                                    grouped = TRUE,
                                    relaxed = TRUE,
                                    parallelize = FALSE,
                                    solver = "eigen",
                                    method = 7L)) 
# 18s for n=p=500 with parallelize = FALSE else returns some warnings
# TODO fix warnings when parallelize = TRUE
# Warning messages:
# 1: In (function (..., deparse.level = 1)  :
#         number of rows of result is not a multiple of vector length (arg 6)
       
R2 = round(cor(as.vector(rowMeans(sim.gaus.mc$beta_true)), as.vector(rowMeans(L0glm.chol.gaus$coefficients))), digits = 4)^2 # Calculate R2
title = paste("Relaxed nonnegative L0 group GLM\nR2=", R2)
plot(as.vector(rowMeans(L0glm.chol.gaus$coefficients)), 
     as.vector(rowMeans(sim.gaus.mc$beta_true)), pch=16, 
     xlab="estimated coefficients",
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


### 3.4.2 BEST SUBSET OF GROUP SELECTION (BSGS) USING ABESS ####
system.time(bsgs.abess.gaus <- fit.GLM.mc(sim = sim.gaus.mc, method = 1, 
                                          family = gaussian(identity), plot = T, 
                                          IC = "gic")) 
# 4.2s for n=p=500
# TODO post this as a reproducible example in https://github.com/abess-team/abess/issues/510
# and ask if the solution quality could be improved here somehow
# TODO also compare what solution & timings you get with
# https://skscope.readthedocs.io/en/latest/gallery/GeneralizedLinearModels/poisson-identity-link.html

cor(as.vector(as.vector(sim.gaus.mc$beta_true)), as.vector(bsgs.abess.gaus)) # 0.6662117
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(bsgs.abess.gaus)!=0) # TP = 17
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(bsgs.abess.gaus)!=0) # FP = 7
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(bsgs.abess.gaus)==0) # TN = 443
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(bsgs.abess.gaus)==0) # FN = 33


# glmnet nonnegative mgaussian multivariate group LASSO fit
system.time(glmnet.mc.gaus <- fit.GLM.mc(sim = sim.gaus.mc, method = 2, 
                                    family = gaussian(identity), plot = T, 
                                    IC = "gic"))
# 29s for n=p=500

image(x=1:nrow(sim.gaus.mc$y), y=1:ncol(sim.gaus.mc$y), z=as.matrix(glmnet.mc.gaus)^0.1, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", 
      main="nonnegative multivariate group LASSO\n(red = true peaks, cyan = estimated peaks)")
abline(v=(1:nrow(sim.gaus.mc$X))[as.vector(rowMax(sim.gaus.mc$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim.gaus.mc$X))[as.vector(rowMin(glmnet.mc.gaus)!=0)], col="cyan")


cor(as.vector(as.vector(sim.gaus.mc$beta_true)), as.vector(glmnet.mc.gaus)) # 0.9815128
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(glmnet.mc.gaus)!=0) # 47
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(glmnet.mc.gaus)!=0) # 32
sum(rowMin(sim.gaus.mc$beta_true)==0 & rowMin(glmnet.mc.gaus)==0) # 418
sum(rowMin(sim.gaus.mc$beta_true)!=0 & rowMin(glmnet.mc.gaus)==0) # 3


### 3.4.3 BENCHMARKS TABLE ####

benchmarks_3_4 = rbind(
  data.frame(method="Relaxed_group_L0glm_gaussian", do.call(cbind, get.metrics(beta = as.matrix(L0glm.chol.gaus$coefficients),  # TODO: allow coefs to be either dense or sparse matrix in get.metrics
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
#                              method   p   n   m  deviance deviance_true    loglik    MD_coefs MAD_coefs sMAPE_coefs bias_coefs
# 1      Relaxed_group_L0glm_gaussian 501 500 165   1125238        207933 -231325.7  0.09550816  0.254989  0.09548499  1.0000758
# 3 multivariate_group_LASSO_gaussian 501 500 165  11803607      10895240 -332418.1 -0.22903666  2.471618  0.21089491  1.0481612
# 2               BSGS_abess_gaussian 501 500 165 312127788     311389742 -460225.2 -0.73411430  7.242264  0.12506039  0.6556086
#    R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true support_size_estimated
# 1 0.9999020   1.2232292      2.92006       1.067991     1.0000758   0.9999054                50                     50
# 3 0.9635367   0.4048812     16.23747       1.495982     1.0481612   0.9707416                50                     79
# 2 0.4456978  69.6143367     70.57427       1.393076     0.6556086   0.7048314                50                     24
#   estimated_over_true_support_size TP FP  TN FN precision recall specificity accuracy  F1_score      gic     ebic      bic      aic
# 1                             1.00 50  0 450  0 1.0000000   1.00   1.0000000     1.00 1.0000000 541924.1 549435.2 506029.3 476611.3
# 3                             1.58 47 32 418  3 0.5949367   0.94   0.9288889     0.93 0.7286822 791104.9 803068.7 733930.2 687072.2
# 2                             0.48 17  7 443 33 0.7083333   0.34   0.9844444     0.92 0.4594595 963959.5 968082.0 944258.5 928112.4
#       mbic
# 1 573460.5
# 3 841337.1
# 2 981268.4