# BENCHMARKS L0GLM: SINGLE CHANNEL VS MULTICHANNEL DECONVOLUTION OF BLURRED SPIKE TRAIN ####

library(remotes)
# remotes::install_github("jaredhuling/ordinis")
# remotes::install_github("tomwenseleers/L0glm/L0glm")
# devtools::load_all(".")
library(L0glm)
library(export)
library(Matrix)
library(ordinis)
library(L0Learn) # use at least version 2
library(abess)
library(glmnet)

# 1. SINGLE CHANNEL, SPARSE  ####
n <- p <- 10000 
npeaks <- 1000
sim.gaus <- simul(method = 1, n = n, p = p, s = npeaks/p, sparse = T) 
sim.pois <- simul(method = 2, n = n, p = p, s = npeaks/p, sparse = T) 

## 1.1. L0GLM FIT SINGLE CHANNEL WITH POISSON NOISE ####

## tuning plot (full range)
# full lambda grid range with warm.start = 0.01
# here preset lambda coincides with best lambda
tune_full <- L0glm.tune(X=sim.pois$X, 
                   y=sim.pois$y, 
                   IC = "gic", 
                   nlambda = 10,
                   prop_range = 1, # full range up to lambda_max
                   warm.start = 0, # 0.01, # choose 0 for no warm starts
                   plot = TRUE, 
                   coefs.true = sim.pois$beta_true, 
                   relaxed = FALSE, 
                   family = poisson(identity), 
                   nonnegative = TRUE, 
                   solver = "eigen", 
                   method = 7L)
tune_full$lambda_best # 20.44997 (best lambda tune)
tune_full$lambda_preset # 20.44997 (lambda preset)
min(tune_full$metrics$gic) # 102767

## tuning plot (zoomed in)
# with lambda grid more zoomed in around lambda_preset preset
tune_zoomed <- L0glm.tune(X=sim.pois$X, 
           y=sim.pois$y, 
           IC = "gic", 
           nlambda = 10,
           prop_range = 0.2, # to zoom in around lambda_preset preset
           warm.start = 0, # 0.01, 
           plot = TRUE, 
           coefs.true = sim.pois$beta_true, 
           relaxed = FALSE, 
           family = poisson(identity), 
           nonnegative = TRUE, 
           solver = "eigen", 
           method = 7L)
# here best lambda (dashed line) is very close to lambda preset expected 
# to optimise GIC (dotted line) so tuning in practice unnecessary
tune_zoomed$lambda_best # 11.49116 (best lambda tune)
tune_zoomed$lambda_preset # 20.44997 (lambda preset)
min(tune_zoomed$metrics$gic) # 102133.6

# tuning lambda with different solvers
system.time(tune.chol.sc.pois <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "gic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = FALSE, family = poisson(identity), nonnegative = TRUE, solver = "eigen", method = 7L))
# Best lambda by GIC = 20.4, time = 2.16s for n=p=10000 
system.time(tune.CG.sc.pois <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "gic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = FALSE, family = poisson(identity), nonnegative = TRUE, solver = "eigen", method = 5L))
# Best lambda by GIC = 20.4, time = 8.55s for n=p=10000, but eigen CG solver failed in IRLS iterations 1 to 3, trying Cholesky LLT solver as fallback.
system.time(tune.LSCG.sc.pois <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "gic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = FALSE, family = poisson(identity), nonnegative = TRUE, solver = "eigen", method = 3L))
# Best lambda by GIC = 20.4, time = 17.53s for n=p=10000, but eigen CG solver failed in IRLS iterations 1 to 3, trying Cholesky LLT solver as fallback.
system.time(tune.osqp.sc.pois <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "gic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = FALSE, family = poisson(identity), nonnegative = TRUE, solver = "osqp"))
# Best lambda by GIC = 20.4, time = 78.06s for n=p=10000
system.time(tune.glmnet.sc.pois <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "gic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = FALSE, family = poisson(identity), nonnegative = TRUE, solver = "glmnet"))
# Best lambda by GIC = 20.4, time = 35.92s for n=p=10000
system.time(tune.QR.sc.pois <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "gic", plot = FALSE, coefs.true = sim.pois$beta_true, relaxed = FALSE, family = poisson(identity), nonnegative = TRUE, solver = "eigen", method = 1L))
# Best lambda by GIC = 20.4, time = 182.98s for n=p=10000
# -> all solvers give same optimal lambda value, but Cholesky fastest

# Chosen lambda also matches expected lambda value that
# was expected to optimise GIC, so tuning is in fact not needed
L0glm:::preset.lambda(IC="gic", p, n) # 20.44997 (lambda preset)
tune.chol.sc.pois$lambda_preset # 20.44997 (lambda preset)
tune.chol.sc.pois$lambda_best # 20.44997


#### 1.1.1 L0GLM FIT WITH CHOLESKY SOLVER ####
system.time(L0glm.chol.pois <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                         family = poisson(identity), 
                                         lambda = "gic", # or tune_zoomed$lambda_best
                                         nonnegative = TRUE,
                                         relaxed = FALSE,
                                         # verbose = TRUE,
                                         solver = "eigen",
                                         method = 7L)) 
# 0.25s with n=p=10000

#### 1.1.2 L0GLM FIT WITH CG SOLVER ####
system.time(L0glm.CG.pois <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                       family = poisson(identity), 
                                       lambda = "gic", # or tune_zoomed$lambda_best
                                       nonnegative = TRUE,
                                       relaxed = FALSE,
                                       solver = "eigen",
                                       method = 5L)) 
# 0.4s with n=p=10000

#### 1.1.3 L0GLM FIT WITH LSCG SOLVER ####
system.time(L0glm.LSCG.pois <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                         family = poisson(identity), 
                                         lambda = "gic", # or tune_zoomed$lambda_best
                                         nonnegative = TRUE,
                                         relaxed = FALSE,
                                         solver = "eigen",
                                         method = 3L)) 
# 0.70s with n=p=10000

#### 1.1.4 L0GLM FIT WITH OSQP SOLVER ####
system.time(L0glm.osqp.pois <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                    family = poisson(identity), 
                                    lambda = "gic", # or tune_zoomed$lambda_best
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "osqp")) 
# 5.78s with n=p=10000

#### 1.1.5 L0GLM FIT WITH glmnet::bigGlm SOLVER ####
system.time(L0glm.glmnet.pois <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                      family = poisson(identity), 
                                      lambda = "gic", # or tune_zoomed$lambda_best
                                      nonnegative = TRUE,
                                      relaxed = FALSE,
                                      solver = "glmnet")) 
# 2.09s with n=p=10000

#### 1.1.6 L0GLM FIT WITH QR SOLVER ####
system.time(L0glm.QR.pois <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                  family = poisson(identity), 
                                  lambda = "gic", # or tune_zoomed$lambda_best
                                  nonnegative = TRUE,
                                  relaxed = FALSE,
                                  solver = "eigen",
                                  method = 1L)) 
# 78.55s with n=p=10000

#### 1.1.7 L0GLM FIT WITH NNLS SOLVER (very slow for larger problem size) ####
# system.time(L0glm.nnls.pois <- L0glm.fit(sim.pois$X, sim.pois$y, 
#                                   family = poisson(identity), 
#                                   lambda = "gic", 
#                                   nonnegative = TRUE,
#                                   relaxed = FALSE,
#                                   solver = "nnls")) 
# XXs with n=p=10000
# 
#### 1.1.8 PLAIN NNLS FIT WITHOUT REGULARISATION ####
# system.time(L0glm.nnlsnoreg.pois <- L0glm.fit(sim.pois$X, sim.pois$y, 
#                                     family = poisson(identity), 
#                                     lambda = 0,
#                                     nonnegative = TRUE,
#                                     relaxed = FALSE,
#                                     solver = "nnls",
#                                     maxit = 1)) 
# XXs with n=p=10000



### 1.1.9 BENCHMARKS TABLE ####

X = sim.pois$X
y = sim.pois$y
benchmarks_1_1 = rbind(
  data.frame(method="L0glm_osqp", do.call(cbind, L0glm:::get.metrics(beta = L0glm.osqp.pois$coefficients,  # some metrics
                                                             X = X, y = y, 
                                                             family = poisson(identity), 
                                                             beta_true = sim.pois$beta_true))),
  data.frame(method="L0glm_bigglm", do.call(cbind, L0glm:::get.metrics(beta = L0glm.glmnet.pois$coefficients,  # some metrics
                                                               X = X, y = y, 
                                                               family = poisson(identity), 
                                                               beta_true = sim.pois$beta_true))),
  data.frame(method="L0glm_qr", do.call(cbind, L0glm:::get.metrics(beta = L0glm.QR.pois$coefficients,  # some metrics
                                                           X = X, y = y, 
                                                           family = poisson(identity),
                                                           beta_true = as.vector(sim.pois$beta_true)) )),
  data.frame(method="L0glm_cholesky", do.call(cbind, L0glm:::get.metrics(beta = L0glm.chol.pois$coefficients,  # some metrics
                                                                 X = X, y = y, 
                                                                 family = poisson(identity),
                                                                 beta_true = sim.pois$beta_true) )),
  data.frame(method="L0glm_LSCG", do.call(cbind, L0glm:::get.metrics(beta = L0glm.LSCG.pois$coefficients,  # some metrics
                                                             X = X, y = y, 
                                                             family = poisson(identity),
                                                             beta_true = sim.pois$beta_true) )),
  data.frame(method="L0glm_CG", do.call(cbind, L0glm:::get.metrics(beta = L0glm.CG.pois$coefficients,  # some metrics
                                                           X = X, y = y, 
                                                           family = poisson(identity), 
                                                           beta_true = sim.pois$beta_true) ))
)

benchmarks_1_1 = benchmarks_1_1[order(benchmarks_1_1[,"deviance_true"]),] # sort by deviance relative to ground truth
benchmarks_1_1
#           method     p     n m deviance deviance_true    loglik
# 2   L0glm_bigglm 10000 10000 1 8543.501      2152.794 -40633.23
# 3       L0glm_qr 10000 10000 1 9348.514      2987.828 -41035.74
# 4 L0glm_cholesky 10000 10000 1 9348.514      2987.828 -41035.74
# 5     L0glm_LSCG 10000 10000 1 9348.525      2987.833 -41035.74
# 6       L0glm_CG 10000 10000 1 9348.533      2987.846 -41035.75
# 1     L0glm_osqp 10000 10000 1 9462.469      3075.212 -41092.71
#      MD_coefs MAD_coefs sMAPE_coefs bias_coefs  R2_coefs MD_selcoefs
# 2 -0.02795328  73.68727  0.06511543  0.9858442 0.8615312    242.5099
# 3 -0.02862166  80.26889  0.07160071  0.9928768 0.9093261    262.8916
# 4 -0.02862166  80.26889  0.07160071  0.9928768 0.9093261    262.8916
# 5 -0.02861981  80.26933  0.07160092  0.9928763 0.9093256    262.8917
# 6 -0.02862117  80.26951  0.07160104  0.9928761 0.9093254    262.8916
# 1 -0.02869495  80.65774  0.07307045  0.9914587 0.9086122    270.6447
#   MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs
# 2     491.1486      0.3557314     0.9858442   0.8869727
# 3     529.9963      0.3992165     0.9928768   0.9399447
# 4     529.9963      0.3992165     0.9928768   0.9399447
# 5     530.0008      0.3992185     0.9928763   0.9399442
# 6     530.0025      0.3992198     0.9928761   0.9399440
# 1     526.0850      0.4058344     0.9914587   0.9403085
#   support_size_true support_size_estimated
# 2              1000                   1004
# 3              1000                   1012
# 4              1000                   1012
# 5              1000                   1012
# 6              1000                   1012
# 1              1000                   1012
#   estimated_over_true_support_size  TP  FP   TN  FN precision recall
# 2                            1.004 853 151 8849 147 0.8496016  0.853
# 3                            1.012 844 168 8832 156 0.8339921  0.844
# 4                            1.012 844 168 8832 156 0.8339921  0.844
# 5                            1.012 844 168 8832 156 0.8339921  0.844
# 6                            1.012 844 168 8832 156 0.8339921  0.844
# 1                            1.012 840 172 8828 160 0.8300395  0.840
#   specificity accuracy  F1_score      gic      ebic      bic      aic
# 2   0.9832222   0.9702 0.8512974 101798.2  99760.82 90513.64 83274.46
# 3   0.9813333   0.9676 0.8389662 102766.8 100713.20 91392.34 84095.47
# 4   0.9813333   0.9676 0.8389662 102766.8 100713.20 91392.34 84095.47
# 5   0.9813333   0.9676 0.8389662 102766.8 100713.21 91392.35 84095.48
# 6   0.9813333   0.9676 0.8389662 102766.9 100713.22 91392.35 84095.49
# 1   0.9808889   0.9668 0.8349901 102880.8 100827.16 91506.29 84209.43
#       mbic
# 2 106224.3
# 3 107228.2
# 4 107228.2
# 5 107228.2
# 6 107228.2
# 1 107342.2

# -> L0glm_bigglm gives most accurate solution but more
# downward bias on coefficients & ca 15x slower than cholesky
# eigen solvers all equally good & cholesky the fastest
# osqp solver the worst
# -> Cholesky is the best


## 1.2. L0GLM FIT SINGLE CHANNEL WITH GAUSSIAN NOISE ####

## tuning plot (full range)
# full lambda grid range with warm.start = 0.01
# here preset lambda coincides with best lambda
tune_full <- L0glm.tune(X=sim.gaus$X, 
                   y=sim.gaus$y, 
                   IC = "gic", 
                   nlambda = 10,
                   prop_range = 1, # full range up to lambda_max
                   warm.start = 0, # 0.01, # choose 0 for no warm starts
                   plot = TRUE, 
                   coefs.true = sim.gaus$beta_true, 
                   relaxed = FALSE, 
                   family = gaussian(identity), 
                   nonnegative = TRUE, 
                   solver = "eigen", 
                   method = 7L)
tune_full$lambda_best # 20.44997 (best lambda tune)
tune_full$lambda_preset # 20.44997 (lambda preset)
min(tune_full$metrics$gic) # 82400.61

## tuning plot (zoomed in)
# tuning plots - with lambda grid more zoomed in around lambda_preset preset
tune_zoomed <- L0glm.tune(X=sim.gaus$X, 
                   y=sim.gaus$y, 
                   IC = "gic", 
                   nlambda = 10,
                   prop_range = 0.2, # to zoom in around lambda_preset preset
                   warm.start = 0, # 0.01, 
                   plot = TRUE, 
                   coefs.true = sim.gaus$beta_true, 
                   relaxed = FALSE, 
                   family = gaussian(identity), 
                   nonnegative = TRUE, 
                   solver = "eigen", 
                   method = 7L)
# here best lambda (dashed line) is very close to lambda preset expected 
# to optimise GIC (dotted line) so tuning in fact not needed
tune_zoomed$lambda_best # 60.73341 (best lambda tune)
tune_zoomed$lambda_preset # 20.44997 (lambda preset)
min(tune_zoomed$metrics$gic) # 82203.1

# TODO: check why coefs are going to zero before calculated lambda_max

system.time(tune.chol.sc.gaus <- L0glm.tune(X=sim.gaus$X, y=sim.gaus$y, IC = "gic", plot = TRUE, coefs.true = sim.gaus$beta_true, relaxed = FALSE, family = gaussian(identity), nonnegative = TRUE, solver = "eigen", method = 7L))
# Best lambda by GIC = 20.4, time = 3.00s with n=p=10000
system.time(tune.CG.sc.gaus <- L0glm.tune(X=sim.gaus$X, y=sim.gaus$y, IC = "gic", plot = FALSE, coefs.true = sim.gaus$beta_true, relaxed = FALSE, family = gaussian(identity), nonnegative = TRUE, solver = "eigen", method = 5L))
# Best lambda by GIC = 20.4, time = 18.17s with n=p=10000, warning eigen CG solver failed in IRLS iteration 3, trying Cholesky LLT solver as fallback.
system.time(tune.LSCG.sc.gaus <- L0glm.tune(X=sim.gaus$X, y=sim.gaus$y, IC = "gic", plot = TRUE, coefs.true = sim.gaus$beta_true, relaxed = FALSE, family = gaussian(identity), nonnegative = TRUE, solver = "eigen", method = 3L))
# Best lambda by GIC = 20.4, time = 39.12s with n=p=10000, warning eigen LSCG solver failed in IRLS iteration 3 & 4, trying Cholesky LLT solver as fallback
system.time(tune.osqp.sc.gaus <- L0glm.tune(X=sim.gaus$X, y=sim.gaus$y, IC = "gic", plot = FALSE, coefs.true = sim.gaus$beta_true, relaxed = FALSE, family = gaussian(identity), nonnegative = TRUE, solver = "osqp"))
# Best lambda by GIC = 20.4, time = 168.51s with n=p=10000
system.time(tune.glmnet.sc.gaus <- L0glm.tune(X=sim.gaus$X, y=sim.gaus$y, IC = "gic", plot = FALSE, coefs.true = sim.gaus$beta_true, relaxed = FALSE, family = gaussian(identity), nonnegative = TRUE, solver = "glmnet"))
# Best lambda by GIC = 20.4, time = 44.43s with n=p=10000
system.time(tune.QR.sc.gaus <- L0glm.tune(X=sim.gaus$X, y=sim.gaus$y, IC = "gic", plot = FALSE, coefs.true = sim.gaus$beta_true, relaxed = TRUE, family = gaussian(identity), nonnegative = TRUE, solver = "eigen", method = 1L))
# Best lambda by GIC = 20.4, time = 167.56s with n=p=10000

# Chosen lambda here matches expected lambda value that
# was expected to optimise GIC, so tuning in fact not needed
L0glm:::preset.lambda(IC="gic", p, n) # 20.44997 (lambda preset)
tune.chol.sc.gaus$lambda_preset # 20.44997 (lambda preset)
tune.chol.sc.gaus$lambda_best # 20.44997 (best lambda tune)


#### 1.2.1 L0GLM FIT WITH Cholesky SOLVER ####
system.time(L0glm.chol.gaus <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                         family = gaussian(identity), 
                                         lambda = "gic", # or tune_zoomed$lambda_best
                                         nonnegative = TRUE,
                                         relaxed = FALSE,
                                         solver = "eigen",
                                         method = 7L)) 
# 0.16s for n=p=10000

#### 1.2.2 L0GLM FIT WITH CG SOLVER ####
system.time(L0glm.CG.gaus <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                       family = gaussian(identity), 
                                       lambda = "gic", # or tune_zoomed$lambda_best
                                       nonnegative = TRUE,
                                       relaxed = TRUE,
                                       solver = "eigen",
                                       method = 5L)) 
# 3.38s for n=p=10000


#### 1.2.3 L0GLM FIT WITH LSCG SOLVER ####
system.time(L0glm.LSCG.gaus <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                         family = gaussian(identity), 
                                         lambda = "gic", # or tune_zoomed$lambda_best
                                         nonnegative = TRUE,
                                         relaxed = FALSE,
                                         solver = "eigen",
                                         method = 4L)) 
# 7s for n=p=10000


#### 1.2.4 L0GLM FIT WITH OSQP SOLVER ####
system.time(L0glm.osqp.gaus <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                    family = gaussian(identity), 
                                    lambda = "gic", # or tune_zoomed$lambda_best
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "osqp")) 
# 37s for n=p=10000


#### 1.2.5 L0GLM FIT WITH glmnet::bigGlm SOLVER ####
system.time(L0glm.glmnet.gaus <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                      family = gaussian(identity), 
                                      lambda = "gic", # or tune_zoomed$lambda_best
                                      nonnegative = TRUE,
                                      relaxed = FALSE,
                                      solver = "glmnet")) 
# 6.5s for n=p=10000 

#### 1.2.6 L0GLM FIT WITH QR SOLVER ####
system.time(L0glm.QR.gaus <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                  family = gaussian(identity), 
                                  lambda = "gic", # or tune_zoomed$lambda_best
                                  nonnegative = TRUE,
                                  relaxed = FALSE,
                                  solver = "eigen",
                                  method = 1L)) 
# 4.85s for n=p=10000


#### 1.2.7 L0GLM FIT WITH NNLS SOLVER (slow for large problem sizes) ####
# system.time(L0glm.nnls.gaus <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
#                                     family = gaussian(identity), 
#                                     lambda = "gic", # or tune_zoomed$lambda_best
#                                     nonnegative = TRUE,
#                                     relaxed = FALSE,
#                                     solver = "nnls")) 
# XXs with n=p=10000

 
# #### 1.2.8 PLAIN NNLS FIT WITHOUT REGULARISATION ####
# system.time(L0glm.nnlsnoreg.gaus <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
#                                          family = gaussian(identity), 
#                                          lambda = 0,
#                                          nonnegative = TRUE,
#                                          relaxed = FALSE,
#                                          solver = "nnls",
#                                          maxit = 1)) 
# XXs with n=p=10000



### 1.2.9 BENCHMARKS TABLE ####

X = sim.gaus$X
y = sim.gaus$y
benchmarks_1_2 = rbind(
  data.frame(method="L0glm_osqp", do.call(cbind, L0glm:::get.metrics(beta = L0glm.osqp.gaus$coefficients,  # some metrics
                                                             X = X, y = y, 
                                                             family = gaussian(identity), 
                                                             beta_true = sim.gaus$beta_true))),
  data.frame(method="L0glm_bigglm", do.call(cbind, L0glm:::get.metrics(beta = L0glm.glmnet.gaus$coefficients,  # some metrics
                                                               X = X, y = y, 
                                                               family = gaussian(identity), 
                                                               beta_true = sim.gaus$beta_true))),
  data.frame(method="L0glm_qr", do.call(cbind, L0glm:::get.metrics(beta = L0glm.QR.gaus$coefficients,  # some metrics
                                                           X = X, y = y, 
                                                           family = gaussian(identity),
                                                           beta_true = as.vector(sim.gaus$beta_true)) )),
  data.frame(method="L0glm_cholesky", do.call(cbind, L0glm:::get.metrics(beta = L0glm.chol.gaus$coefficients,  # some metrics
                                                                 X = X, y = y, 
                                                                 family = gaussian(identity),
                                                                 beta_true = sim.gaus$beta_true) )),
  data.frame(method="L0glm_LSCG", do.call(cbind, L0glm:::get.metrics(beta = L0glm.LSCG.gaus$coefficients,  # some metrics
                                                             X = X, y = y, 
                                                             family = gaussian(identity),
                                                             beta_true = sim.gaus$beta_true) )),
  data.frame(method="L0glm_CG", do.call(cbind, L0glm:::get.metrics(beta = L0glm.CG.gaus$coefficients,  # some metrics
                                                           X = X, y = y, 
                                                           family = gaussian(identity), 
                                                           beta_true = sim.gaus$beta_true) )) 
)

benchmarks_1_2 = benchmarks_1_2[order(benchmarks_1_2[,"deviance_true"]),] # sort by deviance relative to ground truth
benchmarks_1_2
#           method     p     n m deviance deviance_true    loglik
# 1     L0glm_osqp 10001 10000 1 206186.4      36010.75 -29321.36
# 4 L0glm_cholesky 10001 10000 1 254232.4      87854.15 -30368.70
# 3       L0glm_qr 10001 10000 1 254232.4      87854.16 -30368.70
# 5     L0glm_LSCG 10001 10000 1 267722.9     100738.49 -30627.22
# 6       L0glm_CG 10001 10000 1 273371.6     105528.80 -30731.62
# 2   L0glm_bigglm 10001 10000 1 697016.3     541087.64 -35411.50
#        MD_coefs MAD_coefs sMAPE_coefs bias_coefs  R2_coefs MD_selcoefs
# 1 -0.0062328918  16.00636  0.02378202   1.000924 0.9898070    33.75786
# 4 -0.0016361433  22.74395  0.03754825   1.016802 0.9860743    61.49606
# 3 -0.0016361433  22.74395  0.03754825   1.016802 0.9860743    61.49606
# 5 -0.0005941144  27.75970  0.03825287   1.026381 0.9659579    63.88428
# 6  0.0008094115  28.16403  0.03838781   1.026418 0.9657681    65.39540
# 2  0.0253012046  99.17600  0.09883913   1.040941 0.8931648   105.35647
#   MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs
# 1     120.5346      0.1579751      1.000924   0.9914452
# 4     153.4596      0.2414769      1.016802   0.9895949
# 3     153.4596      0.2414769      1.016802   0.9895949
# 5     198.4891      0.2462464      1.026381   0.9697854
# 6     201.0648      0.2458639      1.026418   0.9696615
# 2     709.7733      0.6264514      1.040941   0.8997974
#   support_size_true support_size_estimated
# 1              1000                   1037
# 4              1000                   1058
# 3              1000                   1058
# 5              1000                   1058
# 6              1000                   1057
# 2              1000                   1217
#   estimated_over_true_support_size  TP  FP   TN  FN precision recall
# 1                            1.037 963  74 8926  37 0.9286403  0.963
# 4                            1.058 940 118 8882  60 0.8884688  0.940
# 3                            1.058 940 118 8882  60 0.8884688  0.940
# 5                            1.058 939 119 8881  61 0.8875236  0.939
# 6                            1.057 938 119 8881  62 0.8874172  0.938
# 2                            1.217 887 330 8670 113 0.7288414  0.887
#   specificity accuracy  F1_score      gic     ebic      bic      aic
# 1   0.9917778   0.9889 0.9455081 79870.02 77763.60 68203.06 60718.73
# 4   0.9868889   0.9822 0.9135083 82394.16 80245.12 70491.16 62855.41
# 3   0.9868889   0.9822 0.9135083 82394.16 80245.12 70491.16 62855.41
# 5   0.9867778   0.9820 0.9125364 82911.19 80762.16 71008.19 63372.44
# 6   0.9867778   0.9819 0.9120078 83099.54 80952.53 71207.78 63579.24
# 2   0.9633333   0.9557 0.8001804 95731.34 93259.64 82041.20 73259.01
#        mbic
# 1  84445.99
# 4  87062.70
# 3  87062.70
# 5  87579.74
# 6  87763.67
# 2 101100.82

# -> here osqp solver gave best solution & glmnet's bigGlm the worst
# Cholesky comes in at 2nd place though & was by far the fastest
# -> Cholesky performs best for different noise families & is fastest


# 2. MULTICHANNEL ####

# TODO: add option precluster = T or F to L0glm.fit?
# (see ?glm_cluster)

# Simulate MC n=p=500
n <- p <- 500
npeaks <- 50
sim.gaus.mc <- simul(method = 3, n = n, p = p, s = npeaks/p, seed = 1)
sim.pois.mc <- simul(method = 4, n = n, p = p, s = npeaks/p, seed = 2)
dim(sim.pois.mc$y) # 500 timepoints x 168 m/z channels
m <- ncol(sim.pois.mc$y)


## 2.1. MULTICHANNEL WITH POISSON NOISE ####

### 2.1.1 GROUPED L0GLM FIT WITH relaxed = FALSE, grouped = FALSE (not great) ####

## tuning plot (full range)
# here preset lambda would be ca 10x higher than best lambda
tune_full <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, 
           IC = "gic", 
           nlambda = 10,
           prop_range = 1, # to zoom in around lambda_preset preset
           warm.start = 0.01, #  set to 0 for no warm starts to not affect location of optimal lambda
           plot = TRUE, coefs.true = sim.pois.mc$beta_true, 
           relaxed = FALSE, family = poisson(identity), 
           nonnegative = TRUE, 
           parallelize = FALSE, 
           grouped = FALSE,
           solver = "eigen", method = 7L)
tune_full$lambda_best # 0.9753673 (best lambda tune)
tune_full$lambda_preset # 11.35348 (lambda preset)
L0glm:::preset.lambda(IC="gic", p, n) # 11.35348 (lambda preset)

## tuning plot (zoomed in)
# here preset lambda would be ca 5x higher than best lambda
tune_zoomed <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, 
           IC = "gic", 
           nlambda = 10,
           prop_range = 0.2, # to zoom in around lambda_preset preset
           warm.start = 0.01, # set to 0 for no warm starts to not affect location of optimal lambda
           plot = TRUE, coefs.true = sim.pois.mc$beta_true, 
           relaxed = FALSE, family = poisson(identity), 
           nonnegative = TRUE, 
           parallelize = FALSE, 
           grouped = FALSE,
           solver = "eigen", method = 7L)
tune_zoomed$lambda_best # 1.593531 (best lambda tune)
tune_zoomed$lambda_preset # 11.35348 (lambda preset)
L0glm:::preset.lambda(IC="gic", p, n) # 11.35348 (lambda preset)


#### Cholesky ####
system.time(L0glm.chol <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                    family = poisson(identity), 
                                    lambda = "gic", # or tune_zoomed$lambda_best
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallelize = FALSE,
                                    warm.start = 1,
                                    grouped = FALSE,
                                    solver = "eigen",
                                    method = 7L)) 
# 4.01s with parallelize=F, 5.63s with parallelize=T

x_vals <- as.vector(rowMeans(L0glm.chol$coefficients))        # estimated
y_vals <- as.vector(rowMeans(sim.pois.mc$beta_true))               # true
fit <- lm(y_vals ~ 0 + x_vals)
R2p <- round(100 * summary(fit)$r.squared, 2)
main_expr <- bquote(
  atop("true vs estimated coefficients (row means)",
       atop("nonnegative independent L0 Poisson GLMs",
            R^2 == .(sprintf("%.2f", R2p)) * "%")))
plot(x_vals, y_vals, pch = 16, col = "steelblue",
     xlab = "estimated coefficients", 
     ylab = "true coefficients",
     main = main_expr)



#### OSQP ####
system.time(L0glm.osqp <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                    family = poisson(identity), 
                                    lambda = "gic", # or tune_zoomed$lambda_best
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallelize = FALSE,
                                    grouped = FALSE,
                                    solver = "osqp")) 
# time = 8.06s with parallelize=F, 6.3s with parallelize=T





### 2.1.2 GROUPED L0GLM FIT WITH relaxed = FALSE, grouped = TRUE (better) ####

## tuning plot (full range)
# here preset lambda would be ca 10x higher than best lambda
tune_full <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, 
                        IC = "gic", 
                        nlambda = 10,
                        prop_range = 1, # to zoom in around lambda_preset preset
                        warm.start = 0.01, # here set to 0 for no warm starts to not affect location of optimal lambda
                        plot = TRUE, coefs.true = sim.pois.mc$beta_true, 
                        relaxed = FALSE, family = poisson(identity), 
                        nonnegative = TRUE, 
                        parallelize = FALSE, 
                        grouped = TRUE,
                        solver = "eigen", method = 7L)
tune_full$lambda_best # 0.9753657 (best lambda tune)
tune_full$lambda_preset # 11.35348 (lambda preset)
L0glm:::preset.lambda(IC="gic", p, n) # 11.35348 (lambda preset)

## tuning plot (zoomed in)
# here preset lambda would be ca 5x higher than best lambda
tune_zoomed <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, 
                          IC = "gic", 
                          nlambda = 10,
                          prop_range = 0.2, # to zoom in around lambda_preset preset
                          warm.start = 0.01, # set to 0 for no warm starts to not affect location of optimal lambda
                          plot = TRUE, coefs.true = sim.pois.mc$beta_true, 
                          relaxed = FALSE, family = poisson(identity), 
                          nonnegative = TRUE, 
                          parallelize = FALSE, 
                          grouped = TRUE,
                          solver = "eigen", method = 7L)
tune_zoomed$lambda_best # 2.60347 (best lambda tune)
tune_zoomed$lambda_preset # 11.35348 (lambda preset)
L0glm:::preset.lambda(IC="gic", p, n) # 11.35348 (lambda preset)



#### Cholesky ####
system.time(L0glm.chol <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                    family = poisson(identity), 
                                    lambda = "gic", # or tune_zoomed$lambda_best
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallelize = FALSE,
                                    grouped = TRUE,
                                    solver = "eigen",
                                    method = 7L)) 
# 18s


#### OSQP ####
system.time(L0glm.osqp <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                    family = poisson(identity), 
                                    lambda = "gic", # or tune_zoomed$lambda_best
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallelize = FALSE,
                                    grouped = TRUE,
                                    solver = "osqp"))
# 32s




## 2.1.3 GROUPED L0GLM FIT WITH relaxed = TRUE, grouped = TRUE (best) ####

## tuning plot (full range)
# here preset lambda would be ca 10x higher than best lambda
tune_full <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, 
                        IC = "gic", 
                        nlambda = 10,
                        prop_range = 1, # to zoom in around lambda_preset preset
                        warm.start = 0.01, # here set to 0 for no warm starts to not affect location of optimal lambda
                        plot = TRUE, coefs.true = sim.pois.mc$beta_true, 
                        relaxed = TRUE, family = poisson(identity), 
                        nonnegative = TRUE, 
                        parallelize = FALSE, 
                        grouped = TRUE,
                        solver = "eigen", method = 7L)
tune_full$lambda_best # 0.9753657 (best lambda tune)
tune_full$lambda_preset # 11.35348 (lambda preset)
L0glm:::preset.lambda(IC="gic", p, n) # 11.35348 (lambda preset)

## tuning plot (zoomed in)
# here preset lambda would be ca 10x higher than best lambda
tune_zoomed <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, 
                          IC = "gic", 
                          nlambda = 10,
                          prop_range = 0.2, # to zoom in around lambda_preset preset
                          warm.start = 0.01, # set to 0 for no warm starts to not affect location of optimal lambda
                          plot = TRUE, coefs.true = sim.pois.mc$beta_true, 
                          relaxed = TRUE, family = poisson(identity), 
                          nonnegative = TRUE, 
                          parallelize = FALSE, 
                          grouped = TRUE,
                          solver = "eigen", method = 7L)
tune_zoomed$lambda_best # 0.9755287 (best lambda tune)
tune_zoomed$lambda_preset # 11.35348 (lambda preset)
L0glm:::preset.lambda(IC="gic", p, n) # 11.35348 (lambda preset)


#### Cholesky ####
system.time(L0glm.chol <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                    family = poisson(identity), 
                                    lambda = "gic", # or tune_zoomed$lambda_best
                                    nonnegative = TRUE,
                                    relaxed = TRUE,
                                    parallelize = FALSE,
                                    grouped = TRUE,
                                    solver = "eigen",
                                    method = 7L)) 
# 20s

x_vals <- as.vector(rowMeans(L0glm.chol$coefficients))        # estimated
y_vals <- as.vector(rowMeans(sim.pois.mc$beta_true))               # true
fit <- lm(y_vals ~ 0 + x_vals)
R2p <- round(100 * summary(fit)$r.squared, 2)
main_expr <- bquote(
  atop("true vs estimated coefficients  (row means)",
       atop("relaxed nonnegative grouped L0 Poisson GLM",
            R^2 == .(sprintf("%.2f", R2p)) * "%")))
plot(x_vals, y_vals, pch = 16, col = "steelblue",
     xlab = "estimated coefficients", 
     ylab = "true coefficients",
     main = main_expr)



## Fit relaxed = TRUE, grouped = TRUE
#### OSQP ####
system.time(L0glm.osqp <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                    family = poisson(identity), 
                                    lambda = "gic", # or tune_zoomed$lambda_best
                                    nonnegative = TRUE,
                                    relaxed = TRUE,
                                    parallelize = FALSE,
                                    grouped = TRUE,
                                    solver = "osqp")) 
# 37s




## 2.2. MULTICHANNEL WITH GAUSSIAN NOISE ####

### 2.2.1 GROUPED L0GLM FIT WITH relaxed = FALSE, grouped = FALSE (not great) ####

## tuning plot (full range)
# here best lambda is identical to preset lambda
tune_full <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, 
                        IC = "gic", 
                        nlambda = 10,
                        prop_range = 1, # to zoom in around lambda_preset preset
                        warm.start = 0.01, # here set to 0 for no warm starts to not affect location of optimal lambda
                        plot = TRUE, 
                        coefs.true = sim.gaus.mc$beta_true, 
                        relaxed = FALSE, 
                        family = gaussian(identity), 
                        nonnegative = TRUE, 
                        parallelize = FALSE, 
                        grouped = FALSE,
                        solver = "eigen", method = 7L)
tune_full$lambda_best # 11.35348
tune_full$lambda_preset # 11.35348 
L0glm:::preset.lambda(IC="gic", p, n) # 11.35348

## tuning plot (zoomed in)
# here best lambda is identical to preset lambda
tune_zoomed <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, 
                          IC = "gic", 
                          nlambda = 10,
                          prop_range = 0.2, # to zoom in around lambda_preset preset
                          warm.start = 0, # here set to 0 for no warm starts to not affect location of optimal lambda
                          plot = TRUE, 
                          coefs.true = sim.gaus.mc$beta_true, 
                          relaxed = FALSE, 
                          family = gaussian(identity), 
                          nonnegative = TRUE, 
                          parallelize = FALSE, 
                          grouped = FALSE,
                          solver = "eigen", method = 7L)
tune_zoomed$lambda_best # 11.35348
tune_zoomed$lambda_preset # 11.35348 
L0glm:::preset.lambda(IC="gic", p, n) # 11.35348

#### Cholesky ####
system.time(L0glm.chol <- L0glm.fit(sim.gaus.mc$X, sim.gaus.mc$y, 
                                    family = gaussian(identity), 
                                    lambda = "gic", # or tune_zoomed$lambda_best
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallelize = FALSE,
                                    grouped = FALSE,
                                    solver = "eigen",
                                    method = 7L)) 
# 4s

x_vals <- as.vector(rowMeans(L0glm.chol$coefficients))        # estimated
y_vals <- as.vector(rowMeans(sim.gaus.mc$beta_true))               # true
fit <- lm(y_vals ~ 0 + x_vals)
R2p <- round(100 * summary(fit)$r.squared, 2)
main_expr <- bquote(
  atop("true vs estimated coefficients  (row means)",
       atop("nonnegative L0 independent Gaussian GLMs",
            R^2 == .(sprintf("%.2f", R2p)) * "%")))
plot(x_vals, y_vals, pch = 16, col = "steelblue",
     xlab = "estimated coefficients", 
     ylab = "true coefficients",
     main = main_expr)




#### OSQP ####
system.time(L0glm.osqp <- L0glm.fit(sim.gaus.mc$X, sim.gaus.mc$y, 
                                    family = gaussian(identity), 
                                    lambda = "gic", # or tune_zoomed$lambda_best
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallelize = FALSE,
                                    grouped = FALSE,
                                    solver = "osqp")) 
# 16s




### 2.2.2 GROUPED L0GLM FIT WITH relaxed = FALSE, grouped = TRUE (better) ####

## tuning plot (full range)
# here best lambda is identical to preset lambda
tune_full <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, 
                        IC = "gic", 
                        nlambda = 10,
                        prop_range = 1, # to zoom in around lambda_preset preset
                        warm.start = 0.01, # here set to 0 for no warm starts to not affect location of optimal lambda
                        plot = TRUE, 
                        coefs.true = sim.gaus.mc$beta_true, 
                        relaxed = FALSE, 
                        family = gaussian(identity), 
                        nonnegative = TRUE, 
                        parallelize = FALSE, 
                        grouped = TRUE,
                        solver = "eigen", method = 7L)
tune_full$lambda_best # 11.35348
tune_full$lambda_preset # 11.35348 
L0glm:::preset.lambda(IC="gic", p, n) # 11.35348

## tuning plot (zoomed in)
# here best lambda is identical to preset lambda
tune_zoomed <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, 
                          IC = "gic", 
                          nlambda = 10,
                          prop_range = 0.2, # to zoom in around lambda_preset preset
                          warm.start = 0, # here set to 0 for no warm starts to not affect location of optimal lambda
                          plot = TRUE, 
                          coefs.true = sim.gaus.mc$beta_true, 
                          relaxed = FALSE, 
                          family = gaussian(identity), 
                          nonnegative = TRUE, 
                          parallelize = FALSE, 
                          grouped = TRUE,
                          solver = "eigen", method = 7L)
tune_zoomed$lambda_best # 11.35348
tune_zoomed$lambda_preset # 11.35348 
L0glm:::preset.lambda(IC="gic", p, n) # 11.35348

#### Cholesky ####
system.time(L0glm.chol <- L0glm.fit(sim.gaus.mc$X, sim.gaus.mc$y, 
                                    family = gaussian(identity), 
                                    lambda = "gic", # or tune.chol.mc$lambda_best,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallelize = FALSE,
                                    grouped = TRUE,
                                    solver = "eigen",
                                    method = 7L)) 
# 19s for n=p=500

x_vals <- as.vector(rowMeans(L0glm.chol$coefficients))        # estimated
y_vals <- as.vector(rowMeans(sim.gaus.mc$beta_true))               # true
fit <- lm(y_vals ~ 0 + x_vals)
R2p <- round(100 * summary(fit)$r.squared, 2)
main_expr <- bquote(
  atop("true vs estimated coefficients (row means)",
       atop("nonnegative grouped L0 Gaussian GLM",
            R^2 == .(sprintf("%.2f", R2p)) * "%")))
plot(x_vals, y_vals, pch = 16, col = "steelblue",
     xlab = "estimated coefficients", 
     ylab = "true coefficients",
     main = main_expr)


#### OSQP ####
system.time(L0glm.osqp <- L0glm.fit(sim.gaus.mc$X, sim.gaus.mc$y, 
                                    family = gaussian(identity), 
                                    lambda = "gic", # or tune.osqp.mc$lambda_best,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallelize = FALSE,
                                    grouped = TRUE,
                                    solver = "osqp")) 
# 43s for n=p=500



### 2.2.3 GROUPED L0GLM FIT WITH relaxed = TRUE, grouped = TRUE (best) ####

## tuning plot (full range)
# here best lambda is identical to preset lambda
tune_full <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, 
                        IC = "gic", 
                        nlambda = 10,
                        prop_range = 1, # to zoom in around lambda_preset preset
                        warm.start = 0.01, # here set to 0 for no warm starts to not affect location of optimal lambda
                        plot = TRUE, 
                        coefs.true = sim.gaus.mc$beta_true, 
                        relaxed = TRUE, 
                        family = gaussian(identity), 
                        nonnegative = TRUE, 
                        parallelize = FALSE, 
                        grouped = TRUE,
                        solver = "eigen", method = 7L)
tune_full$lambda_best # 11.35348 (best lambda tune)
tune_full$lambda_preset # 11.35348 (lambda preset)
L0glm:::preset.lambda(IC="gic", p, n) # 11.35348 (lambda preset)

## tuning plot (zoomed in)
# here best lambda is again very close to preset lambda
tune_zoomed <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, 
                          IC = "gic", 
                          nlambda = 10,
                          prop_range = 0.2, # to zoom in around lambda_preset preset
                          warm.start = 0.01, # set to 0 for no warm starts to not affect location of optimal lambda
                          plot = TRUE, 
                          coefs.true = sim.gaus.mc$beta_true, 
                          relaxed = TRUE, 
                          family = gaussian(identity), 
                          nonnegative = TRUE, 
                          parallelize = FALSE, 
                          grouped = TRUE,
                          solver = "eigen", method = 7L)
tune_zoomed$lambda_best # 4.82 (best lambda tune)
tune_zoomed$lambda_preset # 11.35348 (lambda preset)
L0glm:::preset.lambda(IC="gic", p, n) # 11.35348 (lambda preset)


#### Cholesky ####
system.time(L0glm.chol <- L0glm.fit(sim.gaus.mc$X, sim.gaus.mc$y, 
                                    family = gaussian(identity), 
                                    lambda = "gic", # or tune_zoomed$lambda_best
                                    nonnegative = TRUE,
                                    relaxed = TRUE,
                                    parallelize = FALSE,
                                    grouped = TRUE,
                                    solver = "eigen",
                                    method = 7L)) 
# 19s

x_vals <- as.vector(rowMeans(L0glm.chol$coefficients))        # estimated
y_vals <- as.vector(rowMeans(sim.gaus.mc$beta_true))               # true
fit <- lm(y_vals ~ 0 + x_vals)
R2p <- round(100 * summary(fit)$r.squared, 2)
main_expr <- bquote(
  atop("true vs estimated coefficients (row means)",
       atop("relaxed nonnegative grouped L0 Gaussian GLM",
            R^2 == .(sprintf("%.2f", R2p)) * "%")))
plot(x_vals, y_vals, pch = 16, col = "steelblue",
     xlab = "estimated coefficients", 
     ylab = "true coefficients",
     main = main_expr)



#### OSQP ####
system.time(L0glm.osqp <- L0glm.fit(sim.gaus.mc$X, sim.gaus.mc$y, 
                                    family = gaussian(identity), 
                                    lambda = "gic", # or tune_zoomed$lambda_best
                                    nonnegative = TRUE,
                                    relaxed = TRUE,
                                    parallelize = FALSE,
                                    grouped = TRUE,
                                    solver = "osqp")) 
# 23s for n=p=500






# 3. BENCHMARKS AGAINST OTHER METHODS ####
# Take our best working solver & method and compare it against others: 
# glmnet (LASSO), ordinis (MCP & SCAD), L0Learn & abess

## 3.1 SINGLE CHANNEL POISSON ####
n <- p <- 500
npeaks <- 50
sim.pois <- simul(method = 2, n = n, p = p, s = npeaks/p, seed = 2)

### 3.1.1 L0GLM VIA CHOLESKY ####
system.time(tune.chol.sc <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "gic", 
                                       plot = FALSE, coefs.true = sim.pois$beta_true, 
                                       relaxed = FALSE, family = poisson(identity), 
                                       nonnegative = TRUE, 
                                       solver = "eigen", method = 7L))
# 0.5s

system.time(L0glm.chol <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                    family = poisson(identity), 
                                    lambda = "gic", # or tune.chol.sc$lambda_best,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "eigen",
                                    method = 7L)) 
# 0.04s
plot(x=L0glm.chol$coefficients, y=sim.pois$beta_true, pch=16, col='blue',
     main="L0 glm identity link Poisson with lambda tuned to optimise GIC",
     xlab="estimated coefficients", ylab="true coefficients")


### 3.1.2 LASSO USING GLMNET ####
system.time(glmnet.sc.pois <- fit.GLM(sim = sim.pois, method = 2, 
                                      family = poisson(identity), plot = TRUE)) 
# 0.18s


### 3.1.3 ADAPTIVE LASSO USING GLMNET ####
system.time(glmnet.al.sc.pois <- fit.GLM(sim = sim.pois, method = 3, 
                                         family = poisson(identity), plot = T)) 
# 0.25s

### 3.1.4 L0LEARN ####
system.time(L0learn.sc.pois <- fit.GLM(sim = sim.pois, method = 4, 
                                       family = poisson(identity), plot = T))  
# 0.19s


### 3.1.5 L0LEARN AFTER LASSO ####
# better solution than regular L0Learn
system.time(L0learn.after.lasso.sc.pois <- fit.GLM(sim = sim.pois, method = 5, family = poisson(identity), plot = T)) 
# 0.36s


### 3.1.6 ABESS ####
system.time(abess.sc.pois <- fit.GLM(sim = sim.pois, method = 6, 
                                     family = poisson(identity), 
                                     plot = T, IC = "gic"))  
# 1.2s
# TODO: NOTE: result is better if I do not use weights 
# CHECK WHY

### 3.1.7 ABESS AFTER LASSO ####
system.time(abess.after.lasso.sc.pois <- fit.GLM(sim = sim.pois, method = 7, 
                                                 family = poisson(identity), 
                                                 plot = T, 
                                                 IC = "gic"))  
# 0.19s
# TODO: NOTE: result is better if I do not use weights - then cor 0.8849495 - CHECK
# CHECK WHY


### 3.1.8 MCP USING ORDINIS ####
system.time(ordinis.mcp.sc.pois <- fit.GLM(sim = sim.pois, method = 8, 
                                           family = poisson(identity), plot = T,
                                           IC = "gic"))  
# 10s

### 3.1.9 SCAD USING ORDINIS ####
system.time(ordinis.scad.sc.pois <- fit.GLM(sim = sim.pois, method = 9, 
                                            family = poisson(identity), 
                                            plot = T, IC = "gic"))  
# 10.5s




### 3.1.10 BENCHMARKS TABLE ####

X <- sim.pois$X
y <- sim.pois$y
# W <- 1/(y+0.1)
IC <- "gic"

benchmarks_3_1 = rbind(
  data.frame(method="L0glm", do.call(cbind, L0glm:::get.metrics(beta = L0glm.chol$coefficients,  # some metrics
                                                        X = X, y = y, 
                                                        family = poisson(identity), 
                                                        beta_true = sim.pois$beta_true))),
  data.frame(method="nnalasso_glmnet", do.call(cbind, L0glm:::get.metrics(beta = glmnet.al.sc.pois,  # some metrics
                                                                  X = cbind(1, X), y = y, 
                                                                  family = poisson(identity), 
                                                                  beta_true = c(0, as.vector(sim.pois$beta_true))) )),
  data.frame(method="nnL0Learn_after_glmnet_nnlasso", do.call(cbind, L0glm:::get.metrics(beta = L0learn.after.lasso.sc.pois,  # some metrics
                                                                                 X = X, y = y, 
                                                                                 family = poisson(identity),
                                                                                 beta_true = as.vector(sim.pois$beta_true)) )),
  data.frame(method="nnabess", do.call(cbind, L0glm:::get.metrics(beta = abess.sc.pois,  # some metrics
                                                          X = X, y = y, 
                                                          family = poisson(identity),
                                                          beta_true = as.vector(sim.pois$beta_true)) )),
  data.frame(method="nnabess_after_glmnet_nnlasso", do.call(cbind, L0glm:::get.metrics(beta = abess.after.lasso.sc.pois,  # some metrics
                                                                               X = X, y = y, 
                                                                               family = poisson(identity),
                                                                               beta_true = as.vector(sim.pois$beta_true)) )),
  data.frame(method="nnL0Learn", do.call(cbind, L0glm:::get.metrics(beta = L0learn.sc.pois,  # some metrics
                                                            X = X, y = y, 
                                                            family = poisson(identity), 
                                                            beta_true = as.vector(sim.pois$beta_true)) )),
  data.frame(method="nnlasso_glmnet", do.call(cbind, L0glm:::get.metrics(beta = glmnet.sc.pois,  # some metrics
                                                                 X = cbind(1, X), y = y, 
                                                                 family = poisson(identity), 
                                                                 beta_true = c(0, as.vector(sim.pois$beta_true))) )),
  data.frame(method="nnMCP_ordinis", do.call(cbind, L0glm:::get.metrics(beta = ordinis.mcp.sc.pois,  # some metrics
                                                                X = cbind(1, X), y = y, 
                                                                family = poisson(identity),
                                                                beta_true = c(0, as.vector(sim.pois$beta_true))) )),
  data.frame(method="nnSCAD_ordinis", do.call(cbind, L0glm:::get.metrics(beta = ordinis.scad.sc.pois,  # some metrics
                                                                 X = cbind(1, X), y = y, 
                                                                 family = poisson(identity), 
                                                                 beta_true = c(0, as.vector(sim.pois$beta_true))) )))

benchmarks_3_1 = benchmarks_3_1[order(benchmarks_3_1[,"deviance_true"]),] # sort by deviance relative to ground truth
benchmarks_3_1
#                           method   p   n m     deviance deviance_true
# 1                          L0glm 500 500 1 4.407922e+02  5.888956e+01
# 2                nnalasso_glmnet 501 500 1 8.910231e+02  5.260599e+02
# 7                 nnlasso_glmnet 501 500 1 9.088932e+02  5.869995e+02
# 9                 nnSCAD_ordinis 501 500 1 1.387974e+03  1.172627e+03
# 8                  nnMCP_ordinis 501 500 1 1.419627e+03  1.204754e+03
# 3 nnL0Learn_after_glmnet_nnlasso 500 500 1 1.417924e+07  1.414749e+07
# 6                      nnL0Learn 500 500 1 1.417924e+07  1.414749e+07
# 5   nnabess_after_glmnet_nnlasso 500 500 1 1.038960e+08  1.038980e+08
# 4                        nnabess 500 500 1 1.052046e+08  1.052121e+08
#          loglik      MD_coefs   MAD_coefs sMAPE_coefs  bias_coefs
# 1     -2139.572     0.4515654    54.97077  0.05471239 0.959293624
# 2     -2364.688    -1.2191931    73.01081  0.09703150 1.010322884
# 7     -2373.623    -6.4230012   168.57499  0.36660963 1.156047013
# 9     -2613.163    -0.2072492   322.92210  0.40784001 1.195603989
# 8     -2628.990    -0.2296764   325.44781  0.40850065 1.194418337
# 3  -7091540.369   -87.2668739   204.21553  0.16016479 0.794622032
# 6  -7091540.369   -87.2668739   204.21553  0.16016479 0.794622032
# 5 -51949911.351 17789.8020892 18214.16146  0.25006342 0.004393986
# 4 -52604208.929 17978.4522010 18430.23369  0.25890768 0.004496484
#    R2_coefs   MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs
# 1 0.9385868    259.072495     306.4615      0.2725754   0.959293624
# 2 0.9564822    134.245569     519.7108      0.5929597   1.010322884
# 7 0.8535633     -5.307428     635.0605      1.3928018   1.156047013
# 9 0.5315166     29.392924    1195.4566      1.5024837   1.195603989
# 8 0.5271784     29.307803    1205.0429      1.5049911   1.194418337
# 3 0.7291883   2088.368783    2088.3688      0.4344569   0.794622032
# 6 0.7291883   2088.368783    2088.3688      0.4344569   0.794622032
# 5 0.1385699 333370.032902  333370.0329      1.9641374   0.004393986
# 4 0.1395601 350083.518163  350083.5182      1.9789938   0.004496484
#   R2_selcoefs support_size_true support_size_estimated
# 1   0.9719381                50                     49
# 2   0.9645622                50                     55
# 7   0.8560704                50                    129
# 9   0.5341256                50                    132
# 8   0.5297660                50                    132
# 3   0.8924053                50                     14
# 6   0.8924053                50                     14
# 5   0.4167017                50                     27
# 4   0.4412216                50                     26
#   estimated_over_true_support_size TP FP  TN FN precision recall
# 1                             0.98 43  6 444  7 0.8775510   0.86
# 2                             1.10 42 13 438  8 0.7636364   0.84
# 7                             2.58 48 81 370  2 0.3720930   0.96
# 9                             2.64 47 85 366  3 0.3560606   0.94
# 8                             2.64 47 85 366  3 0.3560606   0.94
# 3                             0.28 13  1 449 37 0.9285714   0.26
# 6                             0.28 13  1 449 37 0.9285714   0.26
# 5                             0.54 14 13 437 36 0.5185185   0.28
# 4                             0.52 11 15 435 39 0.4230769   0.22
#   specificity  accuracy  F1_score          gic         ebic          bic
# 1   0.9866667 0.9740000 0.8686869 4.835465e+03 4.888176e+03 4.583660e+03
# 2   0.9711752 0.9580838 0.8000000 5.354018e+03 5.413202e+03 5.071179e+03
# 7   0.8203991 0.8343313 0.5363128 6.212316e+03 6.351130e+03 5.548930e+03
# 9   0.8115299 0.8243513 0.5164835 6.725468e+03 6.867510e+03 6.046654e+03
# 8   0.8115299 0.8243513 0.5164835 6.757121e+03 6.899163e+03 6.078308e+03
# 3   0.9977778 0.9240000 0.4062500 1.418324e+07 1.418325e+07 1.418317e+07
# 6   0.9977778 0.9240000 0.4062500 1.418324e+07 1.418325e+07 1.418317e+07
# 5   0.9711111 0.9020000 0.3636364 1.039001e+08 1.039002e+08 1.039000e+08
# 4   0.9666667 0.8920000 0.2894737 1.052087e+08 1.052087e+08 1.052086e+08
#            aic         mbic
# 1 4.377145e+03 5.056835e+03
# 2 4.839375e+03 5.602513e+03
# 7 5.005246e+03 6.795150e+03
# 9 5.490326e+03 7.321857e+03
# 8 5.521979e+03 7.353510e+03
# 3 1.418311e+07 1.418330e+07
# 6 1.418311e+07 1.418330e+07
# 5 1.038999e+08 1.039003e+08
# 4 1.052085e+08 1.052088e+08

# TODO: should also include timings (but L0glm is much faster than all the rest)
# TODO: maybe also include relaxed LASSO?
# -> L0glm beats all the other methods!


## 3.2 SINGLE CHANNEL GAUSSIAN ####
n <- p <- 500
npeaks <- 50
sim.gaus <- simul(method = 1, n = n, p = p, s = npeaks/p, seed = 2)

### 3.2.1 L0GLM VIA CHOLESKY ####

system.time(tune.chol.sc <- L0glm.tune(X=sim.gaus$X, y=sim.gaus$y, IC = "gic", 
                                       plot = FALSE, coefs.true = sim.gaus$beta_true, 
                                       relaxed = FALSE, family = gaussian(identity), 
                                       nonnegative = TRUE, 
                                       solver = "eigen", method = 7L))
# 0.36s

system.time(L0glm.chol <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                    family = gaussian(identity), 
                                    lambda = "gic", # or tune.chol.sc$lambda_best,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "eigen",
                                    method = 7L)) 
# 0.03s for n=p=500

plot(x=L0glm.chol$coefficients, y=sim.gaus$beta_true, pch=16, col='blue',
     main="L0 glm identity link gaussian with lambda tuned to optimise GIC",
     xlab="estimated coefficients", ylab="true coefficients")


### 3.2.2 LASSO USING GLMNET ####
system.time(glmnet.sc.gaus <- fit.GLM(sim = sim.gaus, method = 2, 
                                      family = gaussian(identity), plot = TRUE)) 
# 0.16s for n=p=500 


### 3.2.3 ADAPTIVE LASSO USING GLMNET ####
system.time(glmnet.al.sc.gaus <- fit.GLM(sim = sim.gaus, method = 3, 
                                         family = gaussian(identity), plot = T)) 
# 0.16s for n=p=500 


### 3.2.4 L0LEARN ####

system.time(L0learn.sc.gaus <- fit.GLM(sim = sim.gaus, method = 4, 
                                       family = gaussian(identity), plot = T))  
#  0.14s for n=p=500


### 3.2.5 L0LEARN AFTER LASSO ####
# Better than regular L0Learn & faster
system.time(L0learn.after.lasso.sc.gaus <- fit.GLM(sim = sim.gaus, method = 5, 
                                                   family = gaussian(identity), plot = T))  
# 0.31s for n=p=500 


### 3.2.6 ABESS ####
system.time(abess.sc.gaus <- fit.GLM(sim = sim.gaus, method = 6, 
                                     family = gaussian(identity), 
                                     plot = T, IC = "gic"))  
# 1.29s for n=p=500


### 3.2.7 ABESS AFTER LASSO ####
system.time(abess.after.lasso.sc.gaus <- fit.GLM(sim = sim.gaus, method = 7, 
                                                 family = gaussian(identity), 
                                                 plot = T, IC = "gic"))  
# 0.17s for n=p=500


### 3.2.8 MCP USING ORDINIS ####
system.time(ordinis.mcp.sc.gaus <- fit.GLM(sim = sim.gaus, method = 8, 
                                           family = gaussian(identity), 
                                           plot = T, IC = "gic"))  
# 0.27s for n=p=500


### 3.2.9 SCAD USING ORDINIS ####
system.time(ordinis.scad.sc.gaus <- fit.GLM(sim = sim.gaus, method = 9, 
                                            family = gaussian(identity), 
                                            plot = T, IC = "gic"))  
# 0.28s for n=p=500


### 3.2.10 BENCHMARKS TABLE ####

X <- sim.gaus$X
y <- sim.gaus$y
IC <- "gic"

benchmarks_3_2 = rbind(
  data.frame(method="L0glm", do.call(cbind, L0glm:::get.metrics(beta = L0glm.chol$coefficients,  # some metrics
                                                                     X = X, y = y, 
                                                                     family = gaussian(identity), 
                                                                     beta_true = sim.gaus$beta_true))),
  data.frame(method="nnalasso_glmnet", do.call(cbind, L0glm:::get.metrics(beta = glmnet.al.sc.gaus,  # some metrics
                                                                  X = cbind(1, X), y = y, 
                                                                  family = gaussian(identity), 
                                                                  beta_true = c(0, as.vector(sim.gaus$beta_true))) )),
  data.frame(method="nnL0Learn_after_glmnet_nnlasso", do.call(cbind, L0glm:::get.metrics(beta = L0learn.after.lasso.sc.gaus,  # some metrics
                                                                                 X = X, y = y, 
                                                                                 family = gaussian(identity),
                                                                                 beta_true = as.vector(sim.gaus$beta_true)) )),
  data.frame(method="abess_after_glmnet_nnlasso", do.call(cbind, L0glm:::get.metrics(beta = abess.after.lasso.sc.gaus,  # some metrics
                                                                             X = X, y = y, 
                                                                             family = gaussian(identity),
                                                                             beta_true = as.vector(sim.gaus$beta_true)) )),
  data.frame(method="nnL0Learn", do.call(cbind, L0glm:::get.metrics(beta = L0learn.sc.gaus,  # some metrics
                                                            X = X, y = y, 
                                                            family = gaussian(identity), 
                                                            beta_true = as.vector(sim.gaus$beta_true)) )),
  data.frame(method="abess", do.call(cbind, L0glm:::get.metrics(beta = abess.sc.gaus,  # some metrics
                                                        X = X, y = y, 
                                                        family = gaussian(identity), 
                                                        beta_true = as.vector(sim.gaus$beta_true)) )),
  data.frame(method="nnlasso_glmnet", do.call(cbind, L0glm:::get.metrics(beta = glmnet.sc.gaus,  # some metrics
                                                                 X = cbind(1, X), y = y, 
                                                                 family = gaussian(identity), 
                                                                 beta_true = c(0, as.vector(sim.gaus$beta_true))) )),
  data.frame(method="nnMCP_ordinis", do.call(cbind, L0glm:::get.metrics(beta = ordinis.mcp.sc.gaus,  # some metrics
                                                                X = cbind(1, X), y = y, 
                                                                family = gaussian(identity),
                                                                beta_true = c(0, as.vector(sim.gaus$beta_true))) )),
  data.frame(method="nnSCAD_ordinis", do.call(cbind, L0glm:::get.metrics(beta = ordinis.scad.sc.gaus,  # some metrics
                                                                 X = cbind(1, X), y = y, 
                                                                 family = gaussian(identity), 
                                                                 beta_true = c(0, as.vector(sim.gaus$beta_true))) )))

benchmarks_3_2 = benchmarks_3_2[order(benchmarks_3_2[,"deviance_true"]),] # sort by deviance relative to ground truth
benchmarks_3_2
#                           method   p   n m   deviance deviance_true
# 3 nnL0Learn_after_glmnet_nnlasso 501 500 1   13547.22      6733.678
# 1                          L0glm 501 500 1   31125.10     19760.389
# 5                      nnL0Learn 501 500 1  232710.80    221496.014
# 9                 nnSCAD_ordinis 502 500 1  508208.93    491714.228
# 8                  nnMCP_ordinis 502 500 1 1273567.56   1249035.323
# 4     abess_after_glmnet_nnlasso 501 500 1 3258284.46   3269706.782
# 2                nnalasso_glmnet 502 500 1 3665597.68   3667508.250
# 7                 nnlasso_glmnet 502 500 1 4338167.47   4350827.911
# 6                          abess 501 500 1 7116162.10   7116773.335
#      loglik      MD_coefs MAD_coefs sMAPE_coefs bias_coefs  R2_coefs
# 3 -1535.301  -0.007830595 113.60805   0.2313449  0.9414228 0.8942910
# 1 -1743.260  -0.003500867  38.66343   0.0519593  0.9881805 0.9739436
# 5 -2246.205   0.411679765 146.54041   0.2057775  0.9247740 0.8441756
# 9 -2441.479   0.601437492 131.55853   0.1594143  0.9207443 0.8543947
# 8 -2671.150   1.123362852 260.28405   0.2713568  0.8101574 0.5979304
# 4 -2905.995   0.940147627 176.14620   0.1669200  0.8598522 0.7867479
# 2 -2935.443  -6.096004919  95.93336   0.1134316  0.9764912 0.9301852
# 7 -2977.558 -11.601690959 125.95294   0.3323378  1.0491874 0.9274381
# 6 -3101.287   2.076538493 186.68587   0.1757650  0.8381777 0.7684228
#   MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs
# 3   105.25784     518.9192      1.1392574     0.9414228   0.9093948
# 1    95.23115     276.4989      0.3073010     0.9881805   0.9791516
# 5   473.67692     622.9805      1.0580407     0.9247740   0.9178036
# 9   576.38337     605.9764      0.8547599     0.9207443   0.9355237
# 8   959.01220    1055.8357      1.2915344     0.8101574   0.8073822
# 4   916.22988     928.4196      0.8637499     0.8598522   0.8957397
# 2   223.27211     695.2690      0.6291677     0.9764912   0.9447399
# 7   -39.80944     529.4674      1.4059936     1.0491874   0.9285669
# 6   979.49260     986.7825      0.9142191     0.8381777   0.8803839
#   support_size_true support_size_estimated
# 3                50                     91
# 1                50                     52
# 5                50                     67
# 9                50                     56
# 8                50                     65
# 4                50                     48
# 2                50                     49
# 7                50                    117
# 6                50                     48
#   estimated_over_true_support_size TP FP  TN FN precision recall
# 3                             1.82 44 47 403  6 0.4835165   0.88
# 1                             1.04 45  7 443  5 0.8653846   0.90
# 5                             1.34 34 33 417 16 0.5074627   0.68
# 9                             1.12 34 22 429 16 0.6071429   0.68
# 8                             1.30 24 41 410 26 0.3692308   0.48
# 4                             0.96 29 19 431 21 0.6041667   0.58
# 2                             0.98 37 12 439 13 0.7551020   0.74
# 7                             2.34 49 68 383  1 0.4188034   0.98
# 6                             0.96 28 20 430 22 0.5833333   0.56
#   specificity  accuracy  F1_score      gic     ebic      bic      aic
# 3   0.8955556 0.8940000 0.6241135 4115.459 4214.458 3642.347 3254.603
# 1   0.9844444 0.9760000 0.8823529 4088.448 4145.480 3815.894 3592.519
# 5   0.9266667 0.9020000 0.5811966 5264.695 5337.869 4915.004 4628.410
# 9   0.9512195 0.9241517 0.6415094 5530.523 5591.879 5237.191 4996.958
# 8   0.9090909 0.8662675 0.4173913 6092.112 6163.156 5752.465 5474.301
# 4   0.9577778 0.9200000 0.5918367 6368.490 6421.218 6116.506 5909.990
# 2   0.9733925 0.9500998 0.7474747 6438.924 6492.745 6181.616 5970.885
# 7   0.8492239 0.8622754 0.5868263 7295.687 7422.705 6688.440 6191.116
# 6   0.9555556 0.9160000 0.5714286 6759.074 6811.801 6507.090 6300.574
#       mbic
# 3 4531.124
# 1 4327.907
# 5 5571.926
# 9 5788.074
# 8 6390.329
# 4 6589.877
# 2 6664.846
# 7 7828.864
# 6 6980.460

# TODO: should also include timings (but L0glm is much faster than all the rest)
# deviance relative to ground truth best for nnL0Learn_after_glmnet_nnlasso but
# F1_score & gic much better for L0glm
# -> L0glm is the winner, both in terms of timings & quality of fit (esp support recovery)


## 3.3 MULTICHANNEL POISSON ####
n <- p <- 500
npeaks <- 50
sim.pois.mc <- simul(method = 4, n = n, p = p, s = npeaks/p, seed = 2)

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
# Best lambda by GIC = 0.976, 232s for n=p=500 with parallelize=F

system.time(L0glm.chol.pois <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                         family = poisson(identity), 
                                         lambda = "gic", # or tune.chol.mc$lambda_best
                                         nonnegative = TRUE,
                                         grouped = TRUE,
                                         relaxed = TRUE,
                                         parallelize = FALSE,
                                         solver = "eigen",
                                         method = 7L)) 
# 21s

x_vals <- as.vector(rowMeans(L0glm.chol.pois$coefficients))        # estimated
y_vals <- as.vector(rowMeans(sim.pois.mc$beta_true))               # true
fit <- lm(y_vals ~ 0 + x_vals)
R2p <- round(100 * summary(fit)$r.squared, 2)
main_expr <- bquote(
  atop("true vs estimated coefficients (row means)",
       atop("relaxed nonnegative grouped L0 Poisson GLM",
            R^2 == .(sprintf("%.2f", R2p)) * "%")))
plot(x_vals, y_vals, pch = 16, col = "steelblue",
     xlab = "estimated coefficients", 
     ylab = "true coefficients",
     main = main_expr)


plotmat(L0glm.chol.pois$coefficients,
      xlab="Time index", ylab="Channel", 
      main="nonnegative relaxed L0 group GLM\n(red = true peaks, cyan = estimated peaks)")
# abline(v=(1:nrow(sim.pois.mc$X))[as.vector(rowMax(sim.pois.mc$beta_true)!=0)], col="red")
# abline(v=(1:nrow(sim.pois.mc$X))[as.vector(rowMin(L0glm.chol.pois$coefficients)!=0)], col="cyan")
# TODO Fix bug Error: unable to find an inherited method for function ‘rowQ’ for signature ‘imat = "dgCMatrix", which = "numeric"’


### 3.3.2 BEST SUBSET OF GROUP SELECTION (BSGS) USING ABESS ####
system.time(bsgs.abess.pois <- fit.GLM.mc(sim = sim.pois.mc, method = 1, 
                                             family = poisson(identity), 
                                             plot = T, IC = "gic")) 
# 4.2s for n=p=500
# NOTE: this is now fit as if noise was Gaussian
# TODO also compare what solution & timings you get with
# https://skscope.readthedocs.io/en/latest/gallery/GeneralizedLinearModels/poisson-identity-link.html
# maybe implement this as an additional fitting method in 
# fit.GLM.mc?




### 3.3.3 MULTIVARIATE GROUP LASSO USING GLMNET ####
system.time(glmnet.mc.pois <- fit.GLM.mc(sim = sim.pois.mc, method = 2, 
                                    family = poisson(identity), plot = T, 
                                    IC = "gic"))
# 32.75s for n=p=500

plotmat(glmnet.mc.pois, 
      xlab="Time index", ylab="Channel", 
      main="glmnet nonnegative multivariate group LASSO\n(red = true peaks, cyan = estimated peaks)")
# abline(v=(1:nrow(sim.pois.mc$X))[as.vector(rowMax(sim.pois.mc$beta_true)!=0)], col="red")
# abline(v=(1:nrow(sim.pois.mc$X))[as.vector(rowMin(glmnet.mc.pois)!=0)], col="cyan")




### 3.3.4 BENCHMARKS TABLE ####

benchmarks_3_3 = rbind(
  data.frame(method="relaxed_group_L0glm_poisson", do.call(cbind, L0glm:::get.metrics(beta = as.matrix(L0glm.chol.pois$coefficients),  # TODO: make sure that metrics work with either dense or sparse coefs
                                                                        X = sim.pois.mc$X, y = sim.pois.mc$y,
                                                                        family = poisson(identity), 
                                                                        beta_true = sim.pois.mc$beta_true) ) ),
  data.frame(method="BSGS_abess_poisson", do.call(cbind, L0glm:::get.metrics(beta = bsgs.abess.pois, 
                                                                     X = sim.pois.mc$X, y = sim.pois.mc$y, 
                                                                     family = poisson(identity), 
                                                                     beta_true = sim.pois.mc$beta_true) ) ),
  data.frame(method="multivariate_group_LASSO_poisson", do.call(cbind, L0glm:::get.metrics(beta = glmnet.mc.pois,  
                                                                     X = sim.pois.mc$X, y = sim.pois.mc$y, 
                                                                     family = poisson(identity), 
                                                                     beta_true = sim.pois.mc$beta_true) ) ) 
  )
benchmarks_3_3 = benchmarks_3_3[order(benchmarks_3_3[,"deviance_true"]),] # sort by deviance relative to ground truth
benchmarks_3_3
#                             method   p   n   m     deviance
# 1      relaxed_group_L0glm_poisson 500 500 168     40521.27
# 3 multivariate_group_LASSO_poisson 500 500 168     94770.17
# 2               BSGS_abess_poisson 500 500 168 153465982.78
#   deviance_true      loglik     MD_coefs MAD_coefs sMAPE_coefs
# 1      10894.17   -102211.0  0.005847954 0.3961107  0.02496972
# 3      63018.86   -129335.4 -0.443704074 1.9998690  0.09829645
# 2  153222097.46 -76814941.7 -4.244634011 6.7621233  0.08466689
#   bias_coefs  R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs
# 1  0.9992973 0.9992087    2.240395     7.387358      0.3585562
# 3  1.0471089 0.9818982   -3.782244    24.653687      1.1730644
# 2  0.8237457 0.5912379   99.797468   101.793668      1.0962998
#   bias_selcoefs R2_selcoefs support_size_true support_size_estimated
# 1     0.9992973   0.9994279                50                     48
# 3     1.0471089   0.9824858                50                     66
# 2     0.8237457   0.8375311                50                      9
#   estimated_over_true_support_size TP FP  TN FN precision recall
# 1                             0.96 48  0 450  2 1.0000000   0.96
# 3                             1.32 46 20 430  4 0.6969697   0.92
# 2                             0.18  8  1 449 42 0.8888889   0.16
#   specificity accuracy  F1_score         gic        ebic         bic
# 1   1.0000000    0.996 0.9795918    244238.6    248011.2    226216.6
# 3   0.9555556    0.952 0.7931034    329777.7    336515.0    297592.9
# 2   0.9977778    0.914 0.2711864 153641793.3 153642921.7 153636402.6
#           aic        mbic
# 1    211435.9    260082.4
# 3    271196.8    358072.4
# 2 153631981.5 153646532.4




## 3.4 MULTICHANNEL GAUSSIAN ####
n <- p <- 500
npeaks <- 50
sim.gaus.mc <- simul(method = 3, n = n, p = p, s = npeaks/p, seed = 1)

### 3.4.1 RELAXED GROUP L0GLM VIA CHOLESKY ####

system.time(tune.chol.mc <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, IC = "gic", 
                                       plot = FALSE, coefs.true = sim.gaus.mc$beta_true, 
                                       family = gaussian(identity), 
                                       nonnegative = TRUE, 
                                       grouped = TRUE,
                                       relaxed = TRUE, 
                                       parallelize = FALSE, 
                                       solver = "eigen", 
                                       method = 7L)) 
# Best lambda by GIC = 11.4 = identical to preset lambda
tune.chol.mc$lambda_best # 11.35348
tune.chol.mc$lambda_preset # 11.35348

system.time(L0glm.chol.gaus <- L0glm.fit(sim.gaus.mc$X, sim.gaus.mc$y, 
                                    family = gaussian(identity), 
                                    lambda = "gic", 
                                    nonnegative = TRUE,
                                    grouped = TRUE,
                                    relaxed = TRUE,
                                    parallelize = FALSE,
                                    solver = "eigen",
                                    method = 7L)) 
# 20s for n=p=500 with parallelize = FALSE
       
x_vals <- as.vector(rowMeans(L0glm.chol.gaus$coefficients))        # estimated
y_vals <- as.vector(rowMeans(sim.gaus.mc$beta_true))               # true
fit <- lm(y_vals ~ 0 + x_vals)
R2p <- round(100 * summary(fit)$r.squared, 2)
main_expr <- bquote(
  atop("true vs estimated coefficients (row means)",
       atop("relaxed nonnegative grouped L0 Gaussian GLM",
            R^2 == .(sprintf("%.2f", R2p)) * "%")))
plot(x_vals, y_vals, pch = 16, col = "steelblue",
     xlab = "estimated coefficients", 
     ylab = "true coefficients",
     main = main_expr)

plotmat(L0glm.chol.gaus$coefficients,
      xlab="Time index", ylab="Channel", 
      main="relaxed nonnegative group L0glm\n(red = true peaks, cyan = estimated peaks)")
# abline(v=(1:nrow(sim.gaus.mc$X))[as.vector(rowMax(sim.gaus.mc$beta_true)!=0)], col="red")
# abline(v=(1:nrow(sim.gaus.mc$X))[as.vector(rowMin(L0glm.chol.gaus$coefficients)!=0)], col="cyan")



### 3.4.2 BEST SUBSET OF GROUP SELECTION (BSGS) USING ABESS ####
system.time(bsgs.abess.gaus <- fit.GLM.mc(sim = sim.gaus.mc, method = 1, 
                                          family = gaussian(identity), plot = T, 
                                          IC = "gic")) 
# 5s for n=p=500
# TODO post this as a reproducible example in https://github.com/abess-team/abess/issues/510
# and ask if the solution quality could be improved here somehow



# glmnet nonnegative mgaussian multivariate group LASSO fit
system.time(glmnet.mc.gaus <- fit.GLM.mc(sim = sim.gaus.mc, method = 2, 
                                    family = gaussian(identity), plot = T, 
                                    IC = "gic"))
# 42s for n=p=500

plotmat(glmnet.mc.gaus,
      xlab="Time index", ylab="Channel", 
      main="nonnegative multivariate group LASSO\n(red = true peaks, cyan = estimated peaks)")
# abline(v=(1:nrow(sim.gaus.mc$X))[as.vector(rowMax(sim.gaus.mc$beta_true)!=0)], col="red")
# abline(v=(1:nrow(sim.gaus.mc$X))[as.vector(rowMin(glmnet.mc.gaus)!=0)], col="cyan")




### 3.4.3 BENCHMARKS TABLE ####

benchmarks_3_4 = rbind(
  data.frame(method="relaxed_group_L0glm_gaussian", do.call(cbind, L0glm:::get.metrics(beta = as.matrix(L0glm.chol.gaus$coefficients),  # TODO: allow coefs to be either dense or sparse matrix in get.metrics
                                                                         X = sim.gaus.mc$X, y = sim.gaus.mc$y, 
                                                                         family = gaussian(identity),
                                                                         beta_true = sim.gaus.mc$beta_true)) ),
  data.frame(method="BSGS_abess_gaussian", do.call(cbind, L0glm:::get.metrics(beta = bsgs.abess.gaus,  # some metrics
                                                                      X = sim.gaus.mc$X, y = sim.gaus.mc$y, 
                                                                      family = gaussian(identity),
                                                                      beta_true = sim.gaus.mc$beta_true)) ),
  data.frame(method="multivariate_group_LASSO_gaussian", do.call(cbind, L0glm:::get.metrics(beta = glmnet.mc.gaus,  
                                                                                   X = sim.gaus.mc$X, y = sim.gaus.mc$y, 
                                                                                   family = gaussian(identity), 
                                                                                   beta_true = sim.gaus.mc$beta_true)) ) )

benchmarks_3_4 = benchmarks_3_4[order(benchmarks_3_4[,"deviance_true"]),] # sort by deviance relative to ground truth
benchmarks_3_4
#                              method   p   n   m  deviance deviance_true
# 1      relaxed_group_L0glm_gaussian 501 500 178   1190956      226789.9
# 3 multivariate_group_LASSO_gaussian 501 500 178   9862040     8916586.6
# 2               BSGS_abess_gaussian 501 500 178 235538580   234573516.3
#      loglik   MD_coefs MAD_coefs sMAPE_coefs bias_coefs  R2_coefs
# 1 -241714.0  0.0995758 0.2712892   0.1007760  1.0004618 0.9997597
# 3 -335784.2 -0.1688281 2.5588307   0.2006945  1.0319417 0.9063851
# 2 -476990.7 -0.5410059 5.9083125   0.1270843  0.5640975 0.3437010
#   MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs
# 1  1.32991993     3.085139       1.137114     1.0004618   0.9997793
# 3 -0.04439605    18.778783       1.520152     1.0319417   0.9157136
# 2 54.26261191    55.173881       1.447997     0.5640975   0.5854438
#   support_size_true support_size_estimated
# 1                50                     51
# 3                50                     75
# 2                50                     25
#   estimated_over_true_support_size TP FP  TN FN precision recall
# 1                             1.02 50  1 449  0 0.9803922   1.00
# 3                             1.50 48 27 423  2 0.6400000   0.96
# 2                             0.50 20  5 445 30 0.8000000   0.40
#   specificity accuracy  F1_score       gic      ebic      bic      aic
# 1   0.9977778    0.998 0.9900990  568345.2  576391.1 529894.6 498381.9
# 3   0.9400000    0.942 0.7680000  800528.6  812747.4 742135.2 694278.3
# 2   0.9888889    0.930 0.5333333 1003566.7 1008264.9 981114.5 962713.5
#        mbic
# 1  602127.1
# 3  851831.6
# 2 1023292.8

# -> L0glm is the clear winner!
