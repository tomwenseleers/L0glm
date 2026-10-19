# BENCHMARKS L0GLM: SINGLE CHANNEL VS MULTICHANNEL DECONVOLUTION OF BLURRED SPIKE TRAIN ####
# Historical timings below are machine/build dependent.
# Historical numerical tables are retained for comparison. The Poisson table
# is reproduced by maxit = 25, tol = 1e-10 (the restored fitting defaults).
# More iterations can improve GIC while worsening recovery against the truth.
# Current L0glm.tune fits each lambda independently and verifies its numerical
# zero endpoint. That endpoint is not a global L0 likelihood bound: a zero-mean
# Poisson model has infinite deviance if any count is positive.
# Select an IC with tune = "gic" (formerly IC = "gic"). For cross-validation,
# use tune = "cv", nfolds = 5; the right panel then shows held-out deviance +/- SE.

library(remotes)
# remotes::install_github("jaredhuling/ordinis")
# remotes::install_github("tomwenseleers/L0glm/L0glm")
# devtools::load_all(".")
library(L0glm)
library(export)
library(Matrix)
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
# here preset lambda coincides with best lambda
system.time(tune_full <- L0glm.tune(X=sim.pois$X, 
                   y=sim.pois$y, 
                   tune = "gic",
                   nlambda = 10,
                   prop_range = 1, # full range up to lambda_max
                   plot = TRUE, 
                   coefs.true = sim.pois$beta_true, 
                   relaxed = FALSE, 
                   family = poisson(identity), 
                   nonnegative = TRUE, 
                   solver = "eigen", 
                   method = 1L)) # 4.05s
tune_full$lambda_best # 12.84577 (best lambda tune)
tune_full$lambda_preset # 20.44997 (lambda preset)
min(tune_full$metrics$gic) # 102185

## tuning plot (zoomed in)
# with lambda grid more zoomed in around lambda_preset preset
system.time(tune_zoomed <- L0glm.tune(X=sim.pois$X, 
           y=sim.pois$y, 
           tune = "gic",
           nlambda = 10,
           prop_range = 0.2, # to zoom in around lambda_preset preset
           plot = TRUE, 
           coefs.true = sim.pois$beta_true, 
           relaxed = FALSE, 
           family = poisson(identity), 
           nonnegative = TRUE, 
           solver = "eigen", 
           method = 1L)) # 4.36s
# The dashed line marks the best grid value; the dotted line is the preset.
# They need not coincide, especially for unrelaxed fits.
tune_zoomed$lambda_best # 12.84577 (best lambda tune)
tune_zoomed$lambda_preset # 20.44997 (lambda preset)
min(tune_zoomed$metrics$gic) # 102185

# The preset is an approximate calibration. A finer grid can improve GIC;
# neither the heuristic fitting algorithm nor the coarse grid certifies an optimum.
L0glm:::preset.lambda(IC="gic", p, n) # 20.44997 (lambda preset)
tune_zoomed$lambda_preset # 20.44997 (lambda preset)
tune_zoomed$lambda_best # 12.84577 (best lambda tune)

# note that with the tuning plot above that uses the eigen solver and
# enforced nonnegativity constraints merely by clipping, the support size
# size shrinks again if the level of regularisation (lambda) is chosen too small
# Around the lambda "gic" preset, however, the solution is excellent at a much
# smaller computational cost.
system.time(tune_full_bigglm <- L0glm.tune(X=sim.pois$X, 
                                    y=sim.pois$y, 
                                    tune = "gic",
                                    nlambda = 10,
                                    prop_range = 1, # full range up to lambda_max
                                    plot = TRUE, 
                                    coefs.true = sim.pois$beta_true, 
                                    relaxed = FALSE, 
                                    family = poisson(identity), 
                                    nonnegative = TRUE, 
                                    solver = "bigGlm"))
# 19.50s

system.time(tune_full_osqp <- L0glm.tune(X=sim.pois$X, 
                                           y=sim.pois$y, 
                                           tune = "gic",
                                           nlambda = 10,
                                           prop_range = 1, # full range up to lambda_max
                                           plot = TRUE, 
                                           coefs.true = sim.pois$beta_true, 
                                           relaxed = FALSE, 
                                           family = poisson(identity), 
                                           nonnegative = TRUE, 
                                           solver = "osqp"))
# 30.94s

# 5-fold cross validation is also supported:
system.time(tune_zoomed_cv <- L0glm.tune(X=sim.pois$X, 
                                      y=sim.pois$y, 
                                      tune = "cv",
                                      nfolds = 5,
                                      prop_range = 1, # to zoom in around lambda_preset preset
                                      lambda_min = 2,
                                      plot = TRUE, 
                                      coefs.true = sim.pois$beta_true, 
                                      relaxed = FALSE, 
                                      family = poisson(identity), 
                                      nonnegative = TRUE, 
                                      solver = "eigen", 
                                      method = 1L)) # 10.94s
tune_zoomed_cv$lambda_best # 2


#### 1.1.1 L0GLM FIT WITH CHOLESKY LLT SOLVER ####
system.time(L0glm.chol.pois <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                         family = poisson(identity), 
                                         lambda = "gic", # or tune_zoomed$lambda_best
                                         nonnegative = TRUE,
                                         relaxed = FALSE,
                                         # verbose = TRUE,
                                         solver = "eigen",
                                         method = 1L)) 
# 0.12s

#### 1.1.2 L0GLM FIT WITH CG SOLVER WITH Incomplete Cholesky LLT preconditioner ####
system.time(L0glm.CG.pois <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                       family = poisson(identity), 
                                       lambda = "gic", # or tune_zoomed$lambda_best
                                       nonnegative = TRUE,
                                       relaxed = FALSE,
                                       solver = "eigen",
                                       method = 3L)) 
# 0.33s

#### 1.1.3 L0GLM FIT WITH preconditioned LSCG SOLVER ####
system.time(L0glm.LSCG.pois <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                         family = poisson(identity), 
                                         lambda = "gic", # or tune_zoomed$lambda_best
                                         nonnegative = TRUE,
                                         relaxed = FALSE,
                                         solver = "eigen",
                                         method = 4L)) 
# 0.78s

#### 1.1.4 L0GLM FIT WITH OSQP SOLVER ####
system.time(L0glm.osqp.pois <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                    family = poisson(identity), 
                                    lambda = "gic", # or tune_zoomed$lambda_best
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "osqp")) 
# 4.61s

#### 1.1.5 L0GLM FIT WITH bigGlm::bigGlm SOLVER ####
system.time(L0glm.bigGlm.pois <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                      family = poisson(identity), 
                                      lambda = "gic", # or tune_zoomed$lambda_best
                                      nonnegative = TRUE,
                                      relaxed = FALSE,
                                      solver = "bigGlm")) 
# 1.74s with n=p=10000




### 1.1.9 BENCHMARKS TABLE ####

X = sim.pois$X
y = sim.pois$y
benchmarks_1_1 = rbind(
  data.frame(method="L0glm_osqp", do.call(cbind, L0glm:::get.metrics(beta = L0glm.osqp.pois$coefficients,  # some metrics
                                                             X = X, y = y, 
                                                             family = poisson(identity), 
                                                             beta_true = sim.pois$beta_true))),
  data.frame(method="L0glm_bigglm", do.call(cbind, L0glm:::get.metrics(beta = L0glm.bigGlm.pois$coefficients,  # some metrics
                                                               X = X, y = y, 
                                                               family = poisson(identity), 
                                                               beta_true = sim.pois$beta_true))),
  data.frame(method="L0glm_cholLLT", do.call(cbind, L0glm:::get.metrics(beta = L0glm.chol.pois$coefficients,  # some metrics
                                                                 X = X, y = y, 
                                                                 family = poisson(identity),
                                                                 beta_true = sim.pois$beta_true) )),
  data.frame(method="L0glm_PLSCG", do.call(cbind, L0glm:::get.metrics(beta = L0glm.LSCG.pois$coefficients,  # some metrics
                                                             X = X, y = y, 
                                                             family = poisson(identity),
                                                             beta_true = sim.pois$beta_true) )),
  data.frame(method="L0glm_PCG", do.call(cbind, L0glm:::get.metrics(beta = L0glm.CG.pois$coefficients,  # some metrics
                                                           X = X, y = y, 
                                                           family = poisson(identity), 
                                                           beta_true = sim.pois$beta_true) ))
)

benchmarks_1_1 = benchmarks_1_1[order(benchmarks_1_1[,"deviance_true"]),] # sort by deviance relative to ground truth
benchmarks_1_1
#          method     p     n m deviance deviance_true    loglik    MD_coefs MAD_coefs sMAPE_coefs bias_coefs  R2_coefs
# 4   L0glm_PLSCG 10000 10000 1 9349.103      2988.746 -41036.03 -0.02861098  80.28593  0.07161198  0.9928605 0.9093004
# 5     L0glm_PCG 10000 10000 1 9349.100      2988.750 -41036.03 -0.02861004  80.28585  0.07161198  0.9928607 0.9093003
# 3 L0glm_cholLLT 10000 10000 1 9349.100      2988.750 -41036.03 -0.02861004  80.28585  0.07161198  0.9928607 0.9093003
# 2  L0glm_bigglm 10000 10000 1 9422.418      3019.605 -41072.69 -0.02894092  83.72523  0.07347632  0.9936964 0.9045904
# 1    L0glm_osqp 10000 10000 1 9462.048      3074.622 -41092.50 -0.02870220  80.66691  0.07307096  0.9914590 0.9085972
#   MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true support_size_estimated
# 4    265.3200     528.5209      0.3977446     0.9928605   0.9400566              1000                   1011
# 5    265.3201     528.5201      0.3977446     0.9928607   0.9400565              1000                   1011
# 3    265.3201     528.5201      0.3977446     0.9928607   0.9400565              1000                   1011
# 2    243.4808     581.1132      0.4125746     0.9936964   0.9283964              1000                   1015
# 1    270.6446     526.1757      0.4058395     0.9914590   0.9402930              1000                   1012
#   estimated_over_true_support_size  TP  FP   TN  FN precision recall specificity accuracy  F1_score       MCC      gic
# 4                            1.011 843 168 8832 157 0.8338279  0.843   0.9813333   0.9675 0.8383889 0.8203380 102747.0
# 5                            1.011 843 168 8832 157 0.8338279  0.843   0.9813333   0.9675 0.8383889 0.8203380 102747.0
# 3                            1.011 843 168 8832 157 0.8338279  0.843   0.9813333   0.9675 0.8383889 0.8203380 102747.0
# 2                            1.015 842 173 8827 158 0.8295567  0.842   0.9807778   0.9669 0.8357320 0.8173569 102902.1
# 1                            1.012 840 172 8828 160 0.8300395  0.840   0.9808889   0.9668 0.8349901 0.8165520 102880.4
#       ebic      bic      aic     mbic
# 4 100695.4 91383.71 84094.06 107203.9
# 5 100695.4 91383.71 84094.06 107203.9
# 3 100695.4 91383.71 84094.06 107203.9
# 2 100842.4 91493.87 84175.38 107376.7
# 1 100826.7 91505.87 84209.01 107341.7

# cholesky LLT fastest & best solution quality (deviance_true)
# eigen solvers all equally good & cholesky the fastest
# osqp solver the worst


## 1.2. L0GLM FIT SINGLE CHANNEL WITH GAUSSIAN NOISE ####

## tuning plot (full range)
# here preset lambda coincides with best lambda
system.time(tune_full <- L0glm.tune(X=sim.gaus$X, 
                   y=sim.gaus$y, 
                   tune = "gic",
                   nlambda = 10,
                   prop_range = 1, # full range up to lambda_max
                   plot = TRUE, 
                   coefs.true = sim.gaus$beta_true, 
                   relaxed = FALSE, 
                   family = gaussian(identity), 
                   nonnegative = TRUE, 
                   solver = "eigen", 
                   method = 1L)) # 7.28
tune_full$lambda_best # 20.44997 (best lambda tune)
tune_full$lambda_preset # 20.44997 (lambda preset)
min(tune_full$metrics$gic) # 82342.78



## tuning plot (zoomed in)
# tuning plots - with lambda grid more zoomed in around lambda_preset preset
system.time(tune_zoomed <- L0glm.tune(X=sim.gaus$X, 
                                      y=sim.gaus$y, 
                                      tune = "gic",
                                      nlambda = 10,
                                      prop_range = 0.2, # to zoom in around lambda_preset preset
                                      plot = TRUE, 
                                      coefs.true = sim.gaus$beta_true, 
                                      relaxed = FALSE, 
                                      family = gaussian(identity), 
                                      nonnegative = TRUE, 
                                      solver = "eigen", 
                                      method = 1L)) # 6.74s

# The dashed line marks the best grid value; the dotted line is the preset.
# They need not coincide when dispersion and support are estimated heuristically.
tune_zoomed$lambda_best # 34.03505 (best lambda tune)
tune_zoomed$lambda_preset # 20.44997 (lambda preset)
min(tune_zoomed$metrics$gic) # 81840.43


# Compare the empirical grid choice with the approximate preset calibration.
L0glm:::preset.lambda(IC="gic", p, n) # 20.44997 (lambda preset)
tune_zoomed$lambda_best # 34.03505 (best lambda tune)
tune_zoomed$lambda_preset # 20.44997 (lambda preset)


#### 1.2.1 L0GLM FIT WITH Cholesky LLT SOLVER ####
system.time(L0glm.chol.gaus <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                         family = gaussian(identity), 
                                         lambda = "gic", # or tune_zoomed$lambda_best
                                         nonnegative = TRUE,
                                         relaxed = FALSE,
                                         solver = "eigen",
                                         method = 1L)) 
# 0.28s

#### 1.2.2 L0GLM FIT WITH PCG SOLVER WITH Incomplete Cholesky preconditioner ####
system.time(L0glm.CG.gaus <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                       family = gaussian(identity), 
                                       lambda = "gic", # or tune_zoomed$lambda_best
                                       nonnegative = TRUE,
                                       relaxed = FALSE,
                                       solver = "eigen",
                                       method = 3L)) 
# 0.49s


#### 1.2.3 L0GLM FIT WITH PLSCG SOLVER ####
system.time(L0glm.LSCG.gaus <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                         family = gaussian(identity), 
                                         lambda = "gic", # or tune_zoomed$lambda_best
                                         nonnegative = TRUE,
                                         relaxed = FALSE,
                                         solver = "eigen",
                                         method = 4L)) 
# 4.44s


#### 1.2.4 L0GLM FIT WITH OSQP SOLVER ####
system.time(L0glm.osqp.gaus <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                    family = gaussian(identity), 
                                    lambda = "gic", # or tune_zoomed$lambda_best
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "osqp")) 
# 52.24s


#### 1.2.5 L0GLM FIT WITH bigGlm::bigGlm SOLVER ####
system.time(L0glm.bigGlm.gaus <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                      family = gaussian(identity), 
                                      lambda = "gic", # or tune_zoomed$lambda_best
                                      nonnegative = TRUE,
                                      relaxed = FALSE,
                                      solver = "bigGlm")) 
# 7.04s



### 1.2.9 BENCHMARKS TABLE ####

X = sim.gaus$X
y = sim.gaus$y
benchmarks_1_2 = rbind(
  data.frame(method="L0glm_osqp", do.call(cbind, L0glm:::get.metrics(beta = L0glm.osqp.gaus$coefficients,  # some metrics
                                                             X = X, y = y, 
                                                             family = gaussian(identity), 
                                                             beta_true = sim.gaus$beta_true))),
  data.frame(method="L0glm_bigglm", do.call(cbind, L0glm:::get.metrics(beta = L0glm.bigGlm.gaus$coefficients,  # some metrics
                                                               X = X, y = y, 
                                                               family = gaussian(identity), 
                                                               beta_true = sim.gaus$beta_true))),
  data.frame(method="L0glm_cholLLT", do.call(cbind, L0glm:::get.metrics(beta = L0glm.chol.gaus$coefficients,  # some metrics
                                                                 X = X, y = y, 
                                                                 family = gaussian(identity),
                                                                 beta_true = sim.gaus$beta_true) )),
  data.frame(method="L0glm_PLSCG", do.call(cbind, L0glm:::get.metrics(beta = L0glm.LSCG.gaus$coefficients,  # some metrics
                                                             X = X, y = y, 
                                                             family = gaussian(identity),
                                                             beta_true = sim.gaus$beta_true) )),
  data.frame(method="L0glm_PCG", do.call(cbind, L0glm:::get.metrics(beta = L0glm.CG.gaus$coefficients,  # some metrics
                                                           X = X, y = y, 
                                                           family = gaussian(identity), 
                                                           beta_true = sim.gaus$beta_true) )) 
)

benchmarks_1_2 = benchmarks_1_2[order(benchmarks_1_2[,"deviance_true"]),] # sort by deviance relative to ground truth
benchmarks_1_2
#          method     p     n m deviance deviance_true    loglik      MD_coefs MAD_coefs sMAPE_coefs bias_coefs  R2_coefs
# 1    L0glm_osqp 10000 10000 1 205112.1      34600.50 -29294.24 -0.0056724685  15.73137  0.02272312   1.001199 0.9898175
# 3 L0glm_cholLLT 10000 10000 1 250925.1      84966.57 -30302.23 -0.0024510498  26.77883  0.03974676   1.015738 0.9805885
# 5     L0glm_PCG 10000 10000 1 252476.9      87052.45 -30333.06 -0.0016494480  26.91941  0.03981904   1.015348 0.9801367
# 4   L0glm_PLSCG 10000 10000 1 282882.2     114917.20 -30901.61  0.0001025945  28.74677  0.03839087   1.024859 0.9680017
# 2  L0glm_bigglm 10000 10000 1 676200.0     519661.97 -35258.90  0.0372800564 102.84467  0.10652868   1.043044 0.8920861
#   MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true support_size_estimated
# 1    33.04374     118.7487      0.1517676      1.001199   0.9914453              1000                   1036
# 3    61.99624     190.1353      0.2575024      1.015738   0.9839142              1000                   1062
# 5    65.73558     187.9662      0.2565414      1.015348   0.9839054              1000                   1061
# 4    67.45273     205.0295      0.2444633      1.024859   0.9720554              1000                   1055
# 2    97.11852     722.0053      0.6745914      1.043044   0.8981086              1000                   1256
#   estimated_over_true_support_size  TP  FP   TN  FN precision recall specificity accuracy  F1_score       MCC      gic
# 1                            1.036 965  71 8929  35 0.9314672  0.965   0.9921111   0.9894 0.9479371 0.9422211 79795.10
# 3                            1.062 938 124 8876  62 0.8832392  0.938   0.9862222   0.9814 0.9097963 0.8999435 82342.78
# 5                            1.061 937 124 8876  63 0.8831291  0.937   0.9862222   0.9813 0.9092673 0.8993430 82383.98
# 4                            1.055 937 118 8882  63 0.8881517  0.937   0.9868889   0.9819 0.9119221 0.9022452 83398.39
# 2                            1.256 891 365 8635 109 0.7093949  0.891   0.9594444   0.9526 0.7898936 0.7698698 96223.42
#       ebic      bic      aic      mbic
# 1 77690.73 68139.61 60662.49  84366.68
# 3 80185.65 70395.06 62730.47  87028.98
# 5 80228.88 70447.50 62790.12  87065.77
# 4 81255.46 71529.35 63915.23  88053.73
# 2 93672.60 82095.21 73031.81 101764.86

# -> here osqp solver gave best solution & bigGlm's bigGlm the worst
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
system.time(tune_full <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y,
           tune = "gic",
           nlambda = 10,
           prop_range = 1, # to zoom in around lambda_preset preset
           plot = TRUE, coefs.true = sim.pois.mc$beta_true, 
           relaxed = FALSE, family = poisson(identity), 
           nonnegative = TRUE, 
           parallelize = TRUE, 
           grouped = FALSE,
           solver = "eigen", method = 1L,
           ncores = 4)) # openmp 2.96s
tune_full$lambda_best # 4.005619 (best lambda tune)
tune_full$lambda_preset # 11.35348 (lambda preset)
L0glm:::preset.lambda(IC="gic", p, n) # 11.35348 (lambda preset)

## tuning plot (zoomed in)
# here preset lambda would be ca 5x higher than best lambda
system.time(tune_zoomed <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, 
           tune = "gic",
           nlambda = 10,
           prop_range = 0.2, # to zoom in around lambda_preset preset
           plot = TRUE, coefs.true = sim.pois.mc$beta_true, 
           relaxed = FALSE, family = poisson(identity), 
           nonnegative = TRUE, 
           parallelize = TRUE, 
           grouped = FALSE,
           solver = "eigen", method = 1L,
           ncores = 4)) # openmp 3.31s
tune_zoomed$lambda_best # 8.022493 (best lambda tune)
tune_zoomed$lambda_preset # 11.35348 (lambda preset)
L0glm:::preset.lambda(IC="gic", p, n) # 11.35348 (lambda preset)


#### Cholesky LLT ####
system.time(L0glm.chol <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                    family = poisson(identity), 
                                    lambda = "gic", # or tune_zoomed$lambda_best
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallelize = TRUE,
                                    grouped = FALSE,
                                    solver = "eigen",
                                    method = 1L)) 
# openmp 0.13s

x_vals <- as.vector(rowMeans(L0glm.chol$coefficients))        # estimated
y_vals <- as.vector(rowMeans(sim.pois.mc$beta_true))          # true
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
# time = 7.67s with parallelize=F, 5.2s with parallelize=T





### 2.1.2 GROUPED L0GLM FIT WITH relaxed = FALSE, grouped = TRUE (better) ####

## tuning plot (full range)
# here preset lambda would be ca 10x higher than best lambda
system.time(tune_full <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, 
                        tune = "gic",
                        nlambda = 10,
                        prop_range = 1, # to zoom in around lambda_preset preset
                        plot = TRUE, coefs.true = sim.pois.mc$beta_true, 
                        relaxed = FALSE, family = poisson(identity), 
                        nonnegative = TRUE, 
                        parallelize = TRUE, 
                        grouped = TRUE,
                        solver = "eigen", method = 1L,
                        ncores = 4)) # openmp 4.47s
tune_full$lambda_best # 2.830413 (best lambda tune)
tune_full$lambda_preset # 11.35348 (lambda preset)
L0glm:::preset.lambda(IC="gic", p, n) # 11.35348 (lambda preset)

## tuning plot (zoomed in)
# here preset lambda would be ca 5x higher than best lambda
system.time(tune_zoomed <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, 
                          tune = "gic",
                          nlambda = 10,
                          prop_range = 0.2, # to zoom in around lambda_preset preset
                          plot = TRUE, coefs.true = sim.pois.mc$beta_true, 
                          relaxed = FALSE, family = poisson(identity), 
                          nonnegative = TRUE, 
                          parallelize = TRUE, 
                          grouped = TRUE,
                          solver = "eigen", method = 1L,
                          ncores = 4)) # openmp 5.75s
tune_zoomed$lambda_best # 8.022493 (best lambda tune)
tune_zoomed$lambda_preset # 11.35348 (lambda preset)
L0glm:::preset.lambda(IC="gic", p, n) # 11.35348 (lambda preset)



#### Cholesky LLT ####
system.time(L0glm.chol <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                    family = poisson(identity), 
                                    lambda = "gic", # or tune_zoomed$lambda_best
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallelize = TRUE,
                                    grouped = TRUE,
                                    solver = "eigen",
                                    method = 1L,
                                    ncores = 4)) 
# openmp 0.28s



#### OSQP ####
system.time(L0glm.osqp <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                    family = poisson(identity), 
                                    lambda = "gic", # or tune_zoomed$lambda_best
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallelize = FALSE,
                                    grouped = TRUE,
                                    solver = "osqp"))
# 33.17s




## 2.1.3 GROUPED L0GLM FIT WITH relaxed = TRUE, grouped = TRUE (best) ####

## tuning plot (full range)
# here preset lambda would be ca 10x higher than best lambda
system.time(tune_full <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, 
                        tune = "gic",
                        nlambda = 10,
                        prop_range = 1, # to zoom in around lambda_preset preset
                        plot = TRUE, coefs.true = sim.pois.mc$beta_true, 
                        relaxed = TRUE, family = poisson(identity), 
                        nonnegative = TRUE, 
                        parallelize = TRUE, 
                        grouped = TRUE,
                        solver = "eigen", method = 1L,
                        ncores = 4)) # openmp 5s
tune_full$lambda_best # 2 (best lambda tune)
tune_full$lambda_preset # 11.35348 (lambda preset)
L0glm:::preset.lambda(IC="gic", p, n) # 11.35348 (lambda preset)

## tuning plot (zoomed in)
# here preset lambda would be ca 10x higher than best lambda
system.time(tune_zoomed <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, 
                          tune = "gic",
                          nlambda = 10,
                          prop_range = 0.2, # to zoom in around lambda_preset preset
                          plot = TRUE, coefs.true = sim.pois.mc$beta_true, 
                          relaxed = TRUE, family = poisson(identity), 
                          nonnegative = TRUE, 
                          parallelize = TRUE, 
                          grouped = TRUE,
                          solver = "eigen", method = 1L,
                          ncores = 4)) # openmp 5.66s
tune_zoomed$lambda_best # 8.022493 (best lambda tune)
tune_zoomed$lambda_preset # 11.35348 (lambda preset)
L0glm:::preset.lambda(IC="gic", p, n) # 11.35348 (lambda preset)


#### Cholesky LLT ####
system.time(L0glm.chol <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                    family = poisson(identity), 
                                    lambda = "gic", # or tune_zoomed$lambda_best
                                    nonnegative = TRUE,
                                    relaxed = TRUE,
                                    parallelize = TRUE,
                                    grouped = TRUE,
                                    solver = "eigen",
                                    method = 1L))  # openmp 0.28s

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
# 58s




## 2.2. MULTICHANNEL WITH GAUSSIAN NOISE ####

### 2.2.1 GROUPED L0GLM FIT WITH relaxed = FALSE, grouped = FALSE (not great) ####

## tuning plot (full range)
# here best lambda is identical to preset lambda
system.time(tune_full <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, 
                        tune = "gic",
                        nlambda = 10,
                        prop_range = 1, # to zoom in around lambda_preset preset
                        plot = TRUE, 
                        coefs.true = sim.gaus.mc$beta_true, 
                        relaxed = FALSE, 
                        family = gaussian(identity), 
                        nonnegative = TRUE, 
                        parallelize = FALSE, 
                        grouped = FALSE,
                        solver = "eigen", method = 1L)) # 85s
tune_full$lambda_best # 11.35348
tune_full$lambda_preset # 11.35348 
L0glm:::preset.lambda(IC="gic", p, n) # 11.35348

## tuning plot (zoomed in)
# here best lambda is identical to preset lambda
system.time(tune_zoomed <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, 
                          tune = "gic",
                          nlambda = 10,
                          prop_range = 0.2, # to zoom in around lambda_preset preset
                          plot = TRUE, 
                          coefs.true = sim.gaus.mc$beta_true, 
                          relaxed = FALSE, 
                          family = gaussian(identity), 
                          nonnegative = TRUE, 
                          parallelize = FALSE, 
                          grouped = FALSE,
                          solver = "eigen", method = 1L)) # 89.75s
tune_zoomed$lambda_best # 11.35348
tune_zoomed$lambda_preset # 11.35348 
L0glm:::preset.lambda(IC="gic", p, n) # 11.35348

#### Cholesky LLT ####
system.time(L0glm.chol <- L0glm.fit(sim.gaus.mc$X, sim.gaus.mc$y, 
                                    family = gaussian(identity), 
                                    lambda = "gic", # or tune_zoomed$lambda_best
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallelize = FALSE,
                                    grouped = FALSE,
                                    solver = "eigen",
                                    method = 1L)) 
# 5.94s

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
# 18.54s




### 2.2.2 GROUPED L0GLM FIT WITH relaxed = FALSE, grouped = TRUE (better) ####

## tuning plot (full range)
# here best lambda is identical to preset lambda
system.time(tune_full <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, 
                        tune = "gic",
                        nlambda = 10,
                        prop_range = 1, # to zoom in around lambda_preset preset
                        plot = TRUE, 
                        coefs.true = sim.gaus.mc$beta_true, 
                        relaxed = FALSE, 
                        family = gaussian(identity), 
                        nonnegative = TRUE, 
                        parallelize = FALSE, 
                        grouped = TRUE,
                        solver = "eigen", method = 1L)) # 242.34s
tune_full$lambda_best # 11.35348
tune_full$lambda_preset # 11.35348 
L0glm:::preset.lambda(IC="gic", p, n) # 11.35348

## tuning plot (zoomed in)
# here best lambda is identical to preset lambda
system.time(tune_zoomed <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, 
                          tune = "gic",
                          nlambda = 10,
                          prop_range = 0.2, # to zoom in around lambda_preset preset
                          plot = TRUE, 
                          coefs.true = sim.gaus.mc$beta_true, 
                          relaxed = FALSE, 
                          family = gaussian(identity), 
                          nonnegative = TRUE, 
                          parallelize = FALSE, 
                          grouped = TRUE,
                          solver = "eigen", method = 1L)) # 300.46s
tune_zoomed$lambda_best # 11.35348
tune_zoomed$lambda_preset # 11.35348 
L0glm:::preset.lambda(IC="gic", p, n) # 11.35348

#### Cholesky LLT ####
system.time(L0glm.chol <- L0glm.fit(sim.gaus.mc$X, sim.gaus.mc$y, 
                                    family = gaussian(identity), 
                                    lambda = "gic", # or tune.chol.mc$lambda_best,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallelize = FALSE,
                                    grouped = TRUE,
                                    solver = "eigen",
                                    method = 1L)) # 26.09s

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
# 52.86s



### 2.2.3 GROUPED L0GLM FIT WITH relaxed = TRUE, grouped = TRUE (best) ####

## tuning plot (full range)
# here best lambda is identical to preset lambda
system.time(tune_full <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, 
                        tune = "gic",
                        nlambda = 10,
                        prop_range = 1, # to zoom in around lambda_preset preset
                        plot = TRUE, 
                        coefs.true = sim.gaus.mc$beta_true, 
                        relaxed = TRUE, 
                        family = gaussian(identity), 
                        nonnegative = TRUE, 
                        parallelize = FALSE, 
                        grouped = TRUE,
                        solver = "eigen", method = 1L)) # 277.47s
tune_full$lambda_best # 11.35348 (best lambda tune)
tune_full$lambda_preset # 11.35348 (lambda preset)
L0glm:::preset.lambda(IC="gic", p, n) # 11.35348 (lambda preset)

## tuning plot (zoomed in)
# here best lambda is again very close to preset lambda
system.time(tune_zoomed <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, 
                          tune = "gic",
                          nlambda = 10,
                          prop_range = 0.2, # to zoom in around lambda_preset preset
                          plot = TRUE, 
                          coefs.true = sim.gaus.mc$beta_true, 
                          relaxed = TRUE, 
                          family = gaussian(identity), 
                          nonnegative = TRUE, 
                          parallelize = FALSE, 
                          grouped = TRUE,
                          solver = "eigen", method = 1L)) # 281.16s
tune_zoomed$lambda_best # 3.361205 (best lambda tune)
tune_zoomed$lambda_preset # 11.35348 (lambda preset)
L0glm:::preset.lambda(IC="gic", p, n) # 11.35348 (lambda preset)


#### Cholesky LLT ####
system.time(L0glm.chol <- L0glm.fit(sim.gaus.mc$X, sim.gaus.mc$y, 
                                    family = gaussian(identity), 
                                    lambda = "gic", # or tune_zoomed$lambda_best
                                    nonnegative = TRUE,
                                    relaxed = TRUE,
                                    parallelize = FALSE,
                                    grouped = TRUE,
                                    solver = "eigen",
                                    method = 1L)) 
# 23.99s

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
# 51.46s






# 3. BENCHMARKS AGAINST OTHER METHODS ####
# Take our best working solver & method and compare it against others: 
# glmnet (LASSO variants), L0Learn & abess.
# MCP and SCAD were unacceptably slow for these n = p = 10000 deconvolution
# problems: their ordinis backend requires dense designs. So we omit them here.

# comparison_methods <- c("L0glm", "abess", "L0Learn", "lasso", "relaxed_lasso",
#                         "adaptive_lasso", "adaptive_relaxed_lasso")

comparison_methods <- c("L0glm", "L0Learn", "lasso", "relaxed_lasso",
                        "adaptive_lasso", "adaptive_relaxed_lasso")
# TODO check why abess was unusably slow


## 3.1 SINGLE CHANNEL POISSON ####
n <- p <- 10000
npeaks <- 1000
sim.pois <- simul(method = 2, n = n, p = p, s = npeaks/p, seed = 2)

### 3.1.1 L0GLM VIA CHOLESKY LLT ####
system.time(tune.chol.sc <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, tune = "gic",
                                       plot = FALSE, coefs.true = sim.pois$beta_true, 
                                       relaxed = FALSE, family = poisson(identity), 
                                       nonnegative = TRUE, 
                                       solver = "eigen", method = 1L))
# 2.54s

system.time(L0glm.chol <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                    family = poisson(identity), 
                                    lambda = "gic", # or tune.chol.sc$lambda_best,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "eigen",
                                    method = 1L)) 
# 0.14s

plot(x=L0glm.chol$coefficients, y=sim.pois$beta_true, pch=16, col='blue',
     main="L0 glm identity link Poisson with lambda tuned to optimise GIC",
     xlab="estimated coefficients", ylab="true coefficients")



### 3.1.2 COMPARE SPARSE METHODS ####
# The L0glm candidate uses the same preset as the fit above; other methods
# are rescored with the common Poisson GIC. WLS is labelled in the table.
# glmnet uses 10 lambda values. ABESS needs consecutive support sizes and
# warm starts on this correlated design: 1:2000 here (twice the peak count).
system.time(comparison.sc.pois <- sparseglm(sim.pois, family = poisson(identity),
  method = comparison_methods, nlambda = 10L,
  tune = "gic", nonnegative = TRUE, intercept = FALSE, approximation = "wls",
  support.size = seq_len(min(p, n - 1L, 2 * npeaks)),
  lambda = list(L0glm = L0glm.chol$init.lambda),
  options = list(abess = list(warm.start = TRUE))))
plot(comparison.sc.pois)


# TODO if abess is included it is super slow - check what's wrong there


# BENCHMARK TABLE

benchmarks_3_1 <- L0glm.metrics(comparison.sc.pois)
benchmarks_3_1 <- benchmarks_3_1[order(benchmarks_3_1[,"MCC"], decreasing = TRUE),]
benchmarks_3_1[,!colSums(is.na(benchmarks_3_1))]
#                   method status  scope        approximation tune parameter selected_parameter tuning_score seconds
# 1                  L0glm     ok single                 none  gic    lambda       2.044997e+01     102765.2    0.14
# 6 adaptive_relaxed_lasso     ok single identity_poisson_wls  gic    lambda       9.646484e+04     122003.4    0.41
# 4          relaxed_lasso     ok single identity_poisson_wls  gic    lambda       1.044924e-01     122083.4    0.39
# 5         adaptive_lasso     ok single identity_poisson_wls  gic    lambda       1.245891e+04     143856.3    0.22
# 3                  lasso     ok single identity_poisson_wls  gic    lambda       4.850106e-03     145903.5    0.22
# 2                L0Learn     ok single identity_poisson_wls  gic    lambda       3.296910e-03  891202271.5    0.53
#   path_seconds cv_seconds note     p     n m   df     deviance deviance_true        loglik MAD_coefs RMSE_coefs  R2_coefs
# 1         0.14          0      10000 10000 1  999 9.607885e+03  2.665175e+03     -41167.85  66.33935   486.4906 0.9195383
# 6         0.41          0      10000 10000 1 1900 1.042069e+04  4.229059e+03     -41574.25 186.90603   757.4528 0.8049317
# 4         0.39          0      10000 10000 1 1904 1.041884e+04  4.226631e+03     -41573.32 187.63062   759.1949 0.8040385
# 5         0.22          0      10000 10000 1 2857 1.270289e+04  6.739153e+03     -42715.35 162.67725   626.9684 0.8653472
# 3         0.22          0      10000 10000 1 3046 1.088513e+04  4.980112e+03     -41806.47 163.04241   626.1371 0.8656842
# 2         0.52          0      10000 10000 1   14 8.911293e+08  8.913185e+08 -445600992.61 302.28576  1403.4651 0.3423337
#   support_size_true support_size_estimated  TP   FP   TN  FN precision recall   F1_score       MCC         gic        ebic
# 1              1000                    999 845  154 8846 155 0.8458458  0.845 0.84542271 0.8282571    102765.2    100738.0
# 6              1000                   1900 906  994 8006  94 0.4768421  0.906 0.62482759 0.6083765    122003.4    118147.8
# 4              1000                   1904 905  999 8001  95 0.4753151  0.905 0.62327824 0.6066987    122083.4    118219.6
# 5              1000                   2857 949 1908 7092  51 0.3321666  0.949 0.49209230 0.4894331    143856.3    138058.6
# 3              1000                   3046 950 2096 6904  50 0.3118844  0.950 0.46959960 0.4674395    145903.5    139722.3
# 2              1000                     14  13    1 8999 987 0.9285714  0.013 0.02564103 0.1034134 891202271.5 891202243.1
#            bic          aic        mbic
# 1     91536.83     84333.70    107169.3
# 6    100648.15     86948.50    130379.5
# 4    100683.14     86954.65    130477.1
# 5    111744.65     91144.70    156451.2
# 3    111667.64     89704.94    159331.7
# 2 891202114.16 891202013.22 891202333.2

# L0glm.fit is the clear winner



## 3.2 SINGLE CHANNEL GAUSSIAN ####
n <- p <- 10000
npeaks <- 1000
sim.gaus <- simul(method = 1, n = n, p = p, s = npeaks/p, seed = 2)

### 3.2.1 L0GLM VIA CHOLESKY LLT ####

system.time(tune.chol.sc <- L0glm.tune(X=sim.gaus$X, y=sim.gaus$y, tune = "gic",
                                       plot = FALSE, coefs.true = sim.gaus$beta_true, 
                                       relaxed = FALSE, family = gaussian(identity), 
                                       nonnegative = TRUE, 
                                       solver = "eigen", method = 1L))
# 5.29s

system.time(L0glm.chol <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                    family = gaussian(identity), 
                                    lambda = "gic", # or tune.chol.sc$lambda_best,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "eigen",
                                    method = 1L)) 
# 0.25s for n=p=10000

plot(x=L0glm.chol$coefficients, y=sim.gaus$beta_true, pch=16, col='blue',
     main="L0 glm identity link gaussian with lambda tuned to optimise GIC",
     xlab="estimated coefficients", ylab="true coefficients")


### 3.2.2 COMPARE SPARSE METHODS ####
# The L0glm candidate uses the same preset as the fit above; other methods
# share GIC scoring. ABESS uses unit support increments with warm starts;
# glmnet uses 10 lambda values.
system.time(comparison.sc.gaus <- sparseglm(sim.gaus, family = gaussian(identity),
  method = comparison_methods, nlambda = 10L,
  tune = "gic", nonnegative = TRUE, intercept = FALSE,
  support.size = seq_len(min(p, n - 1L, 2 * npeaks)),
  lambda = list(L0glm = L0glm.chol$init.lambda),
  options = list(abess = list(warm.start = TRUE))))
plot(comparison.sc.gaus)
benchmarks_3_2 <- L0glm.metrics(comparison.sc.gaus)
benchmarks_3_2[,!colSums(is.na(benchmarks_3_2))]

# BENCHMARK TABLE

benchmarks_3_2 <- L0glm.metrics(comparison.sc.gaus)
benchmarks_3_2 <- benchmarks_3_2[order(benchmarks_3_2[,"MCC"], decreasing = TRUE),]
benchmarks_3_2[,!colSums(is.na(benchmarks_3_2))]
#                   method status  scope approximation tune parameter selected_parameter tuning_score seconds path_seconds
# 1                  L0glm     ok single          none  gic    lambda       2.044997e+01     81098.91    0.25         0.25
# 6 adaptive_relaxed_lasso     ok single          none  gic    lambda       3.034891e+11    145973.36    0.28         0.28
# 5         adaptive_lasso     ok single          none  gic    lambda       3.034891e+11    157320.70    0.30         0.30
# 4          relaxed_lasso     ok single          none  gic    lambda       3.543942e+00    148023.42    0.44         0.44
# 3                  lasso     ok single          none  gic    lambda       1.644952e-01    170989.82    0.14         0.13
# 2                L0Learn     ok single          none  gic    lambda       6.545630e-03    190173.27    0.73         0.73
#         cv_seconds note     p     n m   df     deviance deviance_true    loglik MAD_coefs RMSE_coefs  R2_coefs support_size_true
# 1          0      10000 10000 1 1043 2.308267e+05  5.598291e+04 -29884.80  19.06869   216.2561 0.9839711              1000
# 6          0      10000 10000 1  817 2.406940e+08  2.405623e+08 -64632.87 103.74671   552.6626 0.8952481              1000
# 5          0      10000 10000 1  836 7.201108e+08  7.202037e+08 -70112.26 109.22615   537.7683 0.9028029              1000
# 4          0      10000 10000 1 2298 1.429475e+07  1.412862e+07 -50514.70 242.42759   865.5736 0.7513042              1000
# 3          0      10000 10000 1 3457 1.328163e+07  1.313358e+07 -50147.15 209.48527   707.9322 0.8281221              1000
# 2          0      10000 10000 1   13 1.035381e+11  1.035362e+11 -94953.71 300.30677  1388.1245 0.3479385              1000
#   support_size_estimated  TP   FP   TN  FN precision recall   F1_score       MCC       gic      ebic       bic       aic
# 1                   1042 949   93 8907  51 0.9107486  0.949 0.92948090 0.9217078  81098.91  78982.36  69375.98  61855.59
# 6                    816 672  144 8856 328 0.8235294  0.672 0.74008811 0.7188927 145973.36 144315.43 136790.59 130899.74
# 5                    835 674  161 8839 326 0.8071856  0.674 0.73460490 0.7115234 157320.70 155624.21 147924.37 141896.52
# 4                   2297 892 1405 7595 108 0.3883326  0.892 0.54109797 0.5248353 148023.42 143360.12 122194.76 105625.40
# 3                   3456 928 2528 6472  72 0.2685185  0.928 0.41651706 0.4082175 170989.82 163974.59 132134.44 107208.29
# 2                     12  12    0 9000 988 1.0000000  0.012 0.02371542 0.1039855 190173.27 190146.89 190027.15 189933.42
#        mbic
# 1  85696.94
# 6 149575.08
# 5 161006.17
# 4 158154.08
# 3 186229.89
# 2 190230.58

# L0glm.fit is the clear winner



## 3.3 MULTICHANNEL POISSON ####
n <- p <- 10000
npeaks <- 1000
sim.pois.mc <- simul(method = 4, n = n, p = p, s = npeaks/p, seed = 2)

### 3.3.1 RELAXED GROUP L0GLM VIA CHOLESKY LLT ####
system.time(L0glm.chol.pois <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                         family = poisson(identity), 
                                         lambda = "gic",
                                         nonnegative = TRUE,
                                         grouped = TRUE,
                                         relaxed = TRUE,
                                         parallelize = TRUE,
                                         ncores = 4,
                                         solver = "eigen",
                                         method = 1L)) 
# openmp 6.8s, 99s with parallelize = FALSE, 38s with parallelize = TRUE and ncores = 4

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

sum(sim.pois.mc$beta_true!=0) # 73266
sum(L0glm.chol.pois$coefficients!=0) # 73844

plotmat(L0glm.chol.pois$coefficients,
      xlab="Time index", ylab="Channel", 
      main="nonnegative relaxed L0 group GLM\n(red = true peaks, cyan = estimated peaks)")
# abline(v=(1:nrow(sim.pois.mc$X))[as.vector(rowMax(sim.pois.mc$beta_true)!=0)], col="red")
# abline(v=(1:nrow(sim.pois.mc$X))[as.vector(rowMin(L0glm.chol.pois$coefficients)!=0)], col="cyan")
# TODO Fix bug Error: unable to find an inherited method for function ???rowQ??? for signature ???imat = "dgCMatrix", which = "numeric"???


### 3.3.2 COMPARE MULTICHANNEL POISSON METHODS ####
# L0glm fits the exact grouped Poisson model. Grouped ABESS and glmnet use a
# labelled Gaussian approximation with shared row weights 1/(rowMeans(y)+0.1).
# L0Learn uses independent channel-specific WLS; all are scored as Poisson.
system.time(comparison.mc.pois <- sparseglm(sim.pois.mc, family = poisson(identity),
  method = comparison_methods, nlambda = 10L, tune = "gic", nonnegative = TRUE,
  intercept = FALSE, approximation = "wls_grouped",
  support.size = seq_len(min(p, n - 1L, 2 * npeaks)),
  lambda = list(L0glm = L0glm.chol.pois$init.lambda),
  options = list(L0glm = list(relaxed = TRUE), abess = list(warm.start = TRUE))))

# plot(comparison.mc.pois)
# TODO: plot is currently very slow, presumably because it is plotting all coefs, including all the zero ones
# maybe add option to plot row means of true & estimated coefs for the multichannel case?


lasso.mc.pois <- coef(comparison.mc.pois)$lasso

x_vals <- as.vector(rowMeans(lasso.mc.pois))                      # estimated
y_vals <- as.vector(rowMeans(sim.pois.mc$beta_true))              # true
fit <- lm(y_vals ~ 0 + x_vals)
R2p <- round(100 * summary(fit)$r.squared, 2)
main_expr <- bquote(
  atop("true vs estimated coefficients (row means)",
       atop("independent-channel Poisson LASSO",
            R^2 == .(sprintf("%.2f", R2p)) * "%")))
plot(x_vals, y_vals, pch = 16, col = "steelblue",
     xlab = "estimated coefficients", 
     ylab = "true coefficients",
     main = main_expr)

# plotmat(lasso.mc.pois,
#       xlab="Time index", ylab="Channel", 
#       main="nonnegative independent-channel Poisson LASSO\n(red = true peaks, cyan = estimated peaks)")
# abline(v=(1:nrow(sim.pois.mc$X))[as.vector(rowMax(sim.pois.mc$beta_true)!=0)], col="red")
# abline(v=(1:nrow(sim.pois.mc$X))[as.vector(rowMin(lasso.mc.pois)!=0)], col="cyan")




### 3.3.4 BENCHMARKS TABLE ####

system.time(benchmarks_3_3 <- L0glm.metrics(comparison.mc.gaus)) # 4.78s
benchmarks_3_3 <- benchmarks_3_3[order(benchmarks_3_3[,"MCC"], decreasing = TRUE),]
benchmarks_3_3[,!colSums(is.na(benchmarks_3_3))]
#            method status       scope            approximation tune parameter selected_parameter tuning_score seconds
# 1                  L0glm     ok     grouped                     none  gic    lambda       2.044997e+01      5615963   21.12
# 4          relaxed_lasso     ok     grouped identity_poisson_row_wls  gic    lambda       1.242287e+00      6138687  805.60
# 6 adaptive_relaxed_lasso     ok     grouped identity_poisson_row_wls  gic    lambda       1.100621e+05      6191576  724.17
# 5         adaptive_lasso     ok     grouped identity_poisson_row_wls  gic    lambda       3.955426e+04      7348579  204.70
# 3                  lasso     ok     grouped identity_poisson_row_wls  gic    lambda       1.604476e-01      7413616  188.29
# 2                L0Learn     ok independent     identity_poisson_wls  gic    lambda       2.016690e-05     10236143  542.27
#   path_seconds cv_seconds note     p     n   m     df  deviance deviance_true   loglik MAD_coefs RMSE_coefs  R2_coefs
# 1        21.11          0      10000 10000 226  73844  807023.9      167913.4 -2052928 0.4451336   11.69795 0.9968233
# 4       805.57          0      10000 10000 226  94168  914123.0      289338.6 -2106477 0.9253209   18.32347 0.9921854
# 6       724.16          0      10000 10000 226 103963  766704.6      158833.0 -2032768 0.8376681   15.84627 0.9941978
# 5       204.67          0      10000 10000 226 156644  846383.2      217374.9 -2072607 1.2232663   19.87182 0.9911969
# 3       188.28          0      10000 10000 226 158306  877432.0      247170.1 -2088132 1.3251157   21.70045 0.9894781
# 2       542.24          0      10000 10000 226  49057 5934097.6     5471064.8 -4616465 3.4764744   77.19635 0.8702346
#   support_size_true support_size_estimated   TP   FP   TN FN precision recall  F1_score       MCC within_row_TP
# 1              1000                    955  944   11 8989 56 0.9884817  0.944 0.9657289 0.9623313         64452
# 4              1000                   1218  909  309 8691 91 0.7463054  0.909 0.8196573 0.8023123         60980
# 6              1000                   1462  967  495 8505 33 0.6614227  0.967 0.7855402 0.7743982         64291
# 5              1000                   1813  986  827 8173 14 0.5438500  0.986 0.7010309 0.6962279         69676
# 3              1000                   1840  980  860 8140 20 0.5326087  0.980 0.6901408 0.6847591         69464
# 2              1000                   3687 1000 2687 6313  0 0.2712232  1.000 0.4267122 0.4361743         33296
#   within_row_FP within_row_TN within_row_FN within_row_precision within_row_recall within_row_F1_score within_row_MCC
# 1          8559        135959          4374            0.8827711         0.9364484           0.9088179      0.8643182
# 4         16249        123498          4707            0.7895998         0.9283420           0.8533684      0.7819166
# 6         14723        133119          6409            0.8136659         0.9093494           0.8588509      0.7884679
# 5         22885        127728          2547            0.7527576         0.9647342           0.8456646      0.7719728
# 3         23673        126032          2311            0.7458260         0.9678022           0.8424372      0.7676552
# 2           646        152088         39970            0.9809675         0.4544536           0.6211477      0.5898789
#        gic     ebic     bic     aic     mbic
# 1  5615963  5466112 4785984 4253544  5941502
# 4  6138687  5947593 5080274 4401291  6553824
# 6  6191576  5980606 5023071 4273462  6649894
# 5  7348579  7030704 5587959 4458503  8039139
# 3  7413616  7092368 5634316 4492876  8111503
# 2 10236143 10136593 9684761 9331043 10452409

# L0glm.fit is the clear winner


## 3.4 MULTICHANNEL GAUSSIAN ####
n <- p <- 10000
npeaks <- 1000
sim.gaus.mc <- simul(method = 3, n = n, p = p, s = npeaks/p, seed = 1)

### 3.4.1 RELAXED GROUP L0GLM VIA CHOLESKY LLT ####

system.time(L0glm.chol.gaus <- L0glm.fit(sim.gaus.mc$X, sim.gaus.mc$y, 
                                    family = gaussian(identity), 
                                    lambda = "gic", 
                                    nonnegative = TRUE,
                                    grouped = TRUE,
                                    relaxed = TRUE,
                                    parallelize = TRUE,
                                    ncores = 4,
                                    solver = "eigen",
                                    method = 1L)) 
# 142s with parallelize = FALSE, 47s with parallelize = TRUE & ncores = 4

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

# plotmat(L0glm.chol.gaus$coefficients,
#       xlab="Time index", ylab="Channel", 
#       main="relaxed nonnegative group L0glm\n(red = true peaks, cyan = estimated peaks)")
# abline(v=(1:nrow(sim.gaus.mc$X))[as.vector(rowMax(sim.gaus.mc$beta_true)!=0)], col="red")
# abline(v=(1:nrow(sim.gaus.mc$X))[as.vector(rowMin(L0glm.chol.gaus$coefficients)!=0)], col="cyan")



### 3.4.2 COMPARE MULTICHANNEL GAUSSIAN METHODS ####
# Auto grouping uses native grouped fits where supported, independent paths
# otherwise; the table records the scope of every method.
system.time(comparison.mc.gaus <- sparseglm(sim.gaus.mc, family = gaussian(identity),
  method = comparison_methods, nlambda = 10L, tune = "gic", nonnegative = TRUE, intercept = FALSE,
  support.size = seq_len(min(p, n - 1L, 2 * npeaks)),
  lambda = list(L0glm = L0glm.chol.gaus$init.lambda),
  options = list(L0glm = list(relaxed = TRUE), abess = list(warm.start = TRUE))))
# 3053s=50 mins

lasso.mc.gaus <- coef(comparison.mc.gaus)$lasso
x_vals <- as.vector(rowMeans(lasso.mc.gaus))                      # estimated
y_vals <- as.vector(rowMeans(sim.gaus.mc$beta_true))               # true
fit <- lm(y_vals ~ 0 + x_vals)
R2p <- round(100 * summary(fit)$r.squared, 2)
main_expr <- bquote(
  atop("true vs estimated coefficients (row means)",
       atop("multivariate group LASSO",
            R^2 == .(sprintf("%.2f", R2p)) * "%")))
plot(x_vals, y_vals, pch = 16, col = "steelblue",
     xlab = "estimated coefficients", 
     ylab = "true coefficients",
     main = main_expr)

sum(sim.gaus.mc$beta_true!=0)        # 73070
sum(L0glm.chol.gaus$coefficients!=0) # 187365 - TODO: check how within-row sparsity can be improved

plotmat(lasso.mc.gaus,
      xlab="Time index", ylab="Channel", 
      main="nonnegative multivariate group LASSO\n(red = true peaks, cyan = estimated peaks)")
# abline(v=(1:nrow(sim.gaus.mc$X))[as.vector(rowMax(sim.gaus.mc$beta_true)!=0)], col="red")
# abline(v=(1:nrow(sim.gaus.mc$X))[as.vector(rowMin(lasso.mc.gaus)!=0)], col="cyan")

# TODO add plot method giving this sort of plot with multipanel with top showing correct solution
# and bottom plot indicating actual solution with TPs, FPs and FNs indicated in different colours



### 3.4.3 BENCHMARKS TABLE ####

system.time(benchmarks_3_4 <- L0glm.metrics(comparison.mc.gaus)) # 3.49s
benchmarks_3_4 <- benchmarks_3_4[order(benchmarks_3_4[,"MCC"], decreasing = TRUE),]
benchmarks_3_4[,!colSums(is.na(benchmarks_3_4))]
#                   method status       scope approximation tune parameter selected_parameter tuning_score seconds
# 1                  L0glm     ok     grouped          none  gic    lambda       2.044997e+01     16225216  180.80
# 5         adaptive_lasso     ok     grouped          none  gic    lambda       5.422378e+02     20728348   61.79
# 6 adaptive_relaxed_lasso     ok     grouped          none  gic    lambda       5.422378e+02     18419427  103.96
# 4          relaxed_lasso     ok     grouped          none  gic    lambda       3.970177e+00     16939870 1459.41
# 3                  lasso     ok     grouped          none  gic    lambda       1.426808e+00     20034062  178.42
# 2                L0Learn     ok independent          none  gic    lambda       3.888993e-06     16633518 1069.31
#   path_seconds cv_seconds note     p     n   m     df  deviance deviance_true   loglik MAD_coefs RMSE_coefs  R2_coefs
# 1       180.77          0      10000 10000 232 187597  29104315       6321224 -6194432 0.2830137   3.234918 0.9996738
# 5        61.76          0      10000 10000 232 179371 793931184     767698846 -8530109 1.7465920  23.622649 0.9826592
# 6       103.93          0      10000 10000 232 163182 275652065     251281660 -7541180 1.3009148  23.126024 0.9833310
# 4      1459.38          0      10000 10000 232 218471  30859110       9243342 -6236073 1.2541608  23.373356 0.9830773
# 3       178.41          0      10000 10000 232 315674  68711708      45861548 -6789270 1.7883043  27.525177 0.9769680
# 2      1069.29          0      10000 10000 232  88852 734111510     708697966 -7408249 3.0368554  63.963533 0.8796476
#   support_size_true support_size_estimated   TP   FP   TN  FN precision recall  F1_score        MCC within_row_TP
# 1              1000                    991  984    7 8993  16 0.9929364  0.984 0.9884480 0.98718335         67729
# 5              1000                    870  758  112 8888 242 0.8712644  0.758 0.8106952 0.79360839         51642
# 6              1000                    870  758  112 8888 242 0.8712644  0.758 0.8106952 0.79360839         50227
# 4              1000                   1561  958  603 8397  42 0.6137092  0.958 0.7481453 0.73646483         60286
# 3              1000                   1815  989  826 8174  11 0.5449036  0.989 0.7026643 0.69835078         68844
# 2              1000                   9987 1000 8987   13   0 0.1001302  1.000 0.1820333 0.01202632         24418
#   within_row_FP within_row_TN within_row_FN within_row_precision within_row_recall within_row_F1_score within_row_MCC
# 1        118827         37841          3891            0.3630492         0.9456716           0.5246731      0.2247461
# 5        108142         15045          1027            0.3231988         0.9805009           0.4861499      0.1631332
# 6         98067         25120          2442            0.3386988         0.9536350           0.4998632      0.1985036
# 4         86479         66109          9382            0.4107655         0.8653327           0.5570870      0.2924725
# 3        117341         40034          3229            0.3697613         0.9551982           0.5331413      0.2487036
# 2          5811        153119         48652            0.8077674         0.3341727           0.4727635      0.4106558
#        gic     ebic      bic      aic     mbic
# 1 16225216 15844529 14116696 12764058 17052232
# 5 20728348 20364353 18712285 17418959 21519099
# 6 18419427 18088284 16585323 15408725 19138809
# 4 16939870 16496530 14484338 12909088 17902992
# 3 20034062 19393470 16486005 14209888 21425700
# 2 16633518 16453212 15634855 14994202 17025219

# L0glm.fit is the clear winner
# TODO timings can still be improved drastically by porting to Rcpp/OpenMP, 
# maybe make L0nnls function for gaussian case, or port the whole L0glm.fit function to Rcpp/OpenMP?