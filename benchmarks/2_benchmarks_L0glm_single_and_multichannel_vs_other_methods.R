# BENCHMARKS L0GLM: SINGLE CHANNEL VS MULTICHANNEL DECONVOLUTION OF BLURRED SPIKE TRAIN ####
# Historical timings below are machine/build dependent.
# Historical numerical tables are retained for comparison. In particular, the
# Poisson table below used the former 0.9 coefficient-start mixture; removing
# that mixture changes one selected feature on this simulation.
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
library(ordinis)
library(L0Learn) # use at least version 2
library(abess)
library(bigGlm)

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
                   method = 1L)) # 4.19s
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
           method = 1L)) # 3.44s
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
                                      method = 1L)) # 9.76s
tune_zoomed_cv$lambda_best # 1.467145


#### 1.1.1 L0GLM FIT WITH CHOLESKY LLT SOLVER ####
system.time(L0glm.chol.pois <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                         family = poisson(identity), 
                                         lambda = "gic", # or tune_zoomed$lambda_best
                                         nonnegative = TRUE,
                                         relaxed = FALSE,
                                         # verbose = TRUE,
                                         solver = "eigen",
                                         method = 1L)) 
# 0.14s

#### 1.1.2 L0GLM FIT WITH CG SOLVER WITH Incomplete Cholesky LLT preconditioner ####
system.time(L0glm.CG.pois <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                       family = poisson(identity), 
                                       lambda = "gic", # or tune_zoomed$lambda_best
                                       nonnegative = TRUE,
                                       relaxed = FALSE,
                                       solver = "eigen",
                                       method = 3L)) 
# 0.28s

#### 1.1.3 L0GLM FIT WITH preconditioned LSCG SOLVER ####
system.time(L0glm.LSCG.pois <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                         family = poisson(identity), 
                                         lambda = "gic", # or tune_zoomed$lambda_best
                                         nonnegative = TRUE,
                                         relaxed = FALSE,
                                         solver = "eigen",
                                         method = 4L)) 
# 0.80s

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
# 5     L0glm_PCG 10000 10000 1 9349.100      2988.748 -41036.03 -0.02860991  80.28585  0.07161196  0.9928606 0.9093003
# 3 L0glm_cholLLT 10000 10000 1 9349.100      2988.750 -41036.03 -0.02861004  80.28585  0.07161198  0.9928607 0.9093003
# 2  L0glm_bigglm 10000 10000 1 9422.418      3019.605 -41072.69 -0.02894092  83.72523  0.07347632  0.9936964 0.9045904
# 1    L0glm_osqp 10000 10000 1 9462.048      3074.622 -41092.50 -0.02870220  80.66691  0.07307096  0.9914590 0.9085972
#   MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true support_size_estimated
# 4    265.3200     528.5209      0.3977446     0.9928605   0.9400566              1000                   1011
# 5    265.3201     528.5201      0.3977444     0.9928606   0.9400565              1000                   1011
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
                   method = 1L)) # 6.99s
tune_full$lambda_best # 20.44997 (best lambda tune)
tune_full$lambda_preset # 20.44997 (lambda preset)
min(tune_full$metrics$gic) # 82342.78

# TODO: check why coefs are going to zero before calculated lambda_max


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
                                      method = 1L)) # 7.34s

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
# 0.25s

#### 1.2.2 L0GLM FIT WITH PCG SOLVER WITH Incomplete Cholesky preconditioner ####
system.time(L0glm.CG.gaus <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                       family = gaussian(identity), 
                                       lambda = "gic", # or tune_zoomed$lambda_best
                                       nonnegative = TRUE,
                                       relaxed = FALSE,
                                       solver = "eigen",
                                       method = 3L)) 
# 0.31s


#### 1.2.3 L0GLM FIT WITH PLSCG SOLVER ####
system.time(L0glm.LSCG.gaus <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                         family = gaussian(identity), 
                                         lambda = "gic", # or tune_zoomed$lambda_best
                                         nonnegative = TRUE,
                                         relaxed = FALSE,
                                         solver = "eigen",
                                         method = 4L)) 
# 4.53s


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
# 6.84s



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
           ncores = 4)) # openmp 6.35s
tune_full$lambda_best # 0.879464 (best lambda tune)
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
           ncores = 4)) # openmp 7.19s
tune_zoomed$lambda_best # 4.080982 (best lambda tune)
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
# openmp 0.34s

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
                        ncores = 4)) # openmp 10.89s
tune_full$lambda_best # 1.036359 (best lambda tune)
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
                          ncores = 4)) # openmp 10.13s
tune_zoomed$lambda_best # 2.699964 (best lambda tune)
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
# openmp 0.5s



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
                        ncores = 4)) # openmp 11.36s
tune_full$lambda_best # 1.043608 (best lambda tune)
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
                          ncores = 4)) # openmp 10.62s
tune_zoomed$lambda_best # 1.682115 (best lambda tune)
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
                                    method = 1L))  # openmp 0.51s

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
                        solver = "eigen", method = 1L)) # 80s
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
# 5.40s

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
                                    method = 1L)) # 21.27s

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
# 22s

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
# bigGlm (LASSO), ordinis (MCP & SCAD), L0Learn & abess

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
# 2.6s

system.time(L0glm.chol <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                    family = poisson(identity), 
                                    lambda = "gic", # or tune.chol.sc$lambda_best,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "eigen",
                                    method = 1L)) 
# 0.2s

plot(x=L0glm.chol$coefficients, y=sim.pois$beta_true, pch=16, col='blue',
     main="L0 glm identity link Poisson with lambda tuned to optimise GIC",
     xlab="estimated coefficients", ylab="true coefficients")



### 3.1.2 LASSO USING bigGlm ####
system.time(bigGlm.sc.pois <- fit.GLM(sim = sim.pois, method = 2, 
                                      family = poisson(identity), plot = TRUE)) 
# 3.02s


### 3.1.3 ADAPTIVE LASSO USING bigGlm ####
system.time(bigGlm.al.sc.pois <- fit.GLM(sim = sim.pois, method = 3, 
                                         family = poisson(identity), plot = T)) 
# 2.88s
# TODO check what ridge penalty is used to calculate adaptive penalties
# maybe use the same one as used in 1st iteration of L0glm with lambda preset to optimise GIC?

# TODO add relaxed LASSO as additional fitting method in fit.GLM


### 3.1.4 L0LEARN ####
system.time(L0learn.sc.pois <- fit.GLM(sim = sim.pois, method = 4, 
                                       family = poisson(identity), plot = T))  
# 20.04s
# TODO fix warning below
# Warning message:
# `aes_string()` was deprecated in ggplot2 3.0.0.
# ??? Please use tidy evaluation idioms with `aes()`.
# ??? See also `vignette("ggplot2-in-packages")` for more information.
# ??? The deprecated feature was likely used in the L0Learn package.
# Please report the issue at <https://github.com/hazimehh/L0Learn/issues>.
# This warning is displayed once every 8 hours.
# Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated. 



### 3.1.5 L0LEARN AFTER LASSO ####
# better solution than regular L0Learn
system.time(L0learn.after.lasso.sc.pois <- fit.GLM(sim = sim.pois, method = 5, family = poisson(identity), plot = T)) 
# 2.84s


### 3.1.6 ABESS ####
# system.time(abess.sc.pois <- fit.GLM(sim = sim.pois, method = 6, 
#                                      family = poisson(identity), 
#                                      plot = T, IC = "gic"))  
# XXs
# TODO: takes a really long time to run - check if there is a bug in fit.GLM & abess is called correctly?
# TODO: NOTE: result is better if I do not use weights 
# CHECK WHY - maybe a bug in abess?

### 3.1.7 ABESS AFTER LASSO ####
# system.time(abess.after.lasso.sc.pois <- fit.GLM(sim = sim.pois, method = 7, 
#                                                  family = poisson(identity), 
#                                                  plot = T, 
#                                                  IC = "gic"))  
# XXs
# TODO: takes a really long time to run - check if there is a bug in fit.GLM & abess is called correctly?
# TODO: NOTE: result is better if I do not use weights - then cor 0.8849495 - CHECK
# CHECK WHY - maybe a bug in abess?


### 3.1.8 MCP USING ORDINIS ####
# system.time(ordinis.mcp.sc.pois <- fit.GLM(sim = sim.pois, method = 8, 
#                                            family = poisson(identity), plot = T,
#                                            IC = "gic")) 
# # XXs
# # TODO: takes a really long time to run - check if there is a bug in fit.GLM & ordinis is called correctly?
# 
# 
### 3.1.9 SCAD USING ORDINIS ####
# system.time(ordinis.scad.sc.pois <- fit.GLM(sim = sim.pois, method = 9, 
#                                             family = poisson(identity), 
#                                             plot = T, IC = "gic"))  
# # XXs
# # TODO: takes a really long time to run - check if there is a bug in fit.GLM & ordinis is called correctly?




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
  data.frame(method="nnalasso_bigGlm", do.call(cbind, L0glm:::get.metrics(beta = bigGlm.al.sc.pois,  # some metrics
                                                                  X = cbind(1, X), y = y, 
                                                                  family = poisson(identity), 
                                                                  beta_true = c(0, as.vector(sim.pois$beta_true))) )),
  data.frame(method="nnL0Learn_after_bigGlm_nnlasso", do.call(cbind, L0glm:::get.metrics(beta = L0learn.after.lasso.sc.pois,  # some metrics
                                                                                 X = X, y = y, 
                                                                                 family = poisson(identity),
                                                                                 beta_true = as.vector(sim.pois$beta_true)) )),
  # data.frame(method="nnabess", do.call(cbind, L0glm:::get.metrics(beta = abess.sc.pois,  # some metrics
  #                                                         X = X, y = y, 
  #                                                         family = poisson(identity),
  #                                                         beta_true = as.vector(sim.pois$beta_true)) )),
  # data.frame(method="nnabess_after_bigGlm_nnlasso", do.call(cbind, L0glm:::get.metrics(beta = abess.after.lasso.sc.pois,  # some metrics
  #                                                                              X = X, y = y, 
  #                                                                              family = poisson(identity),
  #                                                                              beta_true = as.vector(sim.pois$beta_true)) )),
  data.frame(method="nnL0Learn", do.call(cbind, L0glm:::get.metrics(beta = L0learn.sc.pois,  # some metrics
                                                            X = X, y = y, 
                                                            family = poisson(identity), 
                                                            beta_true = as.vector(sim.pois$beta_true)) )),
  data.frame(method="nnlasso_bigGlm", do.call(cbind, L0glm:::get.metrics(beta = bigGlm.sc.pois,  # some metrics
                                                                 X = cbind(1, X), y = y, 
                                                                 family = poisson(identity), 
                                                                 beta_true = c(0, as.vector(sim.pois$beta_true))) )) # ,
  # data.frame(method="nnMCP_ordinis", do.call(cbind, L0glm:::get.metrics(beta = ordinis.mcp.sc.pois,  # some metrics
  #                                                               X = cbind(1, X), y = y, 
  #                                                               family = poisson(identity),
  #                                                               beta_true = c(0, as.vector(sim.pois$beta_true))) )),
  # data.frame(method="nnSCAD_ordinis", do.call(cbind, L0glm:::get.metrics(beta = ordinis.scad.sc.pois,  # some metrics
  #                                                                X = cbind(1, X), y = y, 
  #                                                                family = poisson(identity), 
  #                                                                beta_true = c(0, as.vector(sim.pois$beta_true))) ))
  )

benchmarks_3_1 = benchmarks_3_1[order(benchmarks_3_1[,"deviance_true"]),] # sort by deviance relative to ground truth
benchmarks_3_1
#                           method     p     n m     deviance deviance_true        loglik     MD_coefs MAD_coefs sMAPE_coefs
# 1                          L0glm 10000 10000 1 9.607885e+03  2.673237e+03     -41167.85   -0.1240172  66.33935  0.06768539
# 5                 nnlasso_bigGlm 10001 10000 1 1.519571e+04  9.481692e+03     -43961.76   -4.4119940 227.28999  0.41171841
# 2                nnalasso_bigGlm 10001 10000 1 1.823586e+04  1.164163e+04     -45481.84   -1.1983929 115.70266  0.12522431
# 3 nnL0Learn_after_bigGlm_nnlasso 10000 10000 1 1.643789e+09  1.643964e+09 -821930748.36 -302.9962940 310.93014  0.19920985
# 4                      nnL0Learn 10000 10000 1 1.721955e+09  1.722136e+09 -861013922.83 -315.9794884 316.36254  0.19980247
#   bias_coefs  R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true
# 1  0.9729390 0.9195383  212.127126     450.6890     0.36722110     0.9729390   0.9394337              1000
# 5  1.1738751 0.8039161   -2.098986     822.1883     1.48255918     1.1738751   0.8054426              1000
# 2  0.9696198 0.8640589  156.721280     798.3723     0.76761330     0.9696198   0.8756978              1000
# 3  0.8593949 0.2625888 7933.847056    7933.8471     0.41970679     0.8593949   0.8842804              1000
# 4  0.9755619 0.2004882 1915.247549    1915.2475     0.02474037     0.9755619   1.0000000              1000
#   support_size_estimated estimated_over_true_support_size  TP   FP   TN  FN precision recall specificity  accuracy
# 1                    999                            0.999 845  154 8846 155 0.8458458  0.845   0.9828889 0.9691000
# 5                   2718                            2.718 956 1762 7239  44 0.3517292  0.956   0.8042440 0.8194181
# 2                   1199                            1.199 834  365 8636 166 0.6955796  0.834   0.9594490 0.9469053
# 3                      5                            0.005   5    0 9000 995 1.0000000  0.005   1.0000000 0.9005000
# 4                      1                            0.001   1    0 9000 999 1.0000000  0.001   1.0000000 0.9001000
#      F1_score        MCC          gic         ebic          bic          aic         mbic
# 1 0.845422711 0.82825713     102765.2     100738.0 9.153683e+04 8.433370e+04     107169.3
# 5 0.514254976 0.51264708     143507.1     137991.5 1.129572e+05 9.335952e+04     155489.3
# 2 0.758526603 0.73276395     115483.4     113050.3 1.020069e+05 9.336167e+04     120769.2
# 3 0.009950249 0.06709882 1643861599.0 1643861588.8 1.643862e+09 1.643862e+09 1643861621.0
# 4 0.001998002 0.03000150 1722027866.1 1722027864.1 1.722028e+09 1.722028e+09 1722027870.5

# TODO: should also include timings (but L0glm is much faster than all the rest)
# TODO: also include relaxed LASSO
# TODO: check why L0Learn results are so bad & check if abess can be made to run on this example
# -> L0glm beats all the other methods!


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
# 5.96s

system.time(L0glm.chol <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                    family = gaussian(identity), 
                                    lambda = "gic", # or tune.chol.sc$lambda_best,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "eigen",
                                    method = 1L)) 
# 0.27s for n=p=10000

plot(x=L0glm.chol$coefficients, y=sim.gaus$beta_true, pch=16, col='blue',
     main="L0 glm identity link gaussian with lambda tuned to optimise GIC",
     xlab="estimated coefficients", ylab="true coefficients")


### 3.2.2 LASSO USING bigGlm ####
system.time(bigGlm.sc.gaus <- fit.GLM(sim = sim.gaus, method = 2, 
                                      family = gaussian(identity), plot = TRUE)) 
# 3.00s 


### 3.2.3 ADAPTIVE LASSO USING bigGlm ####
system.time(bigGlm.al.sc.gaus <- fit.GLM(sim = sim.gaus, method = 3, 
                                         family = gaussian(identity), plot = T)) 
# 2.90s 


### 3.2.4 L0LEARN ####

system.time(L0learn.sc.gaus <- fit.GLM(sim = sim.gaus, method = 4, 
                                       family = gaussian(identity), plot = T))  
# 21.15s


### 3.2.5 L0LEARN AFTER LASSO ####
# Better than regular L0Learn & faster
system.time(L0learn.after.lasso.sc.gaus <- fit.GLM(sim = sim.gaus, method = 5, 
                                                   family = gaussian(identity), plot = T))  
# 2.90s 


### 3.2.6 ABESS ####
# system.time(abess.sc.gaus <- fit.GLM(sim = sim.gaus, method = 6, 
#                                      family = gaussian(identity), 
#                                      plot = T, IC = "gic"))  
# # TODO: takes a really long time to run - check if there is a bug in fit.GLM & abess is called correctly?
# # XXs
# 
# 
### 3.2.7 ABESS AFTER LASSO ####
# system.time(abess.after.lasso.sc.gaus <- fit.GLM(sim = sim.gaus, method = 7, 
#                                                  family = gaussian(identity), 
#                                                  plot = T, IC = "gic"))  
# # TODO: takes a really long time to run - check if there is a bug in fit.GLM & abess is called correctly?
# # XXs


### 3.2.8 MCP USING ORDINIS ####
# system.time(ordinis.mcp.sc.gaus <- fit.GLM(sim = sim.gaus, method = 8, 
#                                            family = gaussian(identity), 
#                                            plot = T, IC = "gic"))  
# # TODO: takes a really long time to run - check if there is a bug in fit.GLM & ordinis is called correctly?
# # 0.27s
# 
# 
### 3.2.9 SCAD USING ORDINIS ####
# system.time(ordinis.scad.sc.gaus <- fit.GLM(sim = sim.gaus, method = 9, 
#                                             family = gaussian(identity), 
#                                             plot = T, IC = "gic"))  
# # TODO: takes a really long time to run - check if there is a bug in fit.GLM & ordinis is called correctly?
# # 0.28s


### 3.2.10 BENCHMARKS TABLE ####

X <- sim.gaus$X
y <- sim.gaus$y
IC <- "gic"

benchmarks_3_2 = rbind(
  data.frame(method="L0glm", do.call(cbind, L0glm:::get.metrics(beta = L0glm.chol$coefficients,  # some metrics
                                                                     X = X, y = y, 
                                                                     family = gaussian(identity), 
                                                                     beta_true = sim.gaus$beta_true))),
  data.frame(method="nnalasso_bigGlm", do.call(cbind, L0glm:::get.metrics(beta = bigGlm.al.sc.gaus,  # some metrics
                                                                  X = cbind(1, X), y = y, 
                                                                  family = gaussian(identity), 
                                                                  beta_true = c(0, as.vector(sim.gaus$beta_true))) )),
  data.frame(method="nnL0Learn_after_bigGlm_nnlasso", do.call(cbind, L0glm:::get.metrics(beta = L0learn.after.lasso.sc.gaus,  # some metrics
                                                                                 X = X, y = y, 
                                                                                 family = gaussian(identity),
                                                                                 beta_true = as.vector(sim.gaus$beta_true)) )),
  # data.frame(method="abess_after_bigGlm_nnlasso", do.call(cbind, L0glm:::get.metrics(beta = abess.after.lasso.sc.gaus,  # some metrics
  #                                                                            X = X, y = y, 
  #                                                                            family = gaussian(identity),
  #                                                                            beta_true = as.vector(sim.gaus$beta_true)) )),
  data.frame(method="nnL0Learn", do.call(cbind, L0glm:::get.metrics(beta = L0learn.sc.gaus,  # some metrics
                                                            X = X, y = y, 
                                                            family = gaussian(identity), 
                                                            beta_true = as.vector(sim.gaus$beta_true)) )),
  # data.frame(method="abess", do.call(cbind, L0glm:::get.metrics(beta = abess.sc.gaus,  # some metrics
  #                                                       X = X, y = y, 
  #                                                       family = gaussian(identity), 
  #                                                       beta_true = as.vector(sim.gaus$beta_true)) )),
  data.frame(method="nnlasso_bigGlm", do.call(cbind, L0glm:::get.metrics(beta = bigGlm.sc.gaus,  # some metrics
                                                                 X = cbind(1, X), y = y, 
                                                                 family = gaussian(identity), 
                                                                 beta_true = c(0, as.vector(sim.gaus$beta_true))) )) # ,
  # data.frame(method="nnMCP_ordinis", do.call(cbind, L0glm:::get.metrics(beta = ordinis.mcp.sc.gaus,  # some metrics
  #                                                               X = cbind(1, X), y = y, 
  #                                                               family = gaussian(identity),
  #                                                               beta_true = c(0, as.vector(sim.gaus$beta_true))) )),
  # data.frame(method="nnSCAD_ordinis", do.call(cbind, L0glm:::get.metrics(beta = ordinis.scad.sc.gaus,  # some metrics
  #                                                                X = cbind(1, X), y = y, 
  #                                                                family = gaussian(identity), 
  #                                                                beta_true = c(0, as.vector(sim.gaus$beta_true))) ))
  )

benchmarks_3_2 = benchmarks_3_2[order(benchmarks_3_2[,"deviance_true"]),] # sort by deviance relative to ground truth
benchmarks_3_2
#                           method     p     n m     deviance deviance_true    loglik      MD_coefs MAD_coefs sMAPE_coefs
# 1                          L0glm 10000 10000 1 2.308267e+05  5.598291e+04 -29884.80  6.277559e-03  19.06869  0.03084378
# 5                 nnlasso_bigGlm 10001 10000 1 8.484388e+07  8.463086e+07 -59419.30 -9.890032e+00 211.32557  0.47852465
# 2                nnalasso_bigGlm 10001 10000 1 1.003780e+09  1.003973e+09 -71772.88 -2.884356e+01 193.22146  0.16273545
# 3 nnL0Learn_after_bigGlm_nnlasso 10000 10000 1 2.208901e+10  2.208724e+10 -87229.41 -1.328707e+02 285.10776  0.18966557
# 4                      nnL0Learn 10000 10000 1 2.215037e+10  2.214849e+10 -87243.28 -1.334420e+02 291.30557  0.19120841
#   bias_coefs  R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true
# 1  1.0033387 0.9839711    52.47761     130.5835      0.1981169     1.0033387   0.9875623              1000
# 5  0.9979581 0.8191980    14.81733     685.6544      1.5986526     0.9979581   0.8297161              1000
# 2  0.8581092 0.7700742   791.32478    1553.8156      0.8523784     0.8581092   0.8439967              1000
# 3  0.6833272 0.5126658  4817.56176    4817.6981      0.7889600     0.6833272   0.7501028              1000
# 4  0.6699811 0.4926224  5027.43283    5027.5701      0.8285608     0.6699811   0.7322389              1000
#   support_size_estimated estimated_over_true_support_size  TP   FP   TN  FN precision recall specificity  accuracy
# 1                   1042                            1.042 949   93 8907  51 0.9107486  0.949   0.9896667 0.9856000
# 5                   2876                            2.876 906 1970 7031  94 0.3150209  0.906   0.7811354 0.7936206
# 2                    701                            0.701 485  216 8785 515 0.6918688  0.485   0.9760027 0.9269073
# 3                    158                            0.158 114   44 8956 886 0.7215190  0.114   0.9951111 0.9070000
# 4                    157                            0.157 109   48 8952 891 0.6942675  0.109   0.9946667 0.9061000
#    F1_score       MCC       gic      ebic       bic       aic      mbic
# 1 0.9294809 0.9217078  81098.91  78982.36  69375.98  61855.59  85696.94
# 5 0.4674923 0.4554082 177673.79 171835.48 145336.75 124592.60 190356.89
# 2 0.5702528 0.5416869 157901.78 156477.21 150011.41 144949.75 160996.51
# 3 0.1968912 0.2624944 177710.37 177387.72 175923.27 174776.83 178411.32
# 4 0.1884183 0.2501767 177717.66 177397.04 175941.80 174802.57 178414.20

# TODO: should also include timings (but L0glm is much faster than all the rest)
# deviance relative to ground truth best for L0glm
# -> L0glm is the winner, both in terms of timings & quality of fit


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
# openmp 23s, 99s with parallelize = FALSE, 38s with parallelize = TRUE and ncores = 4

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


### 3.3.2 BEST SUBSET OF GROUP SELECTION (BSGS) USING ABESS ####
# system.time(bsgs.abess.pois <- fit.GLM.mc(sim = sim.pois.mc, method = 1, 
#                                              family = poisson(identity), 
#                                              plot = T, IC = "gic"))
# # XXs
# # TODO: takes a really long time to run - check if there is a bug in fit.GLM & abess is called correctly?
# # NOTE: this is now fit as if noise was Gaussian
# # TODO also compare what solution & timings one would get with
# # https://skscope.readthedocs.io/en/latest/gallery/GeneralizedLinearModels/poisson-identity-link.html
# # (but tuning optimal support size)
# # maybe implement this as an additional fitting method in 
# # fit.GLM.mc?




### 3.3.3 MULTIVARIATE GROUP LASSO USING bigGlm ####
system.time(bigGlm.mc.pois <- fit.GLM.mc(sim = sim.pois.mc, method = 2, 
                                    family = poisson(identity), plot = T, 
                                    IC = "gic"))
# 527.18s

x_vals <- as.vector(rowMeans(bigGlm.mc.pois))                      # estimated
y_vals <- as.vector(rowMeans(sim.pois.mc$beta_true))               # true
fit <- lm(y_vals ~ 0 + x_vals)
R2p <- round(100 * summary(fit)$r.squared, 2)
main_expr <- bquote(
  atop("true vs estimated coefficients (row means)",
       atop("multivariate weighted group LASSO",
            R^2 == .(sprintf("%.2f", R2p)) * "%")))
plot(x_vals, y_vals, pch = 16, col = "steelblue",
     xlab = "estimated coefficients", 
     ylab = "true coefficients",
     main = main_expr)

plotmat(bigGlm.mc.pois, 
      xlab="Time index", ylab="Channel", 
      main="bigGlm nonnegative multivariate group LASSO\n(red = true peaks, cyan = estimated peaks)")
# abline(v=(1:nrow(sim.pois.mc$X))[as.vector(rowMax(sim.pois.mc$beta_true)!=0)], col="red")
# abline(v=(1:nrow(sim.pois.mc$X))[as.vector(rowMin(bigGlm.mc.pois)!=0)], col="cyan")




### 3.3.4 BENCHMARKS TABLE ####

benchmarks_3_3 = rbind(
  data.frame(method="relaxed_group_L0glm_poisson", do.call(cbind, L0glm:::get.metrics(beta = as.matrix(L0glm.chol.pois$coefficients),  # TODO: make sure that metrics work with either dense or sparse coefs
                                                                        X = sim.pois.mc$X, y = sim.pois.mc$y,
                                                                        family = poisson(identity), 
                                                                        beta_true = sim.pois.mc$beta_true) ) ),
  # data.frame(method="BSGS_abess_poisson", do.call(cbind, L0glm:::get.metrics(beta = bsgs.abess.pois, 
  #                                                                    X = sim.pois.mc$X, y = sim.pois.mc$y, 
  #                                                                    family = poisson(identity), 
  #                                                                    beta_true = sim.pois.mc$beta_true) ) ),
  data.frame(method="multivariate_group_LASSO_poisson", do.call(cbind, L0glm:::get.metrics(beta = bigGlm.mc.pois,  
                                                                     X = sim.pois.mc$X, y = sim.pois.mc$y, 
                                                                     family = poisson(identity), 
                                                                     beta_true = sim.pois.mc$beta_true) ) ) 
  )
benchmarks_3_3 = benchmarks_3_3[order(benchmarks_3_3[,"deviance_true"]),] # sort by deviance relative to ground truth
benchmarks_3_3
#                             method     p     n   m    deviance deviance_true    loglik     MD_coefs MAD_coefs sMAPE_coefs
# 1      relaxed_group_L0glm_poisson 10000 10000 226    817245.2      192387.6  -2058038  0.001080747 0.4451336  0.02006709
# 2 multivariate_group_LASSO_poisson 10000 10000 226 125243175.7   124563443.9 -64271004 -1.717447776 3.6439305  0.08710050
#   bias_coefs  R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true
# 1  0.9970816 0.9968233    3.615294     10.04112      0.3754349     0.9970816   0.9979523              1000
# 2  1.0718804 0.9516476  -22.116457     64.63786      1.3956475     1.0718804   0.9590942              1000
#   support_size_estimated estimated_over_true_support_size  TP  FP   TN  FN precision recall specificity accuracy  F1_score
# 1                    955                            0.955 944  11 8989  56 0.9884817  0.944   0.9987778   0.9933 0.9657289
# 2                    928                            0.928 663 265 8735 337 0.7144397  0.663   0.9705556   0.9398 0.6877593
#         MCC       gic      ebic       bic       aic      mbic
# 1 0.9623313   5626184   5476334   4796205   4263765   5951723
# 2 0.6550590 130635920 130428138 129485073 128746791 131087312

# -> L0glm clearly the best - both timing-wise (ca. 100x faster than bigGlm) & in terms of performance


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

plotmat(L0glm.chol.gaus$coefficients,
      xlab="Time index", ylab="Channel", 
      main="relaxed nonnegative group L0glm\n(red = true peaks, cyan = estimated peaks)")
# abline(v=(1:nrow(sim.gaus.mc$X))[as.vector(rowMax(sim.gaus.mc$beta_true)!=0)], col="red")
# abline(v=(1:nrow(sim.gaus.mc$X))[as.vector(rowMin(L0glm.chol.gaus$coefficients)!=0)], col="cyan")



# ### 3.4.2 BEST SUBSET OF GROUP SELECTION (BSGS) USING ABESS ####
# system.time(bsgs.abess.gaus <- fit.GLM.mc(sim = sim.gaus.mc, method = 1, 
#                                           family = gaussian(identity), plot = T, 
#                                           IC = "gic")) 
# # 5s
# # TODO: takes a really long time to run - check if there is a bug in fit.GLM & abess is called correctly?
# # TODO post this as a reproducible example in https://github.com/abess-team/abess/issues/510
# # and ask if the solution quality could be improved here somehow



# bigGlm nonnegative mgaussian multivariate group LASSO fit
system.time(bigGlm.mc.gaus <- fit.GLM.mc(sim = sim.gaus.mc, method = 2, 
                                    family = gaussian(identity), plot = T, 
                                    IC = "gic"))
# 592.73s

x_vals <- as.vector(rowMeans(bigGlm.mc.gaus))                      # estimated
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
sum(L0glm.chol.gaus$coefficients!=0) # 187365 - TODO: check how overall sparsity can be made to match better

plotmat(bigGlm.mc.gaus,
      xlab="Time index", ylab="Channel", 
      main="nonnegative multivariate group LASSO\n(red = true peaks, cyan = estimated peaks)")
# abline(v=(1:nrow(sim.gaus.mc$X))[as.vector(rowMax(sim.gaus.mc$beta_true)!=0)], col="red")
# abline(v=(1:nrow(sim.gaus.mc$X))[as.vector(rowMin(bigGlm.mc.gaus)!=0)], col="cyan")




### 3.4.3 BENCHMARKS TABLE ####

benchmarks_3_4 = rbind(
  data.frame(method="relaxed_group_L0glm_gaussian", do.call(cbind, L0glm:::get.metrics(beta = as.matrix(L0glm.chol.gaus$coefficients),  # TODO: allow coefs to be either dense or sparse matrix in get.metrics
                                                                         X = sim.gaus.mc$X, y = sim.gaus.mc$y, 
                                                                         family = gaussian(identity),
                                                                         beta_true = sim.gaus.mc$beta_true)) ),
  # data.frame(method="BSGS_abess_gaussian", do.call(cbind, L0glm:::get.metrics(beta = bsgs.abess.gaus,  # some metrics
  #                                                                     X = sim.gaus.mc$X, y = sim.gaus.mc$y, 
  #                                                                     family = gaussian(identity),
  #                                                                     beta_true = sim.gaus.mc$beta_true)) ),
  data.frame(method="multivariate_group_LASSO_gaussian", do.call(cbind, L0glm:::get.metrics(beta = bigGlm.mc.gaus,  
                                                                                   X = sim.gaus.mc$X, y = sim.gaus.mc$y, 
                                                                                   family = gaussian(identity), 
                                                                                   beta_true = sim.gaus.mc$beta_true)) ) )

benchmarks_3_4 = benchmarks_3_4[order(benchmarks_3_4[,"deviance_true"]),] # sort by deviance relative to ground truth
benchmarks_3_4
#                              method     p     n   m   deviance deviance_true   loglik   MD_coefs MAD_coefs sMAPE_coefs
# 1      relaxed_group_L0glm_gaussian 10000 10000 232   29104315       6321224 -6194432  0.1235918 0.2830137   0.1145568
# 2 multivariate_group_LASSO_gaussian 10000 10000 232 3041642625    3013938623 -9761313 -0.8923490 3.0310019   0.1885916
# bias_coefs  R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true
# 1  0.9998545 0.9996738    1.827372     3.207318       1.361459     0.9998545   0.9997703              1000
# 2  1.0662827 0.9513915   -4.935836    25.841975       1.684906     1.0662827   0.9553134              1000
#   support_size_estimated estimated_over_true_support_size  TP  FP   TN  FN precision recall specificity accuracy  F1_score
# 1                    991                            0.991 984   7 8993  16 0.9929364  0.984   0.9992222   0.9977 0.9884480
# 2                   1101                            1.101 782 319 8681 218 0.7102634  0.782   0.9645556   0.9463 0.7444074
#         MCC      gic     ebic      bic      aic     mbic
# 1 0.9871833 16225216 15844529 14116696 12764058 17052232
# 2 0.7155156 24380782 23898699 21710663 19997753 25428070

# -> L0glm is the clear winner!

