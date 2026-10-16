# BENCHMARKS L0GLM: SINGLE CHANNEL VS MULTICHANNEL DECONVOLUTION OF BLURRED SPIKE TRAIN ####
# Historical timings below are machine/build dependent.

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
system.time(tune_full <- L0glm.tune(X=sim.pois$X, 
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
                   method = 1L)) # 3.11s
tune_full$lambda_best # 20.44997 (best lambda tune)
tune_full$lambda_preset # 20.44997 (lambda preset)
min(tune_full$metrics$gic) # 102766.8

## tuning plot (zoomed in)
# with lambda grid more zoomed in around lambda_preset preset
system.time(tune_zoomed <- L0glm.tune(X=sim.pois$X, 
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
           method = 1L)) # 3.37s
# here best lambda (dashed line) is very close to lambda preset expected 
# to optimise GIC (dotted line) so tuning in practice unnecessary
tune_zoomed$lambda_best # 11.49116 (best lambda tune)
tune_zoomed$lambda_preset # 20.44997 (lambda preset)
min(tune_zoomed$metrics$gic) # 102113.2

# Chosen lambda also matches expected lambda value that
# was expected to optimise GIC, so tuning is in fact not needed
L0glm:::preset.lambda(IC="gic", p, n) # 20.44997 (lambda preset)
tune_zoomed$lambda_preset # 20.44997 (lambda preset)
tune_zoomed$lambda_best # 11.49116 (best lambda tune)


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
# 0.33s

#### 1.1.3 L0GLM FIT WITH preconditioned LSCG SOLVER ####
system.time(L0glm.LSCG.pois <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                         family = poisson(identity), 
                                         lambda = "gic", # or tune_zoomed$lambda_best
                                         nonnegative = TRUE,
                                         relaxed = FALSE,
                                         solver = "eigen",
                                         method = 4L)) 
# 0.75s

#### 1.1.4 L0GLM FIT WITH OSQP SOLVER ####
system.time(L0glm.osqp.pois <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                    family = poisson(identity), 
                                    lambda = "gic", # or tune_zoomed$lambda_best
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "osqp")) 
# 4.73s

#### 1.1.5 L0GLM FIT WITH glmnet::bigGlm SOLVER ####
system.time(L0glm.glmnet.pois <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                      family = poisson(identity), 
                                      lambda = "gic", # or tune_zoomed$lambda_best
                                      nonnegative = TRUE,
                                      relaxed = FALSE,
                                      solver = "glmnet")) 
# 1.68s with n=p=10000




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
# 3 L0glm_cholLLT 10000 10000 1 9348.514      2987.828 -41035.74 -0.02862166  80.26889  0.07160071  0.9928768 0.9093261
# 5     L0glm_PCG 10000 10000 1 9348.517      2987.832 -41035.74 -0.02862145  80.26898  0.07160076  0.9928767 0.9093259
# 4   L0glm_PLSCG 10000 10000 1 9348.525      2987.833 -41035.74 -0.02861981  80.26933  0.07160092  0.9928763 0.9093256
# 2  L0glm_bigglm 10000 10000 1 9422.485      3019.769 -41072.72 -0.02894057  83.72414  0.07347620  0.9936962 0.9045913
# 1    L0glm_osqp 10000 10000 1 9462.469      3075.211 -41092.71 -0.02869495  80.65774  0.07307044  0.9914587 0.9086122
#   MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true support_size_estimated
# 3    262.8916     529.9963      0.3992165     0.9928768   0.9399447              1000                   1012
# 5    262.8916     529.9973      0.3992170     0.9928767   0.9399446              1000                   1012
# 4    262.8917     530.0008      0.3992185     0.9928763   0.9399442              1000                   1012
# 2    243.4808     581.1024      0.4125734     0.9936962   0.9283973              1000                   1015
# 1    270.6447     526.0850      0.4058344     0.9914587   0.9403085              1000                   1012
#   estimated_over_true_support_size  TP  FP   TN  FN precision recall specificity accuracy  F1_score       MCC      gic
# 3                            1.012 844 168 8832 156 0.8339921  0.844   0.9813333   0.9676 0.8389662 0.8209730 102766.8
# 5                            1.012 844 168 8832 156 0.8339921  0.844   0.9813333   0.9676 0.8389662 0.8209730 102766.8
# 4                            1.012 844 168 8832 156 0.8339921  0.844   0.9813333   0.9676 0.8389662 0.8209730 102766.8
# 2                            1.015 842 173 8827 158 0.8295567  0.842   0.9807778   0.9669 0.8357320 0.8173569 102902.2
# 1                            1.012 840 172 8828 160 0.8300395  0.840   0.9808889   0.9668 0.8349901 0.8165520 102880.8
#      mebic     ebic      bic      aic     mbic
# 3 123462.2 100713.2 91392.34 84095.47 107228.2
# 5 123462.2 100713.2 91392.34 84095.48 107228.2
# 4 123462.2 100713.2 91392.35 84095.48 107228.2
# 2 123658.9 100842.4 91493.94 84175.44 107376.8
# 1 123576.2 100827.2 91506.29 84209.43 107342.2

# cholesky LLT fastest & best solution quality (deviance_true)
# eigen solvers all equally good & cholesky the fastest
# osqp solver the worst


## 1.2. L0GLM FIT SINGLE CHANNEL WITH GAUSSIAN NOISE ####

## tuning plot (full range)
# full lambda grid range with warm.start = 0.01
# here preset lambda coincides with best lambda
system.time(tune_full <- L0glm.tune(X=sim.gaus$X, 
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
                   method = 1L)) # 3.61s
tune_full$lambda_best # 20.44997 (best lambda tune)
tune_full$lambda_preset # 20.44997 (lambda preset)
min(tune_full$metrics$gic) # 82394.15

# TODO: check why coefs are going to zero before calculated lambda_max


## tuning plot (zoomed in)
# tuning plots - with lambda grid more zoomed in around lambda_preset preset
system.time(tune_zoomed <- L0glm.tune(X=sim.gaus$X, 
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
                   method = 1L)) # 4.09s
# here best lambda (dashed line) is very close to lambda preset expected 
# to optimise GIC (dotted line) so tuning in fact not needed
tune_zoomed$lambda_best # 60.73341 (best lambda tune)
tune_zoomed$lambda_preset # 20.44997 (lambda preset)
min(tune_zoomed$metrics$gic) # 82200.08


# Chosen lambda here closely matches expected lambda value that
# was expected to optimise GIC, so tuning in fact not needed
L0glm:::preset.lambda(IC="gic", p, n) # 20.44997 (lambda preset)
tune_zoomed$lambda_best # 60.73341 (best lambda tune)
tune_zoomed$lambda_preset # 20.44997 (lambda preset)


#### 1.2.1 L0GLM FIT WITH Cholesky LLT SOLVER ####
system.time(L0glm.chol.gaus <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                         family = gaussian(identity), 
                                         lambda = "gic", # or tune_zoomed$lambda_best
                                         nonnegative = TRUE,
                                         relaxed = FALSE,
                                         solver = "eigen",
                                         method = 1L)) 
# 0.17s

#### 1.2.2 L0GLM FIT WITH PCG SOLVER WITH Incomplete Cholesky preconditioner ####
system.time(L0glm.CG.gaus <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                       family = gaussian(identity), 
                                       lambda = "gic", # or tune_zoomed$lambda_best
                                       nonnegative = TRUE,
                                       relaxed = FALSE,
                                       solver = "eigen",
                                       method = 3L)) 
# 0.20s


#### 1.2.3 L0GLM FIT WITH PLSCG SOLVER ####
system.time(L0glm.LSCG.gaus <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                         family = gaussian(identity), 
                                         lambda = "gic", # or tune_zoomed$lambda_best
                                         nonnegative = TRUE,
                                         relaxed = FALSE,
                                         solver = "eigen",
                                         method = 4L)) 
# 4.69s


#### 1.2.4 L0GLM FIT WITH OSQP SOLVER ####
system.time(L0glm.osqp.gaus <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                    family = gaussian(identity), 
                                    lambda = "gic", # or tune_zoomed$lambda_best
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "osqp")) 
# 38.73s


#### 1.2.5 L0GLM FIT WITH glmnet::bigGlm SOLVER ####
system.time(L0glm.glmnet.gaus <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                      family = gaussian(identity), 
                                      lambda = "gic", # or tune_zoomed$lambda_best
                                      nonnegative = TRUE,
                                      relaxed = FALSE,
                                      solver = "glmnet")) 
# 6.64s



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
#          method     p     n m deviance deviance_true    loglik     MD_coefs MAD_coefs sMAPE_coefs bias_coefs  R2_coefs
# 1    L0glm_osqp 10001 10000 1 206186.0      36010.44 -29321.35 -0.006233171  16.00687  0.02378203   1.000925 0.9898064
# 5     L0glm_PCG 10001 10000 1 254086.5      87732.86 -30365.83 -0.001600716  22.75013  0.03754649   1.016804 0.9860733
# 3 L0glm_cholLLT 10001 10000 1 254232.3      87854.06 -30368.70 -0.001636100  22.74396  0.03754825   1.016802 0.9860743
# 4   L0glm_PLSCG 10001 10000 1 271502.5     103340.99 -30697.32 -0.002801993  26.21674  0.03676606   1.025484 0.9707779
# 2  L0glm_bigglm 10001 10000 1 689201.0     533608.73 -35355.12  0.030608637 101.86554  0.10411557   1.041968 0.8918514
#   MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true support_size_estimated
# 1    33.75785     120.5395      0.1579752      1.000925   0.9914446              1000                   1037
# 5    61.49640     153.5181      0.2414602      1.016804   0.9895939              1000                   1058
# 3    61.49606     153.4597      0.2414769      1.016802   0.9895949              1000                   1058
# 4    56.40167     191.8364      0.2383150      1.025484   0.9738614              1000                   1056
# 2   102.88143     719.5249      0.6563000      1.041968   0.8982113              1000                   1239
#   estimated_over_true_support_size  TP  FP   TN  FN precision recall specificity accuracy  F1_score       MCC      gic
# 1                            1.037 963  74 8926  37 0.9286403  0.963   0.9917778   0.9889 0.9455081 0.9395231 79870.00
# 5                            1.058 940 118 8882  60 0.8884688  0.940   0.9868889   0.9822 0.9135083 0.9040423 82388.42
# 3                            1.058 940 118 8882  60 0.8884688  0.940   0.9868889   0.9822 0.9135083 0.9040423 82394.15
# 4                            1.056 942 114 8886  58 0.8920455  0.942   0.9873333   0.9828 0.9163424 0.9071830 83010.48
# 2                            1.239 886 353 8647 114 0.7150928  0.886   0.9607778   0.9533 0.7914247 0.7710423 96068.48
#      mebic     ebic      bic      aic      mbic
# 1 101097.3 77763.58 68203.04 60718.71  84445.97
# 5 104045.2 80239.38 70485.42 62849.67  87056.96
# 3 104050.9 80245.12 70491.16 62855.40  87062.70
# 4 104626.3 80865.51 71129.96 63508.63  87670.21
# 2 121426.7 93552.14 82131.07 73190.25 101534.95

# -> here osqp solver gave best solution & glmnet's bigGlm the worst
# Cholesky comes in at 3d place though & was by far the fastest
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
system.time(system.time(tune_full <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, 
           IC = "gic", 
           nlambda = 10,
           prop_range = 1, # to zoom in around lambda_preset preset
           warm.start = 0.01, #  set to 0 for no warm starts to not affect location of optimal lambda
           plot = TRUE, coefs.true = sim.pois.mc$beta_true, 
           relaxed = FALSE, family = poisson(identity), 
           nonnegative = TRUE, 
           parallelize = TRUE, 
           grouped = FALSE,
           solver = "eigen", method = 1L,
           ncores = 4))) # openmp 2.41s
tune_full$lambda_best # 0.9753673 (best lambda tune)
tune_full$lambda_preset # 11.35348 (lambda preset)
L0glm:::preset.lambda(IC="gic", p, n) # 11.35348 (lambda preset)

## tuning plot (zoomed in)
# here preset lambda would be ca 5x higher than best lambda
system.time(tune_zoomed <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, 
           IC = "gic", 
           nlambda = 10,
           prop_range = 0.2, # to zoom in around lambda_preset preset
           warm.start = 0.01, # set to 0 for no warm starts to not affect location of optimal lambda
           plot = TRUE, coefs.true = sim.pois.mc$beta_true, 
           relaxed = FALSE, family = poisson(identity), 
           nonnegative = TRUE, 
           parallelize = TRUE, 
           grouped = FALSE,
           solver = "eigen", method = 1L,
           ncores = 4)) # openmp 2.25s
tune_zoomed$lambda_best # 1.593531 (best lambda tune)
tune_zoomed$lambda_preset # 11.35348 (lambda preset)
L0glm:::preset.lambda(IC="gic", p, n) # 11.35348 (lambda preset)


#### Cholesky LLT ####
system.time(L0glm.chol <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                    family = poisson(identity), 
                                    lambda = "gic", # or tune_zoomed$lambda_best
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallelize = TRUE,
                                    warm.start = 1,
                                    grouped = FALSE,
                                    solver = "eigen",
                                    method = 1L)) 
# openmp 0.11s 

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
                        IC = "gic", 
                        nlambda = 10,
                        prop_range = 1, # to zoom in around lambda_preset preset
                        warm.start = 0.01, # here set to 0 for no warm starts to not affect location of optimal lambda
                        plot = TRUE, coefs.true = sim.pois.mc$beta_true, 
                        relaxed = FALSE, family = poisson(identity), 
                        nonnegative = TRUE, 
                        parallelize = TRUE, 
                        grouped = TRUE,
                        solver = "eigen", method = 1L,
                        ncores = 4)) # openmp 4.28s
tune_full$lambda_best # 0.9753657 (best lambda tune)
tune_full$lambda_preset # 11.35348 (lambda preset)
L0glm:::preset.lambda(IC="gic", p, n) # 11.35348 (lambda preset)

## tuning plot (zoomed in)
# here preset lambda would be ca 5x higher than best lambda
system.time(tune_zoomed <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, 
                          IC = "gic", 
                          nlambda = 10,
                          prop_range = 0.2, # to zoom in around lambda_preset preset
                          warm.start = 0.01, # set to 0 for no warm starts to not affect location of optimal lambda
                          plot = TRUE, coefs.true = sim.pois.mc$beta_true, 
                          relaxed = FALSE, family = poisson(identity), 
                          nonnegative = TRUE, 
                          parallelize = TRUE, 
                          grouped = TRUE,
                          solver = "eigen", method = 1L,
                          ncores = 4)) # openmp 3.19s
tune_zoomed$lambda_best # 2.60347 (best lambda tune)
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
# openmp 0.25s



#### OSQP ####
system.time(L0glm.osqp <- L0glm.fit(sim.pois.mc$X, sim.pois.mc$y, 
                                    family = poisson(identity), 
                                    lambda = "gic", # or tune_zoomed$lambda_best
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    parallelize = FALSE,
                                    grouped = TRUE,
                                    solver = "osqp"))
# 33.46s




## 2.1.3 GROUPED L0GLM FIT WITH relaxed = TRUE, grouped = TRUE (best) ####

## tuning plot (full range)
# here preset lambda would be ca 10x higher than best lambda
system.time(tune_full <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, 
                        IC = "gic", 
                        nlambda = 10,
                        prop_range = 1, # to zoom in around lambda_preset preset
                        warm.start = 0.01, # here set to 0 for no warm starts to not affect location of optimal lambda
                        plot = TRUE, coefs.true = sim.pois.mc$beta_true, 
                        relaxed = TRUE, family = poisson(identity), 
                        nonnegative = TRUE, 
                        parallelize = TRUE, 
                        grouped = TRUE,
                        solver = "eigen", method = 1L,
                        ncores = 4)) # openmp 4.56s
tune_full$lambda_best # 0.9753657 (best lambda tune)
tune_full$lambda_preset # 11.35348 (lambda preset)
L0glm:::preset.lambda(IC="gic", p, n) # 11.35348 (lambda preset)

## tuning plot (zoomed in)
# here preset lambda would be ca 10x higher than best lambda
system.time(tune_zoomed <- L0glm.tune(X=sim.pois.mc$X, y=sim.pois.mc$y, 
                          IC = "gic", 
                          nlambda = 10,
                          prop_range = 0.2, # to zoom in around lambda_preset preset
                          warm.start = 0.01, # set to 0 for no warm starts to not affect location of optimal lambda
                          plot = TRUE, coefs.true = sim.pois.mc$beta_true, 
                          relaxed = TRUE, family = poisson(identity), 
                          nonnegative = TRUE, 
                          parallelize = TRUE, 
                          grouped = TRUE,
                          solver = "eigen", method = 1L,
                          ncores = 4)) # openmp 3.53s
tune_zoomed$lambda_best # 0.9755287 (best lambda tune)
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
                                    method = 1L))  # openmp 0.39s

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
system.time(tune_full <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, 
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
                        solver = "eigen", method = 1L)) # 39.28s
tune_full$lambda_best # 11.35348
tune_full$lambda_preset # 11.35348 
L0glm:::preset.lambda(IC="gic", p, n) # 11.35348

## tuning plot (zoomed in)
# here best lambda is identical to preset lambda
system.time(tune_zoomed <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, 
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
                          solver = "eigen", method = 1L)) # 51.81s
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
# 3.92s

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
system.time(tune_full <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, 
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
                        solver = "eigen", method = 1L)) # 131.19s
tune_full$lambda_best # 11.35348
tune_full$lambda_preset # 11.35348 
L0glm:::preset.lambda(IC="gic", p, n) # 11.35348

## tuning plot (zoomed in)
# here best lambda is identical to preset lambda
system.time(tune_zoomed <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, 
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
                          solver = "eigen", method = 1L)) # 229s
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
                                    method = 1L)) # 19s

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
# 43s



### 2.2.3 GROUPED L0GLM FIT WITH relaxed = TRUE, grouped = TRUE (best) ####

## tuning plot (full range)
# here best lambda is identical to preset lambda
system.time(tune_full <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, 
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
                        solver = "eigen", method = 1L)) # 146s
tune_full$lambda_best # 11.35348 (best lambda tune)
tune_full$lambda_preset # 11.35348 (lambda preset)
L0glm:::preset.lambda(IC="gic", p, n) # 11.35348 (lambda preset)

## tuning plot (zoomed in)
# here best lambda is again very close to preset lambda
system.time(tune_zoomed <- L0glm.tune(X=sim.gaus.mc$X, y=sim.gaus.mc$y, 
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
                          solver = "eigen", method = 1L)) # 239s
tune_zoomed$lambda_best # 4.82 (best lambda tune)
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
# 23s






# 3. BENCHMARKS AGAINST OTHER METHODS ####
# Take our best working solver & method and compare it against others: 
# glmnet (LASSO), ordinis (MCP & SCAD), L0Learn & abess

## 3.1 SINGLE CHANNEL POISSON ####
n <- p <- 10000
npeaks <- 1000
sim.pois <- simul(method = 2, n = n, p = p, s = npeaks/p, seed = 2)

### 3.1.1 L0GLM VIA CHOLESKY LLT ####
system.time(tune.chol.sc <- L0glm.tune(X=sim.pois$X, y=sim.pois$y, IC = "gic", 
                                       plot = FALSE, coefs.true = sim.pois$beta_true, 
                                       relaxed = FALSE, family = poisson(identity), 
                                       nonnegative = TRUE, 
                                       solver = "eigen", method = 1L))
# 3.53s

system.time(L0glm.chol <- L0glm.fit(sim.pois$X, sim.pois$y, 
                                    family = poisson(identity), 
                                    lambda = "gic", # or tune.chol.sc$lambda_best,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "eigen",
                                    method = 1L)) 
# 0.34s

# Rcpp/OpenMP version
system.time(L0glm.chol <- L0nnpois(sim.pois$X, sim.pois$y, 
                                    family = poisson(identity), 
                                    lambda = "gic", # or tune.chol.sc$lambda_best,
                                    relaxed = FALSE,
                                    solver = "eigen",
                                    method = 1L))
# 0.20s

plot(x=L0glm.chol$coefficients, y=sim.pois$beta_true, pch=16, col='blue',
     main="L0 glm identity link Poisson with lambda tuned to optimise GIC",
     xlab="estimated coefficients", ylab="true coefficients")



### 3.1.2 LASSO USING GLMNET ####
system.time(glmnet.sc.pois <- fit.GLM(sim = sim.pois, method = 2, 
                                      family = poisson(identity), plot = TRUE)) 
# 2.66s


### 3.1.3 ADAPTIVE LASSO USING GLMNET ####
system.time(glmnet.al.sc.pois <- fit.GLM(sim = sim.pois, method = 3, 
                                         family = poisson(identity), plot = T)) 
# 2.48s
# TODO check what ridge penalty is used to calculate adaptive penalties
# maybe use the same one as used in 1st iteration of L0glm with lambda preset to optimise GIC?

# TODO add relaxed LASSO as additional fitting method in fit.GLM


### 3.1.4 L0LEARN ####
system.time(L0learn.sc.pois <- fit.GLM(sim = sim.pois, method = 4, 
                                       family = poisson(identity), plot = T))  
# 10.58s
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
# 3.72s


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
  data.frame(method="nnalasso_glmnet", do.call(cbind, L0glm:::get.metrics(beta = glmnet.al.sc.pois,  # some metrics
                                                                  X = cbind(1, X), y = y, 
                                                                  family = poisson(identity), 
                                                                  beta_true = c(0, as.vector(sim.pois$beta_true))) )),
  data.frame(method="nnL0Learn_after_glmnet_nnlasso", do.call(cbind, L0glm:::get.metrics(beta = L0learn.after.lasso.sc.pois,  # some metrics
                                                                                 X = X, y = y, 
                                                                                 family = poisson(identity),
                                                                                 beta_true = as.vector(sim.pois$beta_true)) )),
  # data.frame(method="nnabess", do.call(cbind, L0glm:::get.metrics(beta = abess.sc.pois,  # some metrics
  #                                                         X = X, y = y, 
  #                                                         family = poisson(identity),
  #                                                         beta_true = as.vector(sim.pois$beta_true)) )),
  # data.frame(method="nnabess_after_glmnet_nnlasso", do.call(cbind, L0glm:::get.metrics(beta = abess.after.lasso.sc.pois,  # some metrics
  #                                                                              X = X, y = y, 
  #                                                                              family = poisson(identity),
  #                                                                              beta_true = as.vector(sim.pois$beta_true)) )),
  data.frame(method="nnL0Learn", do.call(cbind, L0glm:::get.metrics(beta = L0learn.sc.pois,  # some metrics
                                                            X = X, y = y, 
                                                            family = poisson(identity), 
                                                            beta_true = as.vector(sim.pois$beta_true)) )),
  data.frame(method="nnlasso_glmnet", do.call(cbind, L0glm:::get.metrics(beta = glmnet.sc.pois,  # some metrics
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
#                           method     p     n m     deviance deviance_true        loglik     MD_coefs MAD_coefs
# 1                          L0glm 10000 10000 1 9.657865e+03  2.726736e+03     -41192.84   -0.1236613  66.51161
# 5                 nnlasso_glmnet 10001 10000 1 1.519571e+04  9.481692e+03     -43961.76   -4.4119940 227.28999
# 2                nnalasso_glmnet 10001 10000 1 1.823586e+04  1.164163e+04     -45481.84   -1.1983929 115.70266
# 3 nnL0Learn_after_glmnet_nnlasso 10000 10000 1 1.643789e+09  1.643964e+09 -821930748.36 -302.9962940 310.93014
# 4                      nnL0Learn 10000 10000 1 1.721955e+09  1.722136e+09 -861013922.83 -315.9794884 316.36254
#   sMAPE_coefs bias_coefs  R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true
# 1  0.06719598  0.9724896 0.9190199  216.567653     451.3192     0.36012056     0.9724896   0.9391526              1000
# 5  0.41171841  1.1738751 0.8039161   -2.098986     822.1883     1.48255918     1.1738751   0.8054426              1000
# 2  0.12522431  0.9696198 0.8640589  156.721280     798.3723     0.76761330     0.9696198   0.8756978              1000
# 3  0.19920985  0.8593949 0.2625888 7933.847056    7933.8471     0.41970679     0.8593949   0.8842804              1000
# 4  0.19980247  0.9755619 0.2004882 1915.247549    1915.2475     0.02474037     0.9755619   1.0000000              1000
#   support_size_estimated estimated_over_true_support_size  TP   FP   TN  FN precision recall specificity  accuracy
# 1                    994                            0.994 843  151 8849 157 0.8480885  0.843   0.9832222 0.9692000
# 5                   2718                            2.718 956 1762 7239  44 0.3517292  0.956   0.8042440 0.8194181
# 2                   1199                            1.199 834  365 8636 166 0.6955796  0.834   0.9594490 0.9469053
# 3                      5                            0.005   5    0 9000 995 1.0000000  0.005   1.0000000 0.9005000
# 4                      1                            0.001   1    0 9000 999 1.0000000  0.001   1.0000000 0.9001000
#      F1_score          gic         ebic          bic          aic         mbic
# 1 0.845536610     102712.9     100695.8 9.154075e+04 8.437368e+04     107095.0
# 5 0.514254976     143507.1     137991.5 1.129572e+05 9.335952e+04     155489.3
# 2 0.758526603     115483.4     113050.3 1.020069e+05 9.336167e+04     120769.2
# 3 0.009950249 1643861599.0 1643861588.8 1.643862e+09 1.643862e+09 1643861621.0
# 4 0.001998002 1722027866.1 1722027864.1 1.722028e+09 1.722028e+09 1722027870.5

# TODO: should also include timings (but L0glm is much faster than all the rest)
# TODO: also include relaxed LASSO
# TODO: check why L0Learn results are so bad & check if abess can be made to run on this example
# -> L0glm beats all the other methods!


## 3.2 SINGLE CHANNEL GAUSSIAN ####
n <- p <- 10000
npeaks <- 1000
sim.gaus <- simul(method = 1, n = n, p = p, s = npeaks/p, seed = 2)

### 3.2.1 L0GLM VIA CHOLESKY LLT ####

system.time(tune.chol.sc <- L0glm.tune(X=sim.gaus$X, y=sim.gaus$y, IC = "gic", 
                                       plot = FALSE, coefs.true = sim.gaus$beta_true, 
                                       relaxed = FALSE, family = gaussian(identity), 
                                       nonnegative = TRUE, 
                                       solver = "eigen", method = 1L))
# 1.76s

system.time(L0glm.chol <- L0glm.fit(sim.gaus$X, sim.gaus$y, 
                                    family = gaussian(identity), 
                                    lambda = "gic", # or tune.chol.sc$lambda_best,
                                    nonnegative = TRUE,
                                    relaxed = FALSE,
                                    solver = "eigen",
                                    method = 1L)) 
# 0.18s for n=p=10000

plot(x=L0glm.chol$coefficients, y=sim.gaus$beta_true, pch=16, col='blue',
     main="L0 glm identity link gaussian with lambda tuned to optimise GIC",
     xlab="estimated coefficients", ylab="true coefficients")


### 3.2.2 LASSO USING GLMNET ####
system.time(glmnet.sc.gaus <- fit.GLM(sim = sim.gaus, method = 2, 
                                      family = gaussian(identity), plot = TRUE)) 
# 1.81s 


### 3.2.3 ADAPTIVE LASSO USING GLMNET ####
system.time(glmnet.al.sc.gaus <- fit.GLM(sim = sim.gaus, method = 3, 
                                         family = gaussian(identity), plot = T)) 
# 1.70s 


### 3.2.4 L0LEARN ####

system.time(L0learn.sc.gaus <- fit.GLM(sim = sim.gaus, method = 4, 
                                       family = gaussian(identity), plot = T))  
# 7.56s


### 3.2.5 L0LEARN AFTER LASSO ####
# Better than regular L0Learn & faster
system.time(L0learn.after.lasso.sc.gaus <- fit.GLM(sim = sim.gaus, method = 5, 
                                                   family = gaussian(identity), plot = T))  
# 2.84s 


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
  data.frame(method="nnalasso_glmnet", do.call(cbind, L0glm:::get.metrics(beta = glmnet.al.sc.gaus,  # some metrics
                                                                  X = cbind(1, X), y = y, 
                                                                  family = gaussian(identity), 
                                                                  beta_true = c(0, as.vector(sim.gaus$beta_true))) )),
  data.frame(method="nnL0Learn_after_glmnet_nnlasso", do.call(cbind, L0glm:::get.metrics(beta = L0learn.after.lasso.sc.gaus,  # some metrics
                                                                                 X = X, y = y, 
                                                                                 family = gaussian(identity),
                                                                                 beta_true = as.vector(sim.gaus$beta_true)) )),
  # data.frame(method="abess_after_glmnet_nnlasso", do.call(cbind, L0glm:::get.metrics(beta = abess.after.lasso.sc.gaus,  # some metrics
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
  data.frame(method="nnlasso_glmnet", do.call(cbind, L0glm:::get.metrics(beta = glmnet.sc.gaus,  # some metrics
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
# 1                          L0glm 10001 10000 1 2.297457e+05  5.547342e+04 -29862.33  4.274065e-03  19.68919  0.03188896
# 5                 nnlasso_glmnet 10002 10000 1 8.484388e+07  8.463086e+07 -59420.30 -9.890032e+00 211.32557  0.47852465
# 2                nnalasso_glmnet 10002 10000 1 1.003780e+09  1.003973e+09 -71773.88 -2.884356e+01 193.22146  0.16273545
# 3 nnL0Learn_after_glmnet_nnlasso 10001 10000 1 2.208901e+10  2.208724e+10 -87230.41 -1.328707e+02 285.10776  0.18966557
# 4                      nnL0Learn 10001 10000 1 2.215037e+10  2.214849e+10 -87244.28 -1.334420e+02 291.30557  0.19120841
#   bias_coefs  R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true
# 1  1.0029614 0.9833563    57.57769     131.7825      0.2008554     1.0029614   0.9871797              1000
# 5  0.9979581 0.8191980    14.81733     685.6544      1.5986526     0.9979581   0.8297161              1000
# 2  0.8581092 0.7700742   791.32478    1553.8156      0.8523784     0.8581092   0.8439967              1000
# 3  0.6833272 0.5126658  4817.56176    4817.6981      0.7889600     0.6833272   0.7501028              1000
# 4  0.6699811 0.4926224  5027.43283    5027.5701      0.8285608     0.6699811   0.7322389              1000
#   support_size_estimated estimated_over_true_support_size  TP   FP   TN  FN precision recall specificity  accuracy
# 1                   1040                            1.040 945   95 8905  55 0.9086538  0.945   0.9894444 0.9850000
# 5                   2876                            2.876 906 1970 7031  94 0.3150209  0.906   0.7811354 0.7936206
# 2                    701                            0.701 485  216 8785 515 0.6918688  0.485   0.9760027 0.9269073
# 3                    158                            0.158 114   44 8956 886 0.7215190  0.114   0.9951111 0.9070000
# 4                    157                            0.157 109   48 8952 891 0.6942675  0.109   0.9946667 0.9061000
#    F1_score      gic      ebic       bic       aic      mbic
# 1 0.9264706  81013.3  78900.79  69312.61  61806.65  85602.49
# 5 0.4674923 177676.4 171838.05 145338.75 124594.60 190359.46
# 2 0.5702528 157903.9 156479.35 150013.41 144951.75 160998.65
# 3 0.1968912 177712.4 177389.75 175925.27 174778.83 178413.35
# 4 0.1884183 177719.7 177399.07 175943.80 174804.57 178416.23

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
# 99s with parallelize = FALSE, 38s with parallelize = TRUE and ncores = 4

# Rcpp/OpenMP version: 5 to 18x faster!
system.time(L0glm.chol.pois <- L0nnpois(sim.pois.mc$X, sim.pois.mc$y, 
                                         lambda = "gic",
                                         grouped = TRUE,
                                         relaxed = TRUE,
                                         parallelize = TRUE,
                                         ncores = 4)) 
# 19s with parallelize = FALSE, 8s with parallelize = TRUE and ncores = 4, 5.5s with ncores = 8 to 12

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
sum(L0glm.chol.pois$coefficients!=0) # 70891

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




### 3.3.3 MULTIVARIATE GROUP LASSO USING GLMNET ####
system.time(glmnet.mc.pois <- fit.GLM.mc(sim = sim.pois.mc, method = 2, 
                                    family = poisson(identity), plot = T, 
                                    IC = "gic"))
# 643s

x_vals <- as.vector(rowMeans(glmnet.mc.pois))                      # estimated
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
  # data.frame(method="BSGS_abess_poisson", do.call(cbind, L0glm:::get.metrics(beta = bsgs.abess.pois, 
  #                                                                    X = sim.pois.mc$X, y = sim.pois.mc$y, 
  #                                                                    family = poisson(identity), 
  #                                                                    beta_true = sim.pois.mc$beta_true) ) ),
  data.frame(method="multivariate_group_LASSO_poisson", do.call(cbind, L0glm:::get.metrics(beta = glmnet.mc.pois,  
                                                                     X = sim.pois.mc$X, y = sim.pois.mc$y, 
                                                                     family = poisson(identity), 
                                                                     beta_true = sim.pois.mc$beta_true) ) ) 
  )
benchmarks_3_3 = benchmarks_3_3[order(benchmarks_3_3[,"deviance_true"]),] # sort by deviance relative to ground truth
benchmarks_3_3
#                             method     p     n   m  deviance deviance_true    loglik     MD_coefs MAD_coefs sMAPE_coefs
# 1      relaxed_group_L0glm_poisson 10000 10000 226   1068548      498863.1  -2183690  0.000832458 0.4456107  0.01976052
# 2 multivariate_group_LASSO_poisson 10000 10000 226 125243176   124563443.9 -64271004 -1.717447776 3.6439305  0.08710050
#   bias_coefs  R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true
# 1  0.9970806 0.9968227    3.776191     10.45638      0.3444271     0.9970806   0.9979517              1000
# 2  1.0718804 0.9516476  -22.116457     64.63786      1.3956475     1.0718804   0.9590942              1000
#   support_size_estimated estimated_over_true_support_size  TP  FP   TN  FN precision recall specificity accuracy
# 1                    955                            0.955 944  11 8989  56 0.9884817  0.944   0.9987778   0.9933
# 2                    928                            0.928 663 265 8735 337 0.7144397  0.663   0.9705556   0.9398
#    F1_score       gic      ebic       bic       aic      mbic
# 1 0.9657289   5817099   5673241   5020310   4509162   6129619
# 2 0.6877593 130635920 130428138 129485073 128746791 131087312

# -> L0glm clearly the best - both timing-wise (ca. 100x faster than glmnet) & in terms of performance


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
# 142s with parallelize = FALSE, 35s with parallelize = TRUE & ncores = 4
       
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



# glmnet nonnegative mgaussian multivariate group LASSO fit
system.time(glmnet.mc.gaus <- fit.GLM.mc(sim = sim.gaus.mc, method = 2, 
                                    family = gaussian(identity), plot = T, 
                                    IC = "gic"))
# 693s

x_vals <- as.vector(rowMeans(glmnet.mc.gaus))                      # estimated
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
sum(L0glm.chol.gaus$coefficients!=0) # 186610 - TODO: check how overall sparsity can be made to match better

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
  # data.frame(method="BSGS_abess_gaussian", do.call(cbind, L0glm:::get.metrics(beta = bsgs.abess.gaus,  # some metrics
  #                                                                     X = sim.gaus.mc$X, y = sim.gaus.mc$y, 
  #                                                                     family = gaussian(identity),
  #                                                                     beta_true = sim.gaus.mc$beta_true)) ),
  data.frame(method="multivariate_group_LASSO_gaussian", do.call(cbind, L0glm:::get.metrics(beta = glmnet.mc.gaus,  
                                                                                   X = sim.gaus.mc$X, y = sim.gaus.mc$y, 
                                                                                   family = gaussian(identity), 
                                                                                   beta_true = sim.gaus.mc$beta_true)) ) )

benchmarks_3_4 = benchmarks_3_4[order(benchmarks_3_4[,"deviance_true"]),] # sort by deviance relative to ground truth
benchmarks_3_4
#                              method     p     n   m   deviance deviance_true    loglik   MD_coefs MAD_coefs sMAPE_coefs
# 1      relaxed_group_L0glm_gaussian 10001 10000 232   29114886       6329184  -6226370  0.1236105  0.281523   0.1139298
# 2 multivariate_group_LASSO_gaussian 10001 10000 232 3041642625    3013938623 -11619098 -0.8923490  3.031002   0.1885916
#   bias_coefs  R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true
# 1  0.9998454 0.9996843    1.828325     3.208435       1.358797     0.9998454   0.9997756              1000
# 2  1.0662827 0.9513915   -4.935836    25.841975       1.684906     1.0662827   0.9553134              1000
#   support_size_estimated estimated_over_true_support_size  TP  FP   TN  FN precision recall specificity accuracy
# 1                    989                            0.989 984   5 8995  16 0.9949444  0.984   0.9994444   0.9979
# 2                   1101                            1.101 782 319 8681 218 0.7102634  0.782   0.9645556   0.9463
#    F1_score      gic     ebic      bic      aic     mbic
# 1 0.9894419 16268970 15890279 14171491 12825962 17091634
# 2 0.7444074 28091680 27610060 25424104 23712860 29137945

# -> L0glm is the clear winner!

