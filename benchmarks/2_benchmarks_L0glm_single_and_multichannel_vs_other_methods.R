# BENCHMARKS L0GLM: SINGLE CHANNEL VS MULTICHANNEL DECONVOLUTION OF BLURRED SPIKE TRAIN ####
# Historical timings below are machine/build dependent.
# Historical numerical tables are retained for comparison. The Poisson table
# is reproduced by maxit = 25, tol = 1e-10 (the restored fitting defaults).
# More iterations can improve GIC while worsening recovery against the truth.
# Current l0glmTune fits each lambda independently and verifies its numerical
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


# 1. BENCHMARKS OF L0glm AGAINST OTHER METHODS FOR SINGLE & MULTICHANNEL DECONVOLUTION ####

# We take our best working & fastest Eigen Cholesky LLT solver and compare it against others: 
# L0Learn & glmnet (LASSO variants)
# abess, MCP and SCAD were unacceptably slow for these n = p = 10000 deconvolution
# problems, so we omit them here (and for smaller problems they all give a poor solution).

comparison_methods <- c("L0glm", "L0Learn", "lasso", "relaxed_lasso",
                        "adaptive_lasso", "adaptive_relaxed_lasso")


## 1.1 SINGLE CHANNEL POISSON ####
n <- p <- 10000
npeaks <- 1000
sim.pois <- simul(method = 2, n = n, p = p, s = npeaks/p, seed = 2)

system.time(L0glm.chol <- l0glmFit(sim.pois$X, sim.pois$y,
                                         family = poisson(identity), 
                                         lambda = "gic", 
                                         nonnegative = TRUE,
                                         solver = "eigen",
                                         method = 1L)) 
# 0.18s

# The L0glm candidate uses the fastest Eigen Cholesky LLT solver and a lambda GIC preset; 
# other methods are rescored with the common Poisson GIC. WLS is labelled in the table.
# glmnet uses 10 lambda values. 
system.time(comparison.sc.pois <- sparseGLM(sim.pois, family = poisson(identity),
  method = comparison_methods, nlambda = 10L,
  tune = "gic", nonnegative = TRUE, intercept = FALSE, approximation = "wls",
  support.size = seq_len(min(p, n - 1L, 2 * npeaks)),
  lambda = list(L0glm = L0glm.chol$init.lambda),
  options = list(abess = list(warm.start = TRUE)))) # 2.69s
plot(comparison.sc.pois)

benchmarks_1_1 <- l0glmMetrics(comparison.sc.pois)
benchmarks_1_1 <- benchmarks_1_1[order(benchmarks_1_1[,"MCC"], decreasing = TRUE),]
benchmarks_1_1[,!colSums(is.na(benchmarks_1_1))]
#                   method status  scope        approximation tune parameter selected_parameter tuning_score screen seconds
# 1                  L0glm     ok single                 none  gic    lambda       2.044997e+01     102765.2   none    0.17
# 6 adaptive_relaxed_lasso     ok single identity_poisson_wls  gic    lambda       9.646484e+04     122003.4   none    0.45
# 4          relaxed_lasso     ok single identity_poisson_wls  gic    lambda       1.044924e-01     122083.4   none    0.40
# 5         adaptive_lasso     ok single identity_poisson_wls  gic    lambda       1.245891e+04     143856.3   none    0.24
# 3                  lasso     ok single identity_poisson_wls  gic    lambda       4.850106e-03     145903.5   none    0.21
# 2                L0Learn     ok single identity_poisson_wls  gic    lambda       3.296910e-03  891202271.5   none    1.20
#   path_seconds cv_seconds note     p     n m   df     deviance deviance_true        loglik  R2_weighted adj_R2_weighted
# 1         0.17          0      10000 10000 1  999 9.607885e+03  2.665175e+03     -41167.85 9.995248e-01      0.99947208
# 6         0.45          0      10000 10000 1 1900 1.042069e+04  4.229059e+03     -41574.25 9.994833e-01      0.99936212
# 4         0.40          0      10000 10000 1 1904 1.041884e+04  4.226631e+03     -41573.32 9.994836e-01      0.99936213
# 5         0.24          0      10000 10000 1 2857 1.270289e+04  6.739153e+03     -42715.35 9.993702e-01      0.99911834
# 3         0.21          0      10000 10000 1 3046 1.088513e+04  4.980112e+03     -41806.47 9.994635e-01      0.99922856
# 2         1.20          0      10000 10000 1   14 8.911293e+08  8.913185e+08 -445600992.61 3.607428e-06     -0.00139835
#   MAD_coefs RMSE_coefs  R2_coefs support_size_true support_size_estimated  TP   FP   TN  FN precision recall   F1_score
# 1  66.33935   486.4906 0.9195383              1000                    999 845  154 8846 155 0.8458458  0.845 0.84542271
# 6 186.90603   757.4528 0.8049317              1000                   1900 906  994 8006  94 0.4768421  0.906 0.62482759
# 4 187.63062   759.1949 0.8040385              1000                   1904 905  999 8001  95 0.4753151  0.905 0.62327824
# 5 162.67725   626.9684 0.8653472              1000                   2857 949 1908 7092  51 0.3321666  0.949 0.49209230
# 3 163.04241   626.1371 0.8656842              1000                   3046 950 2096 6904  50 0.3118844  0.950 0.46959960
# 2 302.28576  1403.4651 0.3423337              1000                     14  13    1 8999 987 0.9285714  0.013 0.02564103
#         MCC         gic        ebic          bic          aic        mbic
# 1 0.8282571    102765.2    100738.0     91536.83     84333.70    107169.3
# 6 0.6083765    122003.4    118147.8    100648.15     86948.50    130379.5
# 4 0.6066987    122083.4    118219.6    100683.14     86954.65    130477.1
# 5 0.4894331    143856.3    138058.6    111744.65     91144.70    156451.2
# 3 0.4674395    145903.5    139722.3    111667.64     89704.94    159331.7
# 2 0.1034134 891202271.5 891202243.1 891202114.16 891202013.22 891202333.2

# l0glmFit is the clear winner



## 1.2 SINGLE CHANNEL GAUSSIAN ####
n <- p <- 10000
npeaks <- 1000
sim.gaus <- simul(method = 1, n = n, p = p, s = npeaks/p, seed = 2)

system.time(L0glm.chol <- l0glmFit(sim.gaus$X, sim.gaus$y,
                                    family = gaussian(identity), 
                                    lambda = "gic", 
                                    nonnegative = TRUE,
                                    solver = "eigen",
                                    method = 1L)) 
# 0.14s

# The L0glm candidate uses the fastest Eigen Cholesky LLT solver and a lambda GIC preset; 
# other methods share GIC scoring. glmnet uses 10 lambda values. 
system.time(comparison.sc.gaus <- sparseGLM(sim.gaus, family = gaussian(identity),
  method = comparison_methods, nlambda = 10L,
  tune = "gic", nonnegative = TRUE, intercept = FALSE,
  support.size = seq_len(min(p, n - 1L, 2 * npeaks)),
  lambda = list(L0glm = L0glm.chol$init.lambda),
  options = list(abess = list(warm.start = TRUE)))) # 1.86s
plot(comparison.sc.gaus)

benchmarks_1_2 <- l0glmMetrics(comparison.sc.gaus)
benchmarks_1_2 <- benchmarks_1_2[order(benchmarks_1_2[,"MCC"], decreasing = TRUE),]
benchmarks_1_2[,!colSums(is.na(benchmarks_1_2))]
#                   method status  scope approximation tune parameter selected_parameter tuning_score seconds path_seconds
# 1                  L0glm     ok single          none  gic    lambda       2.044997e+01     81098.91    0.25         0.25
# 6 adaptive_relaxed_lasso     ok single          none  gic    lambda       3.034891e+11    145973.36    0.31         0.31
# 5         adaptive_lasso     ok single          none  gic    lambda       3.034891e+11    157320.70    0.23         0.23
# 4          relaxed_lasso     ok single          none  gic    lambda       3.543942e+00    148023.42    0.49         0.49
# 3                  lasso     ok single          none  gic    lambda       1.644952e-01    170989.82    0.15         0.15
# 2                L0Learn     ok single          none  gic    lambda       6.545630e-03    190173.27    0.53         0.53
#   cv_seconds note     p     n m   df     deviance deviance_true    loglik R2_weighted adj_R2_weighted MAD_coefs RMSE_coefs
# 1          0      10000 10000 1 1043 2.308267e+05  5.598291e+04 -29884.80   0.9999986       0.9999984  19.06869   216.2561
# 6          0      10000 10000 1  817 2.406940e+08  2.405623e+08 -64632.87   0.9985318       0.9984014 103.74671   552.6626
# 5          0      10000 10000 1  836 7.201108e+08  7.202037e+08 -70112.26   0.9956074       0.9952072 109.22615   537.7683
# 4          0      10000 10000 1 2298 1.429475e+07  1.412862e+07 -50514.70   0.9999128       0.9998868 242.42759   865.5736
# 3          0      10000 10000 1 3457 1.328163e+07  1.313358e+07 -50147.15   0.9999190       0.9998762 209.48527   707.9322
# 2          0      10000 10000 1   13 1.035381e+11  1.035362e+11 -94953.71   0.3684351       0.3676764 300.30677  1388.1245
#    R2_coefs support_size_true support_size_estimated  TP   FP   TN  FN precision recall   F1_score       MCC       gic
# 1 0.9839711              1000                   1042 949   93 8907  51 0.9107486  0.949 0.92948090 0.9217078  81098.91
# 6 0.8952481              1000                    816 672  144 8856 328 0.8235294  0.672 0.74008811 0.7188927 145973.36
# 5 0.9028029              1000                    835 674  161 8839 326 0.8071856  0.674 0.73460490 0.7115234 157320.70
# 4 0.7513042              1000                   2297 892 1405 7595 108 0.3883326  0.892 0.54109797 0.5248353 148023.42
# 3 0.8281221              1000                   3456 928 2528 6472  72 0.2685185  0.928 0.41651706 0.4082175 170989.82
# 2 0.3479385              1000                     12  12    0 9000 988 1.0000000  0.012 0.02371542 0.1039855 190173.27
#        ebic       bic       aic      mbic
# 1  78982.36  69375.98  61855.59  85696.94
# 6 144315.43 136790.59 130899.74 149575.08
# 5 155624.21 147924.37 141896.52 161006.17
# 4 143360.12 122194.76 105625.40 158154.08
# 3 163974.59 132134.44 107208.29 186229.89
# 2 190146.89 190027.15 189933.42 190230.58

# l0glmFit is the clear winner



## 1.3 MULTICHANNEL POISSON ####
n <- p <- 10000
npeaks <- 1000
sim.pois.mc <- simul(method = 4, n = n, p = p, s = npeaks/p, seed = 2)

system.time(L0glm.chol.pois <- l0glmFit(sim.pois.mc$X, sim.pois.mc$y,
                                         family = poisson(identity), 
                                         lambda = "gic",
                                         nonnegative = TRUE,
                                         grouped = TRUE,
                                         ncores = 4,
                                         solver = "eigen",
                                         method = 1L)) 
# 10.06s with ncores = 4

# L0glm fits the exact grouped Poisson model. glmnet uses a
# labelled Gaussian approximation with shared row weights 1/(rowMeans(y)+0.1).
# L0Learn uses independent channel-specific WLS; all are scored as Poisson.
system.time(comparison.mc.pois <- sparseGLM(sim.pois.mc, family = poisson(identity),
  method = comparison_methods, nlambda = 10L, tune = "gic", nonnegative = TRUE,
  intercept = FALSE, approximation = "wls_grouped",
  support.size = seq_len(min(p, n - 1L, 2 * npeaks)),
  lambda = list(L0glm = L0glm.chol.pois$init.lambda),
  options = list(L0glm = list(), abess = list(warm.start = TRUE))))
# 2690s=44 mins

plot(comparison.mc.pois, channels = "mean")
plot(comparison.mc.pois, which = 1, type = "heatmap") # L0glm
plot(comparison.mc.pois, which = 4, type = "heatmap") # relaxed LASSO

system.time(benchmarks_1_3 <- l0glmMetrics(comparison.mc.pois)) # 7s
benchmarks_1_3 <- benchmarks_1_3[order(benchmarks_1_3[,"MCC"], decreasing = TRUE),]
benchmarks_1_3[,!colSums(is.na(benchmarks_1_3))]
#                   method status       scope            approximation tune parameter selected_parameter tuning_score seconds
# 1                  L0glm     ok     grouped                     none  gic    lambda       2.044997e+01      5615963   24.18
# 4          relaxed_lasso     ok     grouped identity_poisson_row_wls  gic    lambda       1.242287e+00      6138687  805.38
# 6 adaptive_relaxed_lasso     ok     grouped identity_poisson_row_wls  gic    lambda       1.100621e+05      6191576  719.39
# 5         adaptive_lasso     ok     grouped identity_poisson_row_wls  gic    lambda       3.955426e+04      7348579  209.50
# 3                  lasso     ok     grouped identity_poisson_row_wls  gic    lambda       1.604476e-01      7413616  185.72
# 2                L0Learn     ok independent     identity_poisson_wls  gic    lambda       2.016690e-05     10236143  657.78
#   path_seconds cv_seconds note     p     n   m     df  deviance deviance_true   loglik R2_weighted adj_R2_weighted
# 1        24.15          0      10000 10000 226  73844  807023.9      167913.4 -2052928   0.9924739       0.9922197
# 4       805.35          0      10000 10000 226  94168  914123.0      289338.6 -2106477   0.9913598       0.9909842
# 6       719.37          0      10000 10000 226 103963  766704.6      158833.0 -2032768   0.9933573       0.9930370
# 5       209.48          0      10000 10000 226 156644  846383.2      217374.9 -2072607   0.9927158       0.9921733
# 3       185.69          0      10000 10000 226 158306  877432.0      247170.1 -2088132   0.9924869       0.9919210
# 2       657.77          0      10000 10000 226  49057 5934097.6     5471064.8 -4616465   0.9126006       0.9106613
#   MAD_coefs RMSE_coefs  R2_coefs support_size_true support_size_estimated   TP   FP   TN FN precision recall  F1_score
# 1 0.4451336   11.69795 0.9968233              1000                    955  944   11 8989 56 0.9884817  0.944 0.9657289
# 4 0.9253209   18.32347 0.9921854              1000                   1218  909  309 8691 91 0.7463054  0.909 0.8196573
# 6 0.8376681   15.84627 0.9941978              1000                   1462  967  495 8505 33 0.6614227  0.967 0.7855402
# 5 1.2232663   19.87182 0.9911969              1000                   1813  986  827 8173 14 0.5438500  0.986 0.7010309
# 3 1.3251157   21.70045 0.9894781              1000                   1840  980  860 8140 20 0.5326087  0.980 0.6901408
# 2 3.4764744   77.19635 0.8702346              1000                   3687 1000 2687 6313  0 0.2712232  1.000 0.4267122
#         MCC within_row_TP within_row_FP within_row_TN within_row_FN within_row_precision within_row_recall
# 1 0.9623313         64452          8559        135959          4374            0.8827711         0.9364484
# 4 0.8023123         60980         16249        123498          4707            0.7895998         0.9283420
# 6 0.7743982         64291         14723        133119          6409            0.8136659         0.9093494
# 5 0.6962279         69676         22885        127728          2547            0.7527576         0.9647342
# 3 0.6847591         69464         23673        126032          2311            0.7458260         0.9678022
# 2 0.4361743         33296           646        152088         39970            0.9809675         0.4544536
#  within_row_F1_score within_row_MCC      gic     ebic     bic     aic     mbic
# 1           0.9088179      0.8643182  5615963  5466112 4785984 4253544  5941502
# 4           0.8533684      0.7819166  6138687  5947593 5080274 4401291  6553824
# 6           0.8588509      0.7884679  6191576  5980606 5023071 4273462  6649894
# 5           0.8456646      0.7719728  7348579  7030704 5587959 4458503  8039139
# 3           0.8424372      0.7676552  7413616  7092368 5634316 4492876  8111503
# 2           0.6211477      0.5898789 10236143 10136593 9684761 9331043 10452409

# l0glmFit is the clear winner


## 1.4 MULTICHANNEL GAUSSIAN ####
n <- p <- 10000
npeaks <- 1000
sim.gaus.mc <- simul(method = 3, n = n, p = p, s = npeaks/p, seed = 1)

system.time(L0glm.chol.gaus <- l0glmFit(sim.gaus.mc$X, sim.gaus.mc$y,
                                         family = gaussian(identity), 
                                         lambda = "gic", 
                                         nonnegative = TRUE,
                                         grouped = TRUE,
                                         ncores = 4,
                                         solver = "eigen",
                                         method = 1L)) 
# 37.64s with ncores = 4

# Auto grouping uses native grouped fits where supported, independent paths
# otherwise; the table records the scope of every method.
system.time(comparison.mc.gaus <- sparseGLM(sim.gaus.mc, family = gaussian(identity),
  method = comparison_methods, nlambda = 10L, tune = "gic", nonnegative = TRUE, intercept = FALSE,
  support.size = seq_len(min(p, n - 1L, 2 * npeaks)),
  lambda = list(L0glm = L0glm.chol.gaus$init.lambda),
  options = list(L0glm = list(), abess = list(warm.start = TRUE))))
# 3053s=50 mins

plot(comparison.mc.gaus, channels = "mean")
plot(comparison.mc.gaus, which = 1, type = "heatmap") # L0glm
plot(comparison.mc.gaus, which = 4, type = "heatmap") # relaxed LASSO

system.time(benchmarks_1_4 <- l0glmMetrics(comparison.mc.gaus)) # 3.49s
benchmarks_1_4 <- benchmarks_1_4[order(benchmarks_1_4[,"MCC"], decreasing = TRUE),]
benchmarks_1_4[,!colSums(is.na(benchmarks_1_4))]
#                   method status       scope approximation tune parameter selected_parameter tuning_score seconds
# 1                  L0glm     ok     grouped          none  gic    lambda       2.044997e+01     16225216  195.75
# 5         adaptive_lasso     ok     grouped          none  gic    lambda       5.422378e+02     20728348   53.93
# 6 adaptive_relaxed_lasso     ok     grouped          none  gic    lambda       5.422378e+02     18419427   91.89
# 4          relaxed_lasso     ok     grouped          none  gic    lambda       3.970177e+00     16939870 1522.27
# 3                  lasso     ok     grouped          none  gic    lambda       1.426808e+00     20034062  197.91
# 2                L0Learn     ok independent          none  gic    lambda       3.888993e-06     16633518 1457.22
#   path_seconds cv_seconds note     p     n   m     df  deviance deviance_true   loglik R2_weighted adj_R2_weighted
# 1       195.72          0      10000 10000 232 187597  29104315       6321224 -6194432   0.9999165       0.9999092
# 5        53.92          0      10000 10000 232 179371 793931184     767698846 -8530109   0.9977234       0.9975329
# 6        91.88          0      10000 10000 232 163182 275652065     251281660 -7541180   0.9992096       0.9991499
# 4      1522.25          0      10000 10000 232 218471  30859110       9243342 -6236073   0.9999115       0.9999023
# 3       197.86          0      10000 10000 232 315674  68711708      45861548 -6789270   0.9998030       0.9997720
# 2      1457.21          0      10000 10000 232  88852 734111510     708697966 -7408249   0.9978950       0.9978114
#   MAD_coefs RMSE_coefs  R2_coefs support_size_true support_size_estimated   TP   FP   TN  FN precision recall  F1_score
# 1 0.2830137   3.234918 0.9996738              1000                    991  984    7 8993  16 0.9929364  0.984 0.9884480
# 5 1.7465920  23.622649 0.9826592              1000                    870  758  112 8888 242 0.8712644  0.758 0.8106952
# 6 1.3009148  23.126024 0.9833310              1000                    870  758  112 8888 242 0.8712644  0.758 0.8106952
# 4 1.2541608  23.373356 0.9830773              1000                   1561  958  603 8397  42 0.6137092  0.958 0.7481453
# 3 1.7883043  27.525177 0.9769680              1000                   1815  989  826 8174  11 0.5449036  0.989 0.7026643
# 2 3.0368554  63.963533 0.8796476              1000                   9987 1000 8987   13   0 0.1001302  1.000 0.1820333
#          MCC within_row_TP within_row_FP within_row_TN within_row_FN within_row_precision within_row_recall
# 1 0.98718335         67729        118827         37841          3891            0.3630492         0.9456716
# 5 0.79360839         51642        108142         15045          1027            0.3231988         0.9805009
# 6 0.79360839         50227         98067         25120          2442            0.3386988         0.9536350
# 4 0.73646483         60286         86479         66109          9382            0.4107655         0.8653327
# 3 0.69835078         68844        117341         40034          3229            0.3697613         0.9551982
# 2 0.01202632         24418          5811        153119         48652            0.8077674         0.3341727
#   within_row_F1_score within_row_MCC      gic     ebic      bic      aic     mbic
# 1           0.5246731      0.2247461 16225216 15844529 14116696 12764058 17052232
# 5           0.4861499      0.1631332 20728348 20364353 18712285 17418959 21519099
# 6           0.4998632      0.1985036 18419427 18088284 16585323 15408725 19138809
# 4           0.5570870      0.2924725 16939870 16496530 14484338 12909088 17902992
# 3           0.5331413      0.2487036 20034062 19393470 16486005 14209888 21425700
# 2           0.4727635      0.4106558 16633518 16453212 15634855 14994202 17025219

# l0glmFit is the clear winner





# 2. COMPARISON OF L0glm SOLVERS FOR SINGLE & MULTICHANNEL ####

## 2.1 SINGLE CHANNEL, SPARSE  ####
n <- p <- 10000 
npeaks <- 1000
sim.gaus <- simul(method = 1, n = n, p = p, s = npeaks/p, sparse = T) 
sim.pois <- simul(method = 2, n = n, p = p, s = npeaks/p, sparse = T) 

### L0GLM FIT SINGLE CHANNEL WITH POISSON NOISE ####

## tuning plot (full range)
# Compare the GIC preset with a grid using a common pilot start.
system.time(tune_full <- l0glmTune(X=sim.pois$X,
                                    y=sim.pois$y, 
                                    tune = "gic",
                                    nlambda = 10,
                                    start = "pilot", # common pooled start stabilizes the path
                                    prop_range = 1, # full range up to lambda_max
                                    plot = TRUE, 
                                    coefs.true = sim.pois$beta_true, 
                                    family = poisson(identity), 
                                    nonnegative = TRUE, 
                                    solver = "eigen", 
                                    method = 1L))
# A common pilot retains the union of pilot supports and initializes every lambda independently.
tune_full$lambda_best
tune_full$lambda_preset # 20.44997 (lambda preset)
min(tune_full$metrics$gic)

## tuning plot (zoomed in)
# with lambda grid more zoomed in around lambda_preset preset
system.time(tune_zoomed <- l0glmTune(X=sim.pois$X,
                                      y=sim.pois$y, 
                                      tune = "gic",
                                      nlambda = 10,
                                      start = "pilot", # common pooled start stabilizes the path
                                      prop_range = 0.2, # to zoom in around lambda_preset preset
                                      plot = TRUE, 
                                      coefs.true = sim.pois$beta_true, 
                                      family = poisson(identity), 
                                      nonnegative = TRUE, 
                                      solver = "eigen", 
                                      method = 1L))
# Endpoint probes stop once nonzero support is certified; displayed grid fits receive full refinement.
# The dashed line marks the best grid value; the dotted line is the preset.
# They can differ as the selected support and coefficient estimates change along the path.
tune_zoomed$lambda_best
tune_zoomed$lambda_preset # 20.44997 (lambda preset)
min(tune_zoomed$metrics$gic)

# The preset is a useful calibration; the path can refine lambda for correlated predictors.
L0glm:::presetLambda(IC="gic", p, n) # 20.44997 (lambda preset)
tune_zoomed$lambda_preset # 20.44997 (lambda preset)
tune_zoomed$lambda_best

# The pilot gives a common starting point while each lambda selects its own support.
# Discrete coefficient changes reflect changes in the selected sparse model.
system.time(tune_full_bigglm <- l0glmTune(X=sim.pois$X,
                                           y=sim.pois$y, 
                                           tune = "gic",
                                           nlambda = 10,
                                           prop_range = 1, # full range up to lambda_max
                                           plot = TRUE, 
                                           coefs.true = sim.pois$beta_true, 
                                           family = poisson(identity), 
                                           nonnegative = TRUE, 
                                           solver = "bigGlm"))
# 19.50s

system.time(tune_full_osqp <- l0glmTune(X=sim.pois$X,
                                         y=sim.pois$y, 
                                         tune = "gic",
                                         nlambda = 10,
                                         prop_range = 1, # full range up to lambda_max
                                         plot = TRUE, 
                                         coefs.true = sim.pois$beta_true, 

                                         family = poisson(identity), 
                                         nonnegative = TRUE, 
                                         solver = "osqp"))
# 30.94s

# 5-fold cross validation is also supported:
system.time(tune_zoomed_cv <- l0glmTune(X=sim.pois$X,
                                         y=sim.pois$y, 
                                         tune = "cv",
                                         nfolds = 5,
                                         prop_range = 1, # to zoom in around lambda_preset preset
                                         lambda_min = 2,
                                         plot = TRUE, 
                                         coefs.true = sim.pois$beta_true, 
                                         family = poisson(identity), 
                                         nonnegative = TRUE, 
                                         solver = "eigen", 
                                         method = 1L)) # 10.94s
tune_zoomed_cv$lambda_best # 2


#### 1.1.1 L0GLM FIT WITH CHOLESKY LLT SOLVER ####
system.time(L0glm.chol.pois <- l0glmFit(sim.pois$X, sim.pois$y,
                                         family = poisson(identity), 
                                         lambda = "gic", # or tune_zoomed$lambda_best
                                         nonnegative = TRUE,
                                         # verbose = TRUE,
                                         solver = "eigen",
                                         method = 1L)) 
# 0.18s

#### 1.1.2 L0GLM FIT WITH CG SOLVER WITH Incomplete Cholesky LLT preconditioner ####
system.time(L0glm.CG.pois <- l0glmFit(sim.pois$X, sim.pois$y,
                                       family = poisson(identity), 
                                       lambda = "gic", # or tune_zoomed$lambda_best
                                       nonnegative = TRUE,
                                       solver = "eigen",
                                       method = 3L)) 
# 0.31s

#### 1.1.3 L0GLM FIT WITH preconditioned LSCG SOLVER ####
system.time(L0glm.LSCG.pois <- l0glmFit(sim.pois$X, sim.pois$y,
                                         family = poisson(identity), 
                                         lambda = "gic", # or tune_zoomed$lambda_best
                                         nonnegative = TRUE,
                                         solver = "eigen",
                                         method = 4L)) 
# 0.75s

#### 1.1.4 L0GLM FIT WITH OSQP SOLVER ####
system.time(L0glm.osqp.pois <- l0glmFit(sim.pois$X, sim.pois$y,
                                         family = poisson(identity), 
                                         lambda = "gic", # or tune_zoomed$lambda_best
                                         nonnegative = TRUE,
                                         solver = "osqp")) 
# 45.70s

#### 1.1.5 L0GLM FIT WITH bigGlm::bigGlm SOLVER ####
system.time(L0glm.bigGlm.pois <- l0glmFit(sim.pois$X, sim.pois$y,
                                           family = poisson(identity), 
                                           lambda = "gic", # or tune_zoomed$lambda_best
                                           nonnegative = TRUE,
                                           solver = "bigGlm")) 
# 1.87s with n=p=10000




### BENCHMARKS TABLE ####

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
benchmarks_1_1[, !colSums(is.na(benchmarks_1_1)), drop = FALSE]
#          method     p     n m   df deviance deviance_true    loglik R2_weighted adj_R2_weighted MAD_coefs RMSE_coefs
# 4   L0glm_PLSCG 10000 10000 1 1011 9349.100      2988.750 -41036.03   0.9995603       0.9995109  80.28585   553.3250
# 5     L0glm_PCG 10000 10000 1 1011 9349.100      2988.750 -41036.03   0.9995603       0.9995109  80.28585   553.3250
# 3 L0glm_cholLLT 10000 10000 1 1011 9349.100      2988.750 -41036.03   0.9995603       0.9995109  80.28585   553.3250
# 2  L0glm_bigglm 10000 10000 1 1012 9462.048      3074.622 -41092.50   0.9995555       0.9995054  80.66690   555.5264
# 1    L0glm_osqp 10000 10000 1 1012 9462.048      3074.622 -41092.50   0.9995555       0.9995054  80.66690   555.5264
#    R2_coefs support_size_true support_size_estimated  TP  FP   TN  FN precision recall  F1_score      MCC      gic     ebic
# 4 0.9093003              1000                   1011 843 168 8832 157 0.8338279  0.843 0.8383889 0.820338 102747.0 100695.4
# 5 0.9093003              1000                   1011 843 168 8832 157 0.8338279  0.843 0.8383889 0.820338 102747.0 100695.4
# 3 0.9093003              1000                   1011 843 168 8832 157 0.8338279  0.843 0.8383889 0.820338 102747.0 100695.4
# 2 0.9085972              1000                   1012 840 172 8828 160 0.8300395  0.840 0.8349901 0.816552 102880.4 100826.7
# 1 0.9085972              1000                   1012 840 172 8828 160 0.8300395  0.840 0.8349901 0.816552 102880.4 100826.7
#        bic      aic     mbic
# 4 91383.71 84094.06 107203.9
# 5 91383.71 84094.06 107203.9
# 3 91383.71 84094.06 107203.9
# 2 91505.87 84209.01 107341.7
# 1 91505.87 84209.01 107341.7

# cholesky LLT, CG & LSCG tied in terms of solution quality (deviance_true) & LLT fastest
# eigen solvers all equally good & cholesky the fastest
# osqp solver the worst & the slowest


## 1.2. L0GLM FIT SINGLE CHANNEL WITH GAUSSIAN NOISE ####

## tuning plot (full range)
# here preset lambda coincides with best lambda
system.time(tune_full <- l0glmTune(X=sim.gaus$X,
                                    y=sim.gaus$y, 
                                    tune = "gic",
                                    nlambda = 10,
                                    prop_range = 1, # full range up to lambda_max
                                    plot = TRUE, 
                                    coefs.true = sim.gaus$beta_true, 
                                    family = gaussian(identity), 
                                    nonnegative = TRUE, 
                                    solver = "eigen", 
                                    method = 1L)) # 7.28
tune_full$lambda_best # 20.44997 (best lambda tune)
tune_full$lambda_preset # 20.44997 (lambda preset)
min(tune_full$metrics$gic) # 82342.78



## tuning plot (zoomed in)
# tuning plots - with lambda grid more zoomed in around lambda_preset preset
system.time(tune_zoomed <- l0glmTune(X=sim.gaus$X,
                                      y=sim.gaus$y, 
                                      tune = "gic",
                                      nlambda = 10,
                                      prop_range = 0.2, # to zoom in around lambda_preset preset
                                      plot = TRUE, 
                                      coefs.true = sim.gaus$beta_true, 
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
L0glm:::presetLambda(IC="gic", p, n) # 20.44997 (lambda preset)
tune_zoomed$lambda_best # 34.03505 (best lambda tune)
tune_zoomed$lambda_preset # 20.44997 (lambda preset)


#### 1.2.1 L0GLM FIT WITH Cholesky LLT SOLVER ####
system.time(L0glm.chol.gaus <- l0glmFit(sim.gaus$X, sim.gaus$y,
                                         family = gaussian(identity), 
                                         lambda = "gic", # or tune_zoomed$lambda_best
                                         nonnegative = TRUE,
                                         solver = "eigen",
                                         method = 1L)) 
# 0.16s

# Compare default clamping with an exact constrained Eigen ridge polish.
# On this Gaussian example, the polish matches the bounded bigGlm/OSQP solution.
system.time(L0glm.chol.gaus.polished <- l0glmFit(sim.gaus$X, sim.gaus$y,
  family = gaussian(identity), lambda = "gic", nonnegative = TRUE,
  solver = "eigen", method = 1L,
  control.solver = controlSolverGen("eigen", 1L, nonnegative.polish = TRUE)))

#### 1.2.2 L0GLM FIT WITH PCG SOLVER WITH Incomplete Cholesky preconditioner ####
system.time(L0glm.CG.gaus <- l0glmFit(sim.gaus$X, sim.gaus$y,
                                       family = gaussian(identity), 
                                       lambda = "gic", # or tune_zoomed$lambda_best
                                       nonnegative = TRUE,
                                       solver = "eigen",
                                       method = 3L)) 
# 2.05s


#### 1.2.3 L0GLM FIT WITH PLSCG SOLVER ####
system.time(L0glm.LSCG.gaus <- l0glmFit(sim.gaus$X, sim.gaus$y,
                                         family = gaussian(identity), 
                                         lambda = "gic", # or tune_zoomed$lambda_best
                                         nonnegative = TRUE,
                                         solver = "eigen",
                                         method = 4L)) 
# 3.30s


#### 1.2.4 L0GLM FIT WITH OSQP SOLVER ####
system.time(L0glm.osqp.gaus <- l0glmFit(sim.gaus$X, sim.gaus$y,
                                         family = gaussian(identity), 
                                         lambda = "gic", # or tune_zoomed$lambda_best
                                         nonnegative = TRUE,
                                         solver = "osqp")) 
# 115s


#### 1.2.5 L0GLM FIT WITH bigGlm::bigGlm SOLVER ####
system.time(L0glm.bigGlm.gaus <- l0glmFit(sim.gaus$X, sim.gaus$y,
                                           family = gaussian(identity), 
                                           lambda = "gic", # or tune_zoomed$lambda_best
                                           nonnegative = TRUE,
                                           solver = "bigGlm")) 
# 11s



### 1.2.9 BENCHMARKS TABLE ####

X = sim.gaus$X
y = sim.gaus$y
# Score the fitted models with their fitted dispersions and common S3 metrics.
gaussian_solver_fits <- list(
  L0glm_osqp = L0glm.osqp.gaus,
  L0glm_bigglm = L0glm.bigGlm.gaus,
  L0glm_cholLLT = L0glm.chol.gaus,
  L0glm_cholLLT_constrained = L0glm.chol.gaus.polished,
  L0glm_PLSCG = L0glm.LSCG.gaus,
  L0glm_PCG = L0glm.CG.gaus)
benchmarks_1_2 <- do.call(rbind, lapply(names(gaussian_solver_fits), function(method) {
  cbind(method = method, l0glmMetrics(gaussian_solver_fits[[method]],
                                    beta_true = sim.gaus$beta_true))
}))
benchmarks_1_2 <- benchmarks_1_2[order(benchmarks_1_2[, "deviance_true"]), ]
benchmarks_1_2[, !colSums(is.na(benchmarks_1_2)), drop = FALSE]
# Constrained Eigen LLT: coefficient R2 0.9976, F1 0.9668, 976 TP and 43 FP.
# Default Eigen clamping: coefficient R2 0.9806, F1 0.9098, 938 TP and 124 FP.
# The constrained LLT polish matches bigGlm here and took about 0.93 s.
# Clamping remains the package default; compare both strategies for the data at hand.


# 2. MULTICHANNEL ####

# TODO: test option precluster = T or F?
# (see ?glmCluster)

# Simulate MC n=p=500
n <- p <- 500
npeaks <- 50
sim.gaus.mc <- simul(method = 3, n = n, p = p, s = npeaks/p, seed = 1)
sim.pois.mc <- simul(method = 4, n = n, p = p, s = npeaks/p, seed = 2)
dim(sim.pois.mc$y) # 500 timepoints x 168 m/z channels
m <- ncol(sim.pois.mc$y)


## 2.1. MULTICHANNEL WITH POISSON NOISE ####

### 2.1.1 GROUPED L0GLM FIT WITH grouped = FALSE (not great) ####

## tuning plot (full range)
# here preset lambda would be ca 10x higher than best lambda
system.time(tune_full <- l0glmTune(X=sim.pois.mc$X, y=sim.pois.mc$y,
                                    tune = "gic",
                                    nlambda = 10,
                                    prop_range = 1, # to zoom in around lambda_preset preset
                                    plot = TRUE, coefs.true = sim.pois.mc$beta_true, 
                                    family = poisson(identity),
                                    nonnegative = TRUE, 
                                    grouped = FALSE,
                                    solver = "eigen", method = 1L,
                                    ncores = 4)) # historical 2.96s, TODO check why now 25s
tune_full$lambda_best # 4.005619 (best lambda tune)
tune_full$lambda_preset # 11.35348 (lambda preset)
L0glm:::presetLambda(IC="gic", p, n) # 11.35348 (lambda preset)

## tuning plot (zoomed in)
# here preset lambda would be ca 5x higher than best lambda
system.time(tune_zoomed <- l0glmTune(X=sim.pois.mc$X, y=sim.pois.mc$y,
                                      tune = "gic",
                                      nlambda = 10,
                                      prop_range = 0.2, # to zoom in around lambda_preset preset
                                      plot = TRUE, coefs.true = sim.pois.mc$beta_true, 
                                      family = poisson(identity),
                                      nonnegative = TRUE, 
                                      grouped = FALSE,
                                      solver = "eigen", method = 1L,
                                      ncores = 4)) # historical 3.31s
tune_zoomed$lambda_best # 8.022493 (best lambda tune)
tune_zoomed$lambda_preset # 11.35348 (lambda preset)
L0glm:::presetLambda(IC="gic", p, n) # 11.35348 (lambda preset)


#### Cholesky LLT ####
system.time(L0glm.chol <- l0glmFit(sim.pois.mc$X, sim.pois.mc$y,
                                    family = poisson(identity), 
                                    lambda = "gic", # or tune_zoomed$lambda_best
                                    nonnegative = TRUE,
                                    grouped = FALSE,
                                    solver = "eigen",
                                    method = 1L)) 
# 0.29s

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
system.time(L0glm.osqp <- l0glmFit(sim.pois.mc$X, sim.pois.mc$y,
                                    family = poisson(identity), 
                                    lambda = "gic", # or tune_zoomed$lambda_best
                                    nonnegative = TRUE,ncores=1L,
                                    grouped = FALSE,
                                    solver = "osqp")) 
# 13.64s





### 2.1.2 GROUPED L0GLM FIT WITH grouped = TRUE (better) ####

## tuning plot (full range)
# here preset lambda would be ca 10x higher than best lambda
system.time(tune_full <- l0glmTune(X=sim.pois.mc$X, y=sim.pois.mc$y,
                                    tune = "gic",
                                    nlambda = 10,
                                    prop_range = 1, # to zoom in around lambda_preset preset
                                    plot = TRUE, coefs.true = sim.pois.mc$beta_true, 
                                     family = poisson(identity),
                                    nonnegative = TRUE, 
                                    grouped = TRUE,
                                    solver = "eigen", method = 1L,
                                    ncores = 4)) # historical 4.47s, now 23s
tune_full$lambda_best # 2.830413 (best lambda tune)
tune_full$lambda_preset # 11.35348 (lambda preset)
L0glm:::presetLambda(IC="gic", p, n) # 11.35348 (lambda preset)

## tuning plot (zoomed in)
# here preset lambda would be ca 5x higher than best lambda
system.time(tune_zoomed <- l0glmTune(X=sim.pois.mc$X, y=sim.pois.mc$y,
                                      tune = "gic",
                                      nlambda = 10,
                                      prop_range = 0.2, # to zoom in around lambda_preset preset
                                      plot = TRUE, coefs.true = sim.pois.mc$beta_true, 
                                       family = poisson(identity),
                                      nonnegative = TRUE, 
                                      grouped = TRUE,
                                      solver = "eigen", method = 1L,
                                      ncores = 4)) # historical 5.75s, now 19s
tune_zoomed$lambda_best # 8.022493 (best lambda tune)
tune_zoomed$lambda_preset # 11.35348 (lambda preset)
L0glm:::presetLambda(IC="gic", p, n) # 11.35348 (lambda preset)



#### Cholesky LLT ####
system.time(L0glm.chol <- l0glmFit(sim.pois.mc$X, sim.pois.mc$y,
                                    family = poisson(identity), 
                                    lambda = "gic", # or tune_zoomed$lambda_best
                                    nonnegative = TRUE,
                                    grouped = TRUE,
                                    solver = "eigen",
                                    method = 1L,
                                    ncores = 4)) 
# 0.42s



#### OSQP ####
system.time(L0glm.osqp <- l0glmFit(sim.pois.mc$X, sim.pois.mc$y,
                                    family = poisson(identity), 
                                    lambda = "gic", # or tune_zoomed$lambda_best
                                    nonnegative = TRUE,ncores=1L,
                                    grouped = TRUE,
                                    solver = "osqp"))
# 41.95s




## 2.1.3 GROUPED L0GLM FIT WITH grouped = TRUE (best) ####

## tuning plot (full range)
# here preset lambda would be ca 10x higher than best lambda
system.time(tune_full <- l0glmTune(X=sim.pois.mc$X, y=sim.pois.mc$y,
                                    tune = "gic",
                                    nlambda = 10,
                                    prop_range = 1, # to zoom in around lambda_preset preset
                                    plot = TRUE, coefs.true = sim.pois.mc$beta_true, 
                                     family = poisson(identity),
                                    nonnegative = TRUE, 
                                    grouped = TRUE,
                                    solver = "eigen", method = 1L,
                                    ncores = 4)) # openmp 5s
tune_full$lambda_best # 2 (best lambda tune)
tune_full$lambda_preset # 11.35348 (lambda preset)
L0glm:::presetLambda(IC="gic", p, n) # 11.35348 (lambda preset)

## tuning plot (zoomed in)
# here preset lambda would be ca 10x higher than best lambda
system.time(tune_zoomed <- l0glmTune(X=sim.pois.mc$X, y=sim.pois.mc$y,
                                      tune = "gic",
                                      nlambda = 10,
                                      prop_range = 0.2, # to zoom in around lambda_preset preset
                                      plot = TRUE, coefs.true = sim.pois.mc$beta_true, 
                                       family = poisson(identity),
                                      nonnegative = TRUE, 
                                      grouped = TRUE,
                                      solver = "eigen", method = 1L,
                                      ncores = 4)) # openmp 5.66s
tune_zoomed$lambda_best # 8.022493 (best lambda tune)
tune_zoomed$lambda_preset # 11.35348 (lambda preset)
L0glm:::presetLambda(IC="gic", p, n) # 11.35348 (lambda preset)


#### Cholesky LLT ####
system.time(L0glm.chol <- l0glmFit(sim.pois.mc$X, sim.pois.mc$y,
                                    family = poisson(identity), 
                                    lambda = "gic", # or tune_zoomed$lambda_best
                                    nonnegative = TRUE,
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



## Fit grouped = TRUE
#### OSQP ####
system.time(L0glm.osqp <- l0glmFit(sim.pois.mc$X, sim.pois.mc$y,
                                    family = poisson(identity), 
                                    lambda = "gic", # or tune_zoomed$lambda_best
                                    nonnegative = TRUE,ncores=1L,
                                    grouped = TRUE,
                                    solver = "osqp")) 
# 58s




## 2.2. MULTICHANNEL WITH GAUSSIAN NOISE ####

### 2.2.1 GROUPED L0GLM FIT WITH grouped = FALSE (not great) ####

## tuning plot (full range)
# here best lambda is identical to preset lambda
system.time(tune_full <- l0glmTune(X=sim.gaus.mc$X, y=sim.gaus.mc$y,
                                    tune = "gic",
                                    nlambda = 10,
                                    prop_range = 1, # to zoom in around lambda_preset preset
                                    plot = TRUE, 
                                    coefs.true = sim.gaus.mc$beta_true, 
                                    family = gaussian(identity), 
                                    nonnegative = TRUE,ncores=1L,
                                    grouped = FALSE,
                                    solver = "eigen", method = 1L)) # 85s
tune_full$lambda_best # 11.35348
tune_full$lambda_preset # 11.35348 
L0glm:::presetLambda(IC="gic", p, n) # 11.35348

## tuning plot (zoomed in)
# here best lambda is identical to preset lambda
system.time(tune_zoomed <- l0glmTune(X=sim.gaus.mc$X, y=sim.gaus.mc$y,
                                      tune = "gic",
                                      nlambda = 10,
                                      prop_range = 0.2, # to zoom in around lambda_preset preset
                                      plot = TRUE, 
                                      coefs.true = sim.gaus.mc$beta_true, 

                                      family = gaussian(identity), 
                                      nonnegative = TRUE,ncores=1L,
                                      grouped = FALSE,
                                      solver = "eigen", method = 1L)) # 84s
tune_zoomed$lambda_best # 11.35348
tune_zoomed$lambda_preset # 11.35348 
L0glm:::presetLambda(IC="gic", p, n) # 11.35348

#### Cholesky LLT ####
system.time(L0glm.chol <- l0glmFit(sim.gaus.mc$X, sim.gaus.mc$y,
                                    family = gaussian(identity), 
                                    lambda = "gic", # or tune_zoomed$lambda_best
                                    nonnegative = TRUE,ncores=1L,
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
system.time(L0glm.osqp <- l0glmFit(sim.gaus.mc$X, sim.gaus.mc$y,
                                    family = gaussian(identity), 
                                    lambda = "gic", # or tune_zoomed$lambda_best
                                    nonnegative = TRUE,ncores=1L,
                                    grouped = FALSE,
                                    solver = "osqp")) 
# 18.54s




### 2.2.2 GROUPED L0GLM FIT WITH grouped = TRUE (better) ####

## tuning plot (full range)
# here best lambda is identical to preset lambda
system.time(tune_full <- l0glmTune(X=sim.gaus.mc$X, y=sim.gaus.mc$y,
                                    tune = "gic",
                                    nlambda = 10,
                                    prop_range = 1, # to zoom in around lambda_preset preset
                                    plot = TRUE, 
                                    coefs.true = sim.gaus.mc$beta_true, 

                                    family = gaussian(identity), 
                                    nonnegative = TRUE,ncores=1L,
                                    grouped = TRUE,
                                    solver = "eigen", method = 1L)) # 242.34s
tune_full$lambda_best # 11.35348
tune_full$lambda_preset # 11.35348 
L0glm:::presetLambda(IC="gic", p, n) # 11.35348

## tuning plot (zoomed in)
# here best lambda is identical to preset lambda
system.time(tune_zoomed <- l0glmTune(X=sim.gaus.mc$X, y=sim.gaus.mc$y,
                                      tune = "gic",
                                      nlambda = 10,
                                      prop_range = 0.2, # to zoom in around lambda_preset preset
                                      plot = TRUE, 
                                      coefs.true = sim.gaus.mc$beta_true, 

                                      family = gaussian(identity), 
                                      nonnegative = TRUE,ncores=1L,
                                      grouped = TRUE,
                                      solver = "eigen", method = 1L)) # 300.46s
tune_zoomed$lambda_best # 11.35348
tune_zoomed$lambda_preset # 11.35348 
L0glm:::presetLambda(IC="gic", p, n) # 11.35348

#### Cholesky LLT ####
system.time(L0glm.chol <- l0glmFit(sim.gaus.mc$X, sim.gaus.mc$y,
                                    family = gaussian(identity), 
                                    lambda = "gic", # or tune.chol.mc$lambda_best,
                                    nonnegative = TRUE,ncores=1L,
                                    grouped = TRUE,
                                    solver = "eigen",
                                    method = 1L)) # 2.3s

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
system.time(L0glm.osqp <- l0glmFit(sim.gaus.mc$X, sim.gaus.mc$y,
                                    family = gaussian(identity), 
                                    lambda = "gic", # or tune.osqp.mc$lambda_best,
                                    nonnegative = TRUE,ncores=1L,
                                    grouped = TRUE,
                                    solver = "osqp")) 
# 72.58s



### 2.2.3 GROUPED L0GLM FIT WITH grouped = TRUE (best) ####

## tuning plot (full range)
# here best lambda is identical to preset lambda
system.time(tune_full <- l0glmTune(X=sim.gaus.mc$X, y=sim.gaus.mc$y,
                                    tune = "gic",
                                    nlambda = 10,
                                    prop_range = 1, # to zoom in around lambda_preset preset
                                    plot = TRUE, 
                                    coefs.true = sim.gaus.mc$beta_true, 
                                    family = gaussian(identity), 
                                    nonnegative = TRUE,ncores=1L,
                                    grouped = TRUE,
                                    solver = "eigen", method = 1L)) # 32.86s
tune_full$lambda_best # 11.35348 (best lambda tune)
tune_full$lambda_preset # 11.35348 (lambda preset)
L0glm:::presetLambda(IC="gic", p, n) # 11.35348 (lambda preset)

## tuning plot (zoomed in)
# here best lambda is again very close to preset lambda
system.time(tune_zoomed <- l0glmTune(X=sim.gaus.mc$X, y=sim.gaus.mc$y,
                                      tune = "gic",
                                      nlambda = 10,
                                      prop_range = 0.2, # to zoom in around lambda_preset preset
                                      plot = TRUE, 
                                      coefs.true = sim.gaus.mc$beta_true, 
                                      family = gaussian(identity), 
                                      nonnegative = TRUE,ncores=1L,
                                      grouped = TRUE,
                                      solver = "eigen", method = 1L)) # 35.14s
tune_zoomed$lambda_best # 10.59169 (best lambda tune)
tune_zoomed$lambda_preset # 11.35348 (lambda preset)
L0glm:::presetLambda(IC="gic", p, n) # 11.35348 (lambda preset)


#### Cholesky LLT ####
system.time(L0glm.chol <- l0glmFit(sim.gaus.mc$X, sim.gaus.mc$y,
                                    family = gaussian(identity), 
                                    lambda = "gic", # or tune_zoomed$lambda_best
                                    nonnegative = TRUE,ncores=1L,
                                    grouped = TRUE,
                                    solver = "eigen",
                                    method = 1L)) 
# 2.41s

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
system.time(L0glm.osqp <- l0glmFit(sim.gaus.mc$X, sim.gaus.mc$y,
                                    family = gaussian(identity), 
                                    lambda = "gic", # or tune_zoomed$lambda_best
                                    nonnegative = TRUE,ncores=1L,
                                    grouped = TRUE,
                                    solver = "osqp")) 
# 71s
