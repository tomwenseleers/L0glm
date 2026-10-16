# PEAK DECONVOLUTION AFTER PRECLUSTERING AND AGGREGATING CORRELATED MASS TRACES ####

library(L0glm)
library(fastcluster)
library(dynamicTreeCut)
library(PoiClaClu)
library(Matrix)
# PS on Windows best to install OPENBLAS as in https://github.com/david-cortes/R-openblas-in-windows?tab=readme-ov-file
# to have multithreaded BLAS for matrix operations

## 1. TEST WITH SIMULATED MULTICHANNEL DATA ####
### SIMULATE DATA ####
set.seed(123)
n = 10000
L = 41 # length of PSF
psf_fun = function(x, u, w = 2.5) exp(((x - u)^2)/(-2 * (w^2)))
psf_true = psf_fun(1:L, u=(1+L)/2)
plot(psf_true, type="l")
dat=simulate_spike_train(n=n, p=n, 
                         k=1000, # peaks
                         psf = function(x, u, w = 2.5) exp(((x - u)^2)/(-2 * (w^2))),
                         family="poisson", multichannel=T)
Y=dat$y
X_true=dat$X
dim(Y) # 10000 timepoints x 228 mass channels
true_features = which(rowSums(dat$beta_true!=0)!=0)
beta_true = dat$beta_true
length(true_features)
matplot(Y^0.5, type="l", lty=1, col=rgb(0,0,0,0.1))
plot(rowSums(Y)^0.5, type="l", lty=1)

# NOTE: without clustering and aggregating mass channels first & using X_true
# a single grouped L0glm fit using a sparse Cholesky solver takes 106s :
system.time(L0glm.chol_originaldata <- L0glm.fit(X_true, 
                                         Y, 
                                         family = poisson(identity), 
                                         lambda = "gic",
                                         nonnegative = TRUE,
                                         grouped = TRUE,
                                         relaxed = FALSE,
                                         parallelize = FALSE,
                                         solver = "eigen",
                                         method = 1L)) # 105s
selected_features = which(rowSums(L0glm.chol_originaldata$coefficients!=0)!=0)
length(selected_features) # 951
length(intersect(selected_features, true_features)) # 929 true positives
length(selected_features)-length(intersect(selected_features, true_features)) # 22 false positives
beta_est = L0glm.chol_originaldata$coefficients
rsq = summary(lm(as.vector(beta_true) ~ 0 + as.vector(beta_est)))$r.squared
rsq # 0.9912684

plot(x=as.vector(beta_est)[!((as.vector(beta_est)==0)&(as.vector(beta_true)==0))], 
     y=as.vector(beta_true)[!((as.vector(beta_est)==0)&(as.vector(beta_true)==0))], pch=16, col="blue",
     xlab = "estimated coefficients", ylab = "true coefficients",
     main = "True vs. estimated coefficients grouped L0glm\n(no preclustering)")


### CLUSTERING ####

# first estimate PSF using blind deconvolution, using rowSums(Y)
# NOTE: using Yagg1 or Yagg2 below also works but is much slower & solution not better
L <- 41 # length of estimated PSF
g0 <- convkernel("gaussian", L, list(sigma = 6))  # initial guess for PSF
system.time(psf_est <- L0glm_blind_deco(y = rowSums(Y), # note: using Yagg1/Yagg2 is slower & gives same fit
                                        family = poisson(identity), 
                                        baseline = NULL, # or "intercept" or "bspline"
                                        L = L,                 # PSF kernel length 
                                        kernel_est = g0, 
                                        kernel_smoothfun = function(y) log_conc_smooth(y, nknots=round(L/5)),
                                        perc = "mean",         # use centroid centering
                                        IC = "gic",
                                        kernel_true = psf_true/sum(psf_true),
                                        plot = TRUE,
                                        tol = 1e-3)$g_hat) # 5s
X_est <- convmatrix(psf_est, n)

### 1.1. WORKING ON POISSON DISTANCE CLUSTERED DATA, cf. Witten 2011 - BEST? ####
# cluster & aggregate
system.time(Yagg1 <- glm_cluster(Y,
                     family   = poisson(identity),
                     method   = "poissondist",
                     rowbins  = 5)$Y_agg) # aggregate together every bin of rowbins=5 subsequent timepoints for extra performance
# 5s  
dim(Yagg1) # 10000     14
matplot(Yagg1^0.5, type="l", lty=1) # looks good

# fit grouped L0glm model using aggregated Yagg data
system.time(L0glm.chol_aggregated <- L0glm.fit(X_est, Yagg1, 
                                    family = poisson(identity), 
                                    intercept = FALSE,
                                    lambda = "gic",
                                    nonnegative = TRUE,
                                    grouped = TRUE,
                                    relaxed = FALSE,
                                    parallelize = F,
                                    solver = "eigen",
                                    method = 1L)) # 10s
selected_features = which(rowSums(L0glm.chol_aggregated$coefficients!=0)!=0)
length(selected_features) # 1020
length(intersect(selected_features, true_features)) # 925 true positives
length(selected_features)-length(intersect(selected_features, true_features)) # 95 false positives

# fit selected features on all mass channels, now using ungrouped L0glm since features are already grouped
system.time(L0glm.chol_full <- L0glm.fit(X_est[,selected_features,drop=F], 
                                         Y, 
                                         intercept = FALSE,
                                               family = poisson(identity), 
                                               lambda = "gic",
                                               nonnegative = TRUE,
                                               grouped = FALSE,
                                               relaxed = FALSE,
                                               # start = matrix(rowMeans(L0glm.chol_aggregated$coefficients[selected_features,,drop=F]), nrow=length(selected_features), ncol=ncol(Y)),
                                               parallelize = FALSE,
                                               solver = "eigen",
                                               method = 1L)) # 12s with parallelize=F
selected_features_full = selected_features[which(rowSums(L0glm.chol_full$coefficients!=0)!=0)]
length(selected_features_full) # 1019
length(intersect(selected_features_full, true_features)) # 925 true positives
length(selected_features_full)-length(intersect(selected_features_full, true_features)) # 94 false positives
beta_est = Matrix(0, nrow=nrow(beta_true), ncol=ncol(beta_true), sparse=T)
beta_est[selected_features,] = L0glm.chol_full$coefficients
rsq = summary(lm(as.vector(beta_true) ~ 0 + as.vector(beta_est)))$r.squared
rsq # 0.9773037
# NOTE: using X_true performance is only very slightly better, i.e. blind
# deconvolution worked almost as well than if we knew the PSF in advance

plot(x=as.vector(beta_est)[!((as.vector(beta_est)==0)&(as.vector(beta_true)==0))], 
     y=as.vector(beta_true)[!((as.vector(beta_est)==0)&(as.vector(beta_true)==0))], pch=16, col="blue",
     xlab = "estimated coefficients", ylab = "true coefficients",
     main = "True vs. estimated coefficients grouped L0glm\n(with Poisson distance preclustering)")


### 1.2. WORKING ON DEVIANCE DISTANCE CLUSTERED DATA ####
# cluster & aggregate
system.time(Yagg2 <- glm_cluster(Y,
                                 family   = poisson(identity),
                                 method   = "deviance",
                                 rowbins  = 5)$Y_agg) # aggregate together every bin of rowbins=5 subsequent timepoints for extra performance
# 3.05s  
dim(Yagg2) # 10000     4
matplot(Yagg2^0.5, type="l", lty=1)

# fit grouped L0glm model using aggregated Yagg data
system.time(L0glm.chol_aggregated <- L0glm.fit(X_est, Yagg2, 
                                               family = poisson(identity), 
                                               lambda = "gic",
                                               nonnegative = TRUE,
                                               grouped = TRUE,
                                               relaxed = FALSE,
                                               parallelize = FALSE,
                                               solver = "eigen",
                                               method = 1L)) # 2.5s
selected_features = which(rowSums(L0glm.chol_aggregated$coefficients!=0)!=0)
length(selected_features) # 1033
length(intersect(selected_features, true_features)) # 903 true positives
length(selected_features)-length(intersect(selected_features, true_features)) # 130 false positives

# fit selected features on all mass channels, now using ungrouped L0glm since features are already grouped
system.time(L0glm.chol_full <- L0glm.fit(X_est[,selected_features,drop=F], 
                                         Y, 
                                         family = poisson(identity), 
                                         lambda = "gic",
                                         nonnegative = TRUE,
                                         grouped = FALSE,
                                         relaxed = FALSE,
                                         parallelize = FALSE,
                                         solver = "eigen",
                                         method = 1L)) # 10s with parallelize=F
selected_features_full = selected_features[which(rowSums(L0glm.chol_full$coefficients!=0)!=0)]
length(selected_features_full) # 1033
length(intersect(selected_features_full, true_features)) # 903 true positives
length(selected_features_full)-length(intersect(selected_features_full, true_features)) # 130 false positives
beta_est = Matrix(0, nrow=nrow(beta_true), ncol=ncol(beta_true), sparse=T)
beta_est[selected_features,] = L0glm.chol_full$coefficients
rsq = summary(lm(as.vector(beta_true) ~ 0 + as.vector(beta_est)))$r.squared
rsq # 0.9601071
# NOTE: performance only slightly lower than if we use X_true

plot(x=as.vector(beta_est)[!((as.vector(beta_est)==0)&(as.vector(beta_true)==0))], 
     y=as.vector(beta_true)[!((as.vector(beta_est)==0)&(as.vector(beta_true)==0))], pch=16, col="blue",
     xlab = "estimated coefficients", ylab = "true coefficients",
     main = "True vs. estimated coefficients grouped L0glm\n(with deviance distance preclustering)")



## 2. TEST WITH ACTUAL GC/MS SAMPLES ####

Y <- GCMS_data1 # or GCMS_data2 
n <- nrow(Y)
weights=1/(Y+0.1) # approx Poisson 1/variance weights
Yw = Y*sqrt(weights) # to take into account approx 1/variance Poisson weights
L2n_cols = sqrt(colSums(Yw^2)) # L2 norm rows (timepoints)
Yw_colnorm = sweep(Yw, 2, L2n_cols+0.1, FUN="/")  # columns (mz) normalized by per-row L2 norm
sqrt(colSums(Yw_colnorm^2)) # ca 1
plotmat(Y, 
        xlab="Elution time (mins)",
        ylab="Mass (m/z)")
matplot(Y^0.4, type="l", col = rgb(0, 0, 0, 0.1), lty = 1)


### CLUSTERING ####

psf_fun <- function(x, centre, width) exp(-0.5 * ((x - centre) / width)^2)
psf <- psf_fun(0:88, centre=44, width=12) # TO DO use kernel as determine by blind deconvolution algo from alkane ladder
plot(psf, type="l")
X <- convmatrix(kernel=psf, n=nrow(Y))
# image(X)

# first estimate PSF using blind deconvolution, using rowSums(Y)
# NOTE: using Yagg1 or Yagg2 below also works but is much slower & solution not better
L <- 91 # length of estimated PSF
g0 <- convkernel("gaussian", L, list(sigma = 10))  # initial guess for PSF
plot(g0, type="l")
system.time(psf_est <- L0glm_blind_deco(rowSums(Y), # note: using Yagg1/Yagg2 is slower & gives same fit
                                        family = poisson(identity), 
                                        baseline = NULL, 
                                        L = L,                 # PSF kernel length 
                                        kernel_est = g0, 
                                        kernel_smoothfun = function(y) log_conc_smooth(y, nknots=round(L/5)),
                                        perc = "mean",         # use centroid centering
                                        IC = "gic",
                                        plot = TRUE,
                                        tol = 1e-3)$g_hat) # 2s
X_est <- convmatrix(psf_est, n)

matplot(Y^0.4, type="l", lty=1)
lines(250*psf^0.4, col="red", lwd=3)
lines(800*psf_est^0.4, col="blue", lwd=3)

### 2.1. WORKING ON POISSON DISTANCE CLUSTERED DATA, cf. Witten 2011 - BEST SOLUTION ####
# cluster & aggregate
system.time(Yagg1 <- glm_cluster(Y,
                                 family   = poisson(identity),
                                 method   = "poissondist",
                                 rowbins  = 5)$Y_agg) # aggregate together every bin of rowbins=5 subsequent timepoints for extra performance
# 1.28s
dim(Yagg1) # 119 3
matplot(Yagg1^0.5, type="l", lty=1) # looks good & distance in line with PoiClaClu::PoissonDistance

# estimate PSF using blind deconvolution, using Yagg1
L <- 91 # length of estimated PSF
g0 <- convkernel("gaussian", L, list(sigma = 10))  # initial guess for PSF
plot(g0, type="l")
system.time(psf_est1 <- L0glm_blind_deco(Yagg1,
                                        family = poisson(identity), 
                                        baseline = NULL,
                                        L = L,                 # PSF kernel length 
                                        kernel_est = g0, 
                                        kernel_smoothfun = function(y) log_conc_smooth(y, nknots=round(L/5)),
                                        perc = "mean",         # use centroid centering
                                        IC = "gic",
                                        grouped = TRUE,
                                        plot = TRUE,
                                        tol = 1e-3)$g_hat) # 2s
X_est <- convmatrix(psf_est1, n) 
# using this PSF estimate in the deconvolution below did not give a good result
# though (too many peaks - 25 instead of the ca 6 expected)
# It seems that deconvoluting with a PSF that is slightly wider works better.
# So we will use our visually fitted Gaussian PSF
X_est <- X


# fit grouped L0glm model using aggregated Yagg data
system.time(L0glm.chol_aggregated <- L0glm.fit(X_est, Yagg1, 
                                               family = poisson(identity), 
                                               intercept = TRUE,
                                               normalize = TRUE,
                                               lambda = "gic",
                                               nonnegative = TRUE,
                                               grouped = TRUE,
                                               relaxed = FALSE,
                                               parallelize = FALSE,
                                               solver = "eigen",
                                               method = 1L)) # 0.31s
selected_features = which(rowSums(L0glm.chol_aggregated$coefficients!=0)!=0)
length(selected_features) # 7
# plotmat(L0glm.chol_aggregated$coefficients, subsamp=1, xlab="row", ylab="col")
matplot(Yagg1^0.5, type="l", lty=1) # looks good & distance in line with PoiClaClu::PoissonDistance
matplot((0.5*X[,selected_features,drop=F]*max(Yagg1))^0.5, type="l", 
        lty=1, lwd=3, add=TRUE)


# fit selected features on all mass channels, now using ungrouped L0glm since features are already grouped
system.time(L0glm.chol_full <- L0glm.fit(X_est[,selected_features,drop=F], 
                                         Y, 
                                         intercept = TRUE,
                                         normalize = TRUE,
                                         family = poisson(identity), 
                                         lambda = "gic",
                                         nonnegative = TRUE,
                                         grouped = FALSE,
                                         relaxed = FALSE,
                                         # start = matrix(rowMeans(L0glm.chol_aggregated$coefficients[selected_features,,drop=F]), nrow=length(selected_features), ncol=ncol(Y)),
                                         parallelize = FALSE,
                                         solver = "eigen",
                                         method = 1L)) # 5s with parallelize=F
selected_features_full = selected_features[which(rowSums(L0glm.chol_full$coefficients!=0)!=0)]
length(selected_features_full) # intercept+7 with grouped=F or T
dim(L0glm.chol_full$coefficients) # (intercept+7) x 671
plotmat(L0glm.chol_full$coefficients)
matplot(Y^0.5, type="l", lty=1, col=rgb(0,0,0,0.1))
matplot((0.5*X[,selected_features,drop=F]*max(Y))^0.5, type="l", 
        lty=1, lwd=3, add=TRUE)

peakmaxima = rownames(Y)[sapply(selected_features, function (peak) which.max(as.vector(X[,peak,drop=F])))]
for (i in (2:nrow(L0glm.chol_full$coefficients))) { # we leave out peak 1 = intercept/background
  title=paste("Peak",as.character(i-1), " RT", peakmaxima[i-1])
  plot(x=as.numeric(colnames(Y)), y=L0glm.chol_full$coefficients[i,],
       main=title,
       ylab="Ion count",
       xlab="Mass (m/z)",
       type="h", col="blue") }
# peak 1=noise (likely edge artefact)
# peak 2=15-Methylnonacosane 224+225 RT 26.86
# peak 3=13-Methylnonacosane 196+252 RT 26.91
# peak 4=Cyclodecasiloxane, eicosamethyl 73+147+221+281+355++429+503 RT 26.95 + 168 from 11-Methylnonacosane
# peak 5=9-Methylnonacosane 140+308 RT 27.01
# peak 6=noise (low ion intensity) (likely edge artefact) RT 27.107
# peak 7=7-Methylnonacosane 112+336 RT 27.113
# BEST SOLUTION
# see
?GCMS_data1

W = cbind(1,X[,selected_features,drop=F])
H = L0glm.chol_full$coefficients
Y_est = W %*% H
Y_res = Y-Y_est # residual
range(Y_res)
range(Y)
colnames(Y_est) = colnames(Y)
rownames(Y_est) = rownames(Y)
plotmat(Y_est, main="reconstructed signal")
plotmat(Y, main="true signal")
plotmat(Y-Y_est, main="residual")



### 2.2. WORKING ON DEVIANCE DISTANCE CLUSTERED DATA - SLIGHTLY WORSE SOLUTION ####
# cluster & aggregate
system.time(Yagg2 <- glm_cluster(Y,
                                 family   = poisson(identity),
                                 method   = "deviance",
                                 rowbins  = 5)$Y_agg) # aggregate together every bin of rowbins=5 subsequent timepoints for extra performance
dim(Yagg2) # 119     4
matplot(Yagg2, type="l", lty=1)

# fit grouped L0glm model using aggregated Yagg data
system.time(L0glm.chol_aggregated <- L0glm.fit(X, Yagg2, 
                                               intercept = TRUE,
                                               normalize = TRUE,
                                               family = poisson(identity), 
                                               lambda = "gic",
                                               nonnegative = TRUE,
                                               grouped = TRUE,
                                               relaxed = FALSE,
                                               parallelize = FALSE,
                                               solver = "eigen",
                                               method = 1L)) # 0.19s
selected_features = which(rowSums(L0glm.chol_aggregated$coefficients!=0)!=0)
length(selected_features) # 7
# plotmat(L0glm.chol_aggregated$coefficients, subsamp=1, xlab="row", ylab="col")
matplot(Yagg1^0.5, type="l", lty=1) # looks good & distance in line with PoiClaClu::PoissonDistance
matplot((0.5*X[,selected_features,drop=F]*max(Yagg1))^0.5, type="l", 
        lty=1, lwd=3, add=TRUE)

# fit selected features on all mass channels, now using ungrouped L0glm since features are already grouped
system.time(L0glm.chol_full <- L0glm.fit(X[,selected_features,drop=F], 
                                         Y, 
                                         intercept = TRUE,
                                         normalize = TRUE,
                                         family = poisson(identity), 
                                         lambda = "gic",
                                         nonnegative = TRUE,
                                         grouped = FALSE,
                                         relaxed = FALSE,
                                         parallelize = FALSE,
                                         solver = "eigen",
                                         method = 1L)) # 6.8s with parallelize=F
selected_features_full = selected_features[which(rowSums(L0glm.chol_full$coefficients!=0)!=0)]
length(selected_features_full) # intercept+7, with grouped=F or T
dim(L0glm.chol_full$coefficients) # (intercept+7) x 671
plotmat(L0glm.chol_full$coefficients)
matplot(Y^0.5, type="l", lty=1, col=rgb(0,0,0,0.1))
matplot((0.5*X[,selected_features,drop=F]*max(Y))^0.5, type="l", 
        lty=1, lwd=3, add=TRUE)

peakmaxima = rownames(Y)[sapply(selected_features, function (peak) which.max(as.vector(X[,peak,drop=F])))]
for (i in (2:nrow(L0glm.chol_full$coefficients))) { # we leave out peak 1 = intercept/background
  title=paste("Peak",as.character(i-1), " RT", peakmaxima[i-1])
  plot(x=as.numeric(colnames(Y)), y=L0glm.chol_full$coefficients[i,],
       main=title,
       ylab="Ion count",
       xlab="Mass (m/z)",
       type="h", col="blue") }
# peaks 1 & 2=noise (due to edge artefact)
# peak 3=15-Methylnonacosane 224+225 RT 26.87
# peak 4=13-Methylnonacosane 196+252 RT 26.92
# peak 5=Cyclodecasiloxane, eicosamethyl 73+147+221+281+355++429+503 RT 26.94 + 168 from 11-Methylnonacosane
# peak 6=9-Methylnonacosane 140+308 RT 27.01
# peak 7=7-Methylnonacosane 112+336 RT 27.12
# LESS GOOD SOLUTION
# see
?GCMS_data1

W = cbind(1,X[,selected_features,drop=F])
H = L0glm.chol_full$coefficients
Y_est = W %*% H
Y_res = Y-Y_est # residual
range(Y_res)
range(Y)
colnames(Y_est) = colnames(Y)
rownames(Y_est) = rownames(Y)
plotmat(Y_est, subsamp=1, main="reconstructed signal")
plotmat(Y, subsamp=1, main="true signal")
plotmat(Y-Y_est, subsamp=1, main="residual")
