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
dat=simulate_spike_train(n=10000, p=10000, 
                         k=1000, # peaks
                         family="poisson", multichannel=T)
Y=dat$y
X=dat$X
dim(Y) # 10000 timepoints x 228 mass channels
true_features = which(rowSums(dat$beta_true!=0)!=0)
beta_true = dat$beta_true
length(true_features)
matplot(Y^0.5, type="l", lty=1, col=rgb(0,0,0,0.1))
plot(rowSums(Y)^0.5, type="l", lty=1)

# NOTE: without clustering and aggregating mass channels first
# a single grouped L0glm fit using a sparse Cholesky solver takes 106s :
system.time(L0glm.chol_originaldata <- L0glm.fit(X, 
                                         Y, 
                                         family = poisson(identity), 
                                         lambda = "gic",
                                         nonnegative = TRUE,
                                         grouped = TRUE,
                                         relaxed = FALSE,
                                         warm.start = 1,
                                         parallelize = FALSE,
                                         solver = "eigen",
                                         method = 7L)) # 87s
selected_features = which(rowSums(L0glm.chol_originaldata$coefficients!=0)!=0)
length(selected_features) # 950
length(intersect(selected_features, true_features)) # 929 true positives
length(selected_features)-length(intersect(selected_features, true_features)) # 21 false positives
beta_est = L0glm.chol_originaldata$coefficients
rsq = summary(lm(as.vector(beta_true) ~ 0 + as.vector(beta_est)))$r.squared
rsq # 0.9912648

plot(x=as.vector(beta_est)[!((as.vector(beta_est)==0)&(as.vector(beta_true)==0))], 
     y=as.vector(beta_true)[!((as.vector(beta_est)==0)&(as.vector(beta_true)==0))], pch=16, col="blue",
     xlab = "estimated coefficients", ylab = "true coefficients",
     main = "True vs. estimated coefficients grouped L0glm\n(no preclustering)")


### CLUSTERING ####
# for clustering we aggregate together every bin of rowbins=5 subsequent timepoints for extra performance:
rowbins <- 5
Ybinned <- t(glm_aggregate(Y = t(Y),
                           weights = matrix(1, nrow=ncol(Y), ncol=nrow(Y)), # or t(weights) if given as argument
                           family = poisson(identity), 
                           clusterids = rep(1L:ceiling(nrow(Y)/rowbins), each=rowbins)[1:nrow(Y)])$Y_agg)

### 1.1. WORKING ON POISSON DISTANCE CLUSTERED DATA, cf. Witten 2011 - BEST? ####
# cluster & aggregate
system.time(Ybinnedcl1 <- glm_cluster(Ybinned, family=poisson(identity), method="poissondist")) # 5s
cl1 <- Ybinnedcl1$clusterids # obtained clusters
Yagg1 <- glm_aggregate(Y, family=poisson(identity), clusterids = cl1)$Y_agg # aggregate original data within these clusters
dim(Yagg1) # 10000     14
matplot(Yagg1^0.5, type="l", lty=1) # looks good & distance in line with PoiClaClu::PoissonDistance

# fit grouped L0glm model using aggregated Yagg data
system.time(L0glm.chol_aggregated <- L0glm.fit(X, Yagg1, 
                                    family = poisson(identity), 
                                    intercept = FALSE,
                                    lambda = "gic",
                                    nonnegative = TRUE,
                                    grouped = TRUE,
                                    relaxed = FALSE,
                                    warm.start = 1,
                                    parallelize = F,
                                    solver = "eigen",
                                    method = 7L)) # 6s
selected_features = which(rowSums(L0glm.chol_aggregated$coefficients!=0)!=0)
length(selected_features) # 1017
length(intersect(selected_features, true_features)) # 930 true positives
length(selected_features)-length(intersect(selected_features, true_features)) # 87 false positives

# fit selected features on all mass channels, now using ungrouped L0glm since features are already grouped
system.time(L0glm.chol_full <- L0glm.fit(X[,selected_features,drop=F], 
                                         Y, 
                                         intercept = FALSE,
                                               family = poisson(identity), 
                                               lambda = "gic",
                                               nonnegative = TRUE,
                                               grouped = FALSE,
                                               relaxed = FALSE,
                                               # start = matrix(rowMeans(L0glm.chol_aggregated$coefficients[selected_features,,drop=F]), nrow=length(selected_features), ncol=ncol(Y)),
                                               warm.start = 1,
                                               parallelize = FALSE,
                                               solver = "eigen",
                                               method = 7L)) # 9s with parallelize=F
selected_features_full = selected_features[which(rowSums(L0glm.chol_full$coefficients!=0)!=0)]
length(selected_features_full) # 1017
length(intersect(selected_features_full, true_features)) # 930 true positives
length(selected_features_full)-length(intersect(selected_features_full, true_features)) # 87 false positives
beta_est = Matrix(0, nrow=nrow(beta_true), ncol=ncol(beta_true), sparse=T)
beta_est[selected_features,] = L0glm.chol_full$coefficients
rsq = summary(lm(as.vector(beta_true) ~ 0 + as.vector(beta_est)))$r.squared
rsq # 0.9806827
# NOTE: using grouped=TRUE we get 911 true pos, 47 false pos & rsq 0.9850147

plot(x=as.vector(beta_est)[!((as.vector(beta_est)==0)&(as.vector(beta_true)==0))], 
     y=as.vector(beta_true)[!((as.vector(beta_est)==0)&(as.vector(beta_true)==0))], pch=16, col="blue",
     xlab = "estimated coefficients", ylab = "true coefficients",
     main = "True vs. estimated coefficients grouped L0glm\n(with Poisson distance preclustering)")


### 1.2. WORKING ON DEVIANCE DISTANCE CLUSTERED DATA ####
# cluster & aggregate
system.time(Ybinnedcl2 <- glm_cluster(Ybinned, family=poisson(identity), method="deviance")) # 3s
cl2 <- Ybinnedcl2$cluster # obtained clusters
Yagg2 <- glm_aggregate(Y, family=poisson(identity), clusterids = cl2)$Y_agg # aggregate original data within these clusters
dim(Yagg2) # 10000     4
matplot(Yagg2, type="l", lty=1)

# fit grouped L0glm model using aggregated Yagg data
system.time(L0glm.chol_aggregated <- L0glm.fit(X, Yagg2, 
                                               family = poisson(identity), 
                                               lambda = "gic",
                                               nonnegative = TRUE,
                                               grouped = TRUE,
                                               relaxed = FALSE,
                                               warm.start = 1,
                                               parallelize = FALSE,
                                               solver = "eigen",
                                               method = 7L)) # 2.5s
selected_features = which(rowSums(L0glm.chol_aggregated$coefficients!=0)!=0)
length(selected_features) # 1028
length(intersect(selected_features, true_features)) # 917 true positives
length(selected_features)-length(intersect(selected_features, true_features)) # 111 false positives

# fit selected features on all mass channels, now using ungrouped L0glm since features are already grouped
system.time(L0glm.chol_full <- L0glm.fit(X[,selected_features,drop=F], 
                                         Y, 
                                         family = poisson(identity), 
                                         lambda = "gic",
                                         nonnegative = TRUE,
                                         grouped = FALSE,
                                         relaxed = FALSE,
                                         warm.start = 1,
                                         parallelize = FALSE,
                                         solver = "eigen",
                                         method = 7L)) # 10s with parallelize=F
selected_features_full = selected_features[which(rowSums(L0glm.chol_full$coefficients!=0)!=0)]
length(selected_features_full) # 1028
length(intersect(selected_features_full, true_features)) # 917 true positives
length(selected_features_full)-length(intersect(selected_features_full, true_features)) # 111 false positives
beta_est = Matrix(0, nrow=nrow(beta_true), ncol=ncol(beta_true), sparse=T)
beta_est[selected_features,] = L0glm.chol_full$coefficients
rsq = summary(lm(as.vector(beta_true) ~ 0 + as.vector(beta_est)))$r.squared
rsq # 0.9659199
# NOTE: using grouped=TRUE we get 901 true pos, 64 false pos & rsq 0.9719609

plot(x=as.vector(beta_est)[!((as.vector(beta_est)==0)&(as.vector(beta_true)==0))], 
     y=as.vector(beta_true)[!((as.vector(beta_est)==0)&(as.vector(beta_true)==0))], pch=16, col="blue",
     xlab = "estimated coefficients", ylab = "true coefficients",
     main = "True vs. estimated coefficients grouped L0glm\n(with deviance distance preclustering)")



## 2. TEST WITH ACTUAL GC/MS SAMPLES ####

# function to plot matrix as image, also taking sparse matrix as input
plotmat <- function(M, 
                    row_subsamp=round(nrow(M)/min(nrow(M),2000)), 
                    col_subsamp=round(ncol(M)/min(ncol(M),2000)),
                    gamma=0.12, ...) { 
  cols = colorRampPalette(c("black", "black","blue","green","yellow","orange", "red"))(1000)
  M = M[seq(1,nrow(M),by=row_subsamp),
        seq(1,ncol(M),by=col_subsamp), drop=FALSE]
  M=as.matrix(M) # convert sparse to dense matrix
  M[is.na(M)]=0
  M[M<0]=0
  if (!is.null(as.numeric(rownames(M))[[1]])) { xrange = range(as.numeric(rownames(M)))
                                                xvals = seq(min(xrange), max(xrange), length.out = nrow(M)) } else { xvals = 1:nrow(M) }
  if (!is.null(as.numeric(colnames(M))[[1]])) { yrange = range(as.numeric(colnames(M)))
                                                yvals = seq(min(yrange), max(yrange), length.out = ncol(M)) } else { yvals = 1:ncol(M) }
  image(x=xvals,
        y=yvals,
        z=M^gamma, 
        col = cols,
        useRaster=TRUE,
        ...)
}

# load data
Y <- GCMS_data1 # or GCMS_data2 
weights=1/(Y+0.1) # approx Poisson 1/variance weights
Yw = Y*sqrt(weights) # to take into account approx 1/variance Poisson weights
L2n_cols = sqrt(colSums(Yw^2)) # L2 norm rows (timepoints)
Yw_colnorm = sweep(Yw, 2, L2n_cols+0.1, FUN="/")  # columns (mz) normalized by per-row L2 norm
sqrt(colSums(Yw_colnorm^2)) # ca 1
plotmat(Y, 
        xlab="Retention time (mins)",
        ylab="Mass (m/z)")
matplot(Y^0.4, type="l", col = rgb(0, 0, 0, 0.1), lty = 1)

psf_fun <- function(x, centre, width) exp(-0.5 * ((x - centre) / width)^2)
psf <- psf_fun(0:88, centre=44, width=12) # TO DO use kernel as determine by blind deconvolution algo from alkane ladder
X <- convmatrix(kernel=psf, n=nrow(Y))
# image(X)
matplot(Y^0.4, type="l", lty=1)
lines(250*psf^0.4, col="red", lwd=3)

# NOTE: without clustering and aggregating mass channels first
# a single grouped L0glm fit using a sparse Cholesky solver takes 64s and
# returns a fairly high nr of false positive peaks :
system.time(L0glm.chol_originaldata <- L0glm.fit(X, 
                                                 Y, 
                                                 intercept=FALSE,
                                                 family = poisson(identity), 
                                                 lambda = "gic",
                                                 nonnegative = TRUE,
                                                 grouped = TRUE,
                                                 relaxed = FALSE,
                                                 warm.start = 1,
                                                 parallelize = FALSE,
                                                 solver = "eigen",
                                                 method = 7L)) # 57s
selected_features = which(rowSums(L0glm.chol_originaldata$coefficients!=0)!=0)
length(selected_features) # 17 - quite a few false positives
# plotmat(L0glm.chol_originaldata$coefficients, subsamp=1, xlab="row", ylab="col")
matplot(Y^0.5, type="l", lty=1)
matplot((0.5*X[,selected_features,drop=F]*max(Y))^0.5, type="l", 
        lty=1, lwd=5, add=TRUE)

peakmaxima = rownames(Y)[sapply(selected_features, function (peak) which.max(as.vector(X[,peak,drop=F])))]
for (i in (1:length(selected_features))) { 
  title=paste("Peak",as.character(i-1), " RT", peakmaxima[i-1])
  plot(x=as.numeric(colnames(Y)), y=L0glm.chol_originaldata$coefficients[selected_features[i],],
       main=title,
       ylab="Ion count",
       xlab="Mass (m/z)",
       type="h", col="blue") }
# peak 1=noise
# peak 2+3=15-Methylnonacosane 224+225 RT 26.877 (but split)
# peak 4+5+6=13-Methylnonacosane 196+252 RT 26.91 (but split)
# peak 7=Cyclodecasiloxane, eicosamethyl 73+147+221+281+355++429+503 RT 26.94 + 168 from 11-Methylnonacosane
# peak 8+9=9-Methylnonacosane 140+308 RT 27.01 (but split)
# peak 10+11=noise
# peak 12=7-Methylnonacosane 112+336 RT 27.11
# peaks 13-16=noise
# many false positives :-(


### CLUSTERING ####
# for clustering we aggregate together every bin of rowbins subsequent timepoints for extra performance:
rowbins <- 5
Ybinned <- t(glm_aggregate(Y = t(Y), 
                           weights = matrix(1, nrow=ncol(Y), ncol=nrow(Y)), # or t(weights) if given as argument
                           family = poisson(identity), 
                           clusterids = rep(1L:ceiling(nrow(Y)/rowbins), each=rowbins)[1:nrow(Y)])$Y_agg)

### 2.1. WORKING ON POISSON DISTANCE CLUSTERED DATA, cf. Witten 2011 - BEST SOLUTION ####
# cluster & aggregate
system.time(Ybinnedcl1 <- glm_cluster(Ybinned, family=poisson(identity), method="poissondist")) # 1.3s
cl1 <- Ybinnedcl1$cluster # obtained clusters
Yagg1 <- glm_aggregate(Y, family=poisson(identity), clusterids = cl1)$Y_agg # aggregate original data within these clusters
dim(Yagg1) # 119 3
matplot(Yagg1^0.5, type="l", lty=1) # looks good & distance in line with PoiClaClu::PoissonDistance

# fit grouped L0glm model using aggregated Yagg data
system.time(L0glm.chol_aggregated <- L0glm.fit(X, Yagg1, 
                                               family = poisson(identity), 
                                               intercept = FALSE,
                                               lambda = "gic",
                                               nonnegative = TRUE,
                                               grouped = TRUE,
                                               relaxed = FALSE,
                                               warm.start = 1,
                                               parallelize = FALSE,
                                               solver = "eigen",
                                               method = 7L)) # 0.18s
selected_features = which(rowSums(L0glm.chol_aggregated$coefficients!=0)!=0)
length(selected_features) # 6
# plotmat(L0glm.chol_aggregated$coefficients, subsamp=1, xlab="row", ylab="col")
matplot(Yagg1^0.5, type="l", lty=1) # looks good & distance in line with PoiClaClu::PoissonDistance
matplot((0.5*X[,selected_features,drop=F]*max(Yagg1))^0.5, type="l", 
        lty=1, lwd=3, add=TRUE)


# fit selected features on all mass channels, now using ungrouped L0glm since features are already grouped
system.time(L0glm.chol_full <- L0glm.fit(X[,selected_features,drop=F], 
                                         Y, 
                                         intercept = TRUE,
                                         family = poisson(identity), 
                                         lambda = "gic",
                                         nonnegative = TRUE,
                                         grouped = FALSE,
                                         relaxed = FALSE,
                                         # start = matrix(rowMeans(L0glm.chol_aggregated$coefficients[selected_features,,drop=F]), nrow=length(selected_features), ncol=ncol(Y)),
                                         warm.start = 1,
                                         parallelize = FALSE,
                                         solver = "eigen",
                                         method = 7L)) # 5s with parallelize=F
selected_features_full = selected_features[which(rowSums(L0glm.chol_full$coefficients!=0)!=0)]
length(selected_features_full) # intercept+6 with grouped=F or T
dim(L0glm.chol_full$coefficients) # (intercept+6) x 671
image(L0glm.chol_full$coefficients)
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
# peak 1=15-Methylnonacosane 224+225 RT 26.877
# peak 2=13-Methylnonacosane 196+252 RT 26.91
# peak 3=Cyclodecasiloxane, eicosamethyl 73+147+221+281+355++429+503 RT 26.94 + 168 from 11-Methylnonacosane
# peak 4=9-Methylnonacosane 140+308 RT 27.01
# peak 5=7-Methylnonacosane 112+336 RT 27.11
# peak 6=noise (low ion intensity) RT 27.20
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
system.time(Ybinnedcl2 <- glm_cluster(Ybinned, family=poisson(identity), method="deviance")) # 4s
cl2 <- Ybinnedcl2$cluster # obtained clusters
Yagg2 <- glm_aggregate(Y, family=poisson(identity), clusterids = cl2)$Y_agg # aggregate original data within these clusters
dim(Yagg2) # 119     4
matplot(Yagg2, type="l", lty=1)

# fit grouped L0glm model using aggregated Yagg data
system.time(L0glm.chol_aggregated <- L0glm.fit(X, Yagg2, 
                                               intercept = FALSE,
                                               family = poisson(identity), 
                                               lambda = "gic",
                                               nonnegative = TRUE,
                                               grouped = TRUE,
                                               relaxed = FALSE,
                                               warm.start = 1,
                                               parallelize = FALSE,
                                               solver = "eigen",
                                               method = 7L)) # 0.19s
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
                                         family = poisson(identity), 
                                         lambda = "gic",
                                         nonnegative = TRUE,
                                         grouped = FALSE,
                                         relaxed = FALSE,
                                         warm.start = 1,
                                         parallelize = FALSE,
                                         solver = "eigen",
                                         method = 7L)) # 6.8s with parallelize=F
selected_features_full = selected_features[which(rowSums(L0glm.chol_full$coefficients!=0)!=0)]
length(selected_features_full) # intercept+7, with grouped=F or T
dim(L0glm.chol_full$coefficients) # (intercept+7) x 671
image(L0glm.chol_full$coefficients)
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
# peak 1=noise (due to edge artefact)
# peak 2=15-Methylnonacosane 224+225 RT 26.88
# peak 3=13-Methylnonacosane 196+252 RT 26.91
# peak 4=Cyclodecasiloxane, eicosamethyl 73+147+221+281+355++429+503 RT 26.94 + 168 from 11-Methylnonacosane
# peak 5=9-Methylnonacosane 140+308 RT 27.01
# peak 6=noise (very low ion count) (likely edge artefact)
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
