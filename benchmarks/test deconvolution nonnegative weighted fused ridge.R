# DECONVOLUTION USING NONNEGATIVE WEIGHTED FUSED RIDGE REGRESSION WITH 1st & 2ND ORDER DIFFERENCE RIDGE PENALTY
# T. Wenseleers, 28 febr. 2021

# couple of utility functions ####
# make a square n x n banded convolution matrix containing time-shifted copies of convolution kernel K
H = function (n, K) { 
  require(Matrix)
  return(as.matrix(bandSparse(n, 
                              k = seq(-(length(K)-1),0,1), 
                              diag = t(replicate(n, rev(K))), symm=FALSE)))
  }

# linear convolution of vector y with convolution kernel K
linconv = function (y, K) {   K=K/sum(K)
                             # out=convolve(y, K, type="open") # convolution using default R stats function
                             # ntoclip=(length(out)-length(y)) # /2
                             # out=head(out,-ntoclip )
                             # return(out)
                               X = H(length(y), K)
                               out = X %*% as.matrix(y, ncol=1) # convolution by multiplying with banded convolution matrix
                               return(out)
}

# circular shift of a vector
circshift = function(x, n = 0) {
  if (n == 0) x else c(tail(x, -n), head(x, n))
}

# simulate some data to test accuracy of deconvolution ####

# without a weekday effect ####
K = dexp(seq(1, 41, by=1), rate=1/(5.2+2)) # exponential convolution kernel to go from infections to cases
K = K / sum(K)
plot(K, type="l", main="Convolution kernel")

t = (1-length(K)):1000 # time
# inf_real = 200*exp(sin(t*(0.02+t/10000)))
inf_real = 200*exp(sin(t*0.03))-140 # simulated real infections
inf_real[inf_real<0] = 0
cases_real = tail(linconv(inf_real, K), -length(K)) # simulated diagnosed cases
t = tail(t, -length(K))
inf_real = tail(inf_real, -length(K))
cases_obs = as.numeric(rpois(length(cases_real), lambda=cases_real)) # simulated diagnosed cases with poisson noise added
# we extrapolate input signal a little bit to avoid edge effects in deconvolution
sm = smooth.spline(x=t, y=cases_obs, w=1/(cases_obs+0.1))
cases_obs_extrapol = c(cases_obs, predict(sm, x=c((max(t)+1):(max(t)+length(K))))$y) 

library(Hmisc)
subs = 100:800
c(-50:50)[which.min(sapply(X=-50:50, FUN = function(shift) sum((as.vector(inf_real)[subs]-Lag(as.vector(cases_real)[subs], 
                                                                                              shift=shift))^2, na.rm=T)))]
# -6, as expected

# deconvolution using nonnegative fused weighted ridge regression with 1st & 2nd order difference ridge penalty on the coefficients
y_obs = cases_obs_extrapol  
library(Matrix)
lambda1 = 0.005
lambda2 = 3000
lambdas1 = lambda1/norm(abs(diff(y_obs, differences=1)), "2") # 0.00001 # regularization parameter for 1st order difference ridge penalty
lambdas2 = lambda2/norm(abs(diff(y_obs, differences=2)), "2") #  1000000 # regulatization parameter for 2nd order difference ridge penalty

# PS: one could also contemplate using an adaptive penalty, similar to https://www.mghassany.com/bibliography/files/ghassany-master-thesis.pdf
# quick attempt here, but didn't work...
# y_expected = y_obs
# lambdas1 = 0.001*1/c(mean(abs(diff(y_expected, differences=1))),abs(diff(y_expected, differences=1))+1) # 0.00001 # regularization parameter for 1st order difference ridge penalty
# lambdas2 = 1*1/c(rep(mean(abs(diff(y_expected, differences=2))),2),abs(diff(y_expected, differences=2))+1) #  1000000 # regulatization parameter for 2nd order difference ridge penalty

X = H(length(y_obs), K) # banded convolution matrix containing time-shifted copies of the convolution kernel K
X = sweep(X, 1, STATS=rowSums(X), FUN="/") # make all rowSums sum to 1, also at edges (this fixed edge effect a bit)
nc = ncol(X)
# we add 1st & 2nd order fused ridge regularization by row augmenting our banded convolution covariate matrix
Dg1 = bandSparse(nc, k = c(-1,0,1), diag = list(rep(-1,nc)*sqrt(lambdas1), # 1st order difference ridge penalty, i.e. "fused ridge"
                                                rep(1,nc)*sqrt(lambdas1),
                                                rep(0,nc)*sqrt(lambdas1)), symm=FALSE)
Dg2 = bandSparse(nc, k = c(-1,0,1), diag = list(rep(1,nc)*sqrt(lambdas2), # 2nd order difference ridge penalty
                                                rep(-2,nc)*sqrt(lambdas2),
                                                rep(1,nc)*sqrt(lambdas2)), symm=FALSE)
X_augm_fusedridge12 = as.matrix(rbind(X, Dg1, Dg2)) # augmented covariate matrix to use for deconvolution
y_augm12 = c(y_obs, rep(0,nc*2)) # we also augment our observed counts with zeros
w12 = 1/c(y_obs+0.1, rep(1,2*nc)) 
# 1/variance Poisson observation weights, equivalent to doing a 1 iteration IRLS to approximate an identity link Poisson GLM 
# (cf http://bwlewis.github.io/GLM/, https://stats.stackexchange.com/questions/405629/are-glms-just-glorified-wls-regressions), 
# here with nonnegativity constraints by using nnls instead of regular least squares
# one could do several iterations using the expected value of y to get a more accurate results if desired
library(nnls)
system.time(nonnegfusedridgefit12 <- nnls(X_augm_fusedridge12*sqrt(w12), y_augm12*sqrt(w12))) 
# 2.86s
# PS maybe worth checking out this nnls implementation using Forward Active Set Truncation/FAST NNLS
# and using RcppArmadillo, so that it is also optimized for a sparse covariate matrix, 
# https://stackoverflow.com/questions/58006606/rcpp-implementation-of-fast-nonnegative-least-squares
nonnegfusedridgefit12_coef = head(nonnegfusedridgefit12$x, -length(K)) # deconvoluted signal, ie case counts deconvoluted to infections


# plot deconvolution result
plot(t, inf_real, pch=16, type="l", ylab="Count", 
     main="Real infections (green), observed cases (red), infections estimated using\ndeconvolution with nonnegative weighted fused ridge regression\nwith 1st & 2nd order difference ridge penalty (blue)",
     cex.main=0.8, lwd=3, col="green3", ylim=c(0,550)) # theoretical infections without any noise
lines(t, cases_real, col="red", lwd=2) # theoretical case counts without any noise
points(t, cases_obs, col="red", pch=16, cex=0.5) # observed case counts
lines(t, nonnegfusedridgefit12_coef, col="blue", lwd=2) # case counts deconvoluted to infections (blue & green should coincide)
# lines(t, y_expected, col="darkgrey", lwd=3) # model predictions
y_expected = (X_augm_fusedridge12 %*% as.matrix(nonnegfusedridgefit12$x, ncol=1))[1:nc]

# 0, almost perfect fit - as expected

# TO DO: use central difference instead of forward or backward difference? (see diff.zoo)
# TO DO: deal with weekday effect as additional unknown convolution causing varying lags between infection & 
# diagnosis depending on day of the week (e.g. longer if symptom onset is in weekend)
# TO DO: auto-tune regularization parameters using cross-validation 


# deconvolution with L0glm, using fused ridge penalties as above ####

devtools::load_all(".")
family = gaussian(identity) # TODO fit gives error with poisson(identity) or when using weights
# L0glm fit
lambda_preset = preset.lambda(IC = "aic", p=ncol(X), n=nrow(X)) # options aic, bic, gic, ebic, mbic 
system.time(L0glmfit <- L0glm.fit(X=Matrix(X_augm_fusedridge12, sparse=T), 
                                  y=y_augm12, 
                                  weights=w12,
                                  family = family, 
                                  lambda = 0*lambda_preset, # *L0overLqnorm*4*(1/2^(2-q)), # when using q = 1 in adaptive penalties this would have to be lambda_preset*L0overL1norm to be roughly on right scale, with L0overLqnorm calculated from L0glm fit with q=0? 
                                  # ridge.center = start, # ridge.center
                                  nonnegative = TRUE,
                                  start = NULL, # rep(50, ncol(X)), # NULL, # rep(50, ncol(sim$X)), # 10-100 works best, coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                                  solver = "eigen", # note: solver osqp & glmnet better than eigen for very small lambda values)
                                  maxit = 1,
                                  # dispersion = 1,
                                  gamma = 2,
                                  # tol = 0.1,
                                  method = 7L, # 0.27s for n=p=10000 with eigen solver methods 7 & 8 (Cholesky), 0.61s with methods 5 & 6 (CG), 1.27s with methods 3, 4 (LSCG), 219s with methods 1 & 2 (QR), 4.8s for n=p=10000 with osqp solver, 2s for n=p=10000 with glmnet solver
                                  thresh = 1E-5
)) # 0.05s
L0glm_coefs <- head(L0glmfit$coefficients, -length(K))
sum(L0glm_coefs!=0)  # 1000
# plot(beta_true, col="red", type="l")
plot(L0glm_coefs, col=alpha("blue", 0.5), pch=16, type="l")
# this is without any L0 penalty


# plot deconvolution result
plot(t, inf_real, pch=16, type="l", ylab="Count", 
     main="Real infections (green), observed cases (red), infections estimated using\ndeconvolution with nonnegative weighted fused ridge regression\nwith 1st & 2nd order difference ridge penalty (blue)",
     cex.main=0.8, lwd=3, col="green3", ylim=c(0,550)) # theoretical infections without any noise
lines(t, cases_real, col="red", lwd=2) # theoretical case counts without any noise
points(t, cases_obs, col="red", pch=16, cex=0.5) # observed case counts
lines(t, L0glm_coefs, col="blue", lwd=2) # case counts deconvoluted to infections (blue & green should coincide)
# lines(t, y_expected, col="darkgrey", lwd=3) # model predictions
y_expected = (X_augm_fusedridge12 %*% as.matrix(L0glmfit$coefficients, ncol=1))[1:nc]


# do.call(cbind, get.metrics(beta = L0glm_coefs,  # some metrics
#                            X = X, y = y, 
#                            family = family, 
#                            beta_true = beta_true))





# simulated data with a weekday effect
K = dexp(seq(1, 41, by=1), rate=1/(5.2+2)) # exponential convolution kernel to go from infections to cases
K = K / sum(K)
plot(K, type="l", main="Convolution kernel")

# set up week-day effect convolution matrix
t = (1-length(K)):1000 # time
w = rep(factor(c("Saturday","Sunday","Monday","Tuesday","Wednesday","Thursday","Friday"), 
               levels=c("Saturday","Sunday","Monday","Tuesday","Wednesday","Thursday","Friday")),
        ceiling(length(t)/7))[1:length(t)] # weekday
W = model.matrix(~0+., data=data.frame(w=w)) # weekday effect model matrix
K_W = t(rbind( # weekday-specific shifts, i.e. week-day specific convolution kernels
        c(0, 0, 1, 0, 0, 0, 0), # symptomatic on Saturday = diagnosed on Monday
        c(0, 1, 0, 0, 0, 0, 0), # symptomatic on Sunday = diagnosed on Monday
        c(1, 0, 0, 0, 0, 0, 0),
        c(1, 0, 0, 0, 0, 0, 0),
        c(1, 0, 0, 0, 0, 0, 0),
        c(1, 0, 0, 0, 0, 0, 0),
        c(1, 0, 0, 0, 0, 0, 0)
      ))  
H_W = K_W %*% t(W) # convolution matrix to capture weekday-specific shifts
H_W = do.call(rbind, lapply(1:ceiling(length(t)/7), function (i) H_W))[1:length(t),]
dim(H_W)
# PS maybe more efficient using bandSparse to make convolution matrix

# inf_real = 200*exp(sin(t*(0.02+t/10000)))
inf_real = 200*exp(sin(t*0.03)) # simulated real infections
cases_real = H_W %*% as.matrix(inf_real, ncol=1) # add weekday effect
plot(t, cases_real, type="l")
dim( H(length(inf_real), K) ) # 1041 x 1041

cases_real = tail(linconv(inf_real, K), -length(K)) # simulated diagnosed cases without any weekday effect

t = tail(t, -length(K))
inf_real = tail(inf_real, -length(K))
cases_obs = as.numeric(rpois(length(cases_real), lambda=cases_real)) # simulated diagnosed cases with poisson noise added
# we extrapolate input signal a little bit to avoid edge effects in deconvolution
sm = smooth.spline(x=t, y=cases_obs, w=1/(cases_obs+0.1))
cases_obs_extrapol = c(cases_obs, predict(sm, x=c((max(t)+1):(max(t)+length(K))))$y) 

library(Hmisc)
subs = 100:800
c(-50:50)[which.min(sapply(X=-50:50, FUN = function(shift) sum((as.vector(inf_real)[subs]-Lag(as.vector(cases_real)[subs], 
                                                                                              shift=shift))^2, na.rm=T)))]
# -6, as expected

# deconvolution using nonnegative fused weighted ridge regression with 1st & 2nd order difference ridge penalty on the coefficients
y_obs = cases_obs_extrapol  
library(Matrix)
lambdas1 = 0.005/norm(abs(diff(y_obs, differences=1)), "2") # 0.00001 # regularization parameter for 1st order difference ridge penalty
lambdas2 = 3000/norm(abs(diff(y_obs, differences=2)), "2") #  1000000 # regulatization parameter for 2nd order difference ridge penalty

# PS: one could also contemplate using an adaptive penalty, similar to https://www.mghassany.com/bibliography/files/ghassany-master-thesis.pdf
# quick attempt here, but didn't work...
# y_expected = y_obs
# lambdas1 = 0.001*1/c(mean(abs(diff(y_expected, differences=1))),abs(diff(y_expected, differences=1))+1) # 0.00001 # regularization parameter for 1st order difference ridge penalty
# lambdas2 = 1*1/c(rep(mean(abs(diff(y_expected, differences=2))),2),abs(diff(y_expected, differences=2))+1) #  1000000 # regulatization parameter for 2nd order difference ridge penalty

X = H(length(y_obs), K) # banded convolution matrix containing time-shifted copies of the convolution kernel K
X = sweep(X, 1, STATS=rowSums(X), FUN="/") # make all rowSums sum to 1, also at edges (this fixed edge effect a bit)
nc = ncol(X)
# we add 1st & 2nd order fused ridge regularization by row augmenting our banded convolution covariate matrix
Dg1 = bandSparse(nc, k = c(-1,0,1), diag = list(rep(-1,nc)*sqrt(lambdas1), # 1st order difference ridge penalty, i.e. "fused ridge"
                                                rep(1,nc)*sqrt(lambdas1),
                                                rep(0,nc)*sqrt(lambdas1)), symm=FALSE)
Dg2 = bandSparse(nc, k = c(-1,0,1), diag = list(rep(1,nc)*sqrt(lambdas2), # 2nd order difference ridge penalty
                                                rep(-2,nc)*sqrt(lambdas2),
                                                rep(1,nc)*sqrt(lambdas2)), symm=FALSE)
X_augm_fusedridge12 = as.matrix(rbind(X, Dg1, Dg2)) # augmented covariate matrix to use for deconvolution
y_augm12 = c(y_obs, rep(0,nc*2)) # we also augment our observed counts with zeros
w12 = 1/c(y_obs+0.1, rep(1,2*nc)) 
# 1/variance Poisson observation weights, equivalent to doing a 1 iteration IRLS to approximate an identity link Poisson GLM 
# (cf http://bwlewis.github.io/GLM/, https://stats.stackexchange.com/questions/405629/are-glms-just-glorified-wls-regressions), 
# here with nonnegativity constraints by using nnls instead of regular least squares
# one could do several iterations using the expected value of y to get a more accurate results if desired
library(nnls)
nonnegfusedridgefit12 = nnls(X_augm_fusedridge12*sqrt(w12), y_augm12*sqrt(w12))
# PS maybe worth checking out this nnls implementation using Forward Active Set Truncation/FAST NNLS
# and using RcppArmadillo, so that it is also optimized for a sparse covariate matrix, 
# https://stackoverflow.com/questions/58006606/rcpp-implementation-of-fast-nonnegative-least-squares
nonnegfusedridgefit12_coef = head(nonnegfusedridgefit12$x, -length(K)) # deconvoluted signal, ie case counts deconvoluted to infections


# plot deconvolution result
plot(t, inf_real, pch=16, type="l", ylab="Count", 
     main="Real infections (green), observed cases (red), infections estimated using\ndeconvolution with nonnegative weighted fused ridge regression\nwith 1st & 2nd order difference ridge penalty (blue)",
     cex.main=0.8, lwd=3, col="green3", ylim=c(0,550)) # theoretical infections without any noise
lines(t, cases_real, col="red", lwd=2) # theoretical case counts without any noise
points(t, cases_obs, col="red", pch=16, cex=0.5) # observed case counts
lines(t, nonnegfusedridgefit12_coef, col="blue", lwd=2) # case counts deconvoluted to infections (blue & green should coincide)
# lines(t, y_expected, col="darkgrey", lwd=3) # model predictions
y_expected = (X_augm_fusedridge12 %*% as.matrix(nonnegfusedridgefit12_coef, ncol=1))[1:nc]


library(Hmisc)
subs = 100:800
c(-50:50)[which.min(sapply(X=-50:50, FUN = function(shift) sum((as.vector(inf_real)[subs]-Lag(as.vector(nonnegfusedridgefit12_coef)[subs], 
                                                                                              shift=shift))^2, na.rm=T)))]



# attempt at trying to calculate weekday-specific shift using nonnegative matrix factorization ####

plot(t, y_obs)
head(rawcases)
cases = rawcases[rawcases$PROVINCE=="All"&rawcases$REGION=="Belgium"&rawcases$SEX=="All"&rawcases$AGEGROUP=="All",]
plot(cases$CASES, type="l")
y = cases$CASES[1:(floor(length(cases$CASES)/7)*7)]
Y = matrix(y, byrow=T, ncol=7)
# Y_shifted = cbind(circshift(Y[,1],0),
#                   circshift(Y[,2],1),
#                   circshift(Y[,3],2),
#                   circshift(Y[,4],3),
#                   circshift(Y[,5],4),
#                   circshift(Y[,6],5),
#                   circshift(Y[,7],6))
plot(rowSums(Y),type="l")
# plot(rowSums(Y_shifted),type="l")
image(Y^0.1)
# image(Y_shifted^0.1)
library(ALS)
nmf = als(CList=list(matrix(rowSums(Y), ncol=1)), # m x comp matrix (52 x 1) or ncol=7 ??
    PsiList=list(Y), # dataset, matrix m x n (52 x 7)
    S=matrix(rep(c(1,0,0,0,0,0,0),1), ncol=1, byrow=T), # n x comp (7 x 1)
    # WList=list(1/(Y+0.1)), # approximate Poisson noise 1/variance weights
    nonnegS=TRUE, # nonnegativity constraints on S
    nonnegC=TRUE, # nonnegativity constraints on C
    normS=2) # columns of S are normalized to sum to 1
C = nmf$CList[[1]]
S = nmf$S
scaling = colSums(S)
S = S/scaling
C = C*scaling
Y_fitted = C %*% t(S)
image(Y_fitted^0.1)
image(Y^0.1)
plot(y, type="l")
lines(as.vector(t(Y_fitted)), col="red")

plot(C, type="l")
lines(rowSums(Y), col="red")
plot(cases$CASES, col="red", type="l")
image(S)
