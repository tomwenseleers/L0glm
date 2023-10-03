# EXTENSION: TEST L0glm FIT WITH EXTRA FUSED RIDGE PENALTIES WITH DATASET THAT HAS BOTH SMOOTH AND SPARSE PARTS (HERE DECONVOLUTION OF BLURRED SIGNAL WITH PSF TO RECOVER ORIGINAL SMOOTH SIGNAL) ####

# TODO properly add option for a fused ridge penalty in L0glm_1_internal.R

# IDEA: add option allowing for any Lq penalty or 1st or 2nd order fused Lq
# penalty

# EXTRA OPTIONS NEEDED:
# q = specifying what Lq penalty is desired, e.g. 0, 1 or 2
# 2 would correspond to plain ridge, 1 to LASSO (but with less biased coefs due to
# approximation via adaptive penalties) & 0 to best subset
# We would probably always have to do an L0 fit first to be able to put the scale of the L1 or
# L2 penalties on the L0 penalty scale (by multiplying the penalties with the ratio
# of the L0toLq norm) so that we can auto-set our regularisation level
# and use lambda="aic" or "bic" (lambda=2 or log(n)) also for L1 & L2 penalties

# q_fused = specifying what 1st or 2nd order discrete difference fused Lq penalty is desired
# e.g. fused L0 or fused LASSO to reconstruct step-like signal (e.g. CGH data to 
# look at copy number variation across a genome), fused ridge to reconstruct smooth
# signals (e.g. to deconvolute a signal to an underlying smooth temporal signal), etc
# lambda1 = 1st order fused Lq penalty, e.g. "aic" or numerical value
# lambda2 = 2nd order fused Lq penalty, e.g. "aic" or numerical value
# order = integer IDs with correct order of coefficients for fused ridge penalties
# e.g. in simplest case just 1:p if they are all in order

# We could think about order being a named list of vectors specifying the order
# of different dimensions in which a fused penalty would have to be applied
# e.g. x & y for 2D signals or x & y & z for 3D signals or x & y & time for
# a temporal sequence of 2D signals or probe location (1D) for arrayCGH data to
# look at copy number variation across a genome
# these different fused Lq penalties would then have to be all added via row
# augmentation to design matrix X. Probably best to in that case always coerce
# X to be sparse.


devtools::load_all(".")
library(scales)

# 1. UTILITY FUNCTIONS ####
# make a square n x n banded convolution matrix containing time-shifted copies of convolution kernel K
make_banded_X = function (n, K) { 
  K=K/sum(K)
  require(Matrix)
  X = as(bandSparse(n,k = seq(-(length(K)-1),0,1), 
                          diag = t(replicate(n, rev(K))), symm=FALSE), "generalMatrix")
  X <- apply(X, 2, function (col) col/sum(col))
  X <- as(Matrix(X, sparse=T), "generalMatrix")
  return(X)
}

# linear convolution of vector y with convolution kernel K
linconv = function (y, K) {   K=K/sum(K)
# out=convolve(y, K, type="open") # convolution using default R stats function
# ntoclip=(length(out)-length(y)) # /2
# out=head(out,-ntoclip )
# return(out)
X = make_banded_X(n=length(y), K=K)
out = X %*% as.matrix(y, ncol=1) # convolution by multiplying with banded convolution matrix
return(out)
}

# circular shift of a vector
circshift = function(x, n = 0) {
  if (n == 0) x else c(tail(x, -n), head(x, n))
}

# function to augment design matrix X to add 1st & 2nd order ridge penalties
# TODO still think about allowing the fused ridge penalties to be 
# weighted - right now it is implicitly assuming those differences are
# Gaussian distributed whereas in practice they might rather be Poisson distributed
# in which case the weights on the fused penalties would more appropriately
# be taken to be 1/poisson variance = 1/(diff+0.1) with 0.1 used to avoid
# division by zero
row_augX_fusedridge = function (X, lambda1, lambda2, 
                                order=c(NA,1:(ncol(X)-1))) { 
  # NOTE: assuming intercept is in 1st column here - generalize code that 
  # actual order may be different & also take into account unpenalized terms
  # I tried with this but dim(D1) & dim(D2) is not yet correct
  D1 = cbind(0, sqrt(lambda1)*diff(Diagonal(x=diff(order)), 1)) # maybe central difference would be better
  D1[is.na(D1)] = 0
  # note that dim(diff(Diagonal(x=diff(order)), 2))=499 500 with n=p=500
  D2 = cbind(0, sqrt(lambda2)* diff(Diagonal(x=diff(order)), 2))
  D2[is.na(D2)] = 0
  # # note that dim(diff(Diagonal(x=diff(order)), 2))=498 500 with n=p=500
  # # # slightly faster but not general & giving different dim of D1 & D2:
  # p = ncol(X) # is this correct if there is an intercept term?
  # D1 = bandSparse(p, k = c(-1,0,1), diag = list(rep(-1,p)*sqrt(lambda1), # 1st order difference ridge penalty, i.e. "fused ridge"
  #                                                rep(1,p)*sqrt(lambda1),
  #                                                rep(0,p)*sqrt(lambda1)), symm=FALSE)
  # D2 = bandSparse(p, k = c(-1,0,1), diag = list(rep(1,p)*sqrt(lambda2), # 2nd order difference ridge penalty
  #                                                rep(-2,p)*sqrt(lambda2),
  #                                                rep(1,p)*sqrt(lambda2)), symm=FALSE)
  X_augm_fusedridge12 = as(Matrix(rbind(X, D1, D2), sparse=T), "generalMatrix")
  return(X_augm_fusedridge12) }
# NOTE: for a temporal sequence 1:p
# these would be
# order = 1:ncol(X)
# diff(Diagonal(diff(order)), 1)


# Function to create some possible convolution kernels
gaussian_psf <- function(x=c(1:31), m = mean(x), sd = 3) {
  K = exp(-0.5 * ((x - m)/sd)^2)
  K = K / sum(K)
  return(K)
}

exponential_psf <- function(x=c(1:71), rate = 1/14 ) {
  K = dexp(x, rate=rate) # exponential convolution kernel to go from infections to cases
  K = K / sum(K)
  return(K)
}


# SIMULATE DATA (SMOOTH SPARSE SIGNAL CONVOLUTED BY CONVOLUTION KERNEL WITH GAUSSIAN OR POISSON NOISE) ####
# Set parameters
n <- 500 # Number of observations
p <- 500+1 # Number of predictors (nr of timepoints+1 for intercept)
intercept <- 0
error_sd <- 1 # Noise standard deviation for the response
dispersion <- error_sd^2
noise <- "poisson"
if (noise=="gaussian") { family <- gaussian(identity) } else { family <- poisson(identity) }

# convolution kernel (PSF)
Kgaus = gaussian_psf()
plot(Kgaus, type="l", main="Gaussian convolution kernel") 
Kexp = exponential_psf ()
plot(Kexp, type="l", main="exponential convolution kernel") 

# Create a covariate matrix (X) with time-shifted Gaussian peaks plus columns of 1 for intercept
X <- make_banded_X(n=n, K=Kexp)
X <- cbind(1, X)
colnames(X) = c("(Intercept)", paste0("v",1:(p-1)))
image(X)
plot(X[,100], type="l")

# Create true coefficients (betat_true), here a zero-clipped sine wave
freq <- 7
ampl <- 1000
shift <- 0.8
beta_true <- c(intercept, ampl*(sin(seq(0, freq * pi, length.out = p-1))-shift))
beta_true[beta_true < 0] <- 0
plot(intercept+beta_true[-1], type="l")

# Generate zero-noise & observed noisy response (y) based on X and beta, with Gaussian or Poisson noise
y_true <- as.vector(X %*% as.matrix(beta_true, ncol=1)) # = intercept+linconv(y=y, K=Kexp)
# = as.vector(y_true <- intercept+linconv(y=beta_true, K=Kexp))
set.seed(1)
if (noise=="gaussian") { y <- y_true + rnorm(n, mean = 0, sd = error_sd)
                         y[y<0] <- 0 } else {
  y <- as.numeric(rpois(n, lambda = as.vector(y_true))) }
plot(y, type="l", 
     main=paste0("Observed convolved signal with ", noise, " noise (black) &\ntrue original signal (red)"),
     ylim = c(0, (intercept+max(beta_true))*1.1) )
# lines(y_true, col="green3")
lines(intercept+beta_true[-1], col="red")

# Combine into a data frame
# data <- list(X=X, y=y, y_true=y_true, beta_true=beta)

# Settings for Lq penalty - here we use L0 penalty ####
# TODO modify also in AR.weights & L0glm_1_fit_internal.R & make as option in L0glm
q_Lqpen = 0 # norm of penalty to use in Lq penalty
# q=2 would be a ridge penalty (or here fused ridge penalty when applied to c(0,diff(beta)))
L0overLqnorm = ifelse(q_Lqpen==0, 1, sum(beta_true!=0)/(sum(abs(beta_true)^q_Lqpen)^(1/q_Lqpen)))
L0overLqnorm # later used to scale lambda with, can be estimated from L0glm fit with q=0
# NOTE: this is only needed if we would want to set the Lq penalty to q=1 or 2
# drops out when L0 penalty is used
# lambda for the Lq penalty (here L0 penalty)
lambda = preset.lambda(IC = "aic", p=p, n=n) # options aic, bic, gic, ebic, mbic 
# NOTE: for the L0 penalty lambda we can just use a preset lambda=2 expected to optimise AIC
# for q other than 0 (1 or 2) I think lambda would still have to be multiplied
# by 4*(1/2^(2-q_Lq)), but adapt this in code in L0glm_1_fit_internal

# regular L0glm fit without fused ridge penalties & L0 penalty ####
system.time(L0glmfit <- L0glm.fit( X=as(Matrix(X, sparse=T), "generalMatrix"), y=y,
  # X = X_augm_fusedridge12, y = y_augm12,
  weights = 1/(y+0.1), # TODO: using family=poisson(identity) gives a bug due to NAs in the weights that are produced - fix this
  family = gaussian(identity), # family, 
  lambda = 4*lambda, # for q other than 0: *L0overLqnorm*4*(1/2^(2-q_Lqpen)), # when using q = 1 in adaptive penalties this would have to be lambda_preset*L0overL1norm to be roughly on right scale, with L0overLqnorm calculated from L0glm fit with q=0? 
  nonnegative = TRUE,
  # prior.mean = c(min(y),as.vector(y)),
  # start = rep(50, ncol(X)), # rep(mean(y), ncol(X)), # c(min(y),as.vector(y)), #  # c(min(y),as.vector(y)), # rep(50, ncol(X)), # NULL, # rep(50, ncol(sim$X)), # 10-100 works best, coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
  solver = "glmnet", # note: solver osqp & glmnet recommended here - eigen often fails
  maxit = 1, # just doing single iteration here & using 4x higher lambda, doing more iters made solution worse
  dispersion = 1, # use dispersion=1 if poisson weights are provided
  gamma = 2,
  # tol = 0.1,
  method = 7L, # 0.27s for n=p=10000 with eigen solver methods 7 & 8 (Cholesky), 0.61s with methods 5 & 6 (CG), 1.27s with methods 3, 4 (LSCG), 219s with methods 1 & 2 (QR), 4.8s for n=p=10000 with osqp solver, 2s for n=p=10000 with glmnet solver
  thresh = 1
)) # 0.03s
beta_L0glm <- L0glmfit$coefficients
sum(beta_L0glm!=0)  # 90
sum(beta_true!=0) # 116
plot(y, type="l", 
     main="Observed blurred noisy signal (black) and\nestimated deconvoluted (blue) & true original signal (red)",
     ylim = c(0, (intercept+max(beta_true))*1.1) )
# lines(y_true, col="green3")
lines(intercept+beta_true[-1], col="red")
lines(beta_L0glm[1]+beta_L0glm[-1], col="blue", pch=16) # ugh very bad solution...
graph2png(file="..//benchmarks/plots/fused_L0glm/0_deconvolution_ridge_penalty.png", 
          width=7, height=5)

do.call(cbind, get.metrics(beta = beta_L0glm,  # some metrics
                           X = X, y = y, 
                           family = family, 
                           beta_true = beta_true))
# with intercept=10
#        p   n m deviance deviance_true    loglik MD_coefs MAD_coefs sMAPE_coefs bias_coefs  R2_coefs
# [1,] 501 500 1 479.9987       190.948 -1096.493 1.006498  22.35788    0.275295  0.6094027 0.6052805
#      MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true
# [1,]    27.21945     81.16529         0.8141     0.6094027   0.7017686               117
#      support_size_estimated estimated_over_true_support_size TP FP  TN FN precision    recall
# [1,]                    108                        0.9230769 92 16 368 25 0.8518519 0.7863248
#      specificity  accuracy  F1_score      gic     ebic      bic      aic    mbic
# [1,]   0.9583333 0.9181637 0.8177778 3419.556 3535.772 2864.163 2408.985 3907.51
# NOTE: with intercept=0 instead of intercept=10 I get a loglik of -Ing & a ridiculously
# high deviance & Inf values for aic & bic - fix this
# with intercept=0
#        p   n m deviance deviance_true loglik  MD_coefs MAD_coefs sMAPE_coefs bias_coefs  R2_coefs
# [1,] 501 500 1    1e+10         8e+10   -Inf -0.662296  19.33554   0.2127885  0.6300133 0.6326411
#      MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true
# [1,]    27.27991     76.66781      0.5623006     0.6300133    0.747707               116
#      support_size_estimated estimated_over_true_support_size TP FP  TN FN precision    recall
# [1,]                     90                        0.7758621 88  2 383 28 0.9777778 0.7586207
#      specificity  accuracy  F1_score gic ebic bic aic mbic
# [1,]   0.9948052 0.9401198 0.8543689 Inf  Inf Inf Inf  Inf

# CRAP SOLUTION - LET'S SEE WHAT WE GET BY ADDING FUSED RIDGE PENALTIES
# we could use this solution to estimate the following though:
q_Lqpen = 2 # if later on we would like to use an L2 ridge penalty instead of an L0 penalty
beta_hat = beta_L0glm
L0overLqnorm_hat = ifelse(q_Lqpen==0, 1, sum(beta_hat!=0)/(sum(abs(beta_hat)^q_Lqpen)^(1/q_Lqpen)))
L0overLqnorm_hat # 0.05148392, later used to scale lambda with if Lq pen with q=2 is used
# close to ground truth:
L0overLqnorm = ifelse(q_Lqpen==0, 1, sum(beta_true!=0)/(sum(abs(beta_true)^q_Lqpen)^(1/q_Lqpen)))
L0overLqnorm # 0.07448627, later used to scale lambda with if Lq pen with q=2 is used
sd_beta_hat = sd(beta_hat)
sd_beta_hat # 83.72941
sd_beta_true = sd(beta_true)
sd_beta_true # 63.07643

#  for the fused ridge penalty scaling factors:
q_Lqfusedpen = 2
diff1_beta_hat = c(0, diff(beta_hat, 1)) # 1st order discrete difference in beta_hat
L0overLqnorm_diff1beta_hat = ifelse(q_Lqfusedpen==0, 1, sum(diff1_beta_hat!=0)/(sum(abs(diff1_beta_hat)^q_Lqfusedpen)^(1/q_Lqfusedpen)))
L0overLqnorm_diff1beta_hat # 0.06246179, use to scale L2 fused ridge penalty to get lambda on L0 like scale to allow one to use lambda=2 to get fit with best AIC?
# that one was quite far off...
sd_diff1_beta_hat = sd(diff1_beta_hat)
sd_diff1_beta_hat # 94.50928 # far too high...
diff2_beta_hat = c(0, 0, diff(beta_hat, 2)) # 1st order discrete difference in beta_hat
L0overLqnorm_diff2beta_hat = ifelse(q_Lqfusedpen==0, 1, sum(diff2_beta_hat!=0)/(sum(abs(diff2_beta_hat)^q_Lqfusedpen)^(1/q_Lqfusedpen)))
L0overLqnorm_diff2beta_hat # 0.07181416, use to scale L2 fused ridge penalty to get lambda on L0 like scale to allow one to use lambda=2 to get fit with best AIC?
# far too low...
sd_diff2_beta_hat = sd(diff2_beta_hat)
sd_diff2_beta_hat # 81.57859
# much too high...

# these are all still quite far from the ground truth 
# (5-10x too low for L0overLqnorm_diff1beta and L0overLqnorm_diff2beta):
diff1_betatrue = c(0, diff(beta_true, 1)) # 1st order discrete difference in beta_true
L0overLqnorm_diff1beta = ifelse(q_Lqfusedpen==0, 1, sum(diff1_betatrue!=0)/(sum(abs(diff1_betatrue)^q_Lqfusedpen)^(1/q_Lqfusedpen)))
L0overLqnorm_diff1beta # 0.7223459, use to scale L2 fused ridge penalty to get lambda on L0 like scale to allow one to use lambda=2 to get fit with best AIC?
sd_diff1_betatrue = sd(diff1_betatrue)
sd_diff1_betatrue # 7.491239
diff2_betatrue = c(0, 0, diff(beta_true, 2)) # 1st order discrete difference in beta_true
L0overLqnorm_diff2beta = ifelse(q_Lqfusedpen==0, 1, sum(diff2_betatrue!=0)/(sum(abs(diff2_betatrue)^q_Lqfusedpen)^(1/q_Lqfusedpen)))
L0overLqnorm_diff2beta # 0.3803447, use to scale L2 fused ridge penalty to get lambda on L0 like scale to allow one to use lambda=2 to get fit with best AIC?
sd_diff2_betatrue = sd(diff2_betatrue)
sd_diff2_betatrue # 14.69763

# ADD FUSED RIDGE PENALTIES TO COVARIATE MATRIX ####
# we add 1st & 2nd order fused ridge regularization by row augmenting 
# our banded convolution covariate matrix

# manually chosen lambdas
lambda1_manual = 0.005/norm(abs(diff(y, differences=1)), "2") # fused ridge penalty on 1st order discrete difference in coefs
lambda2_manual = 70/norm(abs(diff(y, differences=2)), "2") # fused ridge penalty on 2nd order discrete difference in coefs3
lambda1_manual # 2.365528e-05
lambda2_manual # 0.2009479

# my attempt at sensible automatically set fused ridge penalties ####
# here using sd(diff()) as denominator, because sd is still heavily overestimated at
# this stage, later switching to sd^2
lambda1 = 4*(1/2^(2-q_Lqfusedpen))*L0overLqnorm_diff1beta_hat / sd(diff(y))^2
lambda2 = 4*(1/2^(2-q_Lqfusedpen))*L0overLqnorm_diff2beta_hat / sd(diff(y,2))^2
lambda1
lambda2
# NOTE: in practice L0overLqnorm_diff1beta and L0overLqnorm_diff2beta
# and sd(diff(beta))^2 and sd(diff(beta,2))^2
# would still have to be estimated during the IRLS loop or could
# be estimated based on a first L0glm fit without fused ridge penalties
# as we do here

# By scaling the fused ridge penalties lambdas1 and lambdas2 by 
# the variance of the differences, we are effectively adapting the 
# strength of the penalty to the variability in the differences.

# The 4*(1/2^(2-q_Lqfusedpen)) factor reduces to 1 (1/4 in the actual AR.weights)
# for q=0 in which we get the 1/4 factor of Frommlet & Nuel in the actual
# AR.weights code & for q=2 we would get a factor of 4, so that the net
# multiplier would become 1, as the penalty would then become a plain
# nonadaptive ridge penalty, where no multiplier should be added

# TODO: introduce extra options lambda1 and lambda2 and order
# TODO: note that these fused penalties could be weighted by 1/the variance in
# the difference - now I am implicitly assuming the differences are gaussian
# distributed but poisson here would be more reasonable
# which would suggest a weight of the fused penalty of approx (1/(diff+0.1))
# (with the 0.1 added to avoid division by zero)


X_augm_fusedridge12 = row_augX_fusedridge(X, lambda1, lambda2) # augmented sparse covariate matrix to use for deconvolution
# we also augment our observed counts with zeros
y_augm12 = c(as.vector(y), rep(0,nrow(X_augm_fusedridge12)-nrow(X))) 
w12 = 1/c(as.vector(y)+0.1, rep(1,nrow(X_augm_fusedridge12)-nrow(X))) 
# 1/variance Poisson observation weights if Poisson is approximated by weighted Gaussian
# TODO: allow other values than 1 for weights on the fused ridge penalties,
# e.g. 1/(diff+0.1) if they are Poisson distributed


# L0glm fit with fused ridge penalties here manually added ####
system.time(L0glmfit_fused <- L0glm.fit(# X=X, y=y,
                                  X = X_augm_fusedridge12, y = y_augm12,
                                  weights = w12, # TODO: using family=poisson(identity) gives a bug due to NAs in the weights that are produced - fix this
                                          family = gaussian(identity), # family, 
                                          lambda = lambda, # for q other than 0: *L0overLqnorm*4*(1/2^(2-q_Lqpen)), # when using q = 1 in adaptive penalties this would have to be lambda_preset*L0overL1norm to be roughly on right scale, with L0overLqnorm calculated from L0glm fit with q=0? 
                                          nonnegative = TRUE,
                                          prior.mean = c(min(y),as.vector(y))*0.1,
                                          # start = rep(50, ncol(X)), # rep(mean(y), ncol(X)), # c(min(y),as.vector(y)), #  # c(min(y),as.vector(y)), # rep(50, ncol(X)), # NULL, # rep(50, ncol(sim$X)), # 10-100 works best, coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                                          solver = "eigen", # note: solver osqp & glmnet recommended here - eigen often fails
                                          maxit = 20, # 
                                          dispersion = 1, # use dispersion=1 if poisson weights are provided
                                          gamma = 2,
                                          # tol = 0.1,
                                          method = 7L, # 0.27s for n=p=10000 with eigen solver methods 7 & 8 (Cholesky), 0.61s with methods 5 & 6 (CG), 1.27s with methods 3, 4 (LSCG), 219s with methods 1 & 2 (QR), 4.8s for n=p=10000 with osqp solver, 2s for n=p=10000 with glmnet solver
                                          thresh = 1
)) # 0.04s
beta_L0glm_fused <- as.vector(L0glmfit_fused$coefficients)
sum(beta_L0glm_fused!=0)  # 104
sum(beta_true!=0) # 117
plot(y, type="l", 
     main="Observed blurred noisy signal (black) and\nestimated deconvoluted (blue) & true original signal (red)",
     ylim = c(0, (intercept+max(beta_true))*1.1) )
# lines(y_true, col="green3")
lines(intercept+beta_true[-1], col="red")
lines(beta_L0glm_fused[1]+beta_L0glm_fused[-1], col="blue")
graph2png(file="..//benchmarks/plots/fused_L0glm/1_deconvolution_fused_1st_2nd_order_diff_ridge_penalty.png", 
          width=7, height=5)


do.call(cbind, get.metrics(beta = beta_L0glm_fused,  # some metrics
                           X = X, y = y, 
                           family = family, 
                           beta_true = beta_true))
# with intercept=10
#        p   n m deviance deviance_true    loglik   MD_coefs MAD_coefs sMAPE_coefs bias_coefs
# [1,] 501 500 1 581.0464      122.6247 -1188.016 0.09771086  2.688914   0.1533905   1.029325
#      R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true
# [1,] 0.991153   0.3285446     9.041249      0.5157626      1.029325    0.991153               117
#      support_size_estimated estimated_over_true_support_size  TP FP  TN FN precision recall
# [1,]                    149                         1.273504 117 32 352  0 0.7852349      1
#      specificity  accuracy  F1_score      gic     ebic     bic      aic     mbic
# [1,]   0.9166667 0.9361277 0.8796992 4068.246 4228.582 3302.01 2674.033 4741.442

# with intercept=0
#        p   n m deviance deviance_true    loglik  MD_coefs MAD_coefs sMAPE_coefs bias_coefs
# [1,] 501 500 1 359.5339      49.96214 -592.7018 -1.267349  1.906152  0.05243359   1.028523
#       R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true
# [1,] 0.9946135   -5.162127     7.764084       0.213571      1.028523   0.9946135               116
#      support_size_estimated estimated_over_true_support_size  TP FP  TN FN precision recall
# [1,]                    123                         1.060345 116  7 378  0 0.9430894      1
#      specificity  accuracy  F1_score      gic     ebic    bic      aic     mbic
# [1,]   0.9818182 0.9860279 0.9707113 2582.331 2714.689 1949.8 1431.404 3138.057

# re-estimate fused ridge penalty scalings
q_Lqfusedpen = 2
beta_hat = beta_L0glm_fused
diff1_beta_hat = c(0, diff(beta_hat[beta_hat!=0], 1)) # 1st order discrete difference in beta_hat
L0overLqnorm_diff1beta_hat = ifelse(q_Lqfusedpen==0, 1, sum(diff1_beta_hat!=0)/(sum(abs(diff1_beta_hat)^q_Lqfusedpen)^(1/q_Lqfusedpen)))
L0overLqnorm_diff1beta_hat # 1.056036, use to scale L2 fused ridge penalty to get lambda on L0 like scale to allow one to use lambda=2 to get fit with best AIC?
sd_diff1_beta_hat = sd(diff1_beta_hat)
sd_diff1_beta_hat # 6.563973
diff2_beta_hat = c(0, 0, diff(beta_hat[beta_hat!=0], 2)) # 1st order discrete difference in beta_hat
L0overLqnorm_diff2beta_hat = ifelse(q_Lqfusedpen==0, 1, sum(diff2_beta_hat!=0)/(sum(abs(diff2_beta_hat)^q_Lqfusedpen)^(1/q_Lqfusedpen)))
L0overLqnorm_diff2beta_hat # 0.5491539, use to scale L2 fused ridge penalty to get lambda on L0 like scale to allow one to use lambda=2 to get fit with best AIC?
sd_diff2_beta_hat = sd(diff2_beta_hat)
sd_diff2_beta_hat # 13.02989

# these are all quite close to the ground truth now:
diff1_betatrue = c(0, diff(beta_true, 1)) # 1st order discrete difference in beta_true
L0overLqnorm_diff1beta = ifelse(q_Lqfusedpen==0, 1, sum(diff1_betatrue!=0)/(sum(abs(diff1_betatrue)^q_Lqfusedpen)^(1/q_Lqfusedpen)))
L0overLqnorm_diff1beta # 0.7223459, use to scale L2 fused ridge penalty to get lambda on L0 like scale to allow one to use lambda=2 to get fit with best AIC?
sd_diff1_betatrue = sd(diff1_betatrue)
sd_diff1_betatrue # 7.491239
diff2_betatrue = c(0, 0, diff(beta_true, 2)) # 1st order discrete difference in beta_true
L0overLqnorm_diff2beta = ifelse(q_Lqfusedpen==0, 1, sum(diff2_betatrue!=0)/(sum(abs(diff2_betatrue)^q_Lqfusedpen)^(1/q_Lqfusedpen)))
L0overLqnorm_diff2beta # 0.3803447, use to scale L2 fused ridge penalty to get lambda on L0 like scale to allow one to use lambda=2 to get fit with best AIC?
sd_diff2_betatrue = sd(diff2_betatrue)
sd_diff2_betatrue # 14.69763

# update augmented covariate matrix X
lambda1 = L0overLqnorm_diff1beta_hat / sd(diff1_beta_hat)^2 # sd(diff(y))^2 is estimate of sd(diff(beta,1))^2, could be updated during IRLS loop; alternative could be norm(abs(diff(beta, differences=1)), "2")
lambda2 = L0overLqnorm_diff2beta_hat / sd(diff2_beta_hat)^2 # sd(diff(y,2))^2 is estimate of sd(diff(beta,2))^2, could be updated during IRLS loop; alternative could be norm(abs(diff(beta, differences=2)), "2")
# NOTE: now using sd(diff1_beta_hat)^2 & sd(diff2_beta_hat)^2
# NOW had to leave out 4*(1/2^(2-q_Lqfusedpen))* factor - to check
X_augm_fusedridge12 = row_augX_fusedridge(X, lambda1, lambda2) # augmented sparse covariate matrix to use for deconvolution
y_augm12 = c(as.vector(y), rep(0,nrow(X_augm_fusedridge12)-nrow(X))) 


# new adaptive L0glm fit with fused ridge penalty & updated penalties ####
start = beta_L0glm_fused # use first L0glm estimates as starting values
system.time(L0glmfit_fused_ad <- L0glm.fit(# X=X, y=y,
  X = X_augm_fusedridge12, 
  y = y_augm12,
  weights = w12, # TODO: using family=poisson(identity) gives a bug due to NAs in the weights that are produced - fix this
  family = gaussian(identity), # family, 
  lambda = lambda, # # for q other than 0: *L0overLqnorm*4*(1/2^(2-q_Lqpen)), # when using q = 1 in adaptive penalties this would have to be lambda_preset*L0overL1norm to be roughly on right scale, with L0overLqnorm calculated from L0glm fit with q=0? 
  nonnegative = TRUE,
  # prior.mean = start*0.1,
  warm.start = 0.1, 
  start = start, # use starting values of initial L0glm fit
  solver = "eigen", # note: solver osqp & glmnet recommended here - eigen often fails
  maxit = 20, # 
  dispersion = 1, # use dispersion=1 if poisson weights are provided
  gamma = 2,
  # tol = 0.1,
  method = 7L, # 0.27s for n=p=10000 with eigen solver methods 7 & 8 (Cholesky), 0.61s with methods 5 & 6 (CG), 1.27s with methods 3, 4 (LSCG), 219s with methods 1 & 2 (QR), 4.8s for n=p=10000 with osqp solver, 2s for n=p=10000 with glmnet solver
  thresh = 1
)) # 0.05s
beta_L0glm_fused_ad <- L0glmfit_fused_ad$coefficients
sum(beta_L0glm_fused_ad!=0)  # 134
sum(beta_true!=0) # 117
plot(y, type="l", 
     main="Observed blurred noisy signal (black) and\nestimated deconvoluted (blue) & true original signal (red)",
     ylim = c(0, (intercept+max(beta_true))*1.1) )
# lines(y_true, col="green3")
lines(intercept+beta_true[-1], col="red")
lines(beta_L0glm_fused_ad[1]+beta_L0glm_fused_ad[-1], col="blue")
graph2png(file="..//benchmarks/plots/fused_L0glm/2_deconvolution_adaptive_fused_1st_2nd_order_diff_ridge_penalty.png", 
          width=7, height=5)

do.call(cbind, get.metrics(beta = beta_L0glm_fused_ad,  # some metrics
                           X = X, y = y, 
                           family = family, 
                           beta_true = beta_true))
# with intercept=10
# solution quality 
#        p   n m deviance deviance_true    loglik  MD_coefs MAD_coefs sMAPE_coefs bias_coefs
# [1,] 501 500 1 589.5958      126.2437 -1177.291 -0.178355   2.86895  0.09625587   1.041226
#       R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true
# [1,] 0.9908357  -0.6668349     10.72645       0.359882      1.041226   0.9908357               117
#      support_size_estimated estimated_over_true_support_size  TP FP  TN FN precision recall
# [1,]                    134                         1.145299 117 17 367  0 0.8731343      1
#      specificity  accuracy  F1_score      gic     ebic     bic      aic     mbic
# [1,]   0.9557292 0.9660679 0.9322709 3876.438 4020.633 3187.34 2622.582 4481.863

# with intercept=0
#        p   n m deviance deviance_true    loglik  MD_coefs MAD_coefs sMAPE_coefs bias_coefs
# [1,] 501 500 1 369.3653      56.53202 -593.6249 -1.529516  1.879855  0.03434931   1.044328
#       R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true
# [1,] 0.9964297   -6.439393     7.914348      0.1446135      1.044328   0.9964297               116
#      support_size_estimated estimated_over_true_support_size  TP FP  TN FN precision recall
# [1,]                    119                         1.025862 116  3 382  0 0.9747899      1
#      specificity accuracy F1_score      gic     ebic      bic     aic     mbic
# [1,]   0.9922078 0.994012 0.987234 2538.749 2666.802 1926.788 1425.25 3076.402
# TODO if I use 2*lambda deviance_true increases to 1e-10 - fix this...

# The way I did it now is essentially adding a fixed 1st & 2nd order fused ridge penalty
# it would also be possible to allow q_Lqfusedpen to be set to either 2, 1 or 0
# for a fused ridge, fused LASSO & fused L0 penalty & add this penalty to the
# current adaptive ridge penalty approximating an Lq norm penalty.
# This would require adaptation of the adaptive ridge penalty applied in each
# IRLS iteration.


# TODO check how quality of result would compare with that of
# EpiNow2, DeconRNASeq, MetaboDecon1D, decon, deconvolveR, 
# incidental
# or https://epiforecasts.io/EpiNow2/reference/estimate_infections.html
# or https://search.r-project.org/CRAN/refmans/astrochron/html/deconv.html
# or https://rdrr.io/cran/surveillance/src/R/backprojNP.R
# see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9148632/
# https://cran.r-project.org/web/packages/incidental/index.html


# Nuel gets CIs on all this using a boosting approach (working with bootstrapped datasets)


# SOME UNFINISHED TESTS LOOKING AT CONDITION NUMBER
# what would be a good way to automatically choose between using L0, L1 or L2
# penalized regression?
# for sparse solutions L0 is best
# L1 is usually a poor man's substitute for L0
# L2 may be preferred if solution is non-sparse & under high collinearity
# maybe some criterion in function of the condition number of X or
# the condition number of the moment matrix?

library(RSpectra)
# largest & smallest eigenvalue & condition number
system.time(largest_eigenvalue <- eigs(crossprod(X[,c(1,20,50,100)]), k = 1, which = "LM")$values) # 501.0001
system.time(smallest_eigenvalue <- eigs(crossprod(X[,c(1,20,50,100)]), k = 1, which = "SM")$values) # -3.979039e-13
# Treat very small negative eigenvalues as zero
if (smallest_eigenvalue < 0 && abs(smallest_eigenvalue) <= 1e-3) {
  smallest_eigenvalue <- 1e-3
}
condition_number <- largest_eigenvalue / smallest_eigenvalue
condition_number # 501000.1
1-1/condition_number # 1

XTX <- crossprod(X[,c(1,20,50,100)])
p <- ncol(X)
average_diagonal <- sum(abs(diag(XTX))) / p 
average_off_diagonal <- (sum(abs(XTX)) - sum(abs(diag(XTX)))) / (p  * (p  - 1))
ratio <- average_diagonal / average_off_diagonal
ratio # 41603.66
1-(1/ratio) # 0.9999

XTX <- crossprod(X)
p <- ncol(X)
average_diagonal <- sum(abs(diag(XTX))) / p 
average_off_diagonal <- (sum(abs(XTX)) - sum(abs(diag(XTX)))) / (p  * (p  - 1))
ratio <- average_diagonal / average_off_diagonal
ratio # 174.0154
1-(1/ratio) # 0.994


library(irlba)
# Compute the largest singular value
system.time(largest_singvalue <- irlba(X, nv = 1)$d) # 22.38303 - 22.3830^2=501
# Compute the smallest singular value
system.time(smallest_singvalue <- tail(irlba(X, nu=1, smallest=TRUE, maxit=10, tol=0.1)$d,1)) # 0.03576386 # slow though with large maxit
condition_number <- (largest_singvalue / smallest_singvalue)^2
condition_number
# 293121.4


library(RSpectra)
# largest & smallest eigenvalue & condition number
system.time(largest_eigenvalue <- eigs(crossprod(X_augm_fusedridge12), k = 1, which = "LM")$values) # 501.0001
system.time(smallest_eigenvalue <- eigs(crossprod(X_augm_fusedridge12), k = 1, which = "SM")$values) # -3.979039e-13
# Treat very small negative eigenvalues as zero
if (smallest_eigenvalue < 0 && abs(smallest_eigenvalue) <= 1e-3) {
  smallest_eigenvalue <- 1e-3
}

condition_number <- largest_eigenvalue / smallest_eigenvalue
condition_number # 208408.8

library(irlba)
# Compute the largest singular value
system.time(largest_singvalue <- irlba(X_augm_fusedridge12, nv = 1)$d) 
# Compute the smallest singular value
system.time(smallest_singvalue <- tail(irlba(X_augm_fusedridge12, nu=1, smallest=TRUE, maxit=10, tol=0.1)$d,1)) # 0.03576386 # slow though with large maxit
condition_number <- (largest_singvalue / smallest_singvalue)^2
condition_number
# 193748


X_ortho = simulate_orthogonal(n=500, p=100, n_nonzero=10, 
                              beta_sd=1, 
                              error_sd=1)$X 
X_ortho_augm = row_augX(X_ortho, lambdas=rep(1, ncol(X_ortho)))

XTX <- crossprod(X_ortho)
p <- ncol(X_ortho)
average_diagonal <- sum(abs(diag(XTX))) / p 
average_off_diagonal <- (sum(abs(XTX)) - sum(abs(diag(XTX)))) / (p  * (p  - 1))
ratio <- average_diagonal / average_off_diagonal
ratio # 27.59958
1-(1/ratio) # 0.96

XTX <- crossprod(X_ortho_augm)
p <- ncol(X_ortho_augm)
average_diagonal <- sum(abs(diag(XTX))) / p 
average_off_diagonal <- (sum(abs(XTX)) - sum(abs(diag(XTX)))) / (p  * (p  - 1))
ratio <- average_diagonal / average_off_diagonal
ratio # 27.59958
1-(1/ratio) # 0.96


library(RSpectra)
# largest & smallest eigenvalue & condition number
system.time(largest_eigenvalue <- eigs(crossprod(X_ortho_augm), k = 1, which = "LM")$values) 
system.time(smallest_eigenvalue <- eigs(crossprod(X_ortho_augm), k = 1, which = "SM")$values)
# Treat very small negative eigenvalues as zero
if (smallest_eigenvalue < 0 && abs(smallest_eigenvalue) <= 1e-3) {
  smallest_eigenvalue <- 1e-3
}

condition_number <- largest_eigenvalue / smallest_eigenvalue
condition_number # 5.9
1-(1/condition_number) # 0.8304084, could maybe be used to autoset the q of the Lq penalty??


# library(irlba) # implementation seems buggy
# # Compute the largest singular value
# system.time(largest_singvalue <- irlba(X_ortho_augm, nv = 1)$d) 
# # Compute the smallest singular value
# system.time(smallest_singvalue <- tail(irlba(X_ortho_augm, nu=1, smallest=TRUE, maxit=10, tol=0.1)$d,1)) # 0.03576386 # slow though with large maxit
# condition_number <- (largest_singvalue / smallest_singvalue)^2
# condition_number
# # 1538.424

sv <- svd(X_ortho_augm)
condition_number_svd <- (sv$d[1] / sv$d[length(sv$d)])^2
condition_number_svd # 5.89
1-1/condition_number_svd # 0.83


# Orthogonalize X
Q <- qr.Q(qr(X_ortho_augm))
# Compute condition number using SVD
sv <- svd(Q)
condition_number_svd <- (sv$d[1] / sv$d[length(sv$d)])^2
condition_number_svd # 1
1-1/condition_number_svd # 4E-15
# hmm so the simulation function above that produces X <- matrix(rnorm(n*p), n, p)
# produces uncorrelated columns, but not an orthogonal design matrix
# columns should also have unit length


# TODO check out
# https://hal.science/hal-01853459v1/preview/hal_adaptive_spline.pdf
# https://github.com/Monneret/bandsolve/blob/master/R/bandsolve.R
# https://github.com/goepp/aspline