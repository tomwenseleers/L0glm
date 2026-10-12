# BLIND DECONVOLUTION EXAMPLE: SIMULTANEOUSLY ESTIMATE CONVOLUTION KERNEL/IMPULSE RESPONSE FUNCTION (=PSF) & PEAK LOCATIONS FROM BLURRED SPIKE TRAIN USING L0glm ####

# SOME HELPER FUNCTIONS ####

gauspeak = function(x, u, w, h=1) h*exp(((x-u)^2)/(-2*(w^2)))

bigauspeak = function(x, u, w1, w2, h=1, center=TRUE) { # asymmetrical bigaussian peak with centroids centered at u
  ## analytic offset that moves the centroid onto u
  delta <- if (center) (w2 - w1) * sqrt(2 / pi) else 0
  u0    <- u - delta            # new "mode" position
  out <- numeric(length(x))
  left <- x <= u0
  out[left]  <- h * exp(- (x[left]  - u0)^2 / (2 * w1^2))
  out[!left] <- h * exp(- (x[!left] - u0)^2 / (2 * w2^2))
  return(out)
}

# centroid of a peak
centroid <- function(g) { 
  stopifnot(is.numeric(g), length(g) > 0)
  n       <- length(g)
  idx     <- seq_len(n)
  return(sum(idx * g) / sum(g))
  }



# SPARSE DECONVOLUTION FUNCTION - DECONVOLVE BLURRED SPIKE TRAIN WITH PSF TO ORIGINAL IMPULSES USING L0glm ####

L0_sparse_deco <- function(y,                            # response to be deconvolved (vector or matrix, i.e. univariate or multivariate)
                           family = poisson(identity),   # noise model, e.g. poisson/gaussian/Gamma, we only support identity link families here
                           kernel = NULL, # kernel: either vector with blur kernel or length of blur kernel L & kernelfunc
                           kernelfunc = function(x) gauspeak(x, u = (1+L)/2, w = 6, h = 1),
                           L = 41,
                           X = NULL, # (sparse) deconvolution kernel - can be used to specify matrix with time-varying PSFs & which would then override kernel
                           IC = "gic", # preset L0glm regularisation to approximately optimise this IC, choices aic, bic, gic, mbic or ebic 
                           grouped = TRUE) {  # use grouped feature selection (best subset of group selection), only used if outcome y is multivariate
  dev_res <- family$dev.resids # deviance residuals function
  if (is.vector(y)) y <- as.matrix(y)
  n <- nrow(y)
  p <- ncol(y)
  if (!is.null(kernelfunc)&!is.null(L)) { kernel <- kernelfunc(x=1:L) } else { L <- length(kernel) } # kernel to use for deconvolution
  
  if (is.null(X)) { 
    X <- convmatrix(kernel = kernel, n = n)  # covariate matrix with time-shifted versions of kernel / psf 
  } else { X <- Matrix(X, sparse = TRUE) }
  y_deconv <- L0glm.fit(y = y,   # nonnegative L0-norm penalised GLM
                     X = X,  
                     family = family, 
                     intercept = FALSE,
                     # start = XXX,
                     warm.start = 1,
                     nonnegative = TRUE, # imposing nonnegativity constraints on coefficients
                     lambda = IC, 
                     grouped = grouped,
                     method = 1L, # using sparse Cholesky LLT solver
                     solver = "eigen")$coefficients
  y_hat <- X %*% y_deconv # fitted values # cbind(1, )
  deviance <- sum(dev_res(y = as.numeric(y), 
                          mu = as.numeric(y_hat+1E-10), # to avoid Inf values with poisson
                          wt = rep(1, prod(dim(y)))))
  # if (is.vector(y_deconv)) { y_deconv <- y_deconv[-1] } else { y_deconv <- y_deconv[-1,,drop=FALSE] } # remove intercept again
  return(list(y_deconv = y_deconv, fitted.values = y_hat, deviance = deviance))
}




# ESTIMATE BLUR KERNEL / PSF / IMPULSE RESPONSE FUNCTION BY DECONVOLVING SIGNAL y WITH ESTIMATED UNBLURRED SIGNAL y_deconv USING L0glm ####

L0_deco_estkernel <- function(y,     # vector or matrix (for multivariate case) with response(s) to be deconvolved
                           y_deconv, # vector or matrix (for multivariate case) with estimated deconvolved unblurred signal(s)
                           kernel_est = NULL, # vector with prior estimate of kernel, if available - is always normalised to unit L1 norm before use
                           L = ifelse(!is.null(kernel_est), length(kernel_est), 41), # desired length of  convolution kernel / PSF to be estimated
                           family = poisson(identity),   # noise model, we only support identity link families here
                           kernel_smoothfun = function (y) log_conc_smooth(y), # NULL or function to smooth estimated kernel with (here we use log-concave spline)
                           perc = "mean", # shift PSF to centroid
                           IC = "gic") { # preset L0glm regularisation to approximately optimise this IC, choices aic, bic, gic, mbic or ebic)  
  if (is.vector(y)) y <- as.matrix(y)
  if (is.vector(y_deconv)) y_deconv <- Matrix(y_deconv, sparse=TRUE)
  n <- nrow(y)
  
  kernel_est_scaling <- sum(kernel_est) # scale kernel_est / PSF_est to unit L1 norm (standard for PSF)
  kernel_est <- kernel_est/kernel_est_scaling

  # rs <- rowSums(y_deconv) # L1 norm normalise rows of y_deconv so that with Poisson noise estimated kernel would also have expected Poisson noise
  # y_deconv_scalings <- ifelse(rs > 0, 1/rs, 1)                  # 1/row-sum   (rows that are all-zero keep factor 1)
  # y_deconv <- Diagonal(x = y_deconv_scalings) %*% y_deconv      # left-multiply by a sparse diagonal
  
  # L1 norm normalise rows of y_deconv so that with Poisson noise kernel_raw would also have expected Poisson noise
  cs <- colSums(y_deconv)
  y_deconv_scalings <- ifelse(cs > 0, 1/cs, 1)
  y_deconv <- y_deconv %*% Diagonal(x = y_deconv_scalings)
  
  kernels <- do.call(rbind, lapply(1:ncol(y), # loop over outcome channels for multivariate case
                                   function (col) {
   X <- convmatrix(kernel = y_deconv[, col],
                   n = n,
                   shape = "same",
                   n_cols = L,
                   origin = (L + 1L) %/% 2L)
   kernel <- L0glm.fit(y = y[, col, drop = FALSE],  # nonnegative L0-norm penalised GLM, but now shrinking towards shrink towards prior kernel estimate 
                       X = X, 
                       family = family, 
                       intercept = FALSE,
                       # no.pen = seq(1, ncol(X), by=2), # not penalising every other scan line could be similar to penalising the difference between consecutive coefficients, i.e. to using a fused ridge penalty, but with pre-set lambda that will approx optimise GIC
                       nonnegative = TRUE,
                       # start = XXX,
                       warm.start = 1,
                       # prior.mean = kernel_est*y_deconv_scalings[col], # shrink towards prior kernel estimate (more stable)
                       grouped = FALSE, 
                       lambda = IC,
                       method = 1L, # using sparse Cholesky LLT solver in adaptive ridge IRLS iterations
                       solver = "eigen")$coefficients
    # optional: smooth estimated psf (NOTE: maybe also consider total variation denoiser tvR::denoise1?)
    if (!is.null(kernel_smoothfun)) { kernel <- kernel_smoothfun(y = kernel) }
    # center PSF at given percentile p or with "mean" put centroid in middle
    kernel <- shift_kernel(kernel, 
                           perc = perc, # put centroid in middle 
                           method = "linear")
    return(kernel/sum(kernel)) # normalise PSF to unit L1 norm
  } ))
# average obtained kernels for different channels - we assume they should be identical, but also return the separate kernels
kernel <- colMeans(kernels) 
# plot(kernel, type="l", lwd=5)    
# matplot(t(kernels), type="l", add = TRUE)    

return(list(kernels = kernels, kernel = kernel))
}


# BLIND DECONVOLUTION FUNCTION TO SIMULTANEOUSLY ESTIMATE THE IMPULSE RESPONSE FUNCTION / PSF AND THE PEAK LOCATIONS ####

L0_blind_deco <- function(y,                            # response to be deconvolved (vector or matrix, i.e. univariate or multivariate)
                          family = poisson(identity),   # noise model, we only support identity link families here
                          kernel_est = NULL, # vector with prior estimate of kernel, if available - is always normalised to unit L1 norm before use
                          L = ifelse(!is.null(kernel_est), length(kernel_est), 41), # desired length of  convolution kernel / PSF to be estimated
                          kernelfunc = function(x) gauspeak(x, u = (1+L)/2, w = 6, h = 1),
                          kernel_smoothfun = function (y) log_conc_smooth(y), # NULL or function to smooth estimated kernel with (here we use log-concave spline)
                          p = "mean", # center centroid of estimated kernel in middle of vector?
                          IC = "gic", # preset L0glm regularisation to approximately optimise this IC, choices aic, bic, gic, mbic or ebic 
                          grouped = FALSE,  # use grouped feature selection (best subset of group selection), only used if outcome y is multivariate
                          maxiter = 20) { 
  dev_res <- family$dev.resids # deviance residuals function
  if (is.vector(y)) y <- as.matrix(y)
  n <- nrow(y)
  if (!is.null(kernelfunc)&!is.null(L)) { kernel_est <- kernelfunc(x=1:L) } else { L <- length(kernel_est) } # kernel to use for deconvolution
  
  deviances <- NULL
  a_start <- NULL
  kernel_est <- kernel_est/sum(kernel_est) # normalise kernel to unit L1 norm
  g_est <- g_est_raw <- kernel_est # starting estimate of peak shape function g
  X_a <- convmatrix(kernel = kernel_est,
                    n = n)
  
  for (iter in 1:maxiter) {
    
    plot(g_true, type="l", col="grey",lwd=5, 
         main=paste0("TRUE (grey) & ESTIMATED (red) PSF\niteration ", iter-1))
    lines(g_est, col="red", lwd=1)
    
    # SPARSE DECONVOLUTION WITH CURRENT ESTIMATE OF PSF TO ESTIMATE SPIKE TRAIN USING L0GLM & L0_sparse_deco ####
    a_fit <- L0_sparse_deco(y = y, 
                            family = family, 
                            X = X_a, 
                            IC = IC,
                            grouped = TRUE)
    a_est <- a_fit$y_deconv # 0.28s for n=10000
    # plot(a_est, type="h", col="blue")
    
    # ESTIMATE UPDATED PSF BASED ON ESTIMATED SPIKE TRAIN USING L0GLM ####
    g_fit <- L0_deco_estkernel(y = y, 
                               y_deconv = a_est, 
                               kernel_est = g_est, 
                               family = family, 
                               kernel_smoothfun = kernel_smoothfun, # log-concave spline smooth
                               perc = perc,
                               IC = IC)
    g_est <- g_fit$kernel/sum(g_fit$kernel)   # normalise estimated PSF to unit L1 norm
    # lines(g_est, type="l", col="purple")
    
    X_a = convmatrix(kernel = g_est, n = n)  # update covariate matrix with time-shifted versions of estimated psf
    y_hat = as.vector(X_a %*% a_est) # fitted y
    deviance = sum(dev_res(y = y, mu = y_hat+1E-10, wt = rep(1, length(y))))
    deviances = c(deviances, deviance)
    
  }
  
  # print(deviances)
  # print(deviance) # 6642.775 with aic, 7142.061 with bic, 7422.009 with gic, 7421.458 with ebic, 7776.403 with mbic, less false pos in true vs est peak heights with gic, ebic or mbic
  plot(deviances, type="l", ylab="Deviance", xlab="Iteration", 
       main = "Convergence behaviour of\nblind deconvolution algorithm") # algorithm shows nice convergence
  
  # R2 for true vs estimated spike height (=peak area)
  r2 = summary(lm(a_true ~ 0 + a_est))$r.squared 
  # for n=10000 & npeaks=500: 99.91% with aic, 99.96595% with bic, 99.97961% with gic, 99.97999% with ebic, 99.95534% with mbic
  # for n=10000 & npeaks=1000: 99.069% with gic
  # for n=10000 & npeaks=2500: 95.325% with gic
  
  plot(x=a_est, y=a_true, 
       xlim=c(min(a_est,a_true), max(a_est,a_true)), 
       ylim=c(min(a_est,a_true), max(a_est,a_true)), 
       type="p", col="blue", pch = 16,
       xlab = "Estimated value",
       ylab = "True value",
       main = bquote(
         atop(
           "True vs. estimated spike heights (peak areas)",
           R^2 == .(round(r2 * 100, 3)) * "%") ))
  
  table(a_true!=0, a_est!=0)

  return(list())    
}


# blind_deco = L0_blind_deco(y = y,
#               family = poisson(identity),   # noise model, we only support identity link families here
#               kernel_est = NULL, # vector with prior estimate of kernel, if available - is always normalised to unit L1 norm before use
#               L = ifelse(!is.null(kernel_est), length(kernel_est), 41), # desired length of  convolution kernel / PSF to be estimated
#               kernelfunc = function(x) gauspeak(x, u = (1+L)/2, w = 6, h = 1), # initial width of gaussian psf = 6, 3.3 to 30 all converge to same solution, only can get stuck with too low w0, e.g. < 3.3
#               kernel_smoothfun = function (y) log_conc_smooth(y), # NULL or function to smooth estimated kernel with (here we use log-concave spline)
#               center_kernel = TRUE, # center centroid of estimated kernel in middle of vector?
#               IC = "gic", # preset L0glm regularisation to approximately optimise this IC, choices aic, bic, gic, mbic or ebic 
#               grouped = FALSE,
#               maxiter = 20)
#   
  



# 1. SIMULATE SPIKE TRAIN BLURRED BY BIGAUSSIAN PSF ####

library(Matrix)
library(L0glm)

perc = "mean" # center at centroid, p=0.5 to center at medioid

# SIMULATE
set.seed(123)
n = 10000 
x = 1:n
npeaks = round(n/20) # 500
set.seed(123)
u = sample(x = n, size = npeaks, replace = FALSE) # simulated peak centroids
w = rep(2.5,npeaks) # simulated peak widths
w1 = w/2 # 1.25
w2 = w*2 # 5
h = sample(100:10000, npeaks, replace=FALSE) # runif(npeaks, min=100, max=1000) # unknown peak heights
a_true = rep(0, length(x))
a_true[u] = h
L = 41 # length of PSF
g_true = bigauspeak(1:L, u=(1+L)/2, w1=w1[[1]], w2=w2[[1]], h=1) # actual true bigaussian peak used
centroid(g_true) # 21.00015
# g_true = gauspeak(1:L, u=(1+L)/2, w=w[[1]], h=1) # bigauspeak(1:L, u=(1+L)/2, w1=w1[[1]], w2=w2[[1]], h=1) # actual true bigaussian peak used
# g_true = shift_kernel(g_true,
#                       perc = perc, # "mean", # put centroid in middle, already the case - so not needed here
#                       method = "linear")
scaling_g = sum(g_true) # normalise PSF to unit L1 norm
g_true = g_true/scaling_g
a_true = a_true*scaling_g
# plot(g_true, type="l")

X_sim = convmatrix(kernel = g_true, n = n)
y_true =  as.vector(X_sim %*% a_true) # true noiseless y, equivalent to gsignal::conv(a_true, g_true, shape="same")
y = rpois(n = length(y_true), lambda = y_true) # TODO: allow Poisson, Gaussian or Gamma noise

# example with asymmetrical bigaussian peaks
plot(x,y,type="l", main="Simulated sum of bigaussian peaks with Poisson noise") 
# lines(x[u], h, type="h", col="blue")




# 2. BLIND DECONVOLUTION, ESTIMATE UNKNOWN PSF & PEAK LOCATIONS USING ALTERNATING L0glm FITS (BLOCK COORDINATE DESCENT) ####

# real data to test
# L = 41
# y = rowSums(GCMS_data1) # real example: Cataglyphis iberica queen cuticular hydrocarbon section
# y = attr(GCMS_data1, "alkanes_TIC") # alkane ladder
# n = length(y)
# g_true = NULL

# L = 101 # round(length(y)/5)
# # y = rowSums(GCMS_data2) # more challenging real beer GC/MS examples with triangular acid peak
# # y = Yagg1[,2] # Yagg1[,1]
# # y = attr(GCMS_data2, "alkanes_TIC") # alkane ladder
# # plot(Yagg2[,4] ^ 0.2, type="l")
# y = Yagg3[,3]
# y = Yagg3[,2]
# y = Yagg3
# matplot(y ^ 0.2, type="l", lty=1)
# # n = length(y)
# n = nrow(y)
# g_true = NULL

# INPUT ARGS
w0 = 6 # initial width of gaussian psf, 3.3 to 30 all converge to same solution, only can get stuck with too low w0, e.g. < 3.3
kernelfunc = function(x) gauspeak(x, u = (1+L)/2, w = w0, h = 1)
kernel = kernelfunc(x=1:L)
kernel = shift_kernel(kernel,
                      perc = perc,
                      method = "linear")
kernel = kernel/sum(kernel) # normalise PSF to unit L1 norm
# plot(kernel, type="l")
maxiter = 15
kernel_smoothfun <- function (y) log_conc_smooth(y)
IC = "gic"
family = poisson(identity)

deviances = NULL
a_start = NULL
g_est <- g_est_raw <- kernel # starting estimate of peak shape function g
X_a <- convmatrix(kernel = g_est, n = n)  # covariate matrix with time-shifted versions of psf
                             
for (iter in 1:(maxiter+1)) {
  print(iter)

  if (is.null(g_true)) { plot(g_est, type="l", main=paste0("ESTIMATED PSF\niteration ", iter-1)) } else {   
    plot(g_true, type="l", col="grey",lwd=7, 
           main=paste0("TRUE (grey) & ESTIMATED (red) PSF\niteration ", iter-1))
      lines(g_est, col="red", lwd=3) 
      }

  # SPARSE DECONVOLUTION WITH CURRENT ESTIMATE OF PSF TO ESTIMATE SPIKE TRAIN USING L0GLM & L0_sparse_deco ####
  a_fit <- L0_sparse_deco(y = y, 
                          family = family, 
                          X = X_a, 
                          IC = IC,
                          grouped = TRUE)
  a_est <- a_fit$y_deconv # 0.28s for n=10000
  # plot(a_est, type="h", col="blue")

  # ESTIMATE UPDATED PSF BASED ON ESTIMATED SPIKE TRAIN USING L0GLM ####
  g_fit <- L0_deco_estkernel(y = y, 
                             y_deconv = a_est, 
                             kernel_est = g_est, 
                             family = family, 
                             kernel_smoothfun = kernel_smoothfun, # log-concave spline smooth
                             perc = perc, # "mean",
                             IC = IC)
  g_est <- g_fit$kernel/sum(g_fit$kernel)   
# lines(g_est, type="l", col="purple")

X_a = convmatrix(kernel = g_est, n = n) # update covariate matrix with time-shifted versions of estimated psf
y_hat = as.vector(X_a %*% a_est) # fitted y
deviance = sum(family$dev.resids(y = y, mu = y_hat+1E-10, wt = rep(1, length(y))))
deviances = c(deviances, deviance)

}

print(deviances)
print(deviance) # 7418.995 with gic, less false pos in true vs est peak heights with gic, ebic or mbic
plot(deviances, type="l", ylab="Deviance", xlab="Iteration", log="y", 
     main = "Convergence behaviour of\nblind deconvolution algorithm") # algorithm shows nice convergence





# R2 for true vs estimated spike height (=peak area)
r2 = summary(lm(a_true ~ 0 + as.vector(a_est)))$r.squared 
r2*100
# for n=10000 & npeaks=500: 99.98% with gic
# for n=10000 & npeaks=1000: 99.069% with gic
# for n=10000 & npeaks=2500: 95.325% with gic

plot(x=a_est, y=a_true, 
     xlim=c(min(a_est,a_true), max(a_est,a_true)), 
     ylim=c(min(a_est,a_true), max(a_est,a_true)), 
     type="p", col="blue", pch = 16,
     xlab = "Estimated value",
     ylab = "True value",
     main = bquote(
       atop(
         "True vs. estimated spike heights (peak areas)",
         R^2 == .(round(r2 * 100, 3)) * "%") ))

table(a_true!=0, as.vector(a_est)!=0) # support recovery confusion matrix 
#       FALSE TRUE
# FALSE  9496    4
# TRUE      5  495

plot(x=1:length(a_true), y=a_true, type="h", col="red", 
     ylim=c(-max(a_est,a_true), max(a_est,a_true)), 
     xlab = "Time index",
     ylab = "Peak height (mirrored)",
     main="True (red) & estimated (blue) peak heights")
lines(x=1:length(a_est), y=-a_est, type="h", col="blue")

plot(y, type="l", col="red",
     ylim=c(-max(y,y_hat), max(y,y,y_hat)),
     xlab = "Time index",
     ylab = "Intensity",
     main = "Observed & fitted signal") # observed
lines(-y_hat, type="l", col="blue") # fitted

plot(x=y_hat, y=y, pch=16, log="xy", ylim=c(10, 20000), xlim=c(10, 20000), col="blue",
     main = "Observed vs. estimated signal", xlab = "Estimated", ylab = "Observed")
