# BLIND DECONVOLUTION EXAMPLE: SIMULTANEOUSLY ESTIMATE CONVOLUTION KERNEL/IMPULSE RESPONSE FUNCTION (=PSF) & PEAK LOCATIONS FROM BLURRED SPIKE TRAIN USING L0glm ####

# 1. SIMULATE SPIKE TRAIN BLURRED BY BIGAUSSIAN PSF ####

library(Matrix)
library(L0glm)

set.seed(123)

## ----   point–spread function (PSF) ------------------------------
L          <- 41                       # PSF length
sigma1     <- 1.25                     # left-hand σ
sigma2     <- 5                        # right-hand σ

g_true <- convkernel("bigaussian", # true PSF: bigaussian peak
                     L       = L,
                     params  = list(sigma1 = sigma1,
                                    sigma2 = sigma2),
                     perc    = 0)   # no centering, as this PSF is already centered at centroid by default
sum(g_true) # 1


## ----   latent spike train --------------------------------------
n        <- 10000                      # observation window
npeaks   <- round(n / 20)              # ~5 % non-zero
idx      <- sample(n, npeaks)          # spike locations (1 … n)
amplitudes  <- sample(1000:100000, npeaks)  # spike amplitudes = peak areas after convolution

a_true           <- numeric(n)
a_true[idx]      <- amplitudes         # latent, unblurred signal


## ----   generate blurred observation -----------------------------
X_sim  <- convmatrix(g_true, n)          # sparse Toeplitz conv-matrix
y_true <- as.vector(X_sim %*% a_true)    # noiseless convolution
y      <- rpois(length(y_true), y_true)  # Poisson noise

## ----   quick sanity plots ---------------------------------------
plot(g_true, type = "l", main = "True bigaussian PSF")
plot(y,      type = "l", main = "Blurred spike train (Poisson noise)")


# 2A. DECONVOLUTION WITH MEASURED/KNOWN PSF ####

# see ?L0_sparse_deco
a_est_nonblind   <- L0_sparse_deco(y = y, 
                          family = poisson(identity), 
                          # X = X_sim, # would be equivalent to passing g_true as kernel
                          kernel = g_true, # measured/known PSF
                          IC = "gic")$y_deconv
R2p = 100*summary(lm(a_true ~ 0 + as.vector(a_est_nonblind)))$r.squared 
R2p
# 99.8731%
main_expr <- bquote(
  atop("true vs estimated peak areas",
       atop("Poisson L0 glm deconvolution with known PSF",
            R^2 == .(sprintf("%.2f", R2p)) * "%")))
plot(x=a_est_nonblind, y=a_true, 
     xlim=c(min(a_est_nonblind,a_true), max(a_est_nonblind,a_true)), 
     ylim=c(min(a_est_nonblind,a_true), max(a_est_nonblind,a_true)), 
     type="p", col="blue", pch = 16,
     xlab = "estimated value",
     ylab = "true value",
     main = main_expr)



# 2B. BLIND DECONVOLUTION, ESTIMATE UNKNOWN PSF & PEAK LOCATIONS USING ALTERNATING L0glm FITS (BLOCK COORDINATE DESCENT) ####

# alternates L0_sparse_deco and L0_deco_estkernel
# see ?L0_blind_deco, ?L0_sparse_deco & ?L0_deco_estkernel
fit   <- L0_blind_deco(y = y, 
                          family = poisson(identity), 
                          L = L, # PSF kernel length 
                          kernel_est = convkernel("gaussian", L, list(sigma = 6)), # intial kernel estimate, here we use a Gaussian with sigma=6, sigma=3.5 to 30 all converge to same solution, only can get stuck with too low w0, e.g. sigma=3
                          IC = "gic")
fit$deviance_trace
plot(fit$deviance_trace, type="l", ylab="Deviance", xlab="Iteration")
a_est_blind = fit$a_hat # estimated spikes
g_est_blind = fit$g_hat # estimated PSF

R2p = 100*summary(lm(a_true ~ 0 + as.vector(a_est_blind)))$r.squared 
R2p
# 99.90046%
main_expr <- bquote(
  atop("true vs estimated peak areas",
       atop("Poisson L0 glm blind deconvolution",
            R^2 == .(sprintf("%.2f", R2p)) * "%")))
plot(x=a_est_blind, y=a_true, 
     xlim=c(min(a_est_blind,a_true), max(a_est_blind,a_true)), 
     ylim=c(min(a_est_blind,a_true), max(a_est_blind,a_true)), 
     type="p", col="blue", pch = 16,
     xlab = "estimated value",
     ylab = "true value",
     main = main_expr)
plot(x=1:L, y=g_true, col="red", 
     lwd=7, type="l",
     xlab = "Kernel index",
     ylab = "Density",
     main = "True (red) vs. estimated (blue) bigaussian PSF")
lines(x=1:L, y=g_est_blind, col="blue", lwd=3)



# real datasets - to try....
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
