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


# 1. SIMULATE SPIKE TRAIN BLURRED BY BIGAUSSIAN PSF ####

library(Matrix)
library(L0glm)

set.seed(123)

## ----   latent spike train --------------------------------------
n        <- 10000                      # observation window
npeaks   <- round(n / 20)              # ~5 % non-zero
idx      <- sample(n, npeaks)          # spike locations (1 … n)
heights  <- sample(100:10000, npeaks)  # spike amplitudes

a_true           <- numeric(n)
a_true[idx]      <- heights            # latent, unblurred signal

## ----   point–spread function (PSF) ------------------------------
L          <- 41                       # PSF length
sigma1     <- 1.25                     # left-hand σ
sigma2     <- 5                        # right-hand σ
perc       <- "mean"                   # center at centroid, p=0.5 to center at medioid

g_true   <- bigauspeak(1:L, u=(1+L)/2, w1=sigma1[[1]], w2=sigma2[[1]], h=1) # actual true bigaussian peak used
centroid(g_true) # 21.00015
# g_true   <- g_true/sum(g_true)
# g_true <- convkernel("bigaussian",
#                      L       = L,
#                      params  = list(sigma1 = sigma1,
#                                     sigma2 = sigma2),
#                      perc    = perc)   # centre the centroid

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
                          X = X_sim, 
                          IC = "gic",
                          grouped = TRUE)$y_deconv
r2 = summary(lm(a_true ~ 0 + as.vector(a_est_nonblind)))$r.squared 
r2*100
# 99.98066%

plot(x=a_est_nonblind, y=a_true, 
     xlim=c(min(a_est_nonblind,a_true), max(a_est_nonblind,a_true)), 
     ylim=c(min(a_est_nonblind,a_true), max(a_est_nonblind,a_true)), 
     type="p", col="blue", pch = 16,
     xlab = "estimated value",
     ylab = "true value",
     main = bquote(
       atop(
         "True vs. estimated spike heights (peak areas)\nL0glm deconvolution with known PSF",
         R^2 == .(round(r2 * 100, 3)) * "%") ))


# 2B. BLIND DECONVOLUTION, ESTIMATE UNKNOWN PSF & PEAK LOCATIONS USING ALTERNATING L0glm FITS (BLOCK COORDINATE DESCENT) ####

# alternates L0_sparse_deco and L0_deco_estkernel
# see ?L0_blind_deco, ?L0_sparse_deco & ?L0_deco_estkernel
fit   <- L0_blind_deco(y = y, 
                          family = poisson(identity), 
                          L = 41, # PSF kernel length 
                          IC = "gic",
                          grouped = TRUE)
fit$deviance_trace
plot(fit$deviance_trace, type="l")
a_est_blind = fit$a_hat # estimated spikes
g_est_blind = fit$g_hat # estimated PSF

r2 = summary(lm(a_true ~ 0 + as.vector(a_est_blind)))$r.squared 
r2*100
# 99.98212%

plot(x=a_est_blind, y=a_true, 
     xlim=c(min(a_est_blind,a_true), max(a_est_blind,a_true)), 
     ylim=c(min(a_est_blind,a_true), max(a_est_blind,a_true)), 
     type="p", col="blue", pch = 16,
     xlab = "estimated value",
     ylab = "true value",
     main = bquote(
       atop(
         "True vs. estimated spike heights (peak areas)\nL0glm blind deconvolution",
         R^2 == .(round(r2 * 100, 3)) * "%") ))

matplot(y=cbind(g_true/sum(g_true), 
                g_est_blind/sum(g_est_blind)),
        type="l",
     main = "True vs. estimated bigaussian PSF")




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
