# BLIND DECONVOLUTION EXAMPLE: SIMULTANEOUSLY ESTIMATE CONVOLUTION KERNEL/PSF/IMPULSE RESPONSE FUNCTION & PEAK LOCATIONS FROM BLURRED SPIKE TRAIN USING L0glm ####

# 1. SIMULATE SPIKE TRAIN BLURRED BY BIGAUSSIAN PSF ####

library(Matrix)
library(L0glm)

set.seed(123)

## ----   point–spread function (PSF) ------------------------------
L          <- 51                       # PSF length
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
amplitudes  <- sample(1000:100000, npeaks)  # spike amplitudes = peak areas after convolution with PSF
a_true           <- numeric(n)
a_true[idx]      <- amplitudes         # latent, unblurred signal

## ----   generate blurred noisy observation -----------------------
X_sim  <- convmatrix(g_true, n)          # sparse Toeplitz conv-matrix
y_true <- as.vector(X_sim %*% a_true)    # noiseless convolution
y      <- rpois(length(y_true), y_true)  # with Poisson noise added

## ----   quick sanity plots ---------------------------------------
plot(g_true, type = "l", ylab = "Density", 
     main = "True bigaussian PSF", col="red", lwd=7)
plot(y,      type = "l", main = "Spike train blurred with bigaussian PSF (Poisson noise)")


# 2. DECONVOLUTION WITH MEASURED/KNOWN PSF ####

# see ?L0_sparse_deco
system.time(a_est_nonblind <- L0_sparse_deco(y = y, 
                          family = poisson(identity), 
                          # X = X_sim, # equivalent to passing g_true as kernel
                          kernel = g_true, # measured/known PSF
                          IC = "gic")$y_deconv)
# 0.5s

## ---------- summary statistics ---------------------------------
n_true <- sum(a_true != 0)                      # ground-truth peaks
n_est  <- sum(a_est_nonblind != 0)              # peaks picked by L0glm

tp      <- sum(a_true != 0 & a_est_nonblind != 0)  # true positives
precision <- tp / n_est
recall    <- tp / n_true
F1        <- 2 * precision * recall / (precision + recall)
R2p       <- 100*summary(lm(a_true ~ 0 + as.vector(a_est_nonblind)))$r.squared 

main_expr <- bquote(
  atop("Poisson L0 glm deconvolution with known PSF",
       atop("true vs. estimated peak areas",
            italic(R)^2 == .(sprintf("%.2f", R2p)) * "%" ~ " | " ~
              "peaks:" ~ .(n_true) ~ "(true)," ~ .(n_est) ~ "(est.) | " ~
              "F1 =" ~ .(sprintf("%.2f", 100 * F1)) * "%"))
)
plot(x=a_est_nonblind, y=a_true, 
     xlim=c(min(a_est_nonblind,a_true), max(a_est_nonblind,a_true)), 
     ylim=c(min(a_est_nonblind,a_true), max(a_est_nonblind,a_true)), 
     type="p", col="blue", pch = 16,
     xlab = "estimated value",
     ylab = "true value",
     main = main_expr)


# 3. BLIND DECONVOLUTION, ESTIMATE UNKNOWN PSF & PEAK LOCATIONS USING ALTERNATING L0glm FITS (BLOCK COORDINATE DESCENT) ####
# alternates L0_sparse_deco and L0_deco_estkernel
# see ?L0_blind_deco, ?L0_sparse_deco & ?L0_deco_estkernel

## 3.1 BLIND DECONVOLUTION OF SIMULATED DATA ####

system.time(fit <- L0_blind_deco(y = y, 
                          family = poisson(identity), 
                          L = L, # PSF kernel length 
                          kernel_est = convkernel("gaussian", L, list(sigma = 6)), # intial kernel estimate, here we use a Gaussian with sigma=6, sigma=3.5 to 30 all converge to same solution, only can get stuck with too low w0, e.g. sigma=3
                          kernel_smoothfun = function(y) log_conc_smooth(y, nknots=round(L/5)),
                          IC = "gic",
                          tol = 1e-9))
# 5.32s
fit$deviance_trace
plot(fit$deviance_trace, type="l", log="y",
     xlab="Iteration", ylab="Deviance", 
     main="Poisson L0 glm blind deconvolution convergence behaviour")
a_est_blind = fit$a_hat # estimated spikes
g_est_blind = fit$g_hat # estimated PSF

## ---------- summary statistics ---------------------------------
n_true <- sum(a_true != 0)                   # ground-truth peaks
n_est  <- sum(a_est_blind != 0)              # peaks picked by L0glm blind deconvolution

tp      <- sum(a_true != 0 & a_est_blind != 0)  # true positives
precision <- tp / n_est
recall    <- tp / n_true
F1        <- 2 * precision * recall / (precision + recall)
R2p       <- 100*summary(lm(a_true ~ 0 + as.vector(a_est_blind)))$r.squared 

main_expr <- bquote(
  atop("Poisson L0 glm blind deconvolution",
       atop("true vs. estimated peak areas",
            italic(R)^2 == .(sprintf("%.2f", R2p)) * "%" ~ " | " ~
              "peaks:" ~ .(n_true) ~ "(true)," ~ .(n_est) ~ "(est.) | " ~
              "F1 =" ~ .(sprintf("%.2f", 100 * F1)) * "%"))
)
plot(x=a_est_blind, y=a_true, 
     xlim=c(min(a_est_blind,a_true), max(a_est_blind,a_true)), 
     ylim=c(min(a_est_blind,a_true), max(a_est_blind,a_true)), 
     type="p", col="blue", pch = 16,
     xlab = "estimated value",
     ylab = "true value",
     main = main_expr)

plot(g_true, col="red", 
     lwd=7, type="l",
     ylab = "Density",
     main = "True (red) vs. estimated (blue) bigaussian PSF\nestimated using Poisson L0 glm blind deconvolution")
lines(g_est_blind, col="blue", lwd=3)


## 3.2 BLIND DECONVOLUTION OF REAL GC/MS DATA ####

## 3.2.1 BLIND DECONVOLUTION GCMS_data1 ####
# see ?GCMS_data1 

y_TIC = rowSums(GCMS_data1) # total ion chromatogram (TIC)
y_std = attr(GCMS_data1, "alkanes_TIC") # alkane ladder standard
y_std_subs = y_std[as.numeric(names(y_std))>=(median(as.numeric(rownames(GCMS_data1)))-2)&
                     as.numeric(names(y_std))<=(median(as.numeric(rownames(GCMS_data1)))+2)]
# alkane ladder standard subsetted to time frame of GCMS_data1 section (2 min around middle of window)
plot(y_std_subs, type="l", ylab="Signal")
plot(x=as.numeric(names(y_std)), y=y_std, type="l", col="red", ylab="Signal", xlab="Elution time (min)") # alkane ladder standard
lines(x=as.numeric(rownames(GCMS_data1)), y=y_TIC, col="blue") # TIC

# PSF estimated from alkane ladder standard section
L <- 51
system.time(fit <- L0_blind_deco(y = y_std_subs, 
                                 family = poisson(identity), 
                                 L = L, # PSF kernel length 
                                 kernel_est = convkernel("gaussian", L, list(sigma = 6)), # intial kernel estimate, here we use a Gaussian with sigma=6
                                 kernel_smoothfun = function(y) log_conc_smooth(y, nknots=round(L/5)),
                                 IC = "gic",
                                 tol = 1e-9))
# 6s
fit$deviance_trace
plot(fit$deviance_trace, type="l", log="y",
     xlab="Iteration", ylab="Deviance", 
     main="Poisson L0 glm blind deconvolution convergence behaviour")
a_est_blind = fit$a_hat # estimated spikes
g_est_blind = fit$g_hat # estimated PSF
plot(g_est_blind, type="o")
sum(a_est_blind!=0) # 235

# PSF estimated from TIC of chromatogram section
L <- 51
system.time(fit <- L0_blind_deco(y = y_TIC, 
                                 family = poisson(identity), 
                                 L = L, # PSF kernel length 
                                 kernel_est = convkernel("gaussian", L, list(sigma = 6)), # intial kernel estimate, here we use a Gaussian with sigma=6
                                 kernel_smoothfun = function(y) log_conc_smooth(y, nknots=round(L/5)), # 
                                 IC = "gic",
                                 tol = 1e-9))
# 6s
fit$deviance_trace
plot(fit$deviance_trace, type="l", log="y",
     xlab="Iteration", ylab="Deviance", 
     main="Poisson L0 glm blind deconvolution convergence behaviour")
a_est_blind = fit$a_hat # estimated spikes
g_est_blind = fit$g_hat # estimated PSF
plot(g_est_blind, type="o")
sum(a_est_blind!=0) # 30
plot(a_est_blind, type="h", col="blue")




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
