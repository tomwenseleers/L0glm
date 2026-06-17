## ------------------------------------------------------------------------
## 1  Simulate a blurred spike-train with a BIGAUSSIAN PSF & Poisson noise
## ------------------------------------------------------------------------
library(Matrix)
library(splines)
set.seed(123)
family <- poisson(identity)

## --- true point-spread function (PSF) -----------------------------
L       <- 51                       # kernel length
sigma1  <- 1.25                     # left-hand sigma
sigma2  <- 5                        # right-hand sigma
g_true  <- convkernel("bigaussian", L,
                      params = list(sigma1 = sigma1,
                                    sigma2 = sigma2),
                      perc   = 0)   # already centred

## --- latent sparse spike-train ------------------------------------
n        <- 10000                   # timepoints
npeaks   <- round(n / 10)           # 1000, ~10 % non-zero
idx      <- sample(n, npeaks)
amps     <- sample(1000:100000, npeaks)
a_true   <- numeric(n);  a_true[idx] <- amps   # true spikes

## --- baseline -----------------------------------------------------
baseline <- "intercept" # or "bspline" or NULL
if (!is.null(baseline)) {if (baseline=="intercept") {
  X_baseline <- Matrix::Matrix(1, n, 1, sparse = TRUE)
  a_baseline <- 5000
} else if (baseline=="bspline") {
  nknots <- max(4L, floor(n/500))
  X_baseline <- splines::bs(seq_len(n), df = nknots, degree = 3)
  X_baseline <- Matrix::Matrix(X_baseline, sparse = TRUE)
  a_baseline <- runif(nknots, min = 2000, max = 6000)
}} else { X_baseline <- NULL }

## --- generate blurred, noisy observations -------------------------
X_sim  <- convmatrix(g_true, n)
if (!is.null(X_baseline)) { X_full <- cbind(X_baseline, X_sim)
a_full <- c(a_baseline, a_true)
y_baseline <- as.vector(X_baseline %*% a_baseline) # noiseless baseline
} else {
  X_full <- X_sim
  a_full <- a_true
  y_baseline <- rep(0, n) }

y_true <- as.vector(X_full %*% a_full)    # noiseless convolution
y      <- rpois(n, y_true)                # with Poisson noise


# quick sanity plot
plot(y, type = "l", ylim = c(0,max(y)), main = "Simulated blurred spike train with baseline in red")
lines(y_baseline, col="red")

## ------------------------------------------------------------------------
## 2  Deconvolution with an initial PSF guess
## ------------------------------------------------------------------------
g0 <- convkernel("gaussian", L, list(sigma = 6))  # initial guess for PSF
fit <- L0glm_sparse_deco(y, 
                         kernel = g0,
                         baseline = baseline,
                         family = poisson(identity), 
                         IC = "gic")
a_est <- as.vector(fit$coef_spikes)               # estimated spike train
a_baseline <- fit$coef_baseline                   # estimated baseline coefs
X_baseline <- fit$X_baseline                      # design matrix for baseline 
if (!is.null(X_baseline)) { offset <- X_baseline %*% a_baseline } else { offset <- NULL }

## ------------------------------------------------------------------------
## 3  Update PSF estimate given initial spike train estimate
## ------------------------------------------------------------------------
fit <- L0glm_deco_estkernel(y,
                            y_deconv = a_est,  # estimated spike train
                            kernel_est = g0,   # for more stable & faster convergence
                            offset = offset,   # offset for baseline
                            L = L,
                            family = poisson(identity),
                            IC = "gic")
g_est <- fit$kernel

plot(g_true, col="red",
     lwd=7, type="l",
     ylab = "Density",
     main = paste("True vs. estimated bigaussian PSF",
                  "after 1 Poisson L0glm deconvolution iteration",
                  sep = "\n"))
lines(g_est, col="blue", lwd=3)
