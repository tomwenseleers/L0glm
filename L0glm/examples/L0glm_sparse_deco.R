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
## 2  Deconvolution with the known PSF
## ------------------------------------------------------------------------

fit <- L0glm_sparse_deco(y, kernel = g_true,
                         family = poisson(identity),
                         baseline = baseline, 
                         IC = "gic")
a_est_baseline <- as.vector(fit$coef_baseline) 
a_est <- as.vector(fit$coef_spikes)
sum(a_est!=0) 
# 992 with NULL baseline, 1030 with intercept baseline, 1043 with bspline baseline

## ---------- summary statistics ------------------------------------
n_true <- sum(a_true != 0)
n_est  <- sum(a_est  != 0)
tp     <- sum(a_true != 0 & a_est != 0)
precision <- tp / n_est
recall    <- tp / n_true
F1        <- 2 * precision * recall / (precision + recall)
R2p <- 100 *
  summary(lm(a_true ~ 0 + a_est))$r.squared
R2p
# 99.80369 with NULL baseline, 96.70126 with intercept baseline, 97.15414 with bspline baseline

## ---------- scatter-plot with true vs. estimated peak areas -------
main_expr <- bquote(
  atop("Poisson L0-glm deconvolution with known PSF",
       atop("true vs. estimated peak areas",
            italic(R)^2 == .(sprintf("%.2f", R2p)) * "%" ~ "|" ~
              "peaks:" ~ .(n_true) ~ "(true)," ~ .(n_est) ~ "(est.)" ~ "|" ~
              "F1 =" ~ .(sprintf("%.2f", 100 * F1)) * "%")))

plot(a_est, a_true,
     xlim = range(a_est, a_true),
     ylim = range(a_est, a_true),
     pch  = 16, col = "blue",
     xlab = "estimated value",
     ylab = "true value",
     main = main_expr)
