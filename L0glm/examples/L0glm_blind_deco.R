\donttest{
  
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
plot(y, type = "l", ylim = c(0,max(y)),
     main = "Simulated blurred spike train with baseline in red")
lines(y_baseline, col="red")

## ------------------------------------------------------------------------
## 2  Blind deconvolution starting from intial Gaussian PSF guess
## ------------------------------------------------------------------------
g0 <- convkernel("gaussian", L, list(sigma = 6))  # initial guess for PSF
# = Gaussian with sigma=6, sigma=3.5 to 30 all converge to the same solution

fit <- L0glm_blind_deco(y,
                        family = poisson(identity),
                        L = L,                 # PSF kernel length
                        baseline = baseline,
                        kernel_est = g0,
                        kernel_smoothfun = function(y) {
                          log_conc_smooth(y, nknots = round(L / 5))
                        },
                        perc = "mean",         # use centroid centering
                        kernel_true = g_true,  # for plots
                        IC = "gic",
                        plot = TRUE,
                        tol = 1e-3)

plot(fit$deviance_trace, type="l", log="y",
     xlab="Iteration", ylab="Deviance",
     main="Poisson L0 glm blind deconvolution\nconvergence behaviour")
a_est = as.vector(fit$a_hat)   # estimated spikes
g_est = fit$g_hat              # estimated PSF

## ---------- summary statistics ------------------------------------
n_true <- sum(a_true != 0)
n_est  <- sum(a_est  != 0)
tp     <- sum(a_true != 0 & a_est != 0)
precision <- tp / n_est
recall    <- tp / n_true
F1        <- 2 * precision * recall / (precision + recall)
R2p <- 100 * summary(lm(a_true ~ 0 + a_est))$r.squared

## ---------- scatter-plot with true vs. estimated peak areas -------
main_expr <- bquote(
  atop("Poisson L0 glm blind deconvolution",
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
}
