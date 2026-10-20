# Single-channel Poisson deconvolution: shifted copies of a known PSF.
single <- simul(method=2, n=500, p=500, s=.1, seed=2)
fit_single <- L0glm.fit(single$X, single$y, family=poisson(identity),
  lambda="gic", nonnegative=TRUE, relaxed=TRUE, parallelize=FALSE)
L0glm.metrics(fit_single, beta_true=single$beta_true)[
  c("deviance", "gic", "F1_score", "MCC")]

old_par <- par(mfrow=c(1, 1), mgp=c(2, .7, 0))
matplot(single$x, cbind(single$y, fitted(fit_single)), type="l", lty=1,
  col=c("grey60", "orange"), xlab="Time", ylab="Signal",
  main="Observed and fitted Poisson signal")
legend("topright", c("Observed", "Fitted"), col=c("grey60", "orange"), lty=1)
plot(single$x_beta, as.numeric(single$beta_true), type="h", col="navy",
  xlab="Spike location", ylab="Coefficient", main="Deconvolved spike train",
  ylim=range(0, as.numeric(single$beta_true), as.numeric(coef(fit_single))))
lines(single$x_beta, as.numeric(coef(fit_single)), type="h", col="red")
legend("topright", c("True", "Estimated"), col=c("navy", "red"), lty=1)

# Multichannel Poisson deconvolution: select rows shared across channels.
# Keep six abundant mass channels for a short example.
multi <- simul(method=4, n=500, p=500, s=.05, seed=2)
channels <- head(order(colSums(multi$y), decreasing=TRUE), 6)
for (name in c("y", "y_true", "beta_true"))
  multi[[name]] <- multi[[name]][, channels, drop=FALSE]
fit_multi <- L0glm.fit(multi$X, multi$y, family=poisson(identity),
  lambda="gic", nonnegative=TRUE, grouped=TRUE, relaxed=TRUE,
  parallelize=FALSE)
L0glm.metrics(fit_multi, beta_true=multi$beta_true)[
  c("gic", "F1_score", "MCC", "within_row_F1_score", "within_row_MCC")]
# The within-row scores assess channel sparsity among correctly recovered rows.
dim(coef(fit_multi))                  # original predictors by response channels
par(mfrow=c(2, 1))
plotmat(multi$beta_true, xlab="Spike location", ylab="Channel",
  main="True multichannel coefficients")
plotmat(coef(fit_multi), xlab="Spike location", ylab="Channel",
  main="Estimated multichannel coefficients")
par(old_par)

# High-dimensional Gaussian regression with positive and negative effects.
# Use p <- 1000000L for the original benchmark (X alone needs about 4 GB).
n <- 500L
p <- 20000L
set.seed(123)
X_hd <- matrix(rnorm(n*p), n, p)
support <- sample.int(p, 60L)
beta_hd <- numeric(p)
beta_hd[support] <- runif(60L, -.8, .8)
mean_hd <- as.numeric(X_hd[, support, drop=FALSE] %*% beta_hd[support])
set.seed(1)
y_hd <- mean_hd + rnorm(n)

# screen="auto" activates ISIS for this ultra-wide design. Screened fits
# default to the original-p GIC preset and a relaxed refit.
# Screening threads are independent of BLAS threads (nthreads).
system.time(fit_hd <- L0glm.fit(X_hd, y_hd, screen.nthreads=8, nthreads=1))
fit_hd$screening[c("method", "screen_size", "iterations", "best_iteration", "nthreads")]
fit_hd$screening$history
L0glm.metrics(fit_hd, beta_true=beta_hd)[c("gic", "F1_score", "MCC")]
dim(coef(fit_hd))                    # all p rows, including screened-out zeros
# Explicit alternatives: screen="SIS", "ISIS", or "none".

# Basic controls: intercept constraints, dispersion and binomial fits.
x <- seq(-1, 1, length.out=30)
X <- matrix(x, ncol=1)

# Nonnegativity includes the Gaussian intercept.
constrained <- L0glm.fit(X, -2+3*x, family=gaussian(),
  intercept=TRUE, nonnegative=TRUE, lambda=1e-8,
  dispersion=1, parallelize=FALSE)
as.numeric(constrained$coefficients)

# Estimate Gaussian variance from an automatic or user-supplied start.
set.seed(1)
y_gaussian <- 2*x + rnorm(length(x), sd=.5)
automatic <- L0glm.fit(X, y_gaussian, family=gaussian(),
  lambda="gic", dispersion=NULL, maxit=100, parallelize=FALSE)
initialized <- L0glm.fit(X, y_gaussian, family=gaussian(),
  lambda="gic", dispersion=-.25, maxit=100, parallelize=FALSE)
# dispersion=.25 instead keeps the variance fixed at .25.

# A binomial intercept may remain negative with a nonnegative slope.
y <- round(1000*plogis(-2+x))/1000
logistic <- L0glm.fit(X, y, family=binomial(),
  weights=rep(1000, length(y)), nonnegative=TRUE,
  lambda=1e-8, parallelize=FALSE)
as.numeric(logistic$coefficients)
