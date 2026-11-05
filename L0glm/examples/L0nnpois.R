## Grouped identity-link Poisson deconvolution with nonnegative coefficients.
n <- p <- 1000L
npeaks <- 100L

sim.pois.mc <- simul(method = 4L, n = n, p = p, n_nonzero = npeaks, seed = 2L)
dim(sim.pois.mc$beta_true)

## Multichannel fit. OpenMP is used across outcome channels when available.
system.time(fit.mc <- L0nnpois(
  sim.pois.mc$X,
  sim.pois.mc$y,
  lambda = "gic",
  grouped = TRUE,
  relaxed = TRUE,
  parallelize = TRUE,
  ncores = 4L
))

dim(fit.mc$coefficients)
sum(fit.mc$coefficients != 0)

## A six-coefficient nonnegative spline background, refined on every channel.
B <- splines::bs(seq(0, 1, length.out=n), df=6, intercept=TRUE)
set.seed(3)
alpha <- outer(c(2, 5, 3, 7, 4, 2), rep(1, ncol(sim.pois.mc$y)))
Y_background <- sim.pois.mc$y + matrix(rpois(length(sim.pois.mc$y), B %*% alpha),
  n, ncol(sim.pois.mc$y))
compulsory <- p + seq_len(ncol(B))
fit_background <- L0nnpois(cbind(sim.pois.mc$X, Matrix::Matrix(B, sparse=TRUE)),
  Y_background, always.include=compulsory, lambda="gic", ncores=4L,
  precluster=TRUE, control.precluster=list(k.max=32L, rowbins=5L))
background_hat <- B %*% as.matrix(coef(fit_background)[compulsory, , drop=FALSE])
fit_background$nuisance$seconds

## Single-channel fit uses the same nonnegative model.
system.time(fit.sc <- L0nnpois(
  sim.pois.mc$X,
  sim.pois.mc$y[, 1L],
  lambda = "gic",
  grouped = FALSE,
  relaxed = FALSE,
  parallelize = FALSE
))

sum(fit.sc$coefficients != 0)

## Native validation stopping with an independent repeat measurement.
small <- simul(n = 300L, p = 300L, n_nonzero = 20L, m = 4L, seed = 7L)
set.seed(8)
y.validation <- matrix(rpois(length(small$y_true), small$y_true), 300L, 4L)
validated <- L0nnpois(small$X, small$y, lambda = "gic", ncores = 2L,
  validation = list(X = small$X, y = y.validation, patience = 3L))
validated$validation$deviance
validated$validation$selection$history
