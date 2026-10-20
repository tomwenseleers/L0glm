## Fast nonnegative identity-link Poisson L0 or grouped L0 penalized fit
## All coefficients, including an optional intercept, are always nonnegative.
## No nonnegative, lower or upper argument is needed or accepted.
n <- p <- 1000
npeaks <- 100

sim.pois.mc <- simul(method = 4, n = n, p = p, s = npeaks / p, seed = 2)
dim(sim.pois.mc$beta_true) # 1000 x 176 channels

## Multichannel fit. OpenMP is used across outcome channels when available.
system.time(fit.mc <- L0nnpois(
  sim.pois.mc$X,
  sim.pois.mc$y,
  lambda = "gic",
  grouped = TRUE,
  relaxed = TRUE,
  parallelize = TRUE,
  ncores = round(parallel::detectCores(logical = FALSE)/2)
)) 
# with n <- p <- 1000 and npeaks <- 100   0.56s with 4 cores on my machine
# with n <- p <- 10000 and npeaks <- 1000 9.41s with 4 cores on my machine

dim(fit.mc$coefficients)
sum(fit.mc$coefficients != 0)

\dontrun{
## Compare with the general R fitting path (process-based parallelism):
system.time(fit.mc <- L0glm.fit(
  X = sim.pois.mc$X,
  y = sim.pois.mc$y,
  family = poisson(identity),
  lambda = "gic",
  nonnegative = TRUE,
  openmp = FALSE,
  grouped = TRUE,
  relaxed = TRUE,
  parallelize = TRUE,
  ncores = round(parallel::detectCores(logical = FALSE)/2)
)) 
# with n <- p <- 1000 and npeaks <- 100   11.84s with 4 cores on my machine
# with n <- p <- 10000 and npeaks <- 1000 38.55s with 4 cores on my machine

dim(fit.mc$coefficients)
sum(fit.mc$coefficients != 0)
}

## Single-channel fit uses the same nonnegative model.
system.time(fit.sc <- L0nnpois(
  sim.pois.mc$X,
  sim.pois.mc$y[, 1L],
  lambda = "gic",
  grouped = FALSE,
  relaxed = FALSE,
  parallelize = FALSE
)) 
# with n <- p <- 1000 and npeaks <- 100   0.03s on my machine
# with n <- p <- 10000 and npeaks <- 1000 0.13s on my machine

sum(fit.sc$coefficients != 0)

## Native validation stopping with an independent repeat measurement.
small <- simul(n = 300L, p = 300L, n_nonzero = 20L, m = 4L, seed = 7L)
set.seed(8)
y.validation <- matrix(rpois(length(small$y_true), small$y_true), 300L, 4L)
validated <- L0nnpois(small$X, small$y, lambda = "gic", ncores = 2L,
  validation = list(X = small$X, y = y.validation, patience = 3L))
validated$validation$deviance
validated$validation$selection$history

## Compare with the general R path for a single channel:
system.time(fit.sc <- L0glm.fit(
  X = sim.pois.mc$X,
  y = sim.pois.mc$y[, 1L],
  family = poisson(identity),
  lambda = "gic",
  nonnegative = TRUE,
  openmp = FALSE,
  grouped = FALSE,
  relaxed = FALSE,
  parallelize = FALSE
)) 
# with n <- p <- 1000 and npeaks <- 100   0.04s on my machine
# with n <- p <- 10000 and npeaks <- 1000 0.25s on my machine

sum(fit.sc$coefficients != 0)
