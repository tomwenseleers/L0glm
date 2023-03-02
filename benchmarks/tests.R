library(Matrix)
library(nnls)
library(devtools)
devtools::load_all()
setwd("D:/Documents/GitHub/L0glm/L0glm")
# Simulate data p > n
set.seed(1)
sim <- simulate_spike_train()
sim2 <- simulate_spike_train(n=1000, npeaks=30, peakhrange = c(100, 10000))

X <- sim$X
Xsp <- Matrix(X, sparse=TRUE)
y <- sim$y

coef_true <- sim$a

lambda <- rep(0, ncol(X))
weights <- 1/(y+0.1)
lower <- rep(0, ncol(X))
upper <- rep(Inf, ncol(X))
nonnegative <- TRUE
start <- NULL
verbose <- FALSE

X <- Xsp

system.time(coefs_LSCG <- lmridge(X=Xsp, y, lambda, weights, lower, upper, nonnegative, solver="LSCG"))
lines(coefs_LSCG, col="blue", type="h")

system.time(coefs_osqp <- lmridge(X=Xsp, y, lambda, weights, lower, upper, nonnegative, solver="osqp"))
lines(coefs_osqp, col="green", type="h")

system.time(coefs_nnls <- nnls(X*sqrt(weights), y*sqrt(weights))$x)
lines(coefs_nnls, col="yellow", type="h")

plot(coefs_osqp, coef_true)


# Case I: fit nonnegative identity link Poisson GLM with no penalty
L0glm1 <- L0glm(formula = "y ~ X + .", data = data.frame(X, y = y), solver = "osqp", family = poisson(identity),
                lambda = 0, tune.meth = "none", nonnegative = TRUE, normalize = FALSE,
                control.iwls = list(maxit = 100),
                control.l0 = list(maxit = 1, delta = 1E-2, gamma = 1.8))

plot_L0glm_benchmark(x = sim$x, y = y, fit = L0glm1, a.true = sim$a,
                     main="Ground truth vs L0glm estimates")



# Case I: fit nonnegative identity link Poisson GLM with no penalty
L0glm1 <- L0glm(formula = "y ~ 0 + .", data = data.frame(X, y = y), family = poisson(identity),
                lambda = 0, tune.meth = "none", nonnegative = TRUE, normalize = FALSE,
                control.iwls = list(maxit = 100),
                control.l0 = list(maxit = 1, delta = 1E-2, gamma = 1.8))
plot_L0glm_benchmark(x = sim$x, y = y, fit = L0glm1, a.true = sim$a,
                     main="Ground truth vs L0glm estimates")
