# Simulate data p > n
set.seed(1)
sim <- simulate_spike_train(n=500, npeaks=10, peakhrange = c(100, 10000))
X <- sim$X
Xsp <- Matrix(X, sparse=TRUE)
y <- sim$y
coef_true <- sim$a
lambda <- rep(0, ncol(X))
weights <- 1/(y+0.1)
lower <- rep(0, ncol(X))
upper <- rep(Inf, ncol(X))
nonnegative <- TRUE
solver <- "osqp"
start <- NULL
verbose <- FALSE

X <- Xsp

coefs = lmridge(X=Xsp, y, lambda, weights, lower, upper, nonnegative, solver)

# coefs = nnls(X*sqrt(weights), y*sqrt(weights))$x
lines(coefs, col="blue", type="h")



norm((coef_osqp-coef_true),type="2")
plot(coef_osqp, coef_true)



# Case I: fit nonnegative identity link Poisson GLM with no penalty
L0glm1 <- L0glm(formula = "y ~ X + .", data = data.frame(X, y = y), solver = "osqp", family = poisson(identity),
                lambda = 0, tune.meth = "none", nonnegative = TRUE, normalize = FALSE,
                control.iwls = list(maxit = 100),
                control.l0 = list(maxit = 1, delta = 1E-2, gamma = 1.8))

plot_L0glm_benchmark(x = sim$x, y = y, fit = L0glm1, a.true = sim$a,
                     main="Ground truth vs L0glm estimates")
