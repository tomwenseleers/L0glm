# Simulate data
set.seed(1)
sim <- simulate_spike_train()
X <- sim$X
y <- sim$y
lambda <- 1
pen_weights <- rep(1, ncol(X))
lambdas <- lambda*pen_weights
lower <- rep(0, ncol(X))
upper <- rep(Inf, ncol(X))
nonnegative <- TRUE
start <- NULL
verbose <- FALSE
weights <- 1/(y+0.1) # observation weights




# L0glm
L0glm1 <- L0glm(formula = "y ~ 0 + .", data = data.frame(X, y = y), family = poisson(identity),
                lambda = 1, tune.meth = "none", nonnegative = nonnegative, normalize = FALSE,
                control.iwls = list(maxit = 100),
                control.l0 = list(maxit = 100, delta = 1E-2, gamma = 1.8))
plot_L0glm_benchmark(x = sim$x, y = y, fit = L0glm1, a.true = sim$beta_true,
                     main="Ground truth vs L0glm estimates")




#L0glm2 OSQP
L0glm3 <- L0glm2(formula = "y ~ 0 + .", data = data.frame(X, y = y), family = poisson(identity), 
                 solver = "osqp", lambda = 1, tune.meth = "none", nonnegative = nonnegative, lower = lower,
                 control.l0 = list(maxit = 100, delta = 1E-2, gamma = 1.8),
                 control.iwls = list(maxit = 1),
                 upper = upper, verbose = TRUE)

L0glm3 <- L0glm2(formula = "y ~ 0 + .", data = data.frame(X, y = y), family = poisson(identity), 
                 solver = "osqp", lambda = 1, nonnegative = nonnegative, lower = lower,
                 control.l0 = list(maxit = 100, delta = 1E-2, gamma = 1.8),
                 control.iwls = list(maxit = 100),
                 upper = upper, verbose = verbose)


plot_L0glm_benchmark(x = sim$x, y = y, fit = L0glm2, a.true = sim$beta_true,
                     main="Ground truth vs L0glm estimates")

cor(L0glm2$coefficients, L0glm3$coefficients)
