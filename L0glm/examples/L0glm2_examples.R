# Simulate data
set.seed(1)
n <- 1000
p <- 1000
s <- 0.1 # sparsity as proportion of p that have nonzero coefficients // p*s == #of peaks
k <- round(p*s) # nr of nonzero covariates
sim <- simulate_spike_train(n=n, p=p, k=k, mean_beta = 1000, sd_logbeta = 1, family="poisson")
sim <- simulate_spike_train()
coefs_true <- sim$beta_true

X <- sim$X
y <- sim$y
lambda <- 1
pen_weights <- rep(1, ncol(X))
lambdas <- lambda*pen_weights
lower <- rep(0, ncol(X))
upper <- rep(Inf, ncol(X))
nonnegative <- TRUE
start <- NULL
verbose <- TRUE
weights <- 1/(y+0.1) # observation weights




# L0glm
system.time(L0glm1 <- L0glm(formula = "y ~ 0 + .", data = data.frame(X, y = y), family = poisson(identity),
                lambda = 1, tune.meth = "none", nonnegative = TRUE, normalize = FALSE,
                control.iwls = list(maxit = 100),
                control.l0 = list(maxit = 100, delta = 1E-2, gamma = 1.8))) #Time = 70s

plot_L0glm_benchmark(x = sim$x, y = y, fit = L0glm1, a.true = sim$beta_true,
                     main="Ground truth vs L0glm estimates")

cor(coefs_true, L0glm1$coefficients)



#L0glm2 OSQP
system.time(L0glm_osqp0 <- L0glm2(formula = "y ~ 0 + .", data = data.frame(X, y = y), family = poisson(identity), 
                 solver = "osqp", lambda = 0, tune.meth = "none", nonnegative = TRUE,
                 control.l0 = list(maxit = 100, delta = 1E-2, gamma = 1.8),
                 control.iwls = list(maxit = 100),
                 verbose = TRUE))
plot_L0glm_benchmark(x = sim$x, y = y, fit = L0glm_osqp0, a.true = sim$beta_true,
                     main="Ground truth vs L0glm2 estimates")

cor(coefs_true, L0glm_osqp0$coefficients)



system.time(L0glm_osqp1 <- L0glm2(formula = "y ~ 0 + .", data = data.frame(X, y = y), family = poisson(identity), 
                 solver = "osqp", lambda = 1, nonnegative = nonnegative, lower = lower,
                 control.l0 = list(maxit = 100, delta = 1E-2, gamma = 1.8),
                 control.iwls = list(maxit = 100),
                 upper = upper, verbose = TRUE))

plot_L0glm_benchmark(x = sim$x, y = y, fit = L0glm_osqp1, a.true = sim$beta_true,
                     main="Ground truth vs L0glm2 estimates")

cor(coefs_true, L0glm_osqp1$coefficients)

cor(L0glm_osqp0$coefficients, L0glm_osqp1$coefficients)



#L0glm2 GLMNET

system.time(L0glm_glmnet0 <- L0glm2(formula = "y ~ 0 + .", data = data.frame(X, y = y), family = poisson(identity), 
                                  solver = "glmnet", lambda = 0, nonnegative = nonnegative, lower = lower,
                                  control.l0 = list(maxit = 100, delta = 1E-2, gamma = 1.8),
                                  control.iwls = list(maxit = 100),
                                  upper = upper, verbose = TRUE))

plot_L0glm_benchmark(x = sim$x, y = y, fit = L0glm_glmnet0, a.true = sim$beta_true,
                     main="Ground truth vs L0glm2 estimates")
cor(coefs_true, L0glm_glmnet0$coefficients)


system.time(L0glm_glmnet1 <- L0glm2(formula = "y ~ 0 + .", data = data.frame(X, y = y), family = poisson(identity), 
                                    solver = "glmnet", lambda = 1, nonnegative = TRUE, lower = lower,
                                    control.l0 = list(maxit = 100, delta = 1E-2, gamma = 1.8),
                                    control.iwls = list(maxit = 100),
                                    upper = upper, verbose = TRUE))

plot_L0glm_benchmark(x = sim$x, y = y, fit = L0glm_glmnet1, a.true = sim$beta_true,
                     main="Ground truth vs L0glm2 estimates")

cor(coefs_true, L0glm_glmnet1$coefficients)
cor(L0glm_glmnet0$coefficients, L0glm_glmnet1$coefficients)


#L0glm2 Eigen - QR
system.time(L0glm_qr0 <- L0glm2(formula = "y ~ 0 + .", data = data.frame(X, y = y), family = poisson(identity), 
                                    solver = "glmnet", method = 1, lambda = 0, nonnegative = nonnegative, lower = lower,
                                    control.l0 = list(maxit = 100, delta = 1E-2, gamma = 1.8),
                                    control.iwls = list(maxit = 100),
                                    upper = upper, verbose = TRUE))

plot_L0glm_benchmark(x = sim$x, y = y, fit = L0glm_qr0, a.true = sim$beta_true,
                     main="Ground truth vs L0glm2 estimates")
cor(coefs_true, L0glm_qr0$coefficients)


system.time(L0glm_glmnet1 <- L0glm2(formula = "y ~ 0 + .", data = data.frame(X, y = y), family = poisson(identity), 
                                    solver = "glmnet", lambda = 1, nonnegative = TRUE, lower = lower,
                                    control.l0 = list(maxit = 100, delta = 1E-2, gamma = 1.8),
                                    control.iwls = list(maxit = 100),
                                    upper = upper, verbose = TRUE))

plot_L0glm_benchmark(x = sim$x, y = y, fit = L0glm_glmnet1, a.true = sim$beta_true,
                     main="Ground truth vs L0glm2 estimates")

cor(coefs_true, L0glm_glmnet1$coefficients)
cor(L0glm_glmnet0$coefficients, L0glm_glmnet1$coefficients)

#L0glm2 Eigen - LSCG


#L0glm2 Eigen - CG

