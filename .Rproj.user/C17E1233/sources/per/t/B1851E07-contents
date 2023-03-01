# Simulate data
sim <- simulate_spike_train()
X <- sim$X
y <- sim$y

# Case I: fit nonnegative identity link Poisson GLM with no penalty
L0glm1 <- L0glm(formula = "y ~ 0 + .", data = data.frame(X, y = y), family = poisson(identity),
                lambda = 0, tune.meth = "none", nonnegative = TRUE,
                control.iwls = list(maxit = 100),
                control.l0 = list(maxit = 1, delta = 1E-2, gamma = 1.8))
plot_L0glm_benchmark(x = sim$x, y = y, fit = L0glm1, a.true = sim$a,
                     main="Ground truth vs L0glm estimates")

# Case II: fit nonnegative identity link Poisson GLM with ridge penalty
L0glm2 <- L0glm(formula = "y ~ 0 + .", data = data.frame(X, y = y), family = poisson(identity),
                lambda = 1, tune.meth = "none", nonnegative = TRUE,
                control.iwls = list(maxit = 100),
                control.l0 = list(maxit = 1))
plot_L0glm_benchmark(x = sim$x, y = y, fit = L0glm2, a.true = sim$a,
                     main="Ground truth vs L0glm estimates (with ridge penalty)")

# Case III: fit nonnegative identity link Poisson GLM with adaptive ridge penalty
library(nnls)
a0 <- nnls(A = X*sqrt(1/(y+0.1)),
           b = y*sqrt(1/(y+0.1)))$x
L0glm3 <- L0glm(formula = "y ~ 0 + .", data = data.frame(X, y = y), family = poisson(identity),
                start = a0, lambda = 1, tune.meth = "none", nonnegative = TRUE, normalize = FALSE,
                control.iwls = list(maxit = 100),
                control.l0 = list(maxit = 1))
plot_L0glm_benchmark(x = sim$x, y = y, fit = L0glm3, a.true = sim$a,
                     main="Ground truth vs L0glm estimates (with adaptive ridge penalty)")

# Case IV: fit nonnegative identity link Poisson GLM with L0 penalty and a fixed
#          lambda (no lambda selection)
L0glm4 <- L0glm(formula = "y ~ 0 + .", data = data.frame(X, y = y), family = poisson(identity),
                lambda = 1, tune.meth = "none", nonnegative = TRUE, normalize = FALSE)
plot_L0glm_benchmark(x = sim$x, y = y, fit = L0glm4, a.true = sim$a,
                     main="Ground truth vs L0 penalized L0glm estimates")

\donttest{ # Code below is computationally costly
# Case V: fit nonnegative identity link Poisson GLM with L0 penalty and an
#         optimized lambda using IC on full data
L0glm5 <- L0glm(formula = "y ~ 0 + .", data = data.frame(X, y = y), family = poisson(identity),
                lambda = 10^seq(-3,3, length.out = 51), # Use arbitrary sequence
                tune.crit = "bic", tune.meth = "IC", nonnegative = TRUE, normalize = FALSE)
plot_L0glm_benchmark(x = sim$x, y = y, fit = L0glm5, a.true = sim$a,
                     main="Ground truth vs L0 penalized L0glm estimates (BIC tuning)")

# Case VI: fit nonnegative identity link Poisson GLM with L0 penalty and an
#           optimized lambda using training and validation set
L0glm6 <- L0glm(formula = "y ~ 0 + .", data = data.frame(X, y = y), family = poisson(identity),
                lambda = 10^seq(-3,3, length.out = 51), # Use arbitrary sequence
                tune.crit = "bic", tune.meth = "trainval", nonnegative = TRUE, normalize = FALSE)
plot_L0glm_benchmark(x = sim$x, y = y, fit = L0glm6, a.true = sim$a,
                     main="Ground truth vs L0 penalized L0glm estimates (Training/validation set)")

# Case VII: fit nonnegative identity link Poisson GLM with L0 penalty and an
#           optimized lambda using 3-fold CV
L0glm7 <- L0glm(formula = "y ~ 0 + .", data = data.frame(X, y = y), family = poisson(identity),
                lambda = 10^seq(-3,3, length.out = 51), # Use arbitrary sequence
                tune.crit = "mse", tune.meth = "3-fold", nonnegative = TRUE, normalize = FALSE)
plot_L0glm_benchmark(x = sim$x, y = y, fit = L0glm7, a.true = sim$a,
                     main="Ground truth vs L0 penalized L0glm estimates (3-fold CV)")

# Case VIII: fit nonnegative identity link Poisson GLM with L0 penalty and an
#            optimized lambda using LOOCV
L0glm8 <- L0glm(formula = "y ~ 0 + .", data = data.frame(X, y = y), family = poisson(identity),
                lambda = 10^seq(-3,3, length.out = 51), # Use arbitrary sequence
                tune.meth = "loocv", nonnegative = TRUE, normalize = FALSE)
plot_L0glm_benchmark(x = sim$x, y = y, fit = L0glm8, a.true = sim$a,
                     main="Ground truth vs L0 penalized L0glm estimates (LOOCV)")
}

# Case IX: fit nonnegative identity link Poisson GLM with L0 penalty and a
#          prespecify lambda using BIC
L0glm9 <- L0glm(formula = "y ~ 0 + .", data = data.frame(X, y = y), family = poisson(identity),
                lambda = "aic", # set lambda that minimize the expected BIC
                nonnegative = TRUE, normalize = FALSE)
plot_L0glm_benchmark(x = sim$x, y = y, fit = L0glm9, a.true = sim$a,
                     main="Ground truth vs L0 penalized L0glm estimates (prespecified with BIC)")
