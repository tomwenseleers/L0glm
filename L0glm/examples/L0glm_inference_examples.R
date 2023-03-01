# Simulate data
sim <- simulate_spike_train()
X <- sim$X
y <- sim$y

# Set up the control parameters of the fitting procedure
ctrl.fit <- control.fit.gen()
ctrl.iwls <- control.iwls.gen(maxit = 1)
ctrl.l0 <- control.l0.gen() # No L0 penalty

# Fit nonnegative identity link Poisson GLM with L0 penalty
L0glm_fit <- L0glm(formula = "y ~ 0 + .", data = data.frame(X, y = y),
                   family = poisson(identity), nonnegative = TRUE,
                   normalize = FALSE, lambda = 1, tune.meth = "none",
                   control.fit = ctrl.fit, control.iwls = ctrl.iwls, control.l0 = ctrl.l0)

\donttest{ # Code below is computationally costly
# Perform inference on the coefficients
system.time(L0glm_infer <- L0glm.inference(L0glm_fit, level = 0.95, boot.repl = 1000,
                                           control.fit = ctrl.fit,
                                           control.iwls = ctrl.iwls, control.l0 = ctrl.l0))
# Plot results
plot_L0glm_benchmark(x = sim$x, y = y, fit = L0glm_fit, inference = L0glm_infer, a.true = sim$a,
                     main="Ground truth vs L0glm estimates (with ridge penalty)")
}


