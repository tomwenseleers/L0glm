# BENCHMARKS & TIMINGS OF L0GLM FITS IN FUNCTION OF PROBLEM SIZE & SOLVER USED ####

# 1. SINGLE CHANNEL ####

library(ggplot2)
library(dplyr)
devtools::load_all(".")
# library(L0glm)

# Simulation sizes (with n=p)
sizes <- c(25, 50, 100, 500, 1000, 2000, 5000) # , 10000, 20000)

# Solvers and methods
solvers <- c("osqp", 
             # "nnls", # L0glm fit using nls with regularisation
             # "nnls_noreg", # nnls without regularisation
             "glmnet", 
             "eigen", 
             "eigen", 
             "eigen", 
             "eigen")
methods <- c(1, 
             # 1, 1, 
             1, 
             1, 4, 6, 7)
solver_names <- c("osqp", 
                  # "nnls", 
                  # "nnls_noreg", 
                  "bigGlm/glmnet", 
                  "eigen QR", 
                  "eigen LSCG",  # TODO eigen LSCG & CG failing sometimes - FIX
                  "eigen CG", 
                  "eigen Cholesky LLT")
crit = "gic"


## 1.1 SINGLE CHANNEL - POISSON NOISE ####

family = poisson(identity)

# Initialize an empty data frame to store results
results <- data.frame()

# Loop over sizes, solvers and methods
for (n in sizes) {
  # Simulate data
  sim <- simul(method = 2, n = n, p = n, options=list(sd_noise=1))
  X <- sim$X
  y <- sim$y
  coefs.true <- sim$beta_true
  method <- 1
  for (i in 1:length(methods)) {
    solver = solvers[[i]]
    method = methods[[i]]
    solver_name = solver_names[[i]]
    # Generate regularisation parameter
    lambda_preset <- preset.lambda(IC = crit, p=ncol(X), n=nrow(X))
    
    # Fit model and measure time
    # start = init_coefs(X=X, y=y, mu=y+0.1, family=family,
    #                   lambda=lambda_preset*(1/4)) # TODO: use this as starting coefs in L0glm
    fit.time <- system.time(
      L0glm <- L0glm.fit(
        X=X, y=y,
        start=NULL, # start, 
        family = family, lambda = lambda_preset, 
        nonnegative = TRUE, relaxed = FALSE, 
        solver = solver, method = method, 
        control.solver = do.call(control.solver.gen, 
                                 list(solver = solver, 
                                      method = method)))
    )[3]
    print(solver_name)
    print(n)
    print(fit.time)
    # Calculate metrics
    metrics <- get.metrics(
      beta = L0glm$coefficients, X = X, y = y, 
      weights = 1/(y+0.1), # approx 1/variance weights
      family = gaussian(identity), # treating poisson as gaussian here with adapted weights - TODO fix poisson(identity) so that deviance always comes out at sensible values 
      beta_true = as.vector(coefs.true)
    )
    
    # Save results to the data frame
    results <- rbind(
      results, 
      data.frame(
        Solver = solver_name, 
        Size = n, 
        Time = fit.time,
        Precision = metrics$precision, 
        Recall = metrics$recall, 
        Specificity = metrics$specificity, 
        Accuracy = metrics$accuracy, 
        # BIC = metrics$bic, 
        F1 = metrics$F1_score,
        R2 = metrics$R2_coefs,
        Estimated_over_true_support_size = metrics$estimated_over_true_support_size,
        deviance = metrics$deviance_true/n # mean per obs deviance relative to ground truth
      )
    )
  }
}

# Reshape data for plotting
results_melt <- results %>% 
  tidyr::pivot_longer(c(Time, Precision, Recall, Accuracy, Specificity, 
                        F1, Estimated_over_true_support_size, R2, deviance), 
                      names_to = "metric", values_to = "value") %>% 
  mutate(method = factor(method), solver = factor(solver))
results_melt$Solver = factor(results_melt$Solver, 
                             levels=rev(c("eigen Cholesky LLT", "eigen CG", "eigen LSCG", "eigen QR", "osqp", 
                                          # "nnls", "nnls_noreg", 
                                          "bigGlm/glmnet")),
                             labels=rev(c("L0glm_eigen Cholesky LLT", "L0glm_eigen CG", "L0glm_eigen LSCG", "L0glm_eigen QR", "L0glm_osqp", 
                                          # "L0glm_nnls", "nnls_noreg", 
                                          "L0glm_bigGlm/glmnet")))

# Plot
ggplot(results_melt, aes(x = Size, y = value, colour=interaction(Solver))) +
  geom_line(aes(linewidth=I(1))) +
  facet_wrap(~metric, scales = "free_y") +
  theme_minimal() +
  labs(title = "Benchmarks different L0glm solvers for varying problem sizes\n(deconvolution of blurred spike train with Poisson noise)", color = "method / solver") +
  theme(legend.position = "right") +
  xlab("Problem size (n=p)") +
  scale_x_continuous(trans="log10") +
  scale_color_manual(values=rev(c("green2", "magenta", "red", "orange", "cyan2", 
                                  # "blue", "cyan2", # nnls
                                  "purple"))) 
ggsave(file="..//benchmarks//plots//timings//timings_metrics.png", width=7, height=5, dpi=600)
ggsave(file="../benchmarks/plots/timings/timings_metrics.pdf", width=7, height=5, dpi=600)

ggplot(results_melt[results_melt$metric=="Time",], aes(x = Size, y = value, colour=interaction(Solver))) +
  geom_line(aes(linewidth=I(1))) +
  # facet_wrap(~metric, scales = "free_y") +
  theme_minimal() +
  labs(title = "Timings different L0glm solvers for varying problem sizes\n(deconvolution of blurred spike train with Poisson noise)", color = "method / solver") +
  theme(legend.position = "right") +
  scale_color_manual(values=rev(c("green2", "magenta", "red", "orange", "cyan2", 
                                  # "blue", "cyan2", # nnls
                                  "purple"))) + 
  ylab("Time (s)") + xlab("Problem size (n=p)") +
  scale_y_continuous(trans="log10") + scale_x_continuous(trans="log10")
ggsave(file="..//benchmarks//plots//timings//timings.png", width=7, height=5, dpi=600)
ggsave(file="../benchmarks/plots/timings/timings.pdf", width=7, height=5, dpi=600)

# save.image("~/GitHub/L0glm/benchmarks/environment_benchmarks timings in function of problem size.RData")

# results look good for all solvers & solution much better & faster than plain nnls



## 1.2 SINGLE CHANNEL - GAUSSIAN NOISE ####

family = gaussian(identity)

# Initialize an empty data frame to store results
results <- data.frame()

# Loop over sizes, solvers and methods
for (n in sizes) {
  # Simulate data
  sim <- simul(method = 1, n = n, p = n, options=list(sd_noise=1))
  X <- sim$X
  y <- sim$y
  coefs.true <- sim$beta_true
  method <- 1
  for (i in 1:length(methods)) {
    solver = solvers[[i]]
    method = methods[[i]]
    solver_name = solver_names[[i]]
    # Generate regularisation parameter
    lambda_preset <- preset.lambda(IC = crit, p=ncol(X), n=nrow(X))
    
    # Fit model and measure time
    start = 0.01*y / diag(X) # diagonal initialisation - only works for n=p case though
    # start = 0.001*init_coefs(X=X, y=y, mu=y, family=family, lambda=lambda_preset*(1/4)) # as.vector(coefs.true)
    print(cor(as.vector(start), as.vector(coefs.true)))
    fit.time <- system.time(
      L0glm <- L0glm.fit(
        X=X, y=y, 
        start=start,  # NULL
        warm.start=0.7,
        family = family, lambda = lambda_preset, 
        nonnegative = TRUE, relaxed = FALSE, 
        solver = solver, method = method, 
        dispersion = 1,
        control.solver = do.call(control.solver.gen, 
                                 list(solver = solver, 
                                      method = method)))
    )[3]
    print(solver_name)
    print(n)
    print(fit.time)
    # Calculate metrics
    metrics <- get.metrics(
      beta = L0glm$coefficients, X = X, y = y, 
      family = family, beta_true = as.vector(coefs.true)
    )
    
    # Save results to the data frame
    results <- rbind(
      results, 
      data.frame(
        Solver = solver_name, 
        Size = n, 
        Time = fit.time,
        Precision = metrics$precision, 
        Recall = metrics$recall, 
        Specificity = metrics$specificity, 
        Accuracy = metrics$accuracy, 
        # BIC = metrics$bic, 
        F1 = metrics$F1_score,
        R2 = metrics$R2_coefs,
        Estimated_over_true_support_size = metrics$estimated_over_true_support_size,
        deviance = metrics$deviance_true/n # mean per obs deviance relative to ground truth
      )
    )
  }
}

# Reshape data for plotting
results_melt2 <- results %>% 
  tidyr::pivot_longer(c(Time, Precision, Recall, Accuracy, Specificity, F1,
                        Estimated_over_true_support_size, R2, deviance), names_to = "metric", values_to = "value") %>% 
  mutate(method = factor(method), solver = factor(solver))
results_melt2$Solver = factor(results_melt2$Solver, 
                              levels=rev(c("eigen Cholesky LLT", "eigen CG", "eigen LSCG", "eigen QR", "osqp", "nnls", "nnls_noreg", "bigGlm/glmnet")),
                              labels=rev(c("L0glm_eigen Cholesky LLT", "L0glm_eigen CG", "L0glm_eigen LSCG", "L0glm_eigen QR", "L0glm_osqp", "L0glm_nnls", "nnls_noreg", "L0glm_bigGlm/glmnet")))

# Plot
ggplot(results_melt2, aes(x = Size, y = value, colour=interaction(Solver))) +
  geom_line() +
  facet_wrap(~metric, scales = "free_y") +
  theme_minimal() +
  labs(title = "Benchmarks different L0glm solvers for varying problem sizes\n(deconvolution of blurred spike train with Gaussian noise)", color = "method / solver") +
  theme(legend.position = "right") +
  xlab("Problem size (n=p)") +
  scale_x_continuous(trans="log10") +
  scale_color_manual(values=rev(c("green2", "magenta", "red", "orange", "cyan2", 
                                  # "blue", "cyan2", # nnls
                                  "purple")))
ggsave(file="..//benchmarks//plots//timings//timings_metrics2.png", width=7, height=5, dpi=600)
ggsave(file="../benchmarks/plots/timings/timings_metrics2.pdf", width=7, height=5, dpi=600)

ggplot(results_melt2[results_melt2$metric=="Time",], aes(x = Size, y = value, 
                                                         colour=interaction(Solver))) +
  geom_line() +
  # facet_wrap(~metric, scales = "free_y") +
  theme_minimal() +
  labs(title = "Timings different L0glm solvers for varying problem sizes\n(deconvolution of blurred spike train with Gaussian noise)", color = "method / solver") +
  theme(legend.position = "right") +
  scale_color_manual(values=rev(c("green2", "magenta", "red", "orange", "cyan2", 
                                  # "blue", "cyan2", # nnls
                                  "purple"))) + 
  ylab("Time (s)") + xlab("Problem size (n=p)") +
  scale_y_continuous(trans="log10") + scale_x_continuous(trans="log10")
ggsave(file="..//benchmarks//plots//timings//timings.png", width=7, height=5, dpi=600)
ggsave(file="../benchmarks/plots/timings/timings.pdf", width=7, height=5, dpi=600)

# save.image("~/GitHub/L0glm/benchmarks/environment_benchmarks timings in function of problem size.RData")


# Here Cholesky & osqp solver is clearly the best.
# TODO try to improve solution quality for some of the other solvers??
# maybe init_coefs could use some tuning
# maybe init_coefs needs to be clipped within allowed bounds?


