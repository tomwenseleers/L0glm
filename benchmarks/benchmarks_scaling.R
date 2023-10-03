# BENCHMARKS & TIMINGS OF L0GLM FITS IN FUNCTION OF PROBLEM SIZE & SOLVER USED ####

# 1. SINGLE CHANNEL ####

library(ggplot2)
library(dplyr)
# devtools::load_all(".")
library(L0glm)

# Simulation sizes (with n=p)
sizes <- c(100, 500, 1000, 2000, 5000, 10000, 20000)

# Solvers and methods
solvers <- c("osqp", "glmnet", "eigen", "eigen", "eigen", "eigen")
methods <- c(1, 1, 1, 3, 5, 7)
solver_names <- c("osqp", "bigGlm/glmnet", "eigen QR", "eigen LSCG", "eigen CG", "eigen Cholesky LLT")
crit = "bic"

## 1.1 SINGLE CHANNEL - POISSON NOISE ####

family = poisson(identity)

# Initialize an empty data frame to store results
results <- data.frame()

# Loop over sizes, solvers and methods
for (n in sizes) {
  # Simulate data
  sim <- simul(method = 2, n = n, p = n)
  X <- sim$X
  y <- sim$y
  coefs.true <- sim$beta_true
  method <- 1
  for (i in 1:length(methods)) {
    solver = solvers[[i]]
    method = methods[[i]]
    solver_name = solver_names[[i]]
    # Generate regularisation parameter
    lambda_preset <- preset.lambda(IC = "bic", p=ncol(X), n=nrow(X))
    
    # Fit model and measure time
    if (solver_name %in% c("bigGlm/glmnet", "eigen LSCG")) { rowaug = T } else { rowaug = F }
    fit.time <- system.time(
      L0glm <- L0glm.fit(
        X=X, y=y, family = family, lambda = lambda_preset, 
        nonnegative = TRUE, relaxed = FALSE, 
        solver = solver, method = method, 
        control.solver = do.call(control.solver.gen, 
                                 list(solver = solver, method = method, rowaug = rowaug)))
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
        Estimated_over_true_support_size = metrics$estimated_over_true_support_size
      )
    )
  }
}

# Reshape data for plotting
results_melt <- results %>% 
  tidyr::pivot_longer(c(Time, Precision, Recall, Accuracy, Specificity, F1,Estimated_over_true_support_size, R2), names_to = "metric", values_to = "value") %>% 
  mutate(method = factor(method), solver = factor(solver))
results_melt$Solver = factor(results_melt$Solver, 
                             levels=c("bigGlm/glmnet", "osqp", "eigen Cholesky LLT", "eigen CG", "eigen LSCG", "eigen QR"),
                             labels=c("bigGlm/glmnet", "osqp", "eigen Cholesky LLT", "eigen CG", "eigen LSCG", "eigen QR"))

# Plot
ggplot(results_melt, aes(x = Size, y = value, colour=interaction(Solver))) +
  geom_line() +
  facet_wrap(~metric, scales = "free_y") +
  theme_minimal() +
  labs(color = "Solver/Method") +
  theme(legend.position = "right") +
  xlab("Problem size (n=p)") +
  scale_x_continuous(trans="log10") +
  scale_color_manual(values=c("blue", "magenta", "red", "orange", "green", "cyan3")) 
# ggsave(file="..//benchmarks//plots//timings//timings_metrics.png", width=7, height=5, dpi=600)
# ggsave(file="../benchmarks/plots/timings/timings_metrics.pdf", width=7, height=5, dpi=600)
# LSCG solver could probably be fixed via some thresholding or using other settings for that one
# e.g. better coefficient initialisation

ggplot(results_melt[results_melt$metric=="Time",], aes(x = Size, y = value, colour=interaction(Solver))) +
  geom_line() +
  facet_wrap(~metric, scales = "free_y") +
  theme_minimal() +
  labs(color = "Solver/Method") +
  theme(legend.position = "right") +
  scale_color_manual(values=c("blue", "magenta", "red", "orange", "green", "cyan3")) +
  ylab("Time (s)") + xlab("Problem size (n=p)") +
  scale_y_continuous(trans="log10") + scale_x_continuous(trans="log10")
# ggsave(file="..//benchmarks//plots//timings//timings.png", width=7, height=5, dpi=600)
# ggsave(file="../benchmarks/plots/timings/timings.pdf", width=7, height=5, dpi=600)

# save.image("~/GitHub/L0glm/benchmarks/environment_benchmarks timings in function of problem size.RData")

# LSCG 3 or 4 were both bad for n=p=100 - try other preconditioner or
# other initialisation?



## 1.2 SINGLE CHANNEL - GAUSSIAN NOISE ####

family = gaussian(identity)

# Initialize an empty data frame to store results
results <- data.frame()

# Loop over sizes, solvers and methods
for (n in sizes) {
  # Simulate data
  sim <- simul(method = 1, n = n, p = n)
  X <- sim$X
  y <- sim$y
  coefs.true <- sim$beta_true
  method <- 1
  for (i in 1:length(methods)) {
    solver = solvers[[i]]
    method = methods[[i]]
    solver_name = solver_names[[i]]
    # Generate regularisation parameter
    lambda_preset <- preset.lambda(IC = "bic", p=ncol(X), n=nrow(X))
    
    # Fit model and measure time
    if (solver_name %in% c("bigGlm/glmnet", "eigen LSCG")) { rowaug = T } else { rowaug = F }
    fit.time <- system.time(
      L0glm <- L0glm.fit(
        X=X, y=y, family = family, lambda = lambda_preset, 
        nonnegative = TRUE, relaxed = FALSE, 
        solver = solver, method = method, 
        control.solver = do.call(control.solver.gen, 
                                 list(solver = solver, method = method, rowaug = rowaug)))
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
        Estimated_over_true_support_size = metrics$estimated_over_true_support_size
      )
    )
  }
}

# Reshape data for plotting
results_melt2 <- results %>% 
  tidyr::pivot_longer(c(Time, Precision, Recall, Accuracy, Specificity, F1,Estimated_over_true_support_size, R2), names_to = "metric", values_to = "value") %>% 
  mutate(method = factor(method), solver = factor(solver))
results_melt2$Solver = factor(results_melt2$Solver, 
                             levels=c("bigGlm/glmnet", "osqp", "eigen Cholesky LLT", "eigen CG", "eigen LSCG", "eigen QR"),
                             labels=c("bigGlm/glmnet", "osqp", "eigen Cholesky LLT", "eigen CG", "eigen LSCG", "eigen QR"))

# Plot
ggplot(results_melt2, aes(x = Size, y = value, colour=interaction(Solver))) +
  geom_line() +
  facet_wrap(~metric, scales = "free_y") +
  theme_minimal() +
  labs(color = "Solver/Method") +
  theme(legend.position = "right") +
  xlab("Problem size (n=p)") +
  scale_x_continuous(trans="log10") +
  scale_color_manual(values=c("blue", "magenta", "red", "orange", "green", "cyan3")) 
# ggsave(file="..//benchmarks//plots//timings//timings_metrics2.png", width=7, height=5, dpi=600)
# ggsave(file="../benchmarks/plots/timings/timings_metrics2.pdf", width=7, height=5, dpi=600)
# LSCG solver could probably be fixed via some thresholding or using other settings for that one
# e.g. better coefficient initialisation

ggplot(results_melt2[results_melt2$metric=="Time",], aes(x = Size, y = value, 
                                                         colour=interaction(Solver))) +
  geom_line() +
  facet_wrap(~metric, scales = "free_y") +
  theme_minimal() +
  labs(color = "Solver/Method") +
  theme(legend.position = "right") +
  scale_color_manual(values=c("blue", "magenta", "red", "orange", "green", "cyan3")) +
  ylab("Time (s)") + xlab("Problem size (n=p)") +
  scale_y_continuous(trans="log10") + scale_x_continuous(trans="log10")
# ggsave(file="..//benchmarks//plots//timings//timings.png", width=7, height=5, dpi=600)
# ggsave(file="../benchmarks/plots/timings/timings.pdf", width=7, height=5, dpi=600)

# save.image("~/GitHub/L0glm/benchmarks/environment_benchmarks timings in function of problem size.RData")

# LSCG 3 or 4 were both bad for n=p=100 - try other preconditioner or
# other initialisation?


