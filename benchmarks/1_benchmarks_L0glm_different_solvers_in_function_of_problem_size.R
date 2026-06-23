# BENCHMARKS & TIMINGS OF L0GLM FITS IN FUNCTION OF PROBLEM SIZE & SOLVER USED ####

# devtools::load_all(".")
library(L0glm)
library(dplyr)
library(ggplot2)

# 1. SPARSE CASE ####

## 1.1 SINGLE CHANNEL ####

# Simulation sizes (with n=p)
sizes <- c(100, 500, 1000, 10000)

# Solvers and methods
solvers <- c("osqp", 
             "bigGlm", 
             rep("eigen", 7)
             )
methods <- c(1, 1, 
             3, 4, 5, 6, 7, 8, 9
             )
solver_names <- c("osqp", 
                  "bigGlm", 
                  "eigen LSCG", "eigen PLSCG", 
                  "eigen CG", "eigen PCG diagonal precond", "eigen Cholesky LLT", "eigen Cholesky LDLT",
                  "eigen PCG incomplete Cholesky precond") # QR, QR ColPiv = 1 & 2, BDCSVD & Jacobi SVD = 10 & 11 not benchmarked, since too slow
crit = "gic"


### 1.1.1 SINGLE CHANNEL - POISSON NOISE ####

family = poisson(identity)

# Initialize an empty data frame to store results
results <- data.frame()

# Loop over sizes, solvers and methods
for (n in sizes) {
  # Simulate data
  sim <- simul(method = 2, n = n, p = n, sparse = T, options=list(sd_noise=1))
  X <- sim$X
  y <- sim$y
  coefs.true <- sim$beta_true
  method <- 1
  # Generate regularisation parameter
  lambda_preset <- L0glm:::preset.lambda(IC = crit, p=ncol(X), n=nrow(X))
  start = L0glm:::init_coefs(X=X, y=y, mu=y+0.1, family=family,
                     lambda=lambda_preset) # component-wise initialisation
  for (i in 1:length(methods)) {
    solver = solvers[[i]]
    method = methods[[i]]
    solver_name = solver_names[[i]]
    # Fit model and measure time
    fit.time <- system.time(
      L0glm <- L0glm.fit(
        X=X, y=y,
        start=start, # NULL uses component-wise initialisation above, alt option: rep(1, ncol(X)) 
        warm.start = 1,
        family = family, lambda = lambda_preset, 
        nonnegative = TRUE, 
        relaxed = FALSE, 
        solver = solver, method = method, 
        control.solver = do.call(control.solver.gen, 
                                 list(solver = solver, 
                                      method = method)))
    )[3]
    print(solver_name)
    print(n)
    print(fit.time)
    # Calculate metrics
    metrics <- L0glm:::get.metrics(
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

warnings() # no warnings here

# Reshape data for plotting
results_melt <- results %>% 
  tidyr::pivot_longer(c(Time, Precision, Recall, Accuracy, Specificity, 
                        F1, Estimated_over_true_support_size, R2, deviance), 
                      names_to = "metric", values_to = "value") %>% 
  mutate(method = factor(method), solver = factor(solver))
# results_melt$Solver = factor(results_melt$Solver, 
#                              levels=rev(c("eigen Cholesky LLT", "eigen CG", "eigen LSCG", "eigen QR", "osqp", 
#                                           # "nnls", "nnls_noreg", 
#                                           "bigGlm" 
#                                           )),
#                              labels=rev(c("L0glm_eigen Cholesky LLT", "L0glm_eigen CG", "L0glm_eigen LSCG", "L0glm_eigen QR", "L0glm_osqp", 
#                                           # "L0glm_nnls", "nnls_noreg", 
#                                           "L0glm_bigGlm" 
#                                           )))

# Plot
ggplot(results_melt, aes(x = Size, y = value, colour=interaction(Solver))) +
  geom_line(aes(linewidth=I(1))) +
  facet_wrap(~metric, scales = "free_y") +
  theme_minimal() +
  labs(title = "Benchmarks different L0glm solvers for varying problem sizes\n(deconvolution of blurred spike train with Poisson noise, sparse case)", color = "method / solver") +
  theme(legend.position = "right") +
  xlab("Problem size (n=p)") +
  scale_x_continuous(trans="log10") +
  scale_color_manual(values=rev(c("green4", "magenta", "red", "orange", "cyan2", 
                                  "blue", "yellow2", 
                                  "green", "purple"))) 
ggsave(file="../benchmarks/plots/timings/timings_metrics_sparse_poisson.png", width=7, height=5, dpi=600)
ggsave(file="../benchmarks/plots/timings/timings_metrics_sparse_poisson.pdf", width=7, height=5, dpi=600)

ggplot(results_melt[results_melt$metric=="Time",], aes(x = Size, y = value, colour=interaction(Solver))) +
  geom_line(aes(linewidth=I(1))) +
  # facet_wrap(~metric, scales = "free_y") +
  theme_minimal() +
  labs(title = "Timings different L0glm solvers for varying problem sizes\n(deconvolution of blurred spike train with Poisson noise, sparse case)", color = "method / solver") +
  theme(legend.position = "right") +
  scale_color_manual(values=rev(c("green4", "magenta", "red", "orange", "cyan2", 
                                  "blue", "yellow2", 
                                  "green", "purple"))) +
  ylab("Time (s)") + xlab("Problem size (n=p)") +
  scale_y_continuous(trans="log10") + scale_x_continuous(trans="log10")
ggsave(file="../benchmarks/plots/timings/timings_sparse_poisson.png", width=7, height=5, dpi=600)
ggsave(file="../benchmarks/plots/timings/timings_sparse_poisson.pdf", width=7, height=5, dpi=600)

# save.image("~/GitHub/L0glm/benchmarks/environment_benchmarks timings in function of problem size.RData")

# results look good for all solvers

results_melt[results_melt$metric=="Time"&results_melt$Size==10000,]
# Solver                                 Size metric  value method solver
# 1 osqp                                  10000 Time    4.51  9      eigen 
# 2 bigGlm                                10000 Time    1.61  9      eigen 
# 3 eigen LSCG                            10000 Time   10.2   9      eigen 
# 4 eigen PLSCG                           10000 Time    0.75  9      eigen 
# 5 eigen CG                              10000 Time    4.24  9      eigen 
# 6 eigen PCG diagonal precond            10000 Time    0.460 9      eigen 
# 7 eigen Cholesky LLT                    10000 Time    0.25  9      eigen 
# 8 eigen Cholesky LDLT                   10000 Time    0.25  9      eigen 
# 9 eigen PCG incomplete Cholesky precond 10000 Time    0.280 9      eigen 



### 1.1.2 SINGLE CHANNEL - GAUSSIAN NOISE ####

family = gaussian(identity)

# Initialize an empty data frame to store results
results <- data.frame()

# Loop over sizes, solvers and methods
for (n in sizes) {
  # Simulate data
  sim <- simul(method = 1, n = n, p = n, sparse = T, options=list(sd_noise=1))
  X <- sim$X
  y <- sim$y
  coefs.true <- sim$beta_true
  method <- 1
  # Generate regularisation parameter
  lambda_preset <- L0glm:::preset.lambda(IC = crit, p=ncol(X), n=nrow(X))
  start = L0glm:::init_coefs(X=X, y=y, mu=y, family=family, lambda=lambda_preset) # rescaled component-wise initialisation
  # start = coef(bigGlm::bigGlm(x=X, y=y, family="gaussian", lower=0, path=FALSE))[-1] # this worked less well: initialise with unregularised nonnegative bigGlm fit
  for (i in 1:length(methods)) {
    solver = solvers[[i]]
    method = methods[[i]]
    solver_name = solver_names[[i]]
    
    # Fit model and measure time
    fit.time <- system.time(
      L0glm <- L0glm.fit(
        X=X, y=y, 
        start=NULL, # NULL uses component-wise initialisation above, alt option: rep(1, ncol(X)) 
        warm.start=1,
        family = family, lambda = lambda_preset, 
        nonnegative = TRUE, relaxed = FALSE, 
        solver = solver, method = method, 
        # dispersion = 1,
        control.solver = do.call(control.solver.gen, 
                                 list(solver = solver, 
                                      method = method)))
    )[3]
    print(solver_name)
    print(n)
    print(fit.time)
    # Calculate metrics
    metrics <- L0glm:::get.metrics(
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

warnings()

# Reshape data for plotting
results_melt2 <- results %>% 
  tidyr::pivot_longer(c(Time, Precision, Recall, Accuracy, Specificity, F1,
                        Estimated_over_true_support_size, R2, deviance), names_to = "metric", values_to = "value") %>% 
  mutate(method = factor(method), solver = factor(solver))
# results_melt2$Solver = factor(results_melt2$Solver, 
#                               levels=rev(c("eigen Cholesky LLT", "eigen CG", "eigen LSCG", "eigen QR", "osqp", 
#                                            #"nnls", "nnls_noreg",
#                                            "bigGlm"
#                                            )),
#                               labels=rev(c("L0glm_eigen Cholesky LLT", "L0glm_eigen CG", "L0glm_eigen LSCG", "L0glm_eigen QR", "L0glm_osqp", 
#                                            #"L0glm_nnls", "nnls_noreg", 
#                                            "L0glm_bigGlm"
#                                            )))

# Plot
ggplot(results_melt2, aes(x = Size, y = value, colour=interaction(Solver))) +
  geom_line() +
  facet_wrap(~metric, scales = "free_y") +
  theme_minimal() +
  labs(title = "Benchmarks different L0glm solvers for varying problem sizes\n(deconvolution of blurred spike train with Gaussian noise, sparse case)", color = "method / solver") +
  theme(legend.position = "right") +
  xlab("Problem size (n=p)") +
  scale_x_continuous(trans="log10") +
  scale_color_manual(values=rev(c("green4", "magenta", "red", "orange", "cyan2", 
                                  "blue", "yellow2", 
                                  "green", "purple"))) 
ggsave(file="../benchmarks/plots/timings/timings_metrics_sparse_gaussian.png", width=7, height=5, dpi=600)
ggsave(file="../benchmarks/plots/timings/timings_metrics_sparse_gaussian.pdf", width=7, height=5, dpi=600)

# here bigGlm solver is markedly worse than the rest
# TODO check if this can be improved

ggplot(results_melt2[results_melt2$metric=="Time",], aes(x = Size, y = value, 
                                                         colour=interaction(Solver))) +
  geom_line() +
  # facet_wrap(~metric, scales = "free_y") +
  theme_minimal() +
  labs(title = "Timings different L0glm solvers for varying problem sizes\n(deconvolution of blurred spike train with Gaussian noise, sparse case)", color = "method / solver") +
  theme(legend.position = "right") +
  scale_color_manual(values=rev(c("green4", "magenta", "red", "orange", "cyan2", 
                                  "blue", "yellow2", 
                                  "green", "purple"))) +  
  ylab("Time (s)") + xlab("Problem size (n=p)") +
  scale_y_continuous(trans="log10") + scale_x_continuous(trans="log10")
ggsave(file="../benchmarks/plots/timings/timings_sparse_gaussian.png", width=7, height=5, dpi=600)
ggsave(file="../benchmarks/plots/timings/timings_sparse_gaussian.pdf", width=7, height=5, dpi=600)

# save.image("~/GitHub/L0glm/benchmarks/environment_benchmarks timings in function of problem size.RData")

results_melt2[results_melt2$metric=="Time"&results_melt2$Size==10000,]
# Solver                                 Size metric  value method solver
# 1 osqp                                  10000 Time   51.1   9      eigen 
# 2 bigGlm                                10000 Time    6.65  9      eigen 
# 3 eigen LSCG                            10000 Time    4.52  9      eigen 
# 4 eigen PLSCG                           10000 Time    4.05  9      eigen 
# 5 eigen CG                              10000 Time    2.11  9      eigen 
# 6 eigen PCG diagonal precond            10000 Time    1.68  9      eigen 
# 7 eigen Cholesky LLT                    10000 Time    0.160 9      eigen 
# 8 eigen Cholesky LDLT                   10000 Time    0.150 9      eigen 
# 9 eigen PCG incomplete Cholesky precond 10000 Time    0.230 9      eigen



# 2. DENSE CASE ####

## 2.1 SINGLE CHANNEL ####

# Simulation sizes (with n=p)
sizes <- c(100, 500, 1000, 10000)

# Solvers and methods
solvers <- c("osqp", 
             "bigGlm", 
             rep("eigen", 7)
)
methods <- c(1, 1, 
             3, 4, 5, 6, 7, 8, 9
)
solver_names <- c("osqp", 
                  "bigGlm", 
                  "eigen LSCG", "eigen PLSCG", 
                  "eigen CG", "eigen PCG diagonal precond", "eigen Cholesky LLT", "eigen Cholesky LDLT",
                  "eigen PCG incomplete Cholesky precond") # QR, QR ColPiv = 1 & 2, BDCSVD & Jacobi SVD = 10 & 11 not benchmarked, since too slow
crit = "gic"


### 2.1.1 SINGLE CHANNEL - POISSON NOISE ####

family = poisson(identity)

# Initialize an empty data frame to store results
results <- data.frame()

# Loop over sizes, solvers and methods
for (n in sizes) {
  # Simulate data
  sim <- simul(method = 2, n = n, p = n, sparse = F, options=list(sd_noise=1))
  X <- sim$X
  y <- sim$y
  coefs.true <- sim$beta_true
  method <- 1
  # Generate regularisation parameter
  lambda_preset <- L0glm:::preset.lambda(IC = crit, p=ncol(X), n=nrow(X))
  start = L0glm:::init_coefs(X=X, y=y, mu=y+0.1, family=family,
                     lambda=lambda_preset) # component-wise initialisation
  for (i in 1:length(methods)) {
    solver = solvers[[i]]
    method = methods[[i]]
    solver_name = solver_names[[i]]
    
    # Fit model and measure time
    fit.time <- system.time(
      L0glm <- L0glm.fit(
        X=X, y=y,
        start=start, # NULL uses component-wise initialisation above, alt option: rep(1, ncol(X)) 
        warm.start = 1,
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
    metrics <- L0glm:::get.metrics(
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

warnings()
# no warnings

# Reshape data for plotting
results_melt <- results %>% 
  tidyr::pivot_longer(c(Time, Precision, Recall, Accuracy, Specificity, 
                        F1, Estimated_over_true_support_size, R2, deviance), 
                      names_to = "metric", values_to = "value") %>% 
  mutate(method = factor(method), solver = factor(solver))
# results_melt$Solver = factor(results_melt$Solver, 
#                              levels=rev(c("eigen Cholesky LLT", "eigen CG", "eigen LSCG", "eigen QR", "osqp", 
#                                           # "nnls", "nnls_noreg", 
#                                           "bigGlm"
#                                           )),
#                              labels=rev(c("L0glm_eigen Cholesky LLT", "L0glm_eigen CG", "L0glm_eigen LSCG", "L0glm_eigen QR", "L0glm_osqp", 
#                                           # "L0glm_nnls", "nnls_noreg", 
#                                           "L0glm_bigGlm"
#                                           )))

# Plot
ggplot(results_melt, aes(x = Size, y = value, colour=interaction(Solver))) +
  geom_line(aes(linewidth=I(1))) +
  facet_wrap(~metric, scales = "free_y") +
  theme_minimal() +
  labs(title = "Benchmarks different L0glm solvers for varying problem sizes\n(deconvolution of blurred spike train with Poisson noise, dense case)", color = "method / solver") +
  theme(legend.position = "right") +
  xlab("Problem size (n=p)") +
  scale_x_continuous(trans="log10") +
  scale_color_manual(values=rev(c("green4", "magenta", "red", "orange", "cyan2", 
                                  "blue", "yellow2", 
                                  "green", "purple"))) 
ggsave(file="../benchmarks/plots/timings/timings_metrics_dense_poisson.png", width=7, height=5, dpi=600)
ggsave(file="../benchmarks/plots/timings/timings_metrics_dense_poisson.pdf", width=7, height=5, dpi=600)

ggplot(results_melt[results_melt$metric=="Time",], aes(x = Size, y = value, colour=interaction(Solver))) +
  geom_line(aes(linewidth=I(1))) +
  # facet_wrap(~metric, scales = "free_y") +
  theme_minimal() +
  labs(title = "Timings different L0glm solvers for varying problem sizes\n(deconvolution of blurred spike train with Poisson noise, dense case)", color = "method / solver") +
  theme(legend.position = "right") +
  scale_color_manual(values=rev(c("green4", "magenta", "red", "orange", "cyan2", 
                                  "blue", "yellow2", 
                                  "green", "purple"))) +  
  ylab("Time (s)") + xlab("Problem size (n=p)") +
  scale_y_continuous(trans="log10") + scale_x_continuous(trans="log10")
ggsave(file="../benchmarks/plots/timings/timings_dense_poisson.png", width=7, height=5, dpi=600)
ggsave(file="../benchmarks/plots/timings/timings_dense_poisson.pdf", width=7, height=5, dpi=600)

# save.image("~/GitHub/L0glm/benchmarks/environment_benchmarks timings in function of problem size.RData")

results_melt[results_melt$metric=="Time"&results_melt$Size==10000,]
# Solver                                 Size metric value method solver
# 1 osqp                                  10000 Time    9.18 9      eigen 
# 2 bigGlm                                10000 Time    6.44 9      eigen 
# 3 eigen LSCG                            10000 Time   15.8  9      eigen 
# 4 eigen PLSCG                           10000 Time    6.31 9      eigen 
# 5 eigen CG                              10000 Time   10.2  9      eigen 
# 6 eigen PCG diagonal precond            10000 Time    5.77 9      eigen 
# 7 eigen Cholesky LLT                    10000 Time    6.14 9      eigen 
# 8 eigen Cholesky LDLT                   10000 Time    5.53 9      eigen 
# 9 eigen PCG incomplete Cholesky precond 10000 Time    6.14 9      eigen 



### 2.1.2 SINGLE CHANNEL - GAUSSIAN NOISE ####

family = gaussian(identity)

# Initialize an empty data frame to store results
results <- data.frame()

# Loop over sizes, solvers and methods
for (n in sizes) {
  # Simulate data
  sim <- simul(method = 1, n = n, p = n, sparse = F, options=list(sd_noise=1))
  X <- sim$X
  y <- sim$y
  coefs.true <- sim$beta_true
  method <- 1
  # Generate regularisation parameter
  lambda_preset <- L0glm:::preset.lambda(IC = crit, p=ncol(X), n=nrow(X))
  start = L0glm:::init_coefs(X=X, y=y, mu=y, family=family, lambda=lambda_preset) # rescaled component-wise initialisation
  for (i in 1:length(methods)) {
    solver = solvers[[i]]
    method = methods[[i]]
    solver_name = solver_names[[i]]
    
    # Fit model and measure time
    fit.time <- system.time(
      L0glm <- L0glm.fit(
        X=X, y=y, 
        start=start, # NULL uses component-wise initialisation above, alt option: rep(1, ncol(X)) 
        warm.start=1,
        family = family, lambda = lambda_preset, 
        nonnegative = TRUE, relaxed = FALSE, 
        solver = solver, method = method, 
        # dispersion = 1,
        control.solver = do.call(control.solver.gen, 
                                 list(solver = solver, 
                                      method = method)))
    )[3]
    print(solver_name)
    print(n)
    print(fit.time)
    # Calculate metrics
    metrics <- L0glm:::get.metrics(
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

warnings()
# no warnings

# Reshape data for plotting
results_melt2 <- results %>% 
  tidyr::pivot_longer(c(Time, Precision, Recall, Accuracy, Specificity, F1,
                        Estimated_over_true_support_size, R2, deviance), names_to = "metric", values_to = "value") %>% 
  mutate(method = factor(method), solver = factor(solver))
# results_melt2$Solver = factor(results_melt2$Solver, 
#                               levels=rev(c("eigen Cholesky LLT", "eigen CG", "eigen LSCG", "eigen QR", "osqp", 
#                                            #"nnls", "nnls_noreg", 
#                                            "bigGlm"
#                                            )),
#                               labels=rev(c("L0glm_eigen Cholesky LLT", "L0glm_eigen CG", "L0glm_eigen LSCG", "L0glm_eigen QR", "L0glm_osqp", 
#                                            #"L0glm_nnls", "nnls_noreg", 
#                                            "L0glm_bigGlm"
#                                            )))

# Plot
ggplot(results_melt2, aes(x = Size, y = value, colour=interaction(Solver))) +
  geom_line() +
  facet_wrap(~metric, scales = "free_y") +
  theme_minimal() +
  labs(title = "Benchmarks different L0glm solvers for varying problem sizes\n(deconvolution of blurred spike train with Gaussian noise, dense case)", color = "method / solver") +
  theme(legend.position = "right") +
  xlab("Problem size (n=p)") +
  scale_x_continuous(trans="log10") +
  scale_color_manual(values=rev(c("green4", "magenta", "red", "orange", "cyan2", 
                                  "blue", "yellow2", 
                                  "green", "purple"))) 
ggsave(file="../benchmarks/plots/timings/timings_metrics_dense_gaussian.png", width=7, height=5, dpi=600)
ggsave(file="../benchmarks/plots/timings/timings_metrics_dense_gaussian.pdf", width=7, height=5, dpi=600)

# here bigGlm solver is markedly worse than the rest
# TODO check if this can be improved

ggplot(results_melt2[results_melt2$metric=="Time",], aes(x = Size, y = value, 
                                                         colour=interaction(Solver))) +
  geom_line() +
  # facet_wrap(~metric, scales = "free_y") +
  theme_minimal() +
  labs(title = "Timings different L0glm solvers for varying problem sizes\n(deconvolution of blurred spike train with Gaussian noise, dense case)", color = "method / solver") +
  theme(legend.position = "right") +
  scale_color_manual(values=rev(c("green4", "magenta", "red", "orange", "cyan2", 
                                  "blue", "yellow2", 
                                  "green", "purple"))) +
  ylab("Time (s)") + xlab("Problem size (n=p)") +
  scale_y_continuous(trans="log10") + scale_x_continuous(trans="log10")
ggsave(file="../benchmarks/plots/timings/timings_dense_gaussian.png", width=7, height=5, dpi=600)
ggsave(file="../benchmarks/plots/timings/timings_dense_gaussian.pdf", width=7, height=5, dpi=600)

# save.image("~/GitHub/L0glm/benchmarks/environment_benchmarks timings in function of problem size.RData")

results_melt2[results_melt2$metric=="Time"&results_melt2$Size==10000,]
#  Solver                                 Size metric value method solver
# 1 osqp                                  10000 Time   52.2  9      eigen 
# 2 bigGlm                                10000 Time    9.80 9      eigen 
# 3 eigen LSCG                            10000 Time    8.74 9      eigen 
# 4 eigen PLSCG                           10000 Time    8    9      eigen 
# 5 eigen CG                              10000 Time    5.80 9      eigen 
# 6 eigen PCG diagonal precond            10000 Time    5.61 9      eigen 
# 7 eigen Cholesky LLT                    10000 Time    3.32 9      eigen 
# 8 eigen Cholesky LDLT                   10000 Time    3.78 9      eigen 
# 9 eigen PCG incomplete Cholesky precond 10000 Time    2.87 9      eigen