# 1. SINGLE CHANNEL ####
# Additional libraries
library(ggplot2)
library(dplyr)

# Simulation sizes
sizes <- c(100, 500, 1000, 2000, 5000, 10000, 20000)

# Solvers and methods
solvers <- c("osqp", "glmnet", "eigen", "eigen", "eigen", "eigen")
family = poisson(identity)

# Initialize an empty data frame to store results
results <- data.frame()

# Loop over sizes, solvers and methods
for (n in sizes) {
  # Simulate data
  sim.pois <- simulate(method = 2, n = n, p = n)
  X <- sim.pois$X
  y <- sim.pois$y
  coefs.true <- sim.pois$beta_true
  method <- 1
  for (solver in solvers) {
    # Generate regularisation parameter
    lambda_preset <- preset.lambda(IC = "bic", y = y, X = X, family = family)
    
    # Fit model and measure time
    fit.time <- system.time(
      L0glm <- L0glm.fit(
        X=X, y=y, family = poisson(identity), lambda = lambda_preset, nonnegative = TRUE, relaxed = FALSE, 
        solver = solver, method = method, control.solver = do.call(control.solver.gen, list(solver = solver, method = method, rowaug = F)))
    )[3]
    print(solver)
    print(n)
    print(fit.time)
    # Calculate metrics
    metrics <- get.metrics(
      beta = L0glm$coefficients, X = X, y = y, family = gaussian(identity), beta_true = as.vector(coefs.true)
    )
    
    if(method == 1){
      meth <- "QR"
    } else if(method == 3){
      meth <- "LSCG"
    } else if(method == 5){
      meth <- "CG"
    } else if(method == 7){
      meth <- "Chol"
    }
    # Save results to the data frame
    results <- rbind(
      results, 
      data.frame(
        Solver = ifelse(solver=="eigen", meth, solver), 
        Size = n, 
        Time = fit.time,
        Precision = metrics$precision, 
        Recall = metrics$recall, 
        Specificity = metrics$specificity, 
        Accuracy = metrics$accuracy, 
        BIC = metrics$bic, 
        F1 = metrics$F1_score,
        R2 = metrics$R2_coefs,
        Estimated_over_true_support_size = metrics$estimated_over_true_support_size
      )
    )
    if (solver == "eigen") method = method + 2
  }
}

# Reshape data for plotting
results_melt <- results %>% 
  tidyr::pivot_longer(c(time, precision, recall, accuracy, specificity, BIC, F1, estimated_over_true_support_size, R2), names_to = "metric", values_to = "value") %>% 
  mutate(method = factor(method), solver = factor(solver))

# Plot
ggplot(results_melt, aes(x = size, y = value, color = interaction(solver))) +
  geom_line() +
  facet_wrap(~metric, scales = "free_y") +
  theme_minimal() +
  labs(color = "Solver/Method") +
  theme(legend.position = "bottom")

