# 1. SINGLE CHANNEL ####
# Additional libraries
library(ggplot2)
library(dplyr)
library(export)

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
  tidyr::pivot_longer(c(Time, Precision, Recall, Accuracy, Specificity, BIC, F1,Estimated_over_true_support_size, R2), names_to = "metric", values_to = "value") %>% 
  mutate(method = factor(method), solver = factor(solver))
results_melt$Solver = factor(results_melt$Solver, 
                             levels=c("glmnet", "osqp", "Chol", "CG", "LSCG", "QR"),
                             labels=c("bigGlm/glmnet", "osqp", "eigen Cholesky LLT", "eigen CG", "eigen LSCG", "eigen QR"))

# Plot
ggplot(results_melt, aes(x = Size, y = value, colour=interaction(Solver))) +
  geom_line() +
  facet_wrap(~metric, scales = "free_y") +
  theme_minimal() +
  labs(color = "Solver/Method") +
  theme(legend.position = "right") +
  scale_color_manual(values=c("blue", "magenta", "red", "orange", "green", "cyan3")) 
ggsave(file="..//benchmarks//plots//timings//timings_metrics.png", width=7, height=5, dpi=600)
ggsave(file="../benchmarks/plots/timings/timings_metrics.pdf", width=7, height=5, dpi=600)
# NOTE: BIC makes no sense as that can't be compared across different models
# LSCG solver could probably be fixed via some thresholding or using other settings for that one

ggplot(results_melt[results_melt$metric=="Time",], aes(x = Size, y = value, colour=interaction(Solver))) +
  geom_line() +
  facet_wrap(~metric, scales = "free_y") +
  theme_minimal() +
  labs(color = "Solver/Method") +
  theme(legend.position = "right") +
  scale_color_manual(values=c("blue", "magenta", "red", "orange", "green", "cyan3")) +
  ylab("Time (s)")
ggsave(file="..//benchmarks//plots//timings//timings.png", width=7, height=5, dpi=600)
ggsave(file="../benchmarks/plots/timings/timings.pdf", width=7, height=5, dpi=600)

# R version 4.0.2 (2020-06-22)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19045)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_United Kingdom.1252  LC_CTYPE=English_United Kingdom.1252    LC_MONETARY=English_United Kingdom.1252 LC_NUMERIC=C                           
# [5] LC_TIME=English_United Kingdom.1252    
# 
# attached base packages:
#   [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] glmnet_4.1-4         qlcMatrix_0.9.7      sparsesvd_0.2        slam_0.1-47          dplyr_1.0.7          ggplot2_3.3.5        L0glm_0.1.0          Matrix_1.3-4        
# [9] RevoUtils_11.0.2     RevoUtilsMath_11.0.0
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.7       compiler_4.0.2   pillar_1.6.3     iterators_1.0.12 tools_4.0.2      digest_0.6.28    boot_1.3-28      docopt_0.7.1     gtable_0.3.0     lifecycle_1.0.1 
# [11] tibble_3.1.5     lattice_0.20-45  pkgconfig_2.0.3  rlang_1.0.5      foreach_1.5.0    cli_3.0.1        rstudioapi_0.13  osqp_0.6.0.3     withr_2.4.2      generics_0.1.0  
# [21] vctrs_0.3.8      grid_4.0.2       tidyselect_1.1.1 glue_1.4.2       R6_2.5.1         fansi_0.5.0      survival_3.2-13  farver_2.1.0     purrr_0.3.4      tidyr_1.1.4     
# [31] magrittr_2.0.1   scales_1.1.1     codetools_0.2-18 ellipsis_0.3.2   nnls_1.4         splines_4.0.2    colorspace_2.0-2 shape_1.4.6      labeling_0.4.2   utf8_1.2.2      
# [41] munsell_0.5.0    crayon_1.4.1    

# running Intel MKL using 20 cores
# on machine Intel Xeon CPU E5-2650 v3 @2.3 GHz with 40 logical cores 
