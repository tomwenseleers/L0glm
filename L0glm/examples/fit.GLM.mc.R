\dontrun{
# Simulate multichannel data
sim_data <- simulate(n = 100, p = 10, m = 5, family = "gaussian")

# Fit using method 1 (abess Best Subset of Group Selection)
fit <- fit.GLM.mc(sim = sim_data, method = 1)

# View coefficients
print(fit)

# Fit using method 2 (glmnet Multivariate Group LASSO)
fit <- fit.GLM.mc(sim = sim_data, method = 2)

}
