\dontrun{
# Simulate data
sim_data <- simulate(n = 100, p = 10, family = "gaussian")

# Fit using method 2 (glmnet)
fit <- fit.GLM(sim = sim_data, method = 2)

# View coefficients
print(fit)
}
