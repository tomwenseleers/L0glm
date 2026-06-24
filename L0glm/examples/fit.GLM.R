\dontrun{
# Simulate data
sim.gaus.sc <- simul(method = 1, n = n, p = p, s = npeaks/p, seed = 2)

# Fit using method 2 (glmnet)
fit <- fit.GLM(sim = sim.gaus.sc, method = 2, family = gaussian(identity))

# View coefficients
print(fit)
}
