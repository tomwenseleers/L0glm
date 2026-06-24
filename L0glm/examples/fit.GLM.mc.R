\dontrun{
# Simulate multichannel data
n <- p <- 1000
npeaks <- 100
sim.gaus.mc <- simul(method = 3, n = n, p = p, s = npeaks/p, seed = 1)

# Fit using method 1 (abess Best Subset of Group Selection)
fit <- fit.GLM.mc(sim = sim.gaus.mc, method = 1,
                  family = gaussian(identity))

# View coefficients
print(fit)

# Fit using method 2 (glmnet Multivariate Group LASSO)
fit <- fit.GLM.mc(sim = sim.gaus.mc, method = 2,
                  family = gaussian(identity))

}
