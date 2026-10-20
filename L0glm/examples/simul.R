# Blurred spike trains with a known PSF, sparse peaks and an optional baseline.
single <- simul("deconvolution", n=200, p=200, n_nonzero=20,
  family="poisson", intercept=1, seed=2)
multi <- simul("deconvolution", n=100, p=120, n_nonzero=10,
  family="gaussian", m=4, intercept=c(1, 2, 3, 4), seed=2,
  options=list(mean_beta=100, sd_noise=2))
dim(multi$beta_true)                  # 120 slopes plus an intercept, four channels
if (requireNamespace("mgcv", quietly=TRUE)) {
  nb_spikes <- simul("deconvolution", n=100, p=100, n_nonzero=10,
    family="nb", m=4, intercept=1, options=list(theta=5))
}

# Standard normal predictors, arbitrary n/p, and signed sparse coefficients.
regression <- simul("normal", n=150, p=30, n_nonzero=4, intercept=1.5,
  options=list(beta_range=c(-2, 2), error_sd=.5), seed=10)
fit <- L0glm.fit(regression$X, regression$y, family=regression$family,
  intercept=regression$intercept, nonnegative=FALSE, parallelize=FALSE)
L0glm.metrics(fit, beta_true=regression$beta_true)[c("gic", "F1_score", "MCC")]

# sparseglm reads family, intercept and coefficient truth from the simulation.
comparison <- sparseglm(regression, method=c("L0glm", "lasso"), nlambda=5)
L0glm.metrics(comparison)[c("method", "seconds", "gic", "F1_score", "MCC")]

# ABESS-style correlated predictors and coefficient magnitudes.
abess_data <- simul("abess", n=200, p=100, n_nonzero=5,
  family="binomial", intercept=-1, options=list(cortype=2, rho=.5))

# Best-subset pattern 3: decreasing effects, population SNR=2, validation data.
subset_data <- simul("bestsubset", n=100, p=50, n_nonzero=5, intercept=2,
  options=list(beta.type=3, rho=.5, snr=2, nval=50))
subset_data$dispersion
dim(subset_data$Xval)

# Other standard links and special families use the same normal-design setup.
probit_data <- simul("normal", n=100, p=20, n_nonzero=3,
  family=binomial("probit"), intercept=-.5)
gamma_data <- simul("normal", n=100, p=20, n_nonzero=3,
  family=Gamma("log"), intercept=1, options=list(dispersion=.5))
if (requireNamespace("mgcv", quietly=TRUE) && requireNamespace("survival", quietly=TRUE)) {
  nb_data <- simul("normal", n=150, p=20, n_nonzero=3,
    family=mgcv::nb(), intercept=1, options=list(theta=3))
  multinomial_data <- simul("normal", n=150, p=20, n_nonzero=3,
    family=mgcv::multinom(K=2), intercept=c(.3, 0, -.3))
  cox_data <- simul("normal", n=150, p=20, n_nonzero=3,
    family=mgcv::cox.ph(), intercept=log(2)) # doubles the baseline hazard
  cox_data$intercept                    # FALSE: no estimable Cox intercept
}
