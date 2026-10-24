# L0glm 0.1.0

## Initial version

* Fast sparse generalized linear models with adaptive ridge and IRLS, including
  Gaussian, Poisson, binomial, Gamma and inverse Gaussian families and their links.
* Negative binomial regression with estimated shape, joint multinomial regression
  in centered log-ratio coordinates, and Cox proportional hazards models.
* Conditional boundary polishing for relaxed nonnegative models removes
  likelihood-supported zero coefficients without a coefficient-size cutoff.
  Native Gaussian/Poisson kernels retain dense/sparse storage and Poisson OpenMP;
  other GLMs, negative binomial and Cox use their family-specific objectives.
* Single and multichannel fitting, grouped predictor-row selection, and relaxed
  refitting for accurate coefficient estimates on the selected support.
* Dense and sparse matrices, tall and wide designs, fast Eigen solvers and
  a native OpenMP backend for nonnegative identity-link Poisson models.
* Sure Independence Screening (SIS) and Iterative Sure Independence Screening
  (ISIS), with compiled OpenMP screening and automatic ISIS when
  `p_eligible > 10*n`, excluding mandatory terms and predictors excluded by
  zero starts in every response channel. The same rule uses the training row
  count in each `L0glm.tune()` fit.
  Screening uses physical cores by default, with explicit thread-count overrides.
* Observation weights, offsets, intercepts, mandatory predictors, normalization,
  starting values, nonnegativity constraints and validation-based early stopping.
* Information-criterion presets and lambda-path tuning, five-fold cross-validation,
  one-standard-error models and coefficient-path plots.
* `sparseglm()` comparisons with supported sparse-regression backends, common
  likelihood and support-recovery metrics, prediction and summary methods,
  mirrored deconvolution plots and multichannel heatmaps.
* Reproducible simulations, large-scale examples and a documentation website with
  tutorials, measured benchmarks and a mathematical algorithm appendix.

* Size-adjusted Poisson channel distances use a native OpenMP kernel with
  packed distance storage and scalar, row or channel-specific weights.
