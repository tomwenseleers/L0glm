# L0glm 0.2.0

## Sparse generalized linear models

* L0 and grouped L0 variable selection through adaptive-ridge IRLS, with
  relaxed refitting on the selected support.
* Standard GLM families and links, negative binomial with fixed or estimated
  shape, joint multinomial models in CLR coordinates, and Cox proportional
  hazards models with direct or matrix-free quadratic updates.
* NB1, common-scale Gamma and censored Gamma signal models; quasi families
  expose scaled likelihood kernels and quasi-information criteria.

## Model specification and selection

* Observation weights, offsets, intercepts, mandatory predictors, nonnegative
  coefficients, normalization, starting values and optional validation stopping.
* Shared observation vectors stay compact across channels in fitting, native
  Poisson, tuning and CV. Prediction supports shared or channel-specific offsets.
* Automatic information-criterion presets, lambda-path tuning, n-fold CV,
  one-standard-error models and coefficient-path plots. `preset.lambda()`
  returns AIC, BIC, GIC, EBIC or mBIC penalty values directly.
* Constrained updates for nonnegative compulsory terms and relaxed-fit boundary
  checks use the specified likelihood.

## High-dimensional and multichannel computation

* Dense and sparse designs, tall and wide systems, Eigen solvers and a native
  OpenMP backend for nonnegative identity-link Poisson models.
* Sure Independence Screening (SIS) and Iterative Sure Independence Screening
  (ISIS), with compiled OpenMP screening. Automatic ISIS uses
  `p_eligible > 10*n`, excluding mandatory predictors and rows locked out by
  zero starts. Tuning applies the same rule to each training fold.
* Single and multichannel fitting, shared predictor-row selection and optional
  Poisson outcome aggregation followed by original-channel refitting.
* Size-adjusted Poisson and scalar-family likelihood clustering with native
  distances, pooled-likelihood merge costs and IC-scored cuts.

## Simulation, comparison and interpretation

* `sparseglm()` compares supported sparse-regression methods with common
  likelihood, prediction, coefficient and support-recovery metrics and timings.
* Prediction, fitted-value, summary, likelihood and information-criterion
  methods, including within-row multichannel sparsity metrics.
* Generic simulations, fixed-PSF deconvolution examples, visualization tools,
  reproducible benchmarks and a mathematical algorithm appendix.
