# L0glm 0.2.0.9002

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
  one-standard-error models and coefficient-path plots. `presetLambda()`
  returns AIC, BIC, GIC, EBIC or mBIC penalty values directly.
  Ordinary and special-family IC paths use the parameter counts reported by
  the S3 criterion methods, including multinomial zero-sum constraints.
* Constrained updates for nonnegative compulsory terms and relaxed-fit boundary
  checks use the specified likelihood.
* Optional final support refinement tests additions, deletions and small
  exchanges after SIS/ISIS and original-channel refitting. It supports the
  ordinary and special families, respects starting exclusions and constraints,
  and is available within information-criterion tuning and CV training fits.
  The default `local.search="auto"` applies it to initial supports of at most
  50 coefficient dimensions, preserving one-iteration ridge pilots and
  unpenalized fits. `FALSE` disables refinement.
  Screened searches operate directly on the retained candidate pool while
  preserving original-dimension coefficients and information criteria.
* The shared native quadratic solver offers checked Cholesky solves with
  iterative refinement and a rank-revealing fallback for profile-shape updates.
  Optional feasible starting values accelerate repeated quadratic subproblems.
* Conditional Poisson refits reuse their preceding active sets and cache fixed
  gradient scales; sparse weighted products accelerate local profile designs.

## High-dimensional and multichannel computation

* Dense and sparse designs and tall or wide systems. A native scalar-family
  driver performs initialization, dispersion updates and grouped/relaxed passes.
  All Eigen solvers share their kernels with the R interfaces and parallelize
  channels through OpenMP. OSQP and `bigGlm` use the native driver with their
  R-facing solver interfaces and independent workers for concurrent channels.
* Nonnegative identity-link Poisson LLT fits use the specialized `l0nnPois`
  backend, including native validation stopping.
* Sure Independence Screening (SIS) and Iterative Sure Independence Screening
  (ISIS), with compiled OpenMP screening. Automatic ISIS uses
  `p_eligible > 10*n`, excluding mandatory predictors and rows locked out by
  zero starts. Tuning applies the same rule to each training fold.
* Single and multichannel fitting, shared predictor-row selection and optional
  Poisson outcome aggregation followed by original-channel refitting.
* Size-adjusted Poisson and scalar-family likelihood clustering with native
  distances, pooled-likelihood merge costs and IC-scored cuts.

## Simulation, comparison and interpretation

* `sparseGLM()` compares supported sparse-regression methods with common
  likelihood, prediction, coefficient and support-recovery metrics and timings.
  Shared SIS/ISIS options screen every adapter, with method-specific ISIS
  offsets, original-dimension coefficients and training-only screening for CV.
* Prediction, fitted-value, summary, likelihood and information-criterion
  methods, including within-row multichannel sparsity metrics.
* Generic simulations, fixed-PSF deconvolution examples, visualization tools,
  reproducible benchmarks and a mathematical algorithm appendix.
