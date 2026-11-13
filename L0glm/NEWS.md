# L0glm 0.2.0.9002

## Sparse generalized linear models

* L0 and grouped L0 variable selection through adaptive-ridge IRLS, with
  criterion-checked support refinement. Large multichannel fits compare one
  unpenalized selected-support refit; large single-channel fits retain their
  adaptive penalties. Small supports allow local additions, deletions and swaps.
  Small ISIS candidates also use guarded refits to improve residual screening.
  `local.search=FALSE` retains pure adaptive-ridge selection. Explicit `lambda=0`
  requests an unpenalized model. Internal support refits check information and
  residual degrees of freedom.
  Independent channels retain their own selected supports.
* A common guarded Pearson rule updates ordinary dispersion and
  negative-binomial shape. Model criteria score the fitted noise parameters
  and count free coefficient dimensions. Eligible wide multitask Gaussian
  fits use an adaptive-ridge pilot and default to the BIC preset. A short
  Poisson pilot initializes eligible negative-binomial fits, and constrained
  unpenalized refits can fall back to a projected likelihood descent step.
* Standard GLM families and links, negative binomial with fixed or estimated
  shape, joint multinomial models in CLR coordinates, and Cox proportional
  hazards models with direct or matrix-free quadratic updates.
* NB1, common-scale Gamma and censored Gamma signal models; quasi families
  expose scaled likelihood kernels and quasi-information criteria.

## Model specification and selection

* Predictor-variable groups use norm-based block L0 penalties, with group 0
  unpenalized, dimension-calibrated IC presets, explicit numeric group costs,
  group-aware thresholding, guarded support refits and cross-validation.
  Recognized scalar-family predictor groups use native IRLS with OpenMP channel
  solves and predictions, preserving initialization, penalties and convergence;
  custom/special-family and validation adapters retain their structured driver.
  Dense Gaussian groups reuse bounded weighted cross-products; standard group
  fits cache current/trial predictions and skip fixed-scale updates. Dense ridge
  diagonals are updated directly, preserving sparse and Woodbury solver paths.
* Staged interaction selection generates products among selected main effects:
  `interactions=1L` for pairwise and `2L` for up to three-way terms. Prediction
  reconstructs products from original columns, and CV repeats candidate
  generation independently within every training fold.
* Observation weights, offsets, intercepts, mandatory predictors, nonnegative
  coefficients, normalization, starting values and optional validation stopping.
* Optional conditional posterior summaries and sparse coefficient sampling,
  using a flat prior after an accepted full support refit and a frozen adaptive
  prior otherwise. Nonnegative and multinomial CLR constraints are respected.
  Selected-support standard errors and case-bootstrap intervals provide
  complementary summaries, including reselection frequencies. Ordinary GLMs,
  negative binomial, multinomial and Cox models are supported.
* Shared observation vectors stay compact across channels in fitting, native
  Poisson, tuning and CV. Prediction supports shared or channel-specific offsets.

* Automatic information-criterion presets, lambda-path tuning, n-fold CV,
  full-data minimum and one-standard-error models, and coefficient-path plots.
  Iteration limits are scalar; independent validation data can select the
  best predictive checkpoint during fitting. `presetLambda()`
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
* CG defaults to incomplete Cholesky preconditioning for dense and sparse
  systems. `controlSolverGen("eigen", 3L, precond=1L)` selects diagonal
  preconditioning; `precond=0L` disables it. Dense matrix operators remain dense.
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

## Memory-bounded data access

* Read-only native screening of double-precision bigmemory backing files avoids
  intermediate R design copies, with deterministic OpenMP reductions and exact
  agreement with resident screening across ordinary and special families.
* `l0glmData()` provides lazy subsets of bigmemory, delayed/HDF5 single-cell
  assays and PLINK BED/PGEN designs. VCF input supports an explicit PLINK 2
  conversion cache, including expected dosages.
* Large file-backed Gaussian identity ISIS fits can retain a bounded candidate
  pool, preserving the global ISIS path and selecting cached or small-support
  Gaussian refinements only when original-p GIC improves. Repeated exact
  coefficients avoid a redundant global scan; `working.columns=0L` disables
  this cache. Restored descriptors attach once per fit, and pilot blocks respect
  the working-memory guard.
* All-family SIS/ISIS streams predictor columns; BED screening reduces packed
  two-bit genotypes directly with OpenMP. Retained-design fitting preserves
  original coefficient dimensions and criterion calibration. Training-only
  imputation is reused for prediction and relearned inside CV and bootstraps.
* Ward outcome clustering streams row blocks, including its cut losses and
  aggregation. Explicit memory guards bound working matrices and outputs.
* Native file-backed LASSO through optional biglasso reuses double bigmemory
  designs across training/CV row views. Gaussian, logistic, grouped Gaussian
  and Breslow Cox paths use the common likelihood scoring, IC selection and CV.

## Simulation, comparison and interpretation

* `sparseGLM()` compares supported sparse-regression methods with common
  likelihood, prediction, coefficient and support-recovery metrics and timings.
  Shared SIS/ISIS options screen every adapter, with method-specific ISIS
  offsets, original-dimension coefficients and training-only screening for CV.
* Prediction, fitted-value, summary, likelihood and information-criterion
  methods, including within-row multichannel sparsity metrics.
* Generic simulations, fixed-PSF deconvolution examples, visualization tools,
  reproducible benchmarks and a mathematical algorithm appendix.
