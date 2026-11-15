# L0glm 0.2.0.9002

## Sparse generalized linear models

* `controlSolverGen(nonnegative.method=...)` selects `"clamp"` or
  `"constrained"` working ridge solves across Eigen, OSQP and `bigGlm`.
  Defaults are `"clamp"` for Eigen and `"constrained"` for OSQP/`bigGlm`.
  Bounds, free intercepts, normalization, prior means and numerical guards
  are preserved; Eigen constrained polishing reuses eligible Gram systems.

* Gram and factor workspaces are checked against a configurable memory budget.
  CG/LSCG can use implicit ridge products with a fresh coefficient-space residual
  check; oversized incomplete Cholesky switches to diagonal preconditioning.
  Dispersion pilots support matrix-free solves, and LLT fallbacks are guarded.
  Fixed-weight multichannel fits share a common immutable Gram matrix.
  The same allocation checks cover serial special-family solver adapters and
  adaptive-ridge dispersion pilots.
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
  Estimated-negative-binomial IRLS, shape updates and channel coordination use
  C++/OpenMP for recognized links. Working designs are reused across passes;
  per-channel NB shapes are copied without serializing caller workspaces.
  Exact joint multinomial IRLS coordinates penalties, backtracking, CLR masks
  and convergence natively; independent Fisher-information blocks use OpenMP.
  Wide systems retain observation-space solves and CLR entrywise sparsity.
  Cox initialization and IRLS also use native coordination, cached Breslow
  risk statistics and reusable matrix-free scratch buffers. Small direct
  systems retain LLT, and risk-information columns can run in OpenMP.
* NB1, common-scale Gamma and censored Gamma signal models; quasi families
  expose scaled likelihood kernels and quasi-information criteria.

## Model specification and selection

* Fixed-noise tuning endpoints use certified nonzero probes to avoid completing
  unnecessary local searches. Returned grid fits retain full refinement,
  independent starts and unchanged scoring. Correlated Poisson path examples
  illustrate the common `start="pilot"` option.

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
  All-zero coefficient paths share one plotted line, including multinomial
  contrasts and zero multichannel means, preserving colours and support counts.
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
  R-facing solver interfaces on the main R thread. `ncores` controls native
  outcome-channel threading, with one thread for serial execution. Custom R
  family/link callbacks stay serial, and no process workers are created.
  Serial callbacks and solver bridges use `ncores` as their default BLAS budget;
  explicit `nthreads` overrides it. OpenMP channel fits limit BLAS to one
  thread per channel, restoring the caller's settings on return.
* Native fits share one prepared working design across outcome channels.
  Automatic sparse routing compresses large exact-zero dense designs once,
  retaining the original input/output storage and explicit dense routing.
  Selected designs, fixed-weight crossproducts and sparse symbolic Cholesky
  patterns are reused without changing the fitted objective or support rules.
  Required sparse external-solver boundaries also reuse the packed design.
* CG defaults to incomplete Cholesky preconditioning for dense and sparse
  systems. `controlSolverGen("eigen", 3L, precond=1L)` selects diagonal
  preconditioning; `precond=0L` disables it. Dense matrix operators remain dense.
* LSCG defaults to incomplete Cholesky of the working normal matrix and a
  tolerance of `1e-12`; `controlSolverGen("eigen", 4L, precond=1L)` selects
  diagonal preconditioning. An inner-guess guard chooses zero when its normal
  residual is smaller than that of the previous working solution, protecting
  solves with large changes in adaptive ridge weights. It checks every
  right-hand side using a fresh residual and can reuse its factorization for
  one restart within the iteration limit. Already formed Gram and Woodbury
  systems use CG directly. Inexact or nonfinite solutions use Cholesky LLT in
  both R and native drivers.
* Nonnegative identity-link Poisson LLT fits automatically use a specialized
  internal backend, including native validation stopping. `l0glmFit()` is
  the public interface for every supported family.
* `controlSolverGen(..., nonnegative.method="constrained")` optionally polishes infeasible
  nonnegative Eigen ridge proposals using batched constrained Cholesky solves.
  It reuses existing predictor Gram matrices, including LLT and LDLT paths,
  supports the native GLM and specialized Poisson drivers, and is off by default. Rectangular
  and Woodbury solves retain their usual path. A fresh KKT check and a
  nonincreasing quadratic objective are required for acceptance.
* OSQP uses diagonal equilibration and absolute and relative tolerances of
  `1e-9`. Tall and square OSQP and `bigGlm` ridge steps receive a shared
  constrained Cholesky polish, checked against the box KKT conditions and the
  augmented least-squares objective. Wide designs retain coordinate descent
  without an additional predictor-space Gram allocation.
  Native callbacks to OSQP and `bigGlm` keep their R arguments protected across
  allocation and garbage collection for dense and sparse working designs.
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
