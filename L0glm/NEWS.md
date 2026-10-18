# L0glm 0.2.1

## Development changes

* Added `sparseglm()` to fit, tune and compare L0glm, ABESS, L0Learn, MCP, SCAD,
  LASSO, relaxed LASSO and both adaptive LASSO variants. Every candidate uses
  common likelihood/IC scoring or shared row-fold CV, with training-only
  adaptive ridge pilots. Results retain sparse original-dimension coefficients,
  paths, fold diagnostics, timings and explicit failure/availability statuses.
  `L0glm.metrics()` is now an S3 generic with a comparison-table method.
  Native grouped Gaussian, joint multinomial and Cox adapters are included;
  matrix channel weights are preserved using independent paths where needed.
  Identity-Poisson WLS approximations require an explicit opt-in. Retired the
  old benchmark drivers and `fit.GLM` entry points without compatibility wrappers,
  and migrated their examples to the new interface.

* Fits now inherit from `L0glm_fit`, providing `coef()`, `logLik()`, `nobs()`,
  standard `AIC()`/`BIC()` and S3 `GIC()`/`EBIC()`/`MBIC()` extraction for
  ordinary, native Poisson and special-family models. Original prior weights,
  offsets and dispersion-estimation flags are retained for correct on-demand
  scoring. `L0glm.metrics()` now also supports ordinary models and optional
  `beta_true` recovery scores. The help pages explain likelihood conventions,
  parameter counts, conditional within-row scores and removed diagnostics.

* Recovery metrics now report channel-entry confusion counts, precision, recall,
  F1 and MCC within correctly recovered nonmandatory rows of multichannel fits.
  The pooled scores are conditional on recovering the row; row-support metrics
  still account for missed and spurious rows. Coefficient summaries retain MAE,
  RMSE and through-origin R-squared. Selected-coefficient summaries, signed mean
  error, sMAPE, regression slope, support-size ratio, accuracy and specificity
  were removed. Historical result tables are preserved; retired diagnostics in
  benchmark scripts were removed and the solver benchmark now includes MCC/RMSE.

* Normalized dense Gaussian identity fits reuse fixed weighted cross-products
  for LLT/LDLT solves with at most 2,000 coefficients. The cache follows the
  shrinking active set; effectively wide and sparse systems retain their
  existing solver paths. Adaptive penalties and dispersion updates are unchanged.

* Normalized dense fits avoid repeated prediction-column copies when most
  coefficients remain active. Restricted-domain links retain their previous
  summation order, preserving the difficult Gamma square-root-link fit.

* Eigen CG now defaults to inner relative tolerance `1e-10`, improving agreement
  with LLT in the Gaussian PSF benchmark. Explicit solver tolerances are retained;
  LSCG keeps its default `1e-7`. This does not change LLT or the outer tolerance.

* Restored the fitting defaults `maxit = 25`, `tol = 1e-10`. Increase `maxit`
  explicitly for problems requiring more iterations; returned convergence
  diagnostics distinguish stopping at the iteration limit from convergence.
  Single-channel native fits no longer detect unused OpenMP core counts.
  Dispersion pilots use sparse row augmentation for large dense-stored designs
  containing mostly exact zeros, preserving the same ridge system and main
  solver storage.
* Numeric scalar `start = 1` now denotes ones in normalized coefficient
  coordinates, including a free centered intercept. All other scalars and
  explicit vectors/matrices retain original input units. This applies to both
  fitters, multichannel fits, and independently initialized tuning/CV fits.

* Added `normalize = "auto"` (default), `"none"`, `"l2"`, `"rms"`, and
  `"sd"`. RMS gives L2 norm `sqrt(n)`. Automatic initialization, numerical
  thresholds and convergence use compatible coordinates, while coefficients,
  explicit vector/matrix starts and ridge targets retain the original input units.
  Centering requires an unrestricted intercept: nonnegative intercepts remain
  constrained, and binomial/quasibinomial intercepts retain their exception.
* Wide Eigen and native Poisson systems use typed Woodbury solves, including
  unpenalized nuisance terms. Dense CG uses diagonal preconditioning; sparse CG
  keeps incomplete Cholesky. Dense/sparse design storage is preserved.
* Exact zero starts exclude coefficients separately in each channel and screen
  all-zero rows before fitting. Mandatory terms remain eligible with zero starts.
  Coefficients are returned as sparse matrices in the original dimensions.
  Signed multichannel starts pool absolute magnitudes without cancellation.
* Intercepts and `always.include` terms are unpenalized. GIC uses the original
  selectable predictor count before starting-mask or pilot support restriction.
* Estimated dispersion resumes updating after a support change. A virtually
  interpolating automatic pilot gets one retry with a moderate ridge penalty.
  Early stopping and per-channel/per-pass `convergence` diagnostics are returned.
* `L0glm.tune()` accepts a common numeric `start` or `start = "pilot"`.
  Independent pilot fits supply a shared seed; each requested lambda is fitted
  independently from it. Returned starts and `fit_args` allow direct replay.
  CV builds pilot seeds and preprocessing separately within each training fold.

* Removed the slow NNLS ridge solver and its package dependency. Use Eigen,
  `bigGlm`, or OSQP. Also removed the unused internal block-fitting helper.

* Both fitters accept scalar starts for every coefficient and response channel,
  including an added intercept. Numeric scalar one uses normalized coordinates;
  other scalar starts repeat in original input units.
* CV tuning returns the full-data `fit_1se` alongside `fit_best`, retaining
  both from the existing path without any additional fits.
* The full tuning range starts at `lambda_min = 2` by default; a smaller
  positive bound can be supplied explicitly. Zooming moves each endpoint
  toward the reference on the log scale. Grids retain unique endpoints and,
  when there are at least three points, the reference.
* CV tuning returns `lambda_1se`, the largest evaluated lambda with mean
  deviance within one SE of the minimum (using the SE at that minimum).
  Both panels mark it in blue alongside `lambda_best` and `lambda_preset`.
  Criterion tuning returns `NULL` for this field; `fit_best` is unchanged.
* Tuning plots now use `log10(lambda)` with a reversed axis in both panels, colored coefficient
  paths and a mathematical beta-hat label. Top axes show support sizes (row
  support for multichannel responses) at 16 equally spaced log-scale positions,
  using the nearest fitted models without adding any fits.
* Renamed `L0glm.tune(IC = ...)` to `L0glm.tune(tune = ...)`; the returned
  selection method is also named `tune`. Package and benchmark calls are updated.
* Added `tune = "cv"`, `nfolds = 5`, and optional `foldid` assignments. CV holds
  out rows jointly across channels, initializes every training fit independently,
  and scores weighted held-out deviance. Training penalties are scaled by their
  fraction of total prior weight to retain the full-data lambda scale. The plot
  shows mean deviance +/- one SE, the minimum and the GIC reference preset.
  Fold scores, assignments and training penalties are returned in `cv`.

* `L0glm.tune()` now initializes every lambda independently. A preset fit no
  longer excludes features from all other fits. The upper endpoint is checked
  by fitting it, expanded if necessary, then refined between verified nonzero
  and zero fits to within 5%; it is not a global-optimality bound.
  The fitted path coefficients, numerical endpoint and best model are returned.
* Automatic dispersion initialization is independent of the fitted lambda,
  using the GIC reference penalty for its permissive scale. This prevents a
  large lambda from cancelling itself through an inversely scaled dispersion.
* Tuning scores now include observation weights, offsets and known dispersion.
  Gaussian scales are profiled per channel when unknown. Corrected the extra
  scale charge in log likelihood and kept the predictor count separate from
  estimated scale parameters. Tuning uses a scoring-only prediction floor
  `mu_floor = 1e-8` so that sparse multichannel Poisson fits with zero predictions
  can be compared. `mu_floor = 0` requests exact likelihoods and can give infinite
  criteria; non-finite values are identified in the plot.

* `L0glm.fit(dispersion = NULL)` now uses a permissive weak-ridge pilot and
  bounded Pearson dispersion updates in the existing IRLS loop. The updates
  require adequate residual degrees of freedom and freeze after three changes
  below 1% on unchanged support. No nested fitting is required.
* Positive `dispersion` values remain fixed; negative values initialize at their
  absolute value and are subsequently estimated. Channel-specific mixtures are
  supported through grouped and parallel fits. Relaxed refits report a fresh
  residual scale for estimated channels. Poisson, binomial and known-theta
  negative-binomial dispersion is fixed at 1 on the general path too.
* Fixed Gaussian coefficient initialization with zero observation weights and
  its weighted rescaling; the automatic scale pilot preserves the caller's RNG.

* `L0glm.fit(openmp = TRUE)` now dispatches nonnegative identity-link Poisson
  fits using Eigen Cholesky LLT (`method = 1L`) to `L0nnpois`. Use
  `openmp = FALSE` to select the general R fitting path. Incompatible solver
  controls also retain the R path. `parallelize` and `ncores` control
  outcome-level concurrency in the compiled backend.
  Local Windows PSOCK workers use a clean R session and IPv4 loopback for their
  master connection.
* The default `L0nnpois` penalty is now `lambda = "gic"`.
* Added the Matthews correlation coefficient (`MCC`, or Pearson's phi) for
  support recovery to fit metrics and the README benchmark. Constant support
  vectors give `NA`, because their correlation is undefined.
* Fixed logical configuration flags in the README benchmark and added mirrored
  signal panels, full coefficient scatter plots, and support panels with dark
  blue true positives, red false positives and magenta false negatives.
* Audited package and benchmark solver calls and migrated remaining historical
  Eigen method numbers, including Cholesky LLT from `7L` to `1L`.
* Information-criterion names are now consistently case-insensitive in both
  fitters, tuning and benchmark helpers.
* Added a shared information-criteria help topic with penalty formulas, usage
  guidance, original references and the GIC/SIC/HBIC relationship.
* Fixed metric calculation when true coefficients are not supplied, allowing
  information-criterion tuning on observed data without simulation truth.
* Removed `lower` and `upper` from `L0glm.fit()` and `L0nnpois()`. General box
  constraints are no longer part of either fitting interface.
* `L0glm.fit(nonnegative = TRUE)` constrains all coefficients, including any
  intercept and `always.include` terms. Binomial and quasibinomial intercepts
  remain unrestricted; columns of ones supplied in `X` are recognized too.
  The default `nonnegative = FALSE` allows coefficients of either sign.
* Removed `nonnegative` from `L0nnpois()`: its coefficients, including any
  intercept, are now always nonnegative. Remove this argument from existing
  calls, and use `L0glm.fit(nonnegative = FALSE)` for unrestricted fits.
* Applied the same constraint rules to initialization, grouped and relaxed
  fits, and serial and parallel execution. The NNLS ridge backend now handles
  unrestricted coefficients, including the binomial intercept.
* Updated help, examples and benchmarks to use the simplified interfaces.
