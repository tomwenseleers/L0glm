# L0glm 0.2.1

## Development changes

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
  Local Windows PSOCK workers use `localhost` for their master connection.
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
