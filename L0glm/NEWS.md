# L0glm 0.2.1

## Development changes

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
* Added `lambda = "mebic"` to `L0glm.fit()` and `L0nnpois()`, using
  `2 * log(p) * log(log(n))`, twice the GIC penalty.
* Information-criterion names are now consistently case-insensitive in both
  fitters, tuning and benchmark helpers. MEBIC is listed immediately after GIC.
* Added MEBIC to the fit metrics used for tuning. Added a shared
  information-criteria help topic with penalty formulas, usage guidance,
  original references and the GIC/SIC/HBIC relationship.
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
