# L0glm 0.1.0

## Initial version

* Sparse Poisson log-concave NMF (`spoloconmf`) recovers elution profiles and
  spectra using local Poisson updates and fast spectral fusion, with a
  compound-count criterion, unpenalized spline backgrounds and native OpenMP.

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
* Small constrained solves for nonnegative compulsory terms and smooth spline
  backgrounds, alternating with sparse selection on the original channels.
  Gaussian, Poisson, other scalar GLMs and negative-binomial/Cox fits retain
  their distribution-specific likelihoods, dispersion and validation controls.
* Fixed-PSF and blind deconvolution use positive signal-based starts with
  backgrounds, channel-specific PSF offsets, and synchronized final estimates.
  Automatic identity-link initialization accounts for nonnegative compulsory
  backgrounds and known offsets while preserving original coefficient units.
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

* Native scalar-family deviance and Bregman clustering use explicit pooled
  likelihood merge costs, IC-scored tree cuts and dispersion-aware precision.
  Poisson aggregation preserves shared case weights.

* Optional Poisson channel preclustering screens aggregated responses and
  refines on the original channels; aggregate-only fits retain the requested
  penalty. Starts, weights, offsets and original predictor dimensions are preserved.
* Preclustering-aware lambda paths reuse training-derived memberships while
  estimating every model independently; cross-validation learns memberships and
  size factors within each training fold.
