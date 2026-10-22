# L0glm

<div class="l0-hero">
<div class="l0-hero-copy">
<p class="l0-eyebrow">SPARSE REGRESSION · GENERALIZED LINEAR MODELS</p>
<h1>Fast variable selection.<br><span>Parsimonious models.</span></h1>
<p class="l0-lead">Fast L<sub>0</sub> and grouped L<sub>0</sub> regression. Accurate effects, cleaner supports, and one interface from ordinary GLMs to multinomial and survival models.</p>
<div class="l0-actions"><a class="btn btn-primary" href="articles/L0glm.html">Get started <span aria-hidden="true">↗</span></a><a class="btn btn-outline-light" href="articles/performance.html">Explore the benchmarks</a></div>
<p class="l0-hero-note">R · Dense & sparse designs · Single & multiple responses</p>
</div>
<img class="l0-hero-logo" src="logo-dark.svg" alt="L0glm: selected coefficient rows across response channels inside a navy hexagon" width="240" height="280">
</div>

<div class="l0-stat-grid">
<a href="articles/high-dimensional.html"><strong>1 million</strong><span>predictors in a high dimensional example</span></a>
<!-- BEGIN speed-stat -->
<a href="articles/performance.html#one-million-predictors"><strong>35× faster</strong><span>than five-fold-CV LASSO in that example</span></a>
<!-- END speed-stat -->
<!-- BEGIN multi-stat -->
<a href="articles/performance.html#multichannel-deconvolution"><strong>944 / 1000</strong><span>true peaks recovered in multichannel GC/MS deconvolution</span></a>
<!-- END multi-stat -->
</div>

## Sparse learning across a wide range of models

L0glm combines adaptive ridge with iteratively reweighted least squares, fast
linear solvers and an optional relaxed refit. Small effects drop out; strong
effects retain their magnitude. A shrinking active set concentrates computation
where it matters, producing sparse models with excellent support recovery.

<div class="l0-feature-grid">
<div><span class="l0-feature-number">01 / MODELS</span><h3>One fitting interface</h3><p>Gaussian, Poisson, binomial, Gamma and inverse Gaussian GLMs; negative binomial with estimated shape; joint multinomial regression; and Cox proportional hazards.</p></div>
<div><span class="l0-feature-number">02 / SCALE</span><h3>Built for sparse signals</h3><p>Dense or sparse matrices, tall or wide systems, grouped multichannel selection, and compiled SIS/ISIS screening for ultrahigh dimensions.</p></div>
<div><span class="l0-feature-number">03 / CONTROL</span><h3>Fit the model you need</h3><p>Observation weights, offsets, nonnegative coefficients, unpenalized terms, normalization and original-scale estimates. Native OpenMP accelerates nonnegative identity-link Poisson fits.</p></div>
<div><span class="l0-feature-number">04 / SELECTION</span><h3>Choose your workflow</h3><p>Start with a fast GIC preset, tune an information criterion, or use cross-validation with minimum-deviance and one-standard-error models. Compare methods with <code>sparseglm()</code>.</p></div>
</div>

## Fit a sparse model

Sixty nonzero effects among **one million predictors**, fitted in one call.

```r
library(L0glm)
set.seed(123)
n <- 500L
p <- 1000000L
X <- matrix(rnorm(n*p), n, p)
active <- sample.int(p, 60L)
beta <- numeric(p)
beta[active] <- runif(60L, -.8, .8)
set.seed(1)
y <- as.numeric(X[, active] %*% beta[active]) + rnorm(n)

system.time(fit <- L0glm.fit(X, y, family=gaussian()))
L0glm.metrics(fit, beta_true=beta)[
  c("gic", "TP", "FP", "FN", "TN", "F1_score", "MCC")]
```

<!-- BEGIN million-metrics -->
| Seconds | GIC | TP | FP | FN | TN | F1 | MCC |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 3.49 | 2583 | 40 | 3 | 20 | 999937 | 0.777 | 0.787 |
<!-- END million-metrics -->

<!-- BEGIN million-confusion -->
| Ground truth | Selected | Excluded |
|---|---:|---:|
| Nonzero effects | 40 TP | 20 FN |
| Zero effects | 3 FP | 999937 TN |
<!-- END million-confusion -->

![True versus estimated coefficients in the million-predictor default fit.](million-coefficients.png)

Install with `remotes::install_github("tomwenseleers/L0glm", subdir="L0glm")`.
The dense design occupies about 4 GB. The [getting-started guide](articles/L0glm.html)
explains the example and prediction methods.

### Tuning when you need it

The preset often already gives an excellent sparse model. Detailed tuning over
a lambda grid is also supported; in this example it selects exactly the preset.

```r
tuned <- L0glm.tune(X, y, family=gaussian(), tune="gic",
  nlambda=10L, prop_range=.25, coefs.true=beta, plot=TRUE)
```

![Coefficient paths and GIC in the million-predictor example, with the selected and preset lambda coinciding.](million-path.png)

Use `tune="cv", nfolds=5L` for five-fold validation. See
[model selection](articles/selection.html) for a 100-signal example and the
one-standard-error rule.

## Example: sparse deconvolution

In **GC/MS analysis**, overlapping chromatographic peaks have distinct mass
spectra. Recovering their peak times and spectra together is a natural application
of grouped L₀ regression. This Poisson simulation has **10,000 observations,
10,000 candidate peak locations, 1,000 true peaks and 226 mass channels**.

L0glm recovers **944 of 1,000 true peaks with just 11 false positives**. It also
recovers the sparse mass spectra accurately, with a within-row F1 of **0.909**.
The comparison uses the same Poisson GIC to select and score every method.

<!-- BEGIN multi-home -->
| Method | Seconds | TP | FP | F1 | MCC |
| --- | ---: | ---: | ---: | ---: | ---: |
| L0glm | 34.14 | 944 | 11 | 0.966 | 0.962 |
| Group LASSO | 244.87 | 980 | 860 | 0.690 | 0.685 |
| Relaxed group LASSO | 906.34 | 909 | 309 | 0.820 | 0.802 |
<!-- END multi-home -->

![Multichannel GC/MS peak recovery: observed data, ground truth, L0glm and relaxed group LASSO, retaining all 226 channels. Cyan lines mark true positives, red false positives and magenta false negatives.](deconvolution-multi-large.png)

The coefficient comparison shows how closely L0glm recovers the individual
mass-spectral intensities. Relaxed group LASSO is the strongest alternative in
this benchmark.

![True versus estimated coefficients across all channels for L0glm and relaxed group LASSO.](deconvolution-multi-home-coefficients.png)

L0glm fits the Poisson model directly. The group-LASSO methods use a
weighted-Gaussian approximation with shared observation weights; their timings
include a ten-value tuning path. [Explore the deconvolution examples and full
benchmark tables →](articles/deconvolution.html)

## Choose a starting point

| I would like to… | Start here |
|---|---|
| Fit a model and make predictions | [Get started](articles/L0glm.html) |
| Choose a penalty by GIC or cross-validation | [Model selection](articles/selection.html) |
| Use weights, constraints or held-out early stopping | [Fitting controls](articles/controls.html) |
| Fit counts, classes or survival outcomes | [Negative binomial](articles/negative-binomial.html) · [Multinomial](articles/multinomial.html) · [Cox](articles/cox.html) |
| Search a very large predictor set | [SIS / ISIS](articles/high-dimensional.html) |
| Understand the mathematics | [Algorithm & theory](articles/algorithm.html) |
| Compare speed and support recovery | [Benchmarks](articles/performance.html) |

## Foundations

The method builds on adaptive ridge and broken adaptive ridge, with an IRLS
implementation for generalized linear models. The [algorithm article](articles/algorithm.html)
develops the quadratic updates and their connection to sparse penalties, drawing
on [Frommlet & Nuel (2016)](https://doi.org/10.1371/journal.pone.0148620),
[Li et al. (2021)](https://doi.org/10.1016/j.jspi.2020.12.001) and
[Abergel, Bouaziz & Nuel (2024)](https://doi.org/10.1007/s11222-024-10440-6).
It also explains screening, grouped selection and likelihood-based scoring.
The adaptive updates also have a [hierarchical Bayesian and automatic relevance
determination (ARD) interpretation](articles/algorithm.html#a-hierarchical-bayesian-and-ard-view):
each ridge step estimates a Gaussian-prior MAP, then the prior scales adapt to
the fitted effects. Zero-centered skeptical priors tighten around small effects
and relax around strong ones.

## Citation

An article describing the methods implemented in L0glm is in preparation.
If you use the package before the article is formally published, we would
appreciate a citation to the package:

Wenseleers T, Lanzas L, Vanderaa C (2026).
*L0glm: Fast Sparse Regression for Generalized Linear Models*.
R package version 0.1.0,
<https://tomwenseleers.github.io/L0glm/>.

The [Authors and Citation page](authors.html#citation) also provides the
citation in BibTeX format.

## Acknowledgements

OpenAI Codex, using GPT-6 Astra and GPT-6.1 Sol, was used to optimize and
refactor code and help draft the documentation. The key ideas and methodological
decisions are the authors' own, and the authors take full responsibility for the
work presented here.
