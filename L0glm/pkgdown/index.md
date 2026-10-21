# L0glm

<div class="l0-hero">
<div class="l0-hero-copy">
<p class="l0-eyebrow">SPARSE REGRESSION · GENERALIZED LINEAR MODELS</p>
<h1>Find the signal.<br><span>Keep the essentials.</span></h1>
<p class="l0-lead">Fast L<sub>0</sub> and grouped L<sub>0</sub> regression. Accurate effects, cleaner supports, and one interface from ordinary GLMs to multinomial and survival models.</p>
<div class="l0-actions"><a class="btn btn-primary" href="articles/L0glm.html">Get started <span aria-hidden="true">↗</span></a><a class="btn btn-outline-light" href="articles/performance.html">Explore the benchmarks</a></div>
<p class="l0-hero-note">R · Dense & sparse designs · Single & multiple responses</p>
</div>
<img class="l0-hero-logo" src="logo.png" alt="L0glm: a sparse signal and fitted curve inside a navy hexagon" width="240" height="240">
</div>

<div class="l0-stat-grid">
<a href="articles/high-dimensional.html"><strong>1 million</strong><span>predictors in the high-dimensional example</span></a>
<a href="articles/performance.html#one-million-predictors"><strong>39× faster</strong><span>than five-fold-CV LASSO in that example</span></a>
<a href="articles/performance.html#multichannel-deconvolution"><strong>25 / 25</strong><span>true rows recovered, with zero false positives</span></a>
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

## A sparse model in a few lines

```r
install.packages("remotes")
remotes::install_github("tomwenseleers/L0glm", subdir = "L0glm")

library(L0glm)
set.seed(7)
X <- matrix(rnorm(400 * 30), 400, 30)
y <- 2 * X[, 1] - 1.5 * X[, 8] + rnorm(400)

fit <- L0glm.fit(X, y, family = gaussian(),
  lambda = "gic", relaxed = TRUE)
summary(fit)
predict(fit, X[1:5, ])
```

On Windows, source installation uses [Rtools](https://cran.r-project.org/bin/windows/Rtools/).
The [getting-started guide](articles/L0glm.html) covers fitting, prediction and metrics.

## See the difference

Poisson spike-train deconvolution is one application: recover sharp events from
a blurred, noisy signal. The same fitting machinery also handles sparse
classification, overdispersed counts and survival outcomes.

![Observed and true signals above, mirrored fitted signals and estimated spikes below, with separate panels for L0glm, ABESS, LASSO and relaxed LASSO.](deconvolution-single.png)

[Explore single and multichannel deconvolution →](articles/deconvolution.html)

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
