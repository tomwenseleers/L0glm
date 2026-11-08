# L0glm


<div class="l0-hero">
<div class="l0-hero-copy">
<p class="l0-eyebrow">SPARSE REGRESSION · GENERALIZED LINEAR MODELS</p>
<h1>Fast variable selection.<br><span>Parsimonious models.</span></h1>
<p class="l0-lead">Fast L<sub>0</sub> and grouped L<sub>0</sub> regression to approximate best subset and best subset of groups selection. Accurate effects, cleaner supports, and one interface from ordinary GLMs to negative binomial, multinomial and survival models.</p>
<div class="l0-actions"><a class="btn btn-primary" href="articles/L0glm.html">Get started <span aria-hidden="true">↗</span></a><a class="btn btn-outline-light" href="articles/performance.html">Explore the benchmarks</a></div>
<p class="l0-hero-note">R · Dense & sparse designs · Single & multiple responses</p>
</div>
<img class="l0-hero-logo" src="logo-dark.svg" alt="L0glm: selected coefficient rows across response channels inside a navy hexagon" width="240" height="280">
</div>

<div class="l0-stat-grid">
<a href="articles/high-dimensional.html"><strong>1 million</strong><span>predictors in a high dimensional example</span></a>
<!-- BEGIN speed-stat -->
<a href="articles/performance.html#one-million-predictor-example"><strong>42× faster</strong><span>than five-fold-CV LASSO in that example</span></a>
<!-- END speed-stat -->
<!-- BEGIN multi-stat -->
<a href="articles/performance.html#multichannel-deconvolution"><strong>944 / 1000</strong><span>peaks recovered in simulated GC/MS deconvolution and <b>100× faster</b> than GIC-tuned relaxed group LASSO</span></a>
<!-- END multi-stat -->
</div>


L0glm selects predictors and estimates their effects in a generalized linear
model. For one response it selects individual coefficients. For multiple
responses, grouped selection identifies a common set of predictor rows and
estimates a separate coefficient for each outcome.


## Sparse learning across a wide range of models

L0glm combines adaptive ridge with iteratively reweighted least squares, fast
linear solvers and an optional relaxed refit. Small effects drop out; strong
effects retain their magnitude. A shrinking active set concentrates computation
where it matters, producing sparse models with excellent support recovery.

<div class="l0-feature-grid">
<div><span class="l0-feature-number">01 / MODELS</span><h3>One fitting interface</h3><p>Gaussian, Poisson, binomial, Gamma and inverse Gaussian GLMs; negative binomial with estimated shape; joint multinomial regression; and Cox proportional hazards.</p></div>
<div><span class="l0-feature-number">02 / SCALE</span><h3>Built for sparse signals</h3><p>Dense or sparse matrices, tall or wide systems, grouped multichannel selection, and compiled SIS/ISIS screening for ultrahigh dimensions.</p></div>
<div><span class="l0-feature-number">03 / CONTROL</span><h3>Fit the model you need</h3><p>Observation weights, offsets, nonnegative coefficients, unpenalized terms, normalization and original-scale estimates. Native OpenMP accelerates nonnegative identity-link Poisson fits.</p></div>
<div><span class="l0-feature-number">04 / SELECTION</span><h3>Choose your workflow</h3><p>Start with a fast GIC preset, tune an information criterion, or use cross-validation with minimum-deviance and one-standard-error models. Compare methods with <code>sparseGLM()</code>.</p></div>
</div>

## Fit a sparse model

<!-- BEGIN million-intro -->
This Gaussian simulation has **one million predictors**, 500 observations and 60 nonzero effects. One default L0glm call takes **2.84 seconds** and recovers **40 true effects with 2 false positives**.
<!-- END million-intro -->

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

system.time(fit <- l0glmFit(X, y, family=gaussian()))
l0glmMetrics(fit, beta_true=beta)[
  c("gic", "TP", "FP", "FN", "TN", "F1_score", "MCC")]
```

<!-- BEGIN million-metrics -->
| Seconds | GIC | TP | FP | FN | TN | F1 | MCC |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 2.84 | 2579 | 40 | 2 | 20 | 999938 | 0.784 | 0.797 |
<!-- END million-metrics -->

The recovered coefficients are close to their true values.

<img src="million-coefficients.png" alt="Recovered versus true effects in the million-predictor example." width="100%">

Install with `remotes::install_github("tomwenseleers/L0glm", subdir="L0glm")`.
The dense design occupies about 4 GB. The [getting-started guide](articles/L0glm.html)
explains the example and prediction methods.

## Application: sparse deconvolution

In **GC/MS analysis**, overlapping chromatographic peaks have distinct mass
spectra. Recovering their peak times and spectra together is a natural application
of grouped L₀ regression. This Poisson simulation has **10,000 observations,
10,000 candidate peak locations, 1,000 true peaks and 226 mass channels**.

By regressing the observed signal on a design matrix consisting of shifted
copies of the peak-shape function (the point spread function), L0glm recovers
the underlying peaks and their mass spectra.

<!-- BEGIN multi-speed -->
L0glm recovers **944 of 1,000 true peaks with 11 false positives** in **8.03 seconds**. In this benchmark it is **25.0× faster than GIC-tuned group LASSO** and **306× faster than relaxed group LASSO tuned by five-fold CV**.
<!-- END multi-speed -->
<!-- BEGIN multi-spectral -->
Among correctly located peaks, the recovered mass spectra have a within-row F1 of **0.912** and MCC of **0.869**.
<!-- END multi-spectral -->
A direct `l0glmFit()` call uses `lambda="gic"` and four cores.
The competing methods select lambda from a ten-value grid using GIC or
five-fold cross-validation.

```r
multi <- simul(method=4L, n=10000L, p=10000L,
  n_nonzero=1000L, seed=2L)
fit_multi <- l0glmFit(multi$X, multi$y, family=poisson("identity"),
  lambda="gic", nonnegative=TRUE, ncores=4L)
```

<!-- BEGIN multi-home -->
| Method | Selection | Seconds | GIC | TP | FP | F1 | MCC |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| L0glm | GIC preset | 8.03 | 5606331 | 944 | 11 | 0.966 | 0.962 |
| Group LASSO | GIC path | 201.03 | 7413616 | 980 | 860 | 0.690 | 0.685 |
| Relaxed group LASSO | GIC path | 804.89 | 6138687 | 909 | 309 | 0.820 | 0.802 |
| Group LASSO | 5-fold CV | 593.23 | 9639345 | 998 | 4242 | 0.320 | 0.316 |
| Relaxed group LASSO | 5-fold CV | 2460.38 | 6198648 | 962 | 502 | 0.781 | 0.769 |
<!-- END multi-home -->

![Multichannel GC/MS peak recovery: observed data, ground truth, L0glm and relaxed group LASSO, retaining all 226 channels. Cyan lines mark true positives, red false positives and magenta false negatives.](deconvolution-multi-large.png)

The coefficient panels compare the individual mass-spectral intensities.
They show L0glm and the alternative with the highest peak-level MCC.

![True versus estimated coefficients across all channels for L0glm and relaxed group LASSO.](deconvolution-multi-home-coefficients.png)

L0glm fits the Poisson model directly. The group-LASSO methods use a
weighted-Gaussian approximation with shared observation weights; CV uses four
parallel workers. [Explore the deconvolution examples and full
benchmark tables →](articles/deconvolution.html)

When peaks have different shapes, logconNMF can estimate a response dictionary
and then fit a log-concave profile for each component. GCMSprime handles
chromatogram files, quantification and spectral library search.
[logconNMF deconvolution](https://tomwenseleers.github.io/logconNMF/articles/log-concave-nmf.html) ·
[GCMSprime workflow](https://tomwenseleers.github.io/GCMSprime/articles/gcms-deconvolution.html)

## Choose a starting point

| I would like to… | Start here |
|---|---|
| Fit a model and make predictions | [Get started](articles/L0glm.html) |
| Choose a penalty by GIC or cross-validation | [Model selection](articles/selection.html) |
| Use weights, constraints or held-out early stopping | [Fitting controls](articles/controls.html) |
| Fit standard GLMs | [Regular GLM models](articles/regular-glm.html) |
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
*L0glm: Scalable L₀ and Grouped L₀-Penalized Variable Selection for Generalized Linear Models*.
R package version 0.2.0.9000,
<https://tomwenseleers.github.io/L0glm/>.

The [Authors and Citation page](authors.html#citation) also provides the
citation in BibTeX format.

## Acknowledgements

OpenAI Codex, using GPT-6 Astra and GPT-6.1 Sol, was used to optimize and
refactor code and help draft the documentation. The key ideas and methodological
decisions are the authors' own, and the authors take full responsibility for the
work presented here.
