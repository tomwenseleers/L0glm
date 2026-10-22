These measured examples supply the pkgdown homepage and tutorials.

From the repository root, run:

```r
source("benchmarks/benchmark_website_examples.R")
website_fit_examples()
source("benchmarks/plot_website_examples.R")
website_plot_examples()
```

Then run `python benchmarks/update_website_tables.py` and
`Rscript L0glm/pkgdown/build.R` to publish the tables and build the site.
The fitting function also accepts a subset of cases, for example
`website_fit_examples("million")`. Cached deconvolution fits are reused;
remove the corresponding RDS files from `codex_sandbox/website_large_examples`
to fit those methods again. Cached model objects remain outside Git.

The million-predictor design is dense and occupies about 4 GB. Its default-fit
time is the median of three runs; all three return identical coefficients.
Screening uses the physical cores (eight on this workstation), with one BLAS
thread. Deconvolution has 10,000 observations, 10,000 predictors and 1,000 true
peaks; the multichannel simulation retains all 226 informative mass channels.
All methods use the same Poisson GIC for selection and scoring. The L0glm row
uses one preset fit with a relaxed refit; glmnet uses ten lambda values.

The million-predictor LASSO and relaxed-LASSO comparison is reproduced by
`README.Rmd`, including common random five-fold CV folds. Those CV results use
the package's common-scoring fold loop. Standalone `cv.glmnet` calls use
`benchmarks/cv_glmnet_parallel.R`, which manages a parallel backend and requires
`doParallel`. It limits simultaneous copies of large dense designs by default.

Run `Rscript benchmarks/check_group_irrepresentability.R` after the multichannel
fit to reproduce the standard weighted group-LASSO design diagnostic.
