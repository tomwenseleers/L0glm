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
L0glm uses a direct `L0glm.fit(..., lambda="gic", relaxed=TRUE,
parallelize=TRUE, ncores=4L)` call. Its time is the median of three runs
returning matching coefficients. GIC-tuned LASSO and fully relaxed LASSO
search ten lambda values, scoring each using the original Poisson likelihood.
The five-fold CV comparisons use native `cv.glmnet`, the same ten-value lambda
grid and common random folds, with four parallel workers. CV minimizes weighted
Gaussian prediction error. Every selected model is also scored using Poisson GIC.
CV timings include worker startup and the full-data fit. Simulation and final
metric reporting are excluded. The multichannel glmnet approximation uses
shared row weights `1/(rowMeans(y)+0.1)` and `family="mgaussian"`.
Fully relaxed models use `gamma=0`.

`benchmark_website_deconvolution.R` runs these five workflows, starting from
the same simulated data. Per-method checkpoints allow interrupted runs to
resume. The `*-preset-timings.csv` files retain every preset timing, and
`*-folds.csv` records the five-fold partition, and `*-tuning.csv` records all
candidate scores and the selected path values. The mirrored single-channel
signal panels show the first 1,000 time points; coefficient comparisons retain
all 10,000 predictors and all 226 multichannel outcomes.

Hardware: HP ZBook Fury 17.3-inch G8 Mobile Workstation, 11th-generation Intel
Core i9-11950H at 2.60 GHz, 8 physical cores and 16 logical cores. Benchmarks use
Windows and R 4.5.0 with OpenBLAS and one BLAS thread.

The million-predictor LASSO and relaxed-LASSO comparison is reproduced by
`README.Rmd`, including common random five-fold CV folds. Those CV results use
the package's common-scoring fold loop. Standalone `cv.glmnet` calls use
`benchmarks/cv_glmnet_parallel.R`, which manages a parallel backend and requires
`doParallel`. It limits simultaneous copies of large dense designs by default.

Run `Rscript benchmarks/check_group_irrepresentability.R` after the multichannel
fit to reproduce the standard weighted group-LASSO design diagnostic.
