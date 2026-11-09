These measured examples supply the README and pkgdown tutorials. Install the
checkout of L0glm you want to measure, and run the scripts with that installation
first on the R library path.

The million-predictor runs use L0glm 0.2.0.9002 at fitter revision
[`07af781`](https://github.com/tomwenseleers/L0glm/commit/07af781).
See [`sessionInfo.txt`](sessionInfo.txt) for R and dependency versions.

From the repository root, run:

```r
source("benchmarks/benchmark_website_examples.R")
website_fit_examples()
source("benchmarks/plot_website_examples.R")
website_plot_examples()
```

Then run `Rscript benchmarks/render_readme.R` and
`Rscript L0glm/pkgdown/build.R` to publish the same tables in the README and website.
The README renderer runs `python benchmarks/update_website_tables.py` and reuses
the published figures; it does not rerun the large fitting comparisons.
The fitting function also accepts a subset of cases, for example
`website_fit_examples("million")`. Deconvolution fits are recomputed by default.
Use `reuse=TRUE` only to resume checkpoints made with the same package versions,
simulation and fitting settings. Cached model objects remain outside Git.

The million-predictor design has 1,000 observations and 60 signed nonzero
effects. It is dense and occupies about 8 GB. Its default-fit
time is the median of three runs; all three return identical coefficients.
Screening uses the physical cores (eight on this workstation), with one BLAS
thread. Deconvolution has 10,000 observations, 10,000 predictors and 1,000 true
peaks; the multichannel simulation retains all 226 informative mass channels.
L0glm uses a direct `l0glmFit(..., lambda="gic", relaxed=TRUE,
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
`Rscript benchmarks/benchmark_million_lasso.R`, including common random five-fold
CV folds. Its published results are stored in `million-lasso-metrics.csv`; the
GIC-preset L0glm result comes from `million-metrics.csv`. The LASSO paths set
`screen="none"` and use the package's common-scoring fold loop for CV.
The shared-SIS and shared-ISIS comparisons are available separately in
[`3_benchmarks_sis_vs_SIS_results`](../../3_benchmarks_sis_vs_SIS_results).
Each million-predictor LASSO workflow runs in a fresh R worker, with one BLAS
thread and serial CV folds. Input loading and worker startup are excluded
from its fitting/selection time. This keeps native working allocations bounded;
the data, 40-value grids and shared folds are identical across workflows.
Standalone `cv.glmnet` calls use
`benchmarks/cv_glmnet_parallel.R`, which manages a parallel backend and requires
`doParallel`. It limits simultaneous copies of large dense designs by default.

Run `Rscript benchmarks/check_group_irrepresentability.R` after the multichannel
fit to reproduce the standard weighted group-LASSO design diagnostic.
