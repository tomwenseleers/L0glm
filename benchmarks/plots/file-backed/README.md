# Reproducible file-backed examples

Run the scripts from the repository root with the current package installed.
Large input files and model caches are written beneath the ignored
`codex_sandbox` directory; these CSV outputs retain the measured results.

The matched 80 GB comparison fits free-intercept models to 8,000 rows and
predicts 2,000 held-out rows:

<!-- BEGIN file-lasso-metrics -->
| Method | Selection | Seconds | GIC | TP | FP | FN | F1 | MCC | Test RMSE |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| L0glm | GIC preset | 248.72 | 24652 | 50 | 1 | 10 | 0.901 | 0.904 | 0.993 |
| biglasso LASSO | BIC path | 792.32 | 25933 | 51 | 3 | 9 | 0.895 | 0.896 | 1.074 |
<!-- END file-lasso-metrics -->

* `benchmark_filebacked_normal.R` generates a 10,000-observation by
  1,000,000-predictor Gaussian design (80 GB in double precision) in bounded
  256-column blocks, with 60 signed effects. It fits a free-intercept ISIS
  model with direct native screening and a 10,000-column resident pool on 8,000 training rows and
  evaluates prediction on 2,000 held-out rows. The history and metrics use
  the original predictor count.
* `benchmark_filebacked_lasso.R` reuses that same file for a free-intercept
  comparison: L0glm with a GIC preset and ISIS, versus native biglasso LASSO
  with BIC selection on a ten-value path down to 5% of its native lambda maximum.
  It uses 8,000 training and 2,000 test
  rows. Both fit a free intercept and read the double backing matrix directly.
  L0glm retains four global screening scans and the original predictor count
  for its GIC comparisons.
* `benchmark_filebacked_gwas.R` uses `helpers/hwe_genotypes.cpp` to create a
  portable PLINK BED/BIM/FAM dataset with 4,000 observations and one million
  independent HWE variants (1 GB packed). It simulates six main effects and
  one interaction, compares SIS, ISIS and ISIS with pairwise interactions,
  and reports recovery and prediction on 800 held-out observations.

```r
source("benchmarks/benchmark_filebacked_normal.R")
benchmark_filebacked_normal()
source("benchmarks/benchmark_filebacked_lasso.R")
benchmark_filebacked_lasso()
source("benchmarks/benchmark_filebacked_gwas.R")
benchmark_filebacked_gwas()
```

Generation and fitting times are separate. BLAS uses one thread and screening
uses eight physical cores on an HP ZBook Fury G8 with an Intel i9-11950H,
64 GB RAM, Windows and R 4.5.0. The retained fitting matrix fits in RAM;
screening streams the full source. Windows may cache file reads.
The LASSO comparison uses eight native biglasso threads and retains all
predictors, with no SIS/ISIS prescreen. It excludes file generation and final
recovery/test reporting, and includes the complete native path and common IC
scoring. Five-fold CV is supported by the same row-view interface; the recorded
large-file table uses BIC selection.
The measured BIC minimum is interior (ninth of ten candidates, lambda 0.056899).
`normal-10000-1000000-lasso-path.csv` retains that complete native path.

In the GWAS table TP/FP/FN/F1/MCC use a common evaluation universe containing
all original predictors, the union of generated products and the true product.
An interaction omitted by a main-effect model is a false negative. The column
`interaction_recovered` identifies the correct product explicitly. BED counts
BIM allele 1, and FAM phenotypes encode controls as 1 and cases as 2.
