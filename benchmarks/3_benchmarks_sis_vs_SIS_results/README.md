# One-million-predictor SIS comparison

Gaussian data with 500 observations, 1,000,000 standard-normal predictors,
60 nonzero effects and a retained screening set of 100 predictors. Run
`benchmarks/3_benchmarks_sis_vs_SIS.R` to reproduce the data and comparisons.

| Method | Median elapsed time, three repetitions |
|:--|--:|
| `sis_glm()`, default eight physical cores | 0.23 s |
| `sis_glm(nthreads=1)` | 0.82 s |
| `SIS::SIS(iter=FALSE, standardize=FALSE)` | 2.69 s |

All three methods returned exactly the same complete ranking of all one million
predictors, including the same ordered top 100. Spearman rank correlation was 1
in every repetition. The default native screen was approximately 11.7 times
faster than the configured public SIS call.

The Gaussian SIS initial screen ranks absolute Pearson correlations, so
`standardize=FALSE` preserves its ranking and avoids a redundant large design
copy. The public SIS call also fits a BIC-tuned LASSO on the 100 retained columns;
the comparisons use its initial `sis.ix0` and `ix_list[[1L]]` outputs. Data
generation, scoring agreement and garbage collection before each call are
outside the timed region.

Measured with R 4.5.0, SIS 1.5 and OpenBLAS using one BLAS thread, on an HP ZBook
Fury G8 with Intel Core i9-11950H, 8 physical/16 logical cores. `timings.csv`
contains the individual measurements.
