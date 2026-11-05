# SIS screening and complete-model comparisons

HP ZBook Fury G8 Mobile Workstation, Intel Core i9-11950H, eight physical cores / 16 logical cores; R 4.5.0, SIS 1.5, one BLAS thread.

500 observations, Gaussian identity-link response, 60 nonzero coefficients
and `nsis=100` for one-pass SIS. Three repetitions; data generation is excluded.
The million-predictor data and default L0glm call match the homepage exactly.
L0glm and glmnet omit an intercept; SIS::SIS fits its native free intercept.

## Screening only

| Screening routine | Seconds |
|:--|--:|
| `sis_glm()`, eight physical cores | **0.24** |
| `sis_glm(nthreads=1)` | 0.74 |
| SIS Gaussian screen (`obtain.ix0`) | 2.00 |

All complete million-predictor rankings were identical in every repetition.
The SIS 1.5 internal `obtain.ix0` helper is called without fitting LASSO.
Standardization is disabled because it leaves Gaussian correlation ranks
unchanged. `timings.csv` contains all repetitions.

## Complete selection pipelines

| Predictors | Screen | Complete model | Seconds | GIC | TP | FP | F1 | MCC |
|--:|:--|:--|--:|--:|--:|--:|--:|--:|
| 1,000,000 | SIS | SIS::SIS + BIC LASSO | 7.31 | 4407 | 12 | 86 | 0.152 | 0.156 |
| 1,000,000 | SIS | L0glm.fit, GIC preset | 0.30 | 2793 | 7 | 1 | 0.206 | 0.319 |
| 1,000,000 | SIS | sparseglm LASSO, BIC grid | 1.57 | 4414 | 12 | 87 | 0.151 | 0.156 |
| 1,000,000 | SIS | sparseglm relaxed LASSO, BIC grid | 1.56 | 4240 | 12 | 80 | 0.158 | 0.161 |
| 1,000,000 | ISIS | L0glm.fit, GIC preset | 2.36 | 2583 | 40 | 3 | 0.777 | 0.787 |
| 1,000,000 | ISIS | sparseglm LASSO, BIC grid | 4.65 | 4402 | 9 | 191 | 0.069 | 0.082 |
| 1,000,000 | ISIS | sparseglm relaxed LASSO, BIC grid | 4.70 | 4187 | 9 | 185 | 0.071 | 0.083 |
| 5,000 | ISIS | SIS::SIS + BIC LASSO | 23.04 | 2883 | 42 | 57 | 0.528 | 0.538 |
| 5,000 | ISIS | L0glm.fit, GIC preset | 0.07 | 2180 | 46 | 0 | 0.868 | 0.874 |
| 5,000 | ISIS | sparseglm LASSO, BIC grid | 0.24 | 4581 | 47 | 228 | 0.281 | 0.352 |
| 5,000 | ISIS | sparseglm relaxed LASSO, BIC grid | 0.27 | 2939 | 47 | 71 | 0.528 | 0.552 |

L0glm uses its default GIC preset; LASSO variants use BIC-selected paths because
GIC commonly chooses the null LASSO model in this example. `sparseglm`
uses ten lambda values and includes common path scoring and diagnostics.
SIS::SIS uses its native LASSO path. Million-predictor ISIS uses the default
GIC-based batch and ten rounds. The public SIS-package ISIS comparison uses
5,000 predictors and three rounds; `nsis=100` is an added-candidate batch
size for L0glm and a total support cap for SIS. `models.csv` contains all
repetitions and GIC/BIC scores as well as recovery metrics.
Standalone L0glm.fit selects its best ISIS round by GIC; the LASSO adapters use
BIC for their outer searches. The default million-predictor fit reproduces GIC
2582.75, TP=40, FP=3, FN=20, TN=999937, F1=0.776699 and MCC=0.787489.
The recorded homepage timing is 3.49 s; the current median is 2.36 s.

### Compare all supported adapters

The same homepage data were also fitted with all nine adapters. This table
uses ISIS, the GIC preset for L0glm and BIC-selected paths for the competitors.
Times are one run of each method's screening, fitting and common scoring;
shared input preparation is separate.

| Method | Seconds | GIC | TP | FP | F1 | MCC |
|:--|--:|--:|--:|--:|--:|--:|
| L0glm | 2.20 | 2583 | 40 | 3 | 0.777 | 0.787 |
| abess | 11.32 | 3268 | 9 | 91 | 0.113 | 0.116 |
| L0Learn | 1.50 | 2790 | 11 | 2 | 0.301 | 0.394 |
| MCP | 3.50 | 4040 | 9 | 168 | 0.076 | 0.087 |
| SCAD | 3.29 | 4353 | 9 | 187 | 0.070 | 0.083 |
| lasso | 3.27 | 4402 | 9 | 191 | 0.069 | 0.082 |
| relaxed_lasso | 3.39 | 4187 | 9 | 185 | 0.071 | 0.083 |
| adaptive_lasso | 3.28 | 4316 | 9 | 184 | 0.071 | 0.084 |
| adaptive_relaxed_lasso | 3.39 | 4305 | 9 | 185 | 0.071 | 0.083 |

The [complete SIS and ISIS results](https://github.com/tomwenseleers/L0glm/blob/master/benchmarks/3_benchmarks_sis_vs_SIS_results/all_methods.csv)
and [reproducible calls](https://github.com/tomwenseleers/L0glm/blob/master/benchmarks/3_benchmarks_sparseglm_screening.R)
include every method and its tuning settings.



Run `benchmarks/3_benchmarks_sis_vs_SIS.R` from the checkout to reproduce both
comparisons. A one-pass marginal screen retains only 12 of the 60 true
predictors in the million-predictor simulation. SIS and ISIS use that same
homepage data; a separate 5,000-predictor case compares public SIS-package ISIS.
