# SIS screening and complete-model comparisons

The Gaussian simulation has 1,000 observations, one million standard-normal predictors and 60 signed nonzero effects, matching the homepage exactly. Seeds are 123 for the design/effects and 1 for the noise. All timing medians use three repetitions and exclude simulation. Hardware: HP ZBook Fury G8, Intel i9-11950H, eight physical cores, 16 logical processors, Windows, R 4.5.0 and one BLAS thread.

## Screening only

| Method | Seconds |
| --- | ---: |
| `sisGLM()`, eight physical cores | 0.53 |
| `sisGLM(nthreads=1)` | 1.81 |
| SIS Gaussian screen (`obtain.ix0`) | 4.64 |

All complete predictor rankings agree. SIS's `obtain.ix0` helper excludes its subsequent LASSO fit; Gaussian rankings are invariant to column standardization.

## Complete selection pipelines

| Predictors | Screen | Complete model | Seconds | GIC | TP | FP | F1 | MCC |
| ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1000000 | SIS | SIS::SIS + BIC LASSO | 9.05 | 6351 | 25 | 73 | 0.316 | 0.326 |
| 1000000 | SIS | l0glmFit, GIC preset | 0.63 | 4660 | 25 | 0 | 0.588 | 0.645 |
| 1000000 | SIS | sparseGLM LASSO, BIC grid | 3.44 | 6325 | 25 | 73 | 0.316 | 0.326 |
| 1000000 | SIS | sparseGLM relaxed LASSO, BIC grid | 3.21 | 4952 | 25 | 16 | 0.495 | 0.504 |
| 1000000 | ISIS | l0glmFit, GIC preset | 2.72 | 4198 | 46 | 1 | 0.860 | 0.866 |
| 1000000 | ISIS | sparseGLM LASSO, BIC grid | 7.50 | 7747 | 41 | 329 | 0.191 | 0.275 |
| 1000000 | ISIS | sparseGLM relaxed LASSO, BIC grid | 7.74 | 7504 | 46 | 312 | 0.220 | 0.314 |
| 5000 | ISIS | SIS::SIS + BIC LASSO | 42.97 | 4388 | 47 | 52 | 0.591 | 0.604 |
| 5000 | ISIS | l0glmFit, GIC preset | 0.08 | 3714 | 45 | 0 | 0.857 | 0.865 |
| 5000 | ISIS | sparseGLM LASSO, BIC grid | 0.25 | 5217 | 50 | 121 | 0.433 | 0.485 |
| 5000 | ISIS | sparseGLM relaxed LASSO, BIC grid | 0.29 | 3931 | 46 | 23 | 0.713 | 0.711 |

One-pass SIS retains 100 predictors. Million-predictor ISIS uses the default GIC batch and ten rounds; the SIS-package ISIS comparison uses 5,000 predictors and three rounds. `nsis=100` is a batch size for L0glm and a total support cap for SIS. These comparisons set `local.search=FALSE` for standalone L0glm, isolating screening and adaptive ridge. SIS fits a free intercept; L0glm and glmnet use the simulation's zero intercept.

## All supported adapters with ISIS

| Method | Seconds | GIC | TP | FP | F1 | MCC |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| L0glm | 2.98 | 4194 | 46 | 0 | 0.868 | 0.876 |
| ABESS | 16.39 | 4575 | 43 | 57 | 0.537 | 0.555 |
| L0Learn | 1.48 | 5050 | 13 | 0 | 0.356 | 0.465 |
| MCP | 6.08 | 6859 | 45 | 291 | 0.227 | 0.317 |
| SCAD | 6.21 | 7491 | 47 | 309 | 0.226 | 0.322 |
| LASSO | 5.81 | 7747 | 41 | 329 | 0.191 | 0.275 |
| Relaxed LASSO | 5.83 | 7504 | 46 | 312 | 0.220 | 0.314 |
| Adaptive LASSO | 5.93 | 7373 | 44 | 297 | 0.219 | 0.308 |
| Adaptive relaxed LASSO | 6.16 | 6946 | 43 | 294 | 0.217 | 0.302 |

L0glm uses the GIC preset; competitors select BIC paths on their own screening pools. The adapter table reports one complete fitting/scoring run per method. The CSV files retain every measured run, criterion and setting. Reproduce with [3_benchmarks_sis_vs_SIS.R](../3_benchmarks_sis_vs_SIS.R) and [3_benchmarks_sparseglm_screening.R](../3_benchmarks_sparseglm_screening.R).
