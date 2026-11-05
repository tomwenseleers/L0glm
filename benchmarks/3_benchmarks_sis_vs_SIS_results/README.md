# SIS screening and complete-model comparisons

HP ZBook Fury G8 Mobile Workstation, Intel Core i9-11950H, eight physical cores / 16 logical cores; R 4.5.0, SIS 1.5, one BLAS thread.

500 observations, Gaussian identity-link response, 60 nonzero coefficients
and `nsis=100`. Three repetitions; data generation is excluded. All complete
models use an intercept and are scored against the same original coefficients.

## Screening only

| Screening routine | Seconds |
|:--|--:|
| `sis_glm()`, eight physical cores | **0.26** |
| `sis_glm(nthreads=1)` | 0.83 |
| SIS Gaussian screen (`obtain.ix0`) | 2.16 |

All complete million-predictor rankings were identical in every repetition.
The SIS 1.5 internal `obtain.ix0` helper is called without fitting LASSO.
Standardization is disabled because it leaves Gaussian correlation ranks
unchanged. `timings.csv` contains all repetitions.

## Complete selection pipelines

| Predictors | Screen | Complete model | Seconds | TP | FP | F1 | MCC |
|--:|:--|:--|--:|--:|--:|--:|--:|
| 1,000,000 | SIS | SIS::SIS + BIC LASSO | 2.71 | 12 | 87 | 0.151 | 0.156 |
| 1,000,000 | SIS | L0glm.fit, BIC preset | 1.43 | 12 | 46 | 0.203 | 0.203 |
| 1,000,000 | SIS | sparseglm LASSO, BIC grid | 2.58 | 12 | 87 | 0.151 | 0.156 |
| 1,000,000 | SIS | sparseglm relaxed LASSO, BIC grid | 2.62 | 12 | 76 | 0.162 | 0.165 |
| 1,000,000 | ISIS | L0glm.fit, BIC preset | 3.52 | 7 | 12 | 0.177 | 0.207 |
| 1,000,000 | ISIS | sparseglm LASSO, BIC grid | 5.98 | 10 | 190 | 0.077 | 0.091 |
| 1,000,000 | ISIS | sparseglm relaxed LASSO, BIC grid | 6.65 | 10 | 190 | 0.077 | 0.091 |
| 5,000 | ISIS | SIS::SIS + BIC LASSO | 23.04 | 58 | 39 | 0.739 | 0.757 |
| 5,000 | ISIS | L0glm.fit, BIC preset | 0.09 | 60 | 26 | 0.822 | 0.833 |
| 5,000 | ISIS | sparseglm LASSO, BIC grid | 0.33 | 60 | 182 | 0.397 | 0.489 |
| 5,000 | ISIS | sparseglm relaxed LASSO, BIC grid | 0.39 | 60 | 57 | 0.678 | 0.712 |

L0glm uses a BIC preset; LASSO variants use BIC-selected paths. `sparseglm`
uses ten lambda values and includes common path scoring and diagnostics.
SIS::SIS uses its native LASSO path. Million-predictor ISIS uses the default
GIC-based batch and ten rounds. The public SIS-package ISIS comparison uses
5,000 predictors and three rounds; `nsis=100` is an added-candidate batch
size for L0glm and a total support cap for SIS. `models.csv` contains all
repetitions and GIC/BIC scores as well as recovery metrics.
Standalone L0glm.fit selects its best ISIS round by GIC; sparseglm uses BIC
for its outer search in this benchmark.

Run `benchmarks/3_benchmarks_sis_vs_SIS.R` from the checkout to reproduce both
comparisons. A one-pass marginal screen retains only 12 of the 60 true
predictors in the million-predictor simulation, limiting subsequent recovery;
the ISIS example tests iterative recovery on its own generated design.
