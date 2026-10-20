# One-million-predictor Gaussian benchmark

Run on 2026-10-20 with 500 observations, 1,000,000 iid standard-normal predictors, 60 uniformly distributed signed nonzero coefficients in [-0.8, 0.8], and Gaussian noise SD 1. The design/coefficient seed is 123 and the noise seed is 1. The coefficient L0/L1 ratio is 2.388402, matching the historical script's note.

R 4.5.0 on Windows 11; Intel Core i9-11950H, 64 GB RAM; one BLAS/OpenMP thread. Simulation and final reporting are excluded from the method timings. SIS/ISIS timings include screening and all outer fits; competing methods include their tuning paths and refits. These are individual elapsed-time runs, not replicated timing estimates.

The complete updated benchmark finished in one R process. Every method's original coefficient dimensions, Gaussian likelihood/GIC, F1 and MCC were independently checked. An initial combined run aborted inside ABESS with std::bad_alloc; the script now calls gc() before the native competitor fits. No fitting-package algorithm was changed.

Current results (current.csv):

| method | seconds | GIC | F1 | MCC | selected | TP | FP | FN |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Full L0glm | 43.01 | 2825.64 | 0.000 | NA | 0 | 0 | 0 | 60 |
| SIS + L0glm | 1.92 | 2799.96 | 0.301 | 0.394 | 13 | 11 | 2 | 49 |
| ISIS + L0glm | 17.26 | 2582.75 | 0.777 | 0.787 | 43 | 40 | 3 | 20 |
| ABESS | 53.31 | 2595.51 | 0.763 | 0.785 | 37 | 37 | 0 | 23 |
| Relaxed LASSO | 30.21 | 2825.64 | 0.000 | NA | 0 | 0 | 0 | 60 |

SIS screens 160 predictors once. ISIS admits up to 20 new candidates per round (plus current support), with start=NULL, gamma=2 and relaxed=TRUE; it returns the best GIC round. Full L0glm uses the same GIC preset and fitting options. ABESS searches consecutive sizes 0:150 with warm starts and no intercept. Relaxed LASSO requests ten lambda values; six admissible relaxed models remain after glmnet's near-saturation support restriction. Full L0glm at the GIC preset returns the null fit; GIC selects the null model along the relaxed-LASSO path. MCC is undefined (NA) for these constant selections.

Historical results, rescored using today's GIC convention (historical.csv):

| method | seconds | GIC | F1 | MCC | selected | TP | FP | FN |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ISIS + L0glm | 28.44 | 2627.88 | 0.776 | 0.796 | 38 | 38 | 0 | 22 |
| ABESS | 60.59 | 2772.94 | 0.442 | 0.481 | 26 | 19 | 7 | 41 |
| Relaxed LASSO | 23.33 | 5875.83 | 0.153 | 0.179 | 188 | 19 | 169 | 41 |

The historical values come from the [printed table at commit 0bfde60](https://github.com/tomwenseleers/L0glm/blob/0bfde60181a2d815f173084f3a3fc9c39c87316c/benchmarks/3_benchmarks_high%20dimensional_large%20problem.R#L1036). There were no historical ordinary full-L0glm or one-pass SIS results. Historical MCC was reconstructed from TP, FP, TN and FN; it was not printed in that table. Original printed GIC values remain in GIC_reported. Recalculation from the saved Gaussian RSS lowers them by approximately 2; the historical RSS values were rounded, so the recalculated scores inherit that rounding.

For all comparisons, with k nonzero coefficients and an estimated Gaussian dispersion, scoring uses

$$GIC=n[1+\log(2\pi RSS/n)]+(k+1)\log(p)\log\log(n).$$

The per-parameter GIC penalty is 25.239593; the BIC penalty is 6.214608. The original predictor count is retained after screening. Historical ISIS already used GIC. Historical ABESS used BIC scaled to the GIC penalty. Historical relaxed LASSO selected by BIC using the legacy likelihood adapter, so its selected model is not a like-for-like GIC comparison.

Holding the current fitting options fixed, the original 160-candidate ISIS batch gave GIC 2778.257, F1 0.347 and MCC 0.433. Reducing only that batch to 20 gave GIC 2582.751, F1 0.777 and MCC 0.787. start=1 produced the same coefficients as start=NULL with the 20-candidate batch. Using a BIC fitting penalty with the 160-candidate batch gave GIC 3587.848 and F1 0.231; outer-round selection still used common GIC. Thus changing GIC to BIC does not explain or fix the ISIS discrepancy on this simulation.

Copying the historical manual ISIS workflow into the current package (ten rounds, two unrelaxed gamma=2.1 pruning fits per round, followed by a gamma=2 fit) did not reproduce its old result: GIC 2799.041 and F1 0.286. The new public ISIS workflow with the smaller batch gives a different, better solution on today's implementation. This test does not isolate all historical package-version differences.

The current 20-candidate ISIS solution recovers 40 true features plus 3 false positives, versus 38 true features and no false positives historically. GIC is lower and F1 is similar; MCC is slightly lower because of the extra false positives. Lower GIC and good recovery do not certify a global optimum.

Historical ISIS used a different outer workflow, ABESS searched about 80 sizes and used its default intercept, and historical BLAS/thread settings were not recorded. Accordingly the time ratios in comparison.csv describe the recorded runs rather than a controlled implementation speedup.

The other files retain the matched screening controls, the new ISIS round history, the relaxed-LASSO path with current and legacy-adapter BIC scores, and sessionInfo(). Detailed local diagnostics and sparse coefficient checkpoints are under codex_sandbox/high_dimensional_benchmark_update/ and are intentionally ignored by Git.
