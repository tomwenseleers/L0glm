# OpenMP SIS/ISIS benchmark

The design has 1,000 observations, one million iid standard-normal predictors (8 GB), 60 signed effects in [-0.8, 0.8], problem seed 123 and noise seed 1. HP ZBook Fury G8, Intel i9-11950H, eight physical cores, 16 logical cores and 64 GB RAM; Windows and R 4.5.0. BLAS and fitting use one thread; screening concurrency varies.

SIS retains 288 candidates. Complete ISIS uses the GIC-based batch, up to ten rounds, guarded support refitting and automatic local refinement. Simulation, explicit GC, equality checks and final reporting are outside timed calls.

| Screening threads | SIS seconds (median of 5) | SIS speedup | ISIS + L0glm seconds (median of 3) | ISIS speedup |
|---:|---:|---:|---:|---:|
| 1 | 1.77 | 1.00 | 12.08 | 1.00 |
| 2 | 1.50 | 1.18 | 11.22 | 1.08 |
| 4 | 0.84 | 2.11 | 6.36 | 1.90 |
| 8 | 0.44 | 4.02 | 3.79 | 3.19 |

Every tested score, sign, rank, coefficient and ISIS history agrees exactly across thread counts. The original-dimension sparse coefficients are retained. All fits give GIC 4195.01, F1 0.868 and MCC 0.876.

Workers retain private top-k heaps ordered by score and original column index. Each predictor keeps its summation order, and merging the heaps preserves exact serial selection, including ties. No floating-point reduction spans workers.

Reproduce with [3_benchmarks_sis_openmp_threads.R](../3_benchmarks_sis_openmp_threads.R). `timings.csv`, `summary.csv` and `isis_history.csv` retain the full results.
