# OpenMP SIS/ISIS benchmark

Run on 2026-10-20 using R 4.5.0 on Windows, an Intel Core i9-11950H (8 physical / 16 logical cores) and 64 GB RAM. The design has 500 observations and 1,000,000 iid standard-normal predictors (about 4 GB), with the original 60 signed coefficients in [-0.8, 0.8], problem seed 123 and Gaussian noise SD 1 with seed 1. BLAS and model-fitting threads were held at one; only screening concurrency varied.

Reproduce with [3_benchmarks_sis_openmp_threads.R](../3_benchmarks_sis_openmp_threads.R) after installing the current package. Simulation, explicit garbage collection, equality checks and final metric calculation are outside the timed calls. SIS uses `return_stats=FALSE`, retains 160 candidates and includes null-fit preparation. Complete ISIS uses the automatic 20-candidate batch, up to ten rounds, GIC and relaxed refits. The order of thread counts varies across repetitions.

| Screening threads | SIS median seconds (5 runs) | SIS speedup | Complete ISIS median seconds (3 runs) | ISIS speedup |
| --- | ---: | ---: | ---: | ---: |
| 1 | 0.80 | 1.00 | 8.41 | 1.00 |
| 2 | 0.50 | 1.60 | 6.47 | 1.30 |
| 4 | 0.33 | 2.42 | 4.50 | 1.87 |
| 8 | 0.24 | 3.33 | 3.31 | 2.54 |

All requested thread counts were actually used. Individual measurements are in `timings.csv`; `summary.csv` also records minima and maxima. Speedups use the new one-thread implementation as the baseline, isolating added screening concurrency. The kernel also replaces Rcpp matrix access with read-only numerical views and reuses per-thread working buffers. These are measurements on this machine, not universal speed guarantees.

Every million-predictor score, sign, rank and selected candidate matched exactly across thread counts. Every ISIS coefficient and round history also matched exactly, including the previously validated serial fit. All complete fits returned GIC 2582.75062779836, F1 0.776699029126214 and MCC 0.787488992969348, selecting 43 predictors (40 true positives, 3 false positives, 20 false negatives). The original-dimension sparse coefficient matrix is retained.

The package suite passed 6,711 expectations with no failures, test warnings or skips. Checks cover 52 family/link contexts: all standard Gaussian, binomial, Poisson, Gamma and inverse-Gaussian links; quasi variants; fixed and estimated negative-binomial theta with log, identity and square-root links; multinomial counts/factors; independent channels; and weighted tied Cox risk sets. Dense and sparse scores matched the preceding serial implementation bit for bit in all 416 comparisons (52 contexts x 2 storage types x 4 thread counts), summarized in `family_identity.csv`. Full-statistic and top-only paths, nuisance terms, offsets, zero observation weights, one-sided scoring, tie ordering, invalid-data errors, CV and complete special-family fits are tested. A separately compiled kernel with OpenMP disabled matched all 104 dense/sparse family references and correctly reported one thread.

Each predictor keeps its original summation order. Workers retain private top-k heaps under the total order `(score, -original_column_index)`. A candidate outside its worker's top k has at least k predecessors in that worker, so it cannot be in the global top k. Merging the local heaps therefore gives the exact serial selection, including ties. No floating-point reductions span workers. R allocation, error reporting and interrupt checks remain on the main thread, following [Writing R Extensions](https://www.stat.ethz.ch/R-manual/R-devel/doc/manual/R-exts.html).

Usage after rebuilding/reinstalling the package:

```r
screen <- sis_glm(X, y, nsis = 160, return_stats = FALSE, nthreads = 8)
fit <- L0glm.fit(X, y, screen.nthreads = 8, nthreads = 1)
fit$screening$nthreads
fit$screening$history
```

Screening defaults to one thread to preserve existing resource use. Counts are capped at the available logical processors and predictor count; builds without OpenMP use one thread. `screen.nthreads` is independent of `nthreads` (BLAS), `ncores` and `parallelize` (outcome fitting), and is ignored when screening is disabled. No environment-variable configuration is required. Build logs, raw family comparisons and coefficient checkpoints are in the ignored `codex_sandbox/sis_openmp` directory.
