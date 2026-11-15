# Solver size and nonnegative-strategy benchmarks

This experiment reports **222 fits**, each timed once with no fitting warmup:
192 nonnegative deconvolution fits and 30 dense normal-design Gaussian fits.
Every fit returned finite metrics and converged within the 500-iteration ceiling.
The fitting package is L0glm 0.2.0.9002, including the `nonnegative.method`
interface. The included results were obtained on 15 November 2026 on the
HP ZBook Fury G8 workstation with an Intel Core i9-11950H processor
(8 physical / 16 logical cores; nominal 2.60 GHz). Every fit uses one core and
one BLAS thread.

## Reproduce

Install the current L0glm package and the benchmark dependencies `callr`,
`ggplot2`, `tidyr`, `tidyselect`, `microbenchmark`, `RhpcBLASctl`, `scales`,
`osqp`, and `glmnet`.
From the repository root:

```sh
Rscript benchmarks/1_benchmarks_L0glm_solvers_problem_sizes.R
```

The script writes this directory. `--plots-only` redraws stored results without
fitting. `--resume` reuses completed per-fit checkpoints only when the manifest
matches the cases, settings, helper code, package DLL and R version. `--smoke`
runs smaller validation cases in a separate `smoke/` subdirectory.
Size grids, constraint strategies, iteration ceiling, thread counts and the
per-fit timeout are selectable at the top of the script. The helper file
`benchmarks/solver_problem_sizes_helpers.R` handles data generation and the
single fitting call.

Fresh sequential workers isolate caches and enforce a 900-second per-fit
timeout. Only `l0glmFit()` is timed, using `microbenchmark::get_nanotime()`;
worker startup, namespace loading, simulation, metric calculation and plotting
are excluded. First-call fitting/dispatch costs are included. These are
single-run measurements, so small timing differences should not be overinterpreted.

## Nonnegative deconvolution

The Gaussian and Poisson examples use the same PSF simulation and seeds for
all solvers, with `n=p` in 100, 500, 1,000 and 10,000, ten percent true peaks,
dense or sparse input storage, `start=1`, `lambda="gic"`, `screen="none"`,
`local.search=FALSE`, `maxit=500`, and `tol=1e-10`. Gaussian simulation uses
unit noise and truncates negative observations at zero. The "dense" PSF
design contains exact structural zeros; native compression is allowed.

Both `nonnegative.method="clamp"` and `"constrained"` are requested. All six
engines use the requested strategy. Constrained Eigen LSCG retains its
least-squares solve and applies a checked bound correction when needed, using
a predictor Gram matrix subject to the memory budget. Initialization defaults
and the guarded LLT fallback are unchanged. `Constraint_handling` records the
strategy used. 2 fits reported an automatic LLT fallback; their
times include it.

The package's automatic initialization remains unchanged. These controlled
deconvolution benchmarks explicitly use `start=1`; the dense normal designs
use `start=NULL`.

For this **unit-start** setup, clamping gives lower GIC and higher coefficient
R-squared at n=p=10,000 for both noise families. The following quality values
are from LLT with sparse PSF storage; other engines agree when their effective
constraint strategy is matched.

| Noise | Strategy | GIC | Coefficient R2 | Coefficient RMSE | TP | FP | F1 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| gaussian | clamp | 90,472.03 | 0.995974 | 116.553 | 944 | 91 | 0.9278 |
| gaussian | constrained | 105,662.97 | 0.989072 | 193.965 | 924 | 123 | 0.9028 |
| poisson | clamp | 102,256.13 | 0.921438 | 514.972 | 843 | 156 | 0.8434 |
| poisson | constrained | 102,360.15 | 0.880574 | 634.914 | 839 | 163 | 0.8382 |

Fitting times in seconds at n=p=10,000 with sparse PSF storage:

| Solver | Gaussian clamp | Gaussian constrained | Poisson clamp | Poisson constrained |
| --- | --- | --- | --- | --- |
| OSQP | 60.167 | 77.556 | 1.173 | 23.671 |
| bigGlm | 10.867 | 11.868 | 4.504 | 3.672 |
| Eigen LLT | 0.295 | 0.751 | 0.306 | 0.426 |
| Eigen LDLT | 0.311 | 0.778 | 0.342 | 0.454 |
| Eigen CG-IC | 1.907 | 1.205 | 0.398 | 0.497 |
| Eigen LSCG-IC | 2.840 | 1.723 | 0.462 | 0.596 |

This is an empirical comparison of the complete adaptive-ridge trajectories.
Solving the constrained ridge subproblem more accurately need not improve
the final nonconvex support-selection criterion.

## Dense Gaussian regression with 60 signed effects

The simulation follows the homepage normal-design example: independent
standard-normal predictors, 60 randomly selected effects uniform on
[-0.8,0.8], and independent unit Gaussian errors. Predictor/effect seed is
123 and noise seed is 1. Two sweeps vary p with n=1,000, and n with p=1,000,
using 100, 1,000 and 10,000 (five unique settings).

These fits use `nonnegative=FALSE`, `start=NULL`, the GIC preset, and
`local.search=FALSE`. **ISIS is explicitly enabled whenever p>n**, including
its fitting/screening iterations in the timer; n>=p fits use `screen="none"`.
This benchmark choice is separate from the package's automatic screening
trigger. The change in workflow at p>n can make timing curves nonmonotone.
No duplicate clamp/constrained fits are needed for unrestricted coefficients.

All six engines recover the same support for each setting. Representative
quality metrics from LLT are:

| n | p | Screen | GIC | Coefficient R2 | TP / 60 | FP | F1 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 100 | 1,000 | ISIS | 556.48 | 0.028155 | 2 | 1 | 0.0635 |
| 1,000 | 100 | none | 3,368.93 | 0.994835 | 56 | 0 | 0.9655 |
| 1,000 | 1,000 | none | 3,595.86 | 0.991664 | 54 | 0 | 0.9474 |
| 1,000 | 10,000 | ISIS | 3,898.02 | 0.981424 | 49 | 0 | 0.8991 |
| 10,000 | 1,000 | none | 29,497.08 | 0.999202 | 56 | 0 | 0.9655 |

Fitting times in seconds:

| n | p | OSQP | bigGlm | Eigen LLT | Eigen LDLT | Eigen CG-IC | Eigen LSCG-IC |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 100 | 1,000 | 0.973 | 1.308 | 0.195 | 0.193 | 0.178 | 0.211 |
| 1,000 | 100 | 0.519 | 0.525 | 0.122 | 0.109 | 0.116 | 0.268 |
| 1,000 | 1,000 | 6.494 | 6.590 | 0.314 | 0.388 | 1.306 | 4.400 |
| 1,000 | 10,000 | 2.010 | 2.452 | 0.236 | 0.253 | 0.312 | 1.006 |
| 10,000 | 1,000 | 74.266 | 77.386 | 1.188 | 2.018 | 2.766 | 45.570 |

The n=100, p=1,000 case remains difficult for this 60-effect setup even with
ISIS: every engine selects two true effects and one false positive. This
limitation is retained in the tables and figures.

## Checks and output files

- All 222 configurations have exactly one measured fitting call, finite time,
  GIC, coefficient R-squared and coefficient RMSE, and reported convergence.
- All twelve dense p>n fits record `Screen="ISIS"`.
- True coefficients are identical across solvers and strategies for a given problem.
- Matched-strategy supports agree exactly. The largest relative coefficient
  difference from LLT is reported in `solver_agreement.csv`. Constrained LSCG
  is compared directly with constrained LLT, including rectangular solves.
- Across all 16 constrained LSCG fits, the largest relative coefficient
  difference from LLT is 1.39e-08, and none requires LLT fallback.
  All 206 configurations unaffected by the constraint correction return
  bitwise identical coefficients in the regression check.
- `metrics_all_conditions.csv` contains the complete results; the long version
  contains plotting values. `clamp_vs_constrained.csv` pairs the strategies,
  `dense_normal_results.csv` isolates the normal designs, and
  `solver_agreement.csv` records numerical/support agreement.
- `solver_size_comparison.pdf` combines all twelve vector plots. Separate PDF
  and PNG files provide deconvolution metrics/timings by strategy and the two
  dense normal-design sweeps.
- GIC is evaluated at fitted dispersion with coefficient-count degrees of freedom.
  The tables retain convergence, iteration and warning diagnostics.
