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

Both `nonnegative.method="clamp"` and `"constrained"` are requested. Eigen
LLT, LDLT and CG use the constrained polish on eligible Gram systems.
Rectangular LSCG retains clamping under either setting. `Constraint_handling`
records this distinction. Four n=p=10,000 Gaussian LSCG fits used one
automatic LLT fallback each; their times include that fallback.

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
| OSQP | 62.823 | 77.586 | 2.639 | 54.088 |
| bigGlm | 11.799 | 12.688 | 9.546 | 7.100 |
| Eigen LLT | 0.274 | 0.691 | 0.496 | 0.584 |
| Eigen LDLT | 0.292 | 0.811 | 0.488 | 0.622 |
| Eigen CG-IC | 1.973 | 1.303 | 0.502 | 0.684 |
| Eigen LSCG-IC | 2.913 | 2.854 | 0.572 | 0.561 |

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
| 100 | 1,000 | 1.281 | 1.484 | 0.236 | 0.214 | 0.218 | 0.237 |
| 1,000 | 100 | 0.523 | 0.647 | 0.122 | 0.140 | 0.136 | 0.358 |
| 1,000 | 1,000 | 9.075 | 8.122 | 0.350 | 0.473 | 1.238 | 5.437 |
| 1,000 | 10,000 | 3.071 | 2.866 | 0.230 | 0.243 | 0.256 | 1.027 |
| 10,000 | 1,000 | 80.937 | 82.151 | 1.190 | 1.976 | 2.989 | 45.717 |

The n=100, p=1,000 case remains difficult for this 60-effect setup even with
ISIS: every engine selects two true effects and one false positive. This
limitation is retained in the tables and figures.

## Checks and output files

- All 222 configurations have exactly one measured fitting call, finite time,
  GIC, coefficient R-squared and coefficient RMSE, and reported convergence.
- All twelve dense p>n fits record `Screen="ISIS"`.
- True coefficients are identical across solvers and strategies for a given problem.
- Matched-strategy supports agree exactly. The largest relative coefficient
  difference from LLT is 1.14e-6. LSCG's constrained requests are compared with
  clamping because that is its effective strategy here.
- `metrics_all_conditions.csv` contains the complete results; the long version
  contains plotting values. `clamp_vs_constrained.csv` pairs the strategies,
  `dense_normal_results.csv` isolates the normal designs, and
  `solver_agreement.csv` records numerical/support agreement.
- `solver_size_comparison.pdf` combines all twelve vector plots. Separate PDF
  and PNG files provide deconvolution metrics/timings by strategy and the two
  dense normal-design sweeps.
- GIC is evaluated at fitted dispersion with coefficient-count degrees of freedom.
  The tables retain convergence, iteration and warning diagnostics.
