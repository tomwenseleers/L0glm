# Conditional posterior coverage in deconvolution

Reproduce from the repository root with `Rscript benchmarks/benchmark_posterior_coverage.R`
after installing the current package. Each selected fit uses the GIC preset,
nonnegative coefficients, default automatic support refinement, and 2,000 posterior
coefficient draws. An accepted unpenalized support refit uses a flat prior; retained
adaptive-ridge fits use their frozen inverse-square prior. Noise parameters are held
fixed in the posterior calculation. No prediction or fitting step uses the true coefficients.

Single-channel cases use n=p=1,000 and 100 true peaks with 200 noise replicates.
Three-channel cases use n=p=500 and 50 shared peaks with 100 replicates.
Gaussian noise has SD 5. Lower-signal cases reduce the peak coefficients twentyfold.
The reference restricts the fitting design to true predictors and uses lambda=0.
All 2,600 selected/reference fits and posterior calculations completed successfully.

| Channels | Noise | Coverage among recovered true coefficients | Coverage across all true coefficients | True-predictor reference |
|---:|:---|---:|---:|---:|
| 1 | gaussian estimated | 81.8% | 74.4% | 94.7% |
| 1 | gaussian known | 79.5% | 71.2% | 94.7% |
| 1 | gaussian moderate | 71.1% | 50.8% | 94.6% |
| 1 | poisson high | 74.0% | 62.3% | 94.9% |
| 1 | poisson moderate | 67.1% | 41.0% | 92.6% |
| 3 | gaussian estimated | 91.4% | 90.0% | 95.3% |
| 3 | gaussian known | 87.7% | 84.4% | 95.2% |
| 3 | poisson high | 86.2% | 74.8% | 94.2% |

Coverage is averaged within replicate, then across replicates. Missed true coefficients
count as uncovered in the all-true column. Monte Carlo SEs are recorded in the summary
CSVs; replicate CSVs retain support, prior, sampler, dispersion and coverage diagnostics.
Elapsed times include inference and were collected while independent regression checks
were running; they are not fitting-speed benchmarks.

The true-predictor references check fixed-support calculations. Selection and omitted
neighbouring peaks can substantially reduce coverage in the selected models. Conditional
credible intervals describe uncertainty within the admitted model and do not adjust for
support selection. Gaussian fixed-variance information and truncated-normal sampling are
also checked against analytic and independent rejection-sampling references in the unit tests.
