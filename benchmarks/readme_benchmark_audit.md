# README benchmark: statistical definitions and abess investigation

The executable benchmark is [readme_benchmarks_L0glm_deconvolution.R](readme_benchmarks_L0glm_deconvolution.R).
The independent checks and the small abess reproducer are in
[diagnose_readme_benchmarks.R](diagnose_readme_benchmarks.R).
The investigation used abess 0.4.11 and glmnet 5.0 on 12 October 2026.

## Why the original abess results were unreliable

The coefficient extraction was correct: `coef.abess()` includes the intercept
as its first row, including when its value is zero. The problem was the
meaning of `beta.low = 0` in the Gaussian fitting code. It solves unconstrained
least squares and then truncates the coefficients to the requested bounds.
For correlated predictors, truncation does not solve nonnegative least squares.
The source uses `ldlt().solve(...)` followed by `trunc(beta, ...)` in both the
[Gaussian](https://github.com/cran/abess/blob/0.4.11/src/AlgorithmGLM.h#L405-L427)
and [multivariate Gaussian](https://github.com/cran/abess/blob/0.4.11/src/AlgorithmGLM.h#L981-L1003)
fits. The same clipping pattern was present in the upstream development code
when checked.

A counterexample with no normalization, no intercept, and both predictors
available is:

```r
X <- cbind(rep(1, 4), 1:4)
y <- c(1, 0, 0, 0)
fit <- abess::abess(X, y, family = "gaussian", support.size = 2,
                    normalize = 0, fit.intercept = FALSE,
                    beta.low = 0, lambda = 0, num.threads = 1)
as.vector(fit$beta)  # abess 0.4.11: c(1, 0)
```

Unconstrained least squares gives `(1, -0.3)`. Clipping gives `(1, 0)` and
residual sum of squares 3, which is even worse than the all-zero fit's value 1.
The nonnegative solution is `(0.25, 0)`, with residual sum of squares 0.75.
Its gradient satisfies the nonnegative least-squares optimality conditions:
the first component is zero and the second is positive at its lower bound.
The reproducer confirms the problem for dense and sparse designs and for both
Gaussian families. This establishes an implementation limitation in bounded
fitting; it does not imply that unconstrained abess is generally defective.

Highly overlapping convolution columns can yield large positive and negative
OLS coefficients that cancel in the prediction. Removing the negative terms
destroys that cancellation, explaining the enormous predictions in the
original benchmark. A requested support size of 100 could also leave fewer
than 100 nonzero coefficients after clipping; that discrepancy was not an
indexing error in coefficient extraction.

Abess is omitted from the README comparison. A LASSO support could initialize
its search through `init.active.set`, but this would not change the incorrect
bounded least-squares step. Proper NNLS refitting would require a separate
hybrid procedure or a change to the abess fitting code. The small reproducer
is retained for reference; the README benchmark now compares only L0glm and
glmnet and does not require abess. The installed abess package is unchanged.

## Multichannel weights

For `glmnet(..., family = "mgaussian")`, `weights` must have length `nrow(X)`.
The installed glmnet 5.0 rejects an `n x m` matrix; the
[glmnet entry point](https://github.com/cran/glmnet/blob/master/R/glmnet.R)
checks this length, and its
[multivariate Gaussian driver](https://github.com/cran/glmnet/blob/master/R/mrelnet.R)
uses the same weights for every response column. A response matrix does not
imply support for a matrix of precision weights.

The multichannel glmnet fit uses the proxy
`1 / (rowMeans(Y) + 0.1)`. This is not the elementwise Poisson IRLS weight
`1 / mu[i, j]`: it loses channel-specific variance information and is fixed
rather than updated. Exact elementwise weighting would require a different
multitask solver or a separately implemented algorithm. L0glm uses the actual
Poisson family with channel-specific IRLS weights.

## Noise-free deviance

Let `B0 = beta_true`, `mu0 = X %*% B0 = y_true`, and `muhat = X %*% Bhat`.
With $\widetilde\mu_{ic}=\max(\widehat\mu_{ic},10^{-4})$, the reported quantity is

\[
D_{\rm true}=2\sum_{i,c}\left[
\mu_{0,ic}\log\left(\frac{\mu_{0,ic}}{\widetilde\mu_{ic}}\right)
-\mu_{0,ic}+\widetilde\mu_{ic}\right].
\]

The zero-times-log term is defined as zero. The reference means are neither
rounded nor replaced by observed counts. Predictions in this comparison are
floored at the fixed `cfg$deviance_floor = 1e-4`; all observations are retained.
Observed-data likelihood/GIC and path selection retain `cfg$mu_floor = 1e-8`,
so this reporting convention does not change the fitted coefficients. The diagnostic
script verifies `y_true == X %*% beta_true` and recalculates the sum directly,
independently of R's `poisson()$dev.resids` implementation.

This is twice the sum of Poisson KL divergences using the floored predictions.
For small prediction errors at true means above the floor, its local quadratic approximation is

\[
D_{\rm true}\approx\sum_{i,c}
\frac{(\widehat\mu_{ic}-\mu_{0,ic})^2}{\mu_{0,ic}}.
\]

It therefore need not rank models like coefficient R2, support F1/MCC, or
unweighted squared prediction error. In addition, convolution can make quite
different coefficient vectors yield similar signals.

The original saved multichannel scores were independently reproduced against
`y_true`. Their ranking arose from larger contributions at weak positive
true means for L0glm, despite its smaller coefficient and signal squared
errors. Zero predictions must not be excluded: predicting zero at a positive
true mean is a model error and gives infinite unfloored Poisson KL divergence.
[R's Poisson `dev.resids`](https://github.com/wch/r-source/blob/trunk/src/library/stats/R/family.R)
returns `c(0, Inf)` for `y = c(0, 1)` and
`mu = c(0, 0)`. `glm.fit` sums deviance contributions across observations;
it does not drop an observation merely because its prediction is zero.
For its default log link, R already floors `exp(eta)` at
`.Machine$double.eps`, and its Poisson family requires positive fitted means.
The benchmark's chosen `1e-4` floor makes this external comparison finite;
it is not presented as R's default numerical threshold.

## Likelihood, entry penalties, row penalties, and reference values

For the observed integer counts `Y`, the Poisson log likelihood is

\[
\ell(B;Y)=\sum_{i,c}\{Y_{ic}\log (XB)_{ic}-(XB)_{ic}-\log\Gamma(Y_{ic}+1)\}.
\]

The table now exposes `neg2loglik`, `nonzero_coefficients`, `selected`, and
`lambda_L0`, and reports both

\[
\mathrm{GIC}(B)=-2\ell(B;Y)+\lambda K_{\rm entries}(B),\qquad
\mathrm{GIC}_{\rm rows}(B)=-2\ell(B;Y)+\lambda K_{\rm rows}(B),
\]

where `lambda = log(p) * log(log(n))`, `n = nrow(X)`, and `p = ncol(X)`.
The counts use the stated absolute threshold `1e-5`. These two criteria agree
for one channel. For multiple channels they answer different questions; a row
count is not a count of all fitted scalar parameters.

The previous benchmark used `n * m` inside the reported multichannel penalty
constant and tuned Gaussian comparators using a different Gaussian criterion.
The revised code matches the fitter's `n` convention and uses the same
observed-data, entry-count Poisson GIC for all comparator path selection.
It saves each tuning path to permit checking that the reported candidate is
the one attaining the smallest GIC on that path.

The current grouped L0 implementation runs channelwise adaptive-ridge fits
with shared support-based starting values, followed by relaxation when
requested. It is a support-sharing heuristic; it does not explicitly solve
the joint row-count objective above. Even the adaptive-ridge approximation
does not certify a global L0 minimum. Accordingly, `GIC_rows` is an explicit
evaluation criterion, not a claim that the algorithm exactly optimizes it.

`GIC_beta_true` and `GIC_rows_beta_true` evaluate the true coefficient matrix
against the same **observed** data. They are oracle references, not theoretical
minima: fitting noise or omitting sufficiently weak true coefficients can
lower the realized penalized objective.
In the ideal unfloored model, write $D_0(B)$ for the same deviance expression
without the prediction floor. For either choice of count `K`,

\[
E\{J(B;Y)-J(B_0;Y)\}=D_0(B)+\lambda\{K(B)-K(B_0)\}.
\]

Thus even the expected penalized objective can favor dropping a weak true
component if the saved penalty exceeds its increase in KL divergence.

For a rigorous bound, let `ell_sat = ell(mu = Y; Y)` be the saturated likelihood.
Because `ell(B; Y) <= ell_sat` and both penalties are nonnegative,

\[
-2\ell(B;Y)+\lambda K(B)\ \geq\ -2\ell_{\rm sat}.
\]

The column `objective_lower_bound` reports this bound. It may not be attainable
by any sparse nonnegative coefficient matrix for the given design. The exact
global penalized minimum is not known merely by knowing `beta_true`.
For unfloored noise-free deviance, the theoretical minimum is zero at the
true mean. The reference table also gives its actual numerical value under
the fixed comparison floor; this can be positive at truly zero means.

Reference values are saved in
[readme_benchmark_references.csv](readme_benchmarks/readme_benchmark_references.csv).
All method panels include elapsed fitting/selection time in their headings.
