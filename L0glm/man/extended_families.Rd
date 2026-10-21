\name{extended_families}
\alias{extended_families}
\title{Negative Binomial, Cox and Multinomial Fits}
\description{
Additional likelihood families using existing ridge kernels and family-specific IRLS.
See the separate mathematical help pages for each family.
}
\details{
Use \code{L0glm.fit(X,y,family=mgcv::nb())} for estimated-shape NB,
\code{family=mgcv::nb(theta=2)} for fixed shape, or
\code{family=mgcv::nb(theta=-2)} to estimate from a starting shape of 2.
Each NB channel has its own shape parameter; scale dispersion is one.
Positive weights are prior/case weights. Offsets are included in the link.

Use \code{family=mgcv::cox.ph()} with right-censored
\code{survival::Surv(time,status)} or an ordinary two-column time/status matrix.
Status belongs in the response, not in weights. Cox uses weighted Breslow ties,
partial likelihood, and no intercept. Delayed entry, strata, time-varying
predictors, multiple states and other tie methods are not implemented.
Survival predictions use an estimated Breslow baseline hazard.

Use \code{family=mgcv::multinom(K=2)} for three categories, supplied as a factor
or three-column integer count matrix. All K+1 category coefficients are returned
in centered log-ratio (CLR) coordinates: every coefficient row sums to zero.
A free intercept is added by default. Starts, offsets and prior means also
have K+1 columns and are centered across categories. Common shifts have no
effect on probabilities. All-zero starting rows exclude predictors before CLR
centering; individual category zeroes in other rows remain eligible.
Scalar start=1 centers to zero while keeping all rows eligible; start=0 fits
only mandatory terms. Intercepts and always.include rows remain eligible even
with zero starts. Eligible zero starts receive one weak-ridge IRLS step before
adaptive weighting. Returned sparse coefficients retain their original
predictor and category dimensions, including excluded zero rows. Nonnegative=TRUE is ignored
with a warning. The model and output stay in CLR; K orthonormal contrast
coordinates are used only to solve nonsingular linear systems efficiently.

\code{prior.mean} supplies ridge-penalty centers in original input coefficient
units: a scalar, length-\eqn{p} vector, or \eqn{p} by \eqn{m} matrix, with
any added intercept included in \eqn{p}. NB uses one column per response
channel; Cox uses one column, no intercept, and centers log-relative-hazard
coefficients rather than the baseline hazard. NB shape is not a prior target.
Multinomial uses all K+1 category columns and centers each target row to CLR;
scalar or vector targets shared across categories become zero. \code{NULL}
also uses zero targets. A nonzero target cannot override an excluded starting
coefficient or multinomial row, which remains zero in the full-size sparse
output. Mandatory terms remain unpenalized. Adaptive weights and support counts
still use coefficients relative to zero, rather than deviations from the
targets. Special-family relaxed refits remove the ridge penalty on retained
support, so the targets no longer exert shrinkage during that refit.

Multinomial uses coupled block IRLS with the exact class covariance. IRLS for
logistic models is algebraically equivalent to a Newton/Fisher-scoring step,
expressed as weighted least squares. It is not independent binomial fitting.
The alternative \code{control.family=list(multinomial="bound")} uses a diagonal
upper curvature bound on the zero-sum plane inspired by Krishnapuram et al.; their paper uses an L1
penalty, while the adaptive L0 extension here is separate.

Grouped selection penalizes a predictor row if any channel or CLR category
coefficient is nonzero. Its adaptive ridge weights are shared across a row.
NB ungrouped selection counts individual channel coefficients. Multinomial
ungrouped selection counts individual CLR entries, with zero-sum rows; it still
fits the joint multinomial likelihood. Grouped multinomial selection instead
counts active rows using a category-invariant CLR row norm. Its parameter-count
IC charges K free directions per active row. Cox has one coefficient per
predictor, so grouped and ungrouped selection coincide.

For normalized working coefficients \eqn{b}, the unsmoothed adaptive
weights use a fixed exponent of two:
\eqn{D_{jc}=\lambda/(4\|b_j\|_2^2)} for grouped multichannel NB and
multinomial models, and \eqn{D_{jc}=\lambda/(4|b_{jc}|^2)} for ungrouped
or single-channel fits, including Cox. Scale dispersion is one.
\code{thresh} (default \code{1e-5}) acts on row norms for grouped selection,
removing whole rows, and on individual magnitudes within each response/category
column for ungrouped or single-channel selection. The relative cutoff is
\code{thresh} times the corresponding largest selectable magnitude; an
absolute \code{1e-20} cutoff also applies. Mandatory terms are excluded from
these maxima and never thresholded. Retained grouped rows do not receive
individual-entry hard cutoffs. Multinomial surviving entries are recentered
within their allowed category set, keeping removed entries zero; rows with
fewer than two allowed entries are zero. Cutoffs apply only in selection passes
with positive lambda. \code{thresh=0} retains the absolute cutoff, and
\code{thresh} does not smooth adaptive weights.

With \code{relaxed=TRUE}, special families refit using \code{lambda=0} without
hard thresholding. Grouped refits admit all eligible entries in retained rows;
ungrouped refits admit only entries nonzero at the end of selection. A grouped
multinomial refit can recover category entries within a retained row, whereas
an ungrouped refit holds unselected CLR entries at zero. Both retain zero-sum
rows: s admitted category entries give s-1 free directions. No new predictor
rows enter the refit. Starting exclusions remain binding, and mandatory terms
remain eligible and unpenalized. Estimated NB shapes continue to be updated.
Sparse outputs retain all original coefficient dimensions. Relaxed fitting is
enabled by default for grouped multichannel NB and all grouped multinomial
responses, but not for Cox, ungrouped fits or single-channel NB. Each pass uses
\code{maxit}; \code{group.maxit} is unused. Tuning applies refits to each
candidate and within CV training folds before scoring; pilot fits are unrelaxed.
The reported L0 objective is a target for a heuristic, not a global certificate.
There is no reference category in multinomial fitting or selection.

\code{L0glm.metrics} reports exact family likelihoods when mu_floor=0, the
family deviance, the target L0 objective, row GIC and parameter-count ICs.
Object metrics and \code{logLik} use exact likelihoods. External NB scoring
defaults to a mean floor of 1e-8. Cox deviance is -2 times partial log likelihood,
not an ordinary saturated-model deviance. Grouped IC tuning uses row counts;
\code{metrics$tuning_score} is the actual selection score, while \code{gic}
retains the parameter-count definition. CV returns full-data best and one-SE fits.
NB CV keeps theta fixed at its training estimate when scoring validation data.
Cox CV uses full-minus-training partial likelihood, preserving full risk sets.

These loops are serial. Cox caches risk sets and uses matrix-free PCG for
large, wide or sparse unrestricted Eigen systems, retaining LLT for small dense
active systems. Its weak-ridge initializer is also matrix-free. Direct Cox and
block multinomial curvature systems are guarded by control.family$max.working.
Wide direct systems use whitening and the existing Woodbury ridge route. The bound option changes curvature, not the definition of the likelihood.
No new compiled code or dependency imports are required; the mgcv constructors
are optional and reference packages are listed in Suggests.
}
\references{
Krishnapuram, B., Carin, L., Figueiredo, M. A. T. and Hartemink, A. J. (2005).
Sparse multinomial logistic regression: Fast algorithms and generalization bounds.
IEEE TPAMI 27(6), 957--968.

Verweij, P. J. M. and van Houwelingen, H. C. (1993).
Cross-validation in survival analysis. Statistics in Medicine 12, 2305--2314.
}
\seealso{\code{\link{L0glm.fit}}, \code{\link{L0glm.tune}},
\code{\link{L0glm.metrics}}, \code{\link{L0glm_negative_binomial}},
\code{\link{L0glm_multinomial}}, \code{\link{L0glm_cox}}}
