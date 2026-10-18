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
effect on probabilities. Scalar start=1 centers to zero, as does start=0;
zero starts receive one weak-ridge IRLS step before adaptive weighting.
Starting zeroes impose no support restriction. Nonnegative=TRUE is ignored
with a warning. The model and output stay in CLR; K orthonormal contrast
coordinates are used only to solve nonsingular linear systems efficiently.

Multinomial uses coupled block IRLS with the exact class covariance. IRLS for
logistic models is algebraically equivalent to a Newton/Fisher-scoring step,
expressed as weighted least squares. It is not independent binomial fitting.
The alternative \code{control.family=list(multinomial="bound")} uses a diagonal
upper curvature bound on the zero-sum plane inspired by Krishnapuram et al.; their paper uses an L1
penalty, while the adaptive L0 extension here is separate.

Grouped selection penalizes a predictor row if any contrast or channel is
nonzero. Its adaptive ridge weights are shared across a row. Relaxed refits
can recover coefficients within retained rows, while respecting explicit
structural-zero starts for NB and Cox. Multinomial ungrouped L0 counts nonzero
CLR entries, with row sums kept exactly constrained during fitting. A row with
s active entries has s-1 free directions. Grouped multinomial L0 counts active
rows and uses a category-invariant CLR row norm. Its parameter-count IC charges
K free directions per active row. Mandatory terms remain unpenalized.
Hard thresholding is disabled during relaxed fitting.
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
