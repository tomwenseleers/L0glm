\name{L0glm_multinomial}
\alias{L0glm_multinomial}
\title{Multinomial Regression with Grouped or CLR Entrywise L0 Penalties}
\description{
Joint multinomial IRLS with centered log-ratio (CLR) coefficients, an invariant
grouped L0 penalty, and an optional entrywise CLR L0 penalty.
}
\details{
For \eqn{C} categories, use \code{family=mgcv::multinom(K=C-1)}.
Supply a factor/class vector or an integer count matrix with \eqn{C} columns.
Every category must have positive weighted training count. For count vectors,
\eqn{N_i=\sum_cY_{ic}} and
\deqn{p_{ic}=\exp(\eta_{ic})/\sum_h\exp(\eta_{ih}),\qquad
 \eta_i=x_i^TB+o_i,\qquad B\mathbf{1}=0.}
Thus \eqn{\eta_{ic}=\log p_{ic}-C^{-1}\sum_h\log p_{ih}} for centered offsets.
The weighted likelihood includes the multinomial count constants.
A factor is a one-trial count response.

All \eqn{C} coefficient columns are returned in CLR coordinates, in original
predictor units. An unrestricted, unpenalized intercept is added by default.
\code{predict(type="link")} returns centered scores;
\code{predict(type="response")} returns probabilities; and
\code{predict(type="class")} returns class labels.
Starting matrices, prior means and offsets have \eqn{C} columns in response
category order and are centered across categories. A starting row supplied as
all zero excludes that predictor before CLR centering, throughout selection
and relaxed refitting. Individual category zeroes within other starting rows
remain eligible. Intercepts and \code{always.include} rows remain eligible
even with zero starts. Scalar \code{start=1} centers to zero but keeps all rows
eligible; \code{start=0} fits only mandatory terms. A weak-ridge step seeds
eligible zero starts before adaptive weighting. The returned sparse coefficient
matrix retains all original predictor rows and \eqn{C} category columns, with
excluded rows filled with zeroes. Use \code{start=NULL} for automatic initialization.
\code{nonnegative=TRUE} is ignored with a warning.

\code{prior.mean} supplies ridge-penalty targets in original predictor units.
Supply a \eqn{p} by \eqn{C} matrix in response-category order, including any
added intercept row. Targets are transformed internally and centered to CLR:
\deqn{a_{jc}\leftarrow a_{jc}-C^{-1}\sum_{h=1}^C a_{jh}.}
Common row shifts have no effect. A scalar or length-\eqn{p} vector is repeated
across categories and therefore gives zero CLR targets; nonzero targets require
category contrasts. \code{NULL} uses zero targets.
The quadratic ridge term is centered on these targets, whereas adaptive weights
and support counts still use coefficients relative to zero. Targets are not
starting coefficients and do not override starting-row exclusions: a nonzero target
for an excluded row leaves that row fixed at zero in the full-size sparse output.
Unpenalized intercept/mandatory rows are not shrunk toward their targets.
The final relaxed refit removes the ridge penalty on retained support, so the
targets no longer exert shrinkage during that refit.
}
\section{Joint quadratic IRLS}{
The observation score and information on class scores are
\deqn{s_i=w_i(Y_i-N_ip_i),\qquad
 W_i=w_iN_i\{\operatorname{diag}(p_i)-p_ip_i^T\}.}
The class covariance has null direction \eqn{\mathbf{1}}.
An orthonormal basis \eqn{Q} with \eqn{Q^TQ=I} and
\eqn{Q^T\mathbf{1}=0} eliminates that redundant direction in the linear solve.
Writing \eqn{B=\Gamma Q^T} preserves \eqn{\|B_j\|_2=\|\Gamma_j\|_2}.
The exact class covariance, including off-diagonal class coupling, is used.
These are weighted least-squares IRLS steps, algebraically equivalent to
Newton/Fisher scoring for multinomial logit. Independent marginal binomial
fits, or centering their coefficients, do not solve this likelihood.

Tall working systems assemble coefficient-space curvature. Wide systems
whiten observation/class blocks and use the existing Woodbury ridge path.
Dense or stacked working arrays can arise even for sparse predictors.
\code{control.family$max.working=2000L} guards the smaller working dimension.
The optional \code{control.family=list(multinomial="bound")} replaces exact
curvature by the upper bound \eqn{0.5I} on the contrast plane, while retaining
the exact likelihood and constraints. The bound is inspired by Krishnapuram
et al.; their sparse penalty is L1, not the adaptive L0 procedure here.
}
\section{Selection, invariance and degrees of freedom}{
With \code{grouped=TRUE}, the target is
\deqn{-2\ell(B)+\lambda\sum_{j\notin M}
  \mathbf{1}\{\|B_j\|_2>0\},}
where \eqn{M} contains mandatory terms. In normalized working CLR coefficients
\eqn{b}, adaptive ridge diagonals are
\eqn{D_{jc}=\lambda/(4\|b_j\|_2^2)}.
The exponent is fixed at two, with \eqn{\|b_j\|_2^2=\sum_h b_{jh}^2}.
Positive and negative class effects cannot cancel in this norm.
Both the row penalty and the model are invariant to the choice of orthonormal
contrast basis and to permutations of category labels.

With \code{grouped=FALSE}, the charge counts nonzero displayed CLR entries,
and \eqn{D_{jc}=\lambda/(4|b_{jc}|^2)}. This is a different sparsity target from
entrywise ILR sparsity. Each active row is solved on the zero-sum plane of its
allowed class entries; a row with fewer than two entries must be zero.
Adaptive weights are unsmoothed and use coefficients relative to zero, not
deviations from \code{prior.mean}; \code{thresh} does not alter them.
With \code{grouped=TRUE}, hard thresholding acts on
\eqn{r_j=\|b_j\|_2}, not individual category entries or arithmetic row means.
An entire row is removed if
\eqn{r_j<\texttt{thresh}\max_{k\in S}r_k} or \eqn{r_j<10^{-20}},
where \eqn{S} contains selectable, nonmandatory rows. Very small category
entries in a retained row are not individually hard thresholded.
With \code{grouped=FALSE}, each category column instead uses
\eqn{|b_{jc}|<\texttt{thresh}\max_{k\in S}|b_{kc}|} or
\eqn{|b_{jc}|<10^{-20}}. Surviving entries are then recentered within their
allowed set to preserve zero-sum rows, while removed entries remain zero;
rows with fewer than two allowed entries are zero. Mandatory rows are exempt
and excluded from cutoff maxima. Thresholded entries or rows remain excluded
during selection. Cutoffs use normalized working coordinates; selectable
slopes equal original-unit coefficients times the stored column scales.
\code{normalize="none"} gives scale one. The default cutoff is
\code{thresh=1e-5}. Setting \code{thresh=0} removes only the relative cutoff.
These cutoffs apply only during selection with \code{lambda>0}.

With \code{relaxed=TRUE}, a final joint multinomial IRLS pass uses
\code{lambda=0} and disables hard thresholding. A grouped refit admits all
category entries in retained predictor rows, so previously zero category
entries within those rows can become nonzero. An ungrouped refit admits only
the selected nonzero CLR entries and fixes other category entries at zero.
Both preserve zero-sum rows: \eqn{s} admitted entries give \eqn{s-1} free
directions, and a nonzero row needs at least two nonzero entries. No new
predictor rows can enter. Explicit all-zero starting rows remain excluded;
intercept and \code{always.include} rows remain eligible and unpenalized.
Sparse output retains all original predictor rows and category columns.
Grouping defaults to \code{TRUE}, and relaxed fitting defaults to \code{TRUE}
for grouped multinomial models with either factor or count-matrix responses.
Set \code{relaxed=TRUE} explicitly to refit an ungrouped model. Each selection
or relaxed pass has its own \code{maxit} limit; \code{group.maxit} is unused.

\code{metrics$objective} uses the requested L0 charge. \code{GIC_rows} charges
selected rows. Parameter-count ICs instead use \eqn{C-1} free directions per
grouped row, or \eqn{s_j-1} for a row with \eqn{s_j} active CLR entries, plus
mandatory directions. Grouped IC tuning uses the row criterion.
The algorithm is a subset-selection heuristic and does not certify global
optimality. These loops are serial. Check \code{convergence}, especially
under separation or at small regularization values.
}
\examples{
if (requireNamespace("mgcv", quietly = TRUE)) {
  set.seed(202)
  X <- matrix(rnorm(300 * 6), 300, 6)
  eta <- cbind(0, 0.2 + X[, 1] - 0.6 * X[, 2],
               -0.3 - 0.8 * X[, 1] + 0.7 * X[, 3])
  P <- exp(eta - apply(eta, 1, max))
  P <- P / rowSums(P)
  y <- factor(vapply(seq_len(nrow(X)), function(i)
    sample.int(3, 1, prob = P[i, ]), integer(1)))
  fit_multinomial <- L0glm.fit(X, y, family = mgcv::multinom(K = 2),
    lambda = "gic", grouped = TRUE, relaxed = TRUE,
    maxit = 200, tol = 1e-8)
  fit_multinomial$coefficients
  fit_multinomial$convergence
  max(abs(Matrix::rowSums(fit_multinomial$coefficients)))
  head(predict(fit_multinomial))
  L0glm.metrics(fit_multinomial)
}
}
\references{
Krishnapuram, B., Carin, L., Figueiredo, M. A. T. and Hartemink, A. J. (2005).
Sparse multinomial logistic regression: Fast algorithms and generalization bounds.
IEEE TPAMI 27(6), 957--968.
}
\seealso{\code{\link{L0glm.fit}}, \code{\link{L0glm.tune}},
\code{\link{L0glm.metrics}}, \code{\link[mgcv]{multinom}}}
