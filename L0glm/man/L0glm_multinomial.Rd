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
Starts, prior means and offsets have \eqn{C} columns and are centered across
categories. Scalar \code{start=1}, like \code{start=0}, centers to zero.
A weak-ridge step seeds zero starting coefficients before adaptive weighting.
Starting zeroes are not structural masks for multinomial fitting.
\code{nonnegative=TRUE} is ignored with a warning.
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
where \eqn{M} contains mandatory terms. Default adaptive ridge diagonals are
\eqn{D_{jc}=\lambda/(4\sum_hB_{jh}^2)}.
Positive and negative class effects cannot cancel in this norm.
Both the row penalty and the model are invariant to the choice of orthonormal
contrast basis and to permutations of category labels.

With \code{grouped=FALSE}, the charge counts nonzero displayed CLR entries,
and \eqn{D_{jc}=\lambda/(4B_{jc}^2)}. This is a different sparsity target from
entrywise ILR sparsity. Each active row is solved on the zero-sum plane of its
allowed class entries; a row with fewer than two entries must be zero.
Entrywise thresholding recenters surviving entries within their allowed set.
The relative \code{thresh} and absolute \code{1e-20} cutoff apply during
selection; relaxed refitting disables thresholding and removes the penalty.
A grouped relaxed refit can restore class entries within retained rows.

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
