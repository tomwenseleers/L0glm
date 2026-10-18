\name{L0glm_negative_binomial}
\alias{L0glm_negative_binomial}
\title{L0-Penalized Negative Binomial Regression}
\description{
Fit NB2 count models with a fixed or estimated shape parameter using
\code{L0glm.fit} and the \code{mgcv::nb} family constructor. Multiple response
channels may share a selected predictor support.
}
\details{
For observation \eqn{i} in channel \eqn{c},
\deqn{E(Y_{ic})=\mu_{ic},\qquad
  \operatorname{Var}(Y_{ic})=\mu_{ic}+\mu_{ic}^2/\theta_c.}
The scale dispersion is one. Shape \eqn{\theta_c} is a separate parameter,
not an additional multiplicative dispersion. With case weights \eqn{w_{ic}},
the log likelihood is the weighted sum of
\deqn{\log\Gamma(y+\theta)-\log\Gamma(\theta)-\log\Gamma(y+1)
 +\theta\log\theta+y\log\mu-(y+\theta)\log(\theta+\mu).}
The mean satisfies \eqn{g(\mu)=X\beta+o}, where \eqn{o} is an offset.
Log, identity and square-root links are supported through the family object;
identity/square-root fits require valid positive means.

\code{mgcv::nb()} estimates shape starting from one;
\code{mgcv::nb(theta=2)} fixes it at two; and
\code{mgcv::nb(theta=-2)} estimates it starting from two. Each response channel
has its own shape. The supplied family object is not mutated.
An intercept must be requested explicitly with \code{intercept=TRUE}.
}
\section{Adaptive ridge and shape estimation}{
At each IRLS iteration, with \eqn{\dot\mu=d\mu/d\eta}, the working quantities are
\deqn{z=\eta-o+(y-\mu)/\dot\mu,\qquad
 W=w\dot\mu^2/(\mu+\mu^2/\theta).}
The existing ridge solver computes
\deqn{(X^T W X+D)\beta^+=X^T Wz+Da,}
where \eqn{a} is \code{prior.mean} in working coefficient coordinates.
For the default \code{gamma=2}, a selectable coefficient has
\eqn{D_j=\lambda/(4\beta_j^2)}. Grouped selection instead shares
\eqn{D_{jc}=\lambda/(4\sum_h\beta_{jh}^2)} across channels.
Mandatory terms, including intercepts, have zero penalty.
The factor one quarter retains the package's adaptive-ridge calibration;
the algorithm heuristically targets
\deqn{-2\ell(B,\theta)+\lambda\,s(B),}
where \eqn{s} counts nonzero coefficients or selected rows according to
\code{grouped}. This is not a certificate of global L0 optimality.

One safeguarded scalar shape update follows each coefficient update.
For fixed fitted means, write \eqn{U=\partial\ell/\partial\theta} and
\eqn{H=\partial^2\ell/\partial\theta^2}. With \eqn{b=\log\theta},
\deqn{U_b=\theta U,\qquad H_b=\theta^2H+\theta U.}
The proposal \eqn{-U_b/H_b} is used when curvature is negative, otherwise an
ascent direction is used. Log steps are capped at one, bounded by
\code{theta.min} and \code{theta.max}, and backtracked against the exact NB
likelihood. This adds linear work per channel without a nested coefficient
fit. Shape convergence uses \code{theta.tol}. Defaults are
\code{control.family=list(theta.min=1e-6,theta.max=1e7,theta.tol=1e-6)}.

The relative \code{thresh} rule and absolute cutoff \code{1e-20} apply during
selection only. Relaxed fitting removes the penalty on selected support and
does not hard threshold. Explicit starting zeroes remain structural zeroes,
except for intercept/mandatory terms. A grouped refit may recover other entries
within retained rows. Returned coefficients are sparse and in original units.
}
\section{Scoring and tuning}{
\code{L0glm.metrics} reports the exact log likelihood, saturated-model
deviance, target objective, row GIC and parameter-count ICs.
Estimated shapes contribute one degree of freedom per channel to the latter.
\code{logLik(fit)} and \code{fit$metrics} use exact likelihoods.
External scoring defaults to \code{mu_floor=1e-8}; use zero for exact scoring.
CV scores held-out deviance using the shape estimated on the training fold,
without re-estimating it on held-out observations. Check all rows of
\code{fit$convergence}; difficult fits may need an explicit larger \code{maxit}.
These family-specific loops are serial; \code{openmp} and \code{ncores} do not
parallelize them.
}
\examples{
if (requireNamespace("mgcv", quietly = TRUE)) {
  set.seed(201)
  X <- matrix(rnorm(300 * 6), 300, 6)
  mu <- exp(0.4 + 0.8 * X[, 1] - 0.6 * X[, 2])
  y <- rnbinom(300, mu = mu, size = 2)
  fit_nb <- L0glm.fit(X, y, family = mgcv::nb(), intercept = TRUE,
    lambda = "gic", relaxed = TRUE, maxit = 200, tol = 1e-8)
  fit_nb$theta
  fit_nb$coefficients
  fit_nb$convergence
  L0glm.metrics(fit_nb, mu_floor = 0)
  head(predict(fit_nb, X[1:5, , drop = FALSE]))
}
}
\seealso{\code{\link{L0glm.fit}}, \code{\link{L0glm.tune}},
\code{\link{L0glm.metrics}}, \code{\link[mgcv]{nb}}}
