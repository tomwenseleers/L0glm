\name{L0glm_cox}
\alias{L0glm_cox}
\title{L0-Penalized Cox Models with Cached Matrix-Free IRLS}
\description{
Fit right-censored Cox proportional hazards models with Breslow ties,
adaptive L0 selection, an optional relaxed refit, and automatic LLT/PCG routing.
}
\details{
Use \code{family=mgcv::cox.ph()} as family metadata and supply
\code{survival::Surv(time,status)} or a two-column time/status matrix.
Status is zero for censoring and one for an observed event; it belongs in
the response, not in \code{weights}. No intercept is identifiable.
Constant columns must be removed. Delayed entry, strata, time-dependent
covariates, multistate responses and Efron/exact tie methods are not supported.

The hazard is \eqn{h_i(t)=h_0(t)\exp(x_i^T\beta+o_i)}.
For event times \eqn{t}, define the risk set \eqn{R_t=\{r:T_r\ge t\}},
\eqn{d_t=\sum_{i:T_i=t}w_i\delta_i}, and
\eqn{S_t=\sum_{r\in R_t}w_r\exp(\eta_r)}.
The weighted Breslow partial log likelihood is
\deqn{\ell(\beta)=\sum_t\left\{\sum_{i:T_i=t}w_i\delta_i\eta_i
       -d_t\log S_t\right\}.}
Censored observations remain in risk sets. Weights follow the weighted
Breslow case-weight definition, including fractional weights.
}
\section{Quadratic IRLS and matrix-free products}{
Let \eqn{u=\partial\ell/\partial\beta} and
\eqn{H=-\partial^2\ell/\partial\beta\partial\beta^T}.
With ridge diagonal \eqn{D} and center \eqn{a}, a working update solves
\deqn{(H+D)\Delta=u-D(\beta-a),\qquad \beta^+=\beta+\Delta.}
The center \eqn{a} is \code{prior.mean}, supplied in original predictor units
as a scalar, length-\eqn{p} vector or \eqn{p} by \eqn{1} matrix, without an
intercept. Normalization is handled internally; \code{NULL} uses zero centers.
These are targets for log-relative-hazard coefficients, not the baseline hazard.
Nonzero targets cannot reactivate predictors excluded by explicit zero starts;
their rows remain zero in the original sparse coefficient dimensions.
Unpenalized mandatory terms are unaffected by their targets. The final relaxed
refit removes the ridge penalty on retained support.
An exact partial-likelihood line search safeguards the current ridge subproblem.
For selectable terms, default adaptive weights are
\eqn{D_j=\lambda/(4\beta_j^2)}. Mandatory terms have zero ridge penalty.
The heuristic target is \eqn{-2\ell(\beta)+\lambda\|\beta\|_0}, excluding
mandatory terms from the L0 count. This does not certify a global optimum.

For \eqn{e_i=w_i\exp(\eta_i)} and
\eqn{q_i=e_i\sum_{t\le T_i}d_t/S_t}, the observation information satisfies
\deqn{(Wv)_i=q_iv_i-e_i\sum_{t\le T_i}{d_t\over S_t^2}
                       \sum_{r\in R_t}e_rv_r,\qquad H=X^TWX.}
After sorting times once, two cumulative sums evaluate this expression.
Thus a product \eqn{(H+D)v} costs a design multiplication, its transpose,
and linear risk-set work, without forming an observation or coefficient Hessian.
Time ordering, tie groups, event weights and risk membership are cached once
per fit; fitted likelihood, baseline hazard and training scores reuse the cache.
The public data and predictions retain their original row order.

PCG uses diagonal preconditioning
\eqn{M_{jj}=\sum_iq_ix_{ij}^2+D_j}, an inexpensive positive curvature bound.
Before accepting convergence it recomputes the actual residual and requires
\deqn{\|M^{-1/2}(A\Delta-b)\|_2/
       \|M^{-1/2}b\|_2\le\texttt{cox.cg.tol},\qquad A=H+D.}
For small direct systems, exact information is assembled using
\deqn{H=X^T\operatorname{diag}(q)X-\sum_t d_t m_tm_t^T,\qquad
  m_t=\sum_{r\in R_t}e_rx_r/S_t,}
and solved with Cholesky LLT.
}
\section{Initialization and solver controls}{
Automatic initialization retains the weak-ridge pilot
\deqn{(X^T\operatorname{diag}(q_0)X+\kappa I)\beta_0=X^Ts_0,
       \qquad \kappa=10^{-3},}
where the score and weights are evaluated at zero coefficients, including offsets.
The pilot uses matrix-free PCG too. Wide designs solve its equivalent
observation-space Woodbury system without forming that matrix.
The resulting step is safeguarded against the pilot's exact penalized likelihood.
\code{fit$initialization} records its space, iterations and scaled residual.
Explicit starts bypass this pilot. Zero starts lock selectable coefficients
out; intercept/mandatory eligibility conventions are retained.

\code{control.family} accepts \code{cox="auto"}, \code{"direct"}, or
\code{"pcg"}; \code{cox.llt.max=c(dense=1024L,sparse=0L)};
\code{cox.cg.tol=1e-10}; and \code{cox.cg.maxit=2000L}.
Automatic routing uses LLT for dense active systems up to the threshold when
the active count does not exceed the observation count. Sparse and wide
systems use PCG. A scalar \code{cox.llt.max} overrides both storage thresholds.
The default reflects measured systems, not a universal timing crossover.
Automatic LLT also limits the risk-mean temporary to ten million elements.
\code{max.working=2000L} guards direct working dimensions; matrix-free solves
are not restricted by it.

PCG requires unrestricted coefficients and \code{solver="eigen"}.
Automatic routing retains the direct path for constrained/non-Eigen solvers.
Eigen nonnegative solves use the package's clipping approximation; use OSQP
for formal bounds. A failed automatic PCG solve warns and falls back to a
permitted small direct system. Forced PCG or a large system fails explicitly;
it never silently accepts an inaccurate inner solve. Increase
\code{cox.cg.maxit} for difficult near-interpolating designs.
\code{fit$history} records each backend, inner iteration count and residual.
These loops are serial. Check \code{fit$convergence} for both fitting passes.
}
\section{Prediction, scoring and cross-validation}{
Breslow increments estimate \eqn{d\widehat H_0(t)=d_t/S_t}.
\code{predict(type="survival",times=...)} returns
\eqn{\exp\{-\widehat H_0(t)\exp(\eta)\}} using log cumulative hazards.
\code{type="link"} returns \eqn{\eta}; \code{type="risk"} and the default
\code{type="response"} return \eqn{\exp(\eta)}.
New-data offsets must be supplied explicitly when needed.
\code{metrics$deviance} is \eqn{-2\ell}, not saturated-model GLM deviance.
Presets use weighted event count as sample size and selectable predictor count
as dimension. Cox CV evaluates
\eqn{-2\{\ell_{\rm full}(\widehat\beta_{-k})-
\ell_{-k}(\widehat\beta_{-k})\}} per held-out event weight.
This retains full risk sets while estimating coefficients on training data only.
\code{fit_best} and \code{fit_1se} are full-data refits.
}
\examples{
if (requireNamespace("mgcv", quietly = TRUE) &&
    requireNamespace("survival", quietly = TRUE) &&
    requireNamespace("glmnet", quietly = TRUE)) {
  # Same right-censored example data as the glmnet Coxnet vignette.
  data("CoxExample", package = "glmnet")
  X <- CoxExample$x
  y <- survival::Surv(CoxExample$y[, "time"], CoxExample$y[, "status"])
  set.seed(203)
  train <- sample(seq_len(nrow(X)), floor(0.7 * nrow(X)))
  Xtr <- X[train, , drop = FALSE]
  ytr <- y[train, , drop = FALSE]
  foldid <- sample(rep(1:3, length.out = length(train)))
  time_l0 <- system.time(tune_cox <- L0glm.tune(Xtr, ytr,
    family = mgcv::cox.ph(), tune = "cv", foldid = foldid,
    nlambda = 6, plot = FALSE, printbest = FALSE, relaxed = TRUE,
    maxit = 200, tol = 1e-8))[["elapsed"]]
  fit_cox <- tune_cox$fit_best
  # Recent glmnet releases offer an explicit tie-method option.
  ties <- if ("cox.ties" \%in\% names(formals(glmnet::glmnet)))
    list(cox.ties = "breslow") else list()
  time_lasso <- system.time(lasso_cox <- do.call(glmnet::cv.glmnet,
    c(list(x = Xtr, y = ytr, family = "cox", alpha = 1,
           foldid = foldid, grouped = TRUE), ties)))[["elapsed"]]
  eta_l0 <- as.numeric(predict(fit_cox, X[-train, , drop = FALSE], type = "link"))
  eta_lasso <- as.numeric(predict(lasso_cox, X[-train, , drop = FALSE],
                                   s = "lambda.min", type = "link"))
  ytest <- y[-train, , drop = FALSE]
  concordance <- function(eta) survival::concordance(ytest ~ eta, reverse = TRUE)$concordance
  comparison_cox <- data.frame(method = c("L0 CV", "LASSO CV"),
    selected = c(fit_cox$metrics$selected,
                 sum(as.numeric(coef(lasso_cox, s = "lambda.min")) != 0)),
    heldout_C = c(concordance(eta_l0), concordance(eta_lasso)),
    seconds = c(time_l0, time_lasso))
  print(comparison_cox)
  # Penalty values differ in definition and scale between L0 and LASSO.
  fit_cox$convergence
  survival_first_three <- predict(fit_cox, X[1:3, , drop = FALSE],
    type = "survival", times = c(1, 2, 4))
  survival_first_three
}
}
\references{
Simon, N., Friedman, J., Hastie, T. and Tibshirani, R. (2011).
Regularization paths for Cox's proportional hazards model via coordinate descent.
Journal of Statistical Software 39(5), 1--13.
The comparison example uses data from the glmnet
\href{https://glmnet.stanford.edu/articles/Coxnet.html}{Coxnet vignette}.

Verweij, P. J. M. and van Houwelingen, H. C. (1993).
Cross-validation in survival analysis. Statistics in Medicine 12, 2305--2314.
}
\seealso{\code{\link{L0glm.fit}}, \code{\link{L0glm.tune}},
\code{\link{L0glm.metrics}}, \code{\link[mgcv]{cox.ph}}}
