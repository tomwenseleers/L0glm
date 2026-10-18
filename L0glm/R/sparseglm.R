#' Fit, tune and compare sparse regression methods
#'
#' @name sparseglm
#' @md
#' @description
#' Compare sparse single-response, independent multichannel, joint multinomial
#' and right-censored Cox PH models. Select methods by a common information
#' criterion or observation-row cross-validation, collect elapsed timings and
#' extract a comparison table with [L0glm.metrics()].
#'
#' @param X Design matrix without an intercept column, or a list containing
#'   `X`, `y` and optionally `beta_true`, `family`, `weights` and `offset`.
#' @param y Numeric response vector or independent-channel matrix, multinomial
#'   factor/count matrix, or `survival::Surv(time,status)` for right censoring.
#'   Omit when `X` is a data list.
#' @param family Standard GLM family object/constructor, `mgcv::multinom(K)`
#'   or `mgcv::cox.ph()`. `"multinomial"` and `"cox"` construct the latter
#'   families automatically. `NULL` uses the data list's family or Gaussian.
#'   A standard-family name in a data list uses identity link, except binomial
#'   which uses logit. Estimated-theta negative-binomial comparisons are not
#'   currently supported.
#' @param method Name, vector of names or `"all"`: `"L0glm"`, `"abess"`,
#'   `"L0Learn"`, `"MCP"`, `"SCAD"`, `"lasso"`, `"relaxed_lasso"`,
#'   `"adaptive_lasso"`, `"adaptive_relaxed_lasso"`. Names are case-insensitive;
#'   spaces and hyphens may replace underscores. `"all"` includes only installed
#'   methods compatible with the family, constraints, grouping, weights and
#'   offsets. Explicitly requested unavailable methods report their reason.
#' @param tune `"aic"`, `"bic"`, `"gic"`, `"ebic"`, `"mbic"` or `"cv"`.
#'   Every candidate is rescored with the same package likelihood/IC code.
#' @param nlambda Requested number of path values, default 40. Actual backend
#'   grids are retained and may contain fewer distinct models. The glmnet methods
#'   use an explicit log-spaced grid of this length from their native
#'   `lambda_max` to `lambda_max * lambda.min.ratio`; set the ratio through
#'   method-specific `options`. A supplied `lambda` sequence takes precedence.
#' @param lambda Optional nonnegative sequence, or named list of method-specific
#'   sequences. Penalty units differ between methods. ABESS instead searches
#'   `support.size`; its additional ridge penalty is fixed at zero.
#' @param support.size Optional integer candidate sizes for ABESS. The default
#'   spans zero to `min(ncol(X), nrow(X)-1-intercept)` in up to `nlambda`
#'   steps, without a 100-feature cap. L0Learn uses its maximum as a support
#'   limit; its L0 path need not visit every intermediate size. For strongly
#'   correlated designs, supply consecutive sizes (for example `1:2000`) and
#'   retain ABESS's `warm.start=TRUE`; coarse jumps can give poor local fits.
#' @param nfolds Number of CV folds, default 5. Inferred from `foldid` when only
#'   that argument is supplied.
#' @param foldid Optional row-fold labels shared by every method and channel.
#'   For odd/even Venetian-blind validation use `rep(1:2,length.out=nrow(X))`.
#' @param weights,offset Original prior observation weights and offsets. Default
#'   to data-list components, then unit weights and zero offsets. Independent
#'   channels accept vectors or channel matrices. Multinomial and Cox use one
#'   weight per row; multinomial offsets have one column per category and are
#'   row-centered. Backend restrictions appear in the returned status table.
#' @param dispersion Known dispersion, negative free starting value or `NULL`,
#'   as in [L0glm.fit()], for consistent ordinary-GLM likelihood scoring.
#' @param intercept Add an intercept, default `TRUE` for signed models/binomial
#'   and `FALSE` otherwise. Cox always omits the unidentifiable intercept.
#' @param nonnegative Constrain coefficients, retaining a signed binomial
#'   intercept. Ignored with a warning for multinomial CLR coefficients.
#' @param grouped `"auto"`, `TRUE` or `FALSE`. Auto uses native row-group
#'   selection when available and independent channels otherwise. `TRUE`
#'   requires native grouping; `FALSE` uses independent channels. For
#'   multinomial it chooses grouped versus ungrouped *joint* selection;
#'   marginal binomial fits are never substituted. `scope` records the choice.
#' @param options Named backend argument lists, e.g.
#'   `list(L0glm=list(maxit=50), abess=list(c.max=5))`. Shared data, constraints,
#'   tuning and path arguments cannot be overridden. Native Gaussian LASSO
#'   variants and ABESS require common channel weights for grouped fits; auto
#'   otherwise preserves the full weight matrix using independent paths.
#' @param approximation `"none"` or `"wls"`. For identity-link Poisson, WLS
#'   fits every competing method except L0glm with an explicitly labelled
#'   Gaussian approximation using weights `weights/(y+0.1)`. This also avoids
#'   glmnet's invalid zero-mean Poisson null model without an intercept or
#'   positive offset. Scoring still uses original
#'   Poisson likelihood, weights and offsets. Other families are not approximated.
#' @param adaptive.gamma Adaptive-LASSO exponent, default 1.
#' @param adaptive.epsilon Positive adaptive denominator floor, default `1e-8`.
#' @param seed Random seed; the caller's random-number state is restored.
#' @param mu_floor Scoring-only positive-mean floor, default `1e-8`.
#' @param beta_true Optional original-scale true coefficients, excluding or
#'   including an added intercept row. Multinomial truth has one column per
#'   category and is centered into CLR coordinates.
#' @param plot Plot the completed coefficient comparison, without refitting.
#'
#' @details
#' L0glm uses [L0glm.tune()] for its automatic full-data grid. LASSO variants
#' use `glmnet`; ABESS searches support sizes; L0Learn searches L0 penalties.
#' MCP/SCAD use optional `ordinis` for Gaussian/logit/Poisson-log models and
#' optional `ncvreg` for Cox. Original dimensions and scales are retained.
#' Backends unable to constrain a nonnegative Gaussian/Poisson intercept
#' report this limitation; L0glm and L0Learn support that bound.
#'
#' Adaptive weights come from exactly one
#' `L0glm.fit(start=NULL, lambda="gic", maxit=1, grouped=FALSE,
#' relaxed=FALSE, thresh=0)` call. Independent-channel weights are
#' \eqn{\max(|\widetilde\beta_{jk}|,\epsilon)^{-\gamma}}. Grouped and
#' multinomial glmnet paths use pilot row RMS instead of signed averages;
#' multinomial factors are shared across categories even for ungrouped fits.
#' Intercepts are unpenalized. Unscaled penalty factors and ridge coefficients
#' are retained; glmnet rescales its factors internally. The pilot and entire
#' sparse/refitting path are recomputed using training data in each CV fold.
#' Relaxed variants use a fully relaxed support-restricted refit (`gamma=0`),
#' without tuning a mixture of penalized and unpenalized coefficients.
#'
#' CV minimizes original-family held-out deviance per total validation weight,
#' or per weighted trial for multinomial. Cox uses Breslow partial-likelihood
#' subtraction with the full risk sets, per held-out event, as in [L0glm.tune()].
#' Folds share observation rows across channels, and whole multinomial count
#' rows remain together. Fold means and SE use the conventions of [L0glm.tune()].
#' L0glm sum-loss penalties are scaled by the training likelihood-mass fraction;
#' other backend lambdas retain their native units. The selected model uses all
#' data. CV used to select a model is not an independent test of that model.
#'
#' Multinomial comparisons support L0glm, ABESS and all four glmnet variants.
#' ABESS counts use at most one weighted row per observation/category, without
#' expanding individual trials. Reference coefficients are reconstructed and
#' all outputs are centered to CLR. Multinomial MCP, SCAD and L0Learn are
#' unsupported. Cox supports L0glm, ABESS, glmnet variants and MCP/SCAD via
#' ncvreg. ABESS and ncvreg require untied times and unit weights for a common
#' Breslow comparison. ncvreg also requires zero offsets and unrestricted
#' coefficients. Different tie conventions are not silently
#' mixed. L0Learn has no Cox loss. For glmnet, zero times are replaced by
#' order-preserving positive ranks, leaving Cox risk sets and ties unchanged.
#'
#' `seconds` includes the full-data path, adaptive pilot, scoring and CV;
#' `path_seconds` and `cv_seconds` separate these components. Missing packages,
#' unsupported combinations and failures retain explanatory status rows so
#' remaining methods can finish. L0Learn's six-digit lambda reporting is
#' accounted for when matching training-fold paths to the requested grid.
#'
#' @return A `sparseglm` object with named `fits`, `tuning` candidate paths/CV
#'   diagnostics, `timings`, common `foldid`, input data and optional truth.
#'   `coef()` extracts named sparse matrices; [L0glm.metrics()] returns a table
#'   of metrics and timings, including unavailable-method rows.
#'   `print()` gives a compact overview; `summary()` adds essential fit and
#'   recovery scores while retaining the complete table in `$metrics`.
#'   `predict()` and `fitted()` return named per-method prediction lists;
#'   unavailable methods retain `NULL` entries. See [model_predictions].
#' @seealso [L0glm.metrics()], [L0glm.tune()], [model_summary], [model_predictions]
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(360), 60, 6)
#' truth <- c(2,-1,rep(0,4))
#' y <- as.numeric(X %*% truth + rnorm(60))
#' comparison <- sparseglm(X, y, method=c("lasso","relaxed_lasso"),
#'   tune="gic", nlambda=5, beta_true=truth)
#' L0glm.metrics(comparison)[c("method","seconds","gic","F1_score","MCC")]
#' @export
sparseglm <- function(X, y=NULL, family=NULL, method="all", tune="gic",
  nlambda=40L, lambda=NULL, support.size=NULL, nfolds=5L, foldid=NULL,
  weights=NULL, offset=NULL, dispersion=NULL, intercept=NULL, nonnegative=FALSE,
  grouped="auto", options=list(), approximation="none", adaptive.gamma=1,
  adaptive.epsilon=1e-8, seed=1L, mu_floor=1e-8, beta_true=NULL, plot=FALSE) {
  if (is.list(X) && !is.null(X$X) && is.null(y)) sim <- X else sim <- list(X=X, y=y)
  if (!is.null(beta_true)) sim$beta_true <- beta_true
  result <- L0glm_benchmark(sim, family, method, tune, nlambda, lambda, support.size,
    if (missing(nfolds) && !is.null(foldid)) NULL else nfolds, foldid,
    weights, offset, dispersion, intercept, nonnegative, options, approximation,
    adaptive.gamma, adaptive.epsilon, seed, mu_floor, plot, grouped)
  result$call <- match.call()
  result
}
