#' Scalar noise models for nonnegative signals
#'
#' Exact fixed-parameter NB1, common-scale Gamma and censored Gamma models for
#' [L0glm.fit()], [L0glm.tune()], [sis_glm()] and [glm_cluster()]. These adapters
#' supply the likelihood score and positive quadratic curvature to the shared
#' IRLS loop, with its usual selection, line search and relaxed refitting.
#' These adapters are intended for L0glm's score-based IRLS. External fitting
#' engines need to implement their likelihood scores to fit the same models.
#'
#' @param dispersion Known variance-to-mean ratio for `nb1()` (at least one),
#'   or known Gamma dispersion for `Gamma_censored()` (strictly positive).
#' @param scale Known positive common Gamma scale.
#' @param link Link name. Defaults to `"log"`; `"identity"` is useful for
#'   additive nonnegative signals.
#' @param censor Detection limit for Gamma observations. A recorded zero means
#'   the latent intensity is below this limit. Positive observations must be
#'   at least the limit. Use zero for uncensored, strictly positive Gamma data.
#'
#' @details
#' For NB1, with \eqn{\phi>1},
#' \deqn{Y\sim\operatorname{NB}\{a=\mu/(\phi-1),\ p=1/\phi\},
#' E(Y)=\mu,\quad\operatorname{Var}(Y)=\phi\mu.}
#' `dispersion=1` uses the Poisson limit. Its negative log-likelihood score is
#' \eqn{\{\psi(a)-\psi(a+y)+\log\phi\}/(\phi-1)} and its curvature is
#' \eqn{\{\psi_1(a)-\psi_1(a+y)\}/(\phi-1)^2}.
#'
#' `Gamma_scale()` uses \eqn{Y\sim\operatorname{Gamma}(\mu/b,b)} with fixed
#' scale \eqn{b}, giving \eqn{\operatorname{Var}(Y)=b\mu}. For a positive
#' observation its score is \eqn{\{\psi(\mu/b)+\log b-\log y\}/b} and
#' curvature is \eqn{\psi_1(\mu/b)/b^2}. Independent NB1 channels sharing
#' \eqn{\phi}, and Gamma channels sharing \eqn{b}, are closed under sums.
#' This property applies to latent Gamma variables; sums of censored recorded
#' intensities provide an approximation.
#'
#' `Gamma_censored()` uses the standard Gamma mean/dispersion model, with shape
#' \eqn{1/\phi}, scale \eqn{\phi\mu} and variance \eqn{\phi\mu^2}. In either
#' Gamma model a zero contributes \eqn{\log P(Y<d;\mu)} and a positive
#' observation contributes its log density. The native kernel evaluates the
#' incomplete Gamma function and its derivatives. With censoring enabled, a
#' zero mean represents the degenerate zero limit.
#'
#' Each family object contains its fixed noise parameter. Its GLM dispersion
#' slot is one because the score and likelihood include that parameter; pass
#' it to the constructor. A pooled noise parameter can be supplied for
#' multichannel observations. Coefficients model the
#' latent mean. Deviance is twice the likelihood gain over the best separate
#' mean for each observation. `deviance_true` is unavailable for these models
#' because a noiseless mean is not an observed count or censored response.
#'
#' Clustering uses exact likelihood-ratio distances and a Ward hierarchy, with
#' fixed, dynamic or information-criterion cuts. A Bregman hierarchy is available
#' for the regular GLM families documented in [glm_cluster()].
#'
#' @return A family object containing its fixed noise parameter, likelihood
#'   score and native likelihood/deviance functions.
#' @references See [stats::NegBinomial] and [stats::GammaDist] for the
#'   distribution parameterizations.
#' @examples
#' set.seed(6)
#' X <- matrix(rnorm(300 * 12), 300, 12)
#' mu <- exp(0.5 + X[, 1] * 0.5)
#' y <- rnbinom(300, size = mu / 2, mu = mu)
#' fit <- L0glm.fit(X, y, family = nb1(3), intercept = TRUE)
#' logLik(fit)
#' @name signal_families
NULL

#' @rdname signal_families
#' @export
nb1 <- function(dispersion = 1, link = "log") {
  L0glm_scalar_family(5L, dispersion, L0glm_scalar_link(link,substitute(link)), 0, "NB1")
}

#' @rdname signal_families
#' @export
Gamma_scale <- function(scale = 1, link = "log", censor = 0) {
  L0glm_scalar_family(6L, scale, L0glm_scalar_link(link,substitute(link)), censor, "Gamma_scale")
}

#' @rdname signal_families
#' @export
Gamma_censored <- function(dispersion = 1, link = "log", censor = 0.5) {
  L0glm_scalar_family(4L, dispersion, L0glm_scalar_link(link,substitute(link)), censor, "Gamma_censored")
}

L0glm_scalar_link <- function(value, expression) {
  if (is.symbol(expression) && as.character(expression) %in%
      c("identity","log","sqrt","inverse","logit","probit","cloglog","cauchit"))
    as.character(expression) else value
}

L0glm_scalar_family <- function(kind, parameter, link, censor, name) {
  if (!is.numeric(parameter) || length(parameter) != 1L || !is.finite(parameter) ||
      parameter <= 0 || kind == 5L && parameter < 1)
    stop("The fixed noise parameter must be positive (at least one for NB1).", call. = FALSE)
  if (!is.numeric(censor) || length(censor) != 1L || !is.finite(censor) || censor < 0)
    stop("censor must be a finite nonnegative detection limit.", call. = FALSE)
  if (!is.character(link) || length(link) != 1L) stop("Supply a link name.", call. = FALSE)
  lk <- stats::make.link(link)
  terms <- function(y, mu, deviance = FALSE)
    L0glm_scalar_terms_cpp(as.numeric(y), as.numeric(mu), kind, parameter, censor, deviance)
  valid <- function(y) all(is.finite(y) & y >= 0) && if (kind == 5L)
    all(abs(y - round(y)) <= 1e-7) else all(ifelse(y == 0, censor > 0, y >= censor))
  initialize <- substitute({
    if (!VALID(y)) stop("Response outside the scalar family domain.")
    n <- rep.int(1, nobs)
    mustart <- pmax(y + 0.1, 0.1)
  }, list(VALID = valid))
  structure(c(list(family = name, link = link), lk,
    list(variance = if (kind == 4L) function(mu) parameter * mu^2 else function(mu) parameter * mu,
      dev.resids = function(y, mu, wt) {z <- terms(y, mu, TRUE)[, 4L];z[wt == 0] <- 0;wt * z},
      aic = function(y, n, mu, wt, dev) {z <- terms(y, mu)[, 1L];-2 * sum(wt[wt > 0] * z[wt > 0])},
      initialize = initialize, validmu = function(mu) all(is.finite(mu) &
        if (kind == 5L || censor > 0) mu >= 0 else mu > 0),
      valideta = function(eta) all(is.finite(eta)), dispersion = 1,
      scalar = list(kind = kind, parameter = parameter, censor = censor),
      score = function(y, mu) terms(y, mu)[, 2:3, drop = FALSE],
      log_density = function(y, mu) terms(y, mu)[, 1L], validate_response = valid,
      random = function(mu) {
        mu <- as.numeric(mu)
        live <- which(mu>0); y <- numeric(length(mu))
        y[live] <- if(kind==5L) {
          if(parameter==1) stats::rpois(length(live),mu[live]) else
            stats::rnbinom(length(live),mu=mu[live],size=mu[live]/(parameter-1))
        } else if(kind==6L) stats::rgamma(length(live),shape=mu[live]/parameter,scale=parameter) else
          stats::rgamma(length(live),shape=1/parameter,scale=mu[live]*parameter)
        if(kind!=5L && censor>0) y[y<censor] <- 0
        y
      })), class = "family")
}
