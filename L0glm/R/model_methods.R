#' Likelihood, information criteria and coefficients of an L0glm fit
#'
#' @name model_criteria
#' @importFrom stats nobs
#' @md
#' @description
#' Extract the fitted log likelihood, parameter-count information criteria,
#' coefficient matrix or criterion sample size from [L0glm.fit()] and
#' [L0nnpois()] models, including the special families.
#'
#' @param object A model returned by [L0glm.fit()] or [L0nnpois()].
#' @param ... Further arguments. [stats::AIC()] and [stats::BIC()] also accept
#'   several fitted models for comparison.
#' @param mu_floor Nonnegative scoring-only floor for positive-family means.
#'   The default, zero, requests the exact score at the fitted means.
#'
#' @details
#' \code{logLik()} evaluates the unpenalized likelihood at the returned
#' coefficients, using the original prior observation weights and offsets.
#' The adaptive-ridge penalty is not part of this value. No additional
#' likelihood pass is performed when fitting ordinary GLMs or native Poisson
#' models; the score is calculated when requested.
#'
#' The \code{logLik} object's \code{df} attribute is a model parameter count,
#' not the effective degrees of freedom of a ridge smoother. Nonzero
#' coefficients and all mandatory coefficient dimensions are counted, along
#' with one nuisance parameter for each estimated dispersion or NB theta.
#' Known dispersions and fixed NB theta do not add a nuisance parameter.
#' Gaussian estimated dispersion is profiled as weighted RSS divided by the
#' number of positive-weight observations in each channel. Estimated Gamma
#' and inverse-Gaussian dispersion use the same deviance-based
#' \code{family$aic} convention as the existing package metrics, which need
#' not equal the jointly maximized likelihood. Quasi-families return
#' \code{NA} likelihood-based criteria.
#'
#' Multinomial grouped selection counts \eqn{K-1} parameters per active row
#' for \eqn{K} categories. Ungrouped CLR selection counts nonzero entries
#' minus active rows, imposing the zero-sum constraint. Mandatory rows count
#' \eqn{K-1} parameters. Cox uses the Breslow partial log likelihood; its
#' unspecified baseline hazard is not counted in these partial-likelihood
#' criteria.
#'
#' With \eqn{d} the parameter count, the scores are
#' \deqn{\mathrm{AIC}=-2\ell+2d,\qquad
#'       \mathrm{BIC}=-2\ell+d\log n,}
#' \deqn{\mathrm{GIC}=-2\ell+d\log p\log\log n.}
#' \code{EBIC()} and \code{MBIC()} use the penalties documented in
#' \link{information_criteria}. Here \eqn{p} is the original selectable
#' predictor count, excluding mandatory terms and using at least one in
#' logarithms. Ordinary multichannel GLMs count coefficient entries, so
#' \code{GIC()} is distinct from a row-count grouped-selection objective.
#' These scores describe the returned fit; penalized or incompletely
#' converged fits need not maximize likelihood on their selected support.
#'
#' \code{nobs()} reports the sample size used by the criteria. Existing
#' conventions are retained: the number of design rows for ordinary GLMs;
#' positive-weight rows for NB and multinomial; and weighted event count
#' for Cox. Special-family criterion sample sizes are bounded below by two.
#' Independent multichannel responses use the number of observations,
#' rather than observations times channels, in the complexity penalty.
#'
#' The usual [stats::AIC()] and [stats::BIC()] methods obtain their values
#' through \code{logLik()}, including its \code{df} and \code{nobs}
#' attributes. \code{coef()} returns the sparse coefficient matrix on the
#' original scale and with the original coefficient dimensions.
#'
#' @return \code{logLik()} returns an object of class \code{logLik}, with
#'   \code{df} and \code{nobs} attributes. Each criterion returns a numeric
#'   score; \code{coef()} returns a sparse matrix; \code{nobs()} returns the
#'   criterion sample size.
#' @seealso [L0glm.metrics()], \link{model_summary}, \link{information_criteria}, [stats::AIC()]
#' @examples
#' set.seed(7)
#' X <- matrix(rnorm(240), 60, 4)
#' y <- 2 * X[, 1] + rnorm(60)
#' fit <- L0glm.fit(X, y, lambda = "gic", parallelize = FALSE)
#' logLik(fit)
#' c(AIC = AIC(fit), BIC = BIC(fit), GIC = GIC(fit),
#'   EBIC = EBIC(fit), MBIC = MBIC(fit))
#' coef(fit)
NULL

#' @rdname model_criteria
#' @export
logLik.L0glm_fit <- function(object, ..., mu_floor = 0) {
  scores <- L0glm.metrics(object, mu_floor = mu_floor)
  structure(scores$loglik, df = scores$df,
    nobs = nobs.L0glm_fit(object), class = "logLik")
}

#' @rdname model_criteria
#' @export
nobs.L0glm_fit <- function(object, ...) {
  if (inherits(object, "L0glm_extended_fit")) max(2, object$n_ic) else nrow(object$X)
}

#' @rdname model_criteria
#' @export
coef.L0glm_fit <- function(object, ...) object$coefficients

#' @rdname model_criteria
#' @export
GIC <- function(object, ...) UseMethod("GIC")

#' @rdname model_criteria
#' @export
GIC.L0glm_fit <- function(object, ..., mu_floor = 0) {
  L0glm_model_criterion(object, "gic", mu_floor)
}

#' @rdname model_criteria
#' @export
EBIC <- function(object, ...) UseMethod("EBIC")

#' @rdname model_criteria
#' @export
EBIC.L0glm_fit <- function(object, ..., mu_floor = 0) {
  L0glm_model_criterion(object, "ebic", mu_floor)
}

#' @rdname model_criteria
#' @export
MBIC <- function(object, ...) UseMethod("MBIC")

#' @rdname model_criteria
#' @export
MBIC.L0glm_fit <- function(object, ..., mu_floor = 0) {
  L0glm_model_criterion(object, "mbic", mu_floor)
}

L0glm_model_criterion <- function(object, criterion, mu_floor) {
  likelihood <- stats::logLik(object, mu_floor = mu_floor)
  get.IC(as.numeric(likelihood), attr(likelihood, "nobs"), max(1, object$p),
    attr(likelihood, "df"))[[criterion]]
}
