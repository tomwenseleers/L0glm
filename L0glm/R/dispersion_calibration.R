#' Dispersion updates and starting values
#'
#' @name dispersion_calibration
#' @md
#' @description
#' [l0glmFit()] estimates ordinary GLM dispersion and negative-binomial
#' shape within the fit. Available information criteria use the returned
#' noise parameter and count each fitted coefficient and estimated noise
#' parameter.
#'
#' @details
#' For an ordinary GLM with variance function \eqn{V} and weights \eqn{w}, the Pearson
#' target for \eqn{\phi} is
#' \deqn{\frac{1}{N-k}\sum_{i:w_i>0}\frac{w_i(y_i-\mu_i)^2}{V(\mu_i)}.}
#' Here \eqn{N} counts observations with positive weight and \eqn{k} counts
#' the fitted coefficient coordinates, including compulsory terms. Negative
#' binomial with variance \eqn{\mu+\mu^2/\theta} solves the corresponding
#' weighted Pearson equation with target \eqn{\bar w(N-k)}. Its left hand
#' side increases with \eqn{\theta}, so the solver uses a bounded scalar root.
#' Fixed dispersion or shape values stay fixed.
#'
#' During adaptive-ridge selection the updates begin at iteration three
#' once \eqn{k < 0.7N}. A step can change the noise parameter by a factor
#' between one half and two. Three successive changes within one percent
#' of the target pause the updates until the support changes. A relaxed
#' pass uses the full Pearson target with positive residual degrees of
#' freedom. These are residual moment estimates. The coefficient step and
#' noise step need not maximize the joint likelihood.
#'
#' The default \code{control.family=list(dispersion.calibration="auto")}
#' enables the following starting values. Set the entry to \code{FALSE}
#' to use the ordinary initializers. \code{TRUE} enables the same eligible
#' pilots as \code{"auto"}. The Pearson updates run with either choice.
#'
#' - Signed, grouped Gaussian identity-link fits with multiple channels and
#'   between \eqn{n+1} and \eqn{10n} selectable predictors use five
#'   adaptive-ridge steps with shared row norms. The last local ridge
#'   smoother supplies a residual scale for each channel. The subsequent
#'   fit uses the BIC preset if \code{lambda} was omitted. The pilot
#'   needs an in-memory unscreened fit, at least five allowed iterations
#'   and automatic starts and dispersions. It leaves explicit group
#'   specifications to the predictor-group interface.
#' - A five-step Poisson fit with the same link supplies starting values for
#'   automatic negative-binomial log-link fits. For direct nonnegative
#'   identity-link fits it supplies only the starting shape. A fixed or
#'   non-default starting shape, supplied coefficients or prior means,
#'   or preclustering uses the ordinary initializer.
#'
#' A constrained negative-binomial unpenalized refit checks conditional
#' likelihood descent when coefficients touch the boundary. If the
#' clipped working solution stalls, the refit tries a diagonally scaled
#' projected score step within the same support.
#'
#' See the [information criteria][information_criteria] for the role of
#' model-space size. The wide Gaussian BIC preset is motivated by multitask
#' examples where each selected row adds many fitted coefficients. Screened
#' fits retain their usual GIC preset.
#'
#' @seealso [l0glmFit()], [l0glmTune()], [model_criteria]
NULL

L0glm_calibration_mode <- function(mode) {
  if (!(isTRUE(mode) || isFALSE(mode) || identical(mode,"auto")))
    stop('control.family$dispersion.calibration must be TRUE, FALSE or "auto".',call.=FALSE)
  mode
}
