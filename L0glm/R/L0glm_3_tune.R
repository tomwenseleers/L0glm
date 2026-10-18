#' @name L0glm.tune
#' @aliases L0glm.tune
#' @title Tunes the regularisation parameter lambda of an L0glm fit
#'
#' @importFrom qlcMatrix colMax colMin
#' 
#' @description
#' Tunes the lambda parameter for the \code{\link{L0glm.fit}} 
#' function. It uses a sequence of penalty parameter values (\code{lambda}) to fit 
#' L0-penalized GLMs, and selects the best lambda by an information criterion
#' or held-out deviance from cross-validation. It can handle multiple output variables (channels),
#' and can produce a plots of the regularization path, true vs. estimated coefficients (if true coefficients are passed, e.g. from a simulation),
#' estimated support size, deviance, log likelihood, mean absolute deviation, 
#' accuracy, precision, recall and the IC of interest, all in function of the level of
#' regularization (L0 penalty lambda). Some plots are only produced if the user 
#' passes the true coefficients \code{coefs.true} (e.g. from a simulation).
#' Each lambda fit is initialized independently. The preset fit is used only
#' to construct the grid; its selected support is not imposed on other fits.
#' A common supplied start or pilot-derived seed can be reused independently
#' at every lambda, avoiding dependence on the order of the grid.
#' The preset is an approximate adaptive-ridge calibration, not a guarantee
#' that the chosen criterion is minimized there. In particular, unrelaxed
#' estimates retain shrinkage and correlated designs need not follow the
#' orthogonal Gaussian calibration.
#' 
#' @param X 
#' A matrix of predictor variables.
#' @param y 
#' A vector or matrix of response variables.
#' @param family 
#' A description of the error distribution and link function to be
#' used in the model. This should be a call to a family function. (See
#' \code{stats::family} for details of family functions.)
#' @param weights 
#' A matrix of observation weights.
#' @param offset 
#' A matrix to be added to the linear predictor.
#' @param dispersion Passed to \code{L0glm.fit}: \code{NULL} automatically
#' initializes and estimates dispersion, a positive value fixes it, and a negative
#' value supplies its absolute value as an initial estimate. A vector can specify
#' separate values for the response channels. Fixed-scale families retain their
#' family dispersion.
#' @param tune A case-insensitive method for model selection:
#' \code{"aic"}, \code{"bic"}, \code{"gic"},
#' \code{"ebic"}, \code{"mbic"}, or \code{"cv"}. Defaults to \code{"gic"}.
#' Cross-validation minimizes mean held-out deviance. This replaces the former
#' \code{IC} argument.
#' See \link{information_criteria} for formulas, guidance and references.
#' @param nfolds Number of cross-validation folds (default 5), used when
#' \code{tune = "cv"}. Must be at least 2 and no larger than the number of
#' positive-weight observation rows.
#' @param foldid Optional integer fold assignment for each row of \code{X}.
#' Labels must run from 1 to the number of folds. When supplied without
#' \code{nfolds}, the number of folds is inferred. All response channels of a
#' row are held out together. Otherwise balanced random folds are generated;
#' use \code{set.seed()} for reproducibility.
#' @param nlambda Number of lambda values to use for tuning (at least 2; default 10).
#' With at least three points and an interior reference, approximately half
#' lie on either side of \code{lambda_reference}, which is included exactly.
#' A two-point grid contains only the displayed endpoints. If the reference
#' equals the lower bound, all points lie at or above it.
#' @param prop_range Numeric in \((0,\,1]\).  
#' Fraction of each log10 distance from \code{lambda_reference} to the full
#' endpoints \code{lambda_min} and \code{lambda_max} retained in the grid.
#' \code{prop_range = 1} keeps the full range; smaller values zoom toward
#' the reference, with each endpoint moving inward on the log scale.
#' @param lambda_min Finite positive lower endpoint of the full-data lambda
#' grid (default 2, the AIC preset). Set a smaller positive value to explore
#' weaker regularization. This is a practical search bound, not a theoretical
#' minimum for predictive optimality or a guarantee of numerical stability.
#' CV training penalties retain their training-weight scaling, so they can
#' be smaller than this full-data bound. A zoomed grid can start above it.
#' @param start \code{NULL} (default) automatically initializes each lambda
#' fit. Numeric scalar \code{1} supplies ones in each fit's normalized working
#' coordinates. Other numeric starts use original coefficient units, with the
#' structural-zero semantics of \code{L0glm.fit}.
#' Use \code{"pilot"} to build a common seed from ten independent unrelaxed
#' fits at log-spaced penalties between \code{lambda_min} and the GIC preset
#' (or one fit when these bounds coincide). The magnitude is the mean absolute
#' coefficient across pilots, with the sign of its mean; exact cancellation
#' uses the sign of the first largest-magnitude pilot coefficient. This retains
#' the union of pilot supports separately in each channel. Mandatory terms
#' remain eligible. No preceding path fit is used to initialize the next fit.
#' @param coefs.true A vector or matrix of true coefficient values, for comparison 
#' with the estimated coefficients.
#' @param printbest A logical value indicating whether to print the best lambda 
#' value according to the selected criterion or CV (default is TRUE).
#' @param plot A logical value indicating whether to produce the plots mentioned
#' in the description. Both panels use \eqn{\log_{10}(\lambda)} on a reversed
#' horizontal axis, with weaker regularization to the right. Coefficient paths
#' have a separate color per predictor and are labeled \eqn{\hat\beta}.
#' Multichannel paths show coefficient row means. The top axes show support
#' sizes at 16 equally spaced log-scale positions, using the nearest fitted
#' lambda at each position. For multichannel fits these count rows with any
#' nonzero coefficient, independently of the plotted row means. Paths and support
#' counts exclude mandatory terms; no additional models are fitted.
#' @param mu_floor Nonnegative prediction floor used only when scoring
#' positive-mean families, such as Poisson (default \code{1e-8}). It retains all
#' observations and makes scores finite when a fitted mean is zero. Set to zero
#' for exact likelihood scoring, which can give infinite criteria.
#' @param ... Other arguments passed to \code{L0glm.fit}.
#'
#' @details For a zero prior mean, the upper endpoint is expanded until an
#' actual fit has no nonzero penalized coefficients, excluding an intercept
#' and \code{always.include} terms. When the preset fit is nonzero, logarithmic
#' bisection refines the nonzero/zero bracket to within 5\%, so a conservative
#' working-model estimate does not leave a long flat tail in the plot.
#' This is an algorithmic endpoint, not a
#' certified bound for the exact L0 likelihood objective. For example, an
#' all-zero Poisson identity-link model with zero offset has infinite deviance
#' when any observed count is positive. The default scoring floor provides
#' numerical comparisons along the path; it does not change the fitted model.
#' With a nonzero prior mean or an unsuccessful endpoint search, a warning is
#' issued and \code{endpoint_zero} is false. A zoomed grid need not reach the
#' full endpoint. Supplied observation weights, offsets and known dispersions
#' are included in scoring; unknown Gaussian scales are profiled separately
#' for each response channel.
#'
#' With \code{tune = "cv"}, the full-data GIC preset defines the reference
#' grid point. If that preset is nonpositive, the grid uses AIC's positive
#' reference value 2 while \code{lambda_preset} still reports GIC. The reference
#' is raised to \code{lambda_min} if necessary; the reported preset is unchanged.
#' The common
#' grid is constructed from full-data fits; all coefficient initialization,
#' dispersion estimation and grouped selection in each fold use training rows
#' only. Folds are fitted sequentially; \code{parallelize} still controls
#' response-channel parallelism within each fit.
#' With \code{start = "pilot"}, each fold constructs its own seed from its
#' training rows, using training-weight-scaled pilot penalties. Full-data
#' coefficients never initialize a CV training fit. A user-supplied numeric
#' start is reused as supplied; the user is responsible for its provenance.
#' Pilot support restriction is a heuristic and can exclude a useful feature
#' missed by every pilot. It does not guarantee monotone supports or global
#' minimization along the path.
#'
#' The fitter's penalty is on a total-likelihood scale. In each fold lambda
#' is multiplied by the fraction of total prior weight in the training rows,
#' so the returned lambda refers to the full-data penalty. Validation deviance
#' is divided by validation prior weight, then averaged across folds with those
#' weights. Its standard error is the square root of the weighted mean squared
#' deviation of fold scores divided by \code{nfolds - 1}; for equal weights this
#' is \code{sd(fold_scores) / sqrt(nfolds)}. Error bars are plus/minus one SE.
#' Quasi-families can be tuned by CV even when likelihood-based ICs are unavailable.
#' The one-standard-error choice \code{lambda_1se} is the largest evaluated
#' lambda with finite mean deviance no greater than the minimum mean deviance
#' plus the SE at that minimum. Both panels mark this value in blue alongside
#' the minimum-deviance choice and the GIC preset. \code{fit_best} continues to
#' correspond to \code{lambda_best}.
#' \code{fit_1se} contains the full-data fit at \code{lambda_1se}, retained
#' from the same path without any additional fits.
#' 
#' @return
#' \code{L0glm.tune} returns a list containing the following elements:
#' \item{\code{tune}}{chosen selection method.}
#' \item{\code{p}}{original number of selectable predictors, retained despite
#' structural starting zeros and pilot support restriction.}
#' \item{\code{start}}{the common numeric start, including any intercept row,
#' or \code{NULL} for independent automatic starts.}
#' \item{\code{pilot}}{pilot penalties, elapsed seconds and row support sizes,
#' or \code{NULL} when no pilot was requested.}
#' \item{\code{fit_args}}{fitting arguments, including the resolved common
#' start, excluding \code{X}, \code{y} and \code{lambda}. For example,
#' \code{do.call(L0glm.fit, c(list(X = X, y = y, lambda = result$lambda_best),
#' result$fit_args))} reproduces \code{fit_best} with the same package version.}
#' \item{\code{preprocessing}}{training preprocessing of the full-data fit.}
#' \item{\code{metrics}}{a data frame containing metrics evaluating the models fit
#' with different regularization parameters; CV adds \code{cv_deviance} and
#' \code{cv_se}.}
#' \item{\code{lambda_best}}{best lambda by the chosen selection method.}
#' \item{\code{lambda_1se}}{largest evaluated lambda within one SE of the
#' minimum CV mean deviance, using the SE at that minimum; \code{NULL} for
#' information-criterion tuning.}
#' \item{\code{lambda_preset}}{preset lambda value that a priori would be expected 
#' to approximately optimise the chosen criterion, or GIC when tuning by CV.}
#' \item{\code{lambda_reference}}{positive grid reference; normally identical
#' to \code{lambda_preset}, but at least \code{lambda_min}.}
#' \item{\code{lambdas}}{a vector containing all the lambdas that were used.}
#' \item{\code{lambda_min}}{lower endpoint of the full numerical path.}
#' \item{\code{lambda_max}}{upper endpoint of the full numerical path.}
#' \item{\code{lambda_nonzero}}{nonzero lower bracket for the refined endpoint,
#' or \code{NA} when the preset is already zero or refinement is unavailable.}
#' \item{\code{endpoint_zero}}{whether the fitted upper endpoint has zero
#' penalized coefficients.}
#' \item{\code{mu_floor}}{prediction floor used for scoring.}
#' \item{\code{coefficients}}{list of coefficient matrices in lambda order.}
#' \item{\code{fit_best}}{the fitted model at \code{lambda_best}.}
#' \item{\code{fit_1se}}{the full-data fitted model at \code{lambda_1se}, or
#' \code{NULL} for information-criterion tuning.}
#' \item{\code{cv}}{CV results, or \code{NULL} for criterion tuning. Contains
#' \code{lambda}, \code{mean}, \code{se}, the lambda-by-fold \code{deviance}
#' matrix, \code{foldid}, \code{nfolds}, \code{fold_weights}, and the training
#' penalties \code{lambda_train}. \code{fit_best}, \code{fit_1se} and \code{coefficients} are
#' full-data fits, not individual fold models. Also contains each fold's
#' \code{starts}, \code{pilots} and training \code{preprocessing}.}
#' 
#' @export
#' @section Additional-family tuning:
#' Family-specific fits are independently initialized at each lambda.
#' Grouped IC tuning uses row support and ungrouped IC tuning uses entry support.
#' CV fits shape, normalization and pilot starts using training observations only.
#' NB validation deviance uses training-estimated shape; multinomial deviance is per weighted trial.
#' The common lambda grid uses full data, so this is not nested validation.
L0glm.tune <- function(X, y,
                       family,
                       weights      = matrix(1, nrow = nrow(as.matrix(y)),
                                             ncol = ncol(as.matrix(y))),
                       offset       = Matrix::Matrix(0, nrow = nrow(X),
                                                     ncol = ncol(as.matrix(y)),
                                                     sparse = TRUE),
                       dispersion   = NULL,
                       tune         = c("gic", "aic", "bic", "ebic", "mbic", "cv"),
                       nlambda      = 10,
                       prop_range   = 1,
                       coefs.true   = NULL,
                       printbest    = TRUE,
                       plot         = TRUE,
                       mu_floor     = 1e-8,
                       nfolds       = 5,
                       foldid       = NULL,
                       lambda_min   = 2,
                       start        = NULL,
                       ...)
{
  tune <- L0glm_ic_name(if (missing(tune)) tune[1L] else tune, allow_cv = TRUE)
  do_cv <- tune == "cv"
  if (!inherits(family, "family")) family <- family()
  if (!is.null(L0glm_ext_kind(family))) return(L0glm_ext_tune(
    X, y, family, if (missing(weights)) NULL else weights,
    if (missing(offset)) NULL else offset, dispersion, tune, nlambda, prop_range,
    coefs.true, printbest, plot, mu_floor,
    if (missing(nfolds) && !is.null(foldid)) NULL else nfolds,
    foldid, lambda_min, start, list(...)))
  if (!is.numeric(nlambda) || length(nlambda) != 1L || !is.finite(nlambda) ||
      nlambda < 2 || nlambda != as.integer(nlambda))
    stop("nlambda must be an integer of at least 2.", call. = FALSE)
  if (!is.numeric(prop_range) || length(prop_range) != 1L ||
      !is.finite(prop_range) || prop_range <= 0 || prop_range > 1)
    stop("prop_range must be in (0, 1].", call. = FALSE)
  if (!is.numeric(lambda_min) || length(lambda_min) != 1L ||
      !is.finite(lambda_min) || lambda_min <= 0)
    stop("lambda_min must be a finite positive number.", call. = FALSE)
  if (!is.numeric(mu_floor) || length(mu_floor) != 1L ||
      !is.finite(mu_floor) || mu_floor < 0)
    stop("mu_floor must be a finite nonnegative number.", call. = FALSE)
  dots <- list(...)
  if ("IC" %in% names(dots))
    stop("IC was renamed to tune in L0glm.tune; use tune = \"gic\" or tune = \"cv\".",
         call. = FALSE)
  pilot_start <- identical(start, "pilot")
  if (is.character(start) && !pilot_start)
    stop('start must be NULL, numeric starting coefficients, or "pilot".', call. = FALSE)
  y <- as.matrix(y)
  n <- nrow(X); m <- ncol(y)
  if (do_cv) {
    weights <- L0glm_tune_matrix(weights, n, m, 1, "weights")
    offset <- L0glm_tune_matrix(offset, n, m, 0, "offset")
    if (any(!is.finite(weights) | weights < 0))
      stop("CV weights must be finite and nonnegative.", call. = FALSE)
    foldid <- L0glm_cv_folds(weights,
      if (missing(nfolds) && !is.null(foldid)) NULL else nfolds, foldid)
  } else if (!is.null(foldid)) stop("foldid is only used with tune = \"cv\".", call. = FALSE)
  intercept <- if (is.null(dots$intercept)) family$family == "binomial" else dots$intercept
  p <- L0glm_predictor_count(X, intercept, dots$always.include)
  lambda_preset <- preset.lambda(IC = if (do_cv) "gic" else tune, p = max(1L, p), n = n)
  lambda_reference <- if (do_cv && (!is.finite(lambda_preset) || lambda_preset <= 0))
    2 else lambda_preset
  if (!is.finite(lambda_reference) || lambda_reference <= 0)
    stop("The chosen criterion gives a nonpositive preset for these dimensions; choose another criterion.",
         call. = FALSE)
  lambda_reference <- max(lambda_reference, lambda_min)
  args <- c(list(X = X, y = y, family = family, weights = weights,
                 offset = offset, dispersion = dispersion,
                 start = if (pilot_start) NULL else start), dots)
  cv_args <- args
  pilot_bounds <- if (pilot_start) c(lambda_min,
    max(lambda_min, preset.lambda("gic", p = max(1L, p), n = n))) else NULL
  initialization <- if (pilot_start) L0glm_pilot_start(args, pilot_bounds) else
    list(start = start, pilot = NULL)
  args$start <- initialization$start
  fit_lambda <- function(lambda) do.call(L0glm.fit, c(args, list(lambda = lambda)))
  preset_fit <- fit_lambda(lambda_reference)
  included <- preset_fit$always.include
  penalized <- setdiff(seq_len(ncol(preset_fit$X)), included)
  is_zero <- function(fit) all(fit$coefficients[penalized, , drop = FALSE] == 0)

  # A working-model threshold only seeds the search. Verify the endpoint by
  # fitting it, using exactly the same initialization and controls as the path.
  lambda_max <- max(2 * lambda_reference, maxlambda(preset_fit$X,
    weights = preset_fit$weights, z = preset_fit$z, dispersion = preset_fit$dispersion,
    always.include = included, coefs = preset_fit$coefficients, family = family,
    nonnegative = isTRUE(dots$nonnegative)))
  endpoint_fit <- fit_lambda(lambda_max)
  zero_prior <- is.null(dots$prior.mean) || all(dots$prior.mean == 0)
  if (zero_prior) for (attempt in seq_len(30L)) {
    if (is_zero(endpoint_fit) || lambda_max > .Machine$double.xmax / 2) break
    lambda_max <- 2 * lambda_max
    endpoint_fit <- fit_lambda(lambda_max)
  }
  # Refine an overly conservative upper estimate without sharing coefficients
  # across fits. Keep both ends verified under the caller's iteration controls.
  lambda_nonzero <- NA_real_
  if (zero_prior && !is_zero(preset_fit) && is_zero(endpoint_fit)) {
    lambda_nonzero <- lambda_reference
    for (attempt in seq_len(16L)) {
      if (log(lambda_max) - log(lambda_nonzero) <= log(1.05)) break
      midpoint <- exp((log(lambda_max) + log(lambda_nonzero)) / 2)
      candidate <- fit_lambda(midpoint)
      if (is_zero(candidate)) {
        lambda_max <- midpoint
        endpoint_fit <- candidate
      } else lambda_nonzero <- midpoint
    }
  }
  endpoint_zero <- is_zero(endpoint_fit)
  if (!endpoint_zero)
    warning("The upper endpoint did not eliminate all penalized coefficients",
            if (!zero_prior) " (the prior mean is nonzero)" else " within the endpoint search",
            ". lambda_max is not a verified zero endpoint.", call. = FALSE)

  # Zoom each full endpoint toward the reference on the log scale.
  g_IC  <- log10(lambda_reference)
  g_min <- g_IC + prop_range * (log10(lambda_min) - g_IC)
  g_max <- g_IC + prop_range * (log10(lambda_max) - g_IC)
  if (nlambda == 2L || g_min == g_IC) {
    lambda_grid <- 10^seq(g_min, g_max, length.out = nlambda)
  } else {
    n_up <- max(2L, ceiling(nlambda / 2))
    n_down <- nlambda - n_up
    lambda_grid <- 10^c(seq(g_min, g_IC, length.out = n_down + 1L)[seq_len(n_down)],
                        seq(g_IC, g_max, length.out = n_up))
  }
  # Keep exact anchors when reusing fits and preserve both full-range endpoints.
  if (nlambda > 2L || g_min == g_IC)
    lambda_grid[which.min(abs(log10(lambda_grid) - g_IC))] <- lambda_reference
  if (prop_range == 1) {
    lambda_grid[1L] <- lambda_min
    lambda_grid[length(lambda_grid)] <- lambda_max
  }
  if (anyDuplicated(lambda_grid))
    stop("prop_range is too small to produce distinct lambda values.", call. = FALSE)
  cv <- if (do_cv) L0glm_cross_validate(cv_args, preset_fit$X, lambda_grid,
                                       foldid, family, mu_floor, pilot_bounds) else NULL
  if (do_cv && !any(is.finite(cv$mean)))
    stop("No lambda has finite cross-validated deviance; check predictions or use a positive mu_floor.",
         call. = FALSE)
  lambda_1se <- NULL
  if (do_cv) {
    cv_best <- which.min(cv$mean)
    within_1se <- is.finite(cv$mean) & cv$mean <= cv$mean[cv_best] + cv$se[cv_best]
    lambda_1se <- max(lambda_grid[within_1se])
  }
  metric_rows <- vector("list", length(lambda_grid))     # each element numeric
  coefficients <- vector("list", length(lambda_grid))
  fit_best <- NULL
  fit_1se <- NULL
  best_score <- Inf
  
  ## ------------------------------------------------------------------ loop --
  for (i in seq_along(lambda_grid)) {
    lambda <- lambda_grid[i]
    
    fit <- if (lambda == lambda_reference) preset_fit else
      if (lambda == lambda_max) endpoint_fit else fit_lambda(lambda)
    if (do_cv && lambda == lambda_1se) fit_1se <- fit
    
    beta_hat <- fit$coefficients                           # p x m (sparse)
    coefficients[[i]] <- beta_hat
    
    ## quality metrics --------------------------------------------------------
    met <- get.metrics(beta = beta_hat,
                       X = fit$X, y = y, family = family,
                       weights = if (is.null(weights)) matrix(1, n, m) else weights,
                       offset = if (is.null(offset)) matrix(0, n, m) else offset,
                       dispersion = dispersion, mu_floor = mu_floor,
                       beta_true = coefs.true, p = p,
                       always.include = included)
    score <- if (do_cv) cv$mean[i] else met[[tune]]
    if (!is.na(score) && score < best_score) {
      best_score <- score
      fit_best <- fit
    }

    ## numeric row *not* data.frame
    metric_rows[[i]] <- c(lambda = unname(lambda), unlist(met))
  }
  
  ## ----------------------------- build tidy data-frame ----------------------
  metrics <- as.data.frame(do.call(rbind, metric_rows),
                           check.names = FALSE,
                           stringsAsFactors = FALSE)
  if (do_cv) {
    metrics$cv_deviance <- cv$mean
    metrics$cv_se <- cv$se
  }

  ## ------------------------------------------------------------- best lambda -----
  selection <- if (do_cv) metrics$cv_deviance else metrics[[tune]]
  best.index <- which.min(selection)
  if (!length(best.index) || is.null(fit_best))
    stop("No model has an available information criterion; check the family and predictions",
         if (mu_floor == 0) "; a positive mu_floor allows scoring zero predicted means" else "",
         ".",
         call. = FALSE)
  if (is.infinite(best_score) && best_score < 0)
    warning("The likelihood is unbounded at zero residual scale; supply a positive dispersion for meaningful tuning.",
            call. = FALSE)
  lambda_best <- metrics$lambda[best.index]
  
  if (printbest)
    message(sprintf("Best lambda by %s = %.3g",
                    if (do_cv) paste0(cv$nfolds, "-fold CV") else toupper(tune), lambda_best))
  
  ## ----------------------------------------------------------------- return
  result <- list(tune         = tune,
       always.include = included,
       p = p,
       start = initialization$start,
       pilot = initialization$pilot,
       fit_args = args[setdiff(names(args), c("X", "y"))],
       preprocessing = fit_best$preprocessing,
       metrics      = metrics,
       lambda_best  = lambda_best, # 
       lambda_1se   = lambda_1se,
       lambda_preset    = lambda_preset, # lambda preset expected to optimise chosen IC
       lambda_reference = lambda_reference,
       lambdas      = lambda_grid,
       lambda_min   = lambda_min,
       lambda_max   = lambda_max,
       lambda_nonzero = lambda_nonzero,
       endpoint_zero = endpoint_zero,
       mu_floor     = mu_floor,
       coefficients = coefficients,
       fit_best     = fit_best,
       fit_1se      = fit_1se,
       cv           = cv)
  if (plot) L0glm_plot_tuning(result, coefs.true)
  result
}

tune.plots <- function(metrics, coefficients, lambdas, IC = "aic", coefs.true = NULL){
  IC <- L0glm_ic_name(IC)
  bestlambda = lambdas[which.min(metrics[,IC])]
  bestlambda = lambdas[which.min(abs(lambdas-bestlambda))] 
  beta_relaxedgroupL0 = as.matrix(coefficients[which(metrics$lambda==bestlambda),]) 
  IC = toupper(IC) # better for plot titles
  # PLOT TRUE VS ESTIMATED COEFFICIENTS 
  if(!is.null(coefs.true)){
    R2 = round(100*summary(lm(as.vector(coefs.true) ~ 0+as.vector(beta_relaxedgroupL0)))$r.squared,2) # Calculate R2
    title = paste("True vs estimated coefficients\nnonnegative L0 group GLM\nR2=", R2)
    plot(as.vector(beta_relaxedgroupL0), as.vector(coefs.true), pch=16, 
         xlab="estimated coefficients nonnegative L0 group GLM",
         ylab="true coefficients", col="steelblue",
         main=title)
  }
  
  # PLOT ESTIMATED SUPPORT SIZE IN FUNCTION OF REGULARISATION
  if(!is.null(coefs.true)){
    title = paste("Estimated support size in function of regularisation\n(red = estimate at best", IC, ", cyan = true support size)")
  } else {
    title = paste("Estimated support size in function of regularisation\n(red = estimate at best", IC, ")")
  }
  
  
  plot(x = log10(lambdas), y = metrics$support_size_estimated,
       ylab="estimated support size",
       main = title,
       pch=16, type="o")
  abline(v = log10(bestlambda), col="red")
  if(!is.null(coefs.true)) abline(h = metrics$support_size_true, col="cyan")
  
  # PLOT DEVIANCE
  title = paste("Theoretical deviance relative to true underlying signal without noise\n(red=estimate at best", IC,", cyan=best deviance)")
  plot(x=log10(lambdas), y=metrics$deviance, pch=16, type="o", 
       ylim=c(0,max(unlist(metrics$deviance))), 
       ylab="deviance",
       main = title) # theoretical deviance relative to true underlying signal without noise
  lambda_bestdeviance = lambdas[which.min(metrics$deviance)]
  abline(v=log10(bestlambda), col="red")
  abline(v=log10(lambda_bestdeviance), col="cyan")
  par(oma=c(0,2,0,0))
  
  # PLOT LOG LIKELIHOOD IN FUNCTION OF REGULARISATION
  title = paste("Log likelihood in function of regularisation\n(red = estimate at best", IC,", cyan = estimate at best log likelihood)")
  plot(x=log10(lambdas), y=metrics$loglik,
       main=title,
       ylab="log likelihood",
       pch=16, type="o") 
  lambda_bestloglik = lambdas[which.max(metrics$loglik)]
  abline(v=log10(bestlambda), col="red")
  abline(v=log10(lambda_bestloglik), col="cyan")
  
  # PLOT MAD
  # TODO plot MD too?
  title = paste("Mean absolute deviation of coefficients vs true coefs\n(red = estimate at best", IC,")")
  plot(x=log10(lambdas), y=metrics$MAD_coefs, 
       main=title,
       ylab="Mean absolute deviation coefficients",
       pch=16, type="o")
  abline(v=log10(bestlambda), col="red")
  
  # PLOT ACCURACY
  title = paste("Accuracy to recover true support\n(red = estimate at best", IC,", cyan = estimate at best accuracy)")
  plot(x=log10(lambdas), y=metrics$accuracy,
       main=title,
       ylab="Accuracy",
       pch=16, type="o") 
  lambda_best_accuracy = lambdas[which.max(metrics$accuracy)] 
  abline(v=log10(bestlambda), col="red")
  abline(v=log10(lambda_best_accuracy), col="cyan")
  
  # PLOT PRECISION
  title = paste("Precision in function of regularisation\n(red = estimate at best", IC, ", blue = estimate at best precision)")
  plot(x=log10(lambdas), y=metrics$precision,
       main=title,
       ylab="Precision",
       pch=16, type="o") 
  lambda_best_precision = lambdas[which.max(metrics$precision)] 
  abline(v=log10(bestlambda), col="red")
  abline(v=log10(lambda_best_precision), col="cyan")
  
  # PLOT RECALL
  # TODO should it be precision and recall separately or F1 score?
  title = paste("Recall in function of regularisation\n(red = estimate at best", IC, ", blue = estimate at best recall)")
  plot(x=log10(lambdas), y=metrics$recall,
       main=title,
       ylab="Recall",
       pch=16, type="o") 
  lambda_best_recall = lambdas[which.max(metrics$recall)] 
  abline(v=log10(bestlambda), col="red")
  abline(v=log10(lambda_best_recall), col="cyan")
  
  # PLOT SPECIFICITY
  title = paste("Specificity in function of regularisation\n(red = estimate at best", IC, ", blue = estimate at best specificity)")
  plot(x=log10(lambdas), y=metrics$specificity,
       main=title,
       ylab="Specificity",
       pch=16, type="o") 
  lambda_best_specficity = lambdas[which.max(metrics$specificity)] 
  abline(v=log10(bestlambda), col="red")
  abline(v=log10(lambda_best_specficity), col="cyan")
  
  
  # R2
  title = paste("R2 value in function of regularisation\n(red = estimate at best", IC,")")
  plot(x=log10(lambdas), y = metrics$R2_coefs,
       main = title, 
       ylab = "R2", 
       pch = 16, type="o")
  abline(v=log10(bestlambda), col="red")
  
  # TP, FP, TN, FN
  color = c("steelblue", "purple")
  title = paste("Effect of regularisation on True and False Positives\n(red = estimate at best", IC, ")")
  plot(x=log10(lambdas), y = metrics$TP,
       main = title,
       ylab = "Value",
       col = color[1],
       pch = 16,
       type = "o")
  lines(log10(unlist(metrics$lambda)), metrics$FP, type="o", pch=16, col=color[2])
  # lines(log10(metrics$lambda), metrics$TN, type="o", pch=16, col=color[1])
  # lines(log10(metrics$lambda), metrics$FN, type="o", pch=16, col=color[4])
  
  legend("topright", legend=c("TP", "FP"), col=color, pch=16, lty=1)
  abline(v=log10(bestlambda), col="red")
  
  # PLOT IC OF INTEREST
  # TODO should it also plot other ICs?
  title = paste("Effect of regularisation on", IC, "\n(red = estimate at best", IC, ")")
  plot(x = log10(lambdas), y = metrics[[paste0(tolower(IC))]],
       main = title,
       ylab = IC, 
       pch = 16,
       type = "o")
  abline(v=log10(bestlambda), col="red")
}

optimize.threshold <- function(thresh, IC = "bic", coefficients, Xmat, y, family){
  coefficients[abs(coefficients) < thresh[1]] = 0
  IC.val <- get.metrics(beta = coefficients, 
                        X = Xmat, y = y, 
                        family = family)[[paste0(IC)]]
  return(IC.val)
}


