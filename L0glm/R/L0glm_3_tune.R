#' @name L0glm.tune
#' @aliases L0glm.tune
#' @title Tunes the regularisation parameter lambda of an L0glm fit.
#'
#' @importFrom qlcMatrix colMax colMin
#' 
#' @description
#' Tunes the lambda parameter for the \code{\link{L0glm.fit}} 
#' function. It uses a sequence of penalty parameter values (\code{lambda}) to fit 
#' L0-penalized GLMs, and selects the best  lambda value according to a specified 
#' information criterion. It can handle multiple output variables (channels), 
#' and can produce a plots of the regularization path, true vs. estimated coefficients (if true coefficients are passed, e.g. from a simulation),
#' estimated support size, deviance, log likelihood, mean absolute deviation, 
#' accuracy, precision, recall and the IC of interest, all in function of the level of
#' regularization (L0 penalty lambda). Some plots are only produced if the user 
#' passes the true coefficients \code{coefs.true} (e.g. from a simulation).
#' Each lambda fit is initialized independently. The preset fit is used only
#' to construct the grid; its selected support is not imposed on other fits.
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
#' @param IC A case-insensitive information-criterion name for model selection:
#' \code{"aic"}, \code{"bic"}, \code{"gic"}, \code{"mebic"},
#' \code{"ebic"}, or \code{"mbic"}. Defaults to \code{"gic"}.
#' See \link{information_criteria} for formulas, guidance and references.
#' @param nlambda Number of lambda values to use for tuning (at least 2; default 10).
#' If \code{nlambda} is even, half of the values lie above 
#' the lambda value that would be expected to optimise
#' the chosen IC (\code{lambda_preset}) and half below;  
#' if odd, the upper branch receives the extra point so that
#' \code{lambda_preset} is included exactly once.
#' @param prop_range Numeric in \((0,\,1]\).  
#' Proportion of the log10 distance between \code{lambda_preset} and
#' \code{lambda_max} (the verified numerical endpoint) that is retained on each side of
#' \code{lambda_preset}. \code{prop_range = 1} (default) keeps the 
#' entire range, whilst smaller values "zoom-in", yielding grids 
#' that are denser around \code{lambda_preset}. For example, 
#' \code{prop_range = 0.1} preserves only 10% of the original 
#' log10 span above and below the preset lambda value that 
#' a priori would be expected to optimise the chosen IC.
#' @param coefs.true A vector or matrix of true coefficient values, for comparison 
#' with the estimated coefficients.
#' @param printbest A logical value indicating whether to print the best lambda 
#' value according to the selected IC (default is TRUE).
#' @param plot A logical value indicating whether to produce the plots mentioned
#' in the description, all in function of regularization (lambda).
#' @param mu_floor Nonnegative prediction floor used only when scoring
#' positive-mean families, such as Poisson (default \code{1e-8}). It retains all
#' observations and makes scores finite when a fitted mean is zero. Set to zero
#' for exact likelihood scoring, which can give infinite criteria.
#' @param ... Other arguments to pass to the `L0glm.fit` function except for \code{start}.
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
#' @return
#' \code{L0glm.tune} returns a list containing the following elements:
#' \item{\code{IC}}{chosen \code{IC}.}
#' \item{\code{metrics}}{a data frame containing metrics evaluating the models fit
#' with different regularization parameters.}
#' \item{\code{lambda_best}}{best lambda according to the \code{IC} 
#' criterion specified.}
#' \item{\code{lambda_preset}}{preset lambda value that a priori would be expected 
#' to approximately optimise the chosen \code{IC}.}
#' \item{\code{lambdas}}{a vector containing all the lambdas that were used.}
#' \item{\code{lambda_max}}{upper endpoint of the full numerical path.}
#' \item{\code{lambda_nonzero}}{nonzero lower bracket for the refined endpoint,
#' or \code{NA} when the preset is already zero or refinement is unavailable.}
#' \item{\code{endpoint_zero}}{whether the fitted upper endpoint has zero
#' penalized coefficients.}
#' \item{\code{mu_floor}}{prediction floor used for scoring.}
#' \item{\code{coefficients}}{list of coefficient matrices in lambda order.}
#' \item{\code{fit_best}}{the fitted model at \code{lambda_best}.}
#' 
#' @export
L0glm.tune <- function(X, y,
                       family,
                       weights      = matrix(1, nrow = nrow(as.matrix(y)),
                                             ncol = ncol(as.matrix(y))),
                       offset       = Matrix::Matrix(0, nrow = nrow(X),
                                                     ncol = ncol(as.matrix(y)),
                                                     sparse = TRUE),
                       dispersion   = NULL,
                       IC           = c("gic", "mebic", "aic", "bic", "ebic", "mbic"),
                       nlambda      = 10,
                       prop_range   = 1,
                       coefs.true   = NULL,
                       printbest    = TRUE,
                       plot         = TRUE,
                       mu_floor     = 1e-8,
                       ...)
{
  IC <- L0glm_ic_name(if (missing(IC)) IC[1L] else IC)
  if (!inherits(family, "family")) family <- family()
  if (!is.numeric(nlambda) || length(nlambda) != 1L || !is.finite(nlambda) ||
      nlambda < 2 || nlambda != as.integer(nlambda))
    stop("nlambda must be an integer of at least 2.", call. = FALSE)
  if (!is.numeric(prop_range) || length(prop_range) != 1L ||
      !is.finite(prop_range) || prop_range <= 0 || prop_range > 1)
    stop("prop_range must be in (0, 1].", call. = FALSE)
  if (!is.numeric(mu_floor) || length(mu_floor) != 1L ||
      !is.finite(mu_floor) || mu_floor < 0)
    stop("mu_floor must be a finite nonnegative number.", call. = FALSE)
  dots <- list(...)
  if ("start" %in% names(dots))
    stop("L0glm.tune initializes each fit independently; do not supply start.", call. = FALSE)
  y <- as.matrix(y)
  n <- nrow(X); m <- ncol(y)
  intercept <- if (is.null(dots$intercept)) family$family == "binomial" else dots$intercept
  p <- ncol(X) + as.integer(intercept)
  lambda_preset <- preset.lambda(IC = IC, p = p, n = n)
  if (!is.finite(lambda_preset) || lambda_preset <= 0)
    stop("The chosen criterion gives a nonpositive preset for these dimensions; choose another criterion.",
         call. = FALSE)
  args <- c(list(X = X, y = y, family = family, weights = weights,
                 offset = offset, dispersion = dispersion), dots)
  fit_lambda <- function(lambda) do.call(L0glm.fit, c(args, list(lambda = lambda)))
  preset_fit <- fit_lambda(lambda_preset)
  included <- dots$always.include
  if (intercept) included <- c(1L, included + 1L)
  penalized <- setdiff(seq_len(p), included)
  is_zero <- function(fit) all(fit$coefficients[penalized, , drop = FALSE] == 0)

  # A working-model threshold only seeds the search. Verify the endpoint by
  # fitting it, using exactly the same initialization and controls as the path.
  lambda_max <- max(2 * lambda_preset, maxlambda(preset_fit$X,
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
    lambda_nonzero <- lambda_preset
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

  n_up   <- max(2L, ceiling(nlambda / 2)) # includes lambda_preset
  n_down <- nlambda - n_up              # below lambda_preset
  # ---- log10 limits, shrunk by prop_range 
  g_IC  <- log10(lambda_preset)
  g_max <- log10(lambda_max)
  g_up_max <- g_IC + prop_range * (g_max - g_IC)   # new "effective" lambda_max
  log_step <- (g_up_max - g_IC) / (n_up - 1)
  log_up   <- g_IC + log_step * seq(0, n_up - 1)
  log_down <- g_IC - log_step * seq_len(n_down)
  lambda_grid <- sort(10^c(log_up, log_down))
  
  # Avoid rounding the preset or endpoint before deciding which fits to reuse.
  lambda_grid[which.min(abs(log10(lambda_grid) - g_IC))] <- lambda_preset
  if (prop_range == 1) lambda_grid[length(lambda_grid)] <- lambda_max
  coefmat  <- Matrix::Matrix(0, nrow = length(lambda_grid), ncol = p,
                             sparse = TRUE)         # row-means of beta(lambda)
  metric_rows <- vector("list", length(lambda_grid))     # each element numeric
  coefficients <- vector("list", length(lambda_grid))
  fit_best <- NULL
  best_score <- Inf
  
  ## ------------------------------------------------------------------ loop --
  for (i in seq_along(lambda_grid)) {
    lambda <- lambda_grid[i]
    
    fit <- if (lambda == lambda_preset) preset_fit else
      if (lambda == lambda_max) endpoint_fit else fit_lambda(lambda)
    
    beta_hat <- fit$coefficients                           # p x m (sparse)
    coefficients[[i]] <- beta_hat
    
    ## store coefficient path (row mean for multichannel)
    coefmat[i, ] <- if (m > 1) Matrix(Matrix::rowMeans(beta_hat), sparse=TRUE) else beta_hat
    
    ## quality metrics --------------------------------------------------------
    met <- get.metrics(beta = beta_hat,
                       X = fit$X, y = y, family = family,
                       weights = if (is.null(weights)) matrix(1, n, m) else weights,
                       offset = if (is.null(offset)) matrix(0, n, m) else offset,
                       dispersion = dispersion, mu_floor = mu_floor,
                       beta_true = coefs.true)
    if (!is.na(met[[IC]]) && met[[IC]] < best_score) {
      best_score <- met[[IC]]
      fit_best <- fit
    }

    ## numeric row *not* data.frame
    metric_rows[[i]] <- c(lambda = unname(lambda), unlist(met))
  }
  
  ## ----------------------------- build tidy data-frame ----------------------
  metrics <- as.data.frame(do.call(rbind, metric_rows),
                           check.names = FALSE,
                           stringsAsFactors = FALSE)

  ## ------------------------------------------------------------- best lambda -----
  best.index  <- which.min(metrics[[IC]])
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
    message(sprintf("Best lambda by %s = %.3g", toupper(IC), lambda_best))
  
  ## ---------------------------------------------------------------- plotting
  if (plot) {
    old_par <- par(no.readonly = TRUE); on.exit(par(old_par))
    layout(matrix(1:2, nrow = 1))          # side-by-side
    ## 1. coefficient paths
    matplot(x  = log10(metrics$lambda),
            y  = as.matrix(coefmat),
            type = "l", lty = 1, col = "grey40",
            xlab = expression(log[10](lambda)),
            ylab = "Coefficient estimate",
            main = "Regularisation path")
    abline(v = log10(lambda_max), col = "black", lty = 3)
    abline(v = log10(lambda_best), col = "red", lty = 2)
    abline(v = log10(lambda_preset), col = "red", lty = 3)
    if (!is.null(coefs.true)) {
      if (m > 1) coefs.true <- Matrix::rowMeans(coefs.true)
      abline(h = as.vector(coefs.true[coefs.true!=0]),
             col = scales::alpha("blue", 0.1)) 
    }
    
    ## 2. chosen information criterion
    finite <- is.finite(metrics[[IC]])
    plot(log10(metrics$lambda), replace(metrics[[IC]], !finite, NA_real_),
         type = "b", pch = 19,
         log = if (any(finite) && all(metrics[[IC]][finite] > 0)) "y" else "",
         ylim = if (any(finite)) NULL else c(0, 1),
         xlab = expression(log[10](lambda)), ylab = toupper(IC),
         main = paste("Tuning -", toupper(IC)))
    abline(v = log10(lambda_best), col = "red", lty = 2)
    abline(v = log10(lambda_preset), col = "red", lty = 3)
    if (any(!finite)) mtext("Non-finite criterion values omitted", side = 3, cex = 0.75)
  }
  
  ## ----------------------------------------------------------------- return
  list(IC           = IC, # chosen IC
       metrics      = metrics,
       lambda_best  = lambda_best, # 
       lambda_preset    = lambda_preset, # lambda preset expected to optimise chosen IC
       lambdas      = lambda_grid,
       lambda_max   = lambda_max,
       lambda_nonzero = lambda_nonzero,
       endpoint_zero = endpoint_zero,
       mu_floor     = mu_floor,
       coefficients = coefficients,
       fit_best     = fit_best)
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


