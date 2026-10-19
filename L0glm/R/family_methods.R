L0glm_ext_loglik <- function(data, B, state) {
  if (data$kind == "nb") {
    mu <- L0glm_ext_nb_mean(data, B)$mu
    if (any(!is.finite(mu) | mu < 0)) return(-Inf)
    return(sum(vapply(seq_len(data$m), function(j) {
      good <- data$weights[, j] > 0
      sum(data$weights[good, j] * stats::dnbinom(data$y[good, j],
        mu = mu[good, j], size = state$theta[j], log = TRUE))
    }, numeric(1))))
  }
  if (data$kind == "cox") return(L0glm_ext_cox_stats(data, B)$loglik)
  probability <- L0glm_ext_softmax(as.matrix(data$Z %*% B + data$offset))
  constant <- lgamma(rowSums(data$y) + 1) - rowSums(lgamma(data$y + 1))
  sum(data$weights[, 1L] * (constant + rowSums(ifelse(data$y == 0, 0, data$y * probability$logp))))
}

L0glm_ext_result <- function(data, result, passes, lambda, grouped, relaxed, args) {
  beta <- result$B / data$preprocessing$scale
  if (length(data$shifted_intercept)) beta[data$shifted_intercept, ] <-
    beta[data$shifted_intercept, ] - as.numeric(crossprod(data$shift, beta))
  rownames(beta) <- colnames(data$X)
  colnames(beta) <- if (data$kind == "multinomial") data$response$levels else
    if (data$kind == "cox") "linear_predictor" else colnames(data$y)
  convergence <- do.call(rbind, lapply(names(passes), function(stage)
    data.frame(stage, iterations = passes[[stage]]$iterations,
               converged = passes[[stage]]$converged)))
  families <- if (data$kind == "nb") lapply(result$state$theta, function(theta) {
    family <- unserialize(serialize(data$family, NULL))
    family$putTheta(log(theta))
    family
  }) else NULL
  answer <- structure(list(coefficients = Matrix::Matrix(beta, sparse = TRUE),
    family = if (length(families) == 1L) families[[1L]] else data$family,
    families = families, kind = data$kind, theta = result$state$theta,
    theta.estimated = result$state$estimated, dispersion = rep(1, data$m),
    X = data$X, y = data$y, weights = data$weights, offset = data$offset,
    response = data$response, added_intercept = data$added_intercept,
    preprocessing = data$preprocessing, always.include = data$included,
    p = data$p_selectable, n_ic = data$n_ic, init.lambda = lambda,
    grouped = grouped, relaxed = relaxed, convergence = convergence,
    history = lapply(passes, `[[`, "history"), solver = args$solver,
    method = args$method, control.family = data$control),
    class = c("L0glm_extended_fit", "L0glm_fit"))
  if (data$kind == "multinomial") answer$coordinates <- "clr"
  if (data$kind == "cox") {
    answer$cox.cache <- data$cox.cache
    answer$initialization <- data$initialization
    answer$baseline <- L0glm_ext_cox_baseline(answer)
  }
  answer$metrics <- L0glm.metrics(answer, mu_floor = 0)
  answer$objective <- answer$metrics$objective
  answer
}

#' Log likelihood of an extended-family fit
#' @param object A fitted extended-family model.
#' @param ... Passed to the common \code{logLik} method, including
#'   \code{mu_floor = 0} for exact scoring by default.
#' @seealso \link{model_criteria} for likelihood and information-criterion conventions.
#' @export
logLik.L0glm_extended_fit <- function(object, ...) {
  logLik.L0glm_fit(object, ...)
}

#' Essential fit and recovery metrics for an L0glm model
#'
#' @name L0glm.metrics
#' @md
#' @param object A model returned by [L0glm.fit()] or [L0nnpois()], or a
#'   comparison from [sparseglm()].
#' @param newdata New predictor matrix, without an automatically added
#'   intercept, or \code{NULL} for training data.
#' @param y Responses for \code{newdata}; required when it is supplied.
#'   \code{NULL} otherwise uses the training responses. Special-family response
#'   formats are described in [L0glm.fit()].
#' @param weights Optional prior observation weights for scoring. Training
#'   scores default to the original prior weights, rather than IRLS working
#'   weights. New data default to unit weights. Supports the fitting family's
#'   vector or matrix weight formats.
#' @param offset Optional prediction offsets. Training scores use the stored
#'   offsets; new data default to zero. Supports the fitting family's offset
#'   formats.
#' @param mu_floor Nonnegative scoring-only floor for Poisson, NB, Gamma and
#'   inverse-Gaussian means. The default is \code{1e-8}; zero requests exact
#'   scoring. This does not alter the fitted model. \code{logLik()} and the
#'   criterion methods default to zero.
#' @param beta_true Optional true coefficient vector or matrix, on the original
#'   scale and with the same dimensions as \code{coef(object)}, including any
#'   intercept. Multinomial truth must use the displayed CLR coordinates.
#' @param ... Further S3 method arguments; currently unused.
#'
#' @details
#' Row support is defined by at least one exactly nonzero coefficient in a
#' row. Intercepts and other \code{always.include} terms are excluded from
#' support classification. Overall row metrics include TP, FP, TN, FN,
#' precision, recall, F1 and Matthews correlation coefficient (MCC/Phi).
#'
#' For multichannel data, \code{within_row_*} metrics assess coefficient
#' zero/nonzero classification within correctly recovered rows. If
#' \eqn{R=\{j:\|\widehat\beta_{j\cdot}\|_1>0,
#'                 \|\beta^*_{j\cdot}\|_1>0\}}
#' excludes mandatory rows, these counts pool all entries \eqn{(j,k)} with
#' \eqn{j\in R} across channels. They are micro-averaged, rather than an
#' average of per-row scores. Missed and spurious rows remain visible in the
#' overall row metrics and are not included in the conditional scores.
#' \code{within_row_TP}, \code{within_row_FP}, \code{within_row_TN} and
#' \code{within_row_FN} use exact nonzero indicators, and their derived
#' precision, recall, \code{within_row_F1_score} and \code{within_row_MCC}
#' describe recovery of the within-row sparsity pattern. No additional
#' numerical coefficient threshold is applied.
#'
#' Within-row metrics are \code{NA} for a single response or unavailable
#' truth. If no rows are correctly recovered, the multichannel counts are
#' zero and the ratios are \code{NA}. Any ratio with zero denominator is
#' \code{NA}; MCC is \code{NA} if either binary margin is constant.
#'
#' Coefficient accuracy retains mean absolute error (the legacy name
#' \code{MAD_coefs}), \code{RMSE_coefs} and \code{R2_coefs} across all
#' entries, including mandatory terms. For compatibility,
#' \code{R2_coefs} is the through-origin regression R-squared,
#' \eqn{(\widehat\beta^T\beta^*)^2/
#'       (\|\widehat\beta\|_2^2\|\beta^*\|_2^2)}, with matrices vectorized.
#' It measures alignment and can equal one despite multiplicative bias;
#' MAE and RMSE quantify absolute error. Selected-coefficient-only summaries,
#' signed mean error, sMAPE, slope, accuracy, specificity and support-size
#' ratio have been removed.
#'
#' Prediction accuracy includes `R2_weighted` and `adj_R2_weighted` for ordinary
#' GLMs and negative binomial models. Working weights are recomputed at the
#' fitted means as \eqn{w_{ij}=v_{ij}\{d\mu_{ij}/d\eta_{ij}\}^2/V(\mu_{ij})},
#' using original prior weights \eqn{v}, the fitted link, and fitted NB theta.
#' Positive-family means use the fitter's machine-precision lower bound;
#' changing the deviance scoring floor does not change these R-squared scores.
#' These are response-scale residual
#' scores, not working-response or deviance pseudo-R-squared measures.
#'
#' Let \eqn{SSE=\sum_{ij}w_{ij}(y_{ij}-\widehat\mu_{ij})^2} and
#' \eqn{SST=\sum_{ij}w_{ij}(y_{ij}-c_j)^2}. With a unit intercept column,
#' \eqn{c_j} is the channel's working-weighted response mean; without one it is
#' zero, following the centered/uncentered convention of [stats::summary.lm()].
#' Then \eqn{R_w^2=1-SSE/SST} and
#' \eqn{R_{w,adj}^2=1-(SSE/SST)(N-b)/(N-k)}. Here \eqn{N} counts
#' positive-working-weight response entries, \eqn{b} counts channel intercepts,
#' and \eqn{k} counts nonzero coefficients plus mandatory terms even when zero.
#' Multichannel sums pool within-channel residuals and baselines, rather than
#' averaging channel R-squared values. Precision weights do not change \eqn{N}
#' into a sum of weights. Dispersion/theta are not mean-model degrees of freedom.
#' The adjustment uses active-set size, not an estimate of selection or ridge
#' effective degrees of freedom. Negative scores are retained.
#'
#' Prediction R-squared is `NA` for zero total sum of squares, undefined working
#' weights, or Cox/multinomial models, which do not have the scalar independent
#' GLM working weights used by this definition. Adjusted R-squared is also
#' `NA` when residual degrees of freedom are nonpositive. New-data scores use
#' predictions on the new rows and the original fitted active-set size.
#'
#' Likelihood and criterion conventions, including dispersion estimates,
#' sample sizes and grouped parameter counts, are detailed in
#' \link{model_criteria}. Ordinary GLMs also report \code{deviance_true}
#' when truth is supplied: deviance comparing the noiseless true mean with
#' the fitted mean, using the scoring weights and offsets.
#'
#' @return A one-row data frame of fit metrics and information criteria.
#'   Ordinary GLMs include row and coefficient summaries, with \code{NA}
#'   values where truth is unavailable. Special-family scores retain their
#'   model-specific objective summaries; supplying \code{beta_true} also
#'   adds the same recovery metrics. The information criteria are
#'   \code{aic}, \code{bic}, \code{gic}, \code{ebic} and \code{mbic}.
#'   For a \code{sparseglm} collection, one row per requested method
#'   includes timings, status, grouping scope, approximation and tuning values.
#'   Failed methods retain a row with unavailable fit metrics.
#' @seealso \link{model_criteria}, \link{information_criteria}
#' @examples
#' set.seed(12)
#' X <- matrix(runif(300), 60, 5)
#' truth <- rbind(c(3, 0, 2), c(0, 4, 1), matrix(0, 3, 3))
#' y <- matrix(rpois(180, 1 + X %*% truth), 60, 3)
#' fit <- L0nnpois(X, y, offset = rep(1, 60), parallelize = FALSE)
#' scores <- L0glm.metrics(fit, beta_true = truth)
#' scores[c("F1_score", "MCC", "within_row_F1_score", "within_row_MCC")]
#' @export
L0glm.metrics <- function(object, newdata = NULL, y = NULL, weights = NULL,
                          offset = NULL, mu_floor = 1e-8, beta_true = NULL, ...) {
  UseMethod("L0glm.metrics")
}

#' @rdname L0glm.metrics
#' @export
L0glm.metrics.L0glm_fit <- function(object, newdata = NULL, y = NULL, weights = NULL,
                          offset = NULL, mu_floor = 1e-8, beta_true = NULL, ...) {
  if (!inherits(object, "L0glm_fit"))
    stop("object must be a model returned by L0glm.fit or L0nnpois.", call. = FALSE)
  training <- is.null(newdata)
  if (!training && is.null(y))
    stop("Supply y when scoring newdata.", call. = FALSE)
  if (!inherits(object, "L0glm_extended_fit"))
    return(L0glm_metrics_glm(object, newdata, y, weights, offset, mu_floor, beta_true))
  if (is.null(y)) y <- object$y
  n <- if (training) nrow(object$X) else nrow(newdata)
  m <- if (object$kind == "nb") ncol(object$coefficients) else 1L
  if (is.null(weights)) weights <- if (training) object$weights else matrix(1, n, m)
  w <- L0glm_ext_matrix(weights, n, m, 1, "weights")
  if (any(!is.finite(w) | w < 0) || !any(w > 0) ||
      length(mu_floor) != 1L || !is.finite(mu_floor) || mu_floor < 0)
    stop("Invalid scoring weights or mu_floor.", call. = FALSE)
  if (is.null(offset) && training) offset <- object$offset
  if (object$kind == "nb") {
    yy <- as.matrix(y)
    fitted_mean <- predict(object, newdata, offset = offset)
    mu <- pmax(fitted_mean, mu_floor)
    likelihood <- saturated <- 0
    for (j in seq_len(m)) {
      good <- w[, j] > 0
      likelihood <- likelihood + sum(w[good, j] * stats::dnbinom(yy[good, j], mu = mu[good, j], size = object$theta[j], log = TRUE))
      saturated <- saturated + sum(w[good, j] * stats::dnbinom(yy[good, j], mu = yy[good, j], size = object$theta[j], log = TRUE))
    }
    deviance <- 2 * (saturated - likelihood)
    mass <- sum(w)
  } else if (object$kind == "multinomial") {
    yy <- L0glm_ext_multinomial_response(y, object$family, n, object$response$levels)$y
    eta <- predict(object, newdata, type = "link", offset = offset)
    lp <- L0glm_ext_softmax(eta)$logp
    likelihood <- sum(w[, 1L] * (lgamma(rowSums(yy) + 1) - rowSums(lgamma(yy + 1)) + rowSums(ifelse(yy == 0, 0, yy * lp))))
    saturated <- rowSums(ifelse(yy == 0, 0, yy * log(yy/pmax(1, rowSums(yy)))))
    deviance <- 2 * sum(w[, 1L] * (saturated - rowSums(ifelse(yy == 0, 0, yy * lp))))
    mass <- sum(w[, 1L] * rowSums(yy))
  } else {
    same_response <- training && identical(y, object$y)
    response <- if (same_response) object$response else L0glm_ext_cox_response(y, n)
    cache <- if (same_response && identical(as.numeric(w), as.numeric(object$weights))) object$cox.cache else NULL
    eta <- as.numeric(predict(object, newdata, type = "link", offset = offset))
    likelihood <- L0glm_ext_cox_eta(eta, response, w[, 1L], cache)$loglik
    deviance <- -2 * likelihood
    mass <- sum(w[, 1L] * response$status)
  }
  beta <- object$coefficients
  selectable <- setdiff(seq_len(nrow(beta)), object$always.include)
  selected <- sum(Matrix::rowSums(beta[selectable, , drop = FALSE] != 0) > 0)
  elements <- sum(beta[selectable, , drop = FALSE] != 0)
  mandatory <- length(object$always.include) * (ncol(beta) - (object$kind == "multinomial"))
  nuisance <- sum(object$theta.estimated)
  df <- if (object$kind == "multinomial") {
    (if (object$grouped) selected * (ncol(beta) - 1L) else elements - selected) + mandatory
  } else elements + mandatory + nuisance
  IC <- get.IC(likelihood, max(2, object$n_ic), max(1, object$p), df)
  IC <- IC[c("gic", "ebic", "bic", "aic", "mbic")]
  charge <- if (object$grouped) selected else elements
  prediction <- if (object$kind == "nb") L0glm_prediction_metrics(yy, fitted_mean, w,
    if (is.null(object$families)) object$family else object$families, beta,
    object$always.include, L0glm_metric_intercept(object$X, object$always.include)) else
      list(R2_weighted=NA_real_, adj_R2_weighted=NA_real_)
  scores <- data.frame(loglik = likelihood, deviance, mean_deviance = deviance/if (mass > 0) mass else 1,
    df, selected, nonzero = elements, mandatory, nuisance,
    prediction,
    objective = -2 * likelihood + object$init.lambda * charge,
    GIC_rows = -2 * likelihood + preset.lambda("gic", max(1, object$p), max(2, object$n_ic)) * selected,
    IC, row.names = NULL)
  if (!is.null(beta_true)) scores <- cbind(scores,
    as.data.frame(L0glm_recovery_metrics(beta, beta_true, object$always.include, object$p)))
  scores
}

#' @rdname L0glm.metrics
#' @export
L0glm.metrics.default <- function(object, ...) {
  stop("object must be a model returned by L0glm.fit or L0nnpois, or a benchmark result.", call. = FALSE)
}

# Ordinary GLMs retain prior weights separately from the final IRLS weights.
L0glm_metrics_glm <- function(object, newdata, y, weights, offset, mu_floor, beta_true) {
  training <- is.null(newdata)
  X <- if (training) object$X else newdata
  if (!training && object$added_intercept) X <- cbind(Intercept = 1, X)
  if (is.null(y)) y <- object$y
  if (is.null(weights) && training) weights <- object$prior.weights
  if (is.null(offset) && training) offset <- object$offset
  dispersion <- ifelse(object$dispersion.estimated, -object$dispersion, object$dispersion)
  as.data.frame(get.metrics(object$coefficients, X, y, weights = weights,
    offset = offset, family = object$family, beta_true = beta_true,
    p = object$p, dispersion = dispersion, mu_floor = mu_floor,
    always.include = object$always.include))
}
