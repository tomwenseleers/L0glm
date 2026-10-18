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
    method = args$method, control.family = data$control), class = "L0glm_extended_fit")
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

#' Predict from an extended-family fit
#' @param object A fitted extended-family model.
#' @param newdata Predictor matrix without an automatically added intercept.
#' @param type Prediction type; response, link, class, risk or survival.
#'   Multinomial link predictions are CLR scores for every category.
#' @param offset Optional prediction offset, on each linear predictor.
#'   Multinomial offsets have one column per category and are centered rowwise.
#' @param times Evaluation times for Cox survival probabilities.
#' @param ... Unused.
#' @export
predict.L0glm_extended_fit <- function(object, newdata = NULL, type = "response",
                                       offset = NULL, times = NULL, ...) {
  if (is.null(newdata)) {
    X <- object$X
    if (is.null(offset)) offset <- object$offset
  } else {
    X <- newdata
    if (object$added_intercept) X <- cbind(Intercept = 1, X)
  }
  if (ncol(X) != nrow(object$coefficients)) stop("Wrong number of prediction columns.", call. = FALSE)
  coefficients <- object$coefficients
  eta <- as.matrix(X %*% coefficients)
  eta <- eta + L0glm_ext_matrix(offset, nrow(X), ncol(eta), 0, "offset")
  if (object$kind == "multinomial") eta <- eta - rowMeans(eta)
  if (type == "link") return(eta)
  if (object$kind == "nb") {
    if (type != "response") stop("NB prediction type must be response or link.", call. = FALSE)
    return(matrix(object$family$linkinv(eta), nrow(eta), ncol(eta)))
  }
  if (object$kind == "multinomial") {
    probability <- L0glm_ext_softmax(eta)$p
    colnames(probability) <- object$response$levels
    if (type == "class") return(factor(object$response$levels[max.col(probability)], levels = object$response$levels))
    if (type != "response") stop("Multinomial prediction type must be response, class or link.", call. = FALSE)
    return(probability)
  }
  if (type %in% c("risk", "response")) return(as.numeric(exp(eta)))
  if (type != "survival" || is.null(times)) stop("Supply times for Cox survival prediction.", call. = FALSE)
  if (any(!is.finite(times) | times < 0)) stop("Survival times must be finite and nonnegative.", call. = FALSE)
  index <- findInterval(times, object$baseline$time)
  loghazard <- if (is.null(object$baseline$log_cumhaz)) log(object$baseline$cumhaz) else object$baseline$log_cumhaz
  loghazard <- c(-Inf, loghazard)[index + 1L]
  exp(-exp(outer(as.numeric(eta), loghazard, "+")))
}

#' Log likelihood of an extended-family fit
#' @param object A fitted extended-family model.
#' @param ... Unused.
#' @export
logLik.L0glm_extended_fit <- function(object, ...) {
  value <- L0glm.metrics(object, mu_floor = 0)$loglik
  structure(value, df = object$metrics$df, nobs = object$n_ic, class = "logLik")
}

#' Score an extended-family model
#' @param object A fitted extended-family model.
#' @param newdata New predictor matrix, or NULL for training data.
#' @param y Responses for newdata, or NULL for training responses.
#' @param weights Optional prior weights for the scored observations.
#' @param offset Optional prediction offsets.
#' @param mu_floor Scoring-only NB mean floor; zero requests exact likelihood.
#' @export
L0glm.metrics <- function(object, newdata = NULL, y = NULL, weights = NULL,
                          offset = NULL, mu_floor = 1e-8) {
  training <- is.null(newdata)
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
    mu <- pmax(predict(object, newdata, offset = offset), mu_floor)
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
  beta <- as.matrix(object$coefficients)
  selectable <- setdiff(seq_len(nrow(beta)), object$always.include)
  selected <- sum(rowSums(beta[selectable, , drop = FALSE] != 0) > 0)
  elements <- sum(beta[selectable, , drop = FALSE] != 0)
  mandatory <- length(object$always.include) * (ncol(beta) - (object$kind == "multinomial"))
  nuisance <- sum(object$theta.estimated)
  df <- if (object$kind == "multinomial") {
    (if (object$grouped) selected * (ncol(beta) - 1L) else elements - selected) + mandatory
  } else elements + mandatory + nuisance
  IC <- get.IC(likelihood, max(2, object$n_ic), max(1, object$p), df)
  charge <- if (object$grouped) selected else elements
  data.frame(loglik = likelihood, deviance, mean_deviance = deviance/if (mass > 0) mass else 1,
    df, selected, nonzero = elements, mandatory, nuisance,
    objective = -2 * likelihood + object$init.lambda * charge,
    GIC_rows = -2 * likelihood + preset.lambda("gic", max(1, object$p), max(2, object$n_ic)) * selected,
    IC, row.names = NULL)
}
