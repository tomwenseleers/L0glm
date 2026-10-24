#' Predictions and fitted values from sparse regression models
#'
#' @name model_predictions
#' @importFrom stats fitted
#' @md
#' @description
#' S3 prediction methods for [L0glm.fit()], [L0nnpois()] and [sparseglm()],
#' including negative-binomial, multinomial and right-censored Cox models.
#'
#' @param object A fitted model or `sparseglm` comparison.
#' @param newdata Predictor matrix, sparse or dense, or numeric data frame.
#'   Omit the intercept added by `intercept=TRUE`. Supply all original predictor
#'   columns, including columns excluded by starting zeroes. Complete unique
#'   column names are matched to the training names; unnamed data use training
#'   column order. A numeric vector denotes one row, or multiple observations
#'   for a one-predictor model. `NULL` uses the stored training design.
#' @param type Default `"response"`: GLM/NB means, multinomial probabilities,
#'   or Cox relative risks. `"link"` returns linear predictors (row-centered CLR
#'   scores for multinomial). `"probabilities"` is also supported for binomial,
#'   quasibinomial and multinomial. Multinomial `"class"` returns a factor;
#'   exact ties choose the first stored category. Cox additionally supports
#'   `"risk"` and `"survival"`; the latter requires `times`.
#' @param offset Prediction offsets: a scalar, vector with one value per
#'   predicted row, or matrix with that many rows and one or one-per-response columns.
#'   Training predictions default to the stored offsets. New-data predictions
#'   default to zero; supply new offsets explicitly, even if the new rows happen
#'   to equal the training rows. Multinomial offsets are on the category logits;
#'   channel/category columns follow the stored response order. Row-constant
#'   shifts have no effect on multinomial probabilities or CLR scores.
#' @param times Finite nonnegative evaluation times for Cox survival prediction.
#' @param ... For `fitted()`, prediction arguments such as `type`, `offset` and
#'   `times`. Use `predict()` for new data. Further `predict()` arguments are
#'   currently unused.
#'
#' @details
#' Predictions use coefficients on their original scale:
#' \deqn{\eta_{ik}=x_i^T\widehat\beta_k+o_{ik},\qquad
#'       \widehat\mu_{ik}=g^{-1}(\eta_{ik}).}
#' Internal centering/normalization is already undone in the returned
#' coefficients and is not reapplied to new data. Designs retain their dense or
#' sparse storage during multiplication. No model is refitted and no numerical
#' scoring floor is applied to predictions. Binomial outputs are success
#' probabilities, not trial counts.
#'
#' Multinomial probabilities use stable softmax evaluation,
#' \deqn{\widehat p_{ik}=\frac{\exp(\eta_{ik})}{\sum_l\exp(\eta_{il})},}
#' in the stored category order. Count-matrix fits still return probabilities
#' summing to one, rather than expected counts at the training trial totals.
#'
#' Cox risk is \eqn{\exp(\eta_i)}, relative to a zero linear predictor rather
#' than centered at the training mean. Survival predictions are
#' \eqn{\widehat S(t|x_i)=\exp\{-\widehat H_0(t)\exp(\eta_i)\}} using
#' the training Breslow baseline. A comparison model's baseline is calculated
#' on demand if it has not been stored. Log cumulative hazards keep this
#' calculation stable under large common offset shifts.
#'
#' `fitted()` evaluates the training rows with stored offsets, and accepts the
#' same prediction types. Collection methods apply these rules to each selected
#' full-data model; these are not out-of-fold predictions.
#'
#' @return GLM/NB means, binomial probabilities, multinomial probabilities and
#'   all link predictions are matrices, including for one response or one row.
#'   Cox risks are numeric vectors; Cox survival is an observations-by-times
#'   matrix. Multinomial classes are factors with the stored category levels.
#'   `sparseglm` methods return named lists in method order, retaining `NULL`
#'   entries for failed or unavailable methods.
#' @seealso [L0glm.fit()], [L0glm.metrics()], [sparseglm()], [model_summary]
#' @examples
#' set.seed(5)
#' X <- matrix(rnorm(240),60,4)
#' y <- 2*X[,1] + rnorm(60)
#' fit <- L0glm.fit(X,y,lambda="gic",parallelize=FALSE)
#' head(fitted(fit))
#' predict(fit,X[1:3,,drop=FALSE])
#' comparison <- sparseglm(X,y,method=c("lasso","relaxed_lasso"),nlambda=4)
#' predict(comparison,X[1:3,,drop=FALSE])
#' if (requireNamespace("mgcv",quietly=TRUE)) {
#'   labels <- factor(rep(c("A","B","C"),20))
#'   multi <- L0glm.fit(X,labels,family=mgcv::multinom(K=2),lambda="gic")
#'   predict(multi,X[1:3,,drop=FALSE],type="probabilities")
#' }
NULL

#' @rdname model_predictions
#' @method predict L0glm_fit
#' @export
predict.L0glm_fit <- function(object, newdata=NULL, type="response", offset=NULL,
                             times=NULL, ...) {
  L0glm_predict(object, newdata, type, offset, times)
}

#' @rdname model_predictions
#' @method predict L0glm_extended_fit
#' @export
predict.L0glm_extended_fit <- function(object, newdata=NULL, type="response", offset=NULL,
                                      times=NULL, ...) {
  L0glm_predict(object, newdata, type, offset, times)
}

#' @rdname model_predictions
#' @method predict sparseglm
#' @export
predict.sparseglm <- function(object, newdata=NULL, type="response", offset=NULL,
                             times=NULL, ...) {
  lapply(object$fits, function(fit) if (is.null(fit)) NULL else
    stats::predict(fit, newdata=newdata, type=type, offset=offset, times=times, ...))
}

#' @rdname model_predictions
#' @method fitted L0glm_fit
#' @export
fitted.L0glm_fit <- function(object, ...) {
  if ("newdata" %in% names(list(...))) stop("Use predict() for new data.", call.=FALSE)
  stats::predict(object, newdata=NULL, ...)
}

#' @rdname model_predictions
#' @method fitted sparseglm
#' @export
fitted.sparseglm <- function(object, ...) {
  if ("newdata" %in% names(list(...))) stop("Use predict() for new data.", call.=FALSE)
  stats::predict(object, newdata=NULL, ...)
}

L0glm_predict <- function(object, newdata, type, offset, times) {
  kind <- if (is.null(object$kind)) "glm" else object$kind
  choices <- switch(kind, multinomial=c("response","probabilities","link","class"),
    cox=c("response","risk","link","survival"), c("response","link"))
  if (object$family$family %in% c("binomial","quasibinomial")) choices <- c(choices,"probabilities")
  if (!is.character(type) || length(type) != 1L || is.na(type) || !type %in% choices)
    stop("Prediction type must be one of: ", paste(choices,collapse=", "), ".", call.=FALSE)
  eta <- L0glm_prediction_eta(object, newdata, offset)
  if (kind == "multinomial") {
    eta <- eta - rowMeans(eta)
    colnames(eta) <- object$response$levels
  }
  if (type == "link") return(eta)
  if (kind == "multinomial") {
    probability <- if (nrow(eta)) L0glm_ext_softmax(eta)$p else eta
    if (type == "class") return(factor(object$response$levels[
      if (nrow(eta)) max.col(probability,ties.method="first") else integer()],levels=object$response$levels))
    return(probability)
  }
  if (kind != "cox") {
    if (!nrow(eta)) return(eta)
    return(matrix(object$family$linkinv(as.numeric(eta)),
      nrow(eta),ncol(eta),dimnames=dimnames(eta)))
  }
  if (type %in% c("risk","response")) return(as.numeric(exp(eta)))
  if (is.null(times)) stop("Supply times for Cox survival prediction.",call.=FALSE)
  if (!is.numeric(times) || any(!is.finite(times) | times < 0))
    stop("Survival times must be finite and nonnegative.",call.=FALSE)
  baseline <- object$baseline
  if (is.null(baseline)) baseline <- L0glm_ext_cox_baseline(object)
  index <- findInterval(times,baseline$time)
  loghazard <- if (is.null(baseline$log_cumhaz)) log(baseline$cumhaz) else baseline$log_cumhaz
  loghazard <- c(-Inf,loghazard)[index+1L]
  exp(-exp(outer(as.numeric(eta),loghazard,"+")))
}

L0glm_prediction_eta <- function(object, newdata, offset) {
  beta <- stats::coef(object)
  if (is.null(newdata)) {
    X <- object$X
    if (is.null(offset)) offset <- object$offset
  } else {
    X <- newdata
    p <- nrow(beta) - isTRUE(object$added_intercept)
    if (is.data.frame(X)) {
      if (!all(vapply(X,is.numeric,logical(1)))) stop("newdata must contain numeric predictors.",call.=FALSE)
      X <- as.matrix(X)
    }
    if (is.null(dim(X)) && is.numeric(X)) {
      X <- if (p == 1L) matrix(X,ncol=1L) else matrix(X,nrow=1L,dimnames=list(NULL,names(X)))
    }
    if (!(is.matrix(X) || inherits(X,"Matrix")) || ncol(X) != p)
      stop("Wrong number of prediction columns; supply all original predictors without an added intercept.",call.=FALSE)
    expected <- colnames(object$X)
    if (isTRUE(object$added_intercept)) expected <- expected[-1L]
    provided <- colnames(X)
    if (!is.null(expected) && !is.null(provided) && any(!is.na(expected) & nzchar(expected))) {
      complete <- !anyNA(expected) && all(nzchar(expected)) && !anyDuplicated(expected)
      if (complete) {
        if (anyNA(provided) || anyDuplicated(provided) || !setequal(expected,provided))
          stop("Prediction column names must match the original predictor names.",call.=FALSE)
        if (!identical(expected,provided)) X <- X[,match(expected,provided),drop=FALSE]
      } else if (!identical(expected,provided))
        stop("Prediction columns with ambiguous names must follow their training order.",call.=FALSE)
    }
    entries <- if (inherits(X,"Matrix")) {
      if ("x" %in% methods::slotNames(X)) X@x else numeric()
    } else X
    if ((!is.numeric(entries) && !is.logical(entries)) || any(!is.finite(entries)))
      stop("newdata must contain finite numeric predictors.",call.=FALSE)
    if (isTRUE(object$added_intercept)) X <- cbind(Intercept=rep(1,nrow(X)),X)
  }
  if (ncol(X) != nrow(beta)) stop("Wrong number of prediction columns.",call.=FALSE)
  active <- which(Matrix::rowSums(beta != 0) > 0)
  eta <- if (length(active) == nrow(beta)) as.matrix(X %*% beta) else
    as.matrix(X[,active,drop=FALSE] %*% beta[active,,drop=FALSE])
  if (!is.null(offset) && !is.numeric(offset) && !is.logical(offset) && !inherits(offset,"Matrix"))
    stop("Prediction offsets must be numeric.",call.=FALSE)
  off <- L0glm_tune_matrix(offset,nrow(X),ncol(beta),0,"offset")
  if (any(!is.finite(off))) stop("Prediction offsets must be finite.",call.=FALSE)
  eta + off
}
