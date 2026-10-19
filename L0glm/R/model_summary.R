#' Print and summarize sparse regression fits
#'
#' @name model_summary
#' @md
#' @description
#' Concise S3 displays for [L0glm.fit()], [L0nnpois()] and [sparseglm()]
#' results, including negative-binomial, multinomial and Cox models.
#'
#' @param object An `L0glm_fit` model or `sparseglm` comparison.
#' @param x A fitted model, comparison or corresponding summary object.
#' @param ... Reserved for additional method arguments; currently unused.
#' @param mu_floor Nonnegative scoring-only mean floor. Model summaries default
#'   to zero, matching [stats::logLik()] and the criterion methods. Comparison
#'   summaries use the floor stored by [sparseglm()] to match tuning scores.
#' @param beta_true Optional true coefficients for recovery metrics, with the
#'   same dimensions as the returned coefficients, including any intercept.
#'   Comparison summaries default to their saved truth. Use `NULL` to omit it.
#' @param digits Significant digits for printed numeric values.
#' @param max.rows,max.cols Maximum coefficient rows and response columns shown
#'   by `print.summary.L0glm_fit()`, default 20 and 6. Positive integers or `Inf`.
#'   Only nonzero or mandatory coefficient rows are displayed. These limits do
#'   not subset the coefficients retained in the summary.
#'
#' @details
#' Printing a fitted object uses its stored metadata and coefficients; it does
#' not evaluate likelihoods or refit a model. `summary()` calls [L0glm.metrics()]
#' to obtain likelihood, deviance and information criteria, optionally including
#' recovery metrics. The full original-dimension sparse coefficient matrix is
#' retained in the model summary. Coefficients are point estimates; standard
#' errors and Wald tests are not calculated for the selected penalized model.
#'
#' Convergence counts use the last recorded pass for each stage and channel,
#' rather than counting intermediate grouped-selection passes as final failures.
#' Missing backend diagnostics remain unavailable. The printed dispersion is
#' the stored fitting scale; the profiled Gaussian likelihood scale can differ,
#' as described in [model_criteria]. NB summaries also retain theta and its
#' estimation flag. Cox and multinomial do not display an ordinary dispersion.
#'
#' Comparison summaries retain the complete metrics table, including timings,
#' tuning metadata and explanatory notes for unavailable or failed methods.
#' Their printed tables show essential scores, plus recovery measures when
#' truth is available. No fitted models or adaptive pilots are recomputed.
#'
#' @return `summary()` returns a list of class `summary.L0glm_fit` or
#'   `summary.sparseglm`. Both contain `metrics` and model dimensions. Model
#'   summaries also contain the full sparse `coefficients`, selected row
#'   indices, dispersion/theta, and a compact convergence table. `print()`
#'   methods return their input invisibly.
#' @seealso [L0glm.metrics()], [model_criteria], [sparseglm()]
#' @examples
#' set.seed(4)
#' X <- matrix(rnorm(240), 60, 4)
#' y <- 2 * X[, 1] + rnorm(60)
#' fit <- L0glm.fit(X, y, lambda="gic", parallelize=FALSE)
#' print(fit)
#' details <- summary(fit)
#' print(details, max.rows=10)
#' details$metrics
#' comparison <- sparseglm(X, y, method=c("L0glm","lasso"), nlambda=4,
#'                         lambda=list(L0glm=c(10,4,2)))
#' summary(comparison)
NULL

#' @rdname model_summary
#' @method print L0glm_fit
#' @export
print.L0glm_fit <- function(x, digits=max(3L, getOption("digits") - 3L), ...) {
  L0glm_print_model_info(L0glm_model_info(x), digits)
  invisible(x)
}

#' @rdname model_summary
#' @method summary L0glm_fit
#' @export
summary.L0glm_fit <- function(object, ..., mu_floor=0, beta_true=NULL) {
  structure(c(L0glm_model_info(object), list(coefficients=stats::coef(object),
    metrics=L0glm.metrics(object, mu_floor=mu_floor, beta_true=beta_true),
    mu_floor=mu_floor, dispersion=object$dispersion,
    dispersion.estimated=object$dispersion.estimated,
    theta=object$theta, theta.estimated=object$theta.estimated)),
    class="summary.L0glm_fit")
}

#' @rdname model_summary
#' @method print summary.L0glm_fit
#' @export
print.summary.L0glm_fit <- function(x, digits=max(3L, getOption("digits") - 3L),
                                  max.rows=20L, max.cols=6L, ...) {
  limit <- function(value) length(value) == 1L && is.numeric(value) &&
    !is.na(value) && value >= 1 && (is.infinite(value) || value == floor(value))
  if (!limit(max.rows) || !limit(max.cols))
    stop("max.rows and max.cols must be positive integers or Inf.", call.=FALSE)
  L0glm_print_model_info(x, digits)
  if (x$mu_floor > 0) cat("Scoring mean floor:", format(x$mu_floor, digits=digits), "\n")
  cat("\nFit metrics:\n")
  fields <- intersect(c("loglik", "deviance", "adj_R2_weighted", "df", "aic", "bic", "gic", "ebic", "mbic"), names(x$metrics))
  print(x$metrics[fields], digits=digits, row.names=FALSE)
  recovery <- L0glm_summary_recovery(x$metrics)
  if (length(recovery)) {
    cat("\nRecovery metrics:\n")
    print(x$metrics[recovery], digits=digits, row.names=FALSE)
  }
  scale <- if (identical(x$kind, "nb")) x$theta else
    if (!x$kind %in% c("cox", "multinomial")) x$dispersion else NULL
  if (length(scale)) {
    estimated <- if (identical(x$kind, "nb")) x$theta.estimated else x$dispersion.estimated
    cat(if (identical(x$kind, "nb")) "\nNB theta:\n" else "\nFitting dispersion:\n")
    index <- seq_len(min(length(scale), max.cols))
    values <- data.frame(channel=index, value=as.numeric(scale[index]),
      estimated=if (length(estimated)) rep_len(estimated, length(scale))[index] else NA)
    print(values, digits=digits, row.names=FALSE)
    if (length(index) < length(scale)) cat("Remaining channels omitted.\n")
  }
  rows <- sort(unique(c(x$selected_rows, x$mandatory_rows)))
  cat("\nCoefficient estimates:\n")
  if (!length(rows)) cat("No selected coefficients.\n") else {
    shown <- rows[seq_len(min(length(rows), max.rows))]
    columns <- seq_len(min(ncol(x$coefficients), max.cols))
    values <- as.matrix(x$coefficients[shown, columns, drop=FALSE])
    if (is.null(rownames(values))) rownames(values) <- as.character(shown)
    if (is.null(colnames(values))) colnames(values) <- paste0("Response", columns)
    print(values, digits=digits)
    if (length(shown) < length(rows) || length(columns) < ncol(x$coefficients))
      cat(sprintf("Showing %d of %d nonzero or mandatory rows and %d of %d channels.\n",
        length(shown), length(rows), length(columns), ncol(x$coefficients)))
  }
  invisible(x)
}

# Lightweight metadata shared by model and comparison displays.
L0glm_model_info <- function(object) {
  beta <- stats::coef(object)
  mandatory <- object$always.include
  selected <- setdiff(which(Matrix::rowSums(beta != 0) > 0), mandatory)
  list(family=object$family$family, link=object$family$link,
    kind=if (is.null(object$kind)) "glm" else object$kind,
    n=nrow(object$X), p=object$p, m=ncol(beta), nobs=stats::nobs(object),
    selected_rows=selected, mandatory_rows=mandatory,
    nonzero=sum(beta != 0) - sum(beta[mandatory, , drop=FALSE] != 0),
    lambda=object$init.lambda, grouped=object$grouped,
    convergence=L0glm_summary_convergence(object$convergence))
}

L0glm_summary_convergence <- function(convergence) {
  if (is.null(convergence) || !nrow(convergence)) return(NULL)
  if (is.null(convergence$stage)) convergence$stage <- "fit"
  if (is.null(convergence$channel)) convergence$channel <- 1L
  last <- !duplicated(convergence[c("stage", "channel")], fromLast=TRUE)
  stages <- factor(convergence$stage[last], levels=unique(convergence$stage))
  rows <- split(convergence[last, , drop=FALSE], stages, drop=TRUE)
  do.call(rbind, lapply(rows, function(data) data.frame(stage=data$stage[1L],
    channels=nrow(data), converged=sum(data$converged %in% TRUE),
    not_converged=sum(data$converged %in% FALSE), unavailable=sum(is.na(data$converged)),
    max_iterations=if (any(is.finite(data$iterations))) max(data$iterations, na.rm=TRUE) else NA_real_)))
}

L0glm_print_model_info <- function(x, digits) {
  cat("L0glm fit\nFamily:", x$family, "  Link:", x$link, "\n")
  cat(sprintf("Observations: %d   Selectable predictors: %d   Channels: %d\n", x$n, x$p, x$m))
  cat(sprintf("Selected rows: %d   Nonzero selectable coefficients: %d   Mandatory rows: %d\n",
    length(x$selected_rows), x$nonzero, length(x$mandatory_rows)))
  if (length(x$lambda) && !all(is.na(x$lambda)))
    cat("Lambda:", paste(format(x$lambda, digits=digits), collapse=", "), "\n")
  if (x$m > 1L && !is.null(x$grouped)) cat("Grouped selection:", x$grouped, "\n")
  if (x$kind == "cox") cat("Criterion sample size (Cox):", x$nobs, "\n")
  if (is.null(x$convergence)) cat("Convergence diagnostics: unavailable\n") else {
    cat("Convergence (last pass per stage/channel):\n")
    print(x$convergence, row.names=FALSE)
  }
}

L0glm_summary_recovery <- function(metrics) {
  fields <- intersect(c("F1_score", "MCC", "within_row_F1_score", "within_row_MCC",
    "MAD_coefs", "RMSE_coefs", "R2_coefs"), names(metrics))
  fields[vapply(metrics[fields], function(x) any(!is.na(x)), logical(1))]
}
