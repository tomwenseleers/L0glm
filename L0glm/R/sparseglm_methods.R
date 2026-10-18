# S3 extraction and plotting for sparse-method comparisons.

#' @rdname L0glm.metrics
#' @export
L0glm.metrics.sparseglm <- function(object, newdata=NULL, y=NULL, weights=NULL,
  offset=NULL, mu_floor=object$mu_floor, beta_true=object$beta_true, ...) {
  if (!is.null(newdata) && is.null(y)) stop("Supply y when scoring newdata.", call.=FALSE)
  rows <- lapply(seq_len(nrow(object$timings)), function(i) {
    fit <- object$fits[[i]]
    meta <- object$timings[i, , drop=FALSE]
    if (is.null(fit)) {
      recovery <- L0glm_recovery_metrics(Matrix::Matrix(0, ncol(object$data$X) +
        object$data$intercept, object$data$m, sparse=TRUE))
      empty <- c(list(p=ncol(object$data$X), n=nrow(object$data$X), m=object$data$m,
        df=NA_real_, loglik=NA_real_, deviance=NA_real_, deviance_true=NA_real_,
        gic=NA_real_, ebic=NA_real_, bic=NA_real_, aic=NA_real_, mbic=NA_real_), recovery)
      empty$support_size_estimated <- NA_real_
      return(cbind(meta, as.data.frame(empty)))
    }
    scores <- L0glm.metrics(fit, newdata=newdata, y=y, weights=weights,
      offset=offset, mu_floor=mu_floor, beta_true=beta_true)
    if (!"p" %in% names(scores)) scores$p <- fit$p
    if (!"n" %in% names(scores)) scores$n <- if (is.null(newdata)) nrow(fit$X) else nrow(newdata)
    if (!"m" %in% names(scores)) scores$m <- ncol(coef(fit))
    if (!"support_size_estimated" %in% names(scores)) scores <- cbind(scores,
      as.data.frame(L0glm_recovery_metrics(coef(fit), beta_true, fit$always.include, fit$p)))
    cbind(meta, scores)
  })
  columns <- unique(unlist(lapply(rows, names)))
  rows <- lapply(rows, function(row) {
    for (name in setdiff(columns, names(row))) row[[name]] <- NA
    row[columns]
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

#' @rdname sparseglm
#' @param object A \code{sparseglm} result.
#' @param ... Further method arguments; currently unused.
#' @export
coef.sparseglm <- function(object, ...) {
  lapply(object$fits, function(fit) if (is.null(fit)) NULL else stats::coef(fit))
}

#' @rdname sparseglm
#' @method print sparseglm
#' @export
print.sparseglm <- function(x, ...) {
  print(x$timings[c("method", "status", "scope", "tune", "selected_parameter", "seconds")], row.names=FALSE)
  invisible(x)
}

#' @rdname sparseglm
#' @param x A \code{sparseglm} result.
#' @export
plot.sparseglm <- function(x, ...) {
  available <- which(!vapply(x$fits, is.null, logical(1)))
  if (!length(available)) return(invisible(x))
  old <- graphics::par(no.readonly=TRUE); on.exit(graphics::par(old))
  columns <- min(2L, length(available))
  graphics::par(mfrow=c(ceiling(length(available)/columns), columns))
  for (i in available) {
    beta <- as.numeric(stats::coef(x$fits[[i]]))
    title <- sprintf("%s (%.3g s)", names(x$fits)[i], x$timings$seconds[i])
    if (is.null(x$beta_true)) {
      graphics::plot(beta, type="h", xlab="Coefficient entry", ylab=expression(hat(beta)), main=title)
    } else {
      truth <- as.numeric(x$beta_true)
      color <- ifelse(truth == 0 & beta != 0, "red", ifelse(truth != 0 & beta == 0, "magenta", "navy"))
      limits <- range(c(beta, truth))
      graphics::plot(beta, truth, pch=16, col=color, cex=.6, xlim=limits, ylim=limits,
        xlab="Estimated coefficients", ylab="True coefficients", main=title)
      graphics::abline(0, 1, lty=2, col="grey60")
    }
  }
  invisible(x)
}
