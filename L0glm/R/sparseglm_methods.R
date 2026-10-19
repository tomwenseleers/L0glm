# S3 extraction and plotting for sparse-method comparisons.

#' @rdname L0glm.metrics
#' @export
L0glm.metrics.sparseglm <- function(object, newdata=NULL, y=NULL, weights=NULL,
  offset=NULL, mu_floor=object$mu_floor, beta_true=object$beta_true, r2_mu_floor=0.1, ...) {
  if (!is.null(newdata) && is.null(y)) stop("Supply y when scoring newdata.", call.=FALSE)
  rows <- lapply(seq_len(nrow(object$timings)), function(i) {
    fit <- object$fits[[i]]
    meta <- object$timings[i, , drop=FALSE]
    if (is.null(fit)) {
      recovery <- L0glm_recovery_metrics(Matrix::Matrix(0, ncol(object$data$X) +
        object$data$intercept, object$data$m, sparse=TRUE))
      empty <- c(list(p=ncol(object$data$X), n=nrow(object$data$X), m=object$data$m,
        df=NA_real_, loglik=NA_real_, deviance=NA_real_, deviance_true=NA_real_,
        R2_weighted=NA_real_, adj_R2_weighted=NA_real_,
        gic=NA_real_, ebic=NA_real_, bic=NA_real_, aic=NA_real_, mbic=NA_real_), recovery)
      empty$support_size_estimated <- NA_real_
      return(cbind(meta, as.data.frame(empty)))
    }
    scores <- L0glm.metrics(fit, newdata=newdata, y=y, weights=weights,
      offset=offset, mu_floor=mu_floor, beta_true=beta_true, r2_mu_floor=r2_mu_floor)
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
#' @param ... Further arguments passed to [graphics::plot()] for coefficient
#'   plots or [plotmat()] for heatmaps; unused by `coef()`.
#' @export
coef.sparseglm <- function(object, ...) {
  lapply(object$fits, function(fit) if (is.null(fit)) NULL else stats::coef(fit))
}

#' @rdname model_summary
#' @method print sparseglm
#' @export
print.sparseglm <- function(x, digits=max(3L, getOption("digits") - 3L), ...) {
  L0glm_print_comparison_info(x$data$family, nrow(x$data$X), ncol(x$data$X),
    x$data$m, x$data$tune, if (is.null(x$foldid)) NULL else max(x$foldid))
  print(x$timings[c("method", "status", "scope", "parameter", "selected_parameter", "seconds")],
    digits=digits, row.names=FALSE)
  invisible(x)
}

#' @rdname model_summary
#' @method summary sparseglm
#' @export
summary.sparseglm <- function(object, ..., mu_floor=object$mu_floor, beta_true=object$beta_true,
                               r2_mu_floor=0.1) {
  structure(list(call=object$call, family=object$data$family,
    n=nrow(object$data$X), p=ncol(object$data$X), m=object$data$m,
    tune=object$data$tune, nfolds=if (is.null(object$foldid)) NULL else max(object$foldid),
    mu_floor=mu_floor, r2_mu_floor=r2_mu_floor,
    metrics=L0glm.metrics(object, mu_floor=mu_floor, beta_true=beta_true, r2_mu_floor=r2_mu_floor)),
    class="summary.sparseglm")
}

#' @rdname model_summary
#' @method print summary.sparseglm
#' @export
print.summary.sparseglm <- function(x, digits=max(3L, getOption("digits") - 3L), ...) {
  L0glm_print_comparison_info(x$family, x$n, x$p, x$m, x$tune, x$nfolds)
  if (x$mu_floor > 0) cat("Scoring mean floor:", format(x$mu_floor, digits=digits), "\n")
  fields <- unique(c("method", "status", "scope", "parameter", "selected_parameter", "seconds",
    if (x$tune == "cv") c("cv_mean", "cv_se") else x$tune,
    "support_size_estimated", "loglik", "deviance", "adj_R2_weighted", "gic", L0glm_summary_recovery(x$metrics)))
  print(x$metrics[intersect(fields, names(x$metrics))], digits=digits, row.names=FALSE)
  notes <- which(!is.na(x$metrics$note) & nzchar(x$metrics$note))
  if (length(notes)) {
    cat("\nMethod notes:\n")
    for (i in notes) cat(x$metrics$method[i], ": ", x$metrics$note[i], "\n", sep="")
  }
  invisible(x)
}

L0glm_print_comparison_info <- function(family, n, p, m, tune, nfolds) {
  cat("Sparse GLM comparison\nFamily:", family$family, "  Link:", family$link, "\n")
  cat(sprintf("Observations: %d   Predictors: %d   Channels: %d\n", n, p, m))
  cat("Tuning:", toupper(tune))
  if (!is.null(nfolds)) cat(sprintf(" (%d folds)", nfolds))
  cat("\n")
}

#' @rdname sparseglm
#' @param x A \code{sparseglm} result.
#' @param which Numeric method indices, in the order stored in `x$fits`.
#'   Defaults to all available methods for `type="coefs"` and method 1 for
#'   `type="heatmap"`. A vector preserves the requested order; unavailable
#'   selected methods are skipped with a warning.
#' @param type `"coefs"` for coefficient comparisons, or `"heatmap"` for
#'   multichannel ground-truth and estimated coefficient images. Each selected
#'   heatmap method occupies a separate page with two rows. Heatmaps require
#'   `beta_true` and omit an automatically added intercept row. Cyan lines show
#'   true support above, and TP below; FP are red and FN magenta. Lines are
#'   semitransparent. Both images use the same colour scale; signed coefficients
#'   use a diverging palette. Support is determined before image downscaling.
#' @param channels `"all"` for individual coefficient entries or `"mean"` for
#'   row means in coefficient plots. Row-mean colours use support across the
#'   original channels. A coincident `(0,0)` point is drawn only once.
#' @export
plot.sparseglm <- function(x, which=NULL, type=c("coefs", "heatmap"),
                          channels=c("all", "mean"), ...) {
  type <- match.arg(type); channels <- match.arg(channels)
  available <- base::which(!vapply(x$fits, is.null, logical(1)))
  if (is.null(which)) which <- if (type == "heatmap") 1L else available
  if (!length(which)) return(invisible(x))
  if (!is.numeric(which) || anyNA(which) || any(!is.finite(which)) ||
      any(which != floor(which) | which < 1 | which > length(x$fits)))
    stop("which must contain positive whole-number method indices.", call.=FALSE)
  unavailable <- setdiff(which, available)
  if (length(unavailable)) warning("Unavailable methods omitted: ",
    paste(names(x$fits)[unavailable], collapse=", "), call.=FALSE)
  selected <- which[which %in% available]
  if (!length(selected)) return(invisible(x))
  if (type == "heatmap" && (is.null(x$beta_true) || ncol(coef(x$fits[[selected[1L]]])) < 2L))
    stop("Heatmaps require multichannel coefficients and beta_true.", call.=FALSE)
  if (type == "heatmap" && channels != "all")
    stop("channels='mean' is only available for coefficient plots.", call.=FALSE)
  old <- graphics::par(no.readonly=TRUE); on.exit(graphics::par(old))
  columns <- min(2L, length(selected))
  graphics::par(mfrow=if (type == "heatmap") c(2L, 1L) else
    c(ceiling(length(selected)/columns), columns))
  dots <- list(...)
  for (i in selected) {
    beta <- stats::coef(x$fits[[i]])
    title <- sprintf("%s (%.3g s)", names(x$fits)[i], x$timings$seconds[i])
    if (type == "heatmap") {
      L0glm_plot_heatmap(beta, x$beta_true, title, isTRUE(x$data$intercept), dots)
    } else {
      points <- L0glm_plot_coefficients(beta, x$beta_true, channels)
      if (is.null(points$truth)) {
        args <- list(x=points$index, y=points$estimated, type="h",
          xlim=c(1, points$size), ylim=range(c(0, points$estimated)),
          xlab=if (channels == "mean") "Coefficient row" else "Coefficient entry",
          ylab=expression(hat(beta)), main=title)
      } else {
        limits <- range(c(0, points$estimated, points$truth))
        args <- list(x=points$estimated, y=points$truth, pch=16,
          col=points$color, cex=.6, xlim=limits, ylim=limits,
          xlab=if (channels == "mean") "Estimated row means" else "Estimated coefficients",
          ylab=if (channels == "mean") "True row means" else "True coefficients", main=title)
      }
      do.call(graphics::plot.default, utils::modifyList(args, dots, keep.null=TRUE))
      if (!is.null(points$truth)) graphics::abline(0, 1, lty=2, col="grey60")
    }
  }
  invisible(x)
}

# Extract sparse nonzero entries before converting plotting vectors to dense.
L0glm_plot_nonzero <- function(B) {
  if (!inherits(B, "sparseMatrix")) return(base::which(as.vector(B) != 0))
  entries <- summary(B)
  present <- if (is.null(entries$x)) rep(TRUE, nrow(entries)) else entries$x != 0
  entries$i[present] + (entries$j[present] - 1L) * nrow(B)
}

L0glm_plot_coefficients <- function(beta, truth, channels) {
  if (!is.null(truth) && is.null(dim(truth))) truth <- matrix(truth, nrow(beta), ncol(beta))
  if (channels == "mean") {
    present <- Matrix::rowSums(abs(beta)) != 0
    true_support <- if (!is.null(truth)) Matrix::rowSums(abs(truth)) != 0
    beta <- Matrix::rowMeans(beta)
    if (!is.null(truth)) truth <- Matrix::rowMeans(truth)
  }
  size <- length(beta)
  index <- sort(unique(c(L0glm_plot_nonzero(beta),
    if (!is.null(truth)) L0glm_plot_nonzero(truth))))
  if (!is.null(truth) && length(index) < size) {
    gaps <- base::which(index != seq_along(index))
    # Retain exactly one common zero without allocating seq_len(size).
    index <- c(if (length(gaps)) gaps[1L] else length(index) + 1L, index)
  }
  estimated <- as.numeric(beta[index])
  actual <- if (!is.null(truth)) as.numeric(truth[index])
  color <- if (!is.null(truth)) {
    selected <- if (channels == "mean") present[index] else estimated != 0
    positive <- if (channels == "mean") true_support[index] else actual != 0
    ifelse(!positive & selected, "red", ifelse(positive & !selected, "magenta", "navy"))
  }
  list(index=index, estimated=estimated, truth=actual, color=color, size=size)
}

L0glm_plot_heatmap <- function(beta, truth, title, intercept, dots) {
  if (intercept) {
    beta <- beta[-1L, , drop=FALSE]; truth <- truth[-1L, , drop=FALSE]
  }
  present <- base::which(Matrix::rowSums(abs(beta)) != 0)
  positive <- base::which(Matrix::rowSums(abs(truth)) != 0)
  signed <- min(min(beta), min(truth)) < 0
  maximum <- max(max(abs(beta)), max(abs(truth)))
  args <- utils::modifyList(list(maxresolution=2000, gamma=.13, signed=signed,
    xlab="Time index", ylab="Channel"), dots, keep.null=TRUE)
  if (signed && is.null(args$cols))
    args$cols <- grDevices::colorRampPalette(c("navy", "white", "firebrick"))(1000)
  if (is.null(args$zlim)) args$zlim <- (if (signed) c(-1, 1) else c(0, 1)) *
    if (maximum > 0) maximum^args$gamma else 1
  p <- nrow(beta)
  stride <- if (is.null(args$maxresolution)) 1 else ceiling(p/min(p, args$maxresolution))
  image_rows <- ceiling(p/stride)
  position <- function(rows) if (p == 1L) rep(1, length(rows)) else
    1 + (rows - 1) * (image_rows - 1)/(p - 1)
  draw <- function(B, main) {
    rownames(B) <- as.character(seq_len(p))
    do.call(plotmat, c(list(M=B, main=main), args[setdiff(names(args), "main")]))
  }
  cyan <- grDevices::adjustcolor("cyan", alpha.f=.4)
  red <- grDevices::adjustcolor("red", alpha.f=.4)
  magenta <- grDevices::adjustcolor("magenta", alpha.f=.4)
  draw(truth, "Ground truth")
  graphics::abline(v=position(positive), col=cyan)
  graphics::legend("topright", "True support", col=cyan, lty=1,
    bg=grDevices::adjustcolor("white", alpha.f=.85), cex=.75)
  draw(beta, title)
  graphics::abline(v=position(intersect(present, positive)), col=cyan)
  graphics::abline(v=position(setdiff(present, positive)), col=red)
  graphics::abline(v=position(setdiff(positive, present)), col=magenta)
  graphics::legend("topright", c("TP", "FP", "FN"), col=c(cyan, red, magenta),
    lty=1, bg=grDevices::adjustcolor("white", alpha.f=.85), cex=.75)
}
