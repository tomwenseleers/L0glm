# Plot completed fits without changing the lambda grid or refitting models.
L0glm_plot_tuning <- function(object, coefs.true = NULL) {
  metrics <- object$metrics
  cv <- object$cv
  do_cv <- object$tune == "cv"
  has_1se <- do_cv && !is.null(object$lambda_1se)
  has_preset <- is.finite(object$lambda_preset) && object$lambda_preset > 0
  x <- log10(metrics$lambda)
  xlim <- rev(range(x))
  beta <- object$coefficients[[1L]]
  selectable <- setdiff(seq_len(nrow(beta)), object$always.include)
  multichannel <- ncol(beta) > 1L
  coefmat <- t(matrix(vapply(object$coefficients, function(b) Matrix::rowMeans(b[selectable, , drop = FALSE]),
                             numeric(length(selectable))), nrow = length(selectable)))
  colors <- grDevices::hcl(h = (15 + 137.508 * (seq_len(ncol(coefmat)) - 1)) %% 360,
                          c = 70, l = 50)

  # Support is discrete: label the nearest evaluated model, never interpolate
  # counts or infer row support from channel means (which can cancel).
  top_at <- seq(xlim[1L], xlim[2L], length.out = 16L)
  nearest <- vapply(top_at, function(at) which.min(abs(x - at)), integer(1))
  support <- metrics$support_size_estimated[nearest]
  annotate <- function(main) {
    graphics::axis(3, at = top_at, labels = support, las = 2, cex.axis = 0.7)
    graphics::mtext(if (multichannel) "Row support size" else "Support size",
                    side = 3, line = 2.6, cex = 0.8)
    graphics::title(main, line = 4.3)
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  graphics::layout(matrix(1:2, nrow = 1))
  margins <- graphics::par("mar")
  margins[3L] <- max(margins[3L], 6.2)
  graphics::par(mar = margins)

  if (ncol(coefmat)) {
    graphics::matplot(x, coefmat, type = "l", lty = 1, col = colors,
                      xlim = xlim, xlab = expression(log[10](lambda)), ylab = "")
  } else graphics::plot(x, rep(0, length(x)), type = "n", xlim = xlim,
                         xlab = expression(log[10](lambda)), ylab = "")
  graphics::mtext(expression(hat(beta)), side = 2, line = 3, las = 1)
  graphics::abline(v = log10(object$lambda_max), col = "black", lty = 3)
  graphics::abline(v = log10(object$lambda_best), col = "red", lty = 2)
  if (has_1se)
    graphics::abline(v = log10(object$lambda_1se), col = "royalblue3", lty = 4)
  if (has_preset)
    graphics::abline(v = log10(object$lambda_preset), col = "red", lty = 3)
  if (!is.null(coefs.true)) {
    if (multichannel) coefs.true <- Matrix::rowMeans(coefs.true)
    coefs.true <- coefs.true[selectable]
    graphics::abline(h = as.vector(coefs.true[coefs.true != 0]),
                     col = scales::alpha("blue", 0.1))
  }
  annotate("Regularisation path")

  selection <- if (do_cv) metrics$cv_deviance else metrics[[object$tune]]
  finite <- is.finite(selection)
  limits <- if (do_cv) range(cv$mean[finite] - cv$se[finite],
                             cv$mean[finite] + cv$se[finite]) else NULL
  graphics::plot(x, replace(selection, !finite, NA_real_), type = "b", pch = 19,
    xlim = xlim,
    log = if (!do_cv && any(finite) && all(selection[finite] > 0)) "y" else "",
    ylim = if (any(finite)) limits else c(0, 1),
    xlab = expression(log[10](lambda)),
    ylab = if (do_cv) "Mean held-out deviance" else toupper(object$tune))
  if (do_cv) {
    bars <- finite & is.finite(cv$se) & cv$se > 0
    if (any(bars)) {
      lower <- cv$mean[bars] - cv$se[bars]
      upper <- cv$mean[bars] + cv$se[bars]
      width <- abs(graphics::xinch(0.04))
      graphics::segments(x[bars], lower, x[bars], upper, col = "grey40")
      graphics::segments(x[bars] - width, lower, x[bars] + width, lower, col = "grey40")
      graphics::segments(x[bars] - width, upper, x[bars] + width, upper, col = "grey40")
    }
  }
  graphics::abline(v = log10(object$lambda_best), col = "red", lty = 2)
  if (has_1se)
    graphics::abline(v = log10(object$lambda_1se), col = "royalblue3", lty = 4)
  if (has_preset)
    graphics::abline(v = log10(object$lambda_preset), col = "red", lty = 3)
  if (do_cv)
    graphics::legend("top", c("Best CV", if (has_1se) "CV 1-SE", if (has_preset) "GIC preset"),
      col = c("red", if (has_1se) "royalblue3", if (has_preset) "red"),
      lty = c(2, if (has_1se) 4, if (has_preset) 3),
      bty = "o", box.col = NA, bg = "white", cex = 0.8)
  annotate(if (do_cv) paste0(cv$nfolds, "-fold cross-validation") else
             paste("Tuning -", toupper(object$tune)))
  if (any(!finite)) graphics::mtext("Non-finite scores omitted", side = 3,
                                    line = 5.5, cex = 0.7)
  invisible(NULL)
}
