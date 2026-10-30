L0glm_loglik_mean <- function(y, mu, weights, family, dispersion=NULL) {
  m <- ncol(y)
  free_scale <- if(is.null(dispersion)) rep(TRUE,m) else dispersion<0
  -0.5 * sum(vapply(seq_len(m), function(j) {
    keep_j <- weights[, j] > 0 & is.finite(weights[, j])
    yy <- y[keep_j, j]; mm <- mu[keep_j, j]; ww <- weights[keep_j, j]
    if (!length(yy)) return(0)
    if (is.function(family$log_density)) return(-2 * sum(ww * family$log_density(yy, mm)))
    if (family$family == "quasipoisson") {
      scale <- if (is.null(dispersion)) sum(ww * (yy-mm)^2 / pmax(mm, .Machine$double.eps)) /
        max(1, length(yy)-1L) else abs(dispersion[j])
      kernel <- ifelse(yy == 0, 0, yy * log(mm)) - mm - lgamma(yy + 1)
      return(-2 * sum(ww * kernel) / max(scale, .Machine$double.eps))
    }
    phi <- if (is.null(dispersion) || free_scale[j]) NULL else dispersion[j]
    if (family$family == "gaussian") {
      rss <- sum(ww * (yy - mm)^2)
      if (is.null(phi)) phi <- rss / length(yy)
      if (phi == 0) return(-Inf) # unbounded likelihood at an exact interpolant
      return(length(yy) * log(2 * pi * phi) - sum(log(ww)) + rss / phi)
    }
    if (family$family == "poisson") yy <- round(yy)
    if (family$family == "Gamma" && !is.null(phi))
      return(-2 * sum(ww * stats::dgamma(yy, shape = 1 / phi, scale = mm * phi, log = TRUE)))
    if (family$family == "inverse.gaussian" && !is.null(phi))
      return(sum(ww * (log(2 * pi * phi * yy^3) + (yy - mm)^2 / (phi * yy * mm^2))))
    dev <- sum(family$dev.resids(yy, mm, ww))
    value <- family$aic(y = yy, n = ww, mu = mm, wt = ww, dev = dev)
    if (family$family %in% c("Gamma", "inverse.gaussian")) value <- value - 2
    value
  }, numeric(1)))
}
