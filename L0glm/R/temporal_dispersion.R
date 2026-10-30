#' Estimate noise dispersion from temporal differences
#'
#' Estimate a pooled marginal noise scale in ordered response channels, allowing
#' for an optional lag-one covariance approximation and a slowly changing mean.
#'
#' @param y Finite ordered observations-by-channels matrix.
#' @param family Scalar GLM family, or one family per channel, supplying the
#'   variance function without its dispersion multiplier.
#' @param order Positive difference orders below `nrow(y)`, default `c(2,4,8)`.
#' @param mu Optional independent or smooth pilot means. With `NULL`, Poisson
#'   uses raw observations as a linear variance proxy; nonlinear variance
#'   functions use a local smooth pilot.
#' @param weights Positive observation precisions, scalar, row vector or matrix.
#'   Zero or nonfinite weights exclude difference stencils.
#' @param rows Optional logical eligibility mask, in the same formats as weights.
#' @param minimum_mean Minimum local mean used for relative variance estimates.
#'   Gaussian channels use every otherwise eligible stencil.
#' @param method `"robust"` uses squared median absolute differences calibrated
#'   to Gaussian noise; `"moment"` uses second moments.
#' @param correlation `"independent"` or `"lag1"`. The latter needs at least two
#'   difference orders.
#' @param ncores OpenMP workers; `NULL` uses physical cores.
#'
#' @details
#' With \eqn{a_k=(-1)^k{r\choose k}}, independent locally stationary noise has
#' \eqn{\operatorname{Var}(\Delta^rY_i)=\phi\sum_k a_k^2V(\mu_{i+k})}.
#' Weighted stencil variances replace \eqn{V(\mu)} by \eqn{V(\mu)/w}.
#' The native kernel collects valid, locally variance-normalized differences.
#' For locally constant variance and lag-one covariance \eqn{c_1},
#' \deqn{\frac{\operatorname{Var}(\Delta^rY)}{\sum_k a_k^2}
#' =c_0-\frac{2r}{r+1}c_1.}
#' A nonnegative two-endpoint spectral parameterization fits this relationship
#' across orders. The returned `dispersion` is marginal \eqn{c_0}, `lag1` is
#' \eqn{c_1}, and `long_run` is \eqn{c_0+2c_1}. Slowly changing signals are
#' suppressed by high-order differences. Quieter rows or a supplied pilot can
#' reduce contamination from peaks. More extensive serial dependence requires
#' a richer covariance model.
#'
#' @return A list with per-order channel and pooled estimates, summaries
#'   `channel` and `pooled`, estimation method and covariance model.
#' @seealso [splogconNMF()]
#' @examples
#' set.seed(4)
#' y <- matrix(rnorm(3000, sd=3), 500, 6)
#' temporal_dispersion(y, gaussian(), method="moment")$pooled
#' @export
temporal_dispersion <- function(y, family = stats::quasipoisson("identity"),
                               order = c(2L, 4L, 8L), mu = NULL, weights = NULL,
                               rows = NULL, minimum_mean = 20,
                               method = c("robust", "moment"),
                               correlation = c("independent", "lag1"), ncores = NULL) {
  ncores <- L0glm_sis_threads(ncores)
  method <- match.arg(method);correlation <- match.arg(correlation)
  y <- as.matrix(y);storage.mode(y) <- "double"
  n <- nrow(y);m <- ncol(y)
  if (!m || !is.numeric(order) || !length(order) || any(!is.finite(order)) ||
      any(order != as.integer(order) | order < 1 | order >= n))
    stop("order must contain positive integers smaller than nrow(y).")
  order <- sort(unique(as.integer(order)))
  if (correlation == "lag1" && length(order) < 2L)
    stop("The lag-one covariance approximation needs at least two orders.")
  if (length(minimum_mean) != 1L || !is.finite(minimum_mean) || minimum_mean < 0 ||
      length(ncores) != 1L || !is.finite(ncores) || ncores < 1)
    stop("minimum_mean must be nonnegative and ncores must be positive.")
  expand <- function(value, default, name) {
    if (is.null(value)) return(matrix(default, n, m))
    if (length(value) == 1L) return(matrix(value, n, m))
    if (is.null(dim(value)) && length(value) == n) return(matrix(value, n, m))
    if (!identical(dim(value), c(n, m))) stop(name, " must be scalar, row vector or n by m matrix.")
    as.matrix(value)
  }
  w <- expand(weights, 1, "weights")
  selected <- expand(rows, TRUE, "rows")
  families <- if (is.list(family) && is.function(family$variance)) rep(list(family), m) else family
  if (length(families) != m || any(!vapply(families, function(f) is.function(f$variance), TRUE)))
    stop("family must be a scalar GLM family or a list of one family per channel.")
  if (any(w[is.finite(w)] < 0)) stop("weights cannot be negative.")
  supplied <- !is.null(mu)
  level <- if (supplied) expand(mu, 0, "mu") else y
  variance <- matrix(NA_real_, n, m)
  if (!supplied) {
    # Poisson variance is linear: raw y provides an unbiased variance proxy.
    # Other variance functions use a local pilot; an independent mean is preferable.
    span <- min(n, 2L * max(order) + 1L)
    for (j in seq_len(m)) if (!families[[j]]$family %in% c("poisson", "quasipoisson", "gaussian")) {
      h <- stats::filter(y[, j], rep(1 / span, span), sides = 2)
      good <- which(is.finite(h))
      if (length(good) > 1L) level[, j] <- stats::approx(good, h[good], seq_len(n), rule = 2)$y
    }
  }
  for (j in seq_len(m)) variance[, j] <- families[[j]]$variance(level[, j]) / w[, j]
  filter_level <- level
  for (j in seq_len(m)) if (families[[j]]$family == "gaussian")
    filter_level[, j] <- minimum_mean
  eligible <- selected & is.finite(y) & is.finite(level) & is.finite(variance) &
    variance >= 0 & is.finite(w) & w > 0
  eligible[is.na(eligible)] <- FALSE
  residual <- if (supplied) y - level else y
  result <- lapply(order, function(r) {
    a <- (-1)^(0:r) * choose(r, 0:r)
    ans <- difference_scale_cpp(residual, variance, filter_level, eligible,
      a, minimum_mean, as.integer(ncores))
    ans$channel$order <- r
    ans$pooled <- data.frame(order = r, as.list(ans$pooled), row.names = NULL)
    ans
  })
  by_order <- do.call(rbind, lapply(result, `[[`, "channel"))
  pooled <- do.call(rbind, lapply(result, `[[`, "pooled"))
  summarize <- function(d) {
    value <- d[[method]];valid <- is.finite(value) & d$cells > 0
    if (!any(valid)) return(c(dispersion = NA_real_, lag1 = NA_real_,
      correlation = NA_real_, long_run = NA_real_, relative_error = NA_real_))
    if (correlation == "independent") return(c(dispersion = median(value[valid]),
      lag1 = 0, correlation = 0, long_run = median(value[valid]),
      relative_error = NA_real_))
    if (sum(valid) < 2L) return(c(dispersion = NA_real_, lag1 = NA_real_,
      correlation = NA_real_, long_run = NA_real_, relative_error = NA_real_))
    # For lag-one covariance, Var(diff^r)/sum(a^2) = c0 - 2*r/(r+1)*c1.
    # Parameterize by nonnegative spectrum endpoints to ensure a valid covariance.
    r <- d$order[valid];v <- value[valid];weight <- sqrt(d$cells[valid] / max(d$cells[valid]))
    t <- 1 / (2 * (r + 1));A <- cbind(t, 1 - t) * weight;z <- v * weight
    unconstrained <- stats::lm.fit(A, z)$coefficients
    candidates <- rbind(c(0, 0),
      c(max(0, sum(A[, 1] * z) / sum(A[, 1]^2)), 0),
      c(0, max(0, sum(A[, 2] * z) / sum(A[, 2]^2))))
    if (all(is.finite(unconstrained)) && all(unconstrained >= 0)) candidates <- rbind(candidates, unconstrained)
    loss <- apply(candidates, 1L, function(b) sum((z - A %*% b)^2))
    b <- unname(candidates[which.min(loss), ]);c0 <- sum(b) / 2;c1 <- diff(rev(b)) / 4
    c(dispersion = c0, lag1 = c1, correlation = if (c0 > 0) c1 / c0 else NA_real_,
      long_run = b[1L], relative_error = sqrt(min(loss) / sum(z^2)))
  }
  summary <- do.call(rbind, lapply(seq_len(m), function(j) summarize(by_order[by_order$channel == j, ])))
  colnames(summary) <- c("dispersion", "lag1", "correlation", "long_run", "relative_error")
  list(by_order = by_order, pooled_by_order = pooled,
    channel = data.frame(channel = seq_len(m), summary), pooled = summarize(pooled),
    method = method, covariance_model = correlation, supplied_mean = supplied)
}
