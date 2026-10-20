# Only sufficiently wide eligible designs change the default fitting path.
L0glm_screen_method <- function(screen, X, start, intercept, always.include) {
  if (!is.character(screen) || length(screen) != 1L || is.na(screen) ||
      !tolower(screen) %in% c("auto", "none", "sis", "isis"))
    stop('screen must be "auto", "none", "SIS" or "ISIS".', call. = FALSE)
  mode <- tolower(screen)
  if (mode != "auto") return(if (mode == "none") mode else toupper(mode))
  n <- nrow(X); p <- ncol(X)
  if (is.null(n) || is.null(p) || p < max(10000L, 10 * n)) return("none")
  eligible <- rep(TRUE, p)
  if (!is.null(start)) {
    if (!(is.numeric(start) || inherits(start, "dMatrix"))) return("none")
    values <- if (inherits(start, "sparseMatrix")) start@x else start
    if (any(!is.finite(values))) return("none")
    if (length(start) == 1L) eligible[] <- start != 0 else {
      rows <- if (NROW(start) == p) seq_len(p) else
        if (isTRUE(intercept) && NROW(start) == p + 1L) seq_len(p) + 1L else return("none")
      eligible <- if (is.null(dim(start))) start[rows] != 0 else
        as.numeric(Matrix::rowSums(start[rows, , drop = FALSE] != 0)) > 0
    }
  }
  eligible[L0glm_mandatory(X, always.include)] <- FALSE
  if (sum(eligible) >= max(10000L, 10 * n)) "ISIS" else "none"
}

L0glm_screen_fit <- function(X, y, family, weights, offset, nonnegative,
                            method, automatic, nsis, niter, options, nthreads=1L) {
  X <- L0glm_sis_design(X)
  if (is.null(nsis)) {
    n <- nrow(X); p <- L0glm_predictor_count(X, options$intercept, options$always.include)
    nsis <- if (method == "SIS") 2 * floor(n / log(max(n, 3L))) else
      round(n / (log(max(p, 2L)) * log(log(max(n, 3L)))))
    nsis <- max(1L, nsis)
  }
  if (!is.numeric(nsis) || length(nsis) != 1L || !is.finite(nsis) ||
      nsis < 1 || nsis != floor(nsis)) stop("nsis must be a positive integer.", call. = FALSE)
  if (!is.numeric(niter) || length(niter) != 1L || !is.finite(niter) ||
      niter < 1 || niter != floor(niter)) stop("screen.maxit must be a positive integer.", call. = FALSE)
  result <- L0glm_isis(X, y, family, weights, offset, as.integer(min(nsis, ncol(X))),
    FALSE, nonnegative, if (method == "SIS") 1L else niter, options, nthreads)
  fit <- result$fit
  fit$screening <- c(list(method = method, automatic = automatic), result[c(
    "screen_size", "iterations", "best_iteration", "converged", "candidates", "indices", "history", "seconds", "nthreads")])
  fit
}
