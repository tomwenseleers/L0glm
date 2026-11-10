# Automatic ISIS: p_eligible > 10 * n, after starting masks and mandatory terms.
L0glm_screen_method <- function(screen, X, start, intercept, always.include) {
  if (!is.character(screen) || length(screen) != 1L || is.na(screen) ||
      !tolower(screen) %in% c("auto", "none", "sis", "isis"))
    stop('screen must be "auto", "none", "SIS" or "ISIS".', call. = FALSE)
  mode <- tolower(screen)
  if (mode != "auto") return(if (mode == "none") mode else toupper(mode))
  n <- nrow(X); p <- ncol(X)
  if (is.null(n) || is.null(p) || p <= 10 * n) return("none")
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
  if (sum(eligible) > 10 * n) "ISIS" else "none"
}

L0glm_screen_fit <- function(X, y, family, weights, offset, nonnegative,
                            method, automatic, nsis, niter, options, nthreads=NULL,
                            candidates=NULL) {
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
    FALSE, nonnegative, if (method == "SIS") 1L else niter, options, nthreads,
    candidates)
  fit <- result$fit
  fit$screening <- c(list(method = method, automatic = automatic), result[c(
    "screen_size", "iterations", "best_iteration", "converged", "candidates", "indices", "history", "seconds", "nthreads")])
  if(!is.null(result$working.set))fit$screening$working.set <- result$working.set
  if(length(result$nuisance_history))fit$screening$nuisance_history <- result$nuisance_history
  if (!is.null(candidates)) fit$screening$fixed.candidates <- TRUE
  fit
}

# Reuse only the predictor pool. Fit each penalty from its usual initialization.
L0glm_path_screen_reuse <- function(args, fit, screen.path) {
  if (screen.path != "preset" || !identical(fit$screening$method, "ISIS")) return(args)
  if (is.null(args$control.solver))
    args$control.solver <- controlSolverGen(
      solver = if (is.null(args$solver)) "eigen" else args$solver,
      method = if (is.null(args$method)) 1L else args$method)
  args$control.solver$.screen.candidates <- fit$screening$candidates
  args
}
