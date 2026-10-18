L0glm_ext_cox_response <- function(y, n) {
  if (inherits(y, "stratifySurv") || !is.null(attr(y, "strata")))
    stop("Cox strata are not supported.", call. = FALSE)
  if (inherits(y, "Surv") && !identical(attr(y, "type"), "right"))
    stop("Cox fitting supports right censoring only.", call. = FALSE)
  y <- as.matrix(y)
  if (nrow(y) != n || ncol(y) != 2L || any(!is.finite(y)) ||
      any(y[, 1L] < 0) || any(!y[, 2L] %in% c(0, 1)))
    stop("Use Surv(time,status) or a two-column time/status matrix for Cox.", call. = FALSE)
  list(y = unclass(y), time = as.numeric(y[, 1L]), status = as.numeric(y[, 2L]), m = 1L)
}

# Fixed event/risk bookkeeping, shared by every likelihood and working step.
L0glm_ext_cox_cache <- function(response, weights) {
  order <- order(response$time, decreasing = TRUE)
  w <- weights[order]; status <- response$status[order]
  runs <- rle(response$time[order])
  ends <- cumsum(runs$lengths)
  group <- rep(seq_along(ends), runs$lengths)
  list(order = order, weights = w, event = w * status,
    d = as.numeric(rowsum(w * status, group, reorder = FALSE)),
    ends = ends, group = group, time = runs$values)
}

# Weighted Breslow risk sets; censored subjects remain in the denominators.
L0glm_ext_cox_eta <- function(eta, response, weights, cache = NULL, sorted = FALSE) {
  if (is.null(cache)) cache <- L0glm_ext_cox_cache(response, weights)
  if (!sorted) eta <- eta[cache$order]
  if (any(!is.finite(eta))) return(list(loglik = -Inf))
  w <- cache$weights; group <- cache$group; ends <- cache$ends; d <- cache$d
  shift <- max(eta[w > 0])
  e <- w * exp(eta - shift)
  e[w == 0] <- 0
  S0 <- cumsum(e)[ends]
  if (any(S0[d > 0] <= 0))
    stop("Cox risk-set sums underflowed; use less extreme starts.", call. = FALSE)
  increments <- ifelse(S0 > 0, d/S0, 0)
  q <- e * rev(cumsum(rev(increments)))[group]
  event <- cache$event > 0
  list(loglik = sum(cache$event[event] * (eta[event] - shift - log(S0[group[event]]))),
    order = cache$order, eta = eta, e = e, S0 = S0, d = d, group = group,
    ends = ends, time = cache$time, shift = shift, increments = increments,
    q = q, score = cache$event - q)
}

L0glm_ext_cox_stats <- function(data, B) {
  cache <- data$cox.cache
  eta <- as.numeric(data$Z %*% B) + data$offset[cache$order, 1L]
  L0glm_ext_cox_eta(eta, data$response, data$weights[, 1L], cache, sorted = TRUE)
}

L0glm_ext_cox_initial <- function(data) {
  B <- matrix(0, data$p, 1L)
  D <- matrix(1e-3, data$p, 1L)
  proposal <- L0glm_ext_cox_ridge_pilot(data)
  diagnostic <- attr(proposal, "cox.solve"); attr(proposal, "cox.solve") <- NULL
  if (!diagnostic$converged)
    stop(sprintf("The matrix-free Cox ridge initializer did not converge in %d iterations (scaled residual %.3g); increase cox.cg.maxit.",
      diagnostic$iterations, diagnostic$residual), call. = FALSE)
  proposal <- pmax(proposal, data$lower)
  objective <- function(b) -2 * L0glm_ext_cox_stats(data, b)$loglik + sum(D * b^2)
  baseline <- objective(B)
  for (halving in 0:30) {
    candidate <- B + (proposal - B)/2^halving
    value <- objective(candidate)
    if (is.finite(value) && value <= baseline) break
  }
  if (!is.finite(value) || value > baseline)
    stop("No descent step for the Cox ridge initializer.", call. = FALSE)
  diagnostic$step_halvings <- halving
  attr(candidate, "cox.solve") <- diagnostic
  candidate
}

L0glm_ext_whiten <- function(W, eta, score) {
  eig <- eigen(as.matrix(W), symmetric = TRUE)
  keep <- eig$values > max(eig$values) * 1e-13
  if (!any(keep)) return(list(root = matrix(0, 0L, nrow(W)), response = numeric()))
  Q <- t(eig$vectors[, keep, drop = FALSE])
  root <- sqrt(eig$values[keep]) * Q
  list(root = root, response = as.numeric(root %*% eta) +
    as.numeric(Q %*% score)/sqrt(eig$values[keep]))
}

# Exact small-system information: no per-observation p-by-p R loop.
L0glm_ext_cox_information <- function(X, stats) {
  event <- which(stats$d > 0)
  means <- vapply(seq_len(ncol(X)), function(j)
    cumsum(stats$e * as.numeric(X[, j]))[stats$ends[event]]/stats$S0[event],
    numeric(length(event)))
  means <- matrix(means, length(event), ncol(X)) * sqrt(stats$d[event])
  as.matrix(Matrix::crossprod(X, Matrix::Diagonal(x = stats$q) %*% X)) - crossprod(means)
}

L0glm_ext_cox_direct <- function(data, B, D, stats, args) {
  active <- is.finite(D[, 1L])
  answer <- matrix(0, data$p, 1L)
  k <- sum(active)
  if (!k) return(answer)
  X <- data$Z[, active, drop = FALSE]
  if (min(k, nrow(X)) > data$control$max.working)
    stop("Cox working system exceeds control.family$max.working; screen predictors or raise the limit.", call. = FALSE)
  beta <- B[active, 1L]
  d <- D[active, 1L]
  if (k <= nrow(X)) {
    H <- L0glm_ext_cox_information(X, stats)
    score <- as.numeric(Matrix::crossprod(X, stats$score))
    if (args$solver == "eigen") {
      proposal <- eigen_ridge_solve(H + diag(d, k), H %*% beta + score + d * data$prior[active, 1L],
        beta, args$control.solver, 1L)
      if (is.null(proposal)) stop("Cox quadratic solver failed.", call. = FALSE)
      answer[active, 1L] <- pmax(data$lower[active], proposal)
      attr(answer, "cox.solve") <- list(backend = "llt", iterations = 0L, residual = NA_real_, converged = TRUE)
      return(answer)
    }
    working <- L0glm_ext_whiten(H, beta, score)
    A <- working$root
  } else {
    C <- rev(cumsum(rev(ifelse(stats$S0 > 0, stats$d/stats$S0^2, 0))))
    index <- outer(stats$group, stats$group, pmax)
    W <- diag(stats$q) - tcrossprod(stats$e) * matrix(C[index], nrow(X), nrow(X))
    working <- L0glm_ext_whiten(W, stats$eta - data$offset[stats$order, 1L], stats$score)
    A <- working$root %*% X
  }
  answer[active, 1L] <- L0glm_ext_ridge(A, working$response, NULL, d, beta,
    data$lower[active], data$prior[active, 1L], args)
  attr(answer, "cox.solve") <- list(backend = if (args$solver == "eigen") "direct_wide" else args$solver,
    iterations = 0L, residual = NA_real_, converged = TRUE)
  answer
}

L0glm_ext_cox_baseline <- function(object) {
  eta <- as.numeric(predict(object, type = "link"))
  stats <- L0glm_ext_cox_eta(eta, object$response, object$weights[, 1L], object$cox.cache)
  hazard <- exp(log(rev(stats$increments)) - stats$shift)
  log_cumhaz <- log(cumsum(rev(stats$increments))) - stats$shift
  data.frame(time = rev(stats$time), hazard = hazard, cumhaz = exp(log_cumhaz), log_cumhaz = log_cumhaz)
}
