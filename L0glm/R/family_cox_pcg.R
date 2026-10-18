# Exact Cox Hessian products, on the internally time-ordered working design.
L0glm_ext_cox_hvp <- function(v, X, stats, D = 0) {
  u <- as.numeric(X %*% v)
  risk <- cumsum(stats$e * u)[stats$ends]
  a <- ifelse(stats$S0 > 0, stats$d * risk/stats$S0^2, 0)
  back <- rev(cumsum(rev(a)))[stats$group]
  as.numeric(Matrix::crossprod(X, stats$q * u - stats$e * back)) + D * v
}

# Symmetrically precondition the system. Verify the true scaled residual before
# accepting convergence; recursive CG residuals alone can be misleading.
L0glm_ext_cox_cg <- function(operator, rhs, diagonal, tol, maxit) {
  scale <- sqrt(diagonal)
  apply <- function(v) operator(v/scale)/scale
  b <- rhs/scale; value <- numeric(length(b)); r <- b; direction <- r
  norm_b <- sqrt(sum(b^2)); rr <- sum(r^2)
  if (norm_b == 0) return(list(value = value, iterations = 0L, residual = 0, converged = TRUE))
  reason <- "iteration limit"
  for (iteration in seq_len(maxit)) {
    product <- apply(direction); curvature <- sum(direction * product)
    if (!is.finite(curvature) || curvature <= 0) { reason <- "nonpositive curvature"; break }
    alpha <- rr/curvature
    value <- value + alpha * direction; r <- r - alpha * product
    next_rr <- sum(r^2)
    if (sqrt(next_rr) <= tol * norm_b) {
      r <- b - apply(value); next_rr <- sum(r^2)
      if (sqrt(next_rr) <= tol * norm_b)
        return(list(value = value/scale, iterations = iteration,
          residual = sqrt(next_rr)/norm_b, converged = TRUE))
      direction <- r; rr <- next_rr
      next
    }
    direction <- r + (next_rr/rr) * direction; rr <- next_rr
  }
  list(value = value/scale, iterations = iteration,
    residual = sqrt(sum((b - apply(value))^2))/norm_b, converged = FALSE, reason = reason)
}

L0glm_ext_cox_pcg_step <- function(data, B, D, stats) {
  active <- is.finite(D[, 1L]); answer <- matrix(0, data$p, 1L)
  X <- data$Z[, active, drop = FALSE]; d <- D[active, 1L]
  beta <- B[active, 1L]
  rhs <- as.numeric(Matrix::crossprod(X, stats$score)) - d * (beta - data$prior[active, 1L])
  diagonal <- as.numeric(Matrix::crossprod(X^2, stats$q)) + d
  diagonal[diagonal <= 0] <- 1
  if (any(!is.finite(c(rhs, diagonal)))) stop("Nonfinite Cox ridge system; check predictor scaling.", call. = FALSE)
  solve <- L0glm_ext_cox_cg(function(v) L0glm_ext_cox_hvp(v, X, stats, d), rhs,
    diagonal, data$control$cox.cg.tol, data$control$cox.cg.maxit)
  answer[active, 1L] <- beta + solve$value
  solve$value <- NULL; solve$backend <- "pcg"
  attr(answer, "cox.solve") <- solve
  answer
}

# Same weak diagonal-working-weight ridge pilot as before. Use a matrix-free
# Woodbury system when wide, so numerical errors cannot enter its primal nullspace.
L0glm_ext_cox_ridge_pilot <- function(data) {
  stats <- L0glm_ext_cox_stats(data, matrix(0, data$p, 1L))
  X <- data$Z; ridge <- 1e-3; good <- stats$q > 0
  wide <- data$p > sum(good)
  if (wide) {
    X <- X[good, , drop = FALSE]; q <- stats$q[good]; root <- sqrt(q)
    rhs <- stats$score[good]/root
    diagonal <- 1 + q * as.numeric(Matrix::rowSums(X^2))/ridge
    operator <- function(v) v + root * as.numeric(X %*% Matrix::crossprod(X, root * v))/ridge
  } else {
    rhs <- as.numeric(Matrix::crossprod(X, stats$score))
    diagonal <- as.numeric(Matrix::crossprod(X^2, stats$q)) + ridge
    operator <- function(v) as.numeric(Matrix::crossprod(X, stats$q * as.numeric(X %*% v))) + ridge*v
  }
  solve <- L0glm_ext_cox_cg(operator, rhs, diagonal, data$control$cox.cg.tol, data$control$cox.cg.maxit)
  beta <- if (wide) as.numeric(Matrix::crossprod(X, root * solve$value))/ridge else solve$value
  solve$value <- NULL; solve$backend <- "pcg"; solve$space <- if (wide) "observations" else "coefficients"
  solve$dimension <- length(rhs)
  beta <- matrix(beta, data$p, 1L)
  attr(beta, "cox.solve") <- solve
  beta
}

L0glm_ext_cox_working <- function(data, B, D, state, args) {
  active <- is.finite(D[, 1L]); k <- sum(active)
  if (!k) {
    answer <- matrix(0, data$p, 1L)
    attr(answer, "cox.solve") <- list(backend = "empty", iterations = 0L, residual = 0, converged = TRUE)
    return(answer)
  }
  mode <- data$control$cox
  unrestricted <- args$solver == "eigen" && all(data$lower[active] < -1e100)
  if (mode == "pcg" && !unrestricted)
    stop('Cox cox="pcg" requires solver="eigen" and unrestricted coefficients.', call. = FALSE)
  stats <- L0glm_ext_cox_stats(data, B)
  cutoff <- data$control$cox.llt.max
  if (length(cutoff) > 1L) cutoff <- cutoff[if (inherits(data$Z, "sparseMatrix")) "sparse" else "dense"]
  # Avoid the wide eigensystem and large temporary risk-mean arrays in auto mode.
  direct_ok <- k <= data$n && k <= data$control$max.working && k * sum(stats$d > 0) <= 1e7
  small <- k <= cutoff && direct_ok
  if (mode == "direct" || !unrestricted || (mode == "auto" && small))
    return(L0glm_ext_cox_direct(data, B, D, stats, args))
  answer <- L0glm_ext_cox_pcg_step(data, B, D, stats)
  diagnostic <- attr(answer, "cox.solve")
  if (diagnostic$converged) return(answer)
  if (mode == "auto" && direct_ok) {
    warning(sprintf("Cox PCG did not converge in %d iterations (scaled residual %.3g); falling back to the direct solver.",
      diagnostic$iterations, diagnostic$residual), call. = FALSE)
    answer <- L0glm_ext_cox_direct(data, B, D, stats, args)
    attr(answer, "cox.solve")$backend <- "direct_fallback"
    attr(answer, "cox.solve")$iterations <- diagnostic$iterations
    return(answer)
  }
  stop(sprintf("Cox PCG did not converge in %d iterations (scaled residual %.3g; %s). Increase cox.cg.maxit or use cox='direct'.",
    diagnostic$iterations, diagnostic$residual, diagnostic$reason), call. = FALSE)
}
