# Normalize the public NULL / fixed value / negative starting value convention.
L0glm_dispersion <- function(dispersion, family, m = 1L) {
  if (!is.null(dispersion) &&
      (!is.numeric(dispersion) || !length(dispersion) ||
       !(length(dispersion) %in% c(1L, m)) ||
       any(!is.finite(dispersion) | dispersion == 0))) {
    stop("dispersion must be NULL or finite nonzero numeric values of length 1 or ncol(y).",
         call. = FALSE)
  }
  fixed <- family$dispersion
  if (family$family %in% c("poisson", "binomial") ||
      startsWith(family$family, "Negative Binomial(")) fixed <- 1
  if (length(fixed) == 1L && is.finite(fixed) && fixed > 0) {
    if (!is.null(dispersion) && any(dispersion != fixed))
      warning("Dispersion is fixed at ", fixed, " for family ", family$family,
              "; dispersion was ignored.", call. = FALSE)
    return(rep(fixed, m))
  }
  if (is.null(dispersion)) NULL else rep_len(as.numeric(dispersion), m)
}

# Evaluate Pearson residuals at the current mean, not at the previous IRLS weights.
L0glm_pearson_scale <- function(y, mu, weights, family, k) {
  good <- is.finite(weights) & weights > 0
  df <- sum(good) - k
  if (df <= 0) return(NA_real_)
  sum(weights[good] * (y[good] - mu[good])^2 /
        family$variance(mu[good])) / df
}

# One weak ridge pilot of the initial working model. Four probe RHS reuse its
# factorization to estimate tr((I-H)'(I-H)), the residual noise degrees of freedom.
# This is an initialization heuristic: omitted-signal bias is not removed by df.
L0glm_initial_dispersion <- function(X, y, mu, weights, offset, family,
                                    always.include = integer(), preprocessing = NULL) {
  if (!is.null(preprocessing) && preprocessing$mode != "none")
    return(L0glm_canonical_dispersion(X, y, mu, weights, offset, family,
                                     always.include, preprocessing$p))
  L0glm_initial_dispersion_legacy(X, y, mu, weights, offset, family)
}

L0glm_initial_dispersion_legacy <- function(X, y, mu, weights, offset, family) {
  eta <- family$linkfun(mu)
  derivative <- family$mu.eta(eta)
  W <- as.numeric(weights * derivative^2 / family$variance(mu))
  z <- as.numeric(eta - offset + (y - mu) / derivative)
  good <- is.finite(W) & W > 0 & is.finite(z)
  if (!any(good)) return(1)
  A <- Matrix::Diagonal(x = sqrt(W[good])) %*% X[good, , drop = FALSE]
  info <- as.numeric(Matrix::colSums(A^2))
  use <- is.finite(info) & info > 0
  A <- A[, use, drop = FALSE] %*% Matrix::Diagonal(x = 1 / sqrt(info[use]))
  response <- z[good] * sqrt(W[good])
  phi <- mean(response^2)
  if (ncol(A)) {
    old_seed <- if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
      get(".Random.seed", .GlobalEnv) else NULL
    on.exit(if (is.null(old_seed)) rm(".Random.seed", envir = .GlobalEnv) else
      assign(".Random.seed", old_seed, envir = .GlobalEnv))
    set.seed(1701)
    probes <- matrix(sample(c(-1, 1), nrow(A) * 4L, replace = TRUE), nrow(A))
    B <- cbind(response, probes)
    tau <- 1e-6
    dual <- ncol(A) > nrow(A)
    K <- if (dual) Matrix::tcrossprod(A) + Matrix::Diagonal(nrow(A), tau) else
      L0glm_pilot_gram(A, rep(tau, ncol(A)))
    rhs <- if (dual) B else Matrix::crossprod(A, B)
    if (inherits(K, "sparseMatrix")) {
      factor <- Matrix::Cholesky(Matrix::forceSymmetric(K), LDL = FALSE)
      solution <- Matrix::solve(factor, rhs, system = "A")
    } else {
      U <- chol(as.matrix(K))
      solution <- backsolve(U, forwardsolve(t(U), as.matrix(rhs)))
    }
    # The dual formula avoids subtracting almost equal fitted values when p > n.
    residual <- if (dual) tau * solution else B - A %*% solution
    df <- mean(Matrix::colSums(residual[, -1L, drop = FALSE]^2))
    if (is.finite(df) && df > 0) phi <- sum(residual[, 1L]^2) / df
  }
  # Keep the permissive start independent of the fitted penalty. Dividing by
  # that penalty would cancel increasing lambda in lambda * dispersion / 4.
  reference <- max(1, preset.lambda("gic", p = ncol(X), n = nrow(X)), na.rm = TRUE)
  if (!is.finite(phi) || phi <= 0) phi <- .Machine$double.eps
  max(.Machine$double.xmin, phi / reference)
}

# Profile mandatory terms in the pilot as in the fitted model. Conditional
# column norms make this pilot independent of predictor units and, with a free
# intercept, of their origin. The four probes share one matrix factorization.
L0glm_canonical_dispersion <- function(X, y, mu, weights, offset, family,
                                      always.include, p) {
  eta <- family$linkfun(mu)
  derivative <- family$mu.eta(eta)
  W <- as.numeric(weights * derivative^2 / family$variance(mu))
  z <- as.numeric(eta - offset + (y - mu) / derivative)
  good <- is.finite(W) & W > 0 & is.finite(z)
  if (!any(good)) return(1)
  X <- X[good, , drop = FALSE]; W <- W[good]
  A <- Xweights(X, W)
  free <- intersect(always.include, seq_len(ncol(X)))
  if (length(free)) {
    nuisance_qr <- qr(as.matrix(A[, free, drop = FALSE]))
    free <- free[nuisance_qr$pivot[seq_len(nuisance_qr$rank)]]
  }
  use <- setdiff(seq_len(ncol(X)), always.include)
  U <- A[, free, drop = FALSE]
  moments <- L0glm_moments(X, W)
  info <- moments$ss
  if (length(free)) {
    G <- Matrix::crossprod(U)
    C <- Matrix::crossprod(U, A)
    info <- pmax(0, info - colSums(as.matrix(C) * L0glm_nuisance_solve(G, C)))
    if (length(free) == 1L && free %in% L0glm_unit_intercepts(X))
      info <- moments$centered_ss
  }
  use <- use[is.finite(info[use]) & info[use] > 0]
  A <- A[, use, drop = FALSE] %*% Matrix::Diagonal(x = 1 / sqrt(info[use]))
  # Dense designs can be residualized directly; sparse designs keep their
  # stored entries and use the nuisance block in the factorization instead.
  if (length(free) && !inherits(A, "sparseMatrix"))
    A <- A - U %*% L0glm_nuisance_solve(G, Matrix::crossprod(U, A))
  old_seed <- if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
    get(".Random.seed", .GlobalEnv) else NULL
  on.exit(if (is.null(old_seed)) rm(".Random.seed", envir = .GlobalEnv) else
    assign(".Random.seed", old_seed, envir = .GlobalEnv))
  set.seed(1701)
  B <- cbind(z[good] * sqrt(W),
             matrix(sample(c(-1, 1), nrow(A) * 4L, replace = TRUE), nrow(A)))
  if (length(free)) B <- B - U %*%
    L0glm_nuisance_solve(G, Matrix::crossprod(U, B))
  make_solve <- function(K) {
    if (inherits(K, "sparseMatrix")) {
      factor <- Matrix::Cholesky(Matrix::forceSymmetric(K), LDL = FALSE)
      function(rhs) Matrix::solve(factor, rhs, system = "A")
    } else {
      factor <- chol(as.matrix(K))
      function(rhs) backsolve(factor, forwardsolve(t(factor), as.matrix(rhs)))
    }
  }
  target <- min(10, (nrow(A) - length(free)) / 5)
  pilot_scale <- function(tau) {
    if (!ncol(A)) {
      residual_operator <- identity
    } else if (ncol(A) > nrow(A)) {
      K <- Matrix::tcrossprod(A) + Matrix::Diagonal(nrow(A), tau)
      solve <- make_solve(K)
      inverse_U <- if (length(free)) solve(U) else NULL
      residual_operator <- function(rhs) {
        value <- solve(rhs)
        if (length(free)) value <- value - inverse_U %*% L0glm_nuisance_solve(
          Matrix::crossprod(U, inverse_U), Matrix::crossprod(U, value))
        tau * value
      }
    } else {
      design <- cbind(U, A)
      K <- L0glm_pilot_gram(design, c(rep(0, length(free)), rep(tau, ncol(A))))
      solve <- make_solve(K)
      residual_operator <- function(rhs) rhs - design %*% solve(Matrix::crossprod(design, rhs))
    }
    residual <- residual_operator(B)
    df2 <- mean(Matrix::colSums(residual[, -1L, drop = FALSE]^2))
    # For R = I-H, nu = tr(R^2)^2 / tr(R^4) measures the number of
    # independent residual directions. Since 0 <= R <= I, nu >= tr(R^2):
    # ordinary pilots pass this cheap bound without any additional solves.
    effective_df <- df2
    if (!is.finite(df2) || df2 < target) {
      repeated <- residual_operator(residual[, -1L, drop = FALSE])
      df4 <- mean(Matrix::colSums(repeated^2))
      effective_df <- if (is.finite(df4) && df4 > 0) df2^2 / df4 else 0
    }
    list(phi = if (is.finite(df2) && df2 > 0) sum(residual[, 1L]^2) / df2 else NA_real_,
         effective_df = effective_df)
  }
  value <- pilot_scale(1e-6)
  null_scale <- sum(B[, 1L]^2) / max(1, nrow(A) - length(free))
  # Near interpolation, a few residual directions can give a very noisy scale
  # even when its magnitude appears plausible. Retry once with moderate ridge.
  if (!is.finite(value$phi) || value$phi < 1e-6 * null_scale || value$effective_df < target)
    value <- pilot_scale(0.1)
  phi <- value$phi
  reference <- max(1, preset.lambda("gic", p = max(1, p), n = nrow(X)), na.rm = TRUE)
  if (!is.finite(phi) || phi <= 0) phi <- .Machine$double.eps
  max(.Machine$double.xmin, phi / reference)
}

# A large dense-stored PSF design can contain mostly exact zeros. Use the
# existing sparse row-augmentation path for this auxiliary ridge factorization;
# the original design and main solver storage remain unchanged. The system is
# exactly crossprod(A) + diag(penalty); only its construction/storage differ.
L0glm_pilot_gram <- function(A, penalty) {
  if (!inherits(A, "sparseMatrix") && length(A) >= 1e6 &&
      Matrix::nnzero(A) < 0.1 * length(A))
    return(Matrix::crossprod(row_augX(A, penalty)))
  Matrix::crossprod(A) + Matrix::Diagonal(ncol(A), penalty)
}
