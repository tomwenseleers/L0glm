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
L0glm_initial_dispersion <- function(X, y, mu, weights, offset, family) {
  eta <- family$linkfun(mu)
  derivative <- family$mu.eta(eta)
  W <- as.numeric(weights * derivative^2 / family$variance(mu))
  z <- as.numeric(eta - offset + (y - mu) / derivative)
  good <- is.finite(W) & W > 0 & is.finite(z)
  if (!any(good)) return(1)
  A <- Matrix::Diagonal(x = sqrt(W[good])) %*% X[good, , drop = FALSE]
  if (!inherits(A, "sparseMatrix") &&
      Matrix::nnzero(A) <= 0.1 * as.double(nrow(A)) * ncol(A))
    A <- Matrix::Matrix(A, sparse = TRUE)
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
    G <- if (dual) Matrix::tcrossprod(A) else Matrix::crossprod(A)
    K <- G + Matrix::Diagonal(nrow(G), tau)
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
  reference <- max(1, preset.lambda("gic", p = ncol(X), n = nrow(X)))
  if (!is.finite(phi) || phi <= 0) phi <- .Machine$double.eps
  max(.Machine$double.xmin, phi / reference)
}
