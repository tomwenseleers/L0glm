L0glm_ext_multinomial_response <- function(y, family, n, levels = NULL) {
  scoring <- !is.null(levels)
  if (is.null(dim(y))) {
    if (length(y) != n || anyNA(y)) stop("Invalid multinomial class response.", call. = FALSE)
    if (is.null(levels)) levels <- if (is.factor(y)) base::levels(y) else {
      K <- family$nlp + 1L
      if (is.numeric(y) && all(y %in% 0:(K - 1L))) as.character(0:(K - 1L)) else
        if (is.numeric(y) && all(y %in% seq_len(K))) as.character(seq_len(K)) else sort(unique(as.character(y)))
    }
    index <- match(as.character(y), levels)
    if (anyNA(index)) stop("Unknown multinomial class.", call. = FALSE)
    Y <- matrix(0, n, length(levels)); Y[cbind(seq_len(n), index)] <- 1
  } else {
    Y <- as.matrix(y)
    if (nrow(Y) != n || any(!is.finite(Y) | Y < 0 | abs(Y - round(Y)) > 1e-7))
      stop("Multinomial count responses must be nonnegative integers.", call. = FALSE)
    if (is.null(levels)) levels <- if (is.null(colnames(Y))) as.character(seq_len(ncol(Y)) - 1L) else colnames(Y)
  }
  if (ncol(Y) < 2L || ncol(Y) != family$nlp + 1L || length(levels) != ncol(Y))
    stop("The response category count must equal mgcv::multinom(K)$K + 1.", call. = FALSE)
  if (!scoring && any(colSums(Y) == 0)) stop("Every multinomial category must occur in training data.", call. = FALSE)
  colnames(Y) <- levels
  list(y = Y, levels = levels, m = ncol(Y), coordinates = "clr")
}

L0glm_ext_softmax <- function(eta) {
  shifted <- eta - apply(eta, 1L, max)
  normalizer <- log(rowSums(exp(shifted)))
  logp <- shifted - normalizer
  list(p = exp(logp), logp = logp)
}

# Solve normal equations of a coupled IRLS weighted least-squares subproblem.
L0glm_ext_quadratic <- function(H, rhs, B, D, lower, prior, args) {
  active <- is.finite(D)
  answer <- numeric(length(B))
  if (!any(active)) return(answer)
  H <- H[active, active, drop = FALSE]
  d <- D[active]
  if (args$solver == "eigen") {
    coefficient <- eigen_ridge_solve(H + diag(d, length(d)),
      rhs[active] + d * prior[active], B[active], args$control.solver, 1L)
    if (is.null(coefficient)) stop("The coupled IRLS ridge solve failed.", call. = FALSE)
    answer[active] <- pmax(lower[active], coefficient)
  } else {
    eig <- eigen(H, symmetric = TRUE)
    keep <- eig$values > max(eig$values) * 1e-13
    Q <- t(eig$vectors[, keep, drop = FALSE])
    A <- sqrt(eig$values[keep]) * Q
    z <- as.numeric(Q %*% rhs[active])/sqrt(eig$values[keep])
    answer[active] <- L0glm_ext_ridge(A, z, NULL, d, B[active], lower[active], prior[active], args)
  }
  answer
}

L0glm_ext_multinomial_working <- function(data, B, D, state, args) {
  # Model, adaptive penalties and masks are CLR. Only the linear solve uses a
  # full-rank orthonormal basis of the zero-sum class plane.
  vars <- which(rowSums(is.finite(D)) > 1L)
  answer <- matrix(0, data$p, data$m)
  if (!length(vars)) return(answer)
  total <- rowSums(data$y); w <- data$weights[, 1L]; good <- w * total > 0
  X <- data$Z[good, vars, drop = FALSE]; n <- nrow(X); p <- ncol(X)
  C <- data$m; m <- C - 1L; Q <- L0glm_ext_contrasts(C)
  eta <- as.matrix(data$Z %*% B)
  P <- L0glm_ext_softmax(eta + data$offset)$p[good, , drop = FALSE]
  score <- (data$y[good, , drop = FALSE] - total[good] * P) * w[good]
  weight <- w[good] * total[good]
  d <- D[vars, , drop = FALSE]
  shared <- all(is.finite(d)) && all(d == d[, 1L])
  if (shared) {
    start <- as.numeric(B[vars, , drop = FALSE] %*% Q)
    prior <- as.numeric(data$prior[vars, , drop = FALSE] %*% Q)
    penalty <- rep(d[, 1L], m)
    map <- NULL; dimension <- p * m
  } else {
    basis <- L0glm_ext_multinomial_basis(d, B[vars, , drop = FALSE], data$prior[vars, , drop = FALSE])
    map <- basis$map; start <- basis$start; prior <- basis$prior
    penalty <- basis$penalty; dimension <- ncol(map)
  }
  if (min(dimension, n * m) > data$control$max.working)
    stop("Multinomial working system exceeds control.family$max.working.", call. = FALSE)
  bound <- data$control$multinomial == "bound"
  M <- P %*% Q
  if (dimension <= n * m) {
    if (shared) {
      H <- matrix(0, dimension, dimension)
      for (j in seq_len(m)) for (k in seq_len(j)) {
        value <- if (bound) rep(if (j == k) .5 else 0, n) else
          as.numeric(P %*% (Q[, j] * Q[, k])) - M[, j] * M[, k]
        block <- as.matrix(Matrix::crossprod(X, Matrix::Diagonal(x = weight * value) %*% X))
        a <- seq_len(p) + (j - 1L) * p; b <- seq_len(p) + (k - 1L) * p
        H[a, b] <- block; H[b, a] <- t(block)
      }
      gradient <- as.numeric(Matrix::crossprod(X, score %*% Q))
    } else {
      H <- matrix(0, dimension, dimension)
      for (j in seq_len(C)) for (k in seq_len(C)) {
        value <- if (bound) rep(.5 * ((j == k) - 1/C), n) else P[, j] * ((j == k) - P[, k])
        block <- Matrix::crossprod(X, Matrix::Diagonal(x = weight * value) %*% X)
        a <- seq_len(p) + (j - 1L) * p; b <- seq_len(p) + (k - 1L) * p
        H <- H + as.matrix(Matrix::crossprod(map[a, , drop = FALSE], block %*% map[b, , drop = FALSE]))
      }
      gradient <- as.numeric(Matrix::crossprod(map, as.numeric(Matrix::crossprod(X, score))))
    }
    coefficient <- L0glm_ext_quadratic(H, H %*% start + gradient, start, penalty,
      rep(-.Machine$double.xmax, dimension), prior, args)
  } else {
    # Whiten in the C-1 observation contrasts. The ridge kernel then uses
    # Woodbury without ever forming a pC-by-pC coefficient Hessian.
    roots <- vector("list", n); z <- numeric(n * m)
    linear <- eta[good, , drop = FALSE] %*% Q
    for (i in seq_len(n)) {
      v <- pmax(P[i, ], 1e-12); v <- v/sum(v); mean <- as.numeric(crossprod(Q, v))
      V <- if (bound) diag(.5, m) else crossprod(Q, Q * v) - tcrossprod(mean)
      R <- chol(weight[i] * V); roots[[i]] <- R
      rows <- i + (seq_len(m) - 1L) * n
      z[rows] <- as.numeric(R %*% linear[i, ]) + as.numeric(forwardsolve(t(R), as.numeric(score[i, ] %*% Q)))
    }
    indices <- which(matrix(TRUE, m, m), arr.ind = TRUE)
    ii <- unlist(lapply(seq_len(n), function(i) i + (indices[, 1L] - 1L) * n))
    jj <- unlist(lapply(seq_len(n), function(i) i + (indices[, 2L] - 1L) * n))
    root <- Matrix::sparseMatrix(i = ii, j = jj, x = unlist(roots), dims = c(n*m, n*m))
    if (shared) design <- if (inherits(X, "sparseMatrix")) Matrix::kronecker(Matrix::Diagonal(m), X) else kronecker(diag(m), as.matrix(X)) else {
      transform <- Matrix::kronecker(Matrix::Matrix(t(Q), sparse = TRUE), Matrix::Diagonal(p)) %*% map
      design <- if (inherits(X, "sparseMatrix")) Matrix::kronecker(Matrix::Diagonal(m), X) %*% transform else kronecker(diag(m), as.matrix(X)) %*% transform
    }
    coefficient <- L0glm_ext_ridge(root %*% design, z, NULL, penalty, start,
      rep(-.Machine$double.xmax, dimension), prior, args)
  }
  answer[vars, ] <- if (shared) matrix(coefficient, p, m) %*% t(Q) else matrix(map %*% coefficient, p, C)
  answer
}

L0glm_ext_contrasts <- function(C) {
  Q <- stats::contr.helmert(C)
  sweep(Q, 2L, sqrt(colSums(Q^2)), "/")
}

# Entrywise CLR sparsity: remove zero class coefficients, then parameterize the
# remaining zero-sum plane separately for each predictor. Whitening its ridge
# penalty retains the existing diagonal-ridge/Woodbury solver interface.
L0glm_ext_multinomial_basis <- function(D, B, prior) {
  p <- nrow(D); C <- ncol(D); ii <- jj <- xx <- start <- center <- penalty <- numeric()
  count <- 0L
  for (j in seq_len(p)) {
    keep <- which(is.finite(D[j, ])); s <- length(keep)
    if (s < 2L) next
    Q <- L0glm_ext_contrasts(s); d <- D[j, keep]
    unpenalized <- all(d == 0)
    R <- if (unpenalized) diag(s - 1L) else chol(crossprod(Q, Q * d))
    T <- Q %*% backsolve(R, diag(s - 1L))
    entries <- which(matrix(TRUE, s, s - 1L), arr.ind = TRUE)
    ii <- c(ii, j + (keep[entries[, 1L]] - 1L) * p)
    jj <- c(jj, count + entries[, 2L]); xx <- c(xx, as.numeric(T))
    start <- c(start, as.numeric(R %*% crossprod(Q, B[j, keep])))
    center <- c(center, if (unpenalized) numeric(s - 1L) else as.numeric(crossprod(T, d * prior[j, keep])))
    penalty <- c(penalty, rep(if (unpenalized) 0 else 1, s - 1L))
    count <- count + s - 1L
  }
  list(map = Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(p*C, count)),
    start = start, prior = center, penalty = penalty)
}
