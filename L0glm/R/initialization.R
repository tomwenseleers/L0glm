L0glm_nuisance_solve <- function(G, rhs) {
  answer <- qr.coef(qr(as.matrix(G)), as.matrix(rhs))
  answer[is.na(answer)] <- 0
  answer
}

# Original-coordinate pilot for nonnegative identity-link signals above a
# compulsory background or known offset. No ridge shrinkage is applied to this
# pilot: a small initial amplitude otherwise creates a large adaptive penalty.
L0glm_signal_family <- function(family) family$link == "identity" &&
  !family$family %in% c("gaussian", "binomial", "quasibinomial")

L0glm_signal_start <- function(X, y, offset, weights, included = integer(), prior = NULL) {
  values <- if (inherits(X, "sparseMatrix") && "x" %in% methods::slotNames(X)) X@x else X
  if (any(!is.finite(values)) || any(values < 0) ||
      any(!is.finite(offset)) || any(!is.finite(weights)) || any(weights < 0)) return(NULL)
  y <- as.matrix(y); n <- nrow(X); m <- ncol(y); p <- ncol(X)
  selectable <- setdiff(seq_len(p), included)
  if (!length(selectable) || !any(weights > 0)) return(NULL)
  residual <- pmax(y - offset, 0)
  result <- matrix(0, p, m)
  if (length(included)) {
    B <- as.matrix(X[, included, drop = FALSE]); mass <- rowSums(B)
    if (!any(mass > 0)) return(NULL)
    constant <- max(abs(mass - mass[1L])) <= 1e-12 * max(1, mass[1L]) && mass[1L] > 0
    for (j in seq_len(m)) {
      good <- weights[, j] > 0
      level <- if (any(good)) as.numeric(stats::quantile(residual[good, j], 0.05)) else 0
      result[included, j] <- if (constant) level / mass[1L] else
        as.numeric(nuisance_quadratic_cpp(crossprod(B, B * weights[, j]),
          as.numeric(crossprod(B, weights[, j] * level)), rep(TRUE, ncol(B))))
    }
    residual <- pmax(residual - B %*% result[included, , drop = FALSE], 0)
  }
  A <- X[, selectable, drop = FALSE]
  information <- as.matrix(Matrix::crossprod(A^2, weights))
  rhs <- as.matrix(Matrix::crossprod(A, weights * residual))
  floor <- sweep(1 / sqrt(pmax(information, .Machine$double.xmin)), 2L,
    1e-8 * sqrt(colSums(weights * residual^2)), "*")
  beta <- pmax(rhs / pmax(information, .Machine$double.xmin), floor)
  if (!is.null(prior)) beta <- pmax(beta, prior[selectable, , drop = FALSE])
  beta[information == 0] <- 0
  result[selectable, ] <- beta
  result
}

L0glm_initial_coefficients <- function(X, y, mu, weights, offset, prior.mean,
                                        family, lambda, always.include, lower,
                                        preprocessing) {
  if (preprocessing$mode == "none")
    return(init_coefs(X, y, mu, weights, offset, prior.mean, family, lambda))
  p <- ncol(X)
  if (!p) return(numeric())
  eta <- family$linkfun(mu)
  derivative <- family$mu.eta(eta)
  W <- as.numeric(weights * derivative^2 / family$variance(mu))
  z <- if (family$link == "identity") as.numeric(y - offset) else
    as.numeric(eta - offset + (y - mu) / derivative)
  if (is.function(family$score)) {
    working <- L0glm_working_model(family)(y, eta, mu, weights, offset)
    W <- working$weights; z <- working$response
  }
  good <- is.finite(W) & W > 0 & is.finite(z)
  if (!any(good)) return(rep(0, p))
  X <- X[good, , drop = FALSE]; W <- W[good]
  A <- Xweights(X, W); b <- z[good] * sqrt(W)
  moments <- L0glm_moments(X, W)
  info <- moments$ss
  free <- intersect(always.include, which(lower < 0))
  selectable <- setdiff(seq_len(p), free)
  U <- A[, free, drop = FALSE]
  response <- b
  if (length(free)) {
    G <- Matrix::crossprod(U)
    C <- Matrix::crossprod(U, A)
    response <- b - as.numeric(U %*% L0glm_nuisance_solve(G, Matrix::crossprod(U, b)))
    info <- pmax(0, info - colSums(as.matrix(C) * L0glm_nuisance_solve(G, C)))
    if (length(free) == 1L && free %in% L0glm_unit_intercepts(X))
      info <- moments$centered_ss
  }
  rhs <- as.numeric(Matrix::crossprod(A, response))
  penalty <- rep_len(lambda, p) * preprocessing$scale^2
  penalty[always.include] <- 0
  beta <- ifelse(info + penalty > 0,
    (rhs + penalty * as.numeric(prior.mean)) / (info + penalty), 0)
  beta[info == 0 & as.numeric(prior.mean) == 0] <- 0
  beta[free] <- 0
  if (family$family == "gaussian") {
    prediction <- as.numeric(A %*% beta)
    if (length(free)) prediction <- prediction -
      as.numeric(U %*% L0glm_nuisance_solve(G, Matrix::crossprod(U, prediction)))
    factor <- if (sum(prediction^2) > 0) sum(prediction * response) / sum(prediction^2) else 0
    beta[selectable] <- beta[selectable] * factor
  }
  if (length(free)) beta[free] <- as.numeric(L0glm_nuisance_solve(G,
    Matrix::crossprod(U, b - A %*% beta)))
  beta
}

L0glm_coefficient_change <- function(beta, previous, scaling, included,
                                    prediction, previous_prediction, weights) {
  selectable <- setdiff(seq_along(beta), included)
  difference <- (beta - previous) * scaling
  slope_change <- if (length(selectable)) sqrt(sum(difference[selectable]^2)) /
    max(1e-12, sqrt(sum((beta[selectable] * scaling[selectable])^2))) else 0
  if (!length(included)) return(slope_change)
  # Changes of input origin or response location must not hide slope changes.
  prediction_change <- prediction - previous_prediction
  total <- sum(weights)
  centered <- prediction - if (total > 0) sum(weights * prediction) / total else 0
  fitted_change <- sqrt(sum(weights * prediction_change^2)) /
    max(sqrt(total), sqrt(sum(weights * centered^2)), 1e-12)
  max(slope_change, fitted_change)
}
