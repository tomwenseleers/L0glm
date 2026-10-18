L0glm_ext_nb_response <- function(y, family, n) {
  y <- as.matrix(y)
  if (nrow(y) != n || any(!is.finite(y) | y < 0 | abs(y - round(y)) > 1e-7))
    stop("Negative-binomial responses must be nonnegative integer counts.", call. = FALSE)
  theta <- as.numeric(family$getTheta(TRUE))
  list(y = y, m = ncol(y), theta = rep(theta, ncol(y)),
       estimated = rep(isTRUE(family$n.theta > 0), ncol(y)))
}

L0glm_ext_nb_mean <- function(data, B) {
  eta <- as.matrix(data$Z %*% B + data$offset)
  mu <- matrix(data$family$linkinv(eta), nrow(eta), ncol(eta))
  list(eta = eta, mu = mu)
}

L0glm_ext_nb_working <- function(data, B, D, state, args) {
  mean <- L0glm_ext_nb_mean(data, B)
  answer <- B
  for (j in seq_len(data$m)) {
    eta <- mean$eta[, j]
    mu <- mean$mu[, j]
    if (any(!is.finite(mu) | mu <= 0)) {
      mu <- data$y[, j] + 0.1
      eta <- data$family$linkfun(mu)
    }
    derivative <- data$family$mu.eta(eta)
    W <- data$weights[, j] * derivative^2 / (mu + mu^2/state$theta[j])
    good <- is.finite(W) & W > 0 & derivative != 0
    z <- eta - data$offset[, j] + (data$y[, j] - mu)/derivative
    answer[, j] <- L0glm_ext_ridge(data$Z[good, , drop = FALSE], z[good], W[good],
      D[, j], B[, j], data$lower, data$prior[, j], args)
  }
  answer
}

L0glm_ext_theta_derivatives <- function(y, mu, theta, w) {
  U <- sum(w * (digamma(y + theta) - digamma(theta) + log(theta) + 1 -
    log(theta + mu) - (y + theta)/(theta + mu)))
  H <- sum(w * (trigamma(y + theta) - trigamma(theta) + 1/theta -
    2/(theta + mu) + (y + theta)/(theta + mu)^2))
  c(gradient = theta * U, curvature = theta^2 * H + theta * U)
}

L0glm_ext_nb_update <- function(data, B, state) {
  mu <- L0glm_ext_nb_mean(data, B)$mu
  for (j in which(state$estimated)) {
    good <- data$weights[, j] > 0
    y <- data$y[good, j]; fitted <- mu[good, j]; w <- data$weights[good, j]
    if (any(!is.finite(fitted) | fitted <= 0)) next
    a <- log(state$theta[j])
    derivatives <- L0glm_ext_theta_derivatives(y, fitted, exp(a), w)
    step <- if (is.finite(derivatives[2L]) && derivatives[2L] < 0)
      -derivatives[1L]/derivatives[2L] else sign(derivatives[1L])
    step <- max(-1, min(1, step))
    likelihood <- function(value) sum(w * stats::dnbinom(y, mu = fitted, size = exp(value), log = TRUE))
    baseline <- likelihood(a)
    for (halving in 0:25) {
      candidate <- max(log(data$control$theta.min), min(log(data$control$theta.max), a + step))
      if (is.finite(likelihood(candidate)) && likelihood(candidate) >= baseline) {
        state$theta[j] <- exp(candidate)
        break
      }
      step <- step/2
    }
  }
  state
}
