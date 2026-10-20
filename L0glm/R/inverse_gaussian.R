# Stable unit deviance for a log-link inverse Gaussian. The usual
# (y-mu)^2/(y*mu^2) can produce Inf/Inf even when its limit is finite.
L0glm_ig_deviance <- function(y, eta, weights) {
  good <- weights > 0
  if (any(!is.finite(eta))) return(Inf)
  inverse <- exp(-as.numeric(eta))
  if (any(!is.finite(inverse) | inverse <= 0) ||
      any(!is.finite(exp(eta)))) return(Inf)
  sum(weights[good] * ((y[good] * inverse[good] - 1) / sqrt(y[good]))^2)
}

# Backtrack only the coefficient/predictor step: no further linear solve.
# Frozen ridge weights define the penalized-likelihood target for this update.
L0glm_ig_step <- function(beta, previous, eta, eta_previous, y, weights,
                          lambdas, prior.mean, X, scaling, thresh,
                          hard_threshold, always.include, check_decrease = TRUE) {
  penalized <- which(is.finite(lambdas) & lambdas > 0)
  score <- function(b, e) L0glm_ig_deviance(y, e, weights) +
    sum(lambdas[penalized] * (b[penalized] - prior.mean[penalized])^2)
  reference <- if (check_decrease) score(previous, eta_previous) else Inf
  tolerance <- 32 * .Machine$double.eps * max(1, abs(reference))
  proposed <- beta
  eta_proposed <- as.numeric(eta)
  step <- 1
  for (halving in 0:50) {
    value <- score(beta, eta)
    if (is.finite(value) && value <= reference + tolerance)
      return(list(beta=beta, eta=as.numeric(eta), halvings=halving, accepted=TRUE))
    step <- step/2
    beta <- previous + step * (proposed - previous)
    eta <- eta_previous + step * (eta_proposed - eta_previous)
    if (hard_threshold) {
      magnitude <- abs(beta) * scaling
      cutoff <- max(c(0, magnitude[setdiff(seq_along(beta), always.include)])) * thresh
      remove <- which((magnitude < cutoff | magnitude < 1e-20) &
        !seq_along(beta) %in% always.include & beta != 0)
      if (length(remove)) {
        eta <- eta - as.numeric(X[, remove, drop=FALSE] %*% beta[remove])
        beta[remove] <- 0
      }
    }
  }
  list(beta=previous, eta=eta_previous, halvings=51L, accepted=FALSE)
}
