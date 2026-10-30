# Stable working models used both by IRLS and restored validation checkpoints.
L0glm_working_model <- function(family) {
  if (is.function(family$score)) return(function(y, eta, mu, weights, offset) {
    eta <- as.numeric(eta); weights <- as.numeric(weights); offset <- as.numeric(offset)
    d <- family$score(y, mu); derivative <- family$mu.eta(eta)
    W <- weights * d[, 2L] * derivative^2
    W[weights == 0] <- 0
    z <- eta - offset - d[, 1L] / (d[, 2L] * derivative)
    list(weights = W, response = z,
      good = which(W > 0 & is.finite(W) & is.finite(z) & derivative != 0))
  })
  if (family$link == "identity") return(function(y, eta, mu, weights, offset) {
    W <- as.numeric(weights / family$variance(mu))
    W[weights == 0] <- 0
    z <- y - offset
    list(weights = W, response = z, good = which(W > 0 & is.finite(W) & is.finite(z)))
  })
  if (family$family == "inverse.gaussian" && family$link == "log")
    return(function(y, eta, mu, weights, offset) {
      W <- as.numeric(weights) * exp(-as.numeric(eta))
      W[weights == 0] <- 0
      z <- as.numeric(eta - offset) - 1 + as.numeric(y / mu)
      list(weights = W, response = z, good = which(W > 0 & is.finite(W) & is.finite(z)))
    })
  function(y, eta, mu, weights, offset) {
    derivative <- family$mu.eta(eta)
    W <- weights * as.numeric(derivative^2 / family$variance(mu))
    W[weights == 0] <- 0
    z <- eta - offset + (y - mu) / derivative
    list(weights = W, response = z,
      good = which(W > 0 & is.finite(W) & is.finite(z) & derivative != 0))
  }
}
