# Family adapters supply the working model and stable unit deviance; the
# iteration and backtracking rules are shared by all distributions.
L0glm_irls_deviance <- function(y, eta, weights, family) {
  good <- weights > 0
  eta <- as.numeric(eta)[good]; y <- as.numeric(y)[good]; weights <- as.numeric(weights)[good]
  if (any(!is.finite(eta)) || is.function(family$valideta) && !isTRUE(all(family$valideta(eta)))) return(Inf)
  mu <- suppressWarnings(family$linkinv(eta))
  if (any(!is.finite(mu))) return(Inf)
  if (family$family %in% c("poisson", "quasipoisson", "Gamma", "inverse.gaussian")) {
    if (any(mu < 0)) return(Inf)
    mu <- pmax(mu, .Machine$double.eps)
  }
  if (is.function(family$validmu) && !isTRUE(all(family$validmu(mu)))) return(Inf)
  terms <- if (family$family == "inverse.gaussian" && family$link == "log")
    weights * ((y * exp(-eta) - 1) / sqrt(y))^2 else
      suppressWarnings(family$dev.resids(y, mu, weights))
  value <- sum(terms)
  if (is.finite(value)) value else Inf
}

# A callback defines the fixed-nuisance penalized objective for this step.
# Optional eta values reuse the linear predictor during coefficient damping.
L0glm_step <- function(proposal, previous, objective, eta = NULL, eta_previous = NULL,
                      baseline = objective(previous, eta_previous),
                      maxit = 50L, rtol = 1e-10) {
  threshold <- if (is.finite(baseline)) baseline + rtol * max(1, abs(baseline)) else Inf
  beta <- proposal; predictor <- eta; fraction <- 1
  for (halving in 0:maxit) {
    value <- tryCatch(objective(beta, predictor), error = function(e) Inf)
    if (is.finite(value) && value <= threshold)
      return(list(beta = beta, eta = predictor, value = value, halvings = halving, accepted = TRUE))
    fraction <- fraction / 2
    beta <- previous + fraction * (proposal - previous)
    predictor <- if (is.null(eta)) NULL else eta_previous + fraction * (eta - eta_previous)
  }
  list(beta = previous, eta = eta_previous, halvings = maxit + 1L, accepted = FALSE)
}

L0glm_threshold_point <- function(beta, eta, X, scaling, thresh, always.include) {
  magnitude <- abs(beta) * scaling
  cutoff <- max(c(0, magnitude[setdiff(seq_along(beta), always.include)])) * thresh
  remove <- which((magnitude < cutoff | magnitude < 1e-20) &
    !seq_along(beta) %in% always.include & beta != 0)
  if (length(remove)) {
    eta <- eta - as.numeric(X[, remove, drop = FALSE] %*% beta[remove])
    beta[remove] <- 0
  }
  list(beta = beta, eta = eta, changed = length(remove) > 0L)
}
