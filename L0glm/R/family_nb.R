L0glm_ext_nb_response <- function(y, family, n) {
  y <- as.matrix(y)
  if (nrow(y) != n || any(!is.finite(y) | y < 0 | abs(y - round(y)) > 1e-7))
    stop("Negative-binomial responses must be nonnegative integer counts.", call. = FALSE)
  theta <- as.numeric(family$getTheta(TRUE))
  list(y = y, m = ncol(y), theta = rep(theta, ncol(y)),
       estimated = rep(isTRUE(family$n.theta > 0), ncol(y)))
}

# Current/trial predictors are reused only within an immutable group-fit design.
L0glm_ext_nb_state <- function(data,control) {
  state <- list(theta=data$response$theta,estimated=data$response$estimated)
  theta <- control$.theta.start
  if(!is.null(theta)) {
    if(length(theta)!=data$m || any(!is.finite(theta)|theta<=0))
      stop("Invalid NB nuisance starting shape.",call.=FALSE)
    state$theta <- as.numeric(theta)
  }
  state
}

L0glm_ext_eta <- function(data, B) {
  cache <- data$prediction.cache
  eta <- if(!is.null(cache) && identical(B,cache$first.B))cache$first.eta else
    if(!is.null(cache) && identical(B,cache$second.B))cache$second.eta else NULL
  if(is.null(eta)) {
    eta <- as.matrix(data$Z %*% B + L0glm_observation_broadcast(data$offset))
    if(!is.null(cache)) {
      cache$second.B <- cache$first.B;cache$second.eta <- cache$first.eta
      cache$first.B <- B;cache$first.eta <- eta
    }
  }
  eta
}

L0glm_ext_nb_mean <- function(data, B) {
  eta <- L0glm_ext_eta(data,B)
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
    W <- L0glm_observation_column(data$weights,j) * derivative^2 / (mu + mu^2/state$theta[j])
    good <- is.finite(W) & W > 0 & derivative != 0
    z <- eta - L0glm_observation_column(data$offset,j) + (data$y[, j] - mu)/derivative
    answer[, j] <- L0glm_ext_ridge(data$Z[good, , drop = FALSE], z[good], W[good],
      D[, j], B[, j], data$lower, data$prior[, j], args)
  }
  answer
}

# A projected Fisher-diagonal step is a feasible descent direction when the
# clipped working ridge step stalls on a constrained relaxed support.
L0glm_nb_projected_step <- function(data, B, D, state) {
  prediction <- L0glm_ext_nb_mean(data,B)
  if (any(!is.finite(prediction$mu) | prediction$mu <= 0)) return(B)
  proposal <- B
  for (j in seq_len(data$m)) {
    active <- which(is.finite(D[,j]))
    if (!length(active)) next
    X <- data$Z[,active,drop=FALSE]
    mu <- prediction$mu[,j]
    derivative <- data$family$mu.eta(prediction$eta[,j])
    w <- L0glm_observation_column(data$weights,j)/(mu+mu^2/state$theta[j])
    ridge <- D[active,j]
    gradient <- as.numeric(Matrix::crossprod(X,w*derivative*(mu-data$y[,j])))+
      ridge*(B[active,j]-data$prior[active,j])
    diagonal <- as.numeric(Matrix::crossprod(X^2,w*derivative^2))+ridge
    good <- is.finite(diagonal) & diagonal > 0 & is.finite(gradient)
    keep <- active[good]
    proposal[keep,j] <- pmax(data$lower[keep],B[keep,j]-gradient[good]/diagonal[good])
  }
  proposal
}

L0glm_nb_refit_step <- function(data, proposal, B, D, state, objective) {
  old <- objective(B)
  search <- function(trial,maxit) L0glm_step(trial,B,function(candidate,eta) objective(candidate),
    baseline=old,maxit=maxit,rtol=1e-12)
  step <- search(proposal,10L)
  if (!is.finite(old) || (step$accepted && step$value < old-1e-12*max(1,abs(old))))
    return(step)
  fallback <- search(L0glm_nb_projected_step(data,B,D,state),18L)
  fallback$halvings <- fallback$halvings + step$halvings
  if (fallback$accepted) return(fallback)
  list(beta=B,halvings=fallback$halvings,accepted=TRUE)
}
