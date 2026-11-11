# A user's explicit starting coefficients, scale or penalty take precedence.
L0glm_pilot_enabled <- function(control) {
  !isFALSE(control$dispersion.calibration)
}

L0glm_wide_gaussian_pilot <- function(X, y, family, dispersion, start, prior.mean,
                                     grouped, nonnegative, maxit, lambda, control, included) {
  if (!identical(family$family, "gaussian") || !identical(family$link, "identity") ||
      NCOL(y) < 2L || !isTRUE(grouped) || isTRUE(nonnegative) || !isTRUE(maxit >= 5L) ||
      !is.null(dispersion) || !is.null(start) || !is.null(prior.mean) || isTRUE(lambda == 0) ||
      !L0glm_pilot_enabled(control)) return(FALSE)
  p <- ncol(X) - length(L0glm_mandatory(X, included))
  n <- nrow(X)
  p > n && p <= 10*n
}

# Automatic starts may seed zero entries within a retained row. Explicit user
# starts go directly to the fitter and retain their original eligibility masks.
L0glm_pilot_start_rows <- function(B, nonnegative) {
  magnitude <- sqrt(rowMeans(B^2))
  for (j in which(magnitude > 0 & rowSums(B==0) > 0)) {
    direction <- if (nonnegative) 1 else sign(B[j,which.max(abs(B[j,]))])
    B[j,B[j,]==0] <- .1*magnitude[j]*direction
  }
  B
}

L0glm_gaussian_pilot <- function(initial, lambda, tol, options) {
  data <- initial$data
  state <- initial$state
  state$phi <- .1*state$phi
  state$estimated[] <- FALSE
  options$validation <- NULL
  pilot <- L0glm_ext_pass(data, initial$B, data$allowed, state, lambda, TRUE, 0, 5L, tol, options)
  if (!nrow(pilot$history)) return(initial)
  probes <- matrix(L0glm_native_probes(data$n), data$n, 4L)
  norm <- rowSums(pilot$B^2)
  shape <- vapply(seq_len(data$m), function(j) {
    penalty <- lambda*data$m*state$phi[j]/(4*norm)
    L0glm_gaussian_pilot_scale(data, penalty, probes, j)
  }, numeric(2))
  phi <- shape[1L,]
  reliable <- is.finite(phi) & phi > 0
  phi[!reliable] <- initial$state$phi[!reliable]
  initial$B <- L0glm_pilot_start_rows(pilot$B, FALSE)
  initial$data$allowed <- initial$B != 0
  initial$data$allowed[data$included,] <- TRUE
  initial$data$dispersion.pilot <- list(method="adaptive ridge", applied=TRUE, parameter="dispersion",
    initial=initial$state$phi, calibrated=phi, residualDf=shape[2L,], reliable=reliable)
  initial$data$initial.dispersion <- initial$state$phi <- phi
  initial
}

# Residual trace of the frozen ridge smoother, with the compulsory terms
# projected out. Use the same operator for the response and the four probes.
L0glm_gaussian_pilot_scale <- function(data, penalty, probes, channel) {
  w <- L0glm_observation_column(data$weights, channel)
  good <- which(w > 0)
  n <- length(good)
  unavailable <- c(NA_real_, NA_real_)
  if (n < 3L) return(unavailable)
  penalty[data$included] <- 0
  active <- which(is.finite(penalty) & penalty > 0)
  w <- w[good]
  X <- data$Z[good,,drop=FALSE]
  free <- data$included
  Q <- if (length(free)) {
    qr_free <- qr(as.matrix(Xweights(X[,free,drop=FALSE], w)))
    qr.Q(qr_free)[,seq_len(qr_free$rank),drop=FALSE]
  } else matrix(0, n, 0L)
  project <- function(x) if (ncol(Q)) x-Q%*%crossprod(Q,x) else x
  Z <- project(Xweights(X[,active,drop=FALSE], w))
  residual <- tryCatch({
    if (!length(active)) project else if (length(active) > n) {
      Z <- Z %*% Matrix::Diagonal(x=1/sqrt(penalty[active]))
      factor <- chol(as.matrix(Matrix::tcrossprod(Z) + Matrix::Diagonal(n)))
      function(x) backsolve(factor, forwardsolve(t(factor), as.matrix(project(x))))
    } else {
      factor <- chol(as.matrix(Matrix::crossprod(Z) + Matrix::Diagonal(x=penalty[active])))
      function(x) {
        x <- project(x)
        rhs <- as.matrix(Matrix::crossprod(Z,x))
        x - Z%*%backsolve(factor, forwardsolve(t(factor), rhs))
      }
    }
  }, error=function(e) NULL)
  if (is.null(residual)) return(unavailable)
  e <- as.matrix(residual(probes[good,,drop=FALSE]))
  df <- mean(colSums(e^2))
  if (!is.finite(df) || df <= 0) return(unavailable)
  if (df < min(10,n/5)) {
    df4 <- mean(colSums(as.matrix(residual(e))^2))
    if (!is.finite(df4) || df4 <= 0 || df^2/df4 < min(10,n/5)) return(c(NA_real_,df))
  }
  off <- L0glm_observation_column(data$offset,channel,good)
  rhs <- matrix(sqrt(w)*(data$y[good,channel]-off),ncol=1L)
  c(sum(residual(rhs)^2)/df, df)
}

L0glm_nb_pilot <- function(X, y, family, weights, offset, intercept, always.include,
                           nonnegative, normalize, grouped, openmp, control.family) {
  control <- L0glm_ext_control(control.family)
  y <- as.matrix(y)
  pilot <- tryCatch(suppressWarnings(l0glmFit(X, y, weights=weights, offset=offset,
    family=stats::poisson(family$link), intercept=intercept, always.include=always.include,
    nonnegative=nonnegative, normalize=normalize, grouped=grouped, openmp=openmp,
    maxit=5L, thresh=0, lambda="gic",  parallelize=FALSE, screen="none",
    local.search=FALSE)), error=function(e) NULL)
  if (is.null(pilot)) return(NULL)
  B <- as.matrix(coef(pilot))
  mu <- as.matrix(fitted(pilot))
  m <- ncol(B)
  w <- L0glm_observation_input(weights,nrow(X),m,1,"weights")
  initial <- rep(family$getTheta(TRUE),m)
  theta <- vapply(seq_len(m), function(j) {
    wj <- L0glm_observation_column(w,j)
    k <- sum(B[,j]!=0) + sum(B[pilot$always.include,j]==0)
    if (k >= .7*sum(wj>0)) return(NA_real_)
    L0glm_nb_pearson_shape(y[,j],mu[,j],wj,k,control)
  }, numeric(1))
  reliable <- is.finite(theta) & theta > 0
  if (!any(reliable)) return(NULL)
  theta[!reliable] <- initial[!reliable]
  start <- if (family$link=="log") L0glm_pilot_start_rows(B,nonnegative) else NULL
  list(start=start,theta=theta,diagnostic=list(method="Poisson pilot",applied=TRUE,parameter="theta",
    initial=initial,calibrated=theta,reliable=reliable))
}
