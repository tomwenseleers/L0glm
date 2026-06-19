#' Fast nonnegative identity-link Poisson L0 or grouped L0 penalized fit
#'
#' `L0nnpois()` is a fast RcppEigen/OpenMP implementation of
#' [L0glm.fit()] for the nonnegative Poisson model with identity link. It keeps
#' the same argument surface as `L0glm.fit()` for drop-in use in e.g. deconvolution
#' workflows, and internally uses the sparse or dense Eigen Cholesky LLT ridge solver.
#' For problems with many outcomes, it is much faster than `L0glm.fit()`,
#' and also much more memory efficient, as the design matrix  `X` is kept in
#' shared memory.
#'
#' @inheritParams L0glm.fit
#' @param ncores Optional integer specifying the number of OpenMP threads to use
#'   across outcome channels when \code{parallelize=TRUE}. For \code{L0nnpois},
#'   the requested value is capped only at the detected number of logical cores,
#'   not at the number of physical cores.
#'
#' @return A list with the same main components as [L0glm.fit()]:
#'   `coefficients`, `weights`, `dispersion`, `init.lambda`, `lambdas`,
#'   `family`, `z`, `X`, and `y`.
#'
#' @seealso [L0glm.fit()]
#'
#' @example examples/L0nnpois.R
#'
#' @export
L0nnpois <- function(X,
                     y,
                     weights       = NULL,
                     family        = poisson(identity),
                     dispersion    = NULL,
                     start         = NULL,
                     offset        = NULL,
                     lambda        = 1,
                     intercept     = FALSE,
                     always.include        = NULL,
                     nonnegative   = TRUE,
                     normalize     = FALSE,
                     solver        = "eigen",
                     method        = 7L,
                     lower         = NULL,
                     upper         = NULL,
                     parallelize   = TRUE,
                     ncores        = L0glm_default_ncores(parallelize),
                     nthreads      = NULL,
                     maxit         = 25,
                     grouped       = TRUE,
                     group.maxit   = 6,
                     relaxed       = grouped & NCOL(y) > 1,
                     warm.start    = 0.9,
                     prior.mean    = NULL,
                     gamma         = 2,
                     thresh        = 1e-5,
                     tol           = 1e-10,
                     control.solver = do.call("control.solver.gen", list(solver = solver, method = method)),
                     verbose       = FALSE) {
  if (!inherits(family, "family")) family <- family()
  if (!identical(family$family, "poisson") || !identical(family$link, "identity")) {
    stop("L0nnpois is only defined for family = poisson(identity).",
         call. = FALSE)
  }
  if (!identical(solver, "eigen") || as.integer(method[1L]) != 7L) {
    warning("L0nnpois always uses the Eigen Cholesky LLT solver; ",
            "solver and method were ignored.", call. = FALSE)
  }
  if (!is.null(dispersion) && any(as.numeric(dispersion) != 1)) {
    warning("Poisson dispersion is fixed at 1 in L0nnpois; dispersion was ignored.",
            call. = FALSE)
  }
  if (!is.null(nthreads) && verbose) {
    message("L0nnpois uses ncores for OpenMP outcome-level threading; nthreads is ignored.")
  }

  y <- as.matrix(y)
  storage.mode(y) <- "double"
  n <- nrow(X)
  p0 <- ncol(X)
  m <- ncol(y)
  if (!identical(nrow(y), n)) {
    stop("nrow(y) must equal nrow(X).", call. = FALSE)
  }

  issparse <- inherits(X, "sparseMatrix")
  if (issparse) {
    X <- methods::as(Matrix::Matrix(X, sparse = TRUE), "dgCMatrix")
  } else {
    X <- as.matrix(X)
    storage.mode(X) <- "double"
  }

  make_bound <- function(x, default, name) {
    if (is.null(x)) return(rep(default, p0))
    x <- as.numeric(x)
    if (length(x) == 1L) return(rep(x, p0))
    if (length(x) != p0) {
      stop(name, " must have length 1 or ncol(X).", call. = FALSE)
    }
    x
  }
  lower <- make_bound(lower,
                      if (isTRUE(nonnegative)) 0 else -.Machine$double.xmax,
                      "lower")
  upper <- make_bound(upper, .Machine$double.xmax, "upper")
  if (isTRUE(nonnegative)) lower[lower < 0] <- 0
  if (any(lower < 0)) {
    stop("L0nnpois requires nonnegative lower bounds; use nonnegative = TRUE ",
         "or lower = 0.", call. = FALSE)
  }

  if (intercept) {
    X <- if (issparse) {
      methods::as(cbind(Intercept = Matrix::Matrix(1, n, 1, sparse = TRUE), X),
                  "dgCMatrix")
    } else {
      cbind(Intercept = 1, X)
    }
    always.include <- c(1L, always.include + 1L)
    lower <- c(-.Machine$double.xmax, lower)
    upper <- c(.Machine$double.xmax, upper)
  }
  p <- ncol(X)

  if (is.character(lambda)) {
    lambda <- match.arg(tolower(lambda),
                        choices = c("aic", "bic", "gic", "ebic", "mbic"))
    lambda <- preset.lambda(IC = lambda, p = p, n = n)
  }
  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda)) {
    stop("lambda must be a finite scalar or one of 'aic', 'bic', 'gic', ",
         "'ebic', or 'mbic'.", call. = FALSE)
  }
  init.lambda <- lambda

  if (is.null(weights)) {
    weights <- matrix(1, n, m)
  } else {
    weights <- as.matrix(weights)
    storage.mode(weights) <- "double"
    if (nrow(weights) != n) {
      stop("nrow(weights) must equal nrow(X).", call. = FALSE)
    }
    if (m > 1L && ncol(weights) == 1L) {
      weights <- matrix(rep(weights, m), n, m)
    }
    if (!identical(dim(weights), c(n, m))) {
      stop("weights must have dimension n x m, or length n.", call. = FALSE)
    }
  }

  if (is.null(offset)) {
    offset <- matrix(0, n, m)
  } else {
    offset <- as.matrix(offset)
    storage.mode(offset) <- "double"
    if (nrow(offset) != n) {
      stop("nrow(offset) must equal nrow(X).", call. = FALSE)
    }
    if (m > 1L && ncol(offset) == 1L) {
      offset <- matrix(rep(offset, m), n, m)
    }
    if (!identical(dim(offset), c(n, m))) {
      stop("offset must have dimension n x m, or length n.", call. = FALSE)
    }
  }

  if (!is.null(start)) {
    start <- as.matrix(start)
    storage.mode(start) <- "double"
    if (intercept && nrow(start) == p0) {
      start <- rbind(0, start)
    }
    if (m > 1L && ncol(start) == 1L) {
      start <- matrix(rep(start, m), nrow = nrow(start), ncol = m)
    }
    if (!identical(dim(start), c(p, m))) {
      stop("start must have dimension p x m (or length p for univariate y).",
           call. = FALSE)
    }
    start <- warm.start * start +
      (1 - warm.start) * rep(stats::sd(start)^(1 / gamma), p)
  }

  if (is.null(prior.mean)) {
    prior.mean <- matrix(0, p, m)
  } else {
    prior.mean <- as.matrix(prior.mean)
    storage.mode(prior.mean) <- "double"
    if (intercept && nrow(prior.mean) == p0) {
      prior.mean <- rbind(0, prior.mean)
    }
    if (m == 1L) {
      if (length(prior.mean) != p) {
        stop("prior.mean must have length p for univariate y.", call. = FALSE)
      }
      prior.mean <- matrix(prior.mean, p, 1L)
    } else if (ncol(prior.mean) == 1L) {
      prior.mean <- matrix(rep(prior.mean, m), p, m)
    } else if (!identical(dim(prior.mean), c(p, m))) {
      stop("prior.mean must be p x m when y has m columns.", call. = FALSE)
    }
  }

  always.include <- as.integer(always.include)
  always.include <- always.include[!is.na(always.include) &
                                     always.include >= 1L &
                                     always.include <= p]

  ncores <- suppressWarnings(as.integer(ncores[1L]))
  if (is.na(ncores) || ncores < 1L) ncores <- 1L
  if (isTRUE(parallelize) && m > 1L) {
    logical_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
    if (!is.na(logical_cores) && logical_cores >= 1L &&
        ncores > logical_cores) {
      warning("ncores capped at logical core count = ", logical_cores,
              call. = FALSE)
      ncores <- as.integer(logical_cores)
    }
  } else {
    ncores <- 1L
  }

  if (issparse) {
    fit <- L0nnpois_sparse_cpp(
      X, y, weights, offset, prior.mean, start, always.include, lower, upper,
      lambda, normalize, grouped, relaxed, as.integer(group.maxit),
      as.integer(maxit), warm.start, gamma, thresh, tol, ncores)
  } else {
    fit <- L0nnpois_dense_cpp(
      X, y, weights, offset, prior.mean, start, always.include, lower, upper,
      lambda, normalize, grouped, relaxed, as.integer(group.maxit),
      as.integer(maxit), warm.start, gamma, thresh, tol, ncores)
  }

  coefficients <- fit$coefficients
  weights.ret <- fit$weights
  lambdas <- fit$lambdas
  z <- fit$z
  dispersion <- fit$dispersion

  weights.ret[is.infinite(weights.ret)] <- 1e11
  lambdas[is.infinite(lambdas)] <- 1e11

  if (m == 1L) {
    coefficients <- Matrix::Matrix(coefficients[, 1L], ncol = 1L,
                                   sparse = TRUE,
                                   dimnames = list(colnames(X), NULL))
    weights.ret <- as.numeric(weights.ret[, 1L])
    lambdas <- as.numeric(lambdas[, 1L])
    z <- as.numeric(z[, 1L])
    dispersion <- as.numeric(dispersion[1L])
  } else {
    coefficients <- Matrix::Matrix(coefficients, sparse = TRUE,
                                   dimnames = list(colnames(X), colnames(y)))
    weights.ret <- Matrix::Matrix(weights.ret, sparse = TRUE)
    lambdas <- Matrix::Matrix(lambdas, sparse = TRUE)
    z <- Matrix::Matrix(z, sparse = TRUE)
    dispersion <- as.numeric(dispersion)
  }

  list(coefficients = coefficients,
       weights      = weights.ret,
       dispersion   = dispersion,
       init.lambda  = init.lambda,
       lambdas      = lambdas,
       family       = family,
       z            = z,
       X            = if (issparse) X else Matrix::Matrix(X, sparse = FALSE),
       y            = y)
}
