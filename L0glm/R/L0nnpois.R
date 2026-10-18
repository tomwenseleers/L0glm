#' Fast nonnegative identity-link Poisson L0 or grouped L0 penalized fit
#'
#' `L0nnpois()` is a fast RcppEigen/OpenMP implementation of
#' [L0glm.fit()] for the nonnegative Poisson model with identity link. It shares
#' the main fitting arguments with `L0glm.fit()` for use in deconvolution
#' workflows, and internally uses the sparse or dense Eigen Cholesky LLT ridge solver.
#' `L0glm.fit(openmp = TRUE)` selects this backend automatically for eligible
#' fits. Compared with the general R fitting path (`openmp = FALSE`), it keeps
#' the design matrix `X` in shared memory across outcome channels.
#' Dense inputs retain dense storage and sparse inputs retain sparse storage.
#' The dense solver forms weighted crossproducts using R's BLAS library.
#' Eligible wide active systems use the Woodbury identity to factor an
#' observation-space system rather than a predictor-space system.
#' The default penalty is `lambda = "gic"`.
#' Every coefficient, including any intercept, is constrained to be nonnegative.
#' This constraint cannot be disabled. Use `L0glm.fit(nonnegative = FALSE)`
#' if coefficients of either sign are required.
#'
#' @inheritParams L0glm.fit
#' @param dispersion Poisson dispersion is fixed at 1. Accepts \code{NULL}
#'   or 1; other finite nonzero numeric values are ignored with a warning.
#'   Negative values cannot request dispersion estimation for this family.
#' @param intercept Logical; if \code{TRUE}, add an intercept column to \code{X}.
#'   The intercept is constrained to be nonnegative, just like every other
#'   coefficient. Defaults to \code{FALSE}.
#' @param ncores Optional integer specifying the number of OpenMP threads to use
#'   across outcome channels when \code{parallelize=TRUE}. For \code{L0nnpois},
#'   the requested value is capped only at the detected number of logical cores,
#'   not at the number of physical cores.
#' @param control.solver Optional list of solver controls. \code{L0nnpois}
#'   always uses Eigen Cholesky LLT internally, but for sparse matrices
#'   \code{control.solver$NaturalOrdering} can be set to \code{TRUE} or
#'   \code{FALSE} to choose Eigen \code{NaturalOrdering} or
#'   \code{AMDOrdering}. The default is \code{FALSE}. The legacy integer
#'   \code{control.solver$ordering} is also accepted, with \code{1L} meaning
#'   natural ordering and \code{0L} meaning AMD ordering.
#'
#' @return An object of class \code{L0glm_fit}, with the same main components as [L0glm.fit()]:
#'   `coefficients`, `weights`, `dispersion`, `init.lambda`, `lambdas`,
#'   `family`, `z`, `X`, and `y`. Original prior weights and offsets are
#'   retained separately from the working weights for scoring. Methods
#'   \code{coef()}, \code{logLik()}, \code{AIC()}, \code{BIC()}, \code{GIC()},
#'   \code{EBIC()} and \code{MBIC()} are available; see \link{model_criteria}.
#'   \code{print()} and \code{summary()} provide compact model displays;
#'   see \link{model_summary}.
#'
#' @seealso [L0glm.fit()], \link{information_criteria} for case-insensitive
#'   criterion presets, penalty formulas, guidance and original references.
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
                     lambda        = "gic",
                     intercept     = FALSE,
                     always.include        = NULL,
                     normalize     = "auto",
                     solver        = "eigen",
                     method        = 1L,
                     parallelize   = TRUE,
                     ncores        = L0glm_default_ncores(parallelize),
                     nthreads      = NULL,
                     maxit         = 25,
                     grouped       = TRUE,
                     group.maxit   = 6,
                     relaxed       = grouped & NCOL(y) > 1,
                     prior.mean    = NULL,
                     gamma         = 2,
                     thresh         = 1e-5,
                     tol           = 1e-10,
                     control.solver = do.call("control.solver.gen", list(solver = solver, method = method)),
                     verbose       = FALSE) {
  validate_thresh(thresh)
  if (!solver %in% c("eigen", "bigGlm", "bigglm", "glmnet", "osqp"))
    stop("Invalid solver choice; use eigen, bigGlm, or osqp.", call. = FALSE)
  normalization <- L0glm_normalization(normalize, nonnegative = TRUE)
  normalization_code <- match(normalization, c("none", "l2", "rms", "sd")) - 1L
  if (!inherits(family, "family")) family <- family()
  if (!identical(family$family, "poisson") || !identical(family$link, "identity")) {
    stop("L0nnpois is only defined for family = poisson(identity).",
         call. = FALSE)
  }
  if (!identical(solver, "eigen") || as.integer(method[1L]) != 1L) {
    warning("L0nnpois always uses the Eigen Cholesky LLT solver; ",
            "solver and method were ignored.", call. = FALSE)
  }
  dispersion <- L0glm_dispersion(dispersion, family, NCOL(y))
  if (!is.null(nthreads) && verbose) {
    message("L0nnpois uses ncores for OpenMP outcome-level threading; nthreads is ignored.")
  }
  natural_ordering <- FALSE
  if (!is.null(control.solver)) {
    if (!is.null(control.solver$NaturalOrdering)) {
      natural_ordering <- isTRUE(control.solver$NaturalOrdering)
    } else if (!is.null(control.solver$ordering)) {
      natural_ordering <- as.integer(control.solver$ordering[1L]) == 1L
    }
  }

  y <- as.matrix(y)
  storage.mode(y) <- "double"
  if (any(!is.finite(y)) || any(y < 0)) {
    stop("Poisson responses must be finite and nonnegative.", call. = FALSE)
  }
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

  if (intercept) {
    predictor_names <- colnames(X)
    X <- if (issparse) {
      methods::as(cbind(Intercept = Matrix::Matrix(1, n, 1, sparse = TRUE), X),
                  "dgCMatrix")
    } else {
      cbind(Intercept = 1, X)
    }
    colnames(X) <- c("Intercept", if (is.null(predictor_names)) rep("", p0) else predictor_names)
    always.include <- c(1L, always.include + 1L)
  }
  p <- ncol(X)
  always.include <- L0glm_mandatory(X, always.include)
  p_selectable <- p - length(always.include)
  preprocessing <- L0glm_preprocessing(X, normalize, rep(0, p), always.include)
  lambda <- L0glm_resolve_lambda(lambda, p = max(1L, p_selectable), n = n)
  init.lambda <- lambda
  start <- L0glm_unit_start(start, X, preprocessing, m)

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
    if (length(start) == 1L) {
      start <- matrix(as.numeric(start), p, m)
    } else {
      start <- as.matrix(start)
      if (intercept && nrow(start) == p0) {
        start <- rbind(0, start)
      }
      if (m > 1L && ncol(start) == 1L) {
        start <- matrix(rep(start, m), nrow = nrow(start), ncol = m)
      }
    }
    storage.mode(start) <- "double"
    if (!identical(dim(start), c(p, m))) {
      stop("start must be a scalar, have dimension p x m, or have length p for univariate y.",
           call. = FALSE)
    }
    if (any(!is.finite(start)))
      stop("start must contain finite coefficients.", call. = FALSE)
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
  allowed <- if (is.null(start)) NULL else start != 0
  if (!is.null(allowed)) allowed[always.include, ] <- TRUE

  if (isTRUE(parallelize) && m > 1L) {
    ncores <- suppressWarnings(as.integer(ncores[1L]))
    if (is.na(ncores) || ncores < 1L) ncores <- 1L
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
      X, y, weights, offset, prior.mean, start, always.include,
      lambda, normalization_code, preprocessing$scale, grouped, relaxed, as.integer(group.maxit),
      as.integer(maxit), gamma, thresh, tol, ncores,
      natural_ordering)
  } else {
    fit <- L0nnpois_dense_cpp(
      X, y, weights, offset, prior.mean, start, always.include,
      lambda, normalization_code, preprocessing$scale, grouped, relaxed, as.integer(group.maxit),
      as.integer(maxit), gamma, thresh, tol, ncores,
      natural_ordering)
  }

  coefficients <- fit$coefficients
  weights.ret <- fit$weights
  lambdas <- fit$lambdas
  z <- fit$z
  dispersion <- fit$dispersion
  if (!is.null(allowed)) lambdas[!allowed] <- Inf

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
  coefficients <- methods::as(methods::as(coefficients, "generalMatrix"), "CsparseMatrix")

  structure(list(coefficients = coefficients,
       weights      = weights.ret,
       prior.weights = weights,
       offset       = offset,
       dispersion   = dispersion,
       dispersion.estimated = rep(FALSE, m),
       init.lambda  = init.lambda,
       lambdas      = lambdas,
       family       = family,
       p            = p_selectable,
       always.include = always.include,
       preprocessing = preprocessing,
       convergence = fit$convergence,
       added_intercept = intercept,
       grouped      = grouped,
       z            = z,
       X            = if (issparse) X else Matrix::Matrix(X, sparse = FALSE),
       y            = y), class = "L0glm_fit")
}
