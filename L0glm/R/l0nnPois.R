# Specialized nonnegative identity-link Poisson backend.
l0nnPois <- function(X,
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
                     ncores        = L0glm_default_ncores(),
                     nthreads      = NULL,
                     maxit         = 25,
                     grouped       = TRUE,
                     group.maxit   = 6,
                     relaxed       = FALSE,
                     prior.mean    = NULL,
                     thresh         = 1e-5,
                     tol           = 1e-10,
                     control.solver = do.call("controlSolverGen", list(solver = solver, method = method)),
                     verbose       = FALSE,
                     validation    = NULL,
                     precluster    = FALSE,
                     cluster.refit = TRUE,
                     control.precluster = list(),
                     local.search = "auto",
                     control.local = list(), group = NULL, group.multiplier = NULL,
                     interactions = 0L, control.interactions = list()) {
  interactions <- L0glm_interaction_order(interactions)
  if(inherits(X,'L0glm_data'))X <- L0glm_data_matrix(X)
  if(interactions>0L) {
    if(!inherits(family,'family'))family <- family()
    if(family$family!='poisson' || family$link!='identity')
      stop('l0nnPois is only defined for family = poisson(identity).',call.=FALSE)
    supplied <- names(match.call())[-1L]
    args <- mget(supplied,environment(),inherits=FALSE)
    args$family <- family;args$nonnegative <- TRUE;args$screen <- 'none'
    if(missing(lambda))args$lambda <- 'gic'
    return(L0glm_interaction_fit(args,interactions,control.interactions))
  }
  if (!is.null(group)) {
    if(!inherits(family,'family'))family <- family()
    if(family$family!='poisson' || family$link!='identity')
      stop('l0nnPois is only defined for family = poisson(identity).',call.=FALSE)
    args <- mget(names(formals(l0nnPois)),environment(),inherits=FALSE)
    args$relaxed <- relaxed
    return(do.call(L0glm_internal_fit,c(args,list(nonnegative=TRUE,screen='none'))))
  }
  if (!identical(local.search, FALSE)) {
    mode <- L0glm_local_mode(local.search)
    supplied <- setdiff(names(match.call())[-1L], c("local.search", "control.local"))
    if (NCOL(y)==1L) supplied <- setdiff(supplied,"ncores")
    arguments <- mget(supplied, envir=environment(), inherits=FALSE)
    arguments$relaxed <- relaxed
    return(L0glm_local_fit(L0glm_internal_pois, arguments, mode, control.local, nonnegative=TRUE))
  }
  L0glm_precluster_flags(precluster,cluster.refit)
  if (!solver %in% c("eigen", "bigGlm", "bigglm", "glmnet", "osqp"))
    stop("Invalid solver choice; use eigen, bigGlm, or osqp.", call. = FALSE)
  nuisance <- isTRUE(maxit>1L) && !isFALSE(control.solver$.nuisance)
  if(nuisance) {
    if(is.data.frame(X)) X <- as.matrix(X)
    if(!inherits(family,"family")) family <- family()
    nuisance <- isTRUE(intercept) || length(always.include)>0L || length(L0glm_intercepts(X))>0L
    if(!nuisance) control.solver$.nuisance <- FALSE
  }
  if (nuisance && (!precluster || cluster.refit)) {
    if(!inherits(family,"family")) family <- family()
    if(!identical(L0glm_ext_kind(family),"multinomial") &&
       !(family$family %in% c("binomial","quasibinomial") && !length(always.include))) {
      supplied <- names(match.call())[-1L]
      if(!is.null(validation) && !inherits(validation,"L0glm_validation")) {
        args <- mget(supplied,envir=environment(),inherits=FALSE)
        args$family <- family
        args$nonnegative <- TRUE; args$screen <- "none"
        return(L0glm_validated_fit(args))
      }
      options <- mget(setdiff(names(formals(L0glm_nnpois_kernel)),c("X","y")),envir=environment(),inherits=FALSE)
      options$nonnegative <- TRUE
      return(L0glm_nuisance_fit(X,y,options,L0glm_internal_pois,supplied))
    }
  }
  if(precluster) {
    if(!inherits(family,"family")) family<-family()
    options<-mget(setdiff(names(formals(L0glm_nnpois_kernel)),c("X","y")),envir=environment(),inherits=FALSE)
    options$nonnegative<-TRUE
    return(L0glm_precluster_fit(X,y,options,L0glm_internal_pois,names(match.call())[-1L]))
  }
  validate_thresh(thresh)
  normalization <- L0glm_normalization(normalize, nonnegative = TRUE)
  normalization_code <- match(normalization, c("none", "l2", "rms", "sd")) - 1L
  if (!inherits(family, "family")) family <- family()
  if (!identical(family$family, "poisson") || !identical(family$link, "identity")) {
    stop("l0nnPois is only defined for family = poisson(identity).",
         call. = FALSE)
  }
  if (!identical(solver, "eigen") || as.integer(method[1L]) != 1L) {
    warning("l0nnPois always uses the Eigen Cholesky LLT solver; ",
            "solver and method were ignored.", call. = FALSE)
  }
  dispersion <- L0glm_dispersion(dispersion, family, NCOL(y))
  if (!is.null(validation) && !inherits(validation, "L0glm_validation")) {
    args <- mget(names(match.call())[-1L], envir = environment(), inherits = FALSE)
    args$family <- family; args$lambda <- lambda; args$intercept <- intercept
    args$nonnegative <- TRUE
    args$solver <- "eigen"; args$method <- 1L
    # Direct calls keep this specialized backend; l0glmFit owns screening.
    args$screen <- "none"
    args$control.solver <- control.solver[intersect(names(control.solver),
      c("NaturalOrdering", "ordering", "nonnegative.polish"))]
    return(do.call(L0glm_internal_fit, args))
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
  p_selectable <- L0glm_criterion_size(p-length(always.include),control.solver)
  preprocessing <- L0glm_preprocessing(X, normalize, rep(0, p), always.include)
  lambda <- L0glm_resolve_lambda(lambda, p = max(1L, p_selectable), n = n)
  init.lambda <- lambda
  start <- L0glm_unit_start(start, X, preprocessing, m)

  weights <- L0glm_observation_input(weights,n,m,1,"weights")
  offset <- L0glm_observation_input(offset,n,m,0,"offset")
  weights <- as.matrix(weights);storage.mode(weights)<-"double"
  offset <- as.matrix(offset);storage.mode(offset)<-"double"

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
  if (m > 1L) {
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
  thread_state <- L0glm_set_thread_count(if(ncores>1L)1L else nthreads)
  if(!is.null(thread_state))on.exit(L0glm_restore_thread_count(thread_state),add=TRUE)
  if (is.null(start) && maxit > 1L && any(offset != 0))
    start <- L0glm_signal_start(X,y,offset,weights,always.include,prior.mean,ncores)
  allowed <- if (is.null(start)) NULL else start != 0
  if (!is.null(allowed)) allowed[always.include, ] <- TRUE

  if (!is.null(validation)) {
    validation$X <- if (inherits(validation$X, "sparseMatrix"))
      methods::as(validation$X, "dgCMatrix") else as.matrix(validation$X)
  }
  if (issparse) {
    fit <- L0nnpois_sparse_cpp(
      X, y, weights, offset, prior.mean, start, always.include,
      lambda, normalization_code, preprocessing$scale, grouped, relaxed, as.integer(group.maxit),
      as.integer(maxit), thresh, tol, ncores,
      natural_ordering, validation,isTRUE(control.solver$.refit.only),
      isTRUE(control.solver$nonnegative.polish))
  } else {
    fit <- L0nnpois_dense_cpp(
      X, y, weights, offset, prior.mean, start, always.include,
      lambda, normalization_code, preprocessing$scale, grouped, relaxed, as.integer(group.maxit),
      as.integer(maxit), thresh, tol, ncores,
      natural_ordering, validation,isTRUE(control.solver$.refit.only),
      isTRUE(control.solver$nonnegative.polish))
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

  result <- structure(list(coefficients = coefficients,
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
  if (!is.null(validation)) {
    result$validation <- fit$validation
    result$validation$deviance <- L0glm_validation_loss(validation, coefficients)
  }
  result$relaxationSkipped <- fit$relaxation_skipped
  result
}

L0glm_nnpois_kernel <- l0nnPois
L0glm_internal_pois <- L0glm_recorded_fitter(L0glm_nnpois_kernel,'l0nnPois',public=FALSE)
l0nnPois <- L0glm_internal_pois
