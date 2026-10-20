# Family-specific IRLS extension. Ordinary GLMs never enter these helpers.
L0glm_ext_kind <- function(family) {
  name <- tolower(family$family)
  if (identical(name, "negative binomial") && is.function(family$putTheta)) return("nb")
  if (name %in% c("cox.ph", "cox ph")) return("cox")
  if (name %in% c("multinom", "multinomial")) return("multinomial")
  NULL
}

L0glm_ext_matrix <- function(value, n, m, default, name) {
  if (is.null(value)) return(matrix(default, n, m))
  if (length(value) == 1L) return(matrix(as.numeric(value), n, m))
  if (length(value) == n && (is.null(dim(value)) || NCOL(value) == 1L))
    return(matrix(as.numeric(value), n, m))
  if (length(value) != n * m) stop(name, " has incompatible dimensions.", call. = FALSE)
  matrix(as.numeric(value), n, m)
}

L0glm_ext_control <- function(control) {
  defaults <- list(multinomial = "irls", theta.min = 1e-6, theta.max = 1e7,
                   theta.tol = 1e-6, max.working = 2000L, cox = "auto",
                   cox.llt.max = c(dense = 1024L, sparse = 0L), cox.cg.tol = 1e-10, cox.cg.maxit = 2000L)
  if (!is.list(control) || (length(control) && is.null(names(control))) ||
      any(!names(control) %in% names(defaults)))
    stop("Unknown control.family option.", call. = FALSE)
  defaults[names(control)] <- control
  if (length(defaults$multinomial) != 1L || !defaults$multinomial %in% c("irls", "bound"))
    stop('control.family$multinomial must be "irls" or "bound".', call. = FALSE)
  if (length(defaults$cox) != 1L || !defaults$cox %in% c("auto", "direct", "pcg"))
    stop('control.family$cox must be "auto", "direct" or "pcg".', call. = FALSE)
  cut <- defaults$cox.llt.max
  if (!is.numeric(cut) || !length(cut) %in% c(1L, 2L) || any(!is.finite(cut)) ||
      any(cut < 0 | cut != floor(cut)) ||
      (length(cut) == 2L && !setequal(names(cut), c("dense", "sparse"))))
    stop("cox.llt.max must be a nonnegative integer or named dense/sparse pair.", call. = FALSE)
  numeric <- defaults[setdiff(names(defaults), c("multinomial", "cox", "cox.llt.max"))]
  if (any(lengths(numeric) != 1L) || !all(vapply(numeric, is.numeric, logical(1))) ||
      any(!is.finite(unlist(numeric))) ||
      defaults$theta.min <= 0 || defaults$theta.max <= defaults$theta.min ||
      defaults$theta.tol <= 0 || defaults$max.working < 2 ||
      defaults$cox.cg.tol <= 0 || defaults$cox.cg.tol >= 1 ||
      defaults$cox.cg.maxit < 1 || defaults$cox.cg.maxit != floor(defaults$cox.cg.maxit))
    stop("Invalid family control bounds/tolerances.", call. = FALSE)
  defaults
}

L0glm_ext_prepare <- function(X, y, family, weights, offset, intercept,
                              always.include, nonnegative, normalize, start,
                              prior.mean, control.family) {
  kind <- L0glm_ext_kind(family)
  if (!is.matrix(X) && !inherits(X, "Matrix")) X <- as.matrix(X)
  if (length(dim(X)) != 2L || !nrow(X) ||
      any(!is.finite(if (inherits(X, "sparseMatrix")) X@x else X)))
    stop("X must be a finite numeric design matrix.", call. = FALSE)
  n <- nrow(X)
  response <- switch(kind, nb = L0glm_ext_nb_response(y, family, n),
    cox = L0glm_ext_cox_response(y, n),
    multinomial = L0glm_ext_multinomial_response(y, family, n))
  m <- response$m
  multinomial <- kind == "multinomial"
  if (multinomial && nonnegative) {
    warning("nonnegative is ignored for multinomial CLR coefficients.", call. = FALSE)
    nonnegative <- FALSE
  }
  if (kind == "cox" && intercept) {
    warning("A Cox intercept is not identifiable; intercept was omitted.", call. = FALSE)
    intercept <- FALSE
  }
  if (intercept) X <- cbind(Intercept = 1, X)
  p <- ncol(X)
  included <- L0glm_mandatory(X, always.include)
  lower <- rep(if (nonnegative) 0 else -.Machine$double.xmax, p)
  if (kind == "multinomial") lower[L0glm_intercepts(X)] <- -.Machine$double.xmax
  if (kind == "cox" && length(L0glm_intercepts(X)))
    stop("Remove constant/intercept columns from a Cox design.", call. = FALSE)
  preprocessing <- L0glm_preprocessing(X, normalize, lower, included)
  centered <- L0glm_center_working_design(X, preprocessing, lower)
  scale <- preprocessing$scale
  Z <- if (inherits(centered$X, "sparseMatrix"))
    centered$X %*% Matrix::Diagonal(x = 1 / scale) else
      sweep(as.matrix(centered$X), 2L, scale, "/")
  to_working <- function(B) {
    if (length(centered$intercept)) B[centered$intercept, ] <-
      B[centered$intercept, ] + as.numeric(crossprod(centered$shift, B))
    B * scale
  }
  if (multinomial && !is.null(dim(start)) && ncol(start) != m)
    stop("Multinomial starts need one CLR column per category.", call. = FALSE)
  unit <- L0glm_unit_start(start, X, preprocessing, m)
  explicit <- !is.null(start)
  B <- if (explicit) L0glm_ext_matrix(unit, p, m, 0, "start") else matrix(0, p, m)
  if (any(!is.finite(B))) stop("Starting coefficients must be finite.", call. = FALSE)
  allowed <- if (explicit) B != 0 else matrix(TRUE, p, m)
  if (multinomial && explicit) allowed <- matrix(rowSums(allowed) > 0, p, m)
  allowed[included, ] <- TRUE
  if (multinomial) B <- B - rowMeans(B)
  prior <- to_working(L0glm_ext_matrix(prior.mean, p, m, 0, "prior.mean"))
  if (multinomial) prior <- prior - rowMeans(prior)
  if (any(!is.finite(prior))) stop("prior.mean must be finite.", call. = FALSE)
  B <- to_working(B)
  B <- pmax(B, lower * scale)
  w <- L0glm_ext_matrix(weights, n, if (kind == "nb") m else 1L, 1, "weights")
  off <- L0glm_ext_matrix(offset, n, m, 0, "offset")
  if (multinomial) off <- off - rowMeans(off)
  if (any(!is.finite(w) | w < 0) || !any(w > 0) || any(!is.finite(off)))
    stop("Weights must be finite and nonnegative; offsets must be finite.", call. = FALSE)
  if (kind == "nb" && any(colSums(w) == 0))
    stop("Each NB channel needs a positive-weight observation.", call. = FALSE)
  if (kind == "multinomial" && any(colSums(response$y * w[, 1L]) == 0))
    stop("Each multinomial category needs a positive-weight observation.", call. = FALSE)
  data <- list(kind = kind, X = X, Z = Z, response = response, y = response$y,
    weights = w, offset = off, included = included, lower = lower * scale,
    prior = prior, preprocessing = preprocessing, shift = centered$shift,
    shifted_intercept = centered$intercept, control = L0glm_ext_control(control.family),
    family = family, B = B, allowed = allowed, explicit = explicit, m = m,
    added_intercept = intercept,
    p = p, n = n, p_selectable = p - length(included))
  data$n_ic <- if (kind == "cox") sum(response$status * w[, 1L]) else sum(rowSums(w) > 0)
  if (kind == "cox" && data$n_ic == 0) stop("Cox fitting needs an observed event.", call. = FALSE)
  if (kind == "cox") {
    data$cox.cache <- L0glm_ext_cox_cache(response, w[, 1L])
    data$Z <- data$Z[data$cox.cache$order, , drop = FALSE]
  }
  data
}

L0glm_ext_penalty <- function(B, allowed, included, lambda, grouped, gamma) {
  denominator <- if (grouped && ncol(B) > 1L)
    matrix(rowSums(B^2)^(gamma/2), nrow(B), ncol(B)) else abs(B)^gamma
  D <- if (lambda == 0) matrix(0, nrow(B), ncol(B)) else lambda / (4 * denominator)
  D[!allowed] <- Inf
  D[included, ] <- 0
  D
}

L0glm_ext_threshold <- function(B, allowed, included, thresh, grouped) {
  selectable <- setdiff(seq_len(nrow(B)), included)
  if (grouped && ncol(B) > 1L) {
    magnitude <- sqrt(rowSums(B^2))
    cutoff <- max(c(0, magnitude[selectable])) * thresh
    drop <- magnitude < cutoff | magnitude < 1e-20
    drop[included] <- FALSE
    B[drop, ] <- 0
    allowed[drop, ] <- FALSE
  } else for (j in seq_len(ncol(B))) {
    cutoff <- max(c(0, abs(B[selectable, j]))) * thresh
    drop <- abs(B[, j]) < cutoff | abs(B[, j]) < 1e-20
    drop[included] <- FALSE
    B[drop, j] <- 0
    allowed[drop, j] <- FALSE
  }
  allowed[included, ] <- TRUE
  list(B = B, allowed = allowed)
}

L0glm_ext_ridge <- function(X, y, w, D, start, lower, prior, args) {
  active <- is.finite(D)
  beta <- numeric(length(D))
  if (any(active)) beta[active] <- lmridge(X[, active, drop = FALSE], y,
    weights = w, lambdas = D[active], start = start[active], lower = lower[active],
    prior.mean = prior[active], solver = args$solver, method = args$method,
    control.solver = args$control.solver)
  beta
}

L0glm_ext_pass <- function(data, B, allowed, state, lambda, grouped, gamma,
                           thresh, maxit, tol, args, relaxed = FALSE) {
  converged <- FALSE
  history <- list()
  included <- data$included
  validation <- args$validation
  if (!is.null(validation)) {
    initial <- if (relaxed) B else L0glm_ext_working_coefficients(data, validation$baseline)
    initial_state <- if (relaxed) state else validation$baseline_state
    value <- L0glm_validation_loss(validation, L0glm_ext_original_coefficients(data, initial),
      state = initial_state)
    validation_track <- L0glm_validation_track(NULL, validation, value,
      list(beta = initial, allowed = allowed, state = initial_state), 0L)
  }
  for (iteration in seq_len(maxit)) {
    D <- L0glm_ext_penalty(B, allowed, included, lambda, grouped, gamma)
    proposal <- switch(data$kind,
      nb = L0glm_ext_nb_working(data, B, D, state, args),
      cox = L0glm_ext_cox_working(data, B, D, state, args),
      multinomial = L0glm_ext_multinomial_working(data, B, D, state, args))
    solve <- attr(proposal, "cox.solve")
    attr(proposal, "cox.solve") <- NULL
    ridge <- function(candidate) sum(D[is.finite(D)] * (candidate - data$prior)[is.finite(D)]^2)
    objective <- function(candidate) -2 * L0glm_ext_loglik(data, candidate, state) + ridge(candidate)
    step <- L0glm_step(proposal, B, function(candidate, eta) objective(candidate),
      baseline = if (any((proposal != 0) != (B != 0))) Inf else objective(B),
      maxit = 30L, rtol = 1e-10)
    halving <- step$halvings
    candidate <- step$beta
    if (!step$accepted) {
      warning("No finite descent step in ", data$kind, " IRLS iteration ", iteration, ".", call. = FALSE)
      break
    }
    if (!relaxed && lambda > 0) {
      smooth <- candidate; allowed_smooth <- allowed
      selection <- L0glm_ext_threshold(candidate, allowed, included, thresh, grouped)
      candidate <- selection$B
      allowed <- selection$allowed
      if (data$kind == "multinomial") {
        for (j in seq_len(nrow(candidate))) {
          keep <- which(allowed[j, ])
          if (length(keep) < 2L) { candidate[j, ] <- 0; allowed[j, ] <- FALSE } else
            candidate[j, keep] <- candidate[j, keep] - mean(candidate[j, keep])
        }
      }
      if (!identical(candidate, smooth) &&
          !is.finite(suppressWarnings(L0glm_ext_loglik(data, candidate, state)))) {
        candidate <- smooth; allowed <- allowed_smooth
      }
    }
    nuisance_change <- 0
    if (data$kind == "nb") {
      updated <- L0glm_ext_nb_update(data, candidate, state)
      nuisance_change <- max(abs(log(updated$theta/state$theta)))
      state <- updated
    }
    selectable <- setdiff(seq_len(nrow(B)), included)
    slope_change <- sqrt(sum((candidate[selectable, , drop = FALSE] - B[selectable, , drop = FALSE])^2)) /
      max(1, sqrt(sum(candidate[selectable, , drop = FALSE]^2)))
    prediction_change <- sqrt(sum((data$Z %*% (candidate - B))^2)) /
      max(1, sqrt(sum((data$Z %*% candidate)^2)))
    B <- candidate
    history[[iteration]] <- data.frame(iteration, loglik = L0glm_ext_loglik(data, B, state),
      selected = sum(rowSums(B != 0)[selectable] > 0), step_halvings = halving,
      coefficient_change = max(slope_change, prediction_change), nuisance_change)
    if (data$kind == "cox") {
      history[[iteration]]$backend <- solve$backend
      history[[iteration]]$inner_iterations <- solve$iterations
      history[[iteration]]$inner_residual <- solve$residual
    }
    if (!is.null(validation)) {
      value <- L0glm_validation_loss(validation, L0glm_ext_original_coefficients(data, B), state = state)
      history[[iteration]]$validation_deviance <- value
      validation_track <- L0glm_validation_track(validation_track, validation, value,
        list(beta = B, allowed = allowed, state = state), iteration)
      if (validation_track$stopped) break
    }
    if (max(slope_change, prediction_change) < tol && nuisance_change < data$control$theta.tol) {
      converged <- TRUE
      break
    }
  }
  if (!is.null(validation)) {
    B <- validation_track$best$beta
    allowed <- validation_track$best$allowed
    state <- validation_track$best$state
  }
  result <- list(B = B, allowed = allowed, state = state, converged = converged,
       iterations = iteration, history = if (length(history)) do.call(rbind, history) else data.frame())
  if (!is.null(validation)) result$validation <- L0glm_validation_info(validation_track, iteration)
  result
}

L0glm_ext_fit <- function(X, y, weights, family, dispersion, start, offset, lambda,
                           intercept, always.include, nonnegative, normalize, solver,
                           method, nthreads, maxit, grouped, relaxed, prior.mean,
                           gamma, thresh, tol, control.solver, control.family, validation = NULL) {
  if (!is.numeric(maxit) || length(maxit) != 1L || maxit < 1 || !is.finite(maxit) ||
      !is.numeric(tol) || length(tol) != 1L || tol <= 0 || !is.finite(tol))
    stop("maxit and tol must be finite positive scalars.", call. = FALSE)
  if (!is.null(dispersion) && any(dispersion != 1))
    stop("These families have scale dispersion=1; NB theta is supplied in family=mgcv::nb().", call. = FALSE)
  thread_state <- L0glm_set_thread_count(nthreads)
  if (!is.null(thread_state)) on.exit(L0glm_restore_thread_count(thread_state), add = TRUE)
  data <- L0glm_ext_prepare(X, y, family, weights, offset, intercept, always.include,
                            nonnegative, normalize, start, prior.mean, control.family)
  lambda <- L0glm_resolve_lambda(lambda, max(1, data$p_selectable), max(2, data$n_ic))
  args <- list(solver = solver, method = method, control.solver = control.solver, validation = validation)
  if (data$kind == "cox" && data$control$cox == "pcg" &&
      (solver != "eigen" || any(data$lower > -1e100)))
    stop('Cox cox="pcg" requires solver="eigen" and unrestricted coefficients.', call. = FALSE)
  state <- if (data$kind == "nb") list(theta = data$response$theta, estimated = data$response$estimated) else list()
  B <- data$B
  if (!data$explicit) {
    D <- matrix(1e-3, data$p, data$m)
    D[data$included, ] <- 0
    B <- switch(data$kind,
      nb = L0glm_ext_nb_working(data, B, D, state, args),
      cox = L0glm_ext_cox_initial(data),
      multinomial = L0glm_ext_multinomial_working(data, B, D, state, args))
  }
  if (data$kind == "cox") {
    data$initialization <- attr(B, "cox.solve")
    attr(B, "cox.solve") <- NULL
  }
  if (data$kind == "multinomial" && data$explicit) {
    # Seed eligible CLR zeroes without reopening explicitly excluded rows.
    zero <- if (grouped) matrix(rowSums(B^2) == 0, data$p, data$m) else B == 0
    zero <- zero & data$allowed
    if (any(zero)) {
      D <- matrix(1e-3, data$p, data$m)
      D[!data$allowed] <- Inf; D[data$included, ] <- 0
      pilot <- L0glm_ext_multinomial_working(data, B, D, state, args)
      B[zero] <- pilot[zero]
      B <- B - rowMeans(B)
    }
  }
  selection <- L0glm_ext_pass(data, B, data$allowed, state, lambda, grouped, gamma,
                              thresh, as.integer(maxit), tol, args)
  result <- selection
  passes <- list(selection = selection)
  if (relaxed) {
    active <- if (grouped) matrix(rowSums(selection$B != 0) > 0, data$p, data$m) else selection$B != 0
    active <- active & data$allowed
    active[data$included, ] <- TRUE
    result <- L0glm_ext_pass(data, selection$B, active, selection$state, 0,
                              grouped, gamma, 0, as.integer(maxit), tol, args, TRUE)
    passes$relaxed <- result
  }
  L0glm_ext_result(data, result, passes, lambda, grouped, relaxed, args)
}
