# Validation uses original predictor/response units, independently of the
# training normalization. Only the training observations enter ridge solves.
L0glm_validation_context <- function(validation, X, y, family, intercept) {
  kind <- L0glm_ext_kind(family)
  if (is.null(kind)) kind <- "glm"
  fields <- c("X", "y", "weights", "offset", "patience", "minit", "tol", "starts", "rule", "mu_floor")
  if (!is.list(validation) || is.null(validation$X) || is.null(validation$y))
    stop("validation must contain X and y.", call. = FALSE)
  if (is.null(names(validation)) || any(!names(validation) %in% fields) ||
      anyDuplicated(names(validation)))
    stop("Use named validation entries: ", paste(fields, collapse = ", "), ".", call. = FALSE)
  V <- validation$X
  if (is.data.frame(V)) V <- as.matrix(V)
  if (!(is.numeric(V) || inherits(V, "dMatrix")) || length(dim(V)) != 2L ||
      ncol(V) != ncol(X) || !nrow(V))
    stop("Validation X and y must have compatible finite numeric dimensions.", call. = FALSE)
  n <- nrow(V)
  response <- switch(if (is.null(kind)) "glm" else kind,
    glm = list(y = as.matrix(validation$y), m = NCOL(y)),
    nb = L0glm_ext_nb_response(validation$y, family, n),
    cox = L0glm_ext_cox_response(validation$y, n),
    multinomial = L0glm_ext_multinomial_response(validation$y, family, n,
      L0glm_ext_multinomial_response(y, family, nrow(X))$levels))
  v <- response$y
  m <- if (kind %in% c("cox", "multinomial")) response$m else NCOL(y)
  if (!is.numeric(v) || nrow(v) != n || (kind != "cox" && ncol(v) != m) ||
      any(!is.finite(v)))
    stop("Validation X and y must have compatible finite numeric dimensions.", call. = FALSE)
  entries <- if (inherits(V, "sparseMatrix")) V@x else V
  if (any(!is.finite(entries))) stop("Validation predictors must be finite.", call. = FALSE)
  if (!is.null(colnames(X)) && !is.null(colnames(V))) {
    if (anyDuplicated(colnames(X)) || anyDuplicated(colnames(V)) ||
        !setequal(colnames(X), colnames(V)))
      stop("Validation predictor names must match training predictor names.", call. = FALSE)
    V <- V[, colnames(X), drop = FALSE]
  }
  w <- L0glm_ext_matrix(validation$weights, n, if (kind %in% c("cox", "multinomial")) 1L else m, 1, "Validation weights")
  off <- L0glm_ext_matrix(validation$offset, n, m, 0, "Validation offsets")
  if (any(!is.finite(w) | w < 0) || !is.finite(sum(w)) || sum(w) <= 0 || any(!is.finite(off)))
    stop("Invalid validation weights or offsets.", call. = FALSE)
  if (kind == "glm") {
    name <- family$family
    if (is.function(family$validate_response) && !family$validate_response(v))
      stop("Validation responses lie outside the family domain.", call. = FALSE)
    if ((name %in% c("poisson", "quasipoisson") || startsWith(name, "Negative Binomial(")) && any(v < 0))
      stop("Validation responses must be nonnegative for this family.", call. = FALSE)
    if (name %in% c("Gamma", "inverse.gaussian") && any(v <= 0))
      stop("Validation responses must be positive for this family.", call. = FALSE)
    if (name %in% c("binomial", "quasibinomial") && any(v < 0 | v > 1))
      stop("Validation binomial responses must lie between zero and one.", call. = FALSE)
  }
  patience <- if (is.null(validation$patience)) 3L else validation$patience
  minit <- if (is.null(validation$minit)) 3L else validation$minit
  tol <- if (is.null(validation$tol)) 1e-8 else validation$tol
  for (value in list(patience, minit)) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 1 || value != floor(value))
      stop("Validation patience and minit must be positive integers.", call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0)
    stop("Validation tol must be a finite nonnegative number.", call. = FALSE)
  mu_floor <- if (is.null(validation$mu_floor)) 1e-8 else validation$mu_floor
  if (!is.numeric(mu_floor) || length(mu_floor) != 1L || !is.finite(mu_floor) || mu_floor < 0)
    stop("Validation mu_floor must be a finite nonnegative number.", call. = FALSE)
  mass <- if (identical(kind, "cox")) w[, 1L] * response$status else
    if (identical(kind, "multinomial")) w[, 1L] * rowSums(v) else rowSums(w)
  if (!is.finite(sum(mass)) || sum(mass) <= 0)
    stop("Validation needs positive response mass (observed events for Cox).", call. = FALSE)
  if (intercept) V <- cbind(Intercept = 1, V)
  structure(list(X = V, y = v, weights = w, offset = off, mass = mass,
    family = family, kind = kind, m = m, response = response, mu_floor = mu_floor,
    cox.cache = if (identical(kind, "cox")) L0glm_ext_cox_cache(response, w[, 1L]) else NULL,
    patience = patience, minit = minit, tol = tol), class = "L0glm_validation")
}

L0glm_validation_errors <- function(context, beta, shift = NULL, intercept = integer(), state = NULL) {
  beta <- matrix(as.numeric(beta), ncol(context$X), context$m)
  if (length(intercept)) beta[intercept, ] <- beta[intercept, ] - as.numeric(crossprod(shift, beta))
  active <- which(rowSums(beta != 0) > 0)
  eta <- as.matrix(if (length(active) > nrow(beta) / 2) context$X %*% beta else
    context$X[, active, drop = FALSE] %*% beta[active, , drop = FALSE]) + context$offset
  tryCatch(L0glm_validation_terms(context, eta, state), error = function(e) rep(Inf, nrow(eta)))
}

L0glm_validation_terms <- function(context, eta, state) {
  y <- context$y; w <- context$weights; mass <- context$mass
  if (identical(context$kind, "cox")) {
    z <- L0glm_ext_cox_eta(as.numeric(eta), context$response, w[, 1L], context$cox.cache)
    if (!is.finite(z$loglik)) return(rep(Inf, nrow(y)))
    errors <- numeric(nrow(y))
    event <- context$response$status[z$order] > 0
    sorted <- numeric(length(event))
    sorted[event] <- -2 * (z$eta[event] - z$shift - log(z$S0[z$group[event]]))
    errors[z$order] <- sorted
    return(errors)
  }
  if (identical(context$kind, "multinomial")) {
    lp <- L0glm_ext_softmax(eta)$logp
    saturated <- ifelse(y == 0, 0, y * log(y / pmax(1, rowSums(y))))
    terms <- ifelse(y == 0, 0, y * lp)
    return(2 * rowSums(saturated - terms) / pmax(1, rowSums(y)))
  }
  mu <- matrix(context$family$linkinv(eta), nrow(y), ncol(y))
  if (context$family$family %in% c("poisson", "quasipoisson", "Gamma", "inverse.gaussian") ||
      identical(context$kind, "nb") || startsWith(context$family$family, "Negative Binomial("))
    mu <- pmax(mu, context$mu_floor)
  if (identical(context$kind, "nb")) {
    theta <- if (is.null(state$theta)) context$response$theta else state$theta
    terms <- matrix(0, nrow(y), ncol(y))
    for (j in seq_len(ncol(y))) {
      # A changing theta must not change the saturated reference. The Poisson
      # saturated mass bounds every NB mass (NB is a Poisson mixture).
      saturated <- if (context$response$estimated[j]) stats::dpois(y[, j], y[, j], log = TRUE) else
        stats::dnbinom(y[, j], mu = y[, j], size = context$response$theta[j], log = TRUE)
      terms[, j] <- 2 * (saturated -
        stats::dnbinom(y[, j], mu = mu[, j], size = theta[j], log = TRUE))
    }
  } else terms <- matrix(context$family$dev.resids(as.numeric(y), as.numeric(mu),
                                                 rep(1, length(y))), nrow(y), ncol(y))
  terms[w == 0] <- 0
  rowSums(w * terms) / pmax(.Machine$double.xmin, mass)
}

L0glm_validation_loss <- function(context, beta, shift = NULL, intercept = integer(), state = NULL) {
  errors <- L0glm_validation_errors(context, beta, shift, intercept, state)
  keep <- context$mass > 0
  value <- sum(context$mass[keep] * errors[keep]) / sum(context$mass)
  if (is.finite(value)) value else Inf
}

L0glm_subset_validation <- function(context, columns = NULL, channels = NULL) {
  if (is.null(context)) return(NULL)
  if (!is.null(columns)) {
    context$X <- context$X[, columns, drop = FALSE]
    if (!is.null(context$baseline)) context$baseline <- context$baseline[columns, , drop = FALSE]
    if (!is.null(context$refit_baseline)) context$refit_baseline <- context$refit_baseline[columns, , drop = FALSE]
  }
  if (!is.null(channels)) {
    for (name in c("y", "weights", "offset", "baseline", "refit_baseline"))
      if (!is.null(context[[name]])) context[[name]] <- context[[name]][, channels, drop = FALSE]
    context$m <- length(channels)
    context$mass <- rowSums(context$weights)
    context$baseline_dispersion <- context$baseline_dispersion[channels]
  }
  context
}

# Conditional marginal slopes avoid the global prediction calibration used by
# the default Gaussian initializer. They provide a second, training-only start.
L0glm_marginal_start <- function(X, y, weights, offset, intercept, always.include) {
  if (NCOL(y) > 1L) {
    y <- as.matrix(y); m <- ncol(y)
    w <- L0glm_ext_matrix(weights, nrow(X), m, 1, "weights")
    off <- L0glm_ext_matrix(offset, nrow(X), m, 0, "offset")
    return(do.call(cbind, lapply(seq_len(m), function(j) L0glm_marginal_start(
      X, y[, j], w[, j], off[, j], intercept, always.include))))
  }
  if (intercept) X <- cbind(Intercept = 1, X)
  n <- nrow(X); p <- ncol(X)
  w <- rep_len(if (is.null(weights)) 1 else as.numeric(weights), n)
  z <- as.numeric(y) - rep_len(if (is.null(offset)) 0 else as.numeric(offset), n)
  free <- L0glm_mandatory(X, if (intercept) c(1L, always.include + 1L) else always.include)
  information <- L0glm_moments(X, w)$ss
  response <- z
  if (length(free)) {
    U <- as.matrix(X[, free, drop = FALSE])
    G <- crossprod(U, U * w)
    response <- z - as.numeric(U %*% L0glm_nuisance_solve(G, crossprod(U, w * z)))
    C <- as.matrix(Matrix::crossprod(U, X * w))
    information <- pmax(0, information - colSums(C * L0glm_nuisance_solve(G, C)))
  }
  rhs <- as.numeric(Matrix::crossprod(X, w * response))
  beta <- numeric(p)
  use <- information > 0
  beta[use] <- rhs[use] / information[use]
  beta[free] <- 0
  if (length(free)) beta[free] <- as.numeric(L0glm_nuisance_solve(G,
    crossprod(U, w * (z - as.numeric(X %*% beta)))))
  matrix(beta, ncol = 1L)
}

L0glm_validated_fit <- function(args) {
  X <- args$X; y <- args$y
  family <- if (is.null(args$family)) gaussian() else args$family
  if (is.function(family)) family <- family()
  kind <- L0glm_ext_kind(family)
  intercept <- if (is.null(args$intercept)) identical(kind, "multinomial") || family$family == "binomial" else isTRUE(args$intercept)
  if (identical(kind, "cox")) intercept <- FALSE
  supplied <- args$validation
  context <- L0glm_validation_context(supplied, X, y, family, intercept)
  gaussian <- family$family == "gaussian" && family$link == "identity"
  starts <- if (is.null(supplied$starts)) if (gaussian) c("default", "marginal") else "default" else supplied$starts
  if (!is.character(starts) || !length(starts) || anyNA(starts) ||
      any(!starts %in% c("default", "marginal")))
    stop('Validation starts must be "default", "marginal" or both.', call. = FALSE)
  if (!is.null(args$start)) starts <- "supplied"
  if ("marginal" %in% starts && !gaussian)
    stop('The "marginal" start is a Gaussian identity-link initializer; use "default" or supply start.', call. = FALSE)
  batches <- if (is.null(args$nsis)) NA_real_ else args$nsis
  if (!is.null(args$nsis) && (!is.numeric(batches) || !length(batches) ||
      any(!is.finite(batches) | batches < 1 | batches != floor(batches))))
    stop("Validation screening sizes must be positive integers.", call. = FALSE)
  screen <- L0glm_screen_method(if (is.null(args$screen)) "auto" else args$screen,
    X, args$start, intercept, args$always.include)
  if (length(batches) > 1L && screen == "none")
    stop('A validation batch grid requires screen="SIS" or "ISIS".', call. = FALSE)
  rule <- if (is.null(supplied$rule)) "min" else supplied$rule
  if (!is.character(rule) || length(rule) != 1L || is.na(rule) || !rule %in% c("min", "1se"))
    stop('Validation rule must be "min" or "1se".', call. = FALSE)
  models <- diagnostics <- list()
  add <- function(fit_args, key, start) {
    elapsed <- system.time(fit <- do.call(L0glm.fit, fit_args), gcFirst = FALSE)[["elapsed"]]
    diagnostics[[key]] <<- data.frame(start = start,
      nsis = if (is.null(fit$screening)) NA_real_ else fit$screening$screen_size,
      validation_deviance = L0glm_validation_loss(context, coef(fit), state = list(theta = fit$theta)),
      seconds = elapsed, selected = sum(coef(fit) != 0))
    design <- fit$X
    fit$X <- NULL
    models[[key]] <<- fit
    invisible(design)
  }
  null_args <- args
  null_args$validation <- NULL; null_args$start <- 0
  null_args$screen <- "none"; null_args$nsis <- NULL; null_args$relaxed <- FALSE
  null_args$intercept <- intercept
  design <- add(null_args, "null", "null")
  null <- models$null
  context$baseline <- as.matrix(coef(null))
  context$baseline_dispersion <- null$dispersion
  context$baseline_state <- list(theta = null$theta, estimated = null$theta.estimated)
  value <- diagnostics$null$validation_deviance
  models$null$validation <- list(selection = list(best_iteration = 0L, iterations = 0L,
    stopped = FALSE, deviance = value, history = data.frame(iteration = 0L, deviance = value,
      selected = sum(coef(null) != 0), dispersion = NA_real_)), deviance = value)
  marginal <- if ("marginal" %in% starts) L0glm_marginal_start(X, y, args$weights,
    args$offset, intercept, args$always.include) else NULL
  for (batch in unique(batches)) for (name in unique(starts)) {
    fit_args <- args; fit_args$validation <- context
    if (!is.na(batch)) fit_args$nsis <- batch
    if (name == "marginal") fit_args$start <- marginal
    add(fit_args, paste(name, batch, sep = "_"), name)
  }
  table <- do.call(rbind, diagnostics)
  if (!any(is.finite(table$validation_deviance)))
    stop("No candidate has finite validation deviance.", call. = FALSE)
  best <- which.min(table$validation_deviance)
  keep <- context$mass > 0
  errors <- L0glm_validation_errors(context, coef(models[[best]]), state = list(theta = models[[best]]$theta))[keep]
  w <- context$mass[keep] / sum(context$mass); sum_w2 <- sum(w^2)
  se <- if (sum_w2 < 1) sqrt(sum(w^2 * (errors - table$validation_deviance[best])^2) /
                              (1 - sum_w2)) else 0
  eligible <- which(table$validation_deviance <= table$validation_deviance[best] +
    if (rule == "1se") se else 0)
  chosen <- eligible[order(table$selected[eligible], table$nsis[eligible], na.last = TRUE)][1L]
  result <- models[[chosen]]
  result$X <- design
  result$validation$starts <- table
  result$validation$selected_start <- table$start[chosen]
  result$validation$selected_nsis <- table$nsis[chosen]
  result$validation$rule <- rule
  result$validation$scoring <- switch(if (is.null(kind)) "glm" else kind,
    glm = "mean unit deviance", nb = "mean predictive deviance with fixed saturated reference",
    multinomial = "mean deviance per trial", cox = "mean partial deviance per event")
  result$validation$minimum <- table$validation_deviance[best]
  result$validation$se <- se
  result
}

# One stopping rule for ordinary GLMs and joint special-family passes.
L0glm_validation_track <- function(track, context, value, checkpoint, iteration) {
  row <- data.frame(iteration = iteration, deviance = value,
    selected = sum(checkpoint$beta != 0),
    dispersion = if (length(checkpoint$dispersion) == 1L) checkpoint$dispersion else NA_real_)
  if (is.null(track)) return(list(best = checkpoint, value = value,
    best_iteration = iteration, stalls = 0L, stopped = FALSE, history = list(row)))
  track$history[[length(track$history) + 1L]] <- row
  if (is.finite(value) && (!is.finite(track$value) ||
      value < track$value - context$tol * max(1, abs(track$value)))) {
    track$best <- checkpoint; track$value <- value
    track$best_iteration <- iteration; track$stalls <- 0L
  } else track$stalls <- track$stalls + 1L
  track$stopped <- iteration >= context$minit && track$stalls >= context$patience
  track
}

L0glm_validation_info <- function(track, iterations) list(best_iteration = track$best_iteration,
  iterations = iterations, stopped = track$stopped, deviance = track$value,
  history = do.call(rbind, track$history))
