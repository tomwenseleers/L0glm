# Validation uses original predictor/response units, independently of the
# training normalization. Only the training observations enter ridge solves.
L0glm_validation_context <- function(validation, X, y, family, intercept) {
  if (family$family != "gaussian" || family$link != "identity" || NCOL(y) != 1L)
    stop("Validation currently supports one Gaussian identity-link response.", call. = FALSE)
  fields <- c("X", "y", "weights", "offset", "patience", "minit", "tol", "starts", "rule")
  if (!is.list(validation) || is.null(validation$X) || is.null(validation$y))
    stop("validation must contain X and y.", call. = FALSE)
  if (is.null(names(validation)) || any(!names(validation) %in% fields) ||
      anyDuplicated(names(validation)))
    stop("Use named validation entries: ", paste(fields, collapse = ", "), ".", call. = FALSE)
  V <- validation$X
  if (is.data.frame(V)) V <- as.matrix(V)
  v <- validation$y
  if (!(is.numeric(v) || inherits(v, "dMatrix")) || NCOL(v) != 1L ||
      !(is.numeric(V) || inherits(V, "dMatrix")) || length(dim(V)) != 2L ||
      ncol(V) != ncol(X) || nrow(V) != length(v) || !length(v) || any(!is.finite(v)))
    stop("Validation X and y must have compatible finite numeric dimensions.", call. = FALSE)
  v <- as.numeric(v)
  entries <- if (inherits(V, "sparseMatrix")) V@x else V
  if (any(!is.finite(entries))) stop("Validation predictors must be finite.", call. = FALSE)
  if (!is.null(colnames(X)) && !is.null(colnames(V))) {
    if (anyDuplicated(colnames(X)) || anyDuplicated(colnames(V)) ||
        !setequal(colnames(X), colnames(V)))
      stop("Validation predictor names must match training predictor names.", call. = FALSE)
    V <- V[, colnames(X), drop = FALSE]
  }
  w <- if (is.null(validation$weights)) 1 else validation$weights
  off <- if (is.null(validation$offset)) 0 else validation$offset
  if (!is.numeric(w) || !length(w) || !length(w) %in% c(1L, length(v)) ||
      !is.numeric(off) || !length(off) || !length(off) %in% c(1L, length(v)))
    stop("Validation weights and offsets must be scalars or have one value per row.", call. = FALSE)
  w <- rep_len(w, length(v)); off <- rep_len(off, length(v))
  if (any(!is.finite(w) | w < 0) || !is.finite(sum(w)) || sum(w) <= 0 || any(!is.finite(off)))
    stop("Invalid validation weights or offsets.", call. = FALSE)
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
  if (intercept) V <- cbind(Intercept = 1, V)
  structure(list(X = V, y = v, weights = w / sum(w), offset = off,
    patience = patience, minit = minit, tol = tol), class = "L0glm_validation")
}

L0glm_validation_errors <- function(context, beta, shift = NULL, intercept = integer()) {
  beta <- as.numeric(beta)
  if (length(intercept)) beta[intercept] <- beta[intercept] - sum(shift * beta)
  active <- which(beta != 0)
  eta <- if (length(active) > length(beta) / 2) as.numeric(context$X %*% beta) else
    as.numeric(context$X[, active, drop = FALSE] %*% beta[active])
  (context$y - eta - context$offset)^2
}

L0glm_validation_loss <- function(context, beta, shift = NULL, intercept = integer()) {
  errors <- L0glm_validation_errors(context, beta, shift, intercept)
  keep <- context$weights > 0
  value <- sum(context$weights[keep] * errors[keep])
  if (is.finite(value)) value else Inf
}

L0glm_subset_validation <- function(context, columns) {
  if (!is.null(context)) context$X <- context$X[, columns, drop = FALSE]
  context
}

# Conditional marginal slopes avoid the global prediction calibration used by
# the default Gaussian initializer. They provide a second, training-only start.
L0glm_marginal_start <- function(X, y, weights, offset, intercept, always.include) {
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
  intercept <- isTRUE(args$intercept)
  supplied <- args$validation
  context <- L0glm_validation_context(supplied, X, y, family, intercept)
  starts <- if (is.null(supplied$starts)) c("default", "marginal") else supplied$starts
  if (!is.character(starts) || !length(starts) || anyNA(starts) ||
      any(!starts %in% c("default", "marginal")))
    stop('Validation starts must be "default", "marginal" or both.', call. = FALSE)
  if (!is.null(args$start)) starts <- "supplied"
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
    elapsed <- system.time(fit <- do.call(L0glm.fit, fit_args))[["elapsed"]]
    diagnostics[[key]] <<- data.frame(start = start,
      nsis = if (is.null(fit$screening)) NA_real_ else fit$screening$screen_size,
      validation_deviance = L0glm_validation_loss(context, coef(fit)),
      seconds = elapsed, selected = sum(coef(fit) != 0))
    design <- fit$X
    fit$X <- NULL
    models[[key]] <<- fit
    invisible(design)
  }
  null_args <- args
  null_args$validation <- context; null_args$start <- 0
  null_args$screen <- "none"; null_args$nsis <- NULL; null_args$relaxed <- FALSE
  design <- add(null_args, "null", "null")
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
  keep <- context$weights > 0
  errors <- L0glm_validation_errors(context, coef(models[[best]]))[keep]
  w <- context$weights[keep]; sum_w2 <- sum(w^2)
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
  result$validation$minimum <- table$validation_deviance[best]
  result$validation$se <- se
  result
}
