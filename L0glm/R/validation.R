# Validation uses original predictor/response units, independently of the
# training normalization. Only the training observations enter ridge solves.
L0glm_validation_context <- function(validation, X, y, family, intercept) {
  if (family$family != "gaussian" || family$link != "identity" || NCOL(y) != 1L)
    stop("Validation currently supports one Gaussian identity-link response.", call. = FALSE)
  fields <- c("X", "y", "weights", "offset", "patience", "minit", "tol")
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
