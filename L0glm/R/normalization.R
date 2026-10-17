# Column scaling changes solver coordinates, not the specified ridge objective.
L0glm_normalization <- function(normalize, nonnegative = NULL) {
  if (is.logical(normalize) && length(normalize) == 1L && !is.na(normalize))
    return(if (normalize) "l2" else "none")
  if (!is.character(normalize) || length(normalize) != 1L || is.na(normalize) ||
      !normalize %in% c("auto", "none", "l2", "rms", "sd"))
    stop('normalize must be TRUE, FALSE, "auto", "none", "l2", "rms", or "sd".', call. = FALSE)
  if (normalize == "auto" && !is.null(nonnegative))
    return(if (nonnegative) "none" else "rms")
  normalize
}

L0glm_moments <- function(X, weights = NULL) {
  if (inherits(X, "Matrix")) {
    X <- as(X, "generalMatrix")
    if (inherits(X, "sparseMatrix")) X <- as(X, "CsparseMatrix")
  } else if (!is.double(X)) storage.mode(X) <- "double"
  L0glm_column_moments_cpp(X, if (is.null(weights)) NULL else as.numeric(weights))
}

L0glm_column_scale <- function(X, normalize) {
  mode <- L0glm_normalization(normalize, nonnegative = FALSE)
  if (mode == "none") return(rep(1, ncol(X)))
  n <- nrow(X)
  moments <- L0glm_moments(X)
  ss <- moments$ss
  if (mode == "rms") ss <- ss / n
  if (mode == "sd") {
    ss <- if (n > 1L) moments$centered_ss / (n - 1L) else
      rep(0, ncol(X))
  }
  scaling <- sqrt(ss)
  scaling[!is.finite(scaling) | scaling == 0] <- 1
  scaling
}

L0glm_preprocessing <- function(X, normalize, lower, always.include) {
  selectable <- setdiff(seq_len(ncol(X)), always.include)
  nonnegative <- !length(selectable) || all(lower[selectable] >= 0)
  free_intercept <- intersect(L0glm_unit_intercepts(X), which(lower < 0))
  mode <- L0glm_normalization(normalize, nonnegative && !length(free_intercept))
  center <- rep(0, ncol(X)); scale <- rep(1, ncol(X))
  if (mode != "none") {
    moments <- L0glm_moments(X)
    if (length(free_intercept)) center <- moments$mean
    ss <- if (length(free_intercept) || mode == "sd") moments$centered_ss else moments$ss
    divisor <- switch(mode, rms = nrow(X), sd = max(1, nrow(X) - 1L), 1)
    scale <- sqrt(ss / divisor)
    scale[!is.finite(scale) | scale == 0] <- 1
  }
  center[always.include] <- 0
  scale[always.include] <- 1
  list(mode = mode, center = center, scale = scale,
       p = length(selectable), centered = mode != "none" && length(free_intercept) > 0L)
}

L0glm_subset_preprocessing <- function(preprocessing, vars) {
  if (is.null(preprocessing)) return(NULL)
  preprocessing$center <- preprocessing$center[vars]
  preprocessing$scale <- preprocessing$scale[vars]
  preprocessing
}

# Change only the internal working origin. For sparse designs, shift only
# fully stored columns: these can have an arbitrarily large mean/SD ratio,
# and shifting their stored values adds no new entries. The remaining columns
# retain their sparse representation and are profiled with the intercept.
L0glm_center_working_design <- function(X, preprocessing, lower) {
  shift <- rep(0, ncol(X))
  intercept <- intersect(L0glm_unit_intercepts(X), which(lower < 0))
  if (!preprocessing$centered || !length(intercept))
    return(list(X = X, shift = shift, intercept = integer()))
  shift <- preprocessing$center
  if (inherits(X, "sparseMatrix")) {
    X <- as(as(X, "generalMatrix"), "CsparseMatrix")
    shift[diff(X@p) != nrow(X)] <- 0
    for (j in which(shift != 0)) {
      entries <- seq.int(X@p[j] + 1L, X@p[j + 1L])
      X@x[entries] <- X@x[entries] - shift[j]
    }
  } else X <- sweep(as.matrix(X), 2L, shift, "-")
  list(X = X, shift = shift, intercept = intercept[1L])
}
