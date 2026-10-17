# Column scaling changes solver coordinates, not the specified ridge objective.
L0glm_normalization <- function(normalize) {
  if (is.logical(normalize) && length(normalize) == 1L && !is.na(normalize))
    return(if (normalize) "l2" else "none")
  if (!is.character(normalize) || length(normalize) != 1L || is.na(normalize) ||
      !normalize %in% c("none", "l2", "rms", "sd"))
    stop('normalize must be TRUE, FALSE, "none", "l2", "rms", or "sd".', call. = FALSE)
  normalize
}

L0glm_column_scale <- function(X, normalize) {
  mode <- L0glm_normalization(normalize)
  if (mode == "none") return(rep(1, ncol(X)))
  n <- nrow(X)
  ss <- as.numeric(Matrix::colSums(X^2))
  if (mode == "rms") ss <- ss / n
  if (mode == "sd") {
    ss <- if (n > 1L) pmax(0, ss - as.numeric(Matrix::colSums(X))^2 / n) / (n - 1L) else
      rep(0, ncol(X))
  }
  scaling <- sqrt(ss)
  scaling[!is.finite(scaling) | scaling == 0] <- 1
  scaling
}
