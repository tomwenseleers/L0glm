# Recognize an explicit intercept without constructing a dense logical copy of X.
L0glm_intercepts <- function(X) {
  named <- grep("intercept", colnames(X), ignore.case = TRUE)
  if (length(named) || !nrow(X)) return(named)
  candidates <- which(as.numeric(X[1L, , drop = FALSE]) == 1)
  for (j in candidates) if (all(X[, j, drop = FALSE] == 1)) return(j)
  integer()
}

L0glm_unit_intercepts <- function(X) {
  candidates <- L0glm_intercepts(X)
  candidates[vapply(candidates, function(j) all(X[, j, drop = FALSE] == 1), logical(1))]
}

L0glm_mandatory <- function(X, always.include = NULL) {
  sort(unique(c(intersect(as.integer(always.include), seq_len(ncol(X))),
                L0glm_intercepts(X))))
}

L0glm_predictor_count <- function(X, intercept, always.include = NULL) {
  # An added intercept leaves an unnamed constant input column selectable.
  included <- if (intercept) unique(c(always.include,
    grep("intercept", colnames(X), ignore.case = TRUE))) else L0glm_mandatory(X, always.include)
  ncol(X) - length(intersect(included, seq_len(ncol(X))))
}
