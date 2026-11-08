# Equivalent fits can be reached through different saved fitting specifications.
without_fit_specification <- function(fit) {
  fit$engine <- fit$specification <- NULL
  fit
}
