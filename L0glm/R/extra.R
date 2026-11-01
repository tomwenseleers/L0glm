#' @useDynLib L0glm
#' @importFrom Rcpp evalCpp
#' @import Matrix
#' @importFrom graphics layout
#' @importFrom stats as.dist pnorm predict quantile rlnorm var
NULL

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("simulated_spectra"))
}
