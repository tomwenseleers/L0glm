#' Simulated Spectra Dataset
#'
#' A dataset containing simulated mass spectra with a structure similar
#' to as they would be encountered in GC/MS analysis.
#'
#' @format The data consists of a sparse \code{dgCMatrix} matrix with 
#' 2000 rows corresponding to different compounds and 258 columns, corresponding
#' to different mass fragments (m/z 41 to 298).
#' \describe{
#'   The 2000 rows correspond to different compounds and the 258 columns correspond 
#'   to different mass fragments (given as column names). Simulated retention
#'   times (in mins) and retention indices are given as attributes 
#'   attr(simulated_spectra, "RT") and attr(simulated_spectra, "RI").
#' }
#'
#' @source Tom Wenseleers, derived from actually measured GC/MS spectra.
#'
#' @examples
#' \dontrun{
#' data(simulated_spectra)
#' }
"simulated_spectra"