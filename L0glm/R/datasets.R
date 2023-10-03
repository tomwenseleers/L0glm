#' Simulated Spectra Dataset
#'
#' A dataset containing simulated mass spectra with a structure similar
#' to as they would be encountered in GC/MS analysis.
#'
#' @format The data consists of a sparse \code{dgCMatrix} matrix with 
#' 100 rows corresponding to different compounds and 460 columns, corresponding
#' to different mass fragments.
#' \describe{
#'   The 100 rows correspond to different compounds and the 460 columns correspond 
#'   to different mass fragments.
#' }
#'
#' @source Tom Wenseleers, derived from actually measured GC/MS spectra.
#'
#' @examples
#' \dontrun{
#' data(simulated_spectra)
#' }
"simulated_spectra"