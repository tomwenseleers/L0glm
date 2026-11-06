#' Simulated Electron Impact Mass Spectra Dataset
#'
#' A dataset containing simulated electron impact mass spectra with a structure 
#' similar to as they would be encountered in GC/MS analysis.
#'
#' @format The data consists of a sparse \code{dgCMatrix} matrix with 
#' 2000 rows corresponding to different compounds and 258 columns, corresponding
#' to different mass fragments (m/z 41 to 298).
#' \describe{
#'   The 2000 rows correspond to different compounds and the 258 columns correspond 
#'   to different mass fragments (given as column names). Simulated retention
#'   times (in mins) and retention indices are given as attributes 
#'   attr(simulatedSpectra, "RT") and attr(simulatedSpectra, "RI").
#' }
#'
#' @source Tom Wenseleers, derived from actually measured GC/MS spectra.
#' @author
#' Tom Wenseleers \email{tom.wenseleers@kuleuven.be} \\
#' Laboratory of Socioecology and Social Evolution, KU Leuven \\
#' \url{https://bio.kuleuven.be/eeb/tw/research} \\
#' ORCID: \href{https://orcid.org/0000-0002-1434-861X}{0000-0002-1434-861X}
#'
#' @example examples/simulatedSpectra.R
"simulatedSpectra"


