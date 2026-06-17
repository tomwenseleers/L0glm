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
#'   attr(simulated_spectra, "RT") and attr(simulated_spectra, "RI").
#' }
#'
#' @source Tom Wenseleers, derived from actually measured GC/MS spectra.
#' @author
#' Tom Wenseleers \email{tom.wenseleers@kuleuven.be} \\
#' Laboratory of Socioecology and Social Evolution, KU Leuven \\
#' \url{https://bio.kuleuven.be/eeb/tw/research} \\
#' ORCID: \href{https://orcid.org/0000-0002-1434-861X}{0000-0002-1434-861X}
#'
#' @example examples/simulated_spectra.R
"simulated_spectra"


#' GC/MS Chromatogram Section - Cuticular Hydrocarbon Profile of a Cataglyphis iberica Ant Queen
#'
#' A small 0.4 min section of a GC/MS chromatogram (RT 26.80-27.20 min) obtained from 
#' the cuticular hydrocarbon extract of a *Cataglyphis iberica* queen. 
#' Six compounds with near-Gaussian peak shapes co-elute in this segment. 
#' Data were collected on a Shimadzu QP-2010 Ultra 
#' equipped with a 30 m x 0.25 mm x 0.25 um DB-5ms column after 
#' liquid injection of a hexane extract. For details  
#' see Van Oystaeyen et al. Science 2014.
#'
#' @format A numeric matrix with 119 rows (time slices) and 671 columns (m/z values).
#' \describe{
#'   \item{Rows}{Individual full-scan acquisitions (time scans).
#'               Row names give the scan mid-point in minutes
#'               (26.80 <= RT <= 27.20).}
#'   \item{Columns}{Unit-mass channels, m/z 30:700.
#'                 Column names are the integer m/z values.}
#' Additional attribute:
#' \itemize{
#'   \item \code{attr(GCMS_data1, "alkanes_TIC")}: named numeric vector with
#'         the total ion chromatogram of a C9 to C40 alkane ladder standard run
#'         in the same sample batch; names give the retention time in minutes;
#'         the total ion trace for the entire chromatogram is included.
#' }
#' }
#' @details
#' Co-eluting compounds present with retention times and diagnostic ions: 
#'         \tabular{lll}{
#'           \strong{RT (min)} \tab \strong{Compound} \tab \strong{Diagnostic ions} \cr
#'           26.84 \tab 15-Methylnonacosane \tab 224, 225 \cr
#'           26.91 \tab 13-Methylnonacosane \tab 196, 225 \cr
#'           26.94 \tab 11-Methylnonacosane \tab 168, 280 \cr
#'           26.95 \tab Eicosamethylcyclodecasiloxane \tab 73, 147, 221, 281, 355, 429, 503 \cr
#'           27.01 \tab 9-Methylnonacosane \tab 140, 308 \cr
#'           27.11 \tab 7-Methylnonacosane \tab 112, 336
#'         }
#' @source Van Oystaeyen, A. et al. (2014) "Conserved class of queen pheromones
#'        stops social insect workers from reproducing." *Science* **343**,
#'        287-290. \doi{10.1126/science.1244899}
#' @author
#' Tom Wenseleers \email{tom.wenseleers@kuleuven.be} \\
#' Laboratory of Socioecology and Social Evolution, KU Leuven \\
#' \url{https://bio.kuleuven.be/eeb/tw/research} \\
#' ORCID: \href{https://orcid.org/0000-0002-1434-861X}{0000-0002-1434-861X}
#' @example examples/GCMS_data1.R
"GCMS_data1"


#' GC/MS Chromatogram Section - Beer Headspace Volatiles
#'
#' A 6 min segment (RT 19.00-25.00 min) of a GC/MS chromatogram recorded from
#' the SPME-headspace of beer. Data were acquired on a Thermo ISQ single-quad
#' coupled to a 30 m x 0.25 mm x 0.25 um RXi-5Sil MS column using a
#' DVB/Carboxen/PDMS (50/30 um) fibre. Many volatiles
#' co-elute in this interval, including by a prominent, triangular,
#' "shark-fin" octanoic-acid peak whose apex lies at RT 21.50 min. For details 
#' see Schreurs et al. Nature Communications 2024.
#'
#' @format A numeric matrix with 764 rows (time slices) and 518 columns (m/z values).
#' \describe{
#'   \item{Rows}{Individual full-scan acquisitions (time scans).
#'               Row names give the scan mid-point in minutes
#'               (19.00 <= RT <= 25.00).}
#'   \item{Columns}{Unit-mass channels, m/z 33:550.
#'                 Column names are the integer m/z values.}
#' Additional attribute:
#' \itemize{
#'   \item \code{attr(GCMS_data2, "alkanes_TIC")}: named numeric vector with
#'         the total ion chromatogram of a C7 to C22 alkane ladder standard run
#'         in the same sample batch; names give the retention time in minutes;
#'         the total ion trace for the chromatogram starting from a RT of
#'         6 mins is included.
#' }
#' }
#' @source Schreurs, M. et al. (2024) "Predicting and improving complex beer
#'        flavor through machine learning." *Nature Communications* **15**, 2368.
#'        \doi{10.1038/s41467-024-46346-0}
#' @author
#' Tom Wenseleers \email{tom.wenseleers@kuleuven.be} \\
#' Laboratory of Socioecology and Social Evolution, KU Leuven \\
#' \url{https://bio.kuleuven.be/eeb/tw/research} \\
#' ORCID: \href{https://orcid.org/0000-0002-1434-861X}{0000-0002-1434-861X}
#'
#' @example examples/GCMS_data2.R
"GCMS_data2"
