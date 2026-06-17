#' @name convmatrix
#' @title Make banded or rectangular convolution design matrices
#'
#' @description
#' \code{convmatrix()} builds a sparse Toeplitz (convolution) matrix
#' \strong{A} such that
#' \code{as.vector(A \%*\% x)} is numerically identical to
#' \code{gsignal::conv(x, kernel, shape)}.
#'
#' \itemize{
#'   \item If \code{n_cols} is \emph{omitted}, the function returns the
#'         classical \emph{banded} design matrix:
#'         square when \code{shape = "same"} and tall
#'         when \code{shape = "full"}.
#'   \item If \code{n_cols} is supplied, only the first
#'         \code{n_cols} time-shifts of \code{kernel} are kept, yielding an
#'         \code{n x n_cols} rectangular design.  Such matrices are useful
#'         for point-spread-function (PSF) update steps in
#'         blind-deconvolution algorithms.
#' }
#'
#' The matrix is produced in a single pass, so no intermediate (larger)
#' object is ever allocated, keeping memory usage to a minimum.
#'
#' @param kernel Numeric vector (`length L`) - FIR filter / PSF
#'               **in forward-time order**.
#' @param n      Integer \eqn{>} 0 - length of the observed signal (the number
#'               of rows for a *"same"* design).
#' @param shape  `"same"` (default, \eqn{n} rows) or
#'               `"full"` (\eqn{n+L-1} rows).
#' @param n_cols Optional.  Number of columns to keep.
#'               *If omitted* `n_cols <- n` and the result is the usual banded
#'               convolution matrix.
#' @param origin Integer (1 ... \eqn{L}) - element of `kernel` that is aligned
#'               with time 0.  Defaults to the centre for `"same"` and to 1
#'               for `"full"`.
#'
#' @return A sparse matrix of class **`dgCMatrix`** with
#'
#' * `n_rows = n`            (shape "same")   or
#' * `n_rows = n + L - 1`    (shape "full")
#'
#' rows, and `n_cols` columns.  
#'
#' @seealso
#'   \code{\link{convkernel}} to generate the 1-D PSF vectors used here.
#'
#' @example examples/convmatrix.R
#' @importFrom Matrix bandSparse
#' @export
#' @author
#' Tom Wenseleers \email{tom.wenseleers@kuleuven.be} \\
#' Laboratory of Socioecology and Social Evolution, KU Leuven \\
#' \url{https://bio.kuleuven.be/eeb/tw/research} \\
#' ORCID: \href{https://orcid.org/0000-0002-1434-861X}{0000-0002-1434-861X}
convmatrix <- function(kernel,
                       n,
                       shape   = c("same", "full"),
                       n_cols  = NULL,
                       origin  = NULL)
{
  shape  <- match.arg(shape)
  L      <- length(kernel)
  
  ## ----- default column count (banded design) ---------------------------
  if (is.null(n_cols))
    n_cols <- n
  stopifnot(n_cols >= 1L, n_cols == as.integer(n_cols))
  
  ## ----- default origin --------------------------------------------------
  if (is.null(origin))
    origin <- if (shape == "same") (L + 1L) %/% 2L else 1L
  stopifnot(origin >= 1L, origin <= L)
  
  ## ----- number of rows --------------------------------------------------
  n_row <- if (shape == "same") n else n + L - 1L
  
  ## ----- flip kernel so that  A %*% x  reproduces  conv(x, kernel) ------
  k_rev      <- rev(kernel)
  origin_rev <- L - origin + 1L
  
  ## ----- offsets of the diagonals we really need ------------------------
  k_off <- seq_len(L) - origin_rev
  keep  <- (k_off > -n_row) & (k_off < n_cols) & (k_rev != 0)
  
  k_off <- k_off[keep]
  diagv <- k_rev[keep]
  
  if (!length(k_off))
    return(Matrix::Matrix(0, n_row, n_cols, sparse = TRUE))
  
  ## ----- build each kept diagonal ---------------------------------------
  diagList <- lapply(seq_along(k_off), function(j) {
    d   <- k_off[j]
    len <- if (d >= 0) min(n_row, n_cols - d)
    else        min(n_row + d, n_cols)
    rep(diagv[j], len)
  })
  
  Matrix::bandSparse(n_row, n_cols,
                     k = k_off,
                     diagonals = diagList)
}
