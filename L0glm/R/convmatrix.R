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
#' @examples
#' ## 1. synthetic spike-train convolved with a bi-Gaussian PSF ------------
#' if (requireNamespace("gsignal", quietly = TRUE)) {
#'   library(Matrix)
#'   library(gsignal)
#'   
#'   ### ---- synthetic spike train ----------------------------------------
#'   set.seed(1)
#'   n <- 100
#'   a <- numeric(n)
#'   a[sample(n, 20)] <- rpois(20, 10)      # random spikes
#'   
#'   ### ---- bi-Gaussian PSF via kernel() ---------------------------------
#'   L <- 31
#'   g <- convkernel("bigaussian", L,
#'                   params = list(sigma1 = 1, sigma2 = 4))  # default p = 0.5
#'   plot(g, type = "l", main = "Bi-Gaussian PSF", ylab = "density")
#'   
#'   ### ---------- 1.1 FULL convolution -----------------------------------
#'   y_ref_full <- as.vector(conv(a, g, shape = "full"))        # reference
#'   
#'   A_full <- convmatrix(g, n, shape = "full")                 # sparse full conv matrix
#'   all.equal(as.matrix(A_full), convmtx(g, n), tol = 1e-12)   # matches convmtx
#'   
#'   y_full <- as.vector(A_full %*% a)                          # A  %*%  a
#'   
#'   plot(1:(n + L - 1), y_ref_full, type = "l",
#'        main = "FULL convolution", xlab = "time index", ylab = "signal")
#'   lines(1:(n + L - 1), y_full, col = 2, lty = 3)             # identical
#'   
#'   centre_shift <- (L + 1) / 2 - 1                            # PSF peak offset
#'   lines((1:n) + centre_shift - 1, a, type = "h", col = 4)    # spike locations
#'   all.equal(y_full, y_ref_full, tol = 1e-12)                 # matches conv
#'   
#'   ### ---------- 1.2 SAME-size convolution ------------------------------
#'   y_ref_same <- as.vector(conv(a, g, shape = "same"))
#'   
#'   A_same <- convmatrix(g, n, shape = "same")                 # sparse conv matrix
#'   y_same <- as.vector(A_same %*% a)
#'   
#'   plot(1:n, y_ref_same, type = "l",
#'        main = "SAME-length convolution", xlab = "time index", ylab = "signal")
#'   lines(1:n, y_same, col = 2, lty = 3)                       # identical
#'   lines((1:n) - 1, a, type = "h", col = 4)                   # spike locations
#'   all.equal(y_same, y_ref_same, tol = 1e-12)                 # matches conv
#' }
#' 
#' ## 2. rectangular design for PSF-update step in blinv deconv algo -------
#' set.seed(1)
#' n      <- 1000                                  # length of obs convolved signal
#' Lpsf   <- 101                                   # length of PSF to be estimated
#' spikes <- rpois(n, 2)                           # true or estimated spikes
#' 
#' G <- convmatrix(kernel = spikes,
#'                 n      = n,
#'                 n_cols = Lpsf,         # keep just the first Lpsf shifts
#'                 shape  = "same",
#'                 origin = (Lpsf + 1L) %/% 2L)
#' all.equal(            # identical to the first Lpsf columns of the banded design
#'   G,
#'   convmatrix(spikes, n, shape = "same",
#'              origin = (Lpsf + 1L) %/% 2L)[, 1:Lpsf, drop = FALSE], tol = 0)
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
