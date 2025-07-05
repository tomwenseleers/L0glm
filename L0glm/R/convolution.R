#' @name ConvMatrix_Banded
#' @title Make a banded convolution design matrix
#'
#' @description
#' Build a sparse banded Toeplitz matrix 'A' such that
#' 'as.vector(A %*% a)' is numerically identical to
#' the convolution 'gsignal::conv(a, kernel, shape)'.
#'
#' @param kernel Numeric vector (`length L`) with the FIR filter / PSF
#'               **in forward time order**.
#' @param n      Length of the signal to be convolved.
#' @param shape  `"same"` (default) – returns an `n × n` matrix whose product
#'               has length `n`;  
#'               `"full"` – returns an `(n + L − 1) × n` matrix that produces
#'               the full linear convolution of length `n + L − 1`.
#' @param origin *One-based* index (1 … L) of the kernel element regarded as
#'               the “origin” (i.e. the sample aligned with signal time
#'               \(t = 0\)). If omitted, the middle element (`(L+1)/2`) or the
#'               first element (1) are used as defaults with `shape=="same"` or
#'               `shape=="full"`.
#'
#' @return A sparse `n × n` (with shape=="same") or `(n + L − 1) × n`
#'         (with shape=="full") convolution design matrix of class `dgCMatrix`.
#' @examples
#' ## --- synthetic spike-train convolved with a bi-Gaussian PSF ----------
#' if (requireNamespace("gsignal", quietly = TRUE)) {
#'   library(Matrix)
#'   library(gsignal)
#'   
#'   ## ---- synthetic spike train ----------------------------------------
#'   set.seed(1)
#'   n <- 100
#'   a <- numeric(n)
#'   a[sample(n, 20)] <- rpois(20, 10)      # random spikes
#'   
#'   ## ---- bi-Gaussian PSF via kernel() ---------------------------------
#'   L <- 31
#'   g <- kernel("bigaussian", L,
#'               params = list(sigma1 = 1, sigma2 = 4))  # default p = 0.5
#'   plot(g, type = "l", main = "Bi-Gaussian PSF", ylab = "density")
#'   
#'   ## ---------- (1) FULL convolution -----------------------------------
#'   y_ref_full <- as.vector(conv(a, g, shape = "full"))       # reference
#'   
#'   A_full <- ConvMatrix_Banded(g, n, shape = "full")         # sparse design
#'   stopifnot(all.equal(as.matrix(A_full), convmtx(g, n),     # dense check
#'                       tol = 1e-12))
#'   
#'   y_full <- as.vector(A_full %*% a)                         # A  %*%  a
#'   
#'   plot(1:(n + L - 1), y_ref_full, type = "l",
#'        main = "FULL convolution", xlab = "time index", ylab = "signal")
#'   lines(1:(n + L - 1), y_full, col = 2, lty = 3)            # identical
#'   
#'   centre_shift <- (L + 1) / 2 - 1                           # PSF peak offset
#'   lines((1:n) + centre_shift - 1, a, type = "h", col = 4)   # spike locations
#'   stopifnot(all.equal(y_full, y_ref_full, tol = 1e-12))
#'   
#'   ## ---------- (2) SAME-size convolution ------------------------------
#'   y_ref_same <- as.vector(conv(a, g, shape = "same"))
#'   
#'   A_same <- ConvMatrix_Banded(g, n, shape = "same")
#'   y_same <- as.vector(A_same %*% a)
#'   
#'   plot(1:n, y_ref_same, type = "l",
#'        main = "SAME-length convolution", xlab = "time index", ylab = "signal")
#'   lines(1:n, y_same, col = 2, lty = 3)                       # identical
#'   lines((1:n) - 1, a, type = "h", col = 4)                   # spike locations
#'   stopifnot(all.equal(y_same, y_ref_same, tol = 1e-12))
#' }
#' 
#' @importFrom Matrix bandSparse
#' @seealso
#'   \code{\link{ConvMatrix_Rect}}  to make rectangular convolution design matrices,
#'   \code{\link{kernel}}           for generating 1-D PSF vectors.
#' @export
ConvMatrix_Banded <- function(kernel,
                              n,
                              shape  = c("same", "full"),
                              origin = NULL)
{
  shape <- match.arg(shape)
  L     <- length(kernel)
  
  ## -------- 1. choose origin --------------------------------------------
  if (is.null(origin)) {
    origin <- if (shape == "same")
      (L + 1L) %/% 2L
    else
      1L
  } else {
    stopifnot(origin >= 1L, origin <= L)
  }
  
  ## -------- 2. flip kernel *and* origin so A %*% a reproduces conv() -----
  kernel <- rev(kernel)
  origin <- L - origin + 1L                # map to flipped index
  
  ## -------- 3. matrix dimensions ----------------------------------------
  n_row <- if (shape == "same")
    n
  else
    n + L - 1L
  k.off <- seq_len(L) - origin             # diagonal offsets
  
  keep  <- kernel != 0
  k.off <- k.off[keep]
  diagv <- kernel[keep]
  
  if (!length(k.off))
    return(Matrix::Matrix(0, n_row, n, sparse = TRUE))
  
  diagList <- lapply(seq_along(k.off), function(j) {
    len <- min(n_row, n - k.off[j]) - max(0L, -k.off[j])
    rep(diagv[j], len)
  })
  
  Matrix::bandSparse(n_row, n, k         = k.off, diagonals = diagList)
}



#' @name ConvMatrix_Rect
#' @title Make a rectangular convolution design matrix
#'
#' @description
#' The matrix contains the first \code{n_cols} time-shifts of \code{kernel}.
#' It is exactly what you would obtain from
#' \preformatted{
#'   ConvMatrix_Banded(kernel, n, shape = shape, origin = ORIGIN)[ , 1:n_cols]
#' }
#' assuming you pass the **same** \code{ORIGIN}.  Here it is created
#' directly in a single pass, so no extra memory is allocated.
#'
#' @param kernel numeric vector (length \eqn{L}) – FIR filter / PSF
#'               **in forward–time order**.
#' @param n      integer \eqn{>} 0 – number of rows of the *“same”* design
#'               (length of the observed signal).
#' @param n_cols integer (1 … \eqn{L}) – how many columns to keep
#'               (often the target PSF length in a PSF-update step).
#' @param shape  `"same"` ( \eqn{n} rows ) or
#'               `"full"` ( \eqn{n+L-1} rows ).
#' @param origin integer (1 … \eqn{L}) – element of \code{kernel} that is
#'               aligned with time 0.  Defaults to the centre for
#'               `"same"` and to 1 for `"full"`.
#'
#' @return A sparse `n × n_cols` (with shape=="same") or `(n + L − 1) × n_cols`
#'         (with shape=="full") convolution design matrix of class `dgCMatrix`.
#'
#' @examples
#' ## rectangular convolution design (for PSF update in blind deconvolution algo)
#' set.seed(1)
#' n      <- 1e3
#' Lpsf   <- 101
#' spikes <- rpois(n, 2)
#'
#' ORIGIN <- (Lpsf + 1L) %/% 2L
#'
#' G <- ConvMatrix_Rect(spikes, n, Lpsf,
#'                      shape  = "same",
#'                      origin = ORIGIN)
#'
#' ## identical to the first Lpsf columns of ConvMatrix_Banded()
#' stopifnot(all.equal(
#'   G,
#'   ConvMatrix_Banded(spikes, n,
#'                     shape  = "same",
#'                     origin = ORIGIN)[ , 1:Lpsf, drop = FALSE],
#'   tol = 0))
#' @importFrom Matrix bandSparse
#' @seealso
#'   \code{\link{ConvMatrix_Banded}}  to make banded convolution design matrices,
#'   \code{\link{kernel}}             for generating 1-D PSF vectors.
#' @export
ConvMatrix_Rect <- function(kernel,
                            n,
                            n_cols,
                            shape  = c("same", "full"),
                            origin = NULL)
{
  shape  <- match.arg(shape)
  L      <- length(kernel)
  
  ## --- argument checks ---------------------------------------------------
  stopifnot(n       >= 1L,
            n       == as.integer(n),
            n_cols  >= 1L,
            n_cols  == as.integer(n_cols),
            n_cols <= L)
  
  ## --- default origin ----------------------------------------------------
  if (is.null(origin))
    origin <- if (shape == "same")
      (n_cols + 1L) %/% 2L
  else
    1L
  stopifnot(origin >= 1L, origin <= L)
  
  ## --- rows & diagonal offsets ------------------------------------------
  n_row <- if (shape == "same")
    n
  else
    n + L - 1L
  
  # flip kernel so that  A %*% x  reproduces  conv(x, kernel)
  k_rev <- rev(kernel)
  origin_rev <- L - origin + 1L          # index in the flipped kernel
  
  k_off <- seq_len(L) - origin_rev       # diagonal offsets
  keep  <- (k_off > -n_row) & (k_off < n_cols) & (k_rev != 0)
  
  k_off <- k_off[keep]
  diagv <- k_rev[keep]
  
  if (!length(k_off))
    return(Matrix::Matrix(0, n_row, n_cols, sparse = TRUE))
  
  ## --- build each kept diagonal -----------------------------------------
  diagList <- lapply(seq_along(k_off), function(j) {
    d   <- k_off[j]
    len <- if (d >= 0)
      min(n_row, n_cols - d)
    else
      min(n_row + d, n_cols)
    rep(diagv[j], len)
  })
  
  Matrix::bandSparse(n_row, n_cols, k = k_off, diagonals = diagList)
}
 