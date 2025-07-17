#' Plot a (Possibly Sparse) Matrix as a False-Colour Image
#'
#' @description
#' Plot a numeric matrix (dense or sparse) as a raster image with a colour
#' gradient. The function optionally downscales large matrices to respect a
#' maximum resolution, applies a gamma transformation for contrast enhancement,
#' and maps the values to a user-defined or default colour palette.
#'
#' @usage
#' plotmat(M, maxresolution = 2000,
#'         row_subsamp = round(nrow(M)/min(nrow(M), maxresolution)),
#'         col_subsamp = round(ncol(M)/min(ncol(M), maxresolution)),
#'         gamma = 0.12,
#'         cols = grDevices::colorRampPalette(
#'           c("black","black","blue","green","yellow","orange","red"))(1000),
#'         ...)
#'
#' @arguments
#' \describe{
#'   \item{\code{M}}{A numeric matrix or a sparse matrix from the \pkg{Matrix}
#'   package (e.g., \code{dgCMatrix}).}
#'
#'   \item{\code{maxresolution}}{Integer. Maximum number of pixels in either
#'   dimension after downscaling. Default is \code{2000}.} Ignored if NULL.
#'
#'   \item{\code{gamma}}{Numeric > 0. Gamma correction applied to matrix values
#'   prior to colour mapping. Values are transformed as
#'   \eqn{v^\gamma}. Smaller values enhance contrast in low-intensity regions.
#'   Default is \code{0.13}.}
#'
#'   \item{\code{cols}}{Character vector of colours to be used by
#'   \code{\link[graphics]{image}}. Defaults to a 1,000-colour palette from
#'   black to red.}
#'
#'   \item{\code{...}}{Additional arguments passed to
#'   \code{\link[graphics]{image}} (e.g., \code{axes = FALSE}, \code{xlab},
#'   \code{ylab}, \code{main}).}
#' }
#'
#' @details
#' \itemize{
#'   \item \strong{Downscaling:} If the matrix is larger than
#'   \code{maxresolution} in either dimension, it is downscaled to fit.
#'
#'   \item \strong{Sparse input:} Sparse matrices are converted to dense form
#'   via \code{\link[base]{as.matrix}}. Missing values are set to zero and
#'   negative values truncated to zero.
#'
#'   \item \strong{Axis scaling:} If row or column names can be coerced to
#'   numeric, they are used for the \code{x} and \code{y} axes; otherwise,
#'   integer indices are used.
#'
#'   \item \strong{Colour mapping:} A gamma transformation
#'   (\eqn{z^\gamma}) is applied before mapping to the colour palette.
#' }
#'
#' @value
#' Invisibly returns the list produced by \code{\link[graphics]{image}}. The
#' main purpose of the function is the graphical output.
#'
#' @author
#' Tom Wenseleers \email{tom.wenseleers@kuleuven.be} \\
#' Laboratory of Socioecology and Social Evolution, KU Leuven \\
#' \url{https://bio.kuleuven.be/eeb/tw/research} \\
#' ORCID: \href{https://orcid.org/0000-0002-1434-861X}{0000-0002-1434-861X}
#'
#' @seealso
#' \code{\link[graphics]{image}}, \code{\link[grDevices]{colorRampPalette}},
#' \code{\link[Matrix]{Matrix}}
#'
#' @examples
#' ## Dense matrix example: plot a GC/MS chromatogram
#' M = simulate_GCMS()$D # simulated GC/MS chromatogram
#' plotmat(M = M, 
#'         gamma = 0.13, 
#'         xlab = "Elution time (mins)", 
#'         ylab = "Mass (m/z)",
#'         main = "GC/MS chromatogram")
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics image
#' @export
plotmat <- function(M,
                    maxresolution = 2000,
                    gamma = 0.13,
                    cols = grDevices::colorRampPalette(
                      c("black", "black", "blue", "green",
                        "yellow", "orange", "red"))(1000),
                    ...) {
  
  # Downscale
  if (!is.null(maxresolution)) { M2 <- downscale(M, maxresolution = maxresolution) } else { M2 <- M }
  M2[is.na(M2)] <- 0
  M2[M2 < 0] <- 0
  
  # Axis coordinates
  if (!is.null(as.numeric(rownames(M))[[1]])) {
    xrange <- range(as.numeric(rownames(M)))
    xvals <- seq(min(xrange), max(xrange), length.out = nrow(M2))
  } else {
    xvals <- seq_len(nrow(M2))
  }
  if (!is.null(as.numeric(colnames(M))[[1]])) {
    yrange <- range(as.numeric(colnames(M)))
    yvals <- seq(min(yrange), max(yrange), length.out = ncol(M2))
  } else {
    yvals <- seq_len(ncol(M2))
  }
  
  # Plot
  graphics::image(
    x = xvals,
    y = yvals,
    z = M2 ^ gamma,
    col = cols,
    useRaster = TRUE,
    ...
  )
}


# ---------------------------------------------------------------------------
# Helper: down‑scale a matrix to ≤ maxresolution × maxresolution by
# block‑averaging.  Works for dense or sparse input.
# ---------------------------------------------------------------------------
downscale <- function(M, maxresolution = 2000) {
  
  ## ---- strides (block sizes) ---------------------------------------------
  rstride <- ceiling(nrow(M) / min(nrow(M), maxresolution))
  cstride <- ceiling(ncol(M) / min(ncol(M), maxresolution))
  
  ## ---- to dense if needed -------------------------------------------------
  M <- as.matrix(M)              # sparse → dense once
  
  ## ---- row grouping factor -----------------------------------------------
  row_grp <- (seq_len(nrow(M)) - 1L) %/% rstride + 1L
  rsizes  <- tabulate(row_grp)                 # group sizes (for means)
  
  ## ---- average rows ------------------------------------------------------
  # rowsum() returns ∑ rows per group; divide by group size for the mean
  M_rowmean <- rowsum(M, row_grp, reorder = FALSE)
  M_rowmean <- M_rowmean / rsizes              # vectorised recycling
  
  ## ---- column grouping factor --------------------------------------------
  col_grp <- (seq_len(ncol(M_rowmean)) - 1L) %/% cstride + 1L
  csizes  <- tabulate(col_grp)
  
  ## ---- average columns ----------------------------------------------------
  # repeat trick on the transposed matrix
  M_colmean <- rowsum(t(M_rowmean), col_grp, reorder = FALSE)
  M_ds <- t(M_colmean / csizes)                # back to (rows × cols)
  
  return(M_ds)
}
