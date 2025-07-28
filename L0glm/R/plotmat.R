#' Plot a (Possibly Sparse) Matrix as a False‑Colour Image
#'
#' @description
#' Plot a numeric matrix (dense or sparse) as a raster image with a colour
#' gradient. The function optionally down‑scales large matrices to respect a
#' maximum resolution, applies a gamma transformation for contrast enhancement,
#' and maps the values to a user‑supplied or default colour palette.
#' Internally, the raster is drawn on an index grid (\code{x = 1:nrow}, \code{y = 1:ncol});
#' if numeric row/column names are available, custom axes are added afterwards
#' whose tick labels reflect those values.
#'
#' @usage
#' plotmat(M, maxresolution = 2000, gamma = 0.13, axes = TRUE,
#'         cols = grDevices::colorRampPalette(
#'           c("black","black","blue","green","yellow","orange","red"))(1000),
#'         ...)
#'
#' @param M A numeric matrix or a sparse matrix from the \pkg{Matrix} package
#'   (e.g., \code{dgCMatrix}).
#' @param maxresolution Integer. Maximum number of pixels in either dimension
#'   after down‑scaling. Set to \code{NULL} to disable down‑scaling. Default \code{2000}.
#' @param gamma Positive numeric. Gamma correction applied to matrix values
#'   prior to colour mapping (\eqn{v^\gamma}). Smaller values enhance contrast
#'   in low‑intensity regions. Default \code{0.13}.
#' @param axes Logical. If \code{TRUE} (default), draw custom axes after plotting
#'   on the index grid. If row/column names are numeric they are used to label
#'   ticks; otherwise default index ticks are used. Set \code{FALSE} to draw no axes.
#' @param cols A character vector of colours (e.g., created by
#'   \code{\link[grDevices]{colorRampPalette}}) passed to the plot. Defaults to a
#'   1000‑colour palette ranging from black to red.
#' @param ... Further arguments passed to \code{\link[graphics]{image}}
#'   (e.g., \code{xlab}, \code{ylab}, \code{main}, graphical parameters).
#'
#' @details
#' \itemize{
#'   \item \strong{Down‑scaling:} If the matrix exceeds \code{maxresolution}
#'     in either dimension, it is down‑scaled by block averaging (see \code{downscale()}).
#'   \item \strong{Sparse input:} Sparse matrices are converted to dense form
#'     via \code{\link[base]{as.matrix}}. Missing values are set to zero and
#'     negative values truncated to zero.
#'   \item \strong{Axes:} The image is rendered on the index grid
#'     \code{x = 1:nrow(M2)}, \code{y = 1:ncol(M2)} with \code{axes = FALSE} in
#'     the internal \code{\link[graphics]{image}} call; if \code{axes = TRUE},
#'     axes are added afterwards. When row/column names are numeric, tick marks
#'     are placed at positions corresponding to those values (using linear
#'     mapping from values to indices) and labeled accordingly; otherwise
#'     default index ticks are drawn.
#'   \item \strong{Colour mapping:} After gamma transformation the values are
#'     mapped to \code{cols}.
#' }
#'
#' @return
#' Invisibly returns the list produced by \code{\link[graphics]{image}}. The
#' primary purpose of the function is its graphical side‑effect.
#'
#' @seealso
#' \code{\link[graphics]{image}}, \code{\link[grDevices]{colorRampPalette}},
#' \code{\link[Matrix]{Matrix}}
#'
#' @examples
#' ## Dense matrix example: plot a simulated GC/MS chromatogram
#' M <- simulate_GCMS()$D
#' plotmat(M,
#'         gamma = 0.13,
#'         xlab = "Elution time (mins)",
#'         ylab = "Mass (m/z)",
#'         main = "GC/MS chromatogram")
#'
#' @importFrom grDevices colorRampPalette
#' @export
plotmat <- function(M,
                    maxresolution = 2000,
                    gamma = 0.13,
                    axes = TRUE,
                    cols = grDevices::colorRampPalette(
                      c("black","black","blue","green","yellow","orange","red"))(1000),
                    ...) {
  
  ## ----- Downscale (block means) ------------------------------------------
  M2 <- if (!is.null(maxresolution)) downscale(M, maxresolution) else M
  M2 <- as.matrix(M2)    # ensure dense for plotting
  M2[is.na(M2)] <- 0
  M2[M2 < 0]   <- 0
  
  ## ----- Axis value grids (for labeling only) ------------------------------
  # Row names -> x-axis values (length = nrow(M2))
  if (!is.null(as.numeric(rownames(M))[[1]])) {
    xrange <- range(as.numeric(rownames(M)))
    xvals  <- seq(min(xrange), max(xrange), length.out = nrow(M2))
  } else {
    xvals  <- seq_len(nrow(M2))
  }
  # Column names -> y-axis values (length = ncol(M2))
  if (!is.null(as.numeric(colnames(M))[[1]])) {
    yrange <- range(as.numeric(colnames(M)))
    yvals  <- seq(min(yrange), max(yrange), length.out = ncol(M2))
  } else {
    yvals  <- seq_len(ncol(M2))
  }
  
  ## ----- Plot on index grid, suppress default axes -------------------------
  out <- graphics::image(
    x = seq_len(nrow(M2)),
    y = seq_len(ncol(M2)),
    z = M2 ^ gamma,
    col = cols,
    useRaster = TRUE,
    axes = FALSE,
    ...
  )
  
  ## ----- Custom axes (map values -> index positions) -----------------------
  if (isTRUE(axes)) {
    # choose 'pretty' tick values in value space
    xticks_val <- pretty(xvals, n = 6)
    yticks_val <- pretty(yvals, n = 6)
    
    # restrict to range; map to index coordinates via linear interpolation
    xticks_val <- xticks_val[xticks_val >= min(xvals) & xticks_val <= max(xvals)]
    yticks_val <- yticks_val[yticks_val >= min(yvals) & yticks_val <= max(yvals)]
    
    xticks_at <- approx(x = xvals, y = seq_along(xvals),
                        xout = xticks_val, rule = 2)$y
    yticks_at <- approx(x = yvals, y = seq_along(yvals),
                        xout = yticks_val, rule = 2)$y
    
    graphics::axis(side = 1, at = xticks_at, labels = xticks_val)
    graphics::axis(side = 2, at = yticks_at, labels = yticks_val, las = 1)
    graphics::box()
  }
  
  invisible(out)
}


# ---------------------------------------------------------------------------
# Helper: down‑scale a matrix to ≤ maxresolution × maxresolution by
# block‑averaging. Works for dense or sparse input.
# ---------------------------------------------------------------------------
downscale <- function(M, maxresolution = 2000) {
  rstride <- ceiling(nrow(M) / min(nrow(M), maxresolution))
  cstride <- ceiling(ncol(M) / min(ncol(M), maxresolution))
  
  M <- as.matrix(M)
  
  row_grp <- (seq_len(nrow(M)) - 1L) %/% rstride + 1L
  rsizes  <- tabulate(row_grp)
  M_rowmean <- rowsum(M, row_grp, reorder = FALSE)
  M_rowmean <- M_rowmean / rsizes
  
  col_grp <- (seq_len(ncol(M_rowmean)) - 1L) %/% cstride + 1L
  csizes  <- tabulate(col_grp)
  M_colmean <- rowsum(t(M_rowmean), col_grp, reorder = FALSE)
  M_ds <- t(M_colmean / csizes)
  
  return(M_ds)
}
