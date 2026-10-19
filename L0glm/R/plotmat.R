#' Plot a (Possibly Sparse) Matrix as a False-Colour Image
#'
#' @description
#' Plot a numeric matrix (dense or sparse) as a raster image with a colour
#' gradient. The function optionally down-scales large matrices to respect a
#' maximum resolution, applies a gamma transformation for contrast enhancement,
#' and maps the values to a user-supplied or default colour palette.
#' Internally, the image is drawn on an index grid (\code{x = 1:nrow}, \code{y = 1:ncol});
#' if numeric row/column names are available, custom axes are added afterwards
#' and tick labels reflect those values.
#'
#' @usage
#' plotmat(M, maxresolution = 2000, gamma = 0.13, axes = TRUE,
#'         cols = grDevices::colorRampPalette(
#'           c("black","black","blue","green","yellow","orange","red"))(1000),
#'         signed = FALSE, ...)
#'
#' @param M A numeric matrix or a sparse matrix from the \pkg{Matrix} package
#'   (e.g., \code{dgCMatrix}).
#' @param maxresolution Integer. Maximum number of pixels in either dimension
#'   after down-scaling. Set to \code{NULL} to disable down-scaling. Default \code{2000}.
#' @param gamma Positive numeric. Gamma correction applied to matrix values
#'   prior to colour mapping (\eqn{v^\gamma}). Smaller values enhance contrast
#'   in low-intensity regions. Default \code{0.13}.
#' @param axes Logical. If \code{TRUE} (default), draw custom axes after plotting
#'   on the index grid. If row/column names are fully numeric they are used to
#'   label ticks; otherwise default index ticks are used. Set \code{FALSE} to draw no axes.
#' @param cols A character vector of colours (e.g., created by
#'   \code{\link[grDevices]{colorRampPalette}}). Defaults to a 1000-colour palette.
#' @param signed Preserve negative values and use the signed gamma transform
#'   \eqn{\operatorname{sign}(v)|v|^\gamma}. Default \code{FALSE} truncates
#'   negative values to zero. Supply a diverging palette for signed data.
#' @param ... Further arguments passed to \code{\link[graphics]{image}}
#'   (e.g., \code{xlab}, \code{ylab}, \code{main}, graphical parameters).
#'
#' @details
#' \itemize{
#'   \item \strong{Down-scaling:} If the matrix exceeds \code{maxresolution}
#'     in either dimension, it is down-scaled by block averaging (see \code{downscale()}).
#'   \item \strong{Sparse input:} Sparse matrices are converted to dense form
#'     via \code{\link[base]{as.matrix}}. Missing values are set to zero and
#'     negative values truncated to zero unless \code{signed=TRUE}.
#'   \item \strong{Axes:} The image is rendered on the index grid with
#'     \code{axes = FALSE} in the internal \code{\link[graphics]{image}} call.
#'     When \code{axes = TRUE}, tick locations are computed by mapping "value
#'     space" (from numeric row/column names, if available) to index positions
#'     via linear interpolation; if names are missing or not all numeric,
#'     integer index ticks are used.
#' }
#'
#' @return
#' Invisibly returns the list produced by \code{\link[graphics]{image}}. The
#' primary purpose of the function is its graphical side-effect.
#'
#' @author
#' Tom Wenseleers \email{tom.wenseleers@kuleuven.be} \\
#' Laboratory of Socioecology and Social Evolution, KU Leuven \\
#' \url{https://bio.kuleuven.be/eeb/tw/research} \\
#' ORCID: \href{https://orcid.org/0000-0002-1434-861X}{0000-0002-1434-861X}
#' 
#' @seealso
#' \code{\link[graphics]{image}}, \code{\link[grDevices]{colorRampPalette}}, \code{\link[Matrix]{Matrix}}
#'
#'
#' @example examples/plotmat.R
#' @importFrom grDevices colorRampPalette
#' @export
plotmat <- function(M,
                    maxresolution = 2000,
                    gamma = 0.13,
                    axes = TRUE,
                    cols = grDevices::colorRampPalette(
                      c("black","black","blue","green","yellow","orange","red"))(1000),
                    signed = FALSE,
                    ...) {
  
  # helper: generate axis "value" grid from dimnames if fully numeric, else index
  axis_values_from_names <- function(nms, n_out) {
    if (is.null(nms)) return(seq_len(n_out))
    v <- suppressWarnings(as.numeric(nms))
    if (length(v) == length(nms) && all(!is.na(v))) {
      seq(min(v), max(v), length.out = n_out)
    } else {
      seq_len(n_out)
    }
  }
  
  ## ----- Downscale (block means) ------------------------------------------
  M2 <- if (!is.null(maxresolution)) downscale(M, maxresolution) else M
  M2 <- as.matrix(M2)
  M2[is.na(M2)] <- 0
  if (!signed) M2[M2 < 0] <- 0
  
  ## ----- Axis "value" vectors for labelling only ---------------------------
  xvals <- axis_values_from_names(rownames(M), nrow(M2))
  yvals <- axis_values_from_names(colnames(M), ncol(M2))
  
  ## ----- Plot on index grid, suppress default axes -------------------------
  out <- graphics::image(
    x = seq_len(nrow(M2)),
    y = seq_len(ncol(M2)),
    z = if (signed) sign(M2) * abs(M2)^gamma else M2^gamma,
    col = cols,
    useRaster = TRUE,
    axes = FALSE,
    ...
  )
  
  ## ----- Custom axes (map values -> index positions) -----------------------
  if (isTRUE(axes)) {
    xticks_val <- pretty(xvals, n = 6)
    yticks_val <- pretty(yvals, n = 6)
    
    xticks_val <- xticks_val[xticks_val >= min(xvals) & xticks_val <= max(xvals)]
    yticks_val <- yticks_val[yticks_val >= min(yvals) & yticks_val <= max(yvals)]
    
    xticks_at <- if (length(xvals) == 1L) rep(1, length(xticks_val)) else
      approx(x = xvals, y = seq_along(xvals), xout = xticks_val, rule = 2)$y
    yticks_at <- if (length(yvals) == 1L) rep(1, length(yticks_val)) else
      approx(x = yvals, y = seq_along(yvals), xout = yticks_val, rule = 2)$y
    
    graphics::axis(side = 1, at = xticks_at, labels = xticks_val)
    graphics::axis(side = 2, at = yticks_at, labels = yticks_val, las = 1)
    graphics::box()
  }
  
  invisible(out)
}


# ---------------------------------------------------------------------------
# Helper: down-scale a matrix to <= maxresolution x maxresolution by
# block-averaging. Works for dense or sparse input.
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
