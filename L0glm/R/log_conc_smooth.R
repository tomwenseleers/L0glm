#' @name log_conc_smooth
#' @title Log-concave spline smoothing of a discrete convolution kernel estimate
#'
#' @description
#' Fits a concave quadratic B-spline (via \code{\link[cobs]{cobs}}) to
#' \eqn{\log_{10}(y + \varepsilon)} and back-transforms and normalises the 
#' fit so the result is a **log-concave, non-negative, unit L1 norm** smoothed
#' version of the original convolution kernel estimate \eqn{y}.  
#' A small offset is added internally to avoid taking the
#' logarithm of zero; the final vector is renormalised to unit L1 norm.
#'
#' @param y      numeric vector (length \eqn{n \ge 3}) with the raw,
#'               non-negative data sampled at the integer points
#'               \code{1:length(y)}.
#' @param tau    quantile to be fitted by \code{cobs}; defaults to
#'               \code{0.5} (median fit).
#' @param nknots integer; number of interior knots for the spline.
#'               If \code{NULL} (default) it is set to \code{round(length(y)/5)}.
#'
#' @return A numeric vector of length \code{length(y)}:
#' \itemize{
#'   \item strictly non-negative,
#'   \item log-concave over its support,
#'   \item rescaled to unit L1 norm.
#' }
#'
#' @details
#' The weights passed to \code{cobs()} are proportional to
#' \code{y + offs}, where \code{offs = max(y) * 1e-6}; in practice that gives
#' the most stable fits for sparse or strongly peaked inputs.
#'
#' @seealso
#' \code{\link[cobs]{cobs}} for the underlying constrained B-spline
#' regression.
#' 
#' @author
#' Tom Wenseleers \email{tom.wenseleers@kuleuven.be} \\
#' Laboratory of Socioecology and Social Evolution, KU Leuven \\
#' \url{https://bio.kuleuven.be/eeb/tw/research} \\
#' ORCID: \href{https://orcid.org/0000-0002-1434-861X}{0000-0002-1434-861X}
#'
#' @example examples/log_conc_smooth.R
#' @importFrom cobs cobs
#' @importFrom stats  approx na.pass
#' @export
log_conc_smooth <- function(y, 
                            tau    = 0.5, # fit median
                            nknots = NULL) {
  
  n <- length(y)
  if (n < 3L) stop("Need at least three points (n >= 3).")
  if (!is.vector(y)) y <- as.vector(y)
  if (is.null(nknots)) nknots <- round(length(y)/5)
  
  x_raw <- seq_len(n)
  
  offs <- max(y) * 1e-6
  weights <- y + offs # empirically gave best convergence
  weights <- weights/sum(weights)*length(weights)
  
  ## ---- 1.  primary cobs fit --------------------------------------------
  fit <- cobs::cobs(x = x_raw,
                    y = log10(y + offs),
                    constraint = "concave",
                    w = weights, 
                    knots = seq(min(x_raw), max(x_raw),
                                length.out = nknots),
                    nknots = nknots,
                    knots.add = FALSE, # adding extra knots disabled for stability
                    repeat.delete.add = FALSE, # removing unneeded knots disabled for stability
                    tau = tau, 
                    degree = 2,
                    lambda = 0)
  
  eval_y <- function(x) {
    10^predict(fit, z = x, na.action = stats::na.pass)[, "fit"] + offs
  }
  
  y_fit <- eval_y(x_raw)
  
  ## ---- 2.  post-processing & output -----------------------------------
  y_final <- y_fit - min(y_fit, na.rm = TRUE)
  if (!all(y_final == 0)) y_final <- y_final / sum(y_final)
  
  return(y_final)
}



