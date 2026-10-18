# Internal utilities and benchmark plotting.

# clamp vector of coefficients between vectors with bounds lower and upper
clamp = function (coefs, lower, upper) {
  if (!is.null(coefs)) {
  return(pmax(lower, pmin(coefs, upper))) }
}
  


# AR.weights ####
# adaptive ridge penalty weights to approximate L0 norm penalty:
# INPUT
#   beta:   a vector of coefficents
#   delta:  a small constant
#   gamma:  the power value
# OUTPUT
#   AR.weights: a vector with adaptive weights associated to beta
# REF
#   Zou, 2006; Candes et al., 2008 & Nuel talk
#   de Rooi & Eilers 2011, Osborne et al. 2000 & Frommlet & Nuel 2016 (eqn 5)
AR.weights <- function(beta, gamma = 2, delta = 0, q = 0) {
  # The fitting algorithms explicitly use delta=0 (unsmoothed weights).
  # Positive delta remains available to callers of this internal helper.
  # Exact-zero coefficients are excluded by the fitting algorithm.
  if (delta == 0) 1 / abs(beta)^gamma else 1 / (abs(beta)^gamma + delta^gamma)
}

validate_thresh <- function(thresh) {
  if (!is.numeric(thresh) || length(thresh) != 1L || !is.finite(thresh) || thresh < 0)
    stop("thresh must be a single finite nonnegative number.", call. = FALSE)
  invisible(NULL)
}

####################################
#### EXPORTED UTILITY FUNCTIONS ####
####################################


#' Plot benchmark examples
#'
#' @importFrom graphics plot legend lines segments par
#'
#' @description
#'
#' Function to plot the comparison between the fitted model and the known true
#' model. It used in the \pkg{L0glm} examples and should not be expected
#' to useful in another use.
#'
#' @param x
#' a vector of length \code{n} with time labels.
#' @param y
#' a vector of length \code{n} with signal.
#' @param coefs XXX
#' @param fitted XXX
#' @param inference
#' the output of \code{L0glm.inference}
#' @param a.true
#' the true coefficients
#' @param ...
#' further arguments passed to \code{plot}
#'
#' @return
#'
#' The function returns \code{NULL}, but plots \code{y} with the expected
#' coefficients and the fitted values from \code{fit} with the estimated
#' coefficients in a way that gives an intuitive view on the quality of the fit.
#'
#' @export
plot_L0glm_benchmark <- function(x, y, coefs, fitted, inference = NULL, a.true, ...){
  a.fit <- coefs
  transf <- function(x) return(x)
  plot(x, transf(y), type="l", ylab="Signal", xlab="Time", col = "grey40",
       ylim=c(-transf(max(y)), transf(max(y))), ...)
  lines(x,-transf(y), col = "grey40")
  lines(transf(a.true), type = "h", col = "red", lwd = 2)
  lines(x, -transf(a.fit), type="h", col="blue", lwd = 2)
  lines(x, -transf(fitted), col="orange2", lwd = 1)
  if(!is.null(inference)){
    # Change color of significant covariates
    sign <- inference$p.value < 0.05
    lines(x[sign], -transf(a.fit[sign]), type="h", col="green4", lwd=2)
    # Draw CIs
    segments(x0 = x, x1 = x, y0 = -transf(inference$CI.lower), y1 = -transf(inference$CI.upper), col = "grey40")
    segments(x0 = x-0.35, x1 = x+0.35, y0 = -transf(inference$CI.lower), y1 = -transf(inference$CI.lower), col = "grey40")
    segments(x0 = x-0.35, x1 = x+0.35, y0 = -transf(inference$CI.upper), y1 = -transf(inference$CI.upper), col = "grey40")
    # Legend
    legend("topleft", lty = c(1,1,0,0,0), col = c("grey40", "orange2", "red", "blue", "green2"),
           pch = c("", "", "|", "|", "|"),
           legend = c("Input signal", "Fitted signal", "Ground truth",
                      "Non significant estimates (p-value >= 0.05)", "Significant estimates (p-value < 0.05)"))
  } else {
    legend("topleft", lty = c(1,1,0,0), col = c("grey40", "orange2", "red", "blue"),
           pch = c("", "", "|", "|"),
           legend = c("Input signal", "Fitted signal", "Ground truth",
                      "L0glm estimates"))
  }
}


# is matrix X sparse? ####
issparse <- function(X) {
  any(class(X) %in% c("dgCMatrix", "dsCMatrix", "sparseMatrix"))
}



gauspeak = function(x, u, w, h=1) h*exp(((x-u)^2)/(-2*(w^2)))

bigauspeak = function(x, u, w1, w2, h=1, center=TRUE) { # asymmetrical bigaussian peak with centroids centered at u
  ## analytic offset that moves the centroid onto u
  delta <- if (center) (w2 - w1) * sqrt(2 / pi) else 0
  u0    <- u - delta            # new "mode" position
  out <- numeric(length(x))
  left <- x <= u0
  out[left]  <- h * exp(- (x[left]  - u0)^2 / (2 * w1^2))
  out[!left] <- h * exp(- (x[!left] - u0)^2 / (2 * w2^2))
  return(out)
}

# centroid of a peak
centroid <- function(g) { 
  stopifnot(is.numeric(g), length(g) > 0)
  n       <- length(g)
  idx     <- seq_len(n)
  return(sum(idx * g) / sum(g))
}
