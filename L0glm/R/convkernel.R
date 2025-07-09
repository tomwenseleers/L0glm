#' @name convkernel
#' @title Generate normalised 1-D convolution kernels (point-spread functions)
#'
#' @description
#' Generates a numeric vector of length \emph{L} that represents a
#' point–spread function (PSF) or impulse–response kernel.
#' The kernel is \eqn{L^{1}}-normalised (area = 1).
#' Optionally it is shifted with \code{\link{shift_kernel}} so that a
#' user-chosen percentile \eqn{perc} of the cumulative distribution falls on
#' index \eqn{1 + perc \cdot (L - 1)}; setting \code{perc = 0} leaves the
#' kernel unshifted.  Alternatively, using \code{perc = "mean"} will move the
#' centroid to the middle of the kernel.
#'
#' The helper is useful for simulating blurred spike trains as observed in
#' one-dimensional spectroscopy (GC–MS, NMR), calcium neuro-imaging, or
#' epidemic curves where diagnosed cases are obtained by convolving
#' infections with a delay-distribution kernel.  The available kernel
#' families are:
#'
#' \itemize{
#'   \item \strong{Gaussian / Bigaussian / Exponentially Modified
#'         Gaussian (Gaussian * Exp)} – common peak-shape models in
#'         gas-chromatography / mass-spectrometry.
#'   \item \strong{Lorentzian} – arises in 1-D NMR as the Fourier transform
#'         of an exponentially decaying free-induction decay
#'         (Cauchy–Lorentz distribution).
#'   \item \strong{Exponential} – standard calcium-indicator impulse
#'         response in neuro-imaging studies.
#'   \item \strong{Gamma / Gamma * Exp} – typical
#'         infection-to-diagnosis delay distributions in epidemiology and
#'         sometimes used for chromatographic peak shapes.
#' }
#'
#' All kernels returned by \code{convkernel()} integrate to 1 and are ready
#' for use in convolution or deconvolution pipelines.
#'
#' @param type  One of `"gaussian"`, `"bigaussian"`, `"emg"` (exponentially
#'   modified Gaussian), `"lorentzian"` (Cauchy), `"laplace"` (two-sided
#'   exponential), `"exponential"` (one-sided exponential), `"gamma"`,
#'   `"gammaexp"` (convolution of Gamma with exponential).
#' @param L     Integer (≥ 3).  Kernel length.
#' @param params Named list of shape parameters (see Details).
#' @param perc Numeric in [0, 1] or "mean".  
#'   Position the chosen percentile of the kernel’s CDF at index
#'   `1 + perc·(L – 1)` or `mean` for centroid centering.  
#'   Defaults to `0` for causal families
#'   (`exponential`, `gamma`, `gammaexp`) and to `mean`
#'   for all other types (resulting in centroid centering).
#'
#' @details Default parameters for the different distributions are as
#'   follows:\cr gaussian: `sigma=2.5`\cr bigaussian:  `sigma1=1`, `sigma2=5`\cr
#'   emg (gaussian * exponential):`sigma=2`,`tau=6 ` lorentzian: `gamma=2`\cr
#'   laplace: `tau=3`\cr exponential: `tau=3`\cr gamma: `k=9`, `tau=0.5`\cr
#'   gammaexp (gamma * exponential): `k=9` , ` tau_g=0.7` , `tau_e=0.7` \cr
#'
#' Returned vectors are **L¹-normalised**.  
#' If `align_percentile > 0` the kernel is shifted (by sub-sample
#' interpolation) so that the requested percentile lies at its
#' theoretical position inside the length-*L* vector.
#'
#' @return Numeric vector of length *L*.
#'
#' @examples
#' ## available convolution kernels
#' cols <- 1:8
#' k1 <- convkernel("gaussian",    51, list(sigma = 4))
#' k2 <- convkernel("bigaussian",  51, list(sigma1 = 2, sigma2 = 6))
#' k3 <- convkernel("emg",         51, list(sigma = 2, tau = 6))
#' k4 <- convkernel("lorentzian",  51, list(gamma = 3))
#' k5 <- convkernel("laplace",     51, list(tau = 4))
#' k6 <- convkernel("exponential", 51, list(tau = 5))
#' k7 <- convkernel("gamma",       51, list(k = 9, tau = 0.5))
#' k8 <- convkernel("gammaexp",    51, list(k = 9, tau_g = 0.7, tau_e = 0.7))
#' matplot(cbind(k1,k2,k3,k4,k5,k6,k7,k8), type = "l", lty = 1, col = cols,
#'         xlab = "kernel index", ylab = "density",
#'         main = "Convolution kernels\n(point–spread functions)")
#' legend("topright",
#'        legend = c("Gaussian", "Bigaussian", "EMG", "Lorentzian",
#'                   "Laplace", "Exponential", "Gamma", "Gammaexp"),
#'        lty = 1, cex = 0.8, y.intersp = 0.8,
#'        col = cols, bty = "n")
#' 
#' ## peak broadening as observed at higher concentrations in GC/MS analysis
#' cols=1:3
#' matplot(cbind(
#'   convkernel("bigaussian", 51,
#'              list(sigma1 = 2, sigma2 = 4), perc = 0.001),
#'   convkernel("bigaussian", 51,
#'              list(sigma1 = 2*2, sigma2 = 4*2), perc = 0.001),
#'   convkernel("bigaussian", 51,
#'              list(sigma1 = 2*3, sigma2 = 4*3), perc = 0.001)), 
#'   xlab = "kernel index", ylab = "density",
#'   main = "Bigaussian point spread functions\n(with sigma1=2, 4 and 6 and sigma2=2*sigma1\nand perc = 0.001)",
#'   col = cols, lty = 1,
#'   type="l") 
#' @importFrom gsignal conv
#' @importFrom stats   approx
#' @seealso
#'   \code{\link{shift_kernel}}, \code{\link{convmatrix}}
#' @export
convkernel <- function(type = c("gaussian","bigaussian","emg","lorentzian",
                                "laplace","exponential","gamma","gammaexp"),
                   L,
                   params  = list(),
                   perc = NULL)     
{
  `%||%` <- function(a,b) if (is.null(a)) b else a
  type <- match.arg(type)
  stopifnot(L >= 3L, L == as.integer(L))
  
  ## ---------------- choose default alignment ----------------------------
  if (is.null(perc)) {
    causal  <- type %in% c("exponential","gamma","gammaexp")
    perc <- if (causal) 0 else "mean"
  }
  stopifnot((is.numeric(perc) && 0 <= perc && perc <= 1) || is.character(perc))
  
  idx <- seq_len(L)
  x   <- idx - (L + 1)/2        # symmetric grid for analytic formulas
  
  ## ---- analytic kernels -------------------------------------------------
  kernel <- switch(type,
              
              gaussian = {
                s <- params$sigma %||% 2.5
                1/(sqrt(2*pi)*s) * exp(-x^2/(2*s^2))
              },
              
              bigaussian = {
                s1 <- params$sigma1 %||% 1
                s2 <- params$sigma2 %||% 5
                delta <- (s2 - s1) * sqrt(2/pi)
                x0 <- x + delta
                out <- ifelse(x0 <= 0,
                              exp(-x0^2/(2*s1^2)),
                              exp(-x0^2/(2*s2^2)))
                out / (sqrt(pi/2)*(s1 + s2))
              },
              
              emg = {
                s   <- params$sigma %||% 2
                tau <- params$tau   %||% 6
                λ   <- 1/tau
                pref <- λ/2 * exp(λ*(λ*s^2/2 - x))
                pref * 2*pnorm(-(λ*s^2 - x)/(s*sqrt(2)))
              },
              
              lorentzian = {
                γ <- params$gamma %||% 2
                γ / pi / (x^2 + γ^2)
              },
              
              laplace = {
                tau <- params$tau %||% 3
                1/(2*tau) * exp(-abs(x)/tau)
              },
              
              exponential = {
                tau <- params$tau %||% 3
                out <- numeric(L)
                t   <- idx - 1                       # t = 0,1,2,...
                out <- 1/tau * exp(-t/tau)
                out
              },
              
              gamma = {
                k   <- params$k   %||% 9
                tau <- params$tau %||% 0.5
                t   <- idx - 1
                rate <- 1/tau
                out <- rate^k / gamma(k) * t^(k-1) * exp(-rate*t)
                out
              },
              
              gammaexp = {
                k      <- params$k      %||% 9
                tau_g  <- params$tau_g  %||% 0.7   # Gamma scale
                tau_e  <- params$tau_e  %||% 0.7   # Exp   scale
                
                t  <- idx - 1                     # 0,1,2,...
                if (tau_g==tau_e) {
                  # discrete Gamma density with shape k=k+1 & rate = 1/tau_g = 1/tau_e
                  k <- k+1
                  g <- (t^(k-1) * exp(-t/tau_g)) /
                    (gamma(k) * tau_g^k)
                  
                } else {
                # discrete Gamma density (shape k, rate = 1/tau_g)
                g1 <- (t^(k-1) * exp(-t/tau_g)) /
                  (gamma(k) * tau_g^k)
                # discrete Exponential density (rate = 1/tau_e)
                g2 <- 1/tau_e * exp(-t/tau_e)
                
                # linear convolution
                g  <- gsignal::conv(g1, g2, shape = "full")[1:L]
                }
                g
              }
  )
  
  ## ---- shift kernel ------------------------------------------------------
  kernel <- shift_kernel(kernel, perc, method = "linear")

  ## ---- return kernel normalised to unit L1 norm --------------------------
  return(kernel / sum(kernel)) 
}




#' @rdname convkernel
#'
#' @param kernel Numeric vector – the PSF to be shifted.
#' @param perc   Either `"mean"` (shift the centroid) or a
#'               numeric percentile in [0, 1].
#' @param method `"linear"` (sub-sample interpolation) or
#'               `"integer"` (whole-sample shift).
#' @return shifted vector of the same length as `kernel`
shift_kernel <- function(kernel,
                         perc = "mean",
                         method = c("linear", "integer"))
{
  method <- match.arg(method)
  stopifnot(is.numeric(kernel),
            length(kernel) >= 3L,
            (is.numeric(perc) && 0 <= perc && perc <= 1) || is.character(perc))
  
  L   <- length(kernel)
  idx <- seq_len(L)
  
  ## -------- locate the reference point ---------------------------------
  if (is.character(perc)) {
    
    ## centroid (= mean)
    pos <- sum(idx * kernel) / sum(kernel)
    
  } else if (perc == 0) {
    
    ## 0-percentile  -> left boundary
    pos <- 1
    
  } else if (perc == 1) {
    
    ## 100-percentile -> right boundary
    pos <- L
    
  } else {
    
    ## p-quantile inside (0, 1)  ------------------------------------------
    cdf <- cumsum(kernel) / sum(kernel)
    pos <- stats::approx(cdf, idx, xout = perc,
                         ties = "ordered", rule = 2)$y
  }
  
  ## -------- where we want it -------------------------------------------
  target <- if (is.character(perc)) { (L + 1) / 2 } else {  # centre for centroid
    1 + perc * (L - 1)                                      # 1-based percentile
  }
  
  delta  <- target - pos                           # +ve → shift right
  if (abs(delta) < .Machine$double.eps)
    return(kernel / sum(kernel))                             # already aligned
  
  out <- switch(method, linear = {
    f <- stats::approxfun(idx, kernel, yleft = 0, yright = 0)
    f(idx - delta)
  }, integer = {
    s <- round(delta)                            # integer shift
    if (s > 0)
      c(rep(0, s), head(kernel, L - s))
    else if (s < 0)
      c(tail(kernel, L + s), rep(0, -s))
    else
      kernel
  })
  
  out / sum(out)                                   # keep L¹ = 1
}
