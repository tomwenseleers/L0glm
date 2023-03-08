#' Simulate example data (blurred spike train)
#'
#' @description
#'
#' This function generates synthetic data consisting of a time series made of a
#' spike train convoluted by a Gaussian blur kernel ("point spread function")
#' with Poisson or Gaussian noise added. Data of this form are
#' commonly encountered in the context of various signal deconvolution problems,
#' where one would like to deconvolute an observed blurry and noisy signal with
#' a known blur kernel (psf) assuming a given noise model. For example, in Gas
#' chromatography/Mass spectrometry one would be interested in deconvoluting a
#' chromatogram, which is a noisy version of an unknown number of superimposed
#' blurred peaks, by a given measured peak shape function / blur kernel.
#' The objective is to recover the underlying spike train
#' (chromatogram peaks in the GC/MS example) used to generate the data by
#' regressing the observed noisy signal onto a banded matrix consisting of
#' shifted copies of the Gaussian kernel used to blur the true signal. The
#' simulated error can be specified to be either Poisson or Gaussian.
#'
#' The function is used as a test case in the \code{\link{L0glm}} examples.
#'
#' @param n Number of observations (time points of the response signal).
#' @param p Number of features (shifted copies of the blur kernel). Defaults to
#'   \code{n}. A situation with \code{p > n} implies that superresolution
#'   deconvolution will be used.
#' @param k Number of non-zero coefficients. Smaller values imply greater
#'   sparsity in the coefficients. This value should be smaller than \code{p}.
#' @param mean_beta The non-zero coefficients are drawn from a log-normal
#'   distribution with mean \code{mean_beta} and standard deviation (on a log
#'   scale) of \code{sd_logbeta}.
#' @param sd_logbeta The non-zero coefficients are drawn from a log-normal
#'   distribution with mean \code{mean_beta} and standard deviation (on a log
#'   scale) of \code{sd_logbeta}.
#' @param family Noise model - currently \code{"poisson"} or \code{"gaussian"}
#'   are supported, with \code{"poisson"} default.
#' @param sd_noise Standard deviation of the noise when
#'   \code{family="gaussian"}.
#' @param psf Point spread function to use to blur spike signal with. Should be
#'   a function with as first two expected arguments \code{x} (time) and
#'   \code{u} (the mode of the peak). Defaults to a Gaussian psf
#'   \code{function(x, u, w=2.5) exp(((x-u)^2)/(-2*(w^2)))}.
#' @param seed The seed used to sample the nonzero coefficients and to generate
#'   the noisy response signal.
#' @param Plot Should the generated data be plotted?
#'
#' @return
#' The function returns a list with the following elements:

#' \item{X}{: a \code{n} by \code{p} covariate matrix. The matrix is a banded
#' matrix with shifted versions of the specified point spread function for
#' which a subset (where \code{beta_true>0}) has associated non-zero
#' (true) coefficients.}
#' \item{beta_true}{: a vector of length \code{p} containing the true
#' coefficients.
#' This is the spike train used to generate the data, so that
#' \code{y = beta_true %*% X + e}.}
#' \item{y}{: a vector of length \code{n} containing the noisy response signal
#' given by \code{y = beta_true %*% X + e}.}
#' \item{y_true}{: a vector of length \code{n} containing the response signal
#' given by \code{y_true = beta_true %*% X}.}
#' \item{x}{: a vector of length \code{n} containing the time points.}
#' \item{x_beta}{: a vector of length \code{p} containing the locations of the
#' simulated spikes. Equal to \code{x} when \code{p=n}.}
#'
#' @export
simulate_spike_train <- function(n = 100, p = n, k = 10,
                                  mean_beta = 1000, sd_logbeta = 1,
                                  family = c("poisson","gaussian"),
                                  sd_noise = 1,
                                  psf = function(x, u, w=2.5)
                                    exp(((x-u)^2)/(-2*(w^2))),
                                  seed = 123, Plot = TRUE
){
  set.seed(seed)
  family = match.arg(family)
  x = 1:n
  x_beta = seq(1, n, length.out = p)
  # simulated peak locations:
  u = sample(2:p, k, replace = FALSE) #changed so does not take intercept column
  beta_true = rep(0, p)
  # simulated peak heights assuming log-normal distribution:
  beta_true[u] = exp(rnorm(k, log(mean_beta), sd_logbeta))
  X = do.call(cbind, lapply(x_beta, function (u) psf(x, u)))
  y_true = as.vector(X %*% beta_true)
  if (family == "poisson") { y = rpois(n, y_true) } else {
    y = y_true + rnorm(n, mean = 0, sd = sd_noise)
    y[y < 0] = 0
  }

  X[,1]=rep(1,n)

  if (Plot) {
    par(mfrow = c(1, 1))
    plot(x, y, type = "l",
         main = paste0(
           "Simulated blurred superimposed spike train\n(red=true coefficients, black=observed signal) with ",
           family, " noise"),
         xlab = "Time", ylab = "Signal")
    points(x, y, pch = 16, cex = 0.5)
    lines(x_beta, beta_true, type = "h", col = "red")
  }

  return(list(
    X = X,
    beta_true = beta_true,
    y = y,
    y_true = y_true,
    x = x,
    x_beta = x_beta
  ))
}
