#' Simulate example data (blurred spike train)
#' 
#' @importFrom Matrix Matrix
#' @importFrom qlcMatrix colMax colMin
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
#' 
simulate_spike_train <- function(n = 100, p = n, k = 10,
                                  mean_beta = 10000, sd_logbeta = 1,
                                  family = c("poisson","gaussian"),
                                  sd_noise = 10,
                                  psf = function(x, u, w=2.5)
                                    exp(((x-u)^2)/(-2*(w^2))),
                                  thresh = 1E-5,
                                  sparse = TRUE,
                                  multichannel = FALSE,
                                  seed = 123, Plot = TRUE
){
  require(qlcMatrix)
  # TODO: update documentation above
  # TODO: add intercept if desired 
  
  set.seed(seed)
  family = match.arg(family)
  x = 1:n
  x_beta = seq(1, n, length.out = p)
  # simulated peak locations:
  u = sample(2:p, k, replace = FALSE) # changed so does not take intercept column
  
  if (multichannel) { 
    data(simulated_spectra) # data with 100 simulated example GC/MS mass spectra
    # simulated_spectra = sweep(simulated_spectra, 1, STATS=rowSums(simulated_spectra), FUN="/") # L1 normalize
    if (k>nrow(simulated_spectra)) stop("Maximum k is 100 when multichannel=TRUE")
    m = ncol(simulated_spectra) # 460 mass channels
    S = simulated_spectra[sample(1:nrow(simulated_spectra), k),] # Randomly select k spectra
  } else { m = 1
           S = matrix(1, nrow=k, ncol=m) }
  
  beta_true = Matrix(0, nrow=p, ncol=m, sparse=sparse) # rep(0, p)
  # simulated peak heights (sum over all channels) assuming log-normal distribution :
  peakheights = exp(rnorm(k, log(mean_beta), sd_logbeta))
  beta_true[u, 1:m] = .sparseDiagonal(k, x=peakheights) %*% S # scale spectra by peak height
  
  # remove invariant channels
  invariant = (qlcMatrix::colMax(beta_true)==qlcMatrix::colMin(beta_true)) # zero across all spectra, these channels are removed
  m = m-sum(invariant)
  beta_true = beta_true[,!invariant,drop=F]
  S = S[,!invariant,drop=F]
  
  # covariate matrix X = shifted copies of PSF
  # TODO maybe calculate this as a linear convolution, 
  # cf https://stackoverflow.com/questions/66373113/linear-convolution-in-r-using-bandsparse-matrix-multiplication:
  X = do.call(cbind, lapply(x_beta, function (u) psf(x, u))) 
  X[X < thresh] = 0
  if (sparse) X = as(X,"dgCMatrix")
  
  y_true = Matrix(X %*% beta_true, sparse=sparse)
  beta_true = Matrix(beta_true, sparse=sparse)
  
  if (family == "poisson") { y = apply(y_true, 2, function (col) rpois(n, col)) } else {
    y = y_true + rnorm(prod(dim(y_true)), mean = 0, sd = sd_noise)
    y[y < 0] = 0
  }
  
  # remove invariant channels
  invariant = (qlcMatrix::colMax(y)==qlcMatrix::colMin(y)) # invariant channels are removed
  m = m-sum(invariant)
  beta_true = beta_true[,!invariant,drop=F]
  S = S[,!invariant,drop=F]
  y_true = y_true[,!invariant,drop=F]
  y = y[,!invariant,drop=F]

  # X[,1]=rep(1,n) # intercept column # TODO add option for intercept

  if (Plot) {
    if (!multichannel) {
    # plot for single channel simulation
    par(mfrow = c(1, 1))
    plot(x, y, type = "l",
         main = paste0(
           "Simulated blurred superimposed spike train\n(red=true coefficients, black=observed signal) with ",
           family, " noise"),
         xlab = "Time", ylab = "Signal")
    points(x, y, pch = 16, cex = 0.5)
    lines(x_beta, beta_true, type = "h", col = "red") } else {
      # plot for multichannel simulation
      par(mfrow = c(1,1))
      image(x=1:n, y=1:m, z=as.matrix(y)^0.2, 
            col = topo.colors(255), useRaster=TRUE,
            main = paste0(
              "Simulated blurred multichannel spike train\n(red=true support) with ",
              family, " noise"), 
            xlab="Time", ylab="Channel")
      abline(v=u*n/p, col="red")
    }
  }

  colnames(X) = paste0("x", 1:ncol(X))
  colnames(y) = paste0("y", 1:ncol(y))
  colnames(y_true) = paste0("y", 1:ncol(y_true))
  
  return(list(
    X = X,
    beta_true = beta_true,
    y = as.matrix(y),
    y_true = as.matrix(y_true),
    x = x,
    x_beta = x_beta
  ))
}
