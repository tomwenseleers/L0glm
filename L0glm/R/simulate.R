#' @name simul
#' @aliases simul
#' @title Simulate data.
#' 
#' @description
#' This function simulates different types of data based on the method specified. 
#' Various aspects of the data generation can be controlled using the \code{options} parameter.
#'
#' @param method 
#' A numeric value indicating the type of data to simulate. The possible values are:
#' 1: Single channel (single variable) data with Gaussian noise.
#' 2: Single channel data with Poisson noise.
#' 3: Multichannel data with Gaussian noise.
#' 4: Multichannel data with Poisson noise.
#' 5: Orthogonal data.
#' 6: High dimensional data (p > 40000 suggested).
#' @param n 
#' A numeric value indicating the number of observations to generate.
#' @param p 
#' A numeric value indicating the number of predictors (variables) to generate.
#' @param S 
#' Sparsity of true solution (proportion of total coefficients p that are nonzero).
#' @param sparse
#' Boolean indicating if the returned design matrix should be sparse (of Matrix class dsCMatrix) or dense.
#' @param options 
#' A list of optional parameters to control the data generation process. This can include:
#'   For methods 1 to 4: \code{sd_noise} (standard deviation of the noise), 
#'   \code{mean_beta} (mean of the beta coefficients), \code{sd_logbeta} (standard deviation of the log-beta).
#'   For method 5: \code{n_nonzero} (number of nonzero coefficients), \code{beta_sd} (standard deviation of the beta coefficients), 
#'   \code{error_sd} (standard deviation of the error term).
#'   For method 6: \code{sd_noise} (standard deviation of the noise).
#'   If options are \code{NULL}, default values will be used.
#' @param seed
#' Random seed.
#' 
#' @return A variable containing the simulated data.
#'
#' @importFrom statmod rinvgauss
#' @importFrom stats gaussian plogis poisson qlogis rbinom rgamma rpois runif
#' @importFrom methods as is
#' @importFrom utils flush.console
#'
#' @examples
#' \dontrun{
#' sim_data <- simul(method=1, n=500, p=500)
#' }
#' 
#' @export
simul <- function(method = 1, n = 500, p = 500, s = 0.1, sparse = T, options = NULL, seed = 1){
  if(method == 1){
    family = "gaussian"
    multichannel = FALSE
    k <- round(s*n)
    sd_noise <- ifelse(is.null(options), 5, options$sd_noise)
    mean_beta <- ifelse(is.null(options$mean_beta), 2000, options$mean_beta)
  } else if(method == 2){
    family = "poisson"
    k <- round(s*n)
    multichannel =  FALSE
    mean_beta <- ifelse(is.null(options$mean_beta), 2000, options$mean_beta)
  } else if(method == 3){
    family = "gaussian"
    multichannel = TRUE
    k <- round(s*n)
    if(k > 100) k <- 100
    sd_noise = ifelse(is.null(options$sd_noise), 5, options$sd_noise)
    mean_beta <- ifelse(is.null(options$mean_beta), 10000, options$mean_beta)
  } else if(method == 4){
    family = "poisson"
    multichannel = TRUE
    k <- round(s*n)
    if(k > 100) k <- 100
    sd_noise = ifelse(is.null(options$sd_noise), 10, options$sd_noise)
    mean_beta <- ifelse(is.null(options$mean_beta), 10000, options$mean_beta)
  } else if(method == 6){
    if(p < 40000) warning("This method is suggested for data with higher dimension.")
    m <- ifelse(is.null(options), 50, m)
  }
  
  if(method == 1 | method == 3){
    sim <- simulate_spike_train(n = n, p = p, family = family, 
                                seed = seed,
                                k = k, 
                                multichannel = multichannel,
                                mean_beta = mean_beta,
                                sparse = sparse,
                                sd_logbeta = ifelse(is.null(options$sd_logbeta), 1, options$sd_logbeta),
                                sd_noise = sd_noise)
  } else if(method == 2 | method == 4){
    sim <- simulate_spike_train(n = n, p = p, family = family, 
                                seed = seed,
                                k = k, 
                                multichannel = multichannel,
                                mean_beta = mean_beta,
                                sparse = sparse,
                                sd_logbeta = ifelse(is.null(options$sd_logbeta), 1, options$sd_logbeta))
  } else if(method == 5){
    sim <- simulate_standardnormal(n = n, p = p, 
                               n_nonzero = ifelse(is.null(options$n_nonzero), 24, options$n_nonzero),
                               beta_sd = ifelse(is.null(options$beta_sd), 0.5, options$beta_sd), 
                               error_sd = ifelse(is.null(options$error_sd), 1, options$error_sd))
  } else if(method == 6){
    sim <- simulate_ordinis(n = n, p = p, m = m, sd_noise = ifelse(is.null(options$sd_noise), 3, options$sd_noise))
  }
  
  return(sim)
}





#' @name simulate_spike_train
#' @title Simulate a single or multichannel spike train convoluted by a given peak shape function (PSF)
#' 
#' @importFrom Matrix Matrix
#' @importFrom qlcMatrix colMax colMin
#' @importFrom grDevices topo.colors
#' @importFrom graphics abline matplot points
#'
#' @description
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
#' @param thresh XXX.
#' @param sparse XXX.
#' @param multichannel XXX.
#' @param seed The seed used to sample the nonzero coefficients and to generate
#'   the noisy response signal.
#' @param Plot Should the generated data be plotted?
#' 
#' @return
#' Returns a list with the following elements:
#' \item{\code{X}}{a \code{n} by \code{p} covariate matrix. The matrix is a banded
#' matrix with time-shifted versions of the specified peak shape function
#' (also known as the Point Spread Function or PSF).}
#' \item{\code{beta_true}}{a vector of length \code{p} containing the true
#' coefficients. This is the spike train used to generate the data, so that
#' \code{y = beta_true \%*\% X + e}.}
#' \item{\code{y}}{a vector of length \code{n} containing the noisy response signal
#' given by \code{y = beta_true \%*\% X + e}.}
#' \item{\code{y_true}}{a vector of length \code{n} containing the response signal
#' given by \code{y_true = beta_true \%*\% X}.}
#' \item{\code{x}}{a vector of length \code{n} containing the time points.}
#' \item{\code{x_beta}}{a vector of length \code{p} containing the locations of the
#' simulated spikes. Equal to \code{x} when \code{p=n}.}
#'
#' @export

simulate_spike_train <- function(n = 100, p = n, k = 10,
                                  mean_beta = 10000, sd_logbeta = 1,
                                  # TO DO: add intercept argument
                                  family = c("poisson","gaussian"),
                                  sd_noise = 10,
                                  psf = function(x, u, w=2.5)
                                    exp(((x-u)^2)/(-2*(w^2))),
                                  thresh = 1E-5,
                                  sparse = TRUE,
                                  multichannel = FALSE,
                                  seed = 123, 
                                  Plot = FALSE
){

  # TODO: update documentation above
  # TODO: add intercept if desired 
  
  set.seed(seed)
  family = match.arg(family)
  x = 1:n
  x_beta = seq(1, n, length.out = p)
  # simulated peak locations:
  u = sample(1:p, k, replace = FALSE) 
  
  if (multichannel) { 
    data("simulated_spectra", package = "L0glm", envir = environment()) # data with 100 simulated example GC/MS mass spectra
    simulated_spectra = sweep(simulated_spectra, 1, STATS=rowSums(simulated_spectra), FUN="/") # L1 normalize
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
  if (sparse) X = as(Matrix(X, sparse=T),"generalMatrix")
  
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
           "Simulated blurred superimposed spike train\n(red = true coefficients, black = observed signal) with ",
           family, " noise"),
         xlab = "Time", ylab = "Signal")
    points(x, y, pch = 16, cex = 0.5)
    lines(x_beta, beta_true, type = "h", col = "red") } else {
      # plot for multichannel simulation
      par(mfrow = c(1,1))
      image(x=1:n, y=1:m, z=as.matrix(y)^0.2, 
            col = topo.colors(255), useRaster=TRUE,
            main = paste0(
              "Simulated blurred multichannel spike train\n(red = true support) with ",
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
    x_beta = x_beta,
    family = family
  ))
}



#' Simulate synthetic GC/MS chromatograms with known ground truth
#'
#' `simulate_GCMS()` generates an \eqn{n \times m}{n x m} data matrix \eqn{D}
#' that mimics a GC/MS chromatogram with \eqn{k} compounds whose peaks have
#' Gaussian shape, but slightly varying widths, and whose electron-impact (EI)
#' spectra are drawn from a reference library.  Each compound elutes once as a
#' Gaussian peak with a random centre, width and total-ion height.  Compound
#' peaks are returned in *time-sorted* order.  Either Poisson (counting) noise
#' or additive Gaussian noise can be added.
#'
#' @details
#' \strong{Spectral library}\cr
#' Let \eqn{S} have one compound per \emph{row} and mass channels in
#' \emph{columns}; rows are renormalised to unit sum.\cr\cr
#'
#' \strong{Elution profile}\cr
#' For compound \eqn{j} the underlying profile is
#' \deqn{%
#'   W_{ij} = \beta_j
#'            \exp\!\Bigl[-\tfrac12\bigl((t_i-c_j)/\sigma_j\bigr)^2\Bigr],
#'            \qquad i = 1,\dots,n%
#' }{%
#'   W_ij = beta_j * exp(-0.5 * ((t_i - c_j)/sigma_j)^2),
#'   i = 1,…,n%
#' }
#' where
#' \eqn{c_j \sim \mathrm{Unif}(1,n),                 }{c[j] ~ Unif(1,n),}
#' \eqn{\sigma_j \sim \mathrm{Unif}(a,b),             }{sigma[j] ~ Unif(a,b),}
#' \eqn{\beta_j  \sim \mathrm{LogNormal}(\mu,\sigma).}{beta[j]  ~ LogNormal(mu,sigma).}
#' Here \code{c(a,b) = peak_width_range},
#' \code{mu = mean_beta} and \code{sigma = sd_logbeta}.
#'\cr\cr
#'
#' \strong{Signal and noise}\cr
#' The noiseless signal is
#' \deqn{ D_{\mathrm{true}} = W_{\mathrm{true}}\,H_{\mathrm{true}} }
#' If \code{family = "poisson"} each entry is replaced by a Poisson draw with
#' that mean; otherwise Gaussian noise with standard deviation
#' \code{sd_noise} is added, and negative values are truncated to 0.
#'\cr\cr
#'
#' \strong{Returned order}\cr
#' Peaks are returned in increasing order of their peak centre \eqn{c_j}.
#'\cr\cr
#'
#' \strong{Background}\cr
#' Any background is assumed to have been subtracted in advance (e.g. by a
#' quantile spline under-approximation).
#'
#' @param n  Integer. Number of time points (rows) in the simulated chromatogram.
#' @param k  Integer. True number of compounds/peaks present (true rank).
#' @param S  Matrix whose rows are EI spectra.  If \code{NULL}, the example data
#'           set \code{data(simulated_spectra)} is used.
#' @param peak_width_range Numeric length-2. Uniform range for Gaussian peak
#'           widths (\eqn{\sigma}) in time points.
#' @param mean_beta  Mean of the total-ion peak heights (original scale).
#' @param sd_logbeta Standard deviation on the \emph{log} scale for the
#'           log-normal peak heights.
#' @param family  \code{"poisson"} (default) or \code{"gaussian"}.
#' @param sd_noise s.d. of additive Gaussian noise (ignored for Poisson).
#' @param seed     Integer seed for reproducibility.
#'
#' @return A list with components
#' \describe{
#'   \item{\code{D}}{Observed \eqn{n \times m}{n x m} chromatogram (dense).}
#'   \item{\code{W_true}}{\eqn{n \times k}{n x k} true elution profiles (sparse).}
#'   \item{\code{H_true}}{\eqn{k \times m}{k x m} true spectra (sparse).}
#'   \item{\code{k}}{True number of compounds.}
#'   \item{\code{centres}, \code{widths}, \code{heights}, \code{peak_width_range}}{Peak parameters.}
#'   \item{\code{family}}{Noise model used.}
#'   \item{\code{sd_noise}}{s.d. of noise if Gaussian.}
#'   \item{\code{compound_idx}}{Row indices of \code{S} that were selected.}
#' }
#'
#' @examples
#' data(simulated_spectra)               # 100 × 460 reference library
#' library(Matrix)
#' sim <- simulate_GCMS(
#'          n = 200, k = 50,
#'          S = simulated_spectra,
#'          family = "poisson")
#' image(as.matrix(sim$D))               # chromatogram with noise
#' matplot(sim$W_true, type = "l")       # true elution profiles
#' image(sim$H_true)                     # EI spectra
#'
#' @import Matrix
#' @export
simulate_GCMS <- function(n = 200,
                          k = 50,
                          S = NULL,              # spectral library (n × m)
                          peak_width_range = c(1, 5),
                          mean_beta   = 1e4,
                          sd_logbeta  = 1,
                          family      = c("poisson", "gaussian"),
                          sd_noise    = if (family=="gaussian") { 1 } else { NULL },
                          seed        = 123) {
  
  if (is.null(S)) {
    data(simulated_spectra, envir = environment())
    S <- simulated_spectra
  } else { S <- Matrix(S, sparse = TRUE) }

  family <- match.arg(family)
  set.seed(seed)
  
  ## -------- helper ------------------------------------------------------
  psf_fun <- function(x, centre, width)
    exp(-0.5 * ((x - centre) / width)^2)
  
  ## -------- pick k spectra ---------------------------------------------
  if (k > nrow(S)) stop("k exceeds number of available spectra")
  idx <- sample(nrow(S), k)
  H_true   <- Matrix(S[idx, , drop = FALSE], sparse = TRUE)
  H_true   <- sweep(H_true, 1, rowSums(H_true), "/")             # rows made to sum to 1
  
  ## -------- elution matrix W -------------------------------------------
  tvec    <- seq_len(n)
  centres <- sort(runif(k, 1, n)) # sorted peak centres
  widths  <- runif(k, peak_width_range[1], peak_width_range[2])
  heights <- rlnorm(k, meanlog = log(mean_beta) - 0.5 * sd_logbeta^2,
                    sdlog   = sd_logbeta)
  
  W_true <- Matrix(sapply(seq_len(k), function(j)
    heights[j] * psf_fun(tvec, centres[j], widths[j])), sparse = TRUE)
  dimnames(W_true) <- list(NULL, paste0("peak", idx))
  
  ## -------- noiseless signal -------------------------------------------
  D_true <- W_true %*% H_true
  
  ## -------- add noise ---------------------------------------------------
  if (family == "poisson") {
    D <- matrix(rpois(length(D_true), lambda = as.vector(D_true)), nrow = n)
  } else {
    D <- D_true + rnorm(length(D_true), sd = sd_noise)
    D[D < 0] <- 0
  }
  
  list(D = as.matrix(D),
       D_true = as.matrix(D_true),
       W_true = W_true,
       H_true = H_true,
       k = k,
       centres = centres,
       widths  = widths,
       heights = heights,
       peak_width_range = peak_width_range,
       family = family,
       sd_noise = sd_noise,
       compound_idx = idx)
}



# simulate single channel data with L2 norm renormalized standard normal variables & gaussian noise with SD of 1 ####
# cf. Frommlet & Nuel 2016 p. 9

simulate_standardnormal <- function(n=100, p=100, n_nonzero=24, beta_sd=0.5, error_sd=1, 
                                    seed.problem=1, seed.noise=sample.int(1E4, 1, replace=T)) {
  # Generate X matrix
  set.seed(seed.problem)
  X <- matrix(rnorm(n*p), n, p)
  
  # # we L2 norm normalize the columns of X first
  # X.n <- sqrt(Matrix::colSums(X^2)) 
  # X <- as.matrix(X %*% Matrix::Diagonal(x = 1 / X.n))
  
  # Generate beta vector with all zeros initially
  beta <- rep(0, p)
  
  # Randomly select n_nonzero indices to assign non-zero values
  non_zero_indices <- sample.int(n=p, size=n_nonzero)
  
  # Generate beta vector
  beta[non_zero_indices] <- rnorm(n_nonzero, 0, beta_sd)
  
  # Generate errors
  set.seed(seed.noise)
  epsilon <- rnorm(n, 0, error_sd)
  
  # Generate y vector
  y <- X %*% beta + epsilon
  
  return(list(X=X, y=y, beta_true=beta))
}




# simulate some univariate data with sparse solution for any GLM family ####
simulate_glm <- function(n = 10000, 
                         p = 100, 
                         n_nonzero = 50, 
                         intercept = 1, 
                         family = "gaussian",
                         error_sd = 1,
                         ntot = 10) {
  # n: number of observations
  # p: number of variables (excluding intercept)
  # n_nonzero: nr of variables (excluding intercept term) with nonzero coefficients
  # intercept: coefficient for intercept term
  # family: response distribution ("gaussian", "poisson", "binomial", "gamma", "inverse.gaussian")
  # error_sd: sd of error for Gaussian case
  # ntot: total sample size for the binomial case
  
  # Create matrix of predictor variables
  X <- cbind(1, # for intercept 
             matrix(runif(n * p), nrow = n, ncol = p))
  
  # Set true coefficients 
  beta <- c(intercept, rep(0, p))
  # if (family=="gamma") { beta[sample(2:(p+1), n_nonzero)] <- abs(rnorm(n_nonzero)) } else { 
  # sample from half normal distribution for gamma to make sure that true coefficients are always positive 
  beta[sample(2:(p+1), n_nonzero)] <- rnorm(n_nonzero) 
  # }
  
  # Set link function and inverse link function based on family
  if (family == "gaussian") {
    linkinv <- identity
    linkf <- identity
  } else if (family == "poisson") {
    linkinv <- exp
    linkf <- log
  } else if (family == "binomial") {
    linkinv <- plogis
    linkf <- qlogis
  } else if (family == "gamma") { # log link Gamma
    linkinv <- exp
    linkf <- log
  } else if (family == "inverse.gaussian") { # log link inverse Gaussian
    linkinv <- exp # function(eta) 1 / eta^2
    linkf <- log # function(mu) sqrt(1 / mu)
  } else {
    stop("Invalid family")
  }
  
  # Generate response variable
  eta <- X %*% beta
  mu <- linkinv(eta)
  if (family == "gaussian") {
    y <- rnorm(n, mean = mu, sd = error_sd)
  } else if (family == "poisson") {
    y <- rpois(n, lambda = mu)
  } else if (family == "binomial") {
    y <- rbinom(n, size = ntot, prob = mu)/ntot
  } else if (family == "gamma") {
    # shape = 1 and scale = 1/rate. scale here is mu & rate 1/mu to get the mean mu
    y <-  rgamma(n, shape = 1, scale = mu)
  } else if (family == "inverse.gaussian") {
    y <- statmod::rinvgauss(n, mean = mu, shape = 1)
  }
  
  # Set prior weights based on family
  if (family == "binomial") {
    weights <- rep(ntot, n)
  } else {
    weights <- rep(1, n)
  }
  
  # L2 norm normalize columns of X
  # X.n <- sqrt(Matrix::colSums(X^2)) # Calculate L2 norm of columns
  # X <- as.matrix(X %*% Matrix::Diagonal(x = 1 / X.n))
  # beta <- beta * X.n
  
  list(X = X, y = y, weights = weights, beta_true = beta, eta_true = eta, mu_true = mu)
}

# data = simulate_glm(n = 100, p = 1000000, n_nonzero = 50, intercept = 1, family = "gaussian")


simulate_ordinis <- function(n = 500, p = 50000, m = 50, sd_noise = 3){
  beta_true <- matrix(c(runif(m), rep(0, p - m))) 
  X <- matrix(rnorm(n * p, sd = sd_noise), n, p)
  y <- drop(X %*% beta_true) + rnorm(n)
  
  return(list(X = X, y = y, beta_true = beta_true))
}



