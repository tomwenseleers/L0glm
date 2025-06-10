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



#' Simulate synthetic GC-MS chromatograms with ground-truth signal + background
#'
#' `simulate_GCMS()` draws \eqn{k} compounds from an EI reference library `S`
#' (rows = compounds, columns = m/z values, e.g. 33–550) and generates an \eqn{n \times m}{n
#' × m} chromatogram `D` whose signal is the sum of
#' k Gaussian peaks (`W_true H_true`) and
#' a smooth, monotonically increasing B-spline background
#' (`W_background*bg_scale H_background`).
#' 
#' @details
#' \itemize{
#'   \item A \emph{compound-specific base width} \eqn{\sigma_{0}} is drawn
#'         once when the compound is first selected from the library.
#'   \item In subsequent simulated chromatograms the same compound retains
#'         the same \eqn{\sigma_{0}}, but the \emph{effective} width becomes
#'         \deqn{\sigma \;=\; \sigma_{0}\Bigl(1 +
#'               \texttt{width\_height\_coef}\,h_{\mathrm{norm}}\Bigr),}
#'         so it can still broaden or narrow in proportion to its
#'         (possibly resampled) peak height.
#' }
#'
#' Gaussian peak width and RT can be coupled to peak height; an additional
#' *chromatographic drift* proportional to the injection number `sample_nr` is
#' also included.  Independent random seeds allow you to keep the same compounds
#' while resampling peak heights (or vice-versa).
#'
#' @section Reference library: `S` must contain
#' - column names (character) for the m/z channels (e.g. `"33":"550"`), and
#' - attribute `attr(S, "RT")` for the retention times *(numeric, unit: mins)*.
#'
#' @param n  Integer ≥ 1. Number of time points (rows).
#' @param k  Integer ≥ 1. Number of compounds/peaks.
#' @param S  Spectral library.  If `NULL` the example data
#'   `data(simulated_spectra)` (2000 × 518) is loaded.
#' @param peak_width_range Length-2 numeric. Uniform base range for Gaussian peak width σ
#'   (*before* any height-coupling is applied and on time point scale).
#' @param mean_beta,sd_logbeta Parameters of the log-normal height distribution.
#' @param width_height_coef  Numeric ≥ 0.  Multiplicative widening:
#'   \deqn{\sigma_j \leftarrow \sigma_j\bigl[1+\texttt{coef}\,h_j\bigr]}{} where
#'   \eqn{h_j=(\beta_j-\min\beta)/(\max\beta-\min\beta)}.
#' @param rt_height_coef Numeric.  Positive values push high peaks to slightly
#'   larger RTs: \deqn{c_j \leftarrow c_j
#'   \bigl[1+\texttt{rt\_height\_coef}\,h_j\bigr]}{}
#' @param sample_nr Integer ≥ 1. Injection/order number.
#' @param rt_drift_coef Numeric.  Additive drift in *time points* per injection:
#'           \eqn{c_j \leftarrow c_j + \texttt{rt\_drift\_coef}\,
#'           (\texttt{sample\_nr}-1)}.
#' @param bg_scale  Numeric ≥ 0. Multiplicative factor applied to the whole
#'   B-spline background.
#' @param family  `"poisson"` (default) or `"gaussian"`.
#' @param sd_noise s.d. of additive Gaussian noise (`family="gaussian"`).
#' @param seed_compounds Integer seed controlling which library spectra are
#'   drawn.
#' @param seed_peaks Integer seed for all other random quantities (heights,
#'   widths, noise, …).
#'
#' @return A list with matrices:
#' \describe{
#'   \item{D}{Observed chromatogram (signal + background + noise) (dense matrix).}
#'   \item{D_true}{Noiseless signal `W_true H_true + W_background H_background` (dense matrix).}
#'   \item{W_true, H_true}{Elution profiles and EI spectra of the *k* peaks (sparse matrices).}
#'   \item{W_background, H_background}{Smooth B-spline background bases and spectra (sparse matrices).}
#'   \item{k_background}{Number of B-spline basis columns (`max(2, n/1000)`).}
#'   \item{centres, widths, heights}{Peak parameters after all applied
#'         height/RT/drift adjustments.}
#'   \item{widths_base}{Base compound width before intensity-based broadening.}
#'   \item{compound_idx}{Row indices of `S` that were selected.}
#'   \item{sample_nr, rt_drift_coef, width_height_coef, rt_height_coef}{Echoed.}
#' }
#'
#' @examples
#' data(simulated_spectra)
#' sim <- simulate_GCMS(n = 10000,
#'                      k = 500,
#'                      S = simulated_spectra,
#'                      peak_width_range = c(3, 5),
#'                      mean_beta  = 1e4,
#'                      sd_logbeta = 1,
#'                      width_height_coef = 3,
#'                      rt_height_coef    = 0.02,
#'                      sample_nr = 1,
#'                      rt_drift_coef = 0.001,
#'                      bg_scale = 100,
#'                      family = "poisson")
#' cols <- colorRampPalette(c("black", "black","blue","green","yellow","orange", "red"))(1000) 
#' Dvis <- sim$D[seq(1, 10000, by = 10), , drop = FALSE]^0.12
#' image(x = as.numeric(rownames(Dvis)), 
#'       y = as.numeric(colnames(Dvis)),
#'       z = Dvis,
#'       col    = cols,
#'       xlab   = "Time (min)",
#'       ylab   = "Mass (m/z)",
#'       main   = "Simulated GC/MS chromatogram")
#'
#' @importFrom Matrix Matrix
#' @importFrom splines bs
#' @export
simulate_GCMS <- function(
    n  = 10000,
    k  = 500,
    S  = NULL,
    peak_width_range = c(3, 5),
    mean_beta  = 1e4,
    sd_logbeta = 1,
    width_height_coef = 3,
    rt_height_coef    = 0.02,
    sample_nr         = 1,
    rt_drift_coef     = 0.001,
    bg_scale          = mean_beta/10,
    family   = c("poisson", "gaussian"),
    sd_noise = if (match.arg(family) == "gaussian") 1 else NULL,
    seed_compounds = 123,
    seed_peaks     = 456
) {
  
  ## -------- load / check spectral library --------------------------------
  if (is.null(S)) {
    data(simulated_spectra, envir = environment())
    S <- simulated_spectra
  }
  if (is.null(attr(S, "RT")))
    stop("Attribute 'RT' is missing from spectral library.")
  if (is.null(colnames(S)))
    stop("Library must have m/z column names (e.g. \"33\":\"550\").")
  
  S <- Matrix(S, sparse = TRUE)              # cheap if already sparse
  m  <- ncol(S)
  mz_names <- colnames(S)
  
  ## -------- pick spectra --------------------------------------------------
  set.seed(seed_compounds)
  if (k > nrow(S))
    stop("k exceeds number of available spectra in library.")
  idx <- sample(nrow(S), k)
  H_true <- Matrix(S[idx, , drop = FALSE], sparse = TRUE)
  H_true <- sweep(H_true, 1, rowSums(H_true), "/")     # rows sum to 1
  rownames(H_true) <- paste0("cmp", idx)
  colnames(H_true) <- mz_names        # keep m/z names
  
  ## -------- compound-specific base widths ---------------------------------
  ## use the *same* RNG stream that selected idx so base widths are tied to compounds
  widths_base <- runif(k, peak_width_range[1], peak_width_range[2])
  
  ## -------- sample-specific quantities ------------------------------------
  set.seed(seed_peaks)
  
  rt_lib   <- attr(S, "RT")[idx]             # minutes
  rt_min   <- min(attr(S, "RT"))
  rt_max   <- max(attr(S, "RT"))
  
  # map minutes → scan indices 1..n
  centres_base <- 1 + (rt_lib - rt_min) / (rt_max - rt_min) * (n - 1)
  
  # log-normal heights & normalised heights
  heights <- rlnorm(
    k,
    meanlog = log(mean_beta) - 0.5 * sd_logbeta^2,
    sdlog   = sd_logbeta
  )
  h_norm <- (heights - min(heights)) /
    (max(heights) - min(heights) + .Machine$double.eps)
  
  # width-height coupling
  widths  <- widths_base * (1 + width_height_coef * h_norm)
  
  # final centres with height shift and chromatographic drift
  centres <- centres_base * (1 + rt_height_coef * h_norm) +
    rt_drift_coef * (sample_nr - 1)
  
  # keep in bounds & sort
  centres <- pmin(pmax(centres, 1), n)
  ord     <- order(centres)
  centres <- centres[ord]
  widths <- widths[ord]
  widths_base  <- widths_base[ord]
  heights <- heights[ord]
  idx     <- idx[ord]
  H_true <- H_true[ord, , drop = FALSE]
  
  ## -------- elution matrix W_true ----------------------------------------
  tvec <- seq_len(n)
  psf  <- function(x, c, s) exp(-0.5 * ((x - c) / s)^2)
  
  W_true <- Matrix(
    vapply(seq_len(k), function(j) {
      heights[j] * psf(tvec, centres[j], widths[j])
    }, numeric(n)),
    sparse = TRUE
  )
  dimnames(W_true) <- list(NULL, paste0("peak", idx))
  
  ## -------- B-spline background ------------------------------------------
  k_background <- max(2L, ceiling(n / 1000))
  B            <- splines::bs(
    x         = seq(0, 1, length.out = n),
    df        = k_background,
    degree    = 3,
    intercept = TRUE
  )
  ## force a monotonically increasing envelope
  env          <- seq(0.2, 1.0, length.out = n)
  W_background <- Matrix(B * env, sparse = TRUE) * bg_scale
  colnames(W_background) <- paste0("bg", seq_len(k_background))
  
  ## single background spectrum = noisy mean of the chosen spectra
  spec_bg <- drop(colMeans(S[idx, , drop = FALSE]))
  spec_bg <- pmax(spec_bg + rnorm(m, sd = 0.01 * max(spec_bg)), 0)
  spec_bg <- spec_bg / sum(spec_bg)          # L1-normalise
  
  H_background <- Matrix(
    matrix(rep(spec_bg, k_background),
           nrow = k_background, byrow = TRUE),
    sparse = TRUE)
  colnames(H_background) <- mz_names
  rownames(H_background) <- paste0("bg", seq_len(k_background))
  
  ## -------- noiseless signals --------------------------------------------
  D_true <- W_true %*% H_true + W_background %*% H_background
  rownames(D_true) <- max(rt_lib)*(1:nrow(D_true))/10000
  
  ## -------- add noise -----------------------------------------------------
  family <- match.arg(family)
  if (family == "poisson") {
    D <- matrix(
      rpois(length(D_true), lambda = as.vector(D_true)),
      nrow = n, dimnames = dimnames(D_true)
    )
  } else {
    D <- D_true + rnorm(length(D_true), sd = sd_noise)
    D[D < 0] <- 0
    dimnames(D) <- dimnames(D_true)
  }
  
  ## -------- return --------------------------------------------------------
  list(
    D             = as.matrix(D),
    D_true        = as.matrix(D_true),
    W_true        = W_true,
    H_true        = H_true,
    W_background  = W_background,
    H_background  = H_background,
    k             = k,
    k_background  = k_background,
    centres       = centres,
    widths        = widths,
    widths_base   = widths_base,
    heights       = heights,
    peak_width_range = peak_width_range,
    family        = family,
    sd_noise      = sd_noise,
    compound_idx  = idx,
    sample_nr     = sample_nr,
    width_height_coef = width_height_coef,
    rt_height_coef    = rt_height_coef,
    rt_drift_coef     = rt_drift_coef
  )
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



