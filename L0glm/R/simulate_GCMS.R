#' Simulate synthetic GC-MS chromatograms with ground-truth signal + background
#'
#' `simulate_GCMS()` draws \eqn{k} compounds from an EI reference library `S`
#' (rows = compounds, columns = m/z values, e.g. 41-298) and generates an \eqn{n \times m}{n
#' x m} chromatogram `D` whose signal is the sum of
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
#' - column names (character) for the m/z channels (e.g. `"41":"298"`), and
#' - attribute `attr(S, "RT")` for the retention times *(numeric, unit: mins)*.
#'
#' @param n  Integer >= 1. Number of time points (rows).
#' @param k  Integer >= 1. Number of compounds/peaks.
#' @param S  Spectral library.  If `NULL` the example data
#'   `data(simulated_spectra)` (2000 spectra x 258 m/z values) is loaded.
#' @param peak_width_range Length-2 numeric. Uniform base range for Gaussian peak width sigma
#'   (*before* any height-coupling is applied and on time point scale).
#' @param mean_beta,sd_logbeta Parameters of the log-normal height distribution.
#' @param width_height_coef  Numeric >= 0.  Multiplicative widening:
#'   \deqn{\sigma_j \leftarrow \sigma_j\bigl[1+\texttt{coef}\,h_j\bigr]}{} where
#'   \eqn{h_j=(\beta_j-\min\beta)/(\max\beta-\min\beta)}.
#' @param rt_height_coef Numeric.  Positive values push high peaks to slightly
#'   larger RTs: \deqn{c_j \leftarrow c_j
#'   \bigl[1+\texttt{rt\_height\_coef}\,h_j\bigr]}{}
#' @param sample_nr Integer >= 1. Injection/order number.
#' @param rt_drift_coef Numeric.  Additive drift in *time points* per injection:
#'           \eqn{c_j \leftarrow c_j + \texttt{rt\_drift\_coef}\,
#'           (\texttt{sample\_nr}-1)}.
#' @param bg_scale  Numeric >= 0. Multiplicative factor applied to the whole
#'   B-spline background.
#' @param family  `"poisson"` (default) or `"gaussian"`.
#' @param sd_noise s.d. of additive Gaussian noise (`family="gaussian"`).
#' @param seed_compounds Integer seed controlling which library spectra are
#'   drawn.
#' @param seed_peaks Integer seed for all other random quantities (heights,
#'   widths, noise, ...).
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
#' @example examples/simulate_GCMS.R
#' @importFrom Matrix Matrix
#' @importFrom splines bs
#' @export
#' @author
#' Tom Wenseleers \email{tom.wenseleers@kuleuven.be} \\
#' Laboratory of Socioecology and Social Evolution, KU Leuven \\
#' \url{https://bio.kuleuven.be/eeb/tw/research} \\
#' ORCID: \href{https://orcid.org/0000-0002-1434-861X}{0000-0002-1434-861X}
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
  
  mz_names <- colnames(S)
  m  <- ncol(S)
  
  ## -------- pick spectra --------------------------------------------------
  set.seed(seed_compounds)
  if (k > nrow(S))
    stop("k exceeds number of available spectra in library.")
  idx <- sample(nrow(S), k)
  H_true <- Matrix(S[idx, , drop = FALSE], sparse = TRUE)
  H_true <- sweep(H_true, 1, rowSums(H_true), "/")     # so that rowSums(H_true)=1
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
  
  # map minutes -> scan indices 1..n
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
  k_background <- max(4L, ceiling(n / 500)) # number of knots
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
