#' Simulate a blurred single or multichannel spike train
#'
#' @description
#' Generate a sparse spike train blurred by a known point spread function (PSF),
#' with Poisson, Gaussian or negative binomial noise. Columns of the design are
#' shifted copies of the PSF. This models signal deconvolution, including GC/MS
#' data with sparse mass spectra shared over time.
#' @param n,p Numbers of time points and predictor locations; `p>n` allows
#'   superresolution. The design excludes any intercept column.
#' @param k Number of peaks, capped at `min(n,p)`.
#' @param mean_beta Median total peak height (the log-normal scale parameter).
#' @param sd_logbeta Standard deviation of log peak heights. Heights follow
#'   `exp(rnorm(k, log(mean_beta), sd_logbeta))`; their arithmetic expectation is
#'   `mean_beta*exp(sd_logbeta^2/2)`.
#' @param family `"poisson"`, `"gaussian"`, `"nb"` or `"negative.binomial"`,
#'   or a corresponding family object with an identity link. NB family objects
#'   are supplied by `mgcv::nb(link="identity")`.
#' @param sd_noise Gaussian noise standard deviation, one value or one per channel.
#' @param psf Function of time `x` and peak location `u`. The default Gaussian
#'   has height one and width 2.5. Values must be finite and nonnegative.
#' @param thresh PSF entries smaller than this nonnegative cutoff are zeroed.
#' @param sparse Return a sparse design/coefficient matrix? Defaults to `TRUE`.
#'   The sparse design is constructed directly without a dense `n` by `p` copy.
#' @param multichannel Generate multiple spectral channels? With `m=NULL`, use
#'   the package's [simulated_spectra] library, retaining informative channels.
#' @param seed Random seed; a specified seed preserves the caller's RNG state.
#'   `NULL` uses and advances the current stream.
#' @param Plot Plot the generated signal or multichannel heatmap?
#' @param intercept `FALSE`/`NULL` omits an intercept; `TRUE` adds a baseline of
#'   one; a scalar or channel vector supplies its value, including zero. The
#'   corresponding row is included in `beta_true`, but not in `X`.
#' @param theta Positive NB shape, one value or one per channel. Defaults to 10
#'   or the supplied family's fixed shape. NB variance is `mu+mu^2/theta`.
#' @param m Explicit number of channels. Values greater than one enable
#'   multichannel generation using random sparse spectra, unless `spectra` is
#'   provided. `NULL` selects the library when `multichannel=TRUE`.
#' @param spectra Optional nonnegative matrix of template spectra, with one
#'   spectrum per row. Rows are normalized to sum to one, then `k` rows are
#'   sampled, with replacement only if the library contains fewer than `k` rows.
#'   Channels are retained, including channels with zero generated signal.
#' @param channel.sparsity Probability of a zero entry in random spectra
#'   (`m` supplied); defaults to 0.7. Each peak retains at least one channel.
#' @param clip Clip negative Gaussian observations to zero? Defaults to `FALSE`:
#'   clipping produces censored data, not Gaussian noise. Numbered `simul()`
#'   setups 1/3 explicitly enable clipping to reproduce historical benchmarks.
#' @return A list with `X`, `beta_true`, `y`, `y_true`, `eta_true`, `family`,
#'   `weights`, `intercept`, `intercept_true`, `dispersion`, `theta`, `spectra`,
#'   time grid `x`, and predictor-location grid `x_beta`. Means satisfy
#'   `y_true = X %*% slopes + intercept_true`. `beta_true` includes the intercept
#'   row when requested. Use `intercept=sim$intercept` when fitting `sim$X`.
#'   With `clip=TRUE`, `y_true` describes the latent uncensored signal.
#' @seealso [simul()], [L0glm.fit()], [sparseglm()]
#' @importFrom Matrix Matrix
#' @importFrom qlcMatrix colMax colMin
#' @importFrom grDevices topo.colors
#' @importFrom graphics abline matplot points
#' @md
#' @export
simulate_spike_train <- function(n=100, p=n, k=10, mean_beta=10000, sd_logbeta=1,
  family="poisson", sd_noise=10,
  psf=function(x, u, w=2.5) exp((x-u)^2/(-2*w^2)), thresh=1e-5,
  sparse=TRUE, multichannel=FALSE, seed=123, Plot=FALSE, intercept=FALSE,
  theta=NULL, m=NULL, spectra=NULL, channel.sparsity=.7, clip=FALSE) {
  n <- L0glm_sim_integer(n, "n", 1); p <- L0glm_sim_integer(p, "p", 1)
  k <- min(L0glm_sim_integer(k, "k"), n, p)
  for (name in c("sparse", "multichannel", "Plot", "clip")) L0glm_sim_flag(get(name), name)
  mean_beta <- L0glm_sim_parameter(mean_beta, 1, 10000, "mean_beta")
  sd_logbeta <- L0glm_sim_parameter(sd_logbeta, 1, 1, "sd_logbeta", zero=TRUE)
  thresh <- L0glm_sim_parameter(thresh, 1, 1e-5, "thresh", zero=TRUE)
  if (!is.function(psf)) stop("psf must be a function of x and u.", call.=FALSE)
  if (!is.null(m)) m <- L0glm_sim_integer(m, "m", 1)
  if (is.character(family)) {
    name <- tolower(family)
    if (length(name) != 1L || is.na(name)) stop("Supply one deconvolution family.", call.=FALSE)
    family <- switch(name, poisson=stats::poisson("identity"), gaussian=stats::gaussian(),
      nb=, negative.binomial=, `negative binomial`={
        if (!requireNamespace("mgcv", quietly=TRUE)) stop("Install mgcv for NB data.", call.=FALSE)
        mgcv::nb(link="identity")
      }, stop("Deconvolution supports Gaussian, Poisson or NB noise.", call.=FALSE))
  }
  if (is.function(family)) family <- family()
  spec <- L0glm_sim_family(family, 1L, list())
  if (!spec$kind %in% c("gaussian", "poisson", "nb") || family$link != "identity")
    stop("Deconvolution requires Gaussian, Poisson or NB with an identity link.", call.=FALSE)
  restore <- L0glm_sim_seed(seed); on.exit(restore(), add=TRUE)
  x <- seq_len(n); x_beta <- seq(1, n, length.out=p)
  u <- sample.int(p, k)
  library_spectra <- is.null(spectra) && is.null(m) && multichannel
  if (!is.null(spectra) || library_spectra) {
    if (library_spectra) {
      utils::data("simulated_spectra", package="L0glm", envir=environment())
      spectra <- simulated_spectra
    }
    if (!is.numeric(spectra) && !inherits(spectra, "Matrix") || length(dim(spectra)) != 2L ||
        !nrow(spectra) || !ncol(spectra) || any(!is.finite(spectra) | spectra < 0) ||
        any(Matrix::rowSums(spectra) <= 0))
      stop("spectra must be a nonnegative matrix with positive row totals.", call.=FALSE)
    if (!is.null(m) && m != ncol(spectra)) stop("m must match ncol(spectra).", call.=FALSE)
    spectra <- spectra/Matrix::rowSums(spectra)
    S <- spectra[sample.int(nrow(spectra), k, replace=k > nrow(spectra)), , drop=FALSE]
    m <- ncol(S)
  } else if (!is.null(m) && m > 1L) {
    if (!is.numeric(channel.sparsity) || length(channel.sparsity) != 1L ||
        !is.finite(channel.sparsity) || channel.sparsity < 0 || channel.sparsity > 1)
      stop("channel.sparsity must be between zero and one.", call.=FALSE)
    S <- matrix(stats::rexp(k*m), k, m)
    S[stats::runif(k*m) < channel.sparsity] <- 0
    empty <- which(rowSums(S) == 0)
    if (length(empty)) S[cbind(empty, sample.int(m, length(empty), replace=TRUE))] <- 1
    if (k) S <- S/rowSums(S)
  } else {m <- 1L; S <- matrix(1, k, m)}
  B <- Matrix::Matrix(0, p, m, sparse=sparse)
  heights <- exp(stats::rnorm(k, log(mean_beta), sd_logbeta))
  if (any(!is.finite(heights))) stop("Peak heights overflow; reduce mean_beta or sd_logbeta.", call.=FALSE)
  if (k) B[u, ] <- Matrix::Diagonal(x=heights) %*% S
  # Preserve the library benchmark's channel selection, without dropping every
  # channel in a null-signal simulation. Explicit channel counts stay fixed.
  if (library_spectra && (is.null(intercept) || identical(intercept, FALSE))) {
    keep <- as.logical(qlcMatrix::colMax(B) != qlcMatrix::colMin(B))
    if (any(keep)) {
      if (length(sd_noise) == m && m > 1L) sd_noise <- sd_noise[keep]
      if (length(theta) == m && m > 1L) theta <- theta[keep]
      B <- B[, keep, drop=FALSE]; S <- S[, keep, drop=FALSE]; m <- ncol(B)
    }
  }
  column <- function(u) {
    value <- psf(x, u)
    if (!is.numeric(value) || length(value) != n || any(!is.finite(value) | value < 0))
      stop("psf must return n finite nonnegative values.", call.=FALSE)
    value[value < thresh] <- 0
    if (sparse) {index <- which(value != 0); list(i=index, x=value[index])} else value
  }
  columns <- lapply(x_beta, column)
  if (sparse) {
    counts <- vapply(columns, function(z) length(z$i), integer(1))
    X <- Matrix::sparseMatrix(i=unlist(lapply(columns, `[[`, "i"), use.names=FALSE),
      j=rep(seq_len(p), counts), x=unlist(lapply(columns, `[[`, "x"), use.names=FALSE), dims=c(n, p))
  } else X <- do.call(cbind, columns)
  colnames(X) <- paste0("x", seq_len(p)); rownames(B) <- colnames(X)
  constant <- L0glm_sim_intercept(intercept, m, spec$kind)
  eta <- sweep(as.matrix(X %*% B), 2L, constant$value, "+")
  phi <- L0glm_sim_parameter(sd_noise, m, 10, "sd_noise", zero=TRUE)^2
  shape <- if (spec$kind == "nb") L0glm_sim_nb_shape(family, theta, m) else NULL
  response <- L0glm_sim_response(eta, family, spec$kind, phi, shape, list())
  if (spec$kind == "gaussian" && clip) response$y[response$y < 0] <- 0
  if (library_spectra && !constant$include) {
    keep <- as.logical(qlcMatrix::colMax(response$y) != qlcMatrix::colMin(response$y))
    if (any(keep)) {
      B <- B[, keep, drop=FALSE]; S <- S[, keep, drop=FALSE]
      for (name in c("y", "y_true", "eta_true")) response[[name]] <- response[[name]][, keep, drop=FALSE]
      constant$value <- constant$value[keep]; phi <- phi[keep]
      if (!is.null(shape)) shape <- shape[keep]
    }
  }
  if (constant$include) B <- rbind(Intercept=constant$value, B)
  answer <- c(list(X=X, beta_true=B, intercept=constant$include, intercept_true=constant$value,
    family=family, dispersion=if (spec$kind == "gaussian") phi else rep(1, ncol(B)),
    theta=shape, spectra=S, x=x, x_beta=x_beta), response)
  if (Plot) {
    if (ncol(answer$y) == 1L) {
      graphics::plot(x, answer$y[, 1L], type="l", xlab="Time", ylab="Signal",
        main=paste("Blurred spike train:", family$family))
      graphics::lines(x, answer$y_true[, 1L], col="orange")
      graphics::segments(x_beta[u], constant$value[1L], x_beta[u],
        constant$value[1L]+heights, col="red")
    } else {
      plotmat(answer$y, xlab="Time", ylab="Channel",
        main=paste("Blurred multichannel spike train:", family$family))
      graphics::abline(v=x_beta[u], col=grDevices::adjustcolor("cyan", .4))
    }
  }
  answer
}
