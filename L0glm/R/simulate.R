#' Simulate sparse GLMs and deconvolution problems
#'
#' @description
#' Generate blurred spike trains or regression data with Gaussian predictors.
#' Returned data and ground truth can be passed to [sparseglm()], or to
#' [L0glm.fit()] and [L0glm.tune()] with the returned family and intercept flag.
#'
#' @param method Setup: `"deconvolution"`, `"normal"`, `"abess"`, or
#'   `"bestsubset"`. Numbered setups remain available: 1/2 are single channel
#'   Gaussian/Poisson deconvolution, 3/4 are their multichannel versions, and 5
#'   is independent normal-design Gaussian regression.
#' @param n,p Numbers of observations and predictors, excluding any intercept.
#' @param s Fraction of predictor rows selected if `n_nonzero` is omitted.
#'   Numbered deconvolution setups retain their original `round(s*n)` rule,
#'   capped at `min(n,p)`; named setups use `round(s*p)`.
#'   Numbered setup 5 retains `min(10,p)` unless a nonzero count is supplied.
#' @param sparse Return sparse design and coefficient matrices? `NULL` uses
#'   sparse storage for deconvolution and dense storage for Gaussian designs.
#'   Sparse storage usually wastes memory for a fully populated Gaussian design.
#' @param options Named list of setup-specific arguments; see Details.
#' @param seed Nonnegative integer seed. The caller's RNG state is restored
#'   on exit. `NULL` uses and advances the current random-number stream.
#' @param family An R GLM family object/function, or a family name. Names include
#'   `"gaussian"`, `"binomial"`, `"poisson"`, `"Gamma"`, `"inverse.gaussian"`,
#'   `"nb"`, `"multinomial"` and `"cox"`. Positive-response regression families
#'   default to a log link; deconvolution uses identity links. Special-family
#'   objects from \pkg{mgcv} are accepted. `NULL` uses Poisson for named
#'   deconvolution, Gaussian for regression, or the numbered setup's family.
#' @param n_nonzero Number of nonzero predictor rows, excluding the intercept.
#'   Overrides `s`; a supplied `beta` determines its own support. Deconvolution
#'   caps this number at `min(n,p)`. See `beta.type` below for weak sparsity.
#' @param intercept `FALSE`/`NULL` omits an intercept; `TRUE` uses one.
#'   A finite numeric scalar or one value per channel/category supplies its
#'   value, including zero. Multinomial intercepts are centered across categories.
#'   For Cox this is a log baseline-hazard multiplier, not a regression intercept.
#' @param m Number of response channels, or multinomial/ordinal categories.
#'   A multinomial family object determines the category count. Numbered setups
#'   3/4 retain all informative library mass channels when `m` is omitted.
#'   Explicit `m>1` deconvolution uses random sparse, normalized spectra with
#'   exactly that many channels.
#' @param beta Optional finite slope vector or `p` by `m` matrix, in original
#'   predictor units, for Gaussian-design setups. A vector is repeated across
#'   channels. Multinomial rows are centered to CLR coordinates. Ordinal models
#'   use a single slope vector shared across categories.
#'
#' @details
#' **Deconvolution.** Options are arguments of [simulate_spike_train()], including
#' `mean_beta`, `sd_logbeta`, `sd_noise`, `psf`, `thresh`, `theta`, `clip`,
#' `spectra` and `channel.sparsity`. Gaussian responses are untruncated by default.
#' Only numbered setups 1/3 retain historical zero clipping, unless
#' `options$clip=FALSE` is supplied. With clipping, `y_true` is the latent
#' uncensored mean.
#'
#' **Gaussian predictors.** `"normal"` uses independent standard normal columns
#' and random support with uniform `beta_range` (default `c(-1,1)`). All regression
#' setups accept `rho` and `correlation` (`"independent"`, `"ar1"`,
#' `"equicorrelated"`), or ABESS `cortype=1:3`. AR(1) covariance is
#' \eqn{\Sigma_{ij}=\rho^{|i-j|}}; equicorrelation has off-diagonal entries
#' \eqn{\rho}. Generation and population signal-variance calculations avoid
#' allocating a `p` by `p` covariance matrix.
#'
#' **ABESS setups.** `"abess"` follows predictor and coefficient distributions
#' in `abess::generate.data()`: random support and family-specific uniform ranges
#' proportional to \eqn{\sqrt{2\log(p)/n}}. `support.size`, `y.dim`, `class.num`,
#' `sigma`, `snr`, `weibull.shape` and `uniform.max` are recognized. Non-Gaussian
#' responses are sampled directly from the requested family, without the source's
#' extra latent Gaussian noise or linear-predictor clipping; therefore `y_true`
#' and `beta_true` describe the same specified model. `"mgaussian"` chooses
#' Gaussian multichannel data. `"ordinal"` generates proportional-odds data
#' with ordered `cutpoints` (default equally spaced from -1 to 1); L0glm does
#' not fit ordinal models. These are statistical setups, not upstream RNG replicas.
#'
#' **Best-subset setups.** `"bestsubset"` implements the `bestsubset::sim.xy()`
#' coefficient patterns. `beta.type=1` places unit coefficients at equally spaced
#' locations; 2 uses the first selected rows; 3 decreases coefficients linearly
#' from 10 to 0.5; 4 uses `c(-10,-6,-2,2,6,10)` in the first six rows (requires
#' `p>=6`, ignoring `n_nonzero`); 5 uses unit coefficients in the first selected
#' rows followed by `0.5^(1:(p-n_nonzero))`. Type 5 is weakly sparse: trailing
#' coefficients are also nonzero until numerical underflow. Nonzero `rho`
#' defaults to AR(1) here. `nval` (default zero) adds independent validation data
#' from the same model, for any Gaussian-design setup, without changing the
#' training sample. Multichannel and
#' non-Gaussian outcomes extend the source's Gaussian model. For multinomial
#' data the patterns multiply centered category contrasts.
#'
#' **Noise and responses.** `dispersion` is a scalar or channel vector, default
#' one: Gaussian variance, Gamma squared coefficient of variation, or inverse
#' Gaussian dispersion. Gaussian `error_sd`, `sd_noise` or `sigma` specifies SD
#' instead. Alternatively `snr` sets Gaussian variance to
#' \eqn{\beta^T\Sigma\beta/\mathrm{snr}} for identity-link Gaussian models
#' (defaults 10 for ABESS, 1 for bestsubset). SNR is not applied to non-Gaussian
#' outcomes. `theta` is positive
#' NB shape (default 10, or the supplied fixed family shape); conflicting fixed
#' family and `theta` values are rejected. NB models have
#' \eqn{\mathrm{Var}(Y)=\mu+\mu^2/\theta}. `trials` (default one) sets binomial
#' or multinomial totals; `trialsval` optionally supplies validation totals.
#' Binomial responses are proportions with trial totals
#' returned as `weights`; multinomial responses are count matrices, or factors
#' with `response="factor"` and one trial. Multinomial coefficients use CLR.
#' Quasipoisson uses Poisson/NB draws matching \eqn{\phi\mu}; quasibinomial uses
#' binomial/beta-binomial draws matching \eqn{\phi\mu(1-\mu)/\mathrm{trials}},
#' requiring attainable overdispersion. Unrestricted `quasi()` is rejected
#' because it does not specify a sampling distribution.
#'
#' Cox data have \eqn{H_0(t)=a t^b}: `baseline.rate=a`, `weibull.shape=b`,
#' both default one. Censoring is independently uniform on `[0, uniform.max]`,
#' with default upper bound one. An intercept multiplies `a` by `exp(intercept)`.
#' `y_true` is relative risk, not an event time; `event_time` contains sampled
#' uncensored times. No Cox intercept row is added to `beta_true`. Invalid means
#' for a requested link are rejected rather than clipped; choose a valid link,
#' intercept and coefficients.
#'
#' @return A list with `X` (without an intercept column), `y`, `y_true`
#'   (conditional means, category probabilities, or Cox relative risks),
#'   `eta_true`, `beta_true` (including any identifiable intercept row),
#'   `intercept` (logical fitting flag), `intercept_true`, `family`, and `weights`.
#'   Unit weights have one shared column, so subsetting response channels does
#'   not require subsetting them; binomial weights contain the trial totals.
#'   Regression data include `covariance` (structure and rho), generating
#'   `dispersion`, and NB `theta`. With `nval>0`, `Xval`, `yval`, `yval_true`,
#'   `eta_val` and `weightsval` describe validation data. Deconvolution also
#'   returns time and spike grids. Cox responses are right-censored `Surv` objects.
#' @references
#' ABESS source: \url{https://github.com/abess-team/abess/blob/master/R-package/R/generate.data.R}.
#' Hastie, Tibshirani and Tibshirani, best-subset source:
#' \url{https://github.com/ryantibs/best-subset/blob/master/bestsubset/R/sim.R}.
#' @seealso [simulate_spike_train()], [L0glm.fit()], [sparseglm()]
#' @importFrom statmod rinvgauss
#' @importFrom stats gaussian plogis poisson qlogis rbinom rgamma rpois runif
#' @importFrom methods as is
#' @importFrom utils flush.console
#' @example examples/simul.R
#' @md
#' @export
simul <- function(method = "deconvolution", n = 500, p = 500, s = 0.1,
                  sparse = NULL, options = list(), seed = 1, family = NULL,
                  n_nonzero = NULL, intercept = FALSE, m = 1L, beta = NULL) {
  n <- L0glm_sim_integer(n, "n", 1); p <- L0glm_sim_integer(p, "p", 1)
  if (is.null(options)) options <- list()
  if (!is.list(options) || (length(options) && (is.null(names(options)) ||
      any(!nzchar(names(options))) || anyDuplicated(names(options)))))
    stop("options must be a uniquely named list.", call.=FALSE)
  if (is.numeric(method)) {
    method <- L0glm_sim_integer(method, "method", 1, 5)
    if (is.null(family)) family <- c("gaussian", "poisson", "gaussian", "poisson", "gaussian")[method]
    if (method %in% c(3L, 4L) && missing(m)) m <- NULL
    if (method < 5L) {
      if (is.null(options$mean_beta)) options$mean_beta <- if (method < 3L) 2000 else 10000
      if (is.null(options$sd_noise) && method %in% c(1L, 3L)) options$sd_noise <- 5
      if (is.null(options$clip) && method %in% c(1L, 3L)) options$clip <- TRUE
      if (is.null(n_nonzero)) n_nonzero <- min(n, p, max(0, round(s*n)))
    } else if (is.null(n_nonzero) && is.null(options$n_nonzero) && is.null(options$support.size))
      n_nonzero <- min(10L, p)
    method <- if (method < 5L) "deconvolution" else "normal"
  }
  method <- match.arg(method, c("deconvolution", "normal", "abess", "bestsubset"))
  if (is.null(sparse)) sparse <- method == "deconvolution"
  L0glm_sim_flag(sparse, "sparse")
  if (is.null(n_nonzero)) n_nonzero <- options$n_nonzero
  if (is.null(n_nonzero)) n_nonzero <- options$support.size
  if (is.null(n_nonzero)) {
    if (!is.numeric(s) || length(s) != 1L || !is.finite(s) || s < 0 || s > 1)
      stop("s must be between zero and one.", call.=FALSE)
    n_nonzero <- round(s*p)
  }
  if (is.null(beta)) n_nonzero <- L0glm_sim_integer(n_nonzero, "n_nonzero", 0,
    if (method == "deconvolution") .Machine$integer.max else p)
  if (method == "deconvolution") {
    L0glm_sim_options(options, c("mean_beta", "sd_logbeta", "sd_noise", "psf", "thresh",
      "theta", "clip", "spectra", "channel.sparsity", "Plot", "n_nonzero", "support.size"))
    if (!is.null(beta)) stop("Use spectra and mean_beta for deconvolution coefficients.", call.=FALSE)
    args <- options[setdiff(names(options), c("n_nonzero", "support.size"))]
    if (missing(m) && !is.null(options$spectra)) m <- ncol(options$spectra)
    if (!is.null(m)) m <- L0glm_sim_integer(m, "m", 1)
    return(do.call(simulate_spike_train, c(list(n=n, p=p, k=n_nonzero,
      family=if (is.null(family)) "poisson" else family, sparse=sparse, seed=seed,
      intercept=intercept, multichannel=is.null(m) || m > 1L, m=m), args)))
  }
  L0glm_sim_options(options, c("n_nonzero", "support.size", "beta_range", "rho",
    "correlation", "cortype", "beta.type", "snr", "sigma", "error_sd", "sd_noise",
    "dispersion", "theta", "trials", "trialsval", "response", "nval", "y.dim", "class.num",
    "cutpoints", "weibull.shape", "uniform.max", "baseline.rate"))
  restore <- L0glm_sim_seed(seed); on.exit(restore(), add=TRUE)
  spec <- L0glm_sim_family(family, m, options, missing(m))
  family <- spec$family; kind <- spec$kind; m <- spec$m
  if (kind == "ordinal") options$class.num <- m
  q <- if (kind == "ordinal") 1L else m
  nval <- L0glm_sim_integer(if (is.null(options$nval)) 0 else options$nval, "nval")
  covariance <- L0glm_sim_covariance(options, method, p)
  X <- L0glm_sim_design(n, p, covariance)
  B <- L0glm_sim_beta(beta, p, q, n_nonzero, method, kind, n, options)
  if (kind == "multinomial") B <- B-rowMeans(B)
  constant <- L0glm_sim_intercept(intercept, q, kind)
  phi <- L0glm_sim_parameter(options$dispersion, q, 1, "dispersion", zero=kind == "gaussian")
  if (kind %in% c("poisson", "binomial", "nb", "multinomial", "cox", "ordinal") && any(phi != 1))
    stop("This family has fixed dispersion=1; choose a quasi-family or NB for overdispersion.", call.=FALSE)
  if (kind == "gaussian") {
    sd <- options$error_sd
    if (is.null(sd)) sd <- options$sd_noise
    if (is.null(sd)) sd <- options$sigma
    snr <- options$snr
    if (is.null(snr) && is.null(options$dispersion) && family$link == "identity" && method != "normal")
      snr <- if (method == "abess") 10 else 1
    if (!is.null(sd)) phi <- L0glm_sim_parameter(sd, q, 1, "noise SD", zero=TRUE)^2 else
      if (!is.null(snr)) {
        if (family$link != "identity") stop("snr requires a Gaussian identity link; use a noise SD or dispersion.", call.=FALSE)
        phi <- L0glm_sim_signal_variance(B, covariance)/L0glm_sim_parameter(snr, q, 1, "snr")
      }
  }
  theta <- if (kind == "nb") L0glm_sim_nb_shape(family, options$theta, q) else NULL
  if (kind == "cox") {
    rate <- L0glm_sim_parameter(options$baseline.rate, 1, 1, "baseline.rate")*exp(constant$value)
    options$baseline.rate <- L0glm_sim_parameter(rate, 1, 1, "effective baseline.rate")
  }
  draw <- function(design, validation=FALSE) {
    eta <- as.matrix(design %*% B)
    if (kind != "cox") eta <- sweep(eta, 2L, constant$value, "+")
    sampling <- options
    if (validation && !is.null(options$trialsval)) sampling$trials <- options$trialsval
    L0glm_sim_response(eta, family, kind, phi, theta, sampling)
  }
  answer <- c(list(X=X, beta_true=if (constant$include) rbind(Intercept=constant$value, B) else B,
    intercept=constant$include, intercept_true=constant$value, family=family,
    covariance=covariance, dispersion=phi, theta=theta), draw(X))
  if (nval) {
    Xval <- L0glm_sim_design(nval, p, covariance)
    val <- draw(Xval, validation=TRUE)
    answer[c("Xval", "yval", "yval_true", "eta_val", "weightsval")] <-
      list(Xval, val$y, val$y_true, val$eta_true, val$weights)
  }
  if (sparse) {
    answer$X <- Matrix::Matrix(answer$X, sparse=TRUE)
    answer$beta_true <- Matrix::Matrix(answer$beta_true, sparse=TRUE)
    if (nval) answer$Xval <- Matrix::Matrix(answer$Xval, sparse=TRUE)
  }
  answer
}

L0glm_sim_integer <- function(x, name, minimum=0, maximum=.Machine$integer.max) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x != floor(x) || x < minimum || x > maximum)
    stop(name, " must be an integer between ", minimum, " and ", maximum, ".", call.=FALSE)
  as.integer(x)
}

L0glm_sim_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) stop(name, " must be TRUE or FALSE.", call.=FALSE)
}

L0glm_sim_options <- function(options, allowed) {
  unknown <- setdiff(names(options), allowed)
  if (length(unknown)) stop("Unknown simulation option: ", paste(unknown, collapse=", "), call.=FALSE)
}

L0glm_sim_parameter <- function(x, m, default, name, zero=FALSE) {
  if (is.null(x)) x <- default
  if (!is.numeric(x) || !length(x) %in% c(1L, m) || any(!is.finite(x)) || any(if (zero) x < 0 else x <= 0))
    stop(name, " must be finite ", if (zero) "nonnegative" else "positive",
      " numbers, one or one per channel.", call.=FALSE)
  rep(x, length.out=m)
}

L0glm_sim_seed <- function(seed) {
  if (is.null(seed)) return(function() NULL)
  seed <- L0glm_sim_integer(seed, "seed")
  had_seed <- exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
  old <- if (had_seed) get(".Random.seed", envir=.GlobalEnv) else NULL
  set.seed(seed)
  function() {
    if (had_seed) assign(".Random.seed", old, envir=.GlobalEnv) else
      if (exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE)) rm(".Random.seed", envir=.GlobalEnv)
  }
}

L0glm_sim_nb_shape <- function(family, theta, m) {
  fixed <- isTRUE(family$n.theta == 0)
  value <- L0glm_sim_parameter(theta, m, if (fixed) family$getTheta(TRUE) else 10, "theta")
  if (fixed) {
    if (any(abs(value/family$getTheta(TRUE)-1) > 1e-12))
      stop("theta conflicts with the supplied fixed-shape NB family.", call.=FALSE)
    value[] <- family$getTheta(TRUE)
  }
  value
}

L0glm_sim_family <- function(family, m, options, default_m=FALSE) {
  name <- if (is.null(family)) "gaussian" else if (is.character(family)) tolower(family) else NULL
  if (!is.null(name)) {
    if (length(name) != 1L || is.na(name)) stop("Supply one family.", call.=FALSE)
    if (name == "mgaussian") {
      if (default_m) m <- if (is.null(options$y.dim)) 3L else options$y.dim
      name <- "gaussian"
    }
    if (name %in% c("multinomial", "ordinal") && default_m)
      m <- if (is.null(options$class.num)) 3L else options$class.num
    m <- L0glm_sim_integer(m, "m", if (name %in% c("multinomial", "ordinal")) 2 else 1)
    family <- switch(name, gaussian=stats::gaussian(), binomial=stats::binomial(),
      poisson=stats::poisson(), gamma=stats::Gamma(link="log"),
      inverse.gaussian=stats::inverse.gaussian(link="log"),
      quasipoisson=stats::quasipoisson(), quasibinomial=stats::quasibinomial(),
      nb=, negative.binomial=, `negative binomial`=, multinomial=, cox=, cox.ph={
        if (!requireNamespace("mgcv", quietly=TRUE)) stop("Install mgcv for this family.", call.=FALSE)
        if (name %in% c("cox", "cox.ph")) mgcv::cox.ph() else
          if (name == "multinomial") mgcv::multinom(K=m-1L) else mgcv::nb()
      }, ordinal="ordinal", stop("Unsupported simulation family: ", name, call.=FALSE))
  }
  if (is.function(family)) family <- family()
  if (identical(family, "ordinal")) kind <- "ordinal" else {
    if (!inherits(family, "family")) stop("Supply a GLM family object or name.", call.=FALSE)
    kind <- L0glm_ext_kind(family)
    if (is.null(kind)) kind <- family$family
    if (kind == "multinomial") m <- family$nlp+1L
  }
  m <- L0glm_sim_integer(m, "m", if (kind %in% c("multinomial", "ordinal")) 2 else 1)
  if (!kind %in% c("gaussian", "binomial", "poisson", "Gamma", "inverse.gaussian",
      "quasipoisson", "quasibinomial", "nb", "multinomial", "cox", "ordinal"))
    stop("This family does not specify a supported sampling distribution.", call.=FALSE)
  if (kind == "cox" && m != 1L) stop("Cox has one event-time response; use m=1.", call.=FALSE)
  list(family=family, kind=kind, m=m)
}

L0glm_sim_intercept <- function(intercept, m, kind) {
  include <- !is.null(intercept) && !identical(intercept, FALSE)
  if (isTRUE(intercept)) intercept <- 1
  value <- if (include) intercept else 0
  if (!is.numeric(value) || !length(value) %in% c(1L, m) || any(!is.finite(value)))
    stop("intercept must be FALSE, TRUE, or finite channel/category values.", call.=FALSE)
  value <- rep(value, length.out=m)
  if (kind == "multinomial") value <- value-mean(value)
  list(include=include && kind != "cox", value=value)
}

L0glm_sim_covariance <- function(options, method, p) {
  type <- options$correlation
  if (!is.null(options$cortype)) {
    cortype <- L0glm_sim_integer(options$cortype, "cortype", 1, 3)
    if (!is.null(type)) stop("Supply correlation or cortype, not both.", call.=FALSE)
    type <- c("independent", "ar1", "equicorrelated")[cortype]
  }
  rho <- if (is.null(options$rho)) 0 else options$rho
  if (!is.numeric(rho) || length(rho) != 1L || !is.finite(rho) || abs(rho) >= 1)
    stop("rho must be between -1 and 1.", call.=FALSE)
  if (is.null(type)) type <- if (method == "abess" || rho == 0) "independent" else "ar1"
  type <- match.arg(type, c("independent", "ar1", "equicorrelated"))
  if (type == "equicorrelated" && p > 1 && rho < -1/(p-1))
    stop("rho is outside the valid equicorrelation range.", call.=FALSE)
  if (type == "independent") rho <- 0
  list(type=type, rho=rho)
}

L0glm_sim_design <- function(n, p, covariance) {
  X <- matrix(stats::rnorm(as.double(n)*p), n, p)
  rho <- covariance$rho
  if (rho != 0 && covariance$type == "ar1" && p > 1L)
    for (j in 2:p) X[, j] <- rho*X[, j-1L]+sqrt(1-rho^2)*X[, j]
  if (rho != 0 && covariance$type == "equicorrelated")
    X <- sqrt(1-rho)*X+(sqrt(1+(p-1)*rho)-sqrt(1-rho))*rowMeans(X)
  colnames(X) <- paste0("x", seq_len(p))
  X
}

L0glm_sim_signal_variance <- function(B, covariance) {
  square <- colSums(B^2); rho <- covariance$rho
  if (rho == 0) return(square)
  if (covariance$type == "equicorrelated") return(pmax(0, (1-rho)*square+rho*colSums(B)^2))
  vapply(seq_len(ncol(B)), function(j) max(0,
    2*sum(B[, j]*as.numeric(stats::filter(B[, j], rho, method="recursive")))-square[j]), numeric(1))
}

L0glm_sim_beta <- function(beta, p, m, k, method, kind, n, options) {
  if (!is.null(beta)) {
    if (inherits(beta, "Matrix")) beta <- as.matrix(beta)
    if (!is.numeric(beta) || any(!is.finite(beta)) ||
        !(length(beta) == p && is.null(dim(beta)) || identical(dim(beta), c(p, m))))
      stop("beta must be a finite slope vector or p by m matrix.", call.=FALSE)
    B <- matrix(beta, p, m)
  } else if (method == "bestsubset") {
    type <- if (is.null(options$beta.type)) 1L else L0glm_sim_integer(options$beta.type, "beta.type", 1, 5)
    b <- numeric(p)
    if (type == 4L) {
      if (p < 6L) stop("beta.type=4 requires p>=6.", call.=FALSE)
      b[1:6] <- c(-10, -6, -2, 2, 6, 10)
    } else {
      if (k) b[if (type == 1L) round(seq(1, p, length.out=k)) else seq_len(k)] <-
        if (type == 3L) seq(10, .5, length.out=k) else 1
      if (type == 5L && k < p) b[(k+1L):p] <- .5^seq_len(p-k)
    }
    B <- if (kind == "multinomial") outer(b, seq_len(m)-mean(seq_len(m))) else matrix(b, p, m)
  } else {
    range <- options$beta_range
    if (is.null(range)) {
      a <- sqrt(2*log(max(2, p))/n)
      range <- if (method == "normal") c(-1, 1) else switch(kind,
        gaussian=5*a*c(1, 100), binomial=, quasibinomial=, cox=5*a*c(2, 10),
        multinomial=2.5*a*c(1, 50), ordinal=2.5*a*c(-50, 50),
        Gamma=5*a*c(1, 100), a*c(-10, 10))
    }
    if (!is.numeric(range) || length(range) != 2L || any(!is.finite(range)) || range[1L] > range[2L])
      stop("beta_range must be two ordered finite numbers.", call.=FALSE)
    B <- matrix(0, p, m)
    if (k) {
      active <- sample.int(p, k)
      B[active, ] <- stats::runif(k*m, range[1L], range[2L])
    }
  }
  rownames(B) <- paste0("x", seq_len(p))
  colnames(B) <- if (kind == "multinomial") as.character(seq_len(m)-1L) else paste0("y", seq_len(m))
  B
}

L0glm_sim_response <- function(eta, family, kind, phi, theta, options) {
  n <- nrow(eta); m <- ncol(eta)
  if (kind == "cox") {
    if (!requireNamespace("survival", quietly=TRUE)) stop("Install survival for Cox data.", call.=FALSE)
    shape <- L0glm_sim_parameter(options$weibull.shape, 1, 1, "weibull.shape")
    limit <- L0glm_sim_parameter(options$uniform.max, 1, 1, "uniform.max")
    event <- exp((log(stats::rexp(n))-log(options$baseline.rate)-eta[, 1L])/shape)
    censor <- stats::runif(n, max=limit)
    return(list(y=survival::Surv(pmin(event, censor), event <= censor), y_true=exp(eta),
      eta_true=eta, weights=matrix(1, n, 1L), event_time=event, censor_time=censor,
      baseline.rate=options$baseline.rate, weibull.shape=shape))
  }
  if (kind == "ordinal") {
    cut <- if (is.null(options$cutpoints)) seq(-1, 1, length.out=options$class.num-1L) else options$cutpoints
    if (!is.numeric(cut) || length(cut) != options$class.num-1L || any(!is.finite(cut)) || is.unsorted(cut, strictly=TRUE))
      stop("cutpoints must be strictly increasing, one fewer than the category count.", call.=FALSE)
    cumulative <- matrix(vapply(cut, function(a) stats::plogis(a-eta[, 1L]), numeric(n)), n)
    mu <- cbind(cumulative, 1)-cbind(0, cumulative)
  } else if (kind == "multinomial") mu <- L0glm_ext_softmax(eta)$p else {
    if (is.function(family$valideta) && !family$valideta(as.numeric(eta)))
      stop("Invalid linear predictor for this link; change beta or intercept.", call.=FALSE)
    mu <- matrix(family$linkinv(as.numeric(eta)), n, m)
    valid <- all(is.finite(mu)) && if (kind %in% c("poisson", "quasipoisson", "nb")) all(mu >= 0) else
      if (kind %in% c("Gamma", "inverse.gaussian")) all(mu > 0) else
        if (kind %in% c("binomial", "quasibinomial")) all(mu >= 0 & mu <= 1) else TRUE
    if (!valid) stop("Invalid means for this family/link; change beta, intercept or link. Means are not clipped.", call.=FALSE)
  }
  if (any(!is.finite(mu))) stop("Non-finite simulation means.", call.=FALSE)
  weights <- matrix(1, n, 1L)
  if (kind %in% c("multinomial", "ordinal")) {
    trials <- L0glm_sim_parameter(options$trials, n, 1, "trials")
    if (any(trials != floor(trials))) stop("trials must be positive integers.", call.=FALSE)
    y <- t(vapply(seq_len(n), function(i) as.numeric(stats::rmultinom(1L, trials[i], mu[i, ])), numeric(ncol(mu))))
    format <- if (is.null(options$response)) "counts" else match.arg(options$response, c("counts", "factor"))
    if (kind == "ordinal" || format == "factor") {
      if (any(trials != 1)) stop("Factor responses require trials=1.", call.=FALSE)
      y <- factor(max.col(y), levels=seq_len(ncol(mu)), labels=as.character(seq_len(ncol(mu))-1L), ordered=kind == "ordinal")
    }
    weights <- matrix(1, n, 1L)
  } else if (kind %in% c("binomial", "quasibinomial")) {
    weights <- L0glm_ext_matrix(options$trials, n, m, 1, "trials")
    if (any(!is.finite(weights) | weights < 1 | weights != floor(weights))) stop("trials must be positive integers.", call.=FALSE)
    probability <- mu
    if (kind == "quasibinomial") for (j in which(phi != 1)) {
      if (phi[j] < 1 || any(weights[, j] <= phi[j]))
        stop("Beta-binomial dispersion must be >=1 and smaller than every trial total.", call.=FALSE)
      concentration <- (weights[, j]-phi[j])/(phi[j]-1)
      good <- mu[, j] > 0 & mu[, j] < 1
      probability[good, j] <- stats::rbeta(sum(good), mu[good, j]*concentration[good], (1-mu[good, j])*concentration[good])
    }
    y <- matrix(stats::rbinom(n*m, weights, probability), n, m)/weights
  } else {
    dispersion <- rep(phi, each=n)
    y <- switch(kind, gaussian=stats::rnorm(n*m, mu, sqrt(dispersion)),
      poisson=stats::rpois(n*m, mu), nb=stats::rnbinom(n*m, mu=mu, size=rep(theta, each=n)),
      Gamma=stats::rgamma(n*m, shape=1/dispersion, scale=as.numeric(mu)*dispersion),
      inverse.gaussian=statmod::rinvgauss(n*m, mean=mu, shape=1/dispersion),
      quasipoisson={
        if (any(phi < 1)) stop("Quasipoisson simulation requires dispersion>=1.", call.=FALSE)
        yy <- stats::rpois(n*m, mu); good <- dispersion > 1 & mu > 0
        yy[good] <- stats::rnbinom(sum(good), mu=mu[good], size=mu[good]/(dispersion[good]-1))
        yy
      })
    y <- matrix(y, n, m)
  }
  if (is.numeric(y) && (any(!is.finite(y)) ||
      kind %in% c("Gamma", "inverse.gaussian") && any(y <= 0)))
    stop("Response generation overflowed or underflowed; reduce beta magnitudes or change the link.", call.=FALSE)
  if (is.matrix(y)) colnames(y) <- if (kind == "multinomial") as.character(seq_len(ncol(y))-1L) else paste0("y", seq_len(ncol(y)))
  colnames(mu) <- if (kind %in% c("multinomial", "ordinal")) as.character(seq_len(ncol(mu))-1L) else colnames(y)
  list(y=y, y_true=mu, eta_true=eta, weights=weights)
}
