slcn_correlation <- function(correlation) {
  if(is.character(correlation)) return(match.arg(correlation,c("auto","none","lag1")))
  if(!is.numeric(correlation)||length(correlation)!=1L||!is.finite(correlation)||correlation<0||correlation>.49)
    stop("A supplied lag-one correlation must lie between zero and 0.49.",call.=FALSE)
  correlation
}

slcn_family <- function(family, dispersion = NULL, theta = NULL) {
  if (is.function(family)) family <- family(link = "identity")
  character_family <- is.character(family)
  if (is.character(family)) {
    family <- switch(tolower(family), gaussian = gaussian("identity"), poisson = poisson("identity"),
      quasipoisson = quasipoisson("identity"), gamma = Gamma("identity"),
      nb1 = list(family = "nb1", link = "identity"),
      gamma_scale = list(family = "Gamma_scale", link = "identity"),
      nb = mgcv::nb(link = "identity"), negative.binomial = mgcv::nb(link = "identity"),
      stop("Unsupported noise family."))
  }
  if (!is.list(family) || !identical(family$link, "identity")) stop("Use a supported identity-link family.")
  name <- tolower(family$family)
  if (grepl("^negative[ .]?binomial", name, ignore.case = TRUE)) name <- "nb"
  if (name == "gamma_censored") name <- "gamma"
  if (name == "gamma") name <- "Gamma"
  if (name == "gamma_scale") name <- "Gamma_scale"
  if (!character_family && !is.null(family$scalar) && is.null(dispersion))
    dispersion <- family$scalar$parameter
  kind <- match(name, c("gaussian", "poisson", "quasipoisson", "nb", "Gamma", "nb1", "Gamma_scale")) - 1L
  if (is.na(kind)) stop("Supported families: gaussian, poisson, quasipoisson, nb, Gamma, nb1 and Gamma_scale.")
  valid <- function(x) is.null(x) || (is.numeric(x) && length(x) == 1L && is.finite(x) && x != 0)
  if (!valid(dispersion) || (!is.null(dispersion) && dispersion < 0) || !valid(theta))
    stop("dispersion must be NULL or positive; theta must be NULL or a finite nonzero scalar.")
  if (!is.null(theta) && name != "nb") stop("theta applies to standard negative binomial (nb); use dispersion for nb1.")
  if (name == "nb1" && !is.null(dispersion) && dispersion < 1) stop("nb1 dispersion must be at least one.")
  if (name %in% c("poisson", "nb") && !is.null(dispersion) && dispersion != 1)
    stop("Poisson/NB2 dispersion is one; use quasipoisson for extra Poisson scale, or theta for NB2.")
  if (name == "nb" && is.null(theta) && !is.null(family$n.theta) && family$n.theta == 0L)
    theta <- family$getTheta(TRUE)
  if (name == "nb" && is.null(theta) && is.null(family$n.theta) && is.function(family$variance))
    theta <- get0(".Theta", envir = environment(family$variance), inherits = FALSE)
  list(family = family, name = name, kind = kind,
    dispersion = if (is.null(dispersion)) NULL else abs(dispersion),
    theta = if (is.null(theta)) NULL else abs(theta),
    theta_update = !is.null(theta) && theta < 0,
    dispersion_estimated = is.null(dispersion) && name %in% c("gaussian", "quasipoisson", "Gamma", "nb1", "Gamma_scale"),
    theta_estimated = name == "nb" && (is.null(theta) || theta < 0))
}

slcn_weights <- function(value, n, m, default = 1) {
  if (is.null(value)) return(default)
  if (length(value) == 1L) return(as.numeric(value))
  if (is.null(dim(value)) && length(value) == n) value <- matrix(value, n, m)
  if (!identical(dim(value), c(n, m))) stop("Weights/offset must be scalar, row vector or match y.")
  as.matrix(value)
}

slcn_pilot <- function(Y, span = 17L) {
  n <- nrow(Y);span <- min(n - (n + 1L) %% 2L, as.integer(span))
  if (span %% 2L == 0L) span <- span - 1L
  out <- apply(pmax(Y, 0), 2L, function(y) {
    h <- stats::filter(y, rep(1 / span, span), sides = 2)
    good <- which(is.finite(h))
    if (length(good) < 2L) return(rep(mean(y), n))
    stats::approx(good, as.numeric(h)[good], seq_len(n), rule = 2)$y
  })
  matrix(pmax(out, 1e-8), n, ncol(Y))
}

slcn_noise <- function(Y, spec, W, ncores, correlation, censor) {
  if ((identical(correlation, "none") || is.numeric(correlation)) &&
      (spec$name %in% c("poisson","nb") || !is.null(spec$dispersion)) &&
      (spec$name != "nb" || !is.null(spec$theta))) {
    if (spec$name %in% c("poisson","nb")) spec$dispersion <- 1
    if (is.null(spec$theta)) spec$theta <- 1
    spec$correlation <- if(is.numeric(correlation)) correlation else 0
    spec$correlation_estimate <- if(is.numeric(correlation)) NA_real_ else 0
    spec$information_factor <- 1+2*spec$correlation
    spec$noise_diagnostics <- NULL
    return(list(spec = spec, pilot = NULL))
  }
  level <- slcn_pilot(Y)
  quiet <- rowSums(level) <= median(rowSums(level))
  # A fixed subset of informative quiet entries supports a shared correlation estimate.
  positive <- level[level > 0]
  cut <- if (length(positive)) max(1e-8, min(1000, as.numeric(quantile(positive, .75)))) else 1
  eligible <- matrix(quiet, nrow(Y), ncol(Y)) & level >= cut
  if (spec$name == "gaussian") eligible <- matrix(quiet, nrow(Y), ncol(Y))
  if (length(W) != 1L) eligible <- eligible & W > 0
  if (sum(eligible) < 32L) {
    cut <- max(1e-8, as.numeric(quantile(positive, .1)))
    eligible <- matrix(quiet, nrow(Y), ncol(Y)) & level >= cut
    if (length(W) != 1L) eligible <- eligible & W > 0
  }
  noise <- tryCatch(temporal_dispersion(Y, family = gaussian("identity"), order = c(4L, 8L),
    rows = eligible, minimum_mean = 0, method = "moment", correlation = "lag1", ncores = ncores),
    error = function(e) NULL)
  rho <- if (is.null(noise)) 0 else unname(noise$pooled[["correlation"]])
  if (!is.finite(rho)) rho <- 0
  rho <- max(0, min(.49, rho))
  sigma2 <- if (is.null(noise)) NA_real_ else unname(noise$pooled[["dispersion"]])
  scale <- spec$dispersion
  if (spec$name == "gaussian" && is.null(scale)) scale <- sigma2
  if (spec$name %in% c("poisson", "nb")) scale <- 1
  if (spec$name %in% c("quasipoisson", "nb1", "Gamma_scale") && is.null(scale)) {
    q <- temporal_dispersion(Y, family = quasipoisson("identity"), order = c(4L, 8L), rows = eligible,
      minimum_mean = cut, correlation = "lag1", method = "moment", ncores = ncores)
    scale <- unname(q$pooled[["dispersion"]])
  }
  if (spec$name == "Gamma" && is.null(scale)) {
    q <- temporal_dispersion(Y, family = Gamma("identity"), order = c(4L, 8L), mu = level, rows = eligible,
      minimum_mean = cut, correlation = "lag1", method = "moment", ncores = ncores)
    scale <- unname(q$pooled[["dispersion"]])
  }
  if (is.null(scale) || !is.finite(scale) || scale <= 0) scale <- if (spec$name == "gaussian") max(1e-8, mean(Y^2) * .01) else 1
  spec$dispersion <- max(1e-8, min(1e20, scale))
  if (spec$name == "nb1") spec$dispersion <- max(1, spec$dispersion)
  if (spec$name == "nb" && is.null(spec$theta)) {
    q <- temporal_dispersion(Y, family = quasipoisson("identity"), order = c(4L, 8L), rows = eligible,
      minimum_mean = cut, correlation = "lag1", method = "moment", ncores = ncores)
    phi <- unname(q$pooled[["dispersion"]])
    reference <- sum(level[eligible]^2) / sum(level[eligible])
    if (!is.finite(reference)) reference <- max(1, mean(level))
    if (!is.finite(phi) || phi <= 1) phi <- 1 + reference / 1e6
    spec$theta <- max(1e-6, min(1e8, reference / (phi - 1)))
  }
  if (is.null(spec$theta)) spec$theta <- 1
  selected <- if (is.numeric(correlation)) max(0, min(.49, correlation)) else if (correlation == "none") 0 else rho
  spec$correlation <- selected
  spec$correlation_estimate <- rho
  spec$information_factor <- 1 + 2 * selected
  spec$noise_diagnostics <- list(raw = if (is.null(noise)) NULL else noise$pooled,
    quiet_entries = sum(eligible), minimum_mean = cut)
  list(spec = spec, pilot = level)
}

