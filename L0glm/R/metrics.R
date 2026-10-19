#' @importFrom stats coef cor is.empty.model lm lm.fit logLik model.matrix model.response rnorm sd
#' @importFrom utils data
#' @importFrom qlcMatrix rowMax

# Essential fit, coefficient and support metrics used by tuning and benchmarks.
get.metrics <- function(beta, X, y, weights = NULL, offset = NULL,
                        family = gaussian(identity), beta_true = NULL, p = NULL,
                        dispersion = NULL, mu_floor = .Machine$double.eps,
                        always.include = NULL, r2_mu_floor = 0.1) {
  if (!inherits(family, "family")) family <- family()
  y <- as.matrix(y)
  n <- nrow(X); m <- ncol(y)
  if (nrow(y) != n || length(beta) != ncol(X) * m)
    stop("X, y and beta have incompatible dimensions.", call. = FALSE)
  beta <- L0glm_metric_coefficients(beta, ncol(X), m)
  if (length(mu_floor) != 1L || !is.finite(mu_floor) || mu_floor < 0)
    stop("mu_floor must be a finite nonnegative number.", call. = FALSE)
  weights <- L0glm_ext_matrix(weights, n, m, 1, "weights")
  offset <- L0glm_ext_matrix(offset, n, m, 0, "offset")
  if (any(!is.finite(weights) | weights < 0))
    stop("Scoring weights must be finite and nonnegative.", call. = FALSE)
  always.include <- L0glm_mandatory(X, always.include)
  if (is.null(p)) p <- ncol(X) - length(always.include)
  recovery <- L0glm_recovery_metrics(beta, beta_true, always.include, p)

  dispersion <- L0glm_dispersion(dispersion, family, m)
  free_scale <- if (is.null(dispersion)) rep(TRUE, m) else dispersion < 0
  scale_parameters <- if (family$family %in% c("gaussian", "Gamma", "inverse.gaussian",
    "quasibinomial", "quasipoisson")) sum(free_scale) else 0L
  edf <- sum(beta != 0) + scale_parameters + sum(beta[always.include, , drop = FALSE] == 0)

  predict_mean <- function(B) {
    active <- which(Matrix::rowSums(abs(B)) != 0)
    eta <- as.matrix(X[, active, drop = FALSE] %*% B[active, , drop = FALSE] + offset)
    matrix(family$linkinv(eta), n, m)
  }
  mu <- predict_mean(beta)
  prediction <- L0glm_prediction_metrics(y, mu, weights, family, beta, always.include,
    L0glm_metric_intercept(X, always.include), r2_mu_floor)
  if (family$family %in% c("poisson", "quasipoisson", "Gamma", "inverse.gaussian"))
    mu <- pmax(mu, mu_floor)
  keep <- which(weights > 0)
  deviance <- sum(family$dev.resids(as.vector(y), as.vector(mu), as.vector(weights))[keep])
  deviance_true <- NA_real_
  if (!is.null(beta_true)) {
    truth <- L0glm_metric_coefficients(beta_true, nrow(beta), m)
    mu_true <- predict_mean(truth)
    deviance_true <- sum(family$dev.resids(as.vector(mu_true), as.vector(mu), as.vector(weights))[keep])
  }

  # calculate some information criteria, cf. Zhang & Chen 2010 ####
  # see also https://github.com/abess-team/abess/blob/master/R-package/R/deviance.abess.R
  # get -2 * log-likelihood & log likelihood from family$aic slot & observed deviance
  # Score each response channel at its own scale. Positive supplied dispersions
  # are known; otherwise Gaussian likelihood profiles phi = weighted RSS / N.
  # family$aic includes an extra 2 for estimated scale in the continuous families.
  min2LL <- sum(vapply(seq_len(m), function(j) {
    keep_j <- weights[, j] > 0 & is.finite(weights[, j])
    yy <- y[keep_j, j]; mm <- mu[keep_j, j]; ww <- weights[keep_j, j]
    if (!length(yy)) return(0)
    phi <- if (is.null(dispersion) || free_scale[j]) NULL else dispersion[j]
    if (family$family == "gaussian") {
      rss <- sum(ww * (yy - mm)^2)
      if (is.null(phi)) phi <- rss / length(yy)
      if (phi == 0) return(-Inf) # unbounded likelihood at an exact interpolant
      return(length(yy) * log(2 * pi * phi) - sum(log(ww)) + rss / phi)
    }
    if (family$family == "poisson") yy <- round(yy)
    if (family$family == "Gamma" && !is.null(phi))
      return(-2 * sum(ww * stats::dgamma(yy, shape = 1 / phi, scale = mm * phi, log = TRUE)))
    if (family$family == "inverse.gaussian" && !is.null(phi))
      return(sum(ww * (log(2 * pi * phi * yy^3) + (yy - mm)^2 / (phi * yy * mm^2))))
    dev <- sum(family$dev.resids(yy, mm, ww))
    value <- family$aic(y = yy, n = ww, mu = mm, wt = ww, dev = dev)
    if (family$family %in% c("Gamma", "inverse.gaussian")) value <- value - 2
    value
  }, numeric(1)))
  loglik = -0.5 * min2LL # plain log-likelihood
  IC = get.IC(loglik = loglik, n = n, p = max(1L, p), edf = edf)
  
  c(list(p = p, n = n, m = m, df = edf,
    deviance = deviance, deviance_true = deviance_true, loglik = loglik), prediction, recovery,
    as.list(IC[c("gic", "ebic", "bic", "aic", "mbic")]))
}

# Response-scale prediction scores using working weights at the fitted means.
# Count active mean parameters, not dispersion/theta or selection/shrinkage df.
L0glm_prediction_metrics <- function(y, mu, weights, family, beta, mandatory=integer(),
                                     intercept=FALSE, r2_mu_floor=0.1) {
  if (length(r2_mu_floor) != 1L || !is.numeric(r2_mu_floor) ||
      !is.finite(r2_mu_floor) || r2_mu_floor < 0)
    stop("r2_mu_floor must be a finite nonnegative number.", call. = FALSE)
  result <- list(R2_weighted=NA_real_, adj_R2_weighted=NA_real_)
  m <- ncol(y)
  families <- if (inherits(family, "family")) rep(list(family), m) else family
  working <- matrix(0, nrow(y), m)
  for (j in seq_len(m)) {
    keep <- weights[, j] > 0
    f <- families[[j]]
    weight_mu <- mu[keep, j]
    count <- f$family %in% c("poisson", "quasipoisson") ||
      grepl("^negative binomial", f$family, ignore.case=TRUE)
    if (f$family %in% c("poisson", "quasipoisson", "Gamma", "inverse.gaussian") ||
        count)
      weight_mu <- pmax(weight_mu, .Machine$double.eps)
    # Boundary count means must not receive machine-epsilon inverse-variance
    # weights. Stabilize the weight calculation, never the fitted residuals.
    if (count) weight_mu <- pmax(weight_mu, r2_mu_floor)
    derivative <- f$mu.eta(f$linkfun(weight_mu))
    variance <- f$variance(weight_mu)
    value <- weights[keep, j] * derivative^2/variance
    value[variance == 0 & derivative == 0] <- 0
    if (any(!is.finite(value) | value < 0)) return(result)
    working[keep, j] <- value
  }
  observed <- colSums(working > 0)
  N <- sum(observed)
  if (!N) return(result)
  # A common rescaling leaves the score unchanged and avoids huge weights.
  working <- working/max(working)
  active <- observed > 0
  keep <- working > 0
  response <- y; response[!keep] <- 0
  center <- if (intercept) colSums(working*response)/pmax(colSums(working), .Machine$double.xmin) else rep(0, m)
  rss <- sum(working[keep] * (y[keep] - mu[keep])^2)
  tss <- sum(working[keep] * sweep(y, 2L, center)[keep]^2)
  if (!is.finite(rss) || !is.finite(tss) || tss <= 0) return(result)
  result$R2_weighted <- 1 - rss/tss
  df <- as.numeric(Matrix::colSums(beta != 0) +
    Matrix::colSums(beta[mandatory, , drop=FALSE] == 0))
  k <- sum(df[active]); baseline <- as.integer(intercept) * sum(active)
  if (N > k && N > baseline)
    result$adj_R2_weighted <- 1 - (rss/tss) * (N - baseline)/(N - k)
  result
}

L0glm_metric_intercept <- function(X, mandatory) {
  any(vapply(mandatory, function(j) all(X[, j, drop=FALSE] == 1), logical(1)))
}

L0glm_metric_coefficients <- function(beta, p, m) {
  if (is.null(dim(beta))) beta <- matrix(beta, p, m)
  Matrix::Matrix(beta, sparse = TRUE)
}

# Binary support counts; undefined ratios and constant-margin MCC return NA.
support_scores <- function(TP, FP, TN, FN) {
  ratio <- function(x, denominator) if (denominator > 0) x / denominator else NA_real_
  list(TP = TP, FP = FP, TN = TN, FN = FN,
    precision = ratio(TP, TP + FP), recall = ratio(TP, TP + FN),
    F1_score = ratio(2 * TP, 2 * TP + FP + FN), MCC = support_mcc(TP, FP, TN, FN))
}

L0glm_recovery_metrics <- function(beta, beta_true = NULL, always.include = integer(),
                                    p = nrow(beta) - length(always.include)) {
  if (!is.null(beta_true) && length(beta_true) != length(beta))
    stop("beta_true must have the same dimensions as the fitted coefficients.", call. = FALSE)
  if (!is.null(dim(beta_true)) && !identical(dim(beta_true), dim(beta)))
    stop("beta_true must have the same dimensions as the fitted coefficients.", call. = FALSE)
  estimated <- setdiff(which(Matrix::rowSums(abs(beta)) != 0), always.include)
  counts <- as.list(stats::setNames(rep(NA_real_, 8L),
    c("TP", "FP", "TN", "FN", "precision", "recall", "F1_score", "MCC")))
  within <- counts
  names(within) <- paste0("within_row_", names(within))
  support_size_true <- MAD_coefs <- RMSE_coefs <- R2_coefs <- NA_real_
  if (!is.null(beta_true)) {
    truth <- L0glm_metric_coefficients(beta_true, nrow(beta), ncol(beta))
    true_rows <- setdiff(which(Matrix::rowSums(abs(truth)) != 0), always.include)
    recovered <- intersect(estimated, true_rows)
    support_size_true <- length(true_rows)
    counts <- support_scores(length(recovered), length(setdiff(estimated, true_rows)),
      p - length(union(estimated, true_rows)), length(setdiff(true_rows, estimated)))
    error <- beta - truth
    MAD_coefs <- sum(abs(error)) / length(beta)
    RMSE_coefs <- sqrt(sum(error^2) / length(beta))
    ss_estimated <- sum(beta^2); ss_true <- sum(truth^2)
    if (ss_true > 0) R2_coefs <- if (ss_estimated > 0)
      min(1, sum(beta * truth)^2 / (ss_estimated * ss_true)) else 0
    if (ncol(beta) > 1L) {
      selected <- beta[recovered, , drop = FALSE] != 0
      present <- truth[recovered, , drop = FALSE] != 0
      TP <- sum(selected & present)
      FP <- sum(selected) - TP; FN <- sum(present) - TP
      within <- support_scores(TP, FP, length(selected) - TP - FP - FN, FN)
      names(within) <- paste0("within_row_", names(within))
    }
  }
  c(list(MAD_coefs = MAD_coefs, RMSE_coefs = RMSE_coefs, R2_coefs = R2_coefs,
      support_size_true = support_size_true, support_size_estimated = length(estimated)),
    counts, within)
}

# Matthews correlation coefficient (Pearson's phi) of binary row support.
# Cast before multiplication to avoid integer overflow on large designs.
# A constant truth or prediction has undefined correlation: return NA.
support_mcc <- function(TP, FP, TN, FN) {
  counts <- as.double(c(TP, FP, TN, FN))
  TP <- counts[1L]; FP <- counts[2L]; TN <- counts[3L]; FN <- counts[4L]
  margins <- c(TP + FP, TP + FN, TN + FP, TN + FN)
  if (any(!is.finite(margins)) || any(margins <= 0)) return(NA_real_)
  (TP * TN - FP * FN) / prod(sqrt(margins))
}


# return information criteria based on log likelihood, n, p & edf
# for a good discussion on information criteria see
# https://www.sciencedirect.com/science/article/pii/S016794731100171X
get.IC = function(loglik, n, p, edf, 
                  gamma_ebic = max(0, 1 - log(n)/(2*log(p))), # reduces to BIC when p <= sqrt(n)
                  Cn_hbic = log(log(n))) {
  min2LL = -2 * loglik
  
  # p=nr of vars, n=nr of cases, edf=effective degrees of freedom
  aic = min2LL + 2 * edf # Akaike information criterion
  # aicc = ifelse(n > edf, min2LL + 2*edf*n/(n-edf), NA) # small sample AIC
  bic = min2LL + log(n) * edf # Bayesian information criterion
  gic = min2LL + log(p) * log(log(n)) * edf # generalized information criterion GIC = SIC in https://www.pnas.org/doi/10.1073/pnas.2014241117
   
  ebic = min2LL + (log(n) + 2 * gamma_ebic * log(p)) * edf # linearized EBIC; Chen and Chen (2008), Biometrika 95, 759-771
  hbic = min2LL + Cn_hbic * log(p) * edf # HBIC, recommended for ultra-high dimensional settings, cf https://arxiv.org/pdf/1311.4981 & https://www.annualreviews.org/docserver/fulltext/statistics/7/1/annurev-statistics-030718-105038.pdf
  mbic = min2LL + log(n * (p ^ 2) / 16) * edf # multiplicity-adjusted mBIC; Bogdan et al. (2008), doi:10.1002/qre.936; distinct from mBIC2
  # TODO maybe still add mbic2, cf Frommlet paper
  return(data.frame(aic = aic, bic = bic, ebic = ebic, hbic = hbic, gic = gic, mbic = mbic))
}


# Shared, case-insensitive names for fitting, scoring and tuning interfaces.
L0glm_ic_name <- function(IC, allow_cv = FALSE) {
  choices <- c("aic", "bic", "gic", "ebic", "mbic")
  if (allow_cv) choices <- c(choices, "cv")
  if (!is.character(IC) || length(IC) != 1L || is.na(IC) ||
      !tolower(IC) %in% choices) {
    stop(if (allow_cv) "tune must be one of " else "Information criterion must be one of ",
         paste(shQuote(choices), collapse = ", "),
         " (case-insensitive).", call. = FALSE)
  }
  tolower(IC)
}

L0glm_resolve_lambda <- function(lambda, p, n) {
  if (is.character(lambda)) lambda <- preset.lambda(lambda, p = p, n = n)
  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda)) {
    stop("lambda must be a finite numeric scalar or an information-criterion name; ",
         "see help('information_criteria', package = 'L0glm').", call. = FALSE)
  }
  lambda
}


# preset.lambda ####
# Per-coefficient penalty on the -2*logLik scale; see information_criteria.
preset.lambda <- function(IC="gic", p, n,
                          gamma_ebic = max(0, 1 - log(n)/(2*log(p)))) { # reduces to BIC when p <= sqrt(n)
  IC <- L0glm_ic_name(IC)
  lambda <- switch(IC,
                   gic = log(p)*log(log(n)),
                   aic = 2,
                   bic = log(n),
                   ebic = log(n) + 2*gamma_ebic*log(p),
                   mbic = log(n*(p^2)/16))
  return(lambda)
}
