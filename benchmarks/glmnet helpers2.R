# ============== helpers (unchanged logic for -2 log L on ML scale) ==============

.min2LL_from_family <- function(y, mu, wt, fam, dev) {
  eps <- .Machine$double.eps
  if (fam$family %in% c("poisson","quasipoisson","Gamma","inverse.gaussian")) mu <- pmax(mu, eps)
  if (fam$family == "binomial") mu <- pmin(pmax(mu, eps), 1 - eps)
  n_arg <- if (fam$family == "binomial") wt else rep.int(1, length(wt))
  fam$aic(y = y, n = n_arg, mu = mu, wt = wt, dev = dev)  # = -2 log L
}

# Unpenalised base-R GLM refit
.refit_glm <- function(Xs, y, wt, offset, fam, intercept) {
  form <- if (intercept) y ~ Xs else y ~ Xs - 1
  fit  <- glm(form, weights = wt, offset = offset, family = fam,
              x = FALSE, y = FALSE, model = FALSE)
  mm   <- model.matrix(form)
  rk   <- qr(mm)$rank - (if (intercept) 1L else 0L)
  mu   <- fit$fitted.values
  dev  <- sum(fam$dev.resids(y = y, mu = mu, wt = wt))
  m2LL <- .min2LL_from_family(y, mu, wt, fam, dev)
  cf   <- coef(fit); beta0 <- if (intercept) unname(cf[1L]) else 0
  bs   <- if (intercept) unname(cf[-1L]) else unname(cf)
  list(min2LL = m2LL, rank = rk, beta0 = beta0, beta = bs)
}

# Tiny-ridge refit via glmnet (exact ridge for GLMs; supports sparse X and bounds)
.refit_ridge_glmnet <- function(Xs, y, wt, offset, family_name,
                                intercept, lower.limits, upper.limits, ridge_lambda = 1e-6) {
  fit <- glmnet::glmnet(x = Xs, y = y, weights = wt, offset = offset,
                        family = family_name, alpha = 0, lambda = ridge_lambda,
                        intercept = intercept, lower.limits = lower.limits,
                        upper.limits = upper.limits, standardize = FALSE)
  b     <- as.numeric(coef(fit, s = ridge_lambda))
  beta0 <- if (intercept) b[1L] else 0
  beta  <- if (intercept) b[-1L] else b[-1L]
  fam   <- switch(family_name, gaussian=gaussian(), binomial=binomial(), poisson=poisson())
  eta   <- (if (intercept) beta0 else 0) + drop(Xs %*% beta) + if (is.null(offset)) 0 else offset
  mu    <- fam$linkinv(eta)
  if (is.null(wt)) wt <- rep(1, length(y))
  dev   <- sum(fam$dev.resids(y = y, mu = mu, wt = wt))
  m2LL  <- .min2LL_from_family(y, mu, wt, fam, dev)
  rk    <- qr(Xs)$rank  # numeric rank of the design (safest for df)
  list(min2LL = m2LL, rank = rk, beta0 = beta0, beta = beta)
}

# NEW: Unpenalised refit via glmnet::bigGlm (inherits glmnet; supports sparse X, bounds, offsets)
.refit_bigGlm <- function(Xs, y, wt, offset, family_name,
                          intercept, lower.limits, upper.limits, path = FALSE) {
  # family_name: "gaussian" | "binomial" | "poisson"
  # Most glmnet args “that make sense” are supported, including bounds and offsets.
  fit <- glmnet::bigGlm(x = Xs, y = y, weights = wt, offset = offset,
                        family = family_name, intercept = intercept,
                        lower.limits = lower.limits, upper.limits = upper.limits,
                        path = path, standardize = FALSE)
  b     <- as.numeric(coef(fit))   # for bigGlm, coef() gives the λ=0 solution
  beta0 <- if (intercept) b[1L] else 0
  beta  <- if (intercept) b[-1L] else b[-1L]
  fam   <- switch(family_name, gaussian=gaussian(), binomial=binomial(), poisson=poisson())
  eta   <- (if (intercept) beta0 else 0) + drop(Xs %*% beta) + if (is.null(offset)) 0 else offset
  mu    <- fam$linkinv(eta)
  if (is.null(wt)) wt <- rep(1, length(y))
  dev   <- sum(fam$dev.resids(y = y, mu = mu, wt = wt))
  m2LL  <- .min2LL_from_family(y, mu, wt, fam, dev)
  rk    <- qr(Xs)$rank
  list(min2LL = m2LL, rank = rk, beta0 = beta0, beta = beta)
}

.IC_values <- function(min2LL, k, n, p, gamma_ebic, Cn_hbic) {
  c(
    aic  = min2LL + 2 * k,
    bic  = min2LL + log(n) * k,
    ebic = min2LL + (log(n) + 2 * gamma_ebic * log(p)) * k,
    hbic = min2LL + (Cn_hbic * log(p)) * k,                 # deviance-scale HBIC
    gic  = min2LL + (log(p) * log(log(n))) * k,
    mbic = min2LL + log(n * p^2 / 16) * k
  )
}

# ============================= main scorer ==================================

score_ic_glmnet <- function(fit, X, y,
                            family = gaussian(),             # base-R family object
                            family_name = c("gaussian","binomial","poisson")[1],
                            weights = NULL, offset = NULL,
                            intercept = FALSE,
                            refit_method = c("glm","bigGlm","ridge_glmnet"),
                            bigGlm_path = FALSE,             # use path warm-start in bigGlm
                            ridge_lambda = 1e-6,             # for ridge_glmnet
                            lower.limits = -Inf, upper.limits = Inf,
                            p_total = ncol(X),
                            gamma_ebic = NULL,
                            Cn_hbic = log(log(NROW(y))),
                            min_resid_df = 3L,
                            return_full_coef = TRUE,
                            fallback_to_ridge = TRUE) {      # try ridge if bigGlm/glm fails
  n <- NROW(y); p <- p_total
  if (is.null(weights)) weights <- rep(1, n)
  if (is.null(offset))  offset  <- rep(0, n)
  if (is.null(gamma_ebic)) gamma_ebic <- max(0, 1 - log(n)/(2*log(p)))
  refit_method <- match.arg(refit_method)
  
  B <- fit$beta
  L <- if (inherits(B, "dgCMatrix")) ncol(B) else NCOL(B)
  
  IC      <- matrix(NA_real_, L, 6, dimnames = list(NULL, c("aic","bic","ebic","hbic","gic","mbic")))
  kvec    <- integer(L)
  rankv   <- integer(L)
  ok      <- logical(L)
  active  <- vector("list", L)
  beta0   <- numeric(L)
  coef_full <- if (return_full_coef) vector("list", L) else NULL
  
  for (j in seq_len(L)) {
    Sj <- which(B[, j] != 0)
    active[[j]] <- Sj
    if (length(Sj) == 0L) {
      # Null model likelihood
      mu0  <- if (intercept) rep(mean(y), n) else rep(0, n)
      dev0 <- sum(family$dev.resids(y = y, mu = mu0, wt = weights))
      m2LL <- .min2LL_from_family(y, mu0, weights, family, dev0)
      k <- 0L; rk <- 0L; okj <- TRUE; b0 <- if (intercept) mean(y) else 0; bs <- numeric(0)
    } else {
      Xs <- X[, Sj, drop = FALSE]
      tmp <- try(
        switch(refit_method,
               glm          = .refit_glm(Xs, y, weights, offset, family, intercept),
               bigGlm       = .refit_bigGlm(Xs, y, weights, offset, family_name,
                                            intercept, lower.limits, upper.limits,
                                            path = bigGlm_path),
               ridge_glmnet = .refit_ridge_glmnet(Xs, y, weights, offset, family_name,
                                                  intercept, lower.limits, upper.limits, ridge_lambda)),
        silent = TRUE
      )
      if (inherits(tmp, "try-error") && fallback_to_ridge) {
        tmp <- .refit_ridge_glmnet(Xs, y, weights, offset, family_name,
                                   intercept, lower.limits, upper.limits, ridge_lambda)
      }
      rk <- tmp$rank
      resid_df <- n - (rk + as.integer(intercept))
      if (is.na(resid_df) || resid_df < min_resid_df) {
        m2LL <- Inf; k <- rk; okj <- FALSE; b0 <- NA_real_; bs <- rep(NA_real_, length(Sj))
      } else {
        m2LL <- tmp$min2LL; k <- rk; okj <- TRUE; b0 <- tmp$beta0; bs <- tmp$beta
      }
    }
    beta0[j] <- b0; rankv[j] <- rk; kvec[j] <- k; ok[j] <- okj
    IC[j, ]  <- .IC_values(m2LL, k, n, p, gamma_ebic, Cn_hbic)
    if (return_full_coef) {
      if (!okj || length(Sj) == 0L) {
        coef_full[[j]] <- Matrix::sparseVector(i = integer(0), x = numeric(0), length = ncol(X))
      } else {
        coef_full[[j]] <- Matrix::sparseVector(i = Sj, x = bs, length = ncol(X))
      }
    }
  }
  
  ans <- data.frame(lambda = fit$lambda, k = kvec, rank = rankv, ok = ok,
                    IC, beta0 = beta0, active = I(active), check.names = FALSE)
  if (return_full_coef) ans$coef_full <- I(coef_full)
  ans
}