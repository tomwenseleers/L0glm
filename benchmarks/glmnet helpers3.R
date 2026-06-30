# ===================== main scorer with optional TRIMMING ======================

score_ic_glmnet <- function(fit, X, y,
                            family = gaussian(), family_name = c("gaussian","binomial","poisson")[1],
                            weights = NULL, offset = NULL,
                            intercept = FALSE,
                            refit_method = c("glm","bigGlm","ridge_glmnet"),
                            bigGlm_path = TRUE, ridge_lambda = 1e-6,
                            lower.limits = -Inf, upper.limits = Inf,
                            p_total = ncol(X),
                            gamma_ebic = NULL, Cn_hbic = log(log(NROW(y))),
                            min_resid_df = 3L,
                            # --- trimming controls ---
                            trim = TRUE,                        # turn trimming on/off
                            trim_choose = c("ebic","bic","hbic","aic","gic","mbic")[1],
                            trim_exact_refit = FALSE,           # FALSE: approximate (no second refit); TRUE: exact refit
                            t_cap = NULL,                       # optional hard cap on t; default uses residual-df guard
                            return_full_coef = TRUE,
                            fallback_to_ridge = TRUE) {
  
  n <- NROW(y); p <- p_total
  if (is.null(weights)) weights <- rep(1, n)
  if (is.null(offset))  offset  <- rep(0, n)
  if (is.null(gamma_ebic)) gamma_ebic <- max(0, 1 - log(n)/(2*log(p)))
  refit_method <- match.arg(refit_method)
  # trim_choose  <- match.arg(trim_choose)
  
  B <- fit$beta
  L <- if (inherits(B, "dgCMatrix")) ncol(B) else NCOL(B)
  
  IC_raw   <- matrix(NA_real_, L, 6, dimnames=list(NULL,c("aic","bic","ebic","hbic","gic","mbic")))
  k_raw    <- integer(L); rank_raw <- integer(L); ok_raw <- logical(L)
  active   <- vector("list", L)
  beta0_raw <- numeric(L)
  coef_full_raw <- if (return_full_coef) vector("list", L) else NULL
  
  # trimmed outputs
  IC_trim   <- matrix(NA_real_, L, 6, dimnames=dimnames(IC_raw))
  k_trim    <- integer(L); rank_trim <- integer(L); ok_trim <- logical(L)
  t_star    <- integer(L)
  beta0_trim <- numeric(L)
  coef_full_trim <- if (return_full_coef) vector("list", L) else NULL
  
  for (j in seq_len(L)) {
    Sj <- which(B[, j] != 0)
    active[[j]] <- Sj
    
    # ---------- raw-support refit ----------
    if (length(Sj) == 0L) {
      mu0  <- if (intercept) rep(mean(y), n) else rep(0, n)
      dev0 <- sum(family$dev.resids(y = y, mu = mu0, wt = weights))
      m2LL0 <- .min2LL_from_family(y, mu0, weights, family, dev0)
      k_raw[j] <- 0L; rank_raw[j] <- 0L; ok_raw[j] <- TRUE
      IC_raw[j,] <- .IC_values(m2LL0, 0L, n, p, gamma_ebic, Cn_hbic)
      beta0_raw[j] <- if (intercept) mean(y) else 0
      if (return_full_coef) coef_full_raw[[j]] <- Matrix::sparseVector(i=integer(0), x=numeric(0), length=ncol(X))
      # nothing to trim
      ok_trim[j] <- TRUE; t_star[j] <- 0L
      IC_trim[j,] <- IC_raw[j,]; k_trim[j] <- 0L; rank_trim[j] <- 0L
      beta0_trim[j] <- beta0_raw[j]
      if (return_full_coef) coef_full_trim[[j]] <- coef_full_raw[[j]]
      next
    }
    
    Xs <- X[, Sj, drop=FALSE]
    
    # one refit on the full active set (raw)
    tmp_raw <- try(
      switch(refit_method,
             glm          = .refit_glm(Xs, y, weights, offset, family, intercept),
             bigGlm       = .refit_bigGlm(Xs, y, weights, offset, family_name, intercept,
                                          lower.limits, upper.limits, path = bigGlm_path),
             ridge_glmnet = .refit_ridge_glmnet(Xs, y, weights, offset, family_name,
                                                intercept, lower.limits, upper.limits, ridge_lambda)),
      silent = TRUE)
    
    if (inherits(tmp_raw,"try-error") && fallback_to_ridge) {
      tmp_raw <- .refit_ridge_glmnet(Xs, y, weights, offset, family_name,
                                     intercept, lower.limits, upper.limits, ridge_lambda)
    }
    
    rk_raw <- tmp_raw$rank
    resid_df_raw <- n - (rk_raw + as.integer(intercept))
    if (is.na(resid_df_raw) || resid_df_raw < min_resid_df) {
      # invalidate raw model
      IC_raw[j,] <- .IC_values(Inf, rk_raw, n, p, gamma_ebic, Cn_hbic)
      k_raw[j] <- rk_raw; rank_raw[j] <- rk_raw; ok_raw[j] <- FALSE
      beta0_raw[j] <- NA_real_
      if (return_full_coef) coef_full_raw[[j]] <- Matrix::sparseVector(i=integer(0), x=numeric(0), length=ncol(X))
      # still attempt trimming (it might reduce rank)
      beta_full <- rep(NA_real_, length(Sj))
    } else {
      IC_raw[j,] <- .IC_values(tmp_raw$min2LL, rk_raw, n, p, gamma_ebic, Cn_hbic)
      k_raw[j] <- rk_raw; rank_raw[j] <- rk_raw; ok_raw[j] <- TRUE
      beta0_raw[j] <- tmp_raw$beta0
      beta_full <- tmp_raw$beta
      if (return_full_coef) {
        coef_full_raw[[j]] <- Matrix::sparseVector(i = Sj, x = beta_full, length = ncol(X))
      }
    }
    
    # ---------- trimming over t ----------
    if (!trim || all(!is.finite(beta_full))) {
      # just copy raw into trimmed
      ok_trim[j] <- ok_raw[j]; t_star[j] <- length(Sj)
      IC_trim[j,] <- IC_raw[j,]; k_trim[j] <- k_raw[j]; rank_trim[j] <- rank_raw[j]
      beta0_trim[j] <- beta0_raw[j]
      if (return_full_coef) coef_full_trim[[j]] <- coef_full_raw[[j]]
      next
    }
    
    ord <- order(abs(beta_full), decreasing=TRUE)
    # cap t by residual-df guard and any user cap
    t_guard <- max(0L, n - as.integer(intercept) - min_resid_df)
    t_max   <- min(length(Sj), t_guard, if (is.null(t_cap)) length(Sj) else t_cap)
    
    best <- list(crit = Inf, j = NA_integer_, IC = rep(NA_real_,6),
                 k = NA_integer_, rk = NA_integer_,
                 beta0 = NA_real_, coef = NULL, idx = NULL)
    
    for (t in 0:t_max) {
      if (t == 0L) {
        # null under trimming
        mu0  <- if (intercept) rep(mean(y), n) else rep(0, n)
        dev0 <- sum(family$dev.resids(y = y, mu = mu0, wt = weights))
        m2LL <- .min2LL_from_family(y, mu0, weights, family, dev0)
        rk   <- 0L; k    <- 0L
        ICv  <- .IC_values(m2LL, k, n, p, gamma_ebic, Cn_hbic)
        crit <- ICv[[trim_choose]]
        beta0_t <- if (intercept) mean(y) else 0
        coef_t  <- Matrix::sparseVector(i = integer(0), x = numeric(0), length = ncol(X))
      } else {
        keep <- ord[seq_len(t)]
        Xst  <- Xs[, keep, drop=FALSE]
        if (!trim_exact_refit) {
          # approximate: zero small coefs, reuse beta0_raw
          b_keep <- beta_full[keep]
          eta    <- (if (intercept) beta0_raw[j] else 0) + drop(Xst %*% b_keep) + offset
          mu     <- family$linkinv(eta)
          dev    <- sum(family$dev.resids(y = y, mu = mu, wt = weights))
          m2LL   <- .min2LL_from_family(y, mu, weights, family, dev)
          rk     <- qr(Xst)$rank
          k      <- rk
          ICv    <- .IC_values(m2LL, k, n, p, gamma_ebic, Cn_hbic)
          crit   <- ICv[[trim_choose]]
          beta0_t <- beta0_raw[j]
          coef_t  <- Matrix::sparseVector(i = Sj[keep], x = b_keep, length = ncol(X))
        } else {
          # exact: refit on top-t via chosen engine
          tmp_t <- try(
            switch(refit_method,
                   glm          = .refit_glm(Xst, y, weights, offset, family, intercept),
                   bigGlm       = .refit_bigGlm(Xst, y, weights, offset, family_name, intercept,
                                                lower.limits, upper.limits, path = bigGlm_path),
                   ridge_glmnet = .refit_ridge_glmnet(Xst, y, weights, offset, family_name,
                                                      intercept, lower.limits, upper.limits, ridge_lambda)),
            silent = TRUE)
          if (inherits(tmp_t,"try-error") && fallback_to_ridge) {
            tmp_t <- .refit_ridge_glmnet(Xst, y, weights, offset, family_name,
                                         intercept, lower.limits, upper.limits, ridge_lambda)
          }
          rk     <- tmp_t$rank
          # guard
          if (n - (rk + as.integer(intercept)) < min_resid_df) next
          k      <- rk
          ICv    <- .IC_values(tmp_t$min2LL, k, n, p, gamma_ebic, Cn_hbic)
          crit   <- ICv[[trim_choose]]
          beta0_t <- tmp_t$beta0
          coef_t  <- Matrix::sparseVector(i = Sj[keep], x = tmp_t$beta, length = ncol(X))
        }
      }
      if (crit < best$crit) best <- list(crit=crit, j=t, IC=ICv, k=k, rk=rk,
                                         beta0=beta0_t, coef=coef_t, idx = if (t==0) integer(0) else Sj[keep])
    }
    
    # record trimmed winner
    t_star[j]     <- best$j
    IC_trim[j, ]  <- best$IC
    k_trim[j]     <- best$k
    rank_trim[j]  <- best$rk
    ok_trim[j]    <- is.finite(best$crit)
    beta0_trim[j] <- best$beta0
    if (return_full_coef) coef_full_trim[[j]] <- best$coef
  }
  
  # assemble result
  out <- data.frame(
    lambda = fit$lambda,
    # raw-support refit summary
    k_raw = k_raw, rank_raw = rank_raw, ok_raw = ok_raw,
    aic_raw = IC_raw[, "aic"], bic_raw = IC_raw[, "bic"], ebic_raw = IC_raw[, "ebic"],
    hbic_raw = IC_raw[, "hbic"], gic_raw = IC_raw[, "gic"], mbic_raw = IC_raw[, "mbic"],
    beta0_raw = beta0_raw,
    # trimmed summary
    t_star = t_star, k_trim = k_trim, rank_trim = rank_trim, ok_trim = ok_trim,
    aic_trim = IC_trim[, "aic"], bic_trim = IC_trim[, "bic"], ebic_trim = IC_trim[, "ebic"],
    hbic_trim = IC_trim[, "hbic"], gic_trim = IC_trim[, "gic"], mbic_trim = IC_trim[, "mbic"],
    beta0_trim = beta0_trim,
    active = I(active),
    check.names = FALSE
  )
  if (return_full_coef) {
    out$coef_full_raw  <- I(coef_full_raw)
    out$coef_full_trim <- I(coef_full_trim)
  }
  out
}
