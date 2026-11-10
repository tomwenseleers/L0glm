# R counterpart to scale_step.h for callbacks, validation and extended families.
L0glm_scale_estimable <- function(n, k, stable, iteration, relaxed) {
  n > k && (relaxed || (stable < 3L && iteration >= 3L && k < .7*n))
}

L0glm_scale_step <- function(value, target, stable, changed, relaxed) {
  if (!is.finite(target) || target <= 0)
    return(list(value = value, stable = stable, updated = FALSE))
  stable <- if (!changed && abs(target/value-1) < .01) stable+1L else 0L
  list(value = if (relaxed) target else min(2*value, max(value/2, target)),
       stable = stable, updated = TRUE)
}

L0glm_pearson_step <- function(phi, stable, y, mu, w, family, k, iteration, changed, relaxed) {
  if (changed) stable <- 0L
  target <- if (L0glm_scale_estimable(sum(w>0), k, stable, iteration, relaxed))
    L0glm_pearson_scale(y, mu, w, family, k) else NA_real_
  L0glm_scale_step(phi, target, stable, changed, relaxed)
}

L0glm_nb_pearson_shape <- function(y, mu, w, k, control) {
  good <- w > 0
  y <- y[good]; mu <- mu[good]; w <- w[good]
  if (length(y) <= k || any(!is.finite(mu) | mu <= 0)) return(NA_real_)
  residual <- w*(y-mu)^2
  mean_squared <- mu^2
  target <- mean(w)*(length(y)-k)
  equation <- function(log_theta) sum(residual/(mu+mean_squared/exp(log_theta)))/target-1
  limits <- log(c(control$theta.min, control$theta.max))
  lo <- equation(limits[1L]); hi <- equation(limits[2L])
  if (is.na(lo) || is.na(hi)) return(NA_real_)
  if (hi <= 0) return(control$theta.max)
  if (lo >= 0) return(control$theta.min)
  exp(stats::uniroot(equation, limits, tol = 1e-7)$root)
}

L0glm_nb_scale_step <- function(data, B, previous, state, iteration, relaxed) {
  if (!any(state$estimated) || isTRUE(data$fixed.dispersion)) return(state)
  if (is.null(state$stable)) state$stable <- integer(data$m)
  mu <- L0glm_ext_nb_mean(data, B)$mu
  for (j in which(state$estimated)) {
    changed <- any((B[,j]!=0) != (previous[,j]!=0))
    stable <- if (changed) 0L else state$stable[j]
    w <- as.numeric(L0glm_observation_column(data$weights, j))
    k <- if (is.null(data$noise.df)) sum(B[,j]!=0 | seq_len(nrow(B)) %in% data$included) else data$noise.df[j]
    target <- if (L0glm_scale_estimable(sum(w>0), k, stable, iteration, relaxed))
      L0glm_nb_pearson_shape(data$y[,j], mu[,j], w, k, data$control) else NA_real_
    step <- L0glm_scale_step(state$theta[j], target, stable, changed, relaxed)
    state$theta[j] <- step$value
    state$stable[j] <- step$stable
  }
  state
}
