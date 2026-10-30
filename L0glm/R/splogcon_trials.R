# Joint local fusion refinement: the candidate profile and overlapping spectra
# can adjust together, while every outside contribution stays in the offset.
slcn_enhance_trials <- function(engine, ctx) {
  original_trial <- engine$spolo_trial
  original_coherent <- engine$spolo_coherent
  engine$spolo_coherent <- function(model, d, ncores = 8L, sweeps = 8L, background = FALSE) {
    proposal <- original_coherent(model, d, ncores, sweeps, background)
    delta <- engine$spolo_delta_cpp(d$Y, model$mu, proposal$mu) + d$lambda_rank * (ncol(proposal$H) - ncol(model$H))
    if (is.finite(delta) && delta <= 1e-6) return(proposal)
    model$seconds <- proposal$seconds;model$coherent <- c(proposal$coherent, list(rejected = TRUE))
    model
  }
  engine$spolo_trial <- function(st, proposal, sweeps = 4L, ncores = 8L) {
    donors <- proposal$group
    left <- min(st$lo[donors], proposal$left);right <- max(st$hi[donors], proposal$right)
    live <- which(colSums(st$T) > 0)
    neighbours <- setdiff(live[st$lo[live] <= right & st$hi[live] >= left], donors)
    if (length(neighbours) > 16L) {
      overlap <- pmax(0, pmin(st$hi[neighbours], right) - pmax(st$lo[neighbours], left) + 1L)
      neighbours <- neighbours[order(overlap * sqrt(colSums(st$T)[neighbours]), decreasing = TRUE)[1:16]]
    }
    expanded <- proposal
    expanded$group <- c(donors, neighbours)
    expanded$E <- cbind(proposal$E, matrix(0, nrow(st$T), length(neighbours)))
    complete <- expanded
    complete$E[, seq_along(donors)] <- st$T[, donors, drop = FALSE]
    complete$kind <- paste0(proposal$kind, "_full")
    best <- NULL
    variants <- if (identical(expanded$E, complete$E)) list(expanded) else list(expanded, complete)
    for (variant in variants) {
      trial <- original_trial(st, variant, sweeps, ncores)
      if (is.null(trial)) next
      a <- engine$spolo_local_cpp(st$handle, trial$left, trial$right)
      focus <- ncol(trial$H)
      for (cycle in seq_len(if (ctx$spec$kind == 0L) 8L else 2L)) {
        rows <- which(trial$H[, focus] > 0);rows <- min(rows):max(rows)
        old_h <- trial$H[rows, focus];old_s <- trial$S[, focus]
        raw <- slcn_profile(a$Y[rows, , drop = FALSE], trial$mu[rows, , drop = FALSE],
          old_h, old_s, ctx$spec$kind, ctx$spec$dispersion, ctx$spec$theta, ctx$W,
          ctx$n, trial$left + min(rows) - 2L, 10L, as.integer(ncores), ctx$censor)
        mass <- sum(raw);if (!is.finite(mass) || mass <= 0) break
        proposed_h <- engine$spolo_lc(raw)$h;proposed_s <- old_s * mass
        rest <- trial$mu[rows, , drop = FALSE] - old_h %o% old_s
        accepted <- FALSE
        for (fraction in c(1, .5, .25)) {
          h <- engine$spolo_lc((1 - fraction) * old_h + fraction * proposed_h)$h
          s <- (1 - fraction) * old_s + fraction * proposed_s
          H <- trial$H;S <- trial$S;mu <- trial$mu
          H[rows, focus] <- h;S[, focus] <- s;mu[rows, ] <- pmax(rest, 0) + h %o% s
          active <- which(colSums(S) > 0)
          fitted <- engine$spolo_spectral_cpp(a$Y, mu, H[, active, drop = FALSE],
            S[, active, drop = FALSE], 3L, ncores)
          S[] <- 0;S[, active] <- fitted$S
          delta <- engine$spolo_delta_cpp(a$Y, trial$mu, fitted$mu) +
            st$lambda * (sum(colSums(S) > 0) - sum(colSums(trial$S) > 0))
          if (is.finite(delta) && delta < -1e-7 && sum(S[, focus]) > 0) {
            trial$H <- H;trial$S <- S;trial$mu <- fitted$mu;accepted <- TRUE;break
          }
        }
        if (!accepted) break
      }
      # Reconsider residual donors after joint reallocation of their signal.
      for (j in order(colSums(trial$S))) {
        if (j == focus || sum(trial$S[, j]) <= 0) next
        mu <- pmax(trial$mu - trial$H[, j] %o% trial$S[, j], 0)
        delta <- engine$spolo_delta_cpp(a$Y, trial$mu, mu) - st$lambda
        if (is.finite(delta) && delta < 0) {trial$S[, j] <- 0;trial$mu <- mu}
      }
      trial$delta <- engine$spolo_delta_cpp(a$Y, a$mu, trial$mu) +
        st$lambda * (sum(colSums(trial$S) > 0) - length(trial$group))
      if (is.finite(trial$delta) && (is.null(best) || trial$delta < best$delta)) best <- trial
    }
    best
  }
  invisible(engine)
}
