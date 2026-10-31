# Reunite complementary ion sets when their broad elution profiles coincide.
# Candidate construction is algebraic; only the original-channel loss is scored.
slcn_fuse_profiles <- function(model, Y, psf_width) {
  model$members <- lapply(seq_len(ncol(model$H)), identity)
  audit <- list()
  repeat {
    H <- model$H;K <- ncol(H)
    if (K < 2L) break
    width <- vapply(seq_len(K), function(k) {
      h <- H@x[seq.int(H@p[k] + 1L, H@p[k + 1L])]
      sum(h > max(h) / 2)
    }, 0L)
    broad <- which(width >= 4 * psf_width)
    if (!length(broad)) break
    norms <- colSums(H^2)
    cosine <- as.matrix(crossprod(H[, broad, drop = FALSE], H)) /
      sqrt(norms[broad] %o% norms)
    pairs <- which(cosine >= .95, arr.ind = TRUE)
    if (!nrow(pairs)) break
    a <- broad[pairs[, 1L]];b <- pairs[, 2L]
    pairs <- unique(cbind(pmin(a, b), pmax(a, b)))
    pairs <- pairs[pairs[, 1L] < pairs[, 2L], , drop = FALSE]
    pairs <- pairs[order(pairs[, 1L], pairs[, 2L]), , drop = FALSE]
    best <- NULL;best_delta <- -1e-6
    for (i in seq_len(nrow(pairs))) {
      group <- pairs[i, ];a <- group[1L];b <- group[2L]
      left <- min(H@i[H@p[group] + 1L]) + 1L
      right <- max(H@i[H@p[group + 1L]]) + 1L;rows <- seq.int(left, right)
      profiles <- as.matrix(H[rows, group, drop = FALSE])
      shape <- spolo_lc(as.numeric(profiles %*% rowSums(model$S[group, , drop = FALSE])))
      spectrum <- colSums(model$S[group, , drop = FALSE])
      before <- model$mu[rows, , drop = FALSE]
      after <- pmax(before - profiles %*% model$S[group, , drop = FALSE], 0) + shape$h %o% spectrum
      delta <- slcn_delta(Y[rows, , drop = FALSE], before, after, model$spec$kind,
        model$spec$dispersion, model$spec$theta, as.numeric(model$weights),
        nrow(Y), left - 1L, model$control$censor) - model$lambda
      if (!is.finite(delta) || delta >= best_delta) next
      best <- list(group = group, rows = rows, h = shape$h, S = spectrum, mu = after,
        gap = shape$gap, cosine = sum(H[, a] * H[, b]) / sqrt(norms[a] * norms[b]))
      best_delta <- delta
    }
    if (is.null(best)) break
    keep <- setdiff(seq_len(K), best$group)
    members <- sort(unique(unlist(model$members[best$group], use.names = FALSE)))
    h <- sparseMatrix(i = best$rows, j = rep(1L, length(best$rows)), x = best$h,
      dims = c(nrow(Y), 1L))
    model$H <- cbind(H[, keep, drop = FALSE], h)
    model$S <- rbind(model$S[keep, , drop = FALSE], best$S)
    model$mu[best$rows, ] <- best$mu
    model$members <- c(model$members[keep], list(members))
    audit[[length(audit) + 1L]] <- data.frame(source = paste(members, collapse = ";"),
      profile_cosine = best$cosine, LC_gap = best$gap, delta = best_delta)
  }
  model$profile_fusion <- list(merges = if (length(audit)) do.call(rbind, audit) else data.frame(),
    fitting_calls_in_fusion = 0L, minimum_cosine = .95, minimum_FWHM = 4 * psf_width)
  if (length(audit)) model$state_ids <- model$source_ids <- NULL
  model
}

slcn_fusion_spectra <- function(f, Y, H, S, fixed, rows, sweeps = 3L) {
  n <- nrow(Y);m <- ncol(Y);W <- f$weights
  if (length(W) > 1L) W <- W[rows, , drop = FALSE]
  if (ncol(H) <= 64L && nrow(H) <= 2000L && f$spec$kind %in% 0:2) {
    W <- matrix(rep_len(W, length(fixed)), nrow(fixed), m)
    W <- sweep(W, 2L, rep_len(f$spec$dispersion, m), "/")
    z <- if (f$spec$kind == 0L) nuisance_gaussian_cpp(as.matrix(H), Y[rows, , drop = FALSE],
      fixed, W, S, rep(1L, ncol(H)), as.integer(f$ncores)) else
      nuisance_poisson_cpp(as.matrix(H), Y[rows, , drop = FALSE], fixed, W, S,
        80L, 1e-8, as.integer(f$ncores))
    return(list(S = z$coefficients, mu = fixed + z$background))
  }
  z <- slcn_spectral(Y[rows, , drop = FALSE], as.matrix(H %*% S) + fixed,
    as(Matrix(H, sparse = TRUE), "dgCMatrix"), t(S), f$spec$kind,
    f$spec$dispersion, f$spec$theta, as.numeric(f$weights), n, min(rows) - 1L,
    as.integer(sweeps), as.integer(f$ncores), 20L, 1e-8, f$control$censor)
  list(S = t(z$S), mu = z$mu)
}

# A final local likelihood pass; neighboring modes and the background stay fixed.
slcn_fusion_refit <- function(f, Y, block, free) {
  if (!length(block)) return(f)
  rows <- which(rowSums(as.matrix(f$H[, block, drop = FALSE])) > 0)
  rows <- seq.int(max(1L, min(rows) - 2L), min(nrow(Y), max(rows) + 2L))
  h <- as.matrix(f$H[rows, block, drop = FALSE]);s <- f$S[block, , drop = FALSE]
  modes <- apply(h, 2L, which.max)
  radius <- pmax(1L, round(colSums(sweep(h, 2L, apply(h, 2L, max), "/") > .5) / 4))
  fixed <- pmax(f$mu[rows, , drop = FALSE] - h %*% s, 0)
  z <- slcn_fusion_spectra(f, Y, h, s, fixed, rows, 6L);s <- z$S;mean <- z$mu
  for (cycle in 1:4) {
    old <- mean
    for (k in order(rowSums(s) * apply(h, 2L, max), decreasing = TRUE)) {
      if (sum(s[k, ]) <= 0) next
      positive <- which(h[, k] > 0)
      at <- seq.int(max(1L, min(positive) - 2L), min(nrow(h), max(positive) + 2L))
      value <- h[at, k];inside <- which(value > 0)
      if (length(inside) < 2L) next
      eta <- numeric(length(value));eta[inside] <- log(value[inside])
      left <- seq_len(min(inside) - 1L);right <- seq.int(max(inside), length(value))
      eta[left] <- eta[min(inside)] + (left - min(inside)) * diff(eta[head(inside, 2L)])
      eta[right] <- eta[max(inside)] + (right - max(inside)) * diff(eta[tail(inside, 2L)])
      if (any(!is.finite(eta) | abs(eta) > 600)) {at <- positive;value <- h[at, k]} else value <- exp(eta)
      if (length(at) < 3L) next
      before <- mean[at, , drop = FALSE];rest <- pmax(before - h[at, k] %o% s[k, ], 0)
      guarded <- !block[k] %in% free
      lo <- if (guarded) max(0L, modes[k] - min(at) - radius[k]) else -1L
      hi <- if (guarded) min(length(at) - 1L, modes[k] - min(at) + radius[k]) else -1L
      first <- min(rows) + min(at) - 2L
      z <- slcn_log_profile(Y[rows[at], , drop = FALSE], rest + value %o% s[k, ], value,
        s[k, ], f$spec$kind, f$spec$dispersion, f$spec$theta, as.numeric(f$weights),
        nrow(Y), first, 2L, f$control$censor, 128L, as.integer(f$ncores), as.integer(lo), as.integer(hi))
      if (guarded && max(z$profile[seq.int(lo, hi) + 1L]) < max(z$profile) * (1 - 1e-8)) next
      area <- sum(z$profile);profile <- z$profile / area;spectrum <- s[k, ] * area
      after <- rest + profile %o% spectrum
      delta <- slcn_delta(Y[rows[at], , drop = FALSE], before, after, f$spec$kind,
        f$spec$dispersion, f$spec$theta, as.numeric(f$weights), nrow(Y), first, f$control$censor)
      if (is.finite(delta) && delta <= 1e-7) {h[at, k] <- profile;s[k, ] <- spectrum;mean[at, ] <- after}
    }
    z <- slcn_fusion_spectra(f, Y, h, s, fixed, rows);s <- z$S;mean <- z$mu
    delta <- slcn_delta(Y[rows, , drop = FALSE], old, mean, f$spec$kind,
      f$spec$dispersion, f$spec$theta, as.numeric(f$weights), nrow(Y), min(rows) - 1L, f$control$censor)
    if (cycle >= 3L && -delta < 1e-6) break
  }
  f$H[, block] <- 0;f$H[rows, block] <- Matrix(h, sparse = TRUE);f$H <- drop0(f$H)
  f$S[block, ] <- s;f$mu[rows, ] <- mean
  live <- rowSums(f$S) > 0
  f$H <- f$H[, live, drop = FALSE];f$S <- f$S[live, , drop = FALSE];f$members <- f$members[live]
  f
}

slcn_refine_profiles <- function(f, Y, verbose = FALSE, psf_width = f$control$psf_width) {
  audit <- list();begun <- proc.time()[["elapsed"]]
  accept <- function(label, expression) {
    seconds <- system.time(candidate <- force(expression))[["elapsed"]]
    delta <- if (identical(f$mu, candidate$mu) && ncol(candidate$H) == ncol(f$H)) 0 else
      slcn_delta(Y, f$mu, candidate$mu, f$spec$kind, f$spec$dispersion,
        f$spec$theta, as.numeric(f$weights), nrow(Y), 0L, f$control$censor) +
        f$lambda * (ncol(candidate$H) - ncol(f$H))
    ok <- is.finite(delta) && delta <= 1e-6
    audit[[length(audit) + 1L]] <<- data.frame(stage = label, seconds = seconds,
      delta = delta, accepted = ok, compounds = ncol(candidate$H))
    if (ok) f <<- candidate
    if (verbose) {cat(label, round(seconds, 2), "s; change", signif(delta, 7),
      if (ok) "accepted;" else "rejected;", ncol(f$H), "compounds\n");flush.console()}
  }
  accept("broad_profile_fusion", slcn_fuse_profiles(f, Y, psf_width))
  if (is.null(f$profile_fusion$merges) || !nrow(f$profile_fusion$merges)) {
    f$history <- rbind(f$history, do.call(rbind, audit))
    f$profile_fusion$seconds <- proc.time()[["elapsed"]] - begun
    return(f)
  }
  if (!is.null(f$profile_fusion$merges) && nrow(f$profile_fusion$merges)) {
    # Disjoint elution regions are refined separately, keeping updates local.
    free <- which(lengths(f$members) > 1L)
    ranges <- t(vapply(free, function(k) range(which(f$H[, k] > 0)), integer(2L)))
    ranges <- ranges[order(ranges[, 1L]), , drop = FALSE];regions <- list(ranges[1L, ])
    if (nrow(ranges) > 1L) for (i in 2:nrow(ranges)) {
      j <- length(regions)
      if (ranges[i, 1L] <= regions[[j]][2L]) regions[[j]][2L] <- max(regions[[j]][2L], ranges[i, 2L]) else
        regions[[j + 1L]] <- ranges[i, ]
    }
    accept("fusion_refit", {
      candidate <- f
      for (region in regions) {
        block <- which(colSums(candidate$H[seq.int(region[1L], region[2L]), , drop = FALSE]) > .02)
        candidate <- slcn_fusion_refit(candidate, Y, block, which(lengths(candidate$members) > 1L))
      }
      candidate
    })
  }
  f$history <- rbind(f$history, do.call(rbind, audit))
  f$rank <- ncol(f$H)
  f$loglik <- slcn_loglik(Y, f$mu, f$spec$kind, f$spec$dispersion,
    f$spec$theta, as.numeric(f$weights), as.integer(f$ncores), f$control$censor)
  f$objective <- -2 * f$loglik + f$lambda * f$rank
  f$complete_delta <- f$objective - f$initial_objective
  f$ledger_error <- f$complete_delta - sum(f$history$delta[f$history$accepted])
  f$profile_fusion$seconds <- proc.time()[["elapsed"]] - begun
  f
}
