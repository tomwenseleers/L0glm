# Hard log-concavity projection, shared by PSF calibration and NMF.
spolo_lc <- function(h, min_relative = 1e-12) {
  slcn_logconcave_cpp(as.numeric(h), min_relative)
}


slcn_context <- function(Y, B, spec, W, offset, ncores, censor) {
  spolo_lc <- spolo_lc
  ctx <- new.env(parent = emptyenv())
  ctx$Y <- Y;ctx$B <- B;ctx$W <- as.numeric(W);ctx$offset <- offset
  ctx$spec <- spec;ctx$n <- nrow(Y);ctx$m <- ncol(Y);ctx$censor <- censor;ctx$threads <- ncores
  ctx$background_fit <- function(y, fixed, alpha) {
    W <- matrix(rep_len(ctx$W, length(y)), nrow(y), ncol(y))
    if (ctx$spec$kind == 0L) return(nuisance_gaussian_cpp(B, y, fixed,
      W, alpha, rep(1L, ncol(B)), as.integer(ncores)))
    if (ctx$spec$kind %in% 1:2) return(nuisance_poisson_cpp(B, y, pmax(fixed, 0), W, alpha,
      40L, 1e-8, as.integer(ncores)))
    slcn_background(y, pmax(fixed, 0), B, alpha, ctx$spec$kind,
      ctx$spec$dispersion, ctx$spec$theta, ctx$W, as.integer(ncores), censor = ctx$censor)
  }


spolo_state <- function(model, d, lambda_rank) {
    st <- new.env(parent = emptyenv())
    st$n <- nrow(d$Y)
    st$m <- ncol(d$Y)
    st$B <- d$B
    st$alpha <- model$alpha
    st$T <- t(model$S)
    st$lambda <- lambda_rank
    st$h <- vector("list", ncol(model$H))
    st$lo <- st$hi <- integer(ncol(model$H))
    st$history <- list()
    st$fit_calls <- c(profile = 0L, background = 0L, residual = 0L)
    for (j in seq_len(ncol(model$H))) {
        rows <- which(model$H[, j] > 0)
        stopifnot(all(diff(rows) == 1L))
        st$lo[j] <- min(rows)
        st$hi[j] <- max(rows)
        st$h[[j]] <- as.numeric(model$H[rows, j])
    }
    st$handle <- spolo_state_cpp(d$Y, model$mu)
    st
}

spolo_record <- function(st, type, delta, accepted, id = NA_integer_) {
    st$history[[length(st$history) + 1L]] <- data.frame(type = type, id = id, delta = delta, accepted = accepted)
}

spolo_export <- function(st) {
    keep <- which(colSums(st$T) > 0)
    i <- as.integer(unlist(lapply(keep, function(k) seq.int(st$lo[k], st$hi[k])), use.names = FALSE))
    j <- rep(seq_along(keep), lengths(st$h[keep]))
    x <- as.numeric(unlist(st$h[keep], use.names = FALSE))
    H <- Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(st$n, length(keep)))
    list(H = H, S = t(st$T[, keep, drop = FALSE]), alpha = st$alpha, mu = spolo_mean_cpp(st$handle), 
        state_ids = keep, history = if (length(st$history)) do.call(rbind, st$history) else data.frame(), 
        fit_calls = st$fit_calls)
}

spolo_spectral <- function(st, ncores = 8L) {
    for (j in order(colSums(st$T), decreasing = TRUE)) {
        old <- st$T[, j]
        if (sum(old) <= 0) 
            next
        a <- spolo_local_cpp(st$handle, st$lo[j], st$hi[j])
        fit <- spolo_spectral_cpp(a$Y, a$mu, matrix(st$h[[j]], ncol = 1L), matrix(old, ncol = 1L), 1L, 
            ncores)
        delta <- fit$change - if (sum(fit$S) == 0) 
            st$lambda
        else 0
        ok <- is.finite(delta) && delta <= 1e-07
        if (ok) {
            spolo_commit_cpp(st$handle, st$lo[j], fit$mu)
            st$T[, j] <- fit$S[, 1L]
        }
        spolo_record(st, "spectrum", delta, ok, j)
    }
}



spolo_js <- function(p, q) {
    a <- (p + q)/2
    good <- p > 0
    goodq <- q > 0
    0.5 * sum(p[good] * log(p[good]/a[good])) + 0.5 * sum(q[goodq] * log(q[goodq]/a[goodq]))
}

spolo_fuse <- function(st, radius = 0.04, max_gap = 0.25) {
    accepted <- 0L
    visit <- order(colSums(st$T) * vapply(st$h, max, 0), decreasing = TRUE)
    for (original in visit) {
        anchor <- original
        repeat {
            area <- colSums(st$T)
            if (area[anchor] <= 0) 
                break
            p <- st$T[, anchor]/area[anchor]
            eligible <- setdiff(which(area > 0 & st$lo <= st$hi[anchor] & st$hi >= st$lo[anchor]), anchor)
            distance <- vapply(eligible, function(j) spolo_js(p, st$T[, j]/area[j]), 0)
            eligible <- eligible[order(distance)]
            eligible <- eligible[sort(distance) <= radius]
            if (!length(eligible)) 
                break
            merged <- FALSE
            for (other in eligible) {
                group <- c(anchor, other)
                left <- min(st$lo[group])
                right <- max(st$hi[group])
                n <- right - left + 1L
                H <- matrix(0, n, 2L)
                for (z in 1:2) {
                  k <- group[z]
                  H[seq.int(st$lo[k], st$hi[k]) - left + 1L, z] <- st$h[[k]]
                }
                shape <- spolo_lc(as.numeric(H %*% area[group]))
                if (shape$gap > max_gap) 
                  next
                s <- rowSums(st$T[, group, drop = FALSE])
                a <- spolo_local_cpp(st$handle, left, right)
                R <- pmax(a$mu - H %*% t(st$T[, group, drop = FALSE]), 0)
                mu <- R + shape$h %o% s
                delta <- spolo_delta_cpp(a$Y, a$mu, mu) - st$lambda
                ok <- is.finite(delta) && delta < -1e-06
                spolo_record(st, "fusion", delta, ok, anchor)
                if (ok) {
                  spolo_commit_cpp(st$handle, left, mu)
                  st$T[, group] <- 0
                  st$T <- cbind(st$T, s)
                  st$lo <- c(st$lo, left)
                  st$hi <- c(st$hi, right)
                  st$h[[length(st$h) + 1L]] <- shape$h
                  accepted <- accepted + 1L
                  anchor <- ncol(st$T)
                  merged <- TRUE
                  break
                }
            }
            if (!merged) 
                break
        }
    }
    accepted
}

spolo_prune <- function(st) {
    for (j in which(colSums(st$T) > 0)) {
        a <- spolo_local_cpp(st$handle, st$lo[j], st$hi[j])
        mu <- pmax(a$mu - st$h[[j]] %o% st$T[, j], 0)
        delta <- spolo_delta_cpp(a$Y, a$mu, mu) - st$lambda
        ok <- is.finite(delta) && delta < -1e-06
        if (ok) {
            spolo_commit_cpp(st$handle, st$lo[j], mu)
            st$T[, j] <- 0
        }
        spolo_record(st, "deletion", delta, ok, j)
    }
}

spolo_algebraic_fusion <- function(model, d, lambda_rank, radius = 0.08, max_gap = 0.5, ncores = 8L) {
    st <- spolo_state(model, d, lambda_rank)
    initial_area <- colSums(st$T)
    visit <- order(initial_area * vapply(st$h, max, 0), decreasing = TRUE)
    proposals <- 0L
    fusion_seconds <- system.time(for (anchor in visit) {
        area <- colSums(st$T)
        if (area[anchor] <= max(initial_area) * 1e-12) 
            next
        p <- st$T[, anchor]/area[anchor]
        group <- spolo_connected_cpp(st$T, area, p, st$lo, st$hi, anchor, radius)
        if (length(group) < 3L) 
            next
        p <- apply(sweep(st$T[, group, drop = FALSE], 2, area[group], "/"), 1, median)
        if (sum(p) <= 0) 
            next
        p <- p/sum(p)
        group <- spolo_connected_cpp(st$T, area, p, st$lo, st$hi, anchor, radius)
        if (length(group) < 3L) 
            next
        extracted <- spolo_extract_cpp(st$T[, group, drop = FALSE], area[group], p, 0.1, radius/5, ncores)
        if (sum(extracted$mass) <= 0) 
            next
        left <- min(st$lo[group])
        right <- max(st$hi[group])
        H <- matrix(0, right - left + 1L, length(group))
        for (z in seq_along(group)) {
            k <- group[z]
            H[seq.int(st$lo[k], st$hi[k]) - left + 1L, z] <- st$h[[k]]
        }
        raw <- as.numeric(H %*% extracted$mass)
        shape <- spolo_lc(raw)
        if (shape$gap > max_gap) 
            next
        s <- rowSums(extracted$E)
        local <- spolo_local_cpp(st$handle, left, right)
        mu <- pmax(local$mu - H %*% t(extracted$E), 0) + shape$h %o% s
        spolo_commit_cpp(st$handle, left, mu)
        st$T[, group] <- pmax(st$T[, group, drop = FALSE] - extracted$E, 0)
        st$T <- cbind(st$T, s)
        st$lo <- c(st$lo, left)
        st$hi <- c(st$hi, right)
        st$h[[length(st$h) + 1L]] <- shape$h
        proposals <- proposals + 1L
    })[["elapsed"]]
    list(model = spolo_export(st), seconds = fusion_seconds, proposals = proposals, radius = radius)
}

spolo_cluster <- function(H, S, psf_width, radius = 0.04, fractions = c(0.05, 0.1, 0.2), ncores = 8L) {
    timer <- proc.time()[["elapsed"]]
    T <- t(S)
    area <- colSums(T)
    K <- ncol(H)
    lo <- hi <- integer(K)
    height <- numeric(K)
    for (j in seq_len(K)) {
        at <- which(H[, j] > 0)
        lo[j] <- min(at)
        hi[j] <- max(at)
        height[j] <- max(H[at, j])
    }
    proposals <- list()
    owners <- integer()
    tested <- 0L
    for (anchor in order(area * height, decreasing = TRUE)) {
        p <- T[, anchor]/area[anchor]
        family <- spolo_connected_cpp(T, area, p, lo, hi, anchor, radius)
        broad <- 1/height[anchor] >= 3 * psf_width || (length(family) >= 3L && max(hi[family]) - min(lo[family]) + 
            1L >= 8 * psf_width)
        if (!broad) 
            next
        if (any(owners %in% family)) 
            next
        owners <- c(owners, anchor)
        tested <- tested + 1L
        extraction <- spolo_extract_cpp(T, area, p, 0.1, radius/5, ncores)
        groups <- list(family)
        kinds <- "whole_spectrum"
        for (fraction in fractions) {
            mass <- extraction$mass
            mass[mass/area < fraction] <- 0
            group <- spolo_connected_cpp(extraction$E, mass, p, lo, hi, anchor, radius)
            if (length(group) >= 2L && !any(vapply(groups, function(g) identical(g, group), FALSE))) {
                groups[[length(groups) + 1L]] <- group
                kinds <- c(kinds, paste0("contained_", fraction))
            }
        }
        alternatives <- list()
        for (j in seq_along(groups)) {
            group <- groups[[j]]
            if (length(group) < 2L) 
                next
            E <- extraction$E[, group, drop = FALSE]
            raw <- as.numeric(H[, group, drop = FALSE] %*% colSums(E))
            at <- which(raw > 0)
            if (length(at) < 2L) 
                next
            shape <- spolo_lc(raw[min(at):max(at)])
            alternatives[[length(alternatives) + 1L]] <- list(anchor = anchor, group = group, E = E, 
                left = min(at), right = max(at), h = shape$h, S = rowSums(E), gap = shape$gap, kind = kinds[j])
        }
        if (length(alternatives)) 
            proposals[[length(proposals) + 1L]] <- alternatives
    }
    list(proposals = proposals, seconds = proc.time()[["elapsed"]] - timer, anchors_tested = tested, 
        fitting_calls = 0L)
}

spolo_trial <- function(st, proposal, sweeps = 4L, ncores = 8L) {
    group <- proposal$group
    if (any(colSums(st$T[, group, drop = FALSE]) <= 0)) 
        return(NULL)
    E <- pmin(proposal$E, st$T[, group, drop = FALSE])
    left <- min(st$lo[group], proposal$left)
    right <- max(st$hi[group], proposal$right)
    n <- right - left + 1L
    g <- length(group)
    H <- matrix(0, n, g + 1L)
    for (j in seq_along(group)) {
        k <- group[j]
        H[seq.int(st$lo[k], st$hi[k]) - left + 1L, j] <- st$h[[k]]
    }
    H[seq.int(proposal$left, proposal$right) - left + 1L, g + 1L] <- proposal$h
    a <- spolo_local_cpp(st$handle, left, right)
    offset <- pmax(a$mu - H[, seq_len(g), drop = FALSE] %*% t(st$T[, group, drop = FALSE]), 0)
    S <- cbind(pmax(st$T[, group, drop = FALSE] - E, 0), rowSums(E))
    if (sum(S[, g + 1L]) <= 0) 
        return(NULL)
    mu <- offset + H %*% t(S)
    active <- which(colSums(S) > 0)
    fitted <- spolo_spectral_cpp(a$Y, mu, H[, active, drop = FALSE], S[, active, drop = FALSE], sweeps, 
        ncores)
    S[] <- 0
    S[, active] <- fitted$S
    mu <- fitted$mu
    for (j in order(colSums(S))) {
        if (sum(S[, j]) <= 0) 
            next
        next_mu <- pmax(mu - H[, j] %o% S[, j], 0)
        change <- spolo_delta_cpp(a$Y, mu, next_mu) - st$lambda
        if (is.finite(change) && change < 0) {
            S[, j] <- 0
            mu <- next_mu
        }
    }
    if (sum(S[, g + 1L]) <= 0) 
        return(NULL)
    nextK <- sum(colSums(S) > 0)
    delta <- spolo_delta_cpp(a$Y, a$mu, mu) + st$lambda * (nextK - g)
    list(group = group, left = left, right = right, H = H, S = S, mu = mu, delta = delta, gap = proposal$gap, 
        kind = proposal$kind, anchor = proposal$anchor)
}

spolo_complete <- function(model, d, ncores = 8L, sweeps = 4L, fractions = c(0.05, 0.1,
    0.2)) {
    total <- proc.time()[["elapsed"]]
    psf_width <- d$psf_width
    clustered <- spolo_cluster(model$H, model$S, psf_width, fractions = fractions, ncores = ncores)
    st <- spolo_state(model, d, d$lambda_rank)
    audit <- list()
    accepted <- 0L
    refit0 <- proc.time()[["elapsed"]]
    for (family in clustered$proposals) {
        best <- NULL
        for (proposal in family) {
            trial <- spolo_trial(st, proposal, sweeps, ncores)
            if (is.null(trial)) 
                next
            audit[[length(audit) + 1L]] <- data.frame(anchor = trial$anchor, kind = trial$kind, donors = length(trial$group), 
                gap = trial$gap, delta = trial$delta, accepted = FALSE)
            trial$audit_id <- length(audit)
            if (is.finite(trial$delta) && trial$delta < -1e-06 && (is.null(best) || trial$delta < best$delta)) 
                best <- trial
        }
        if (is.null(best)) 
            next
        g <- length(best$group)
        spolo_commit_cpp(st$handle, best$left, best$mu)
        st$T[, best$group] <- best$S[, seq_len(g), drop = FALSE]
        st$T <- cbind(st$T, best$S[, g + 1L])
        at <- which(best$H[, g + 1L] > 0)
        st$lo <- c(st$lo, best$left + min(at) - 1L)
        st$hi <- c(st$hi, best$left + max(at) - 1L)
        st$h[[length(st$h) + 1L]] <- best$H[min(at):max(at), g + 1L]
        spolo_record(st, "completion", best$delta, TRUE, best$anchor)
        audit[[best$audit_id]]$accepted <- TRUE
        accepted <- accepted + 1L
    }
    out <- spolo_export(st)
    out$seconds <- c(clustering = clustered$seconds, local_refit = proc.time()[["elapsed"]] - refit0, 
        total = proc.time()[["elapsed"]] - total)
    out$completion <- list(accepted = accepted, anchors = clustered$anchors_tested, fitting_calls_in_clustering = 0L, 
        psf_width = psf_width, sweeps = sweeps, fractions = fractions, audit = if (length(audit)) do.call(rbind,
            audit) else data.frame())
    out
}

spolo_coherent <- function(model, d, ncores = 8L, sweeps = 8L, background = FALSE) {
    t0 <- proc.time()[["elapsed"]]
    algebraic <- spolo_algebraic_fusion(model, d, d$lambda_rank, radius = 0.04, max_gap = Inf, ncores = ncores)
    out <- algebraic$model
    fit0 <- proc.time()[["elapsed"]]
    polished <- spolo_spectral_cpp(d$Y, out$mu, out$H, t(out$S), sweeps, ncores)
    out$S <- t(polished$S)
    out$mu <- polished$mu
    st <- spolo_state(out, d, d$lambda_rank)
    spolo_prune(st)
    if (background) 
        spolo_background(st, d, ncores)
    out <- spolo_export(st)
    out$seconds <- c(clustering = algebraic$seconds, local_refit = proc.time()[["elapsed"]] - fit0, total = proc.time()[["elapsed"]] - 
        t0)
    out$coherent <- list(proposals = algebraic$proposals, sweeps = sweeps, fitting_calls_in_clustering = 0L, 
        background = background)
    out
}

spolo_fast <- function(model, d, ncores = 8L, sweeps = 8L, block_sweeps = 12L, background = FALSE) {
    a <- spolo_coherent(model, d, ncores, sweeps, background)
    b <- spolo_complete(a, d, ncores, block_sweeps)
    b$pipeline <- list(coherent = a$coherent, completion = b$completion, seconds = rbind(coherent = a$seconds, 
        completion = b$seconds))
    b$seconds <- a$seconds + b$seconds
    b
}

spolo_broad_polish <- function(model, d, ncores = 8L, prior_strength = 0) {
    t0 <- proc.time()[["elapsed"]]
    width <- d$psf_width
    st <- spolo_state(model, d, d$lambda_rank)
    height <- vapply(st$h, max, 0)
    broad <- which(1/height >= 3 * width)
    broad <- broad[order(colSums(st$T)[broad] * height[broad], decreasing = TRUE)]
    for (j in broad) spolo_profile(st, j, prior_strength, ncores, pad = width)
    out <- spolo_export(st)
    out$seconds <- proc.time()[["elapsed"]] - t0
    out$polish <- list(candidates = length(broad), accepted = sum(out$history$accepted), prior_strength = prior_strength)
    out
}

spolo_sweep <- function(model, d, cycles, ncores = 8L) {
    t0 <- proc.time()[["elapsed"]]
    st <- spolo_state(model, d, d$lambda_rank)
    fusions <- spolo_fuse(st)
    for (cycle in seq_len(cycles)) {
        spolo_spectral(st, ncores)
        visit <- order(colSums(st$T), decreasing = TRUE)
        for (z in seq_along(visit)) {
            spolo_profile(st, visit[z], prior_strength = d$lambda, ncores = ncores)
        }
        spolo_background(st, d, ncores)
        spolo_prune(st)
        fusions <- fusions + spolo_fuse(st)
    }
    out <- spolo_export(st)
    out$seconds <- proc.time()[["elapsed"]] - t0
    out$control <- list(cycles = cycles, prior_strength = d$lambda, accepted_fusions = fusions, ncores = ncores)
    out
}


  first <- function(y) if (is.null(attr(y, "row_start"))) 0L else as.integer(attr(y, "row_start")) - 1L
  spolo_state_cpp <- slcn_state
  spolo_local_cpp <- slcn_local
  spolo_commit_cpp <- slcn_commit
  spolo_mean_cpp <- slcn_mean
  spolo_delta_cpp <- function(y, old_mu, new_mu) slcn_delta(y, old_mu, new_mu,
    ctx$spec$kind, ctx$spec$dispersion, ctx$spec$theta, ctx$W, ctx$n, first(y), ctx$censor)
  spolo_spectral_cpp <- function(y, mean, profiles, spectra, sweeps = 4L, nthreads = ncores,
                                    maxit = 12L, tol = 1e-8) {
    if (!inherits(profiles, "dgCMatrix")) profiles <- as(Matrix(profiles, sparse = TRUE), "dgCMatrix")
    if (ctx$spec$kind == 0L && ncol(profiles) <= 32L && nrow(profiles) <= 2000L) {
      signal <- pmax(mean - as.matrix(profiles %*% t(spectra)), 0)
      weights <- if (length(ctx$W) == 1L) matrix(ctx$W, nrow(y), ncol(y)) else
        matrix(ctx$W, ctx$n, ctx$m)[first(y) + seq_len(nrow(y)), , drop = FALSE]
      a <- nuisance_gaussian_cpp(as.matrix(profiles), y, signal,
        weights, t(spectra), rep(1L, ncol(profiles)), as.integer(nthreads))
      fitted <- signal + a$background
      delta <- spolo_delta_cpp(y, mean, fitted)
      if (is.finite(delta) && delta <= 1e-7)
        return(list(mu = fitted, S = t(a$coefficients), change = delta, kkt = a$kkt))
    }
    slcn_spectral(y, mean, profiles, spectra, ctx$spec$kind,
      ctx$spec$dispersion, ctx$spec$theta, ctx$W, ctx$n, first(y),
      as.integer(sweeps), as.integer(nthreads), as.integer(maxit), tol, ctx$censor)
  }
  spolo_background <- function(st, d, ncores = 8L) {
    if (!ncol(st$B)) return(invisible(NULL))
    old <- spolo_mean_cpp(st$handle)
    fixed <- pmax(old - st$B %*% st$alpha, 0)
    fit <- ctx$background_fit(d$Y, fixed, st$alpha)
    fit <- list(mu = fixed + fit$background, S = t(fit$coefficients))
    delta <- spolo_delta_cpp(d$Y, old, fit$mu)
    ok <- is.finite(delta) && delta <= 1e-7
    if (ok) {spolo_commit_cpp(st$handle, 1L, fit$mu);st$alpha <- t(fit$S)}
    st$fit_calls["background"] <- st$fit_calls["background"] + 1L
    spolo_record(st, "background", delta, ok)
  }
  spolo_profile <- function(st, j, prior_strength = 0, ncores = 8L, pad = 3L) {
    s <- st$T[, j];if (sum(s) <= 0) return(invisible(FALSE))
    left <- max(1L, st$lo[j] - pad);right <- min(st$n, st$hi[j] + pad)
    h <- numeric(right - left + 1L);h[seq.int(st$lo[j], st$hi[j]) - left + 1L] <- st$h[[j]]
    a <- spolo_local_cpp(st$handle, left, right)
    residual <- pmax(a$mu - h %o% s, 0)
    if (ctx$spec$kind == 1L) {
      W <- if (length(ctx$W)==1L) matrix(ctx$W,length(h),length(s)) else
        matrix(ctx$W,ctx$n,ctx$m)[left:right,,drop=FALSE]
      fit <- L0glm.fit(matrix(s,ncol=1L),t(a$Y),family=poisson("identity"),
        weights=t(W),offset=t(residual),start=matrix(pmax(h,max(h)*1e-8),nrow=1L),
        prior.mean=if(prior_strength>0) matrix(h,nrow=1L) else NULL,
        lambda=prior_strength,nonnegative=TRUE,intercept=FALSE,normalize="none",
        grouped=FALSE,relaxed=FALSE,maxit=4L,thresh=0,tol=1e-7,ncores=ncores,
        control.solver=list(.nuisance=FALSE))
      raw <- as.numeric(coef(fit));st$fit_calls["profile"] <- st$fit_calls["profile"]+1L
    } else raw <- slcn_profile(a$Y, a$mu, h, s, ctx$spec$kind,
        ctx$spec$dispersion, ctx$spec$theta, ctx$W, ctx$n, left - 1L,
        8L, as.integer(ncores), ctx$censor, prior_strength)
    total <- sum(raw)
    if (!is.finite(total) || total <= 0) return(invisible(FALSE))
    shape <- spolo_lc(raw)$h;next_s <- s * total
    candidate_mu <- residual + shape %o% next_s
    polished <- spolo_spectral_cpp(a$Y, candidate_mu, matrix(shape, ncol = 1L), matrix(next_s, ncol = 1L),
      if(ctx$spec$kind==1L) 1L else 2L, ncores)
    delta <- spolo_delta_cpp(a$Y, a$mu, polished$mu) - if(sum(polished$S)==0) st$lambda else 0
    ok <- is.finite(delta) && delta < -1e-6
    if (ok) {spolo_commit_cpp(st$handle, left, polished$mu);st$lo[j] <- left;st$hi[j] <- right;
      st$h[[j]] <- shape;st$T[, j] <- polished$S[, 1L]}
    spolo_record(st, "profile", delta, ok, j)
    invisible(ok)
  }
  list(engine = environment(), context = ctx)
}
