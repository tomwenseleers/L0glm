# Original-channel likelihood refinement with hard log-concave profiles.
slcn_refine_fit <-
function (f, Y, verbose = FALSE)
{
    begun <- proc.time()[["elapsed"]]
    audit <- list()
    criterion <- function(z) -2 * slcn_loglik(Y, z$mu, z$spec$kind,
        z$spec$dispersion, z$spec$theta, as.numeric(z$weights),
        as.integer(z$ncores), z$control$censor) + z$lambda *
        ncol(z$H)
    accept <- function(candidate, label, seconds) {
        delta <- criterion(candidate) - criterion(f)
        ok <- is.finite(delta) && delta <= 1e-06
        audit[[length(audit) + 1L]] <<- data.frame(stage = label,
            seconds, delta, accepted = ok, compounds = ncol(candidate$H))
        if (verbose) {
            cat(label, round(seconds, 2), "s; change", signif(delta,
                7), if (ok)
                "accepted;"
            else "rejected;", ncol(candidate$H), "compounds\n")
            flush.console()
        }
        if (ok)
            f <<- candidate
        ok
    }
    f$members <- lapply(seq_len(ncol(f$H)), identity)
    seconds <- system.time(z <- slcn_refine_block(f, Y, cycles = f$control$profile_cycles,
        profile_steps = 1L, adaptive = FALSE, joint = FALSE,
        background = TRUE))[["elapsed"]]
    accept(z, "smooth_profiles", seconds)
    seconds <- system.time(z <- slcn_joint_fixed(f, Y))[["elapsed"]]
    accept(z, "joint_spectra", seconds)
    for (pairs in c(FALSE, TRUE)) for (pass in seq_len(f$control$fusion_passes)) {
        proposals <- slcn_anchor_proposals(f, pairs)
        if (!length(proposals))
            break
        changed <- FALSE
        for (proposal in head(proposals, 2L)) {
            seconds <- system.time({
                z <- slcn_merge_anchor(f, Y, proposal)
                if (!pairs)
                  z <- slcn_restore_atoms(z, Y, f, proposal$donors,
                    criterion)
                z <- slcn_joint_fixed(z, Y)
                if (pairs && identical(proposal$shape, "anchor") &&
                    criterion(z) >= criterion(f) - 1e-6) {
                    proposal$shape <- "hull"
                    alternative <- slcn_merge_anchor(f, Y, proposal, guarded = TRUE)
                    alternative <- slcn_joint_fixed(alternative, Y)
                    if (criterion(alternative) < criterion(z)) z <- alternative
                }
            })[["elapsed"]]
            if (accept(z, if (pairs)
                "contained_tail_fusion"
            else "broad_anchor_fusion", seconds)) {
                changed <- TRUE
                break
            }
        }
        if (!changed)
            break
    }
    f$rank <- ncol(f$H)
    f$loglik <- (f$lambda * f$rank - criterion(f))/2
    f$objective <- criterion(f)
    f$ic_df <- ncol(Y) * f$rank + f$nuisance_df
    f$history <- rbind(f$history, do.call(rbind, audit))
    f$complete_delta <- f$objective - f$initial_objective
    f$ledger_error <- f$complete_delta - sum(f$history$delta[f$history$accepted])
    f$seconds["refinement"] <- proc.time()[["elapsed"]] - begun
    slcn_temporal_order(f)
}
slcn_profile_fit <- function(y, mean, h, s, family, phi, theta, weights,
                             full_n, first, maxit, censor, q, nthreads, adaptive) {
  slcn_log_profile(y, mean, h, s, family, phi, theta, weights, full_n,
    first, maxit, censor, q, nthreads, mode = if (adaptive) 1L else 0L)
}
slcn_positive_tail <-
function (h)
{
    at <- which(h > 0)
    if (!length(at))
        return(NULL)
    if (length(at) == 1L)
        return(pmax(h, max(h) * 1e-12))
    eta <- numeric(length(h))
    eta[at] <- log(h[at])
    left <- seq_len(min(at) - 1L)
    right <- seq.int(max(at), length(h))
    eta[left] <- eta[min(at)] + (left - min(at)) * diff(eta[head(at,
        2L)])
    eta[right] <- eta[max(at)] + (right - max(at)) * diff(eta[tail(at,
        2L)])
    if (any(!is.finite(eta) | abs(eta) > 600))
        NULL
    else exp(eta)
}
slcn_joint_spectra <-
function (f, Y, H, S, fixed, rows = seq_len(nrow(Y)), maxit = 80L,
    joint = TRUE, sweeps = 3L)
{
    n <- nrow(Y)
    m <- ncol(Y)
    W <- f$weights
    if (length(W) > 1L)
        W <- W[rows, , drop = FALSE]
    phi <- rep_len(f$spec$dispersion, m)
    mean <- as.matrix(H %*% S) + fixed
    if (joint && ncol(H) <= 64L && nrow(H) <= 2000L && f$spec$kind %in%
        0:2) {
        W <- matrix(rep_len(W, length(mean)), nrow(mean), m)
        W <- sweep(W, 2L, phi, "/")
        z <- if (f$spec$kind == 0L)
            nuisance_gaussian_cpp(as.matrix(H), Y[rows, , drop = FALSE],
                fixed, W, S, rep(1L, ncol(H)), as.integer(f$ncores))
        else nuisance_poisson_cpp(as.matrix(H), Y[rows, , drop = FALSE],
            fixed, W, S, as.integer(maxit), 1e-08, as.integer(f$ncores))
        return(list(S = z$coefficients, mu = fixed + z$background))
    }
    z <- slcn_spectral(Y[rows, , drop = FALSE], mean, as(Matrix(H,
        sparse = TRUE), "dgCMatrix"), t(S), f$spec$kind, f$spec$dispersion,
        f$spec$theta, as.numeric(f$weights), n, min(rows) - 1L,
        as.integer(sweeps), as.integer(f$ncores), 20L, 1e-08,
        f$control$censor)
    list(S = t(z$S), mu = z$mu)
}
slcn_refine_block <-
function (f, Y, block = seq_len(ncol(f$H)), cycles = 2L, profile_steps = 2L,
    adaptive = TRUE, joint = TRUE, pad = 2L, background = FALSE, free = NULL)
{
    if (!length(block))
        return(f)
    rows <- which(rowSums(as.matrix(f$H[, block, drop = FALSE])) >
        0)
    rows <- seq.int(max(1L, min(rows) - pad), min(nrow(Y), max(rows) +
        pad))
    if (background)
        rows <- seq_len(nrow(Y))
    h <- as.matrix(f$H[rows, block, drop = FALSE])
    s <- f$S[block, , drop = FALSE]
    if (!is.null(free)) {
        modes <- apply(h, 2L, which.max)
        radius <- pmax(1L, round(colSums(sweep(h, 2L, apply(h, 2L, max), "/") > .5) / 4))
    }
    fixed <- pmax(f$mu[rows, , drop = FALSE] - h %*% s, 0)
    A <- f$alpha
    B <- f$background[rows, , drop = FALSE]
    K <- ncol(h)
    q <- ncol(B)
    if (background)
        fixed <- pmax(fixed - B %*% A, 0)
    refit <- function(s, A, sweeps = 3L) {
        z <- slcn_joint_spectra(f, Y, if (background)
            cbind(h, B)
        else h, if (background)
            rbind(s, A)
        else s, fixed, rows, joint = joint, sweeps = sweeps)
        list(S = z$S[seq_len(K), , drop = FALSE], A = if (background) z$S[K +
            seq_len(q), , drop = FALSE] else A, mu = z$mu)
    }
    z <- refit(s, A, 6L)
    s <- z$S
    A <- z$A
    mean <- z$mu
    for (cycle in seq_len(cycles)) {
        old <- mean
        for (k in order(rowSums(s) * apply(h, 2, max), decreasing = TRUE)) {
            if (sum(s[k, ]) <= 0)
                next
            at <- which(h[, k] > 0)
            at <- seq.int(max(1L, min(at) - pad), min(nrow(h),
                max(at) + pad))
            positive <- slcn_positive_tail(h[at, k])
            if (is.null(positive)) {
                at <- which(h[, k] > 0)
                positive <- h[at, k]
            }
            if (length(at) < 3L)
                next
            before <- mean[at, , drop = FALSE]
            rest <- pmax(before - h[at, k] %o% s[k, ], 0)
            first <- min(rows) + min(at) - 2L
            guarded <- !is.null(free) && !block[k] %in% free
            lo <- if (guarded) max(0L, modes[k] - min(at) - radius[k]) else -1L
            hi <- if (guarded) min(length(at) - 1L, modes[k] - min(at) + radius[k]) else -1L
            if (guarded && lo > hi) next
            fit <- if (is.null(free)) slcn_profile_fit(
                Y[rows[at], , drop = FALSE], rest + positive %o% s[k, ], positive,
                s[k, ], f$spec$kind, f$spec$dispersion, f$spec$theta,
                as.numeric(f$weights), nrow(Y), first, as.integer(profile_steps),
                f$control$censor, if (adaptive) 48L else 12L,
                as.integer(f$ncores), adaptive) else slcn_log_profile(
                Y[rows[at], , drop = FALSE], rest + positive %o% s[k, ], positive,
                s[k, ], f$spec$kind, f$spec$dispersion, f$spec$theta,
                as.numeric(f$weights), nrow(Y), first, as.integer(profile_steps),
                f$control$censor, 128L, as.integer(f$ncores),
                as.integer(lo), as.integer(hi), mode = 2L)
            if (guarded && max(fit$profile[seq.int(lo, hi) + 1L]) < max(fit$profile) * (1 - 1e-8)) next
            a <- sum(fit$profile)
            hh <- fit$profile/a
            ss <- s[k, ] * a
            proposed <- rest + hh %o% ss
            delta <- slcn_delta(Y[rows[at], , drop = FALSE],
                before, proposed, f$spec$kind, f$spec$dispersion,
                f$spec$theta, as.numeric(f$weights), nrow(Y),
                first, f$control$censor)
            if (is.finite(delta) && delta <= 1e-07) {
                h[at, k] <- hh
                s[k, ] <- ss
                mean[at, ] <- proposed
            }
        }
        z <- refit(s, A)
        s <- z$S
        A <- z$A
        mean <- z$mu
        delta <- slcn_delta(Y[rows, , drop = FALSE], old, mean,
            f$spec$kind, f$spec$dispersion, f$spec$theta, as.numeric(f$weights),
            nrow(Y), min(rows) - 1L, f$control$censor)
        W <- if (length(f$weights) == 1L)
            f$weights
        else f$weights[rows, , drop = FALSE]
        ll <- slcn_loglik(Y[rows, , drop = FALSE], old, f$spec$kind,
            f$spec$dispersion, f$spec$theta, as.numeric(W), as.integer(f$ncores),
            f$control$censor)
        if (cycle >= 3L && -delta < 1e-07 *
            max(1, abs(ll)))
            break
    }
    f$H[, block] <- 0
    f$H[rows, block] <- Matrix(h, sparse = TRUE)
    f$S[block, ] <- s
    f$alpha <- A
    f$mu[rows, ] <- mean
    live <- rowSums(f$S) > 0
    f$H <- f$H[, live, drop = FALSE]
    f$S <- f$S[live, , drop = FALSE]
    if (!is.null(f$members))
        f$members <- f$members[live]
    f
}
slcn_joint_fixed <-
function (f, Y, maxit = 100L)
{
    K <- ncol(f$H)
    q <- ncol(f$background)
    H <- cbind(f$H, Matrix(f$background, sparse = TRUE))
    S <- rbind(f$S, f$alpha)
    if (!ncol(H))
        return(f)
    fixed <- matrix(rep_len(f$offset, length(Y)), nrow(Y), ncol(Y))
    z <- slcn_joint_spectra(f, Y, H, S, fixed, maxit = maxit,
        joint = TRUE, sweeps = 8L)
    f$S <- z$S[seq_len(K), , drop = FALSE]
    f$alpha <- z$S[K + seq_len(q), , drop = FALSE]
    f$mu <- z$mu
    live <- rowSums(f$S) > 0
    f$H <- f$H[, live, drop = FALSE]
    f$S <- f$S[live, , drop = FALSE]
    if (!is.null(f$members))
        f$members <- f$members[live]
    f
}
slcn_profile_geometry <-
function (f)
{
    H <- f$H
    K <- ncol(H)
    width <- left <- right <- apex <- integer(K)
    for (k in seq_len(K)) {
        at <- seq.int(H@p[k] + 1L, H@p[k + 1L])
        value <- H@x[at]
        ix <- H@i[at] + 1L
        significant <- value > max(value) * 1e-06
        width[k] <- sum(value > max(value) * 0.5)
        left[k] <- min(ix[significant])
        right[k] <- max(ix[significant])
        apex[k] <- ix[which.max(value)]
    }
    list(width = width, left = left, right = right, apex = apex)
}
slcn_connected <-
function (a)
{
    result <- list()
    remaining <- seq_len(nrow(a))
    while (length(remaining)) {
        group <- front <- remaining[1L]
        while (length(front)) {
            next_group <- setdiff(which(colSums(a[front, , drop = FALSE]) >
                0), group)
            group <- c(group, next_group)
            front <- next_group
        }
        result[[length(result) + 1L]] <- sort(group)
        remaining <- setdiff(remaining, group)
    }
    result
}
slcn_anchor_proposals <-
function (f, pairs = FALSE)
{
    K <- ncol(f$H)
    if (K < 2L)
        return(list())
    g <- slcn_profile_geometry(f)
    S <- f$S
    cosine <- tcrossprod(S)/sqrt(rowSums(S^2) %o% rowSums(S^2))
    result <- list()
    broad <- which(g$width > max(4 * median(g$width), 30))
    for (anchor in broad[order(g$width[broad], decreasing = TRUE)]) {
        support <- seq.int(g$left[anchor], g$right[anchor])
        mass <- as.numeric(colSums(f$H[support, , drop = FALSE]))
        pool <- which(seq_len(K) != anchor & mass >= 0.9)
        if (pairs) {
            pool <- which(seq_len(K) != anchor & cosine[anchor,
                ] >= 0.995 & g$left <= g$right[anchor] + 4 *
                median(g$width) & g$right >= g$left[anchor] -
                4 * median(g$width))
            for (other in pool) {
                shape <- if (mass[other] >= 0.9 &&
                  g$width[other] * 3 < g$width[anchor])
                  "anchor"
                else "hull"
                if (shape == "hull") {
                  raw <- as.numeric(f$H[, c(anchor, other), drop = FALSE] %*%
                    rowSums(S[c(anchor, other), , drop = FALSE]))
                  at <- which(raw > 0)
                  projected <- spolo_lc(raw[min(at):max(at)])
                  if (projected$gap > 0.25)
                    next
                }
                result[[length(result) + 1L]] <- list(anchor = anchor,
                  group = c(anchor, other), donors = other, strength = cosine[anchor,
                    other], shape = shape)
            }
        }
        else if (length(pool) >= 3L) {
            for (component in slcn_connected(cosine[pool, pool,
                drop = FALSE] >= 0.985)) {
                donors <- pool[component]
                if (length(donors) < 3L)
                  next
                centre <- colSums(S[donors, , drop = FALSE])
                compatibility <- sum(centre * S[anchor, ])/sqrt(sum(centre^2) *
                  sum(S[anchor, ]^2))
                if (compatibility < 0.4 || which.max(centre) !=
                  which.max(S[anchor, ]))
                  next
                result[[length(result) + 1L]] <- list(anchor = anchor,
                  group = sort(c(anchor, donors)), donors = donors,
                  strength = sum(S[donors, ]))
            }
        }
    }
    if (length(result))
        result <- result[order(vapply(result, function(z) z$strength,
            0), decreasing = TRUE)]
    result
}
slcn_merge_anchor <-
function (f, Y, proposal, cycles = 8L, guarded = FALSE)
{
    group <- proposal$group
    anchor <- proposal$anchor
    keep <- setdiff(seq_len(ncol(f$H)), group)
    h <- f$H[, anchor, drop = FALSE]
    if (identical(proposal$shape, "hull")) {
        raw <- as.numeric(f$H[, group, drop = FALSE] %*% rowSums(f$S[group,
            , drop = FALSE]))
        at <- which(raw > 0)
        at <- min(at):max(at)
        h <- sparseMatrix(i = at, j = rep(1L, length(at)), x = spolo_lc(raw[at])$h,
            dims = c(nrow(Y), 1L))
    }
    f$H <- cbind(f$H[, keep, drop = FALSE], h)
    f$S <- rbind(f$S[keep, , drop = FALSE], colSums(f$S[group,
        , drop = FALSE]))
    f$members <- c(f$members[keep], list(sort(unique(unlist(f$members[group])))))
    f$mu <- as.matrix(f$H %*% f$S) + f$background %*% f$alpha + f$offset
    new <- ncol(f$H)
    support <- which(f$H[, new] > max(f$H[, new]) * 1e-06)
    block <- which(colSums(f$H[support, , drop = FALSE]) > 0.025)
    if (guarded && nrow(Y) <= 2000L && ncol(f$H) <= 64L)
        block <- seq_len(ncol(f$H))
    slcn_refine_block(f, Y, block, cycles = cycles, profile_steps = 4L,
        pad = 2L, background = guarded, free = if (guarded) new else NULL)
}
slcn_restore_atoms <-
function (f, Y, original, pool, criterion)
{
    for (pass in seq_len(min(3L, length(pool)))) {
        best <- NULL
        best_score <- criterion(f)
        best_id <- NA_integer_
        for (k in pool) {
            candidate <- f
            h <- original$H[, k, drop = FALSE]
            z <- slcn_joint_spectra(f, Y, h, matrix(0, 1L, ncol(Y)),
                f$mu, maxit = 60L)
            candidate$H <- cbind(f$H, h)
            candidate$S <- rbind(f$S, z$S)
            candidate$mu <- z$mu
            candidate$members <- c(f$members, original$members[k])
            score <- criterion(candidate)
            if (is.finite(score) && score < best_score - 1e-06) {
                best <- candidate
                best_score <- score
                best_id <- k
            }
        }
        if (is.null(best))
            break
        f <- best
        pool <- setdiff(pool, best_id)
        support <- which(f$H[, ncol(f$H)] > max(f$H[, ncol(f$H)]) *
            1e-05)
        block <- which(colSums(f$H[support, , drop = FALSE]) >
            0.025)
        z <- slcn_refine_block(f, Y, block, cycles = 2L, profile_steps = 2L)
        if (criterion(z) < criterion(f))
            f <- z
    }
    f
}
