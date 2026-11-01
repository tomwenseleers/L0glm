slcn_design <- function(X, n, name) {
  X <- drop0(as(as(Matrix(X, sparse = TRUE), "generalMatrix"), "dgCMatrix"))
  if (nrow(X) != n || any(!is.finite(X@x) | X@x <= 0) || any(colSums(X) <= 0))
    stop(name, " must have finite nonnegative columns with positive mass.")
  for (j in seq_len(ncol(X))) {
    z <- seq.int(X@p[j] + 1L, X@p[j + 1L]);h <- X@x[z];i <- X@i[z]
    if (any(diff(i) != 1L) || any(diff(log(h), differences = 2L) > 1e-6))
      stop(name, " columns need contiguous log-concave positive support.")
  }
  X
}

slcn_start <- function(start, n, m, B, O) {
  if (!is.list(start) || is.null(start$H) || is.null(start$S)) stop("start must contain H and S.")
  H <- slcn_design(start$H, n, "start$H");K <- ncol(H);S <- as.matrix(start$S)
  A <- if (is.null(start$alpha)) matrix(0, ncol(B), m) else as.matrix(start$alpha)
  if (!identical(dim(S), c(K, m)) || !identical(dim(A), c(ncol(B), m)) ||
      any(!is.finite(S) | S < 0) || any(!is.finite(A) | A < 0)) stop("Invalid start spectra/background dimensions or values.")
  sums <- colSums(H);H <- H %*% Diagonal(x = 1 / sums);S <- S * as.numeric(sums)
  live <- rowSums(S) > 0;H <- H[, live, drop = FALSE];S <- S[live, , drop = FALSE]
  list(H = H, S = S, alpha = A, mu = as.matrix(H %*% S) + B %*% A + O)
}

slcn_temporal_order <- function(model) {
  apex <- vapply(seq_len(ncol(model$H)), function(j) {
    at <- seq.int(model$H@p[j] + 1L, model$H@p[j + 1L])
    model$H@i[at[which.max(model$H@x[at])]] + 1L
  }, 0L)
  index <- order(apex)
  model$H <- model$H[, index, drop = FALSE]
  model$S <- model$S[index, , drop = FALSE]
  model$apex_scan <- apex[index]
  if (length(model$state_ids) == length(index)) model$state_ids <- model$state_ids[index]
  if (length(model$source_ids) == length(index)) model$source_ids <- model$source_ids[index]
  if (length(model$members) == length(index)) model$members <- model$members[index]
  labels <- if (length(index)) paste0("compound", seq_along(index)) else character()
  colnames(model$H) <- rownames(model$S) <- labels
  model
}

slcn_initial <- function(X, Y, B, spec, W, offset, scalar, ncores, screen_factor = 1) {
  p <- ncol(X);q <- ncol(B);pilot <- slcn_pilot(Y)
  variance <- switch(spec$name, gaussian = matrix(rep(rep_len(spec$dispersion, ncol(Y)), each = nrow(Y)), nrow(Y), ncol(Y)),
    poisson = pilot, quasipoisson = sweep(pilot, 2, rep_len(spec$dispersion, ncol(Y)), "*"),
    nb = pilot + sweep(pilot^2, 2, rep_len(spec$theta, ncol(Y)), "/"),
    nb1 = sweep(pilot, 2, rep_len(spec$dispersion, ncol(Y)), "*"),
    Gamma_scale = sweep(pilot, 2, rep_len(spec$dispersion, ncol(Y)), "*"),
    Gamma = sweep(pilot^2, 2, rep_len(spec$dispersion, ncol(Y)), "*"))
  working_weights <- matrix(pmin(1e12, pmax(1e-12, pilot / pmax(variance, 1e-20))), nrow(Y), ncol(Y)) * W
  design <- if (q) cbind(X, Matrix(B, sparse = TRUE)) else X
  fit <- L0glm.fit(design, pmax(Y, 0), family = poisson("identity"),
    weights = working_weights, offset = offset, lambda = scalar * screen_factor,
    nonnegative = TRUE, intercept = FALSE, always.include = p + seq_len(q),
    relaxed = TRUE, grouped = TRUE, group.maxit = if (spec$kind == 1L) 6L else 2L,
    normalize = "none", ncores = ncores)
  beta <- as.matrix(coef(fit));keep <- which(rowSums(beta[seq_len(p), , drop = FALSE]) > 0)
  H <- X[, keep, drop = FALSE];S <- beta[keep, , drop = FALSE]
  alpha <- beta[p + seq_len(q), , drop = FALSE]
  if (q && spec$kind != 0L) alpha <- pmax(alpha, 1e-8)
  mu <- as.matrix(H %*% S) + B %*% alpha + offset
  list(H = H, S = S, alpha = alpha, mu = mu, source_ids = keep)
}

slcn_run <- function(X, Y, B, spec, W, offset, scalar, ncores, control, start = NULL,
                     verbose = FALSE, censor = .5, rank_channels = ncol(Y), selection = TRUE) {
  objects <- slcn_context(Y, B, spec, W, offset, ncores, censor)
  engine <- objects$engine;ctx <- objects$context
  if (isTRUE(control$joint_trials)) slcn_enhance_trials(engine, ctx)
  penalty <- scalar * rank_channels
  d <- list(X = X, Y = Y, B = B, lambda = scalar, lambda_rank = penalty, psf_width = control$psf_width)
  elapsed <- numeric();history <- list()
  timer <- function(label, expression) {
    elapsed[label] <<- system.time(value <- force(expression))[["elapsed"]]
    value
  }
  current <- if (is.null(start)) timer("initial", slcn_initial(X, Y, B, spec, W, offset,
    scalar, ncores, control$screen_factor)) else start
  score <- function(model) slcn_loglik(Y, model$mu, ctx$spec$kind,
    ctx$spec$dispersion, ctx$spec$theta, ctx$W, ncores, censor)
  initial_rank <- ncol(current$H);initial_ll <- score(current)
  if (!is.finite(initial_ll)) stop("The initial mean has no finite family score; supply a positive background or offset.")
  history[[1L]] <- data.frame(stage = "initial", seconds = if (length(elapsed)) sum(elapsed) else 0,
    delta = 0, accepted = TRUE, compounds = ncol(current$H))
  step <- function(label, expression) {
    seconds <- system.time(candidate <- force(expression))[["elapsed"]]
    delta <- engine$spolo_delta_cpp(Y, current$mu, candidate$mu) + penalty * (ncol(candidate$H) - ncol(current$H))
    accepted <- is.finite(delta) && delta <= 1e-6
    history[[length(history) + 1L]] <<- data.frame(stage = label, seconds = seconds,
      delta = delta, accepted = accepted, compounds = ncol(candidate$H))
    if (accepted) current <<- candidate
    elapsed[label] <<- seconds
    if (verbose) {cat(label, round(seconds, 2), "s; change", signif(delta, 6),
      if (accepted) "accepted;" else "rejected;", ncol(current$H), "compounds\n");flush.console()}
  }
  fixed_refit <- function(model, sweeps = control$native_sweeps, reopen = FALSE) {
    K <- ncol(model$H);q <- ncol(B)
    if (q && ctx$spec$kind %in% 1:6) {
      beta <- if (reopen && K) pmax(model$S, .1 * rowMeans(model$S)) else model$S
      mean <- as.matrix(model$H %*% beta) + B %*% model$alpha + offset
      fit <- engine$spolo_spectral_cpp(Y, mean, model$H, t(beta), sweeps, ncores)
      model$S <- t(fit$S)
      fixed <- as.matrix(model$H %*% model$S) + offset
      background_fit <- ctx$background_fit(Y, fixed, model$alpha)
      model$alpha <- background_fit$coefficients
      fit <- engine$spolo_spectral_cpp(Y, fixed + background_fit$background,
        model$H, t(model$S), min(2L, sweeps), ncores)
      model$S <- t(fit$S);model$mu <- fit$mu
      live <- rowSums(model$S) > 0
      model$H <- model$H[, live, drop = FALSE];model$S <- model$S[live, , drop = FALSE]
      return(model)
    }
    design <- as(cbind(model$H, Matrix(B, sparse = TRUE)), "dgCMatrix")
    beta <- rbind(if (reopen && K) pmax(model$S, .1 * rowMeans(model$S)) else model$S,
      pmax(model$alpha, if (ctx$spec$name %in% c("Gamma", "Gamma_scale")) 1e-8 else 0))
    mean <- as.matrix(design %*% beta) + offset
    fit <- engine$spolo_spectral_cpp(Y, mean, design, t(beta), sweeps, ncores)
    result <- model
    result$S <- t(fit$S[, seq_len(K), drop = FALSE])
    result$alpha <- t(fit$S[, K + seq_len(q), drop = FALSE]);result$mu <- fit$mu
    live <- rowSums(result$S) > 0
    result$H <- result$H[, live, drop = FALSE];result$S <- result$S[live, , drop = FALSE]
    result
  }
  if (!is.null(start) || ctx$spec$kind != 1L)
    step("family_refit", fixed_refit(current, reopen = !is.null(start)))
  theta_update <- function(label) {
    if (ctx$spec$name != "nb" || !isTRUE(ctx$spec$theta_update)) return(invisible(NULL))
    old_theta <- ctx$spec$theta;old_ll <- score(current);accepted <- FALSE
    seconds <- system.time({
      for (iteration in 1:3) {
        derivative <- slcn_nb_score(Y, current$mu, ctx$spec$theta, ctx$W, ncores)
        move <- if (is.finite(derivative[2]) && derivative[2] < 0) -derivative[1] / derivative[2] else sign(derivative[1])
        move <- max(-1, min(1, move));reference <- score(current);theta0 <- ctx$spec$theta
        for (half in 0:12) {
          ctx$spec$theta <- max(1e-6, min(1e8, theta0 * exp(move / 2^half)))
          likelihood <- score(current)
          if (is.finite(likelihood) && likelihood >= reference) {accepted <- TRUE;break}
          ctx$spec$theta <- theta0
        }
        if (abs(log(ctx$spec$theta / theta0)) < 1e-4) break
      }
    })[["elapsed"]]
    delta <- -2 * (score(current) - old_ll)
    history[[length(history) + 1L]] <<- data.frame(stage = label, seconds = seconds,
      delta = delta, accepted = accepted, compounds = ncol(current$H))
    elapsed[label] <<- seconds
    if (verbose) {cat(label, "theta", signif(old_theta, 5), "->", signif(ctx$spec$theta, 5),
      round(seconds, 2), "s\n");flush.console()}
  }
  theta_update("pooled_theta_0")
  if (ncol(current$H) && selection) {
    step("fusion_1", engine$spolo_fast(current, d, ncores, control$native_sweeps, control$block_sweeps))
    for (pass in seq_len(control$broad_repeats)) {
      step(paste0("broad_", pass), engine$spolo_broad_polish(current, d, ncores))
      step(paste0("fusion_", pass + 1L), engine$spolo_fast(current, d, ncores,
        control$native_sweeps, control$block_sweeps, background = TRUE))
      theta_update(paste0("pooled_theta_", pass))
    }
    if (control$full_sweeps) step("local_sweep", engine$spolo_sweep(current, d, control$full_sweeps, ncores))
    step("completion", engine$spolo_complete(current, d, ncores, control$block_sweeps))
    for (pass in seq_len(control$final_polishes)) {
      step(paste0("final_fixed_", pass), fixed_refit(current))
      step(paste0("final_fusion_", pass), {
        state <- engine$spolo_state(current, d, penalty)
        engine$spolo_fuse(state);engine$spolo_prune(state)
        engine$spolo_export(state)
      })
    }
    theta_update("pooled_theta_final")
  }
  current$history <- do.call(rbind, history);current$spec <- ctx$spec
  current$seconds <- elapsed;current$loglik <- score(current)
  current$objective <- -2 * current$loglik + penalty * ncol(current$H)
  current$initial_objective <- -2 * initial_ll + penalty * initial_rank
  current$complete_delta <- current$objective - current$initial_objective
  current$ledger_error <- current$complete_delta - sum(current$history$delta[current$history$accepted])
  current
}

#' Sparse log-concave nonnegative matrix factorization
#'
#' Recover a parsimonious set of compounds from multichannel signals, with
#' nonnegative spectra, hard log-concave elution profiles and a smooth background.
#' Fast spectral fusion and local native updates accommodate broad, asymmetric
#' and shark-fin peaks while retaining all original outcome channels.
#'
#' @param X Nonnegative scans-by-candidates matrix, typically shifted copies of
#'   a point-spread function from [convmatrix()] or [L0glm_blind_deco()]. Columns
#'   must have contiguous log-concave positive support. Dense and sparse inputs
#'   are accepted and stored sparsely; columns are normalized to unit sum.
#' @param y Numeric scans-by-channels matrix or a single-channel vector.
#' @param family Identity-link noise family: `gaussian`, `poisson`,
#'   `quasipoisson`, `mgcv::nb`, `Gamma`, [nb1()], [Gamma_scale()] or
#'   [Gamma_censored()]. The corresponding names, including `"nb"`, `"nb1"`
#'   and `"Gamma_scale"`, are also accepted. Count likelihoods require integer
#'   counts; quasipoisson accepts nonnegative intensities. Gaussian observations
#'   may be negative; their fitted latent signals remain nonnegative.
#' @param background Optional nonnegative scans-by-basis matrix, for example
#'   `splines::bs(time, df=6, intercept=TRUE)`. Coefficients are nonnegative and
#'   unpenalized. `NULL` fits without a background.
#' @param lambda `"gic"` (default), `"aic"`, `"bic"`, `"ebic"`, `"mbic"`, or
#'   a nonnegative numeric compound penalty. See Details for its scale.
#' @param dispersion Known common Gaussian variance, quasipoisson scale,
#'   Gamma dispersion, NB1 variance-to-mean ratio or common Gamma scale.
#'   `NULL` estimates a pooled value from temporal differences and freezes it.
#'   Parameters supplied in a [signal_families] object are used when this
#'   argument is omitted. Poisson and NB2 use unit dispersion.
#' @param theta NB2 shape, shared by all original channels. `NULL` estimates and
#'   freezes a temporal moment estimate; a positive value fixes it; a negative
#'   value initializes at its absolute value and enables pooled likelihood
#'   updates. A fixed shape in the supplied NB family is respected.
#' @param precluster `NULL` or `FALSE` uses every channel throughout. `TRUE`
#'   uses up to 32 groups; an integer supplies another maximum. Size-adjusted
#'   Poisson clustering, aggregation, selection and fusion operate in that
#'   reduced space, followed by refinement on all original channels.
#' @param weights Nonnegative likelihood weights: a scalar, row vector or
#'   matrix matching `y`. `NULL` uses one. Zero weights omit observations.
#' @param offset Known nonnegative mean contribution, in the same supported
#'   formats as `weights`. `NULL` uses zero. It is included on the response
#'   scale because the link is identity.
#' @param correlation `"auto"` or `"lag1"` estimates an optional lag-one working
#'   adjustment to the selection penalty. `"none"` uses independent observations;
#'   a number between zero and 0.49 supplies the lag-one correlation.
#' @param ncores OpenMP workers. `NULL` uses the detected physical cores.
#' @param control Named list: `broad_repeats=1L`, `full_sweeps=1L`,
#'   `final_polishes=1L`, `native_sweeps=6L`, `block_sweeps=6L`,
#'   `psf_width=NULL`, `screen_factor=1`, `censor=0.5`, `linkage="ward"`,
#'   `joint_trials=NULL`, `profile_cycles=10L`, `fusion_passes=5L`,
#'   `profile_fusion=TRUE`. `joint_trials` uses direct Poisson fusion for Poisson data
#'   and joint local refinement for the other families. `psf_width` supplies
#'   the typical FWHM in scans. When omitted, local updates use the middle design
#'   column and broad-profile fusion uses the median of up to 21 representative
#'   design columns. `screen_factor` scales the initial selection
#'   penalty. `censor` is the Gamma detection limit: recorded zero intensities
#'   are left-censored below this value; positive values must be at least the
#'   limit. Zero requests an uncensored Gamma model. A supplied Gamma family
#'   constructor provides its own detection limit unless overridden here.
#'   `linkage` may be `"ward"` or `"bregman"`. `joint_trials` allows overlapping
#'   compounds to adjust together during local fusion refinement.
#'   `profile_cycles` bounds the original-channel shape/background refinement;
#'   zero disables it. `fusion_passes` bounds broad-anchor and compatible-tail
#'   fusion after that refinement.
#'   For a narrow, spectrally matching fragment contained in a broad profile,
#'   a guarded whole-profile refit can also update the shared background.
#'   Neighbouring modes stay within one quarter of their initial FWHM, and
#'   the merge is retained only when the original-channel criterion decreases.
#'   `profile_fusion` enables a final conservative fusion pass on all original
#'   channels, followed by one local likelihood-based profile/spectrum refit.
#' @param start Optional list `H`, `S` and `alpha`, with original dimensions.
#'   `H` has scans in rows and compounds in columns; `S` has compounds in rows
#'   and all original channels in columns. `alpha` contains background
#'   coefficients and defaults to zero. Profiles are normalized with the
#'   inverse change applied to spectra. A restart bypasses initial screening
#'   and preclustering; nonnegative spectral entries can reopen during refitting.
#' @param verbose Print stage timings, objective changes and compound counts.
#' @param object,x A fitted `splogconNMF` object.
#' @param ... Further method arguments.
#' @param k Complexity multiplier for [stats::AIC()].
#'
#' @details
#' **Generalized matrix factorization.** The model is an identity-link,
#' nonnegative generalized matrix factorization (GMF), with a smooth background,
#' log-concave profiles and automatic rank selection. See Kidzinski et al.
#' (2022) and Castiglione et al. (2026). The function retains the NMF name;
#' [the mathematical article](../articles/log-concave-nmf.html#generalized-matrix-factorization)
#' relates its mean model to the GMF formulation.
#'
#' **Deviance matrix factorization.** This is a rank-penalized, nonnegative
#' deviance matrix factorization on the mean scale, with hard log-concavity
#' constraints and an additive background. For standard exponential families
#' with fixed noise parameters, the loss is a weighted Bregman divergence:
#' Gaussian gives squared error, Poisson gives generalized KL, and uncensored
#' Gamma gives Itakura--Saito. See Wang and Carvalho (2023), and Li and Tyagi
#' (2023). The [mathematical article](../articles/log-concave-nmf.html#deviance-matrix-factorization-and-bregman-nmf)
#' gives the derivation and the negative-binomial case.
#'
#' **Model and selection criterion.** With scans in rows, write
#' \deqn{\mu=HS+B A+O,\qquad H,S,A\geq0.}
#' Each profile has unit sum, contiguous positive support and
#' \eqn{\Delta^2\log H_{\cdot r}\leq0} on that support. The objective is
#' \deqn{J=-2\sum_{ij}w_{ij}\ell(y_{ij};\mu_{ij},\nu)+\lambda_R K,}
#' where \eqn{K} is the number of compounds and \eqn{\nu} contains the shared
#' noise parameters. An IC preset uses \eqn{\lambda_R=m\kappa\lambda_{IC}},
#' with \eqn{m} original channels. A numeric `lambda` supplies the total
#' compound penalty before multiplication by \eqn{\kappa}. The initial
#' grouped fit uses \eqn{\lambda_R/m}. Each compound contributes one selection
#' count; spectral entries and profile complexity have no separate penalty.
#' This is a compound-count model-selection convention. The IC methods use
#' \eqn{mK} plus the number of estimated shared noise parameters as their
#' complexity count; it is distinct from counting every free profile ordinate.
#'
#' **Initialization and local updates.** A variance-matched grouped Poisson
#' working fit supplies candidate components, using [L0glm.fit()] with relaxed
#' refitting and compulsory background columns. The subsequent updates use
#' the specified family's exact score and likelihood. For a spectral
#' coefficient or profile ordinate, the conditional mean is
#' \eqn{\mu_i=r_i+x_i b}, with the other components fixed in \eqn{r_i}.
#' If \eqn{g_i} and \eqn{I_i>0} are the negative log-likelihood score and a
#' positive working curvature, a nonnegative quadratic proposal is
#' \deqn{b^+=\max\left\{0,b-\frac{\sum_i w_i x_i g_i}
#' {\sum_i w_i x_i^2 I_i}\right\}.}
#' Backtracking checks the actual loss, with an explicit zero-boundary check.
#' Poisson and quasipoisson use a specialized bracketed scalar likelihood solve.
#' Background terms use small constrained quadratic solves, and all channels
#' run independently through OpenMP.
#' Poisson profile proposals use four warm native IRLS iterations. During a
#' full profile sweep, an adaptive ridge tether centred on the preceding
#' profile stabilizes weak tails, with local penalty weights proportional to
#' the inverse squared profile ordinate. The final acceptance test uses the
#' complete likelihood-and-rank criterion.
#'
#' Spectral clustering proposes fusions without fitting a model. It compares
#' normalized spectra through Jensen--Shannon information loss and joins
#' overlapping compatible profiles. A least concave majorant of the log profile
#' enforces hard log concavity. Local refinement reallocates signal between
#' overlapping compounds and checks the actual criterion. Every complete stage
#' is accepted against its full working data, and final refinement and scoring
#' use all original observations. `history` records each accepted or rejected
#' stage; `ledger_error` checks the accumulated criterion change.
#'
#' **Broad-profile fusion.** After refitting the original channels, an additional
#' algebraic pass compares elution profiles. Candidate pairs must have cosine
#' similarity at least 0.95 and at least one profile with FWHM at least four
#' times the typical PSF FWHM. This can reunite complementary ion sets assigned
#' to nearly proportional broad profiles. The area-weighted profile is projected
#' to hard log concavity and the spectra are added. Each merge must improve
#' the original-channel likelihood-and-compound criterion with the noise scale,
#' weights, offsets and background fixed. Candidate fusion requires no refits.
#' One final local likelihood pass updates the profiles and spectra in affected
#' elution regions. It uses Fisher scoring in log-profile coordinates with
#' piecewise-linear interpolation and linear concavity constraints; the noise
#' likelihood is evaluated on the original intensity scale. Backtracking
#' preserves likelihood improvement. Neighboring profile modes stay within
#' one quarter of their starting FWHM, while merged profiles can move freely.
#' In log-profile coordinates \eqn{\eta=\log H_{\cdot r}}, the score and working
#' information are \eqn{g_t=H_{tr}\sum_c w_{tc}S_{rc}g_{\mu,tc}} and
#' \eqn{I_t=H_{tr}^2\sum_c w_{tc}S_{rc}^2I_{\mu,tc}}. With interpolation matrix
#' \eqn{C}, the quadratic step minimizes
#' \eqn{\sum_t I_t\{(Cb)_t-\eta_t+g_t/I_t\}^2}, subject to nonincreasing
#' slopes of \eqn{Cb}. Zero-information rows receive numerical stabilization.
#' `profile_fusion` records the merge audit and source component groups. Source
#' IDs refer to components entering original-channel profile refinement, or the
#' final pass when that refinement is disabled. Groups are retained across
#' passes; `source_components` follows the returned temporal ordering.
#'
#' **Noise models.** Gaussian variance is \eqn{\phi}; quasipoisson variance is
#' \eqn{\phi\mu}; NB2 variance is \eqn{\mu+\mu^2/\theta}; standard Gamma
#' variance is \eqn{\phi\mu^2}. NB1 has variance \eqn{\phi\mu}, while
#' common-scale Gamma has shape \eqn{\mu/b} and scale \eqn{b}, hence variance
#' \eqn{b\mu}. Their exact scores and densities are given in [signal_families].
#' A censored Gamma zero contributes a log CDF at the detection limit, and a
#' positive observation contributes its log density. For quasipoisson,
#' \eqn{\ell_Q=\{y\log\mu-\mu-\log\Gamma(y+1)\}/\phi}; `logLik()` marks
#' this scaled kernel with `quasi=TRUE`, and the regular IC method names return
#' quasi-information criteria. Compare candidate ranks under the same family
#' and frozen scale.
#'
#' For an updated NB2 shape, the pooled score is
#' \deqn{U_\theta=\sum_{ij}w_{ij}\left\{\psi(y_{ij}+\theta)-\psi(\theta)
#' +\log\theta+1-\log(\theta+\mu_{ij})
#' -\frac{y_{ij}+\theta}{\theta+\mu_{ij}}\right\}.}
#' Bounded Newton proposals in \eqn{\log\theta} are backtracked against the
#' pooled likelihood. One common shape is retained across original channels.
#'
#' **Temporal calibration and aggregation.** [temporal_dispersion()] uses
#' fourth and eighth differences in quieter parts of a smoothed pilot signal.
#' This estimates the marginal noise scale; an optional lag-one working
#' correction multiplies the rank penalty by
#' \eqn{\kappa=1+2\widehat\rho_1}. It leaves the fitted variance parameter on
#' its marginal scale. This local covariance approximation is recorded in
#' `spec`; use `correlation="none"` for independent simulated data.
#'
#' Preclustering uses [glm_cluster()] with `size_fac=colSums(pmax(y,0))`,
#' removing differences in total ion intensity when comparing channel shapes.
#' Poisson and common-parameter NB1 latent signals add exactly; independent
#' common-scale Gamma latent signals also add. Standard NB2 and Gamma use
#' moment-matched aggregate parameters. Censored Gamma sums and differing
#' channel weights give working approximations in the reduced space. The
#' returned model always contains all original channels, their offsets and
#' their original-channel likelihood.
#' The compound penalty uses the original channel count throughout aggregation.
#' Selection and initial fusion run in the reduced space. Lifting the spectra
#' to the original channels is followed by a spectral/background refit and
#' local shape refinement, with no repetition of the reduced-space search.
#' Worked full-channel fits to [GCMS_data1] and [GCMS_data2], with traces,
#' heatmaps, profiles and spectra, are shown in the
#' [GC/MS deconvolution tutorial](../articles/gcms-deconvolution.html) and
#' the [GCMSdeco()] help page. [write_GCMS_cdf()] exports these matrices for
#' the complete file-based workflow.
#'
#' **Smooth local profiles.** Original-channel refinement fits
#' \eqn{\eta(t)=a+bt-\tfrac12\sum_j c_j\{(t-u_j)_+^2-(t-v_j)_+^2\}},
#' with \eqn{c_j\geq0}. Thus \eqn{\eta''(t)\leq0} and
#' \eqn{H(t)=\exp\{\eta(t)\}} is log-concave. Fisher scoring uses the
#' original observation model and the chain-rule score/information given above.
#' Small constrained quadratic solves and likelihood backtracking update each
#' local profile. Joint spectral updates separate overlapping compounds;
#' compatible broad-profile and tail fusions are accepted under the same
#' original-channel likelihood and compound-count penalty.
#'
#' @return An object of class `splogconNMF`, containing sparse unit-sum profiles
#'   `H`, spectra `S`, background coefficients `alpha`, the background basis,
#'   fitted means `mu`, `rank`, resolved noise settings `spec`, original data,
#'   offsets and weights, `objective`, `loglik`, stage ledger `history`,
#'   `ledger_error`, aggregation and `profile_fusion` diagnostics, controls and
#'   timings `seconds`.
#'   `coef()` returns spectra and `fitted()` the full mean matrix. `logLik()`,
#'   `AIC()`, `BIC()`, `GIC()`, `EBIC()` and `MBIC()` score the returned model.
#'   Columns of `H` and the matching rows of `S` are ordered by increasing apex
#'   scan, also returned in `apex_scan`. Their labels identify the same compounds.
#' @seealso [GCMSdeco()], [L0glm_blind_deco()], [L0glm.fit()], [glm_cluster()]
#' @references Wang, L. and Carvalho, L. E. (2023). Deviance matrix factorization.
#'   Electronic Journal of Statistics 17(2), 3762--3810. \doi{10.1214/23-EJS2174}.
#'   Li, X. and Tyagi, A. (2023). Block-active ADMM to minimize NMF with Bregman
#'   divergences. Sensors 23(16), 7229. \doi{10.3390/s23167229}.
#'   Kidzinski, L., Hui, F. K. C., Warton, D. I. and Hastie, T. J. (2022).
#'   Generalized matrix factorization: efficient algorithms for fitting
#'   generalized linear latent variable models to large data arrays.
#'   JMLR 23(291), 1--29. \url{https://www.jmlr.org/papers/v23/20-1104.html}.
#'   Castiglione, C., Segers, A., Clement, L. and Risso, D. (2026).
#'   Stochastic gradient descent estimation of generalized matrix factorization
#'   models with application to single-cell RNA sequencing data.
#'   Biostatistics 27(1), kxag010. \doi{10.1093/biostatistics/kxag010}.
#' @example examples/splogconNMF.R
#' @md
#' @export
splogconNMF <- function(X, y, family = poisson("identity"), background = NULL,
                       lambda = "gic", dispersion = NULL, theta = NULL,
                       precluster = NULL, weights = NULL, offset = NULL,
                       correlation = c("auto", "none", "lag1"), ncores = NULL,
                       control = list(), start = NULL, verbose = FALSE) {
  begun <- proc.time()[["elapsed"]]
  correlation <- slcn_correlation(correlation)
  if (is.null(ncores)) {
    ncores <- parallel::detectCores(logical = FALSE)
    if (!is.finite(ncores)) ncores <- 1L
  }
  if (length(ncores) != 1L || !is.finite(ncores) || ncores < 1 || ncores != as.integer(ncores)) stop("ncores must be a positive integer.")
  spec <- slcn_family(family, dispersion, theta)
  Y <- as.matrix(y)
  if (!is.numeric(Y)) stop("y must be a numeric observation matrix.", call. = FALSE)
  storage.mode(Y) <- "double";n <- nrow(Y);m <- ncol(Y)
  if (n < 9L || !m || any(!is.finite(Y))) stop("y must be finite with at least nine scans and one channel.")
  if (spec$kind != 0L && any(Y < 0)) stop("This noise family requires nonnegative observations.")
  if (spec$name %in% c("poisson", "nb", "nb1") && any(abs(Y - round(Y)) > 1e-7)) stop("Poisson and NB observations must be integer counts.")
  X <- slcn_design(X, n, "X")
  if (!ncol(X)) stop("X must have at least one candidate profile.")
  X <- X %*% Diagonal(x = 1 / colSums(X))
  B <- if (is.null(background)) matrix(0, n, 0L) else as.matrix(background)
  if (nrow(B) != n || any(!is.finite(B) | B < 0)) stop("background must be a finite nonnegative basis.")
  W <- slcn_weights(weights, n, m);O <- slcn_weights(offset, n, m, default = 0)
  if (any(!is.finite(W) | W < 0) || !any(W > 0) || any(!is.finite(O) | O < 0)) stop("Nonnegative finite weights/offset and some positive weight are required.")
  defaults <- list(broad_repeats = 1L, full_sweeps = 1L, final_polishes = 1L,
    native_sweeps = 6L, block_sweeps = 6L, psf_width = NULL, screen_factor = 1,
    censor = .5, linkage = "ward", joint_trials = NULL, profile_fusion = TRUE,
    profile_cycles = 10L, fusion_passes = 5L)
  if (!is.list(control) || length(control) && (is.null(names(control)) ||
      any(!nzchar(names(control))) || anyDuplicated(names(control)))) stop("Use a named control list.")
  if (is.list(family) && !is.null(family$scalar) && !"censor" %in% names(control))
    defaults$censor <- family$scalar$censor
  if (length(setdiff(names(control), names(defaults)))) stop("Unknown control entry.")
  control <- utils::modifyList(defaults, control, keep.null = TRUE)
  supplied_width <- !is.null(control$psf_width)
  if (is.null(control$joint_trials)) control$joint_trials <- spec$kind != 1L
  control$linkage <- match.arg(control$linkage, c("ward", "bregman"))
  if (is.null(control$psf_width)) {
    h <- as.numeric(X[, max(1L, ncol(X) %/% 2L)]);control$psf_width <- sum(h >= max(h) / 2)
  }
  for (key in c("broad_repeats", "full_sweeps", "final_polishes", "native_sweeps", "block_sweeps", "psf_width", "profile_cycles", "fusion_passes")) {
    value <- control[[key]];minimum <- if (key %in% c("broad_repeats", "full_sweeps", "final_polishes", "profile_cycles", "fusion_passes")) 0 else 1
    if (length(value) != 1L || !is.finite(value) || value < minimum || value != as.integer(value)) stop("Invalid integer control: ", key)
  }
  if (!is.numeric(control$screen_factor) || length(control$screen_factor)!=1L ||
      !is.finite(control$screen_factor) || control$screen_factor <= 0 ||
      !is.numeric(control$censor) || length(control$censor)!=1L ||
      !is.finite(control$censor) || control$censor < 0) stop("Use a positive screen_factor and nonnegative censor.")
  if (!is.logical(control$joint_trials) || length(control$joint_trials)!=1L || is.na(control$joint_trials))
    stop("joint_trials must be TRUE or FALSE.")
  if (!is.logical(control$profile_fusion) || length(control$profile_fusion)!=1L || is.na(control$profile_fusion))
    stop("profile_fusion must be TRUE or FALSE.")
  fusion_width <- control$psf_width
  if (control$profile_fusion && !supplied_width) {
    columns <- unique(round(seq(1, ncol(X), length.out = min(ncol(X), 21L))))
    fusion_width <- stats::median(vapply(columns, function(j) {
      h <- X@x[seq.int(X@p[j] + 1L, X@p[j + 1L])];sum(h >= max(h) / 2)
    }, 0L))
  }
  if (spec$name %in% c("Gamma", "Gamma_scale") && any(ifelse(Y == 0, control$censor == 0, Y < control$censor)))
    stop("Represent observations below the Gamma detection limit by zero.")
  if (!is.null(start)) start <- slcn_start(start, n, m, B, O)
  if (identical(precluster, TRUE)) precluster <- 32L
  if (identical(precluster, FALSE)) precluster <- NULL
  if (!is.null(precluster) && (length(precluster) != 1L || !is.finite(precluster) || precluster < 1 || precluster != as.integer(precluster))) stop("precluster must be NULL, logical or a positive integer.")
  if (length(lambda)!=1L || !(is.character(lambda) || is.numeric(lambda)))
    stop("lambda must be an IC name or one nonnegative compound penalty.")
  nuisance <- system.time(noise <- slcn_noise(Y, spec, W, ncores, correlation, control$censor))[["elapsed"]]
  spec <- noise$spec
  scalar <- if (is.character(lambda)) preset.lambda(lambda, p = ncol(X), n = n) else lambda / m
  scalar <- scalar * spec$information_factor
  if (!is.finite(scalar) || scalar < 0) stop("lambda must be an IC preset or nonnegative compound penalty.")
  clustering_seconds <- aggregate_seconds <- 0;aggregation <- NULL
  if (is.null(start) && !is.null(precluster) && !identical(precluster, FALSE) && precluster < m) {
    clustering_seconds <- system.time(cl <- glm_cluster(pmax(Y, 0), family = poisson("identity"),
      weights = W, method = "poissondist", size_fac = colSums(pmax(Y, 0)),
      rowbins = if (all(W == 1)) 5L else 1L,
      k = as.integer(precluster), linkage = control$linkage, nthreads = ncores))[["elapsed"]]
    group <- cl$clusterids;g <- max(group)
    P <- sparseMatrix(i = seq_len(m), j = group, x = 1, dims = c(m, g))
    aggregate_Y <- as.matrix(Y %*% P)
    fractions <- colSums(pmax(Y, 0)) / as.numeric(colSums(pmax(Y, 0) %*% P))[group]
    fractions[!is.finite(fractions)] <- 0
    sumsq <- as.numeric(rowsum(fractions^2, group, reorder = TRUE))
    sizes <- tabulate(group, nbins = g)
    aggregated_spec <- spec;aggregated_spec$theta_update <- FALSE
    if (spec$name == "gaussian") aggregated_spec$dispersion <- spec$dispersion * sizes
    if (spec$name == "nb") aggregated_spec$theta <- spec$theta / pmax(sumsq, 1e-12)
    if (spec$name == "Gamma") aggregated_spec$dispersion <- spec$dispersion * pmax(sumsq, 1e-12)
    aggregate_W <- if (length(W) == 1L) W else as.matrix(W %*% P) / rep(sizes, each = n)
    aggregate_O <- if (length(O) == 1L) matrix(O, n, g) * rep(sizes, each = n) else as.matrix(O %*% P)
    aggregate_seconds <- system.time(a <- slcn_run(X, aggregate_Y, B, aggregated_spec,
      aggregate_W, aggregate_O, scalar, ncores, control, verbose = verbose, censor = control$censor,
      rank_channels = m))[["elapsed"]]
    S <- sweep(a$S[, group, drop = FALSE], 2, fractions, "*")
    alpha <- sweep(a$alpha[, group, drop = FALSE], 2, fractions, "*")
    start <- list(H = a$H, S = S, alpha = alpha,
      mu = as.matrix(a$H %*% S) + B %*% alpha + O)
    aggregation <- list(group = group, fractions = fractions, channels = g,
      moment_matched = spec$name %in% c("nb", "Gamma"),
      censored_approximation = spec$name %in% c("Gamma", "Gamma_scale"),
      weighted_approximation = length(W) != 1L, model = a,
      sumsq_fractions = sumsq, original_channels = m)
  }
  fit_seconds <- system.time(fit <- slcn_run(X, Y, B, spec, W, O, scalar, ncores, control,
    start, verbose, control$censor, selection = is.null(aggregation)))[["elapsed"]]
  fit$background <- B;fit$offset <- O;fit$weights <- W;fit$y <- Y
  fit$family <- family;fit$n <- n;fit$p <- ncol(X);fit$lambda <- scalar * m
  fit$rank <- ncol(fit$H)
  fit$nuisance_df <- as.integer(spec$dispersion_estimated) + as.integer(spec$theta_estimated)
  fit$ic_df <- m * fit$rank + fit$nuisance_df
  fit$aggregation <- aggregation;fit$control <- control;fit$ncores <- ncores
  fit <- slcn_temporal_order(fit)
  refinement_seconds <- system.time(
    if (ncol(fit$H) && control$profile_cycles) fit <- slcn_refine_fit(fit, Y, verbose))[["elapsed"]]
  if (control$profile_fusion) refinement_seconds <- refinement_seconds + system.time(
    fit <- slcn_refine_profiles(fit, Y, verbose, fusion_width))[["elapsed"]]
  fit$rank <- ncol(fit$H);fit$ic_df <- m * fit$rank + fit$nuisance_df
  fit$seconds <- c(noise = nuisance, clustering = clustering_seconds, aggregate = aggregate_seconds,
    original = fit_seconds, refinement = refinement_seconds, total = proc.time()[["elapsed"]] - begun)
  if (!is.null(fit$profile_fusion$merges) && nrow(fit$profile_fusion$merges))
    fit <- slcn_temporal_order(fit)
  if (!is.null(fit$profile_fusion)) fit$profile_fusion$source_components <- fit$members
  fit$members <- NULL
  if (!is.null(fit$aggregation$model))
    fit$aggregation$model <- slcn_temporal_order(fit$aggregation$model)
  colnames(fit$S) <- colnames(Y);dimnames(fit$mu) <- dimnames(Y)
  class(fit) <- "splogconNMF"
  fit
}

#' @rdname splogconNMF
#' @export
fitted.splogconNMF <- function(object, ...) {
  if(!is.null(object$mu)) return(object$mu)
  as.matrix(object$H %*% object$S)+object$background %*% object$alpha+object$offset
}
#' @rdname splogconNMF
#' @export
nobs.splogconNMF <- function(object, ...) object$n
#' @rdname splogconNMF
#' @export
coef.splogconNMF <- function(object, ...) object$S
slcn_df <- function(object) ncol(object$S)*ncol(object$H) +
  as.integer(object$spec$dispersion_estimated) + as.integer(object$spec$theta_estimated)
#' @rdname splogconNMF
#' @export
logLik.splogconNMF <- function(object, ...) structure(object$loglik, class = "logLik",
  df = slcn_df(object), nobs = object$n, quasi = object$spec$name == "quasipoisson")
#' @rdname splogconNMF
#' @export
AIC.splogconNMF <- function(object, ..., k = 2) -2 * object$loglik +
  k * object$spec$information_factor * slcn_df(object)
#' @rdname splogconNMF
#' @export
BIC.splogconNMF <- function(object, ...) AIC(object, k = log(object$n))
#' @rdname splogconNMF
#' @export
GIC.splogconNMF <- function(object, ...) AIC(object,
  k = preset.lambda("gic", p = object$p, n = object$n))
#' @rdname splogconNMF
#' @export
EBIC.splogconNMF <- function(object, ...) AIC(object,
  k = preset.lambda("ebic", p = object$p, n = object$n))
#' @rdname splogconNMF
#' @export
MBIC.splogconNMF <- function(object, ...) AIC(object,
  k = preset.lambda("mbic", p = object$p, n = object$n))
#' @rdname splogconNMF
#' @method print splogconNMF
#' @export
print.splogconNMF <- function(x, ...) {
  cat("Sparse log-concave NMF (", x$spec$name, ")\n", ncol(x$H), " compounds; ",
    ncol(x$S), " channels; ", round(x$seconds[["total"]], 2), " s\n",
    "dispersion = ", signif(x$spec$dispersion[1L], 5), "; theta = ", signif(x$spec$theta[1L], 5),
    "; lag-one correlation = ", signif(x$spec$correlation, 4), "\n",
    "GIC = ", format(GIC(x), digits = 9), "\n", sep = "")
  invisible(x)
}
