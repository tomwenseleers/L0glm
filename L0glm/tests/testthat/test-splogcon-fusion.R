fusion_test_model <- function(kind = 1L, distinct = FALSE, precision = 1) {
  n <- 160L;m <- 4L
  profile <- function(centre) {
    h <- numeric(n);rows <- seq.int(centre - 25L, centre + 25L)
    h[rows] <- exp(-(rows - centre)^2 / 50);h / sum(h)
  }
  H <- Matrix::Matrix(cbind(profile(70L), profile(if (distinct) 72L else 70L)), sparse = TRUE)
  S <- rbind(c(80000, 40000, 0, 0), c(0, 0, 100000, 30000))
  B <- matrix(1, n, 1L);A <- matrix(3, 1L, m);O <- matrix(.2, n, m)
  mu <- as.matrix(H %*% S) + B %*% A + O
  W <- matrix(precision, n, m);W[1:4, ] <- 0;W[, 2L] <- W[, 2L] / 2
  phi <- switch(as.character(kind), `0` = 1, `2` = 10, `4` = 1, `5` = 5, `6` = 5, 1)
  model <- list(H = H, S = S, background = B, alpha = A, offset = O, mu = mu,
    weights = W, spec = list(kind = kind, dispersion = phi, theta = 20), lambda = 1000,
    control = list(censor = .5, psf_width = 2L), ncores = 1L, state_ids = 1:2)
  Y <- if (kind %in% c(1L, 3L, 5L)) round(mu) else mu
  list(f = model, Y = Y)
}

test_that("broad complementary fragments merge under each original-channel loss", {
  for (kind in 0:6) {
    d <- fusion_test_model(kind)
    z <- slcn_fuse_profiles(d$f, d$Y, 2L)
    expect_equal(dim(z$H), c(160L, 1L));expect_equal(dim(z$S), c(1L, 4L))
    expect_equal(z$mu, d$f$mu, tolerance = 1e-12)
    expect_equal(as.numeric(z$S), as.numeric(colSums(d$f$S)))
    expect_equal(z$profile_fusion$merges$delta, -1000, tolerance = 1e-8)
    expect_identical(z$members, list(1:2));expect_null(z$state_ids)
    expect_identical(z$alpha, d$f$alpha);expect_identical(z$offset, d$f$offset)
    expect_identical(z$weights, d$f$weights)
    expect_identical(z$profile_fusion$fitting_calls_in_fusion, 0L)
    # Close but resolvably shifted profiles must retain their distinct spectra.
    d <- fusion_test_model(kind, distinct = TRUE, precision = 1000)
    z <- slcn_fuse_profiles(d$f, d$Y, 2L)
    expect_identical(z$H, d$f$H);expect_identical(z$S, d$f$S)
    expect_equal(nrow(z$profile_fusion$merges), 0L)
  }
})

test_that("ordinary narrow profiles are protected and empty models are supported", {
  d <- fusion_test_model()
  z <- slcn_fuse_profiles(d$f, d$Y, 10L)
  expect_identical(z$H, d$f$H);expect_identical(z$S, d$f$S)
  expect_equal(nrow(z$profile_fusion$merges), 0L)
  d$f$H <- d$f$H[, FALSE, drop = FALSE];d$f$S <- d$f$S[FALSE, , drop = FALSE]
  z <- slcn_fuse_profiles(d$f, d$Y, 2L)
  expect_equal(dim(z$H), c(160L, 0L));expect_equal(dim(z$S), c(0L, 4L))
})

test_that("the final profile pass retains lineage from contained-tail refinement", {
  d <- fusion_test_model()
  d$f$members <- list(c(2L, 4L), c(7L, 9L))
  z <- slcn_fuse_profiles(d$f, d$Y, 2L)
  expect_identical(z$members, list(c(2L, 4L, 7L, 9L)))
  expect_identical(z$profile_fusion$merges$source, "2;4;7;9")
})

test_that("guarded tail refinement retains dimensions and every original noise model", {
  for (kind in 0:6) {
    d <- fusion_test_model(kind)
    f <- d$f
    f$members <- list(1L, 2L)
    before <- -2 * slcn_loglik(d$Y, f$mu, kind, f$spec$dispersion,
      20, as.numeric(f$weights), 1L, .5) + f$lambda * 2
    proposal <- list(group = 1:2, anchor = 1L, donors = 2L, shape = "hull")
    z <- slcn_merge_anchor(f, d$Y, proposal, cycles = 2L, guarded = TRUE)
    after <- -2 * slcn_loglik(d$Y, z$mu, kind, f$spec$dispersion,
      20, as.numeric(f$weights), 1L, .5) + f$lambda
    expect_lte(after, before + 1e-6)
    expect_equal(dim(z$H), c(160L, 1L))
    expect_equal(dim(z$S), c(1L, 4L))
    expect_equal(z$mu, as.matrix(z$H %*% z$S) + z$background %*% z$alpha + z$offset,
      tolerance = 1e-8)
    expect_identical(z$offset, f$offset)
    expect_identical(z$weights, f$weights)
    expect_equal(as.numeric(Matrix::colSums(z$H)), 1, tolerance = 1e-10)
    h <- as.numeric(z$H[, 1L]);h <- h[h > 0]
    expect_lte(max(diff(log(h), differences = 2L)), 1e-6)
  }
})

test_that("the public pass uses a representative PSF width and respects overrides", {
  d <- fusion_test_model();n <- nrow(d$Y)
  X <- convmatrix(convkernel("gaussian", 9L, list(sigma=.5)), n)
  wide <- numeric(n);rows <- 40:120
  wide[rows] <- exp(-(rows-80)^2/(2*12^2));X[,80L] <- wide/sum(wide)
  start <- d$f[c("H","S","alpha")]
  a <- list(X=X,y=d$Y,family=poisson("identity"),background=d$f$background,
    offset=d$f$offset,weights=d$f$weights,lambda=1000,correlation="none",ncores=2L,start=start)
  controls <- list(broad_repeats=0L,full_sweeps=0L,final_polishes=0L)
  z <- do.call(splogconNMF,c(a,list(control=controls)))
  expect_equal(z$profile_fusion$minimum_FWHM,4)
  expect_gt(z$control$psf_width,1)
  expect_true("broad_profile_fusion" %in% z$history$stage)
  expect_equal(length(z$profile_fusion$source_components),ncol(z$H))
  controls$psf_width <- 29L
  b <- do.call(splogconNMF,c(a,list(control=controls)))
  expect_equal(b$profile_fusion$minimum_FWHM,116)
  expect_lte(ncol(z$H),ncol(b$H))
  controls$profile_fusion <- FALSE
  off <- do.call(splogconNMF,c(a,list(control=controls)))
  expect_null(off$profile_fusion)
  expect_false("broad_profile_fusion" %in% off$history$stage)
})

test_that("native constrained log profiles improve the original likelihood in parallel", {
  d <- fusion_test_model();rows <- 45:95
  h <- as.numeric(d$f$H[rows, 1L]);s <- colSums(d$f$S)
  mean <- .2 + matrix(3, length(rows), 4L) + h %o% s
  for (kind in 0:6) {
    d <- fusion_test_model(kind)
    call <- function(threads) slcn_log_profile(d$Y[rows, , drop = FALSE], mean, h, s,
      kind, d$f$spec$dispersion, 20, as.numeric(d$f$weights), 160L, 44L, 3L, .5,
      48L, threads, 18L, 32L)
    a <- call(1L);b <- call(2L)
    expect_identical(a, b)
    expect_true(all(is.finite(a$profile) & a$profile > 0))
    expect_lte(max(diff(log(a$profile), differences = 2L)), 1e-7)
    expect_true(all(a$history <= 1e-7))
    after <- mean - h %o% s + a$profile %o% s
    change <- slcn_delta(d$Y[rows, , drop = FALSE], mean, after, kind,
      d$f$spec$dispersion, 20, as.numeric(d$f$weights), 160L, 44L, .5)
    expect_lte(change, 1e-6)
  }
  expect_error(slcn_log_profile(d$Y[rows, ], mean, h * 0, s, 1L, 1, 1,
    1, length(rows)), "Positive finite")
})

test_that("the final local refit preserves constraints, background and the criterion ledger", {
  for (kind in 0:6) {
    d <- fusion_test_model(kind);f <- d$f
    f$initial_objective <- -2 * slcn_loglik(d$Y, f$mu, kind, f$spec$dispersion,
      20, as.numeric(f$weights), 1L, .5) + 2000
    f$history <- data.frame(stage = "initial", seconds = 0, delta = 0,
      accepted = TRUE, compounds = 2L)
    z <- slcn_refine_profiles(f, d$Y)
    expect_equal(ncol(z$H), 1L)
    expect_equal(z$mu, as.matrix(z$H %*% z$S) + f$background %*% f$alpha + f$offset,
      tolerance = 1e-8)
    expect_identical(z$alpha, f$alpha);expect_identical(z$offset, f$offset)
    expect_equal(as.numeric(Matrix::colSums(z$H)), 1, tolerance = 1e-10)
    h <- as.numeric(z$H[, 1L]);rows <- which(h > 0)
    expect_true(all(diff(rows) == 1L))
    expect_lte(max(diff(log(h[rows]), differences = 2L)), 1e-6)
    expect_true(all(z$history$delta[z$history$accepted] <= 1e-6))
    expect_lte(z$objective, f$initial_objective)
    expect_lt(abs(z$ledger_error), 1e-3)
  }
})
