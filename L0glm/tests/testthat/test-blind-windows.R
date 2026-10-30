test_that("time-varying PSF interpolation preserves convolution and log concavity", {
  build <- getFromNamespace(".psf_build_sparse_matrix", "L0glm")
  n <- 160L; L <- 31L
  for (g in list(convkernel("gaussian", L, list(sigma = 2)),
                 log_conc_smooth(convkernel("bigaussian", L, list(sigma1 = 1, sigma2 = 3))))) {
    X <- build(seq_len(n), c(30, 130), seq_len(L) - (L + 1L) %/% 2L,
      rbind(g, g), log_floor = 0, drop_tol = 1e-12, normalize = FALSE, smooth = TRUE)
    expect_equal(as.matrix(X), as.matrix(convmatrix(g, n)), tolerance = 1e-11)
  }
  kernels <- rbind(convkernel("gaussian", L, list(sigma = 1.5)),
                   convkernel("gaussian", L, list(sigma = 4)))
  X <- build(seq_len(n), c(30, 130), seq_len(L) - (L + 1L) %/% 2L,
    kernels, log_floor = 0, drop_tol = 1e-12, normalize = FALSE, smooth = TRUE)
  expect_s4_class(X, "dgCMatrix")
  for (j in seq_len(n)) {
    h <- as.numeric(X[, j]); h <- h[h > 0]
    expect_true(all(diff(log(h), differences = 2L) < 1e-7))
  }
})

test_that("window fits learn varying PSFs and refit every original channel", {
  set.seed(211)
  n <- 1200L; L <- 31L
  centers <- c(200.5, 600.5, 1000.5)
  kernels <- t(vapply(c(1.5, 2.5, 3.5), function(s)
    convkernel("gaussian", L, list(sigma = s)), numeric(L)))
  X <- L0glm:::.psf_build_sparse_matrix(seq_len(n), centers, seq_len(L) - 16L,
    kernels, log_floor = 0, drop_tol = 1e-12, normalize = FALSE, smooth = TRUE)
  A <- matrix(0, n, 3)
  at <- seq.int(45, n - 45L, by = 40L)
  A[at, ] <- matrix(runif(length(at) * 3, 10000, 80000), length(at), 3)
  mu <- as.matrix(X %*% A) + 100
  Y <- matrix(rpois(length(mu), mu), n, 3)
  start <- convkernel("gaussian", L, list(sigma = 6))
  fit <- L0glm_blind_deco(Y, kernel_est = start,
    baseline = "intercept", nwindows = 10L, maxiter = 6L)
  expect_equal(dim(fit$g_hat), c(10L, L))
  expect_equal(rowSums(fit$g_hat), rep(1, 10), ignore_attr = TRUE, tolerance = 1e-12)
  expect_equal(dim(fit$a_hat), dim(A))
  expect_identical(fit$X, fit$X_spikes)
  expect_equal(fit$fitted_values, fit$X %*% fit$a_hat + fit$fitted_baseline,
               ignore_attr = TRUE, tolerance = 1e-10)
  expect_equal(fit$deviance, L0glm:::L0glm_deco_deviance(Y, fit$fitted_values, poisson(identity)))
  truth <- t(vapply(fit$psf_timepoints, function(t) {
    h <- exp(L0glm:::.psf_interp_log_kernel(t, centers, log(kernels), smooth = TRUE))
    h / sum(h)
  }, numeric(L)))
  cosine <- rowSums(fit$g_hat * truth) / sqrt(rowSums(fit$g_hat^2) * rowSums(truth^2))
  expect_true(all(cosine > .97))
  expect_true(which.max(fit$g_hat[1L, ]) %in% 15:17)
  expect_lt(sum(fit$g_hat[10L, ]^2), sum(fit$g_hat[1L, ]^2))
  shared <- L0glm_blind_deco(Y, kernel_est = start, baseline = "intercept", maxiter = 6L)
  f1 <- function(a) {selected <- Matrix::rowSums(a) > 0; 2 * sum(selected[at]) / (sum(selected) + length(at))}
  expect_gt(f1(fit$a_hat), f1(shared$a_hat))
  expect_lt(fit$deviance, shared$deviance)
  # The design is ready for the log-concave compound-fusion workflow.
  expect_error(spoloconmf(fit$X, Y, background = matrix(1, n, 1), ncores = 1L,
    control = list(broad_repeats = 0L, full_sweeps = 0L, final_polishes = 0L)), NA)
})

test_that("blind deconvolution retains its best completed IC checkpoint", {
  set.seed(12)
  g <- convkernel("gaussian", 21L, list(sigma = 2))
  a <- numeric(200); a[c(30, 80, 120, 170)] <- 10000
  d <- list(kernel = g, Y = rpois(200, as.numeric(convmatrix(g, 200) %*% a) + 500))
  fit <- L0glm_blind_deco(d$Y, kernel_est = d$kernel,
    baseline = "intercept", maxiter = 4L)
  expect_equal(fit$selection_score, min(fit$selection_trace))
  expect_equal(fit$best_iteration, which.min(fit$selection_trace) - 1L)
  check <- L0glm_sparse_deco(d$Y, kernel = fit$g_hat, baseline = "intercept")
  expect_equal(fit$selection_score, GIC(check$fit, mu_floor = 1e-8))
  expect_equal(fit$a_hat, check$coef_spikes)
})

test_that("local starts and window counts are validated", {
  Y <- matrix(0, 150, 2)
  g <- convkernel("gaussian", 21L)
  expect_error(L0glm_blind_deco(Y, kernel_est = g, nwindows = 0L), "nwindows")
  expect_error(L0glm_blind_deco(Y, kernel_est = g, nwindows = 20L), "at least L")
  expect_error(L0glm_blind_deco(Y, kernel_est = rbind(g, g), nwindows = 3L), "one row per window")
  fit <- suppressWarnings(L0glm_blind_deco(Y, kernel_est = rbind(g, g, g),
    L = 5L, nwindows = 3L, maxiter = 1L, kernel_smoothfun = NULL))
  expect_equal(fit$g_hat, matrix(rep(g, each = 3L), 3), ignore_attr = TRUE)
  expect_true(all(fit$a_hat == 0))
  expect_equal(dim(fit$X), c(150L, 150L))
})

test_that("time-varying Gaussian fits retain their individual channels", {
  set.seed(14)
  g <- convkernel("gaussian", 21L, list(sigma = 2))
  A <- matrix(0, 240, 2); A[c(35, 95, 155, 210), ] <- 20000
  mu <- as.matrix(convmatrix(g, 240) %*% A) + 200
  Y <- mu + matrix(rnorm(length(mu), sd = 5), 240, 2)
  fit <- L0glm_blind_deco(Y, family = gaussian(identity), kernel_est = g,
    baseline = "intercept", nwindows = 2L, maxiter = 2L, parallelize = FALSE)
  expect_equal(dim(fit$a_hat), c(240L, 2L))
  expect_true(all(is.finite(fit$g_hat)))
  expect_true(is.finite(fit$deviance))
  expect_equal(fit$fitted_values, fit$X %*% fit$a_hat + fit$fitted_baseline,
               ignore_attr = TRUE, tolerance = 1e-9)
})
