deco_fixture <- function(m = 3L, n = 160L, baseline = TRUE) {
  set.seed(470)
  kernel <- convkernel("gaussian", 21L, list(sigma = 2), perc = 0)
  X <- convmatrix(kernel, n)
  A <- matrix(0, n, m)
  A[c(25, 65, 105, 140), ] <- matrix(seq(4000, 16000, length.out = 4L * m), 4L, m)
  B <- splines::bs(seq_len(n), df = 6L, intercept = TRUE)
  alpha <- matrix(rep(seq(300, 1200, length.out = m), each = 6L), 6L, m)
  background <- if (baseline) B %*% alpha else matrix(0, n, m)
  mu <- as.matrix(X %*% A) + background
  Y <- matrix(rpois(n * m, mu), n, m)
  list(Y = Y, A = A, kernel = kernel, X = X, background = background)
}

test_that("fixed deconvolution recovers peaks with fitted backgrounds", {
  d <- deco_fixture()
  for (baseline in c("intercept", "bspline")) {
    fit <- L0glm_sparse_deco(d$Y, kernel = d$kernel, baseline = baseline, nknots = 6L)
    expect_s4_class(fit$coef_spikes, "sparseMatrix")
    expect_equal(dim(fit$coef_spikes), dim(d$A))
    selected <- as.numeric(Matrix::rowSums(fit$coef_spikes)) > 0
    truth <- rowSums(d$A) > 0
    expect_true(all(selected[truth]))
    expect_lte(sum(selected & !truth), 2L)
    expect_true(all(fit$coef_spikes >= 0))
    expect_equal(fit$fitted_values,
      fit$X_spikes %*% fit$coef_spikes + fit$fitted_baseline + fit$offset,
      tolerance = 1e-10, ignore_attr = TRUE)
    expect_equal(fit$deviance, L0glm:::L0glm_deco_deviance(d$Y, fit$fitted_values, poisson(identity)))
    expect_lt(sum((as.matrix(fit$coef_spikes) - d$A)^2) / sum(d$A^2), 0.02)
  }
})

test_that("single-channel and custom rectangular designs preserve dimensions", {
  d <- deco_fixture(m = 1L, baseline = FALSE)
  fit <- L0glm_sparse_deco(as.numeric(d$Y), kernel = d$kernel)
  expect_equal(dim(fit$coef_spikes), c(160L, 1L))
  expect_null(fit$coef_baseline)
  expect_equal(as.matrix(fit$fitted_baseline), matrix(0, 160L, 1L))
  X <- as.matrix(d$X[, c(25, 65, 105, 140)])
  rectangular <- L0glm_sparse_deco(d$Y, X = X,
    kernelfunc = function(x) stop("X overrides kernel generation"))
  expect_equal(dim(rectangular$coef_spikes), c(4L, 1L))
  expect_identical(rectangular$X_spikes, X)
  expect_equal(rectangular$fitted_values, X %*% rectangular$coef_spikes,
    tolerance = 1e-9, ignore_attr = TRUE)
})

test_that("known offsets fit original counts and broadcast correctly", {
  d <- deco_fixture()
  fit <- L0glm_sparse_deco(d$Y, kernel = d$kernel, offset = d$background)
  expect_equal(fit$offset, d$background)
  expect_equal(which(Matrix::rowSums(fit$coef_spikes) > 0), which(rowSums(d$A) > 0))
  expect_equal(fit$fitted_values, fit$X_spikes %*% fit$coef_spikes + d$background,
    tolerance = 1e-10, ignore_attr = TRUE)
  one <- L0glm_sparse_deco(d$Y[, 1L], kernel = d$kernel, offset = 300)
  vector <- L0glm_sparse_deco(d$Y[, 1L], kernel = d$kernel, offset = rep(300, 160L))
  expect_equal(one$coef_spikes, vector$coef_spikes)
})

test_that("PSF updates use the correct channel-specific offsets", {
  d <- deco_fixture()
  joint <- L0glm_deco_estkernel(d$Y, d$A, offset = d$background,
    kernel_est = d$kernel, kernel_smoothfun = NULL)
  separate <- do.call(rbind, lapply(1:3, function(j)
    L0glm_deco_estkernel(d$Y[, j], d$A[, j], offset = d$background[, j],
      kernel_est = d$kernel, kernel_smoothfun = NULL)$kernel))
  expect_equal(joint$kernels, separate, tolerance = 1e-10)
  expect_equal(joint$kernel, colMeans(separate), tolerance = 1e-10)
  expect_equal(sum(joint$kernel), 1, tolerance = 1e-12)
  expect_true(all(is.finite(joint$kernel) & joint$kernel >= 0))
  expect_gt(sum(joint$kernel * d$kernel) / sqrt(sum(joint$kernel^2) * sum(d$kernel^2)), 0.95)
})

test_that("uninformative channels cannot introduce NaN PSFs", {
  d <- deco_fixture()
  d$Y[, 2L] <- 0; d$A[, 2L] <- 0; d$background[, 2L] <- 0
  fit <- L0glm_deco_estkernel(d$Y, d$A, offset = d$background,
    kernel_est = d$kernel, kernel_smoothfun = NULL)
  expect_identical(as.logical(fit$informative), c(TRUE, FALSE, TRUE))
  expect_equal(fit$kernel, colMeans(fit$kernels[c(1, 3), ]))
  expect_true(all(is.finite(fit$kernels)))
  expect_warning(empty <- L0glm_deco_estkernel(matrix(0, 20, 2), matrix(0, 20, 2),
    kernel_est = d$kernel), "No nonzero spikes")
  expect_equal(empty$kernel, d$kernel)
  expect_equal(dim(empty$kernels), c(2L, 21L))
})

test_that("PSF updates preserve the width of broad positive kernels", {
  set.seed(104)
  n <- 500L; x <- -25:25
  spikes <- matrix(0, n, 1L)
  spikes[c(60, 150, 240, 330, 420), 1L] <- c(5, 8, 6, 3, 7) * 1e5
  width <- function(h) diff(range(x[h >= max(h) / 2]))
  for (sigma in c(8, 12)) {
    kernel <- dnorm(x, sd = sigma); kernel <- kernel / sum(kernel)
    mean <- as.numeric(convmatrix(kernel, n) %*% spikes) + 5
    y <- rpois(n, mean)
    fit <- L0glm_deco_estkernel(y, spikes, kernel_est = kernel,
      offset = 5, kernel_smoothfun = NULL)
    expect_equal(width(fit$kernel), width(kernel), tolerance = 2)
    cosine <- sum(fit$kernel * kernel) / sqrt(sum(fit$kernel^2) * sum(kernel^2))
    expect_gt(cosine, .999)
    expect_equal(sum(fit$kernel), 1, tolerance = 1e-12)
  }
})

test_that("blind fits return synchronized PSF, spikes, baseline and deviance", {
  d <- deco_fixture()
  for (m in c(1L, 3L)) {
    fit <- L0glm_blind_deco(d$Y[, seq_len(m), drop = FALSE], kernel_est = d$kernel,
      baseline = "bspline", nknots = 6L, kernel_smoothfun = NULL, maxiter = 2L)
    expect_equal(dim(fit$a_hat), c(160L, m))
    expect_equal(dim(fit$coef_baseline), c(6L, m))
    expect_equal(sum(fit$g_hat), 1, tolerance = 1e-12)
    expect_equal(fit$X_spikes, convmatrix(fit$g_hat, 160L), tolerance = 0)
    expect_equal(fit$fitted_values, fit$X_spikes %*% fit$a_hat + fit$fitted_baseline,
      tolerance = 1e-10, ignore_attr = TRUE)
    expect_equal(fit$fitted_baseline, fit$X_baseline %*% fit$coef_baseline,
      tolerance = 1e-10, ignore_attr = TRUE)
    expect_equal(tail(fit$deviance_trace, 1L), fit$deviance)
    expect_equal(fit$deviance, L0glm:::L0glm_deco_deviance(d$Y[, seq_len(m), drop = FALSE],
      fit$fitted_values, poisson(identity)))
    expect_true(all(is.finite(fit$deviance_trace)))
    expect_gt(sum(Matrix::rowSums(fit$a_hat) > 0), 0)
  }
})

test_that("Gaussian identity fits and default PSF smoothing work", {
  d <- deco_fixture(m = 1L)
  y <- as.numeric(d$X %*% d$A) + as.numeric(d$background)
  fixed <- L0glm_sparse_deco(y, family = gaussian(identity), kernel = d$kernel,
    baseline = "intercept")
  expect_gt(sum(fixed$coef_spikes > 0), 0)
  expect_true(is.finite(fixed$deviance))
  expect_equal(fixed$deviance, sum((y - as.numeric(fixed$fitted_values))^2))
  blind <- L0glm_blind_deco(d$Y, kernel_est = d$kernel, baseline = "intercept", maxiter = 1L)
  expect_true(all(is.finite(blind$g_hat)))
  expect_equal(sum(blind$g_hat), 1, tolerance = 1e-12)
})

test_that("blind deconvolution learns a centered asymmetric PSF", {
  set.seed(123)
  n <- 500L; L <- 51L
  kernel <- convkernel("bigaussian", L, list(sigma1 = 1.25, sigma2 = 5), perc = "mean")
  A <- matrix(0, n, 3L)
  at <- sort(sample(seq.int(25, n - 25), 50L))
  A[at, ] <- matrix(sample(1000:100000, 150L), 50L, 3L)
  Y <- matrix(rpois(n * 3L, as.matrix(convmatrix(kernel, n) %*% A) + 3000), n, 3L)
  fit <- L0glm_blind_deco(Y[, 1L],
    kernel_est = convkernel("gaussian", L, list(sigma = 6)),
    baseline = "intercept", maxiter = 8L)
  cosine <- sum(fit$g_hat * kernel) / sqrt(sum(fit$g_hat^2) * sum(kernel^2))
  estimate <- as.numeric(fit$a_hat); truth <- A[, 1L]
  expect_gt(cosine, 0.99)
  expect_gt(1 - sum((estimate - truth)^2) / sum((truth - mean(truth))^2), 0.95)
})

test_that("invalid deconvolution inputs fail before fitting", {
  d <- deco_fixture(m = 1L)
  expect_error(L0glm_sparse_deco(d$Y, family = poisson(log)), "identity link")
  expect_error(L0glm_sparse_deco(d$Y, kernel = rep(0, 5)), "positive sum")
  expect_error(L0glm_sparse_deco(d$Y, baseline = "unknown"), "arg")
  expect_error(L0glm_sparse_deco(d$Y, X = matrix(1, 3, 2)), "nrow")
  expect_error(L0glm_deco_estkernel(d$Y, d$A[-1, , drop = FALSE]), "same dimensions")
  expect_error(L0glm_deco_estkernel(d$Y, d$A, kernel_est = d$kernel,
    kernel_smoothfun = function(y) y[-1]), "length-L")
  expect_error(L0glm_blind_deco(d$Y, maxiter = 0), "maxiter")
})
