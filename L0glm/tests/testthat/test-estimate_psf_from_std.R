test_that("estimate_psf_from_std returns a normalized sparse PSF matrix", {
  set.seed(1)
  time <- seq(0, 20, length.out = 500)
  centers <- seq(2, 18, length.out = 8)
  baseline <- 25 + 4 * sin(time / max(time) * 2 * pi)
  y <- baseline

  for (k in seq_along(centers)) {
    sigma <- 0.035 + 0.004 * k
    amp <- 500 + 25 * k
    y <- y + amp * exp(-0.5 * ((time - centers[k]) / sigma)^2)
  }

  y <- pmax(y + rnorm(length(y), sd = 2), 0)
  names(y) <- as.character(time)

  X <- estimate_psf_from_std(
    y = y,
    family = stats::gaussian(link = "identity"),
    baseline_df = 6,
    output_timepoints = time,
    n_peaks = length(centers),
    peak_min_distance = 1,
    peak_height = 100,
    psf_length = 31
  )

  expect_s4_class(X, "dgCMatrix")
  expect_equal(dim(X), c(length(time), length(time)))
  expect_equal(length(attr(X, "peak_timepoints")), length(centers))
  expect_equal(ncol(attr(X, "psf_kernels")), 31)
  expect_equal(as.numeric(Matrix::colSums(X)), rep(1, ncol(X)), tolerance = 1e-8)
  expect_equal(rowSums(attr(X, "psf_kernels")), rep(1, length(centers)), tolerance = 1e-8)
})

test_that("estimate_psf_from_std accepts character family aliases", {
  time <- seq(0, 10, length.out = 250)
  centers <- c(2, 4, 6, 8)
  y <- rep(5, length(time))
  for (center in centers) {
    y <- y + 100 * exp(-0.5 * ((time - center) / 0.05)^2)
  }

  X <- estimate_psf_from_std(
    y = y,
    family = "poisson",
    baseline_df = 5,
    output_timepoints = time,
    n_peaks = length(centers),
    peak_min_distance = 1,
    peak_height = 20,
    psf_length = 21
  )

  expect_s4_class(X, "dgCMatrix")
  expect_equal(attr(X, "family")$family, "poisson")
  expect_equal(as.numeric(Matrix::colSums(X)), rep(1, ncol(X)), tolerance = 1e-8)
})

test_that("estimate_psf_from_std infers retention times from numeric names", {
  time <- seq(5, 7, length.out = 120)
  y <- 2 + 50 * exp(-0.5 * ((time - 6) / 0.04)^2)
  names(y) <- as.character(time)

  fit <- estimate_psf_from_std(
    y = y,
    family = "gaussian",
    baseline_df = 5,
    output_timepoints = time,
    n_peaks = 1,
    peak_height = 10,
    psf_length = 15,
    return_details = TRUE
  )

  expect_equal(fit$time, time, tolerance = 1e-12)
  expect_true(abs(fit$peak_timepoints - 6) < 0.05)
})

test_that("standard peak windows retain both sides of neighbouring peaks", {
  time <- seq(0, 10, by = .01)
  sigma <- c(.15, .2)
  y <- 1e6 * (dnorm(time, 3, sigma[1L]) + dnorm(time, 6, sigma[2L]))
  at <- vapply(c(3, 6), function(t) which.min(abs(time - t)), 1L)
  windows <- L0glm:::.psf_peak_windows(time, y, at, 0, .005, .999)
  for (j in seq_along(windows)) {
    w <- windows[[j]]
    expect_lt(min(w$x), time[at[j]])
    expect_gt(max(w$x), time[at[j]])
    fwhm <- diff(range(w$x[w$signal >= max(w$signal) / 2]))
    expect_lt(abs(fwhm - 2 * sqrt(2 * log(2)) * sigma[j]), 2 * median(diff(time)))
  }
})

test_that("automatic standard support recovers broad as well as narrow PSFs", {
  time <- 0:1000
  centers <- c(120, 350, 600, 850)
  sigma <- c(4, 7, 10, 12)
  y <- rep(100, length(time))
  for (j in seq_along(centers)) y <- y + 1e6 * dnorm(time, centers[j], sigma[j])
  fit <- estimate_psf_from_std(y, time = time, n_peaks = 4L,
    peak_height = 1000, peak_min_distance = 100, return_details = TRUE)
  width <- function(x, h) diff(range(x[h >= max(h) / 2]))
  expected <- vapply(sigma, function(s) width(time, dnorm(time, 500, s)), 0)
  estimated <- apply(fit$psf_kernels, 1L, function(h) width(fit$psf_offsets, h))
  expect_gt(fit$psf_length, 51L)
  expect_equal(fit$psf_length %% 2L, 1L)
  expect_lt(max(abs(estimated / expected - 1)), .1)
  expect_s4_class(fit$X, "dgCMatrix")
  expect_equal(as.numeric(Matrix::colSums(fit$X)), rep(1, length(time)), tolerance = 1e-10)
})
