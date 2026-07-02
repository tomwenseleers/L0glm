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
