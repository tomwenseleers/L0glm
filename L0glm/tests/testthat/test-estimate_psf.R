psf_test_data <- function(kind) {
  set.seed(419)
  n <- 260L;m <- 4L;g <- dnorm(-15:15, sd = 2.3);g <- g / sum(g)
  A <- matrix(0, n, m)
  A[c(35, 95, 155, 220), ] <- matrix(runif(16, 1e4, 3e4), 4L, m)
  mu <- as.matrix(convmatrix(g, n) %*% A) + 20
  family <- switch(kind, gaussian = gaussian("identity"), poisson = poisson("identity"),
    quasipoisson = quasipoisson("identity"), nb = mgcv::nb(theta = 20, link = "identity"),
    nb1 = nb1(4, "identity"), Gamma = Gamma_censored(.1, "identity"),
    Gamma_scale = Gamma_scale(4, "identity", censor = .5))
  y <- switch(kind, gaussian = mu + rnorm(length(mu), sd = 4),
    poisson = rpois(length(mu), mu), quasipoisson = 4 * rpois(length(mu), mu / 4),
    nb = rnbinom(length(mu), mu = mu, size = 20),
    nb1 = rnbinom(length(mu), mu = mu, size = mu / 3),
    Gamma = rgamma(length(mu), shape = 10, scale = mu / 10),
    Gamma_scale = rgamma(length(mu), shape = mu / 4, scale = 4))
  Y <- matrix(y, n, m);rownames(Y) <- seq_len(n) / 10
  list(Y = Y, g = g, family = family,
    dispersion = switch(kind, gaussian = 16, quasipoisson = 4, NULL))
}

test_that("standard and blind PSFs fit all signal families with weights and offsets", {
  skip_if_not_installed("mgcv")
  for (kind in c("gaussian", "poisson", "quasipoisson", "nb", "nb1", "Gamma", "Gamma_scale")) {
    d <- psf_test_data(kind)
    for (method in c("standard", "blind")) {
      z <- estimate_psf(d$Y, family = d$family, standard = d$Y, method = method,
        nwindows = 2L, precluster = 2L, dispersion = d$dispersion, ncores = 2L,
        weights = seq(.8, 1.2, length.out = nrow(d$Y)), offset = 3,
        control = list(L = 31L, n_peaks = 4L, kernel_est = d$g,
          maxiter = 2L, baseline_df = 4L))
      expect_s4_class(z$X, "dgCMatrix")
      expect_identical(dim(z$X), c(260L, 260L))
      expect_equal(as.numeric(Matrix::colSums(z$X)), rep(1, 260L), tolerance = 1e-10)
      expect_equal(rowSums(z$psf_kernels), c(1, 1), ignore_attr = TRUE, tolerance = 1e-10)
      cosine <- as.numeric(z$psf_kernels %*% d$g) / sqrt(rowSums(z$psf_kernels^2) * sum(d$g^2))
      expect_gt(min(cosine), .97)
      for (j in c(1L, 40L, 130L, 240L, 260L)) {
        h <- as.numeric(z$X[, j]);h <- h[h > 0]
        expect_lte(max(diff(log(h), differences = 2L)), 1e-6)
      }
      expect_identical(z$aggregation$moment_matched, kind %in% c("nb", "Gamma"))
      expect_false(z$aggregation$weight_approximation)
    }
  }
})

test_that("PSF input formats, attributes and fixed designs retain their meaning", {
  d <- psf_test_data("poisson");Y <- d$Y
  at <- seq_len(nrow(Y)) / 10
  standard <- setNames(rowSums(Y), at)
  attr(Y, "alkanes_TIC") <- standard
  a <- estimate_psf(Y, method = "standard", ncores = 1L, control = list(L = 31L))
  b <- estimate_psf(Y, standard = list(Y = matrix(standard, ncol = 1L), time = at),
    method = "standard", ncores = 2L, control = list(L = 31L))
  expect_equal(a$X, b$X, tolerance = 1e-12)
  expect_equal(a$time, at)
  fixed <- estimate_psf(Y, precluster = FALSE, ncores = 1L,
    control = list(kernel_est = d$g, maxiter = 0L))
  expect_identical(fixed$method, "fixed")
  expect_equal(fixed$g_hat, d$g, tolerance = 1e-12)
  expect_equal(fixed$iterations, 0L)
  expect_equal(fixed$aggregation$original_channels, 4L)
  expect_null(fixed$calibration)
  expect_null(fixed$Y)
  expected <- slcn_project_design_cpp(as(convmatrix(d$g, nrow(Y)), "dgCMatrix"), 1L)
  expect_identical(fixed$X, expected)
  expect_error(estimate_psf(Y, control = list(typo = 1)), "Unknown PSF control")
  expect_error(estimate_psf(Y, control = list(maxiter = 0L), standard = NULL), "Supply kernel_est")
  expect_error(estimate_psf(Y, time = rev(at)), "increasing")
  expect_error(estimate_psf(Y, nwindows = 0L), "nwindows")
  expect_error(estimate_psf(Y, control = list(L = 30L)), "odd")
})

test_that("rectangular convolution matches direct indexing with sparse kernels", {
  for (n in c(7L, 19L, 31L)) for (p in c(3L, 23L, 39L)) {
    set.seed(n + p);g <- runif(19L);g[g < .5] <- 0;origin <- 10L
    z <- convmatrix(g, n = n, n_cols = p, origin = origin)
    index <- outer(seq_len(n), seq_len(p), "-") + origin
    expected <- matrix(0, n, p);keep <- index >= 1L & index <= length(g)
    expected[keep] <- g[index[keep]]
    expect_equal(as.matrix(z), expected, tolerance = 1e-12)
  }
})

test_that("censored aggregation is identified independently of variance closure", {
  Y <- matrix(rep(c(0, 1, 2, 4), 10L), 20L, 2L)
  args <- list(y=Y, family=Gamma_scale(4, "identity", censor=.5), ncores=1L,
    control=list(kernel_est=c(.1,.8,.1), maxiter=0L))
  grouped <- do.call(estimate_psf, args)
  original <- do.call(estimate_psf, c(args, list(precluster=FALSE)))
  expect_true(grouped$aggregation$censored_approximation)
  expect_false(grouped$aggregation$moment_matched)
  expect_false(original$aggregation$censored_approximation)
})
