test_that("multichannel simul methods allow more than 100 peaks", {
  n <- p <- 150L
  npeaks <- 120L

  sim_gaus <- simul(method = 3, n = n, p = p, s = npeaks / p, seed = 2)
  sim_pois <- simul(method = 4, n = n, p = p, s = npeaks / p, seed = 2)

  support_size <- function(sim) {
    sum(Matrix::rowSums(sim$beta_true != 0) > 0)
  }

  expect_equal(support_size(sim_gaus), npeaks)
  expect_equal(support_size(sim_pois), npeaks)
})

test_that("simul caps requested peaks at possible spike locations", {
  n <- p <- 30L
  requested_peaks <- 60L

  sim <- simul(method = 4, n = n, p = p, s = requested_peaks / p, seed = 2)
  support_size <- sum(Matrix::rowSums(sim$beta_true != 0) > 0)

  expect_equal(support_size, n)
})

test_that("named deconvolution supports all noises and explicit channel counts", {
  for (family in c("poisson", "gaussian", "nb")) for (m in c(1L, 5L)) {
    if (family == "nb") skip_if_not_installed("mgcv")
    sim <- simul("deconvolution", n=70, p=90, n_nonzero=8, family=family,
      m=m, intercept=2, seed=81)
    expect_s4_class(sim$X, "sparseMatrix")
    expect_equal(dim(sim$y), c(70L, m))
    expect_equal(dim(sim$beta_true), c(91L, m))
    expect_equal(sum(Matrix::rowSums(abs(sim$beta_true[-1L, , drop=FALSE])) > 0), 8)
    expect_equal(unname(sim$y_true), unname(as.matrix(cbind(1, sim$X) %*% sim$beta_true)))
    expect_true(all(is.finite(sim$y)))
    if (family != "gaussian") expect_true(all(sim$y >= 0 & sim$y == round(sim$y)))
  }
})

test_that("deconvolution preserves explicit spectra and null/one-peak dimensions", {
  templates <- rbind(c(1, 0, 0), c(0, 1, 0))
  sim <- simulate_spike_train(n=50, p=30, k=5, spectra=templates, intercept=c(1, 2, 3))
  expect_equal(dim(sim$y), c(50L, 3L))
  expect_true(all(sim$beta_true[-1L, 3L] == 0))
  expect_true(all(sim$y_true[, 3L] == 3))
  wrapped <- simul("deconvolution", n=50, p=30, n_nonzero=5, options=list(spectra=templates))
  expect_equal(ncol(wrapped$y), 3L)
  for (m in c(1L, 3L)) for (k in c(0L, 1L)) {
    sim <- simulate_spike_train(n=1, p=1, k=k, m=m, mean_beta=1, sd_logbeta=0)
    expect_equal(dim(sim$y), c(1L, m))
    expect_equal(sum(Matrix::rowSums(sim$beta_true != 0) > 0), k)
  }
  sim <- simulate_spike_train(n=5, p=7, k=0, multichannel=TRUE)
  expect_equal(nrow(sim$y), 5L)
  expect_gt(ncol(sim$y), 0L)
  expect_true(all(sim$y == 0))
})

test_that("Gaussian deconvolution is uncensored unless clipping is requested", {
  raw <- simulate_spike_train(n=100, k=0, family="gaussian", seed=8)
  clipped <- simulate_spike_train(n=100, k=0, family="gaussian", clip=TRUE, seed=8)
  expect_true(any(raw$y < 0))
  expect_equal(clipped$y, pmax(raw$y, 0))
  expect_true(all(raw$y_true == 0))
})

test_that("sparse and dense deconvolution have identical data", {
  for (family in c("poisson", "gaussian")) {
    a <- simulate_spike_train(n=60, p=80, k=5, m=3, family=family, sparse=TRUE, seed=2)
    b <- simulate_spike_train(n=60, p=80, k=5, m=3, family=family, sparse=FALSE, seed=2)
    expect_true(is.matrix(b$X))
    for (name in c("X", "beta_true", "y", "y_true"))
      expect_equal(as.matrix(a[[name]]), as.matrix(b[[name]]))
  }
  expect_error(simulate_spike_train(family=poisson()), "identity link")
  expect_error(simulate_spike_train(psf=function(x, u) -x), "nonnegative")
  expect_error(simulate_spike_train(spectra=matrix(0, 2, 3)), "row totals")
  a <- simulate_spike_train(n=40, p=40, k=4, multichannel=TRUE, family="gaussian", intercept=FALSE)
  b <- simulate_spike_train(n=40, p=40, k=4, multichannel=TRUE, family="gaussian", intercept=NULL)
  expect_identical(a$y, b$y)
  expect_equal(dim(a$weights), c(40L, 1L))
})
