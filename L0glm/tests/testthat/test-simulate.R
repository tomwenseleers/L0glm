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
