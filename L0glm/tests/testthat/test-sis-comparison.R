test_that("Gaussian score screening agrees with the initial SIS-package rank", {
  skip_if_not_installed("SIS")
  set.seed(431)
  X <- matrix(rnorm(150 * 400), 150, 400)
  colnames(X) <- paste0("x", seq_len(ncol(X)))
  y <- 2 * X[, 1] - X[, 3] + rnorm(nrow(X))
  actual <- sis_glm(X, y, nsis = 30L, nthreads = 1L)
  for (standardize in c(FALSE, TRUE)) {
    reference <- SIS::SIS(X, y, family = "gaussian", nsis = 30L,
      iter = FALSE, penalty = "lasso", tune = "bic", standardize = standardize,
      parallel = FALSE, seed = 431L)
    expect_identical(actual$indices, as.integer(reference$sis.ix0))
    expect_identical(actual$ranking, as.integer(reference$ix_list[[1L]]))
  }
})
