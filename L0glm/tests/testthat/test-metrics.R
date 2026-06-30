test_that("get.metrics uses row support for multichannel coefficients", {
  beta <- Matrix::Matrix(c(1, 0, -1, 0), nrow = 2, sparse = TRUE)
  beta_true <- Matrix::Matrix(c(0, 2, 0, 0), nrow = 2, sparse = TRUE)
  metrics <- L0glm:::get.metrics(
    beta = beta,
    X = diag(2),
    y = matrix(0, 2, 2),
    beta_true = beta_true
  )

  expect_equal(metrics$support_size_estimated, 1)
  expect_equal(metrics$support_size_true, 1)
})
