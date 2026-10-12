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
  expect_equal(metrics$F1_score, 0)
  expect_equal(metrics$MCC, -1)
})

test_that("MCC is binary Pearson correlation and handles undefined margins", {
  truth <- c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)
  selected <- c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE)
  expect_equal(L0glm:::support_mcc(2, 1, 3, 1), cor(as.numeric(truth), as.numeric(selected)))
  expect_equal(L0glm:::support_mcc(100000L, 50000L, 200000L, 30000L),
               (1e5 * 2e5 - 5e4 * 3e4) / sqrt(1.5e5 * 1.3e5 * 2.5e5 * 2.3e5))
  expect_equal(L0glm:::support_mcc(2, 0, 3, 0), 1)
  expect_true(is.na(L0glm:::support_mcc(0, 0, 5, 2)))
  expect_true(is.na(L0glm:::support_mcc(0, 3, 4, 0)))
})
