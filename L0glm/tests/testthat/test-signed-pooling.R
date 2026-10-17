test_that("signed starts pool magnitudes without channel cancellation", {
  B <- rbind(c(2, -2), c(2, -2 + 1e-12), c(-3, 1), c(0, 0))
  for (sparse in c(FALSE, TRUE)) {
    pooled <- L0glm_group_start(Matrix::Matrix(B, sparse = sparse))
    expect_equal(pooled$vars, 1:3)
    expect_equal(pooled$start[, 1], c(2, 2 - 0.5e-12, -2, 0), tolerance = 1e-13)
    expect_equal(pooled$start[, 1], pooled$start[, 2])
    expect_equal(L0glm_group_start(Matrix::Matrix(abs(B), sparse = sparse), TRUE)$start[, 1],
                 rowMeans(abs(B)))
  }
})

test_that("signed multichannel fits are invariant to response channel order", {
  set.seed(1640)
  X <- matrix(rnorm(80 * 45), 80, 45)
  B <- matrix(0, 45, 3)
  B[1:5, ] <- cbind(seq(1, 2, length.out = 5), seq(-2, -1, length.out = 5), rep(1, 5))
  Y <- X %*% B + matrix(rnorm(240, sd = 0.3), 80, 3)
  for (sparse in c(FALSE, TRUE)) {
    args <- list(X = Matrix::Matrix(X, sparse = sparse), family = gaussian(),
                 dispersion = 0.09, lambda = "gic", grouped = TRUE, relaxed = TRUE,
                 maxit = 50L, parallelize = FALSE)
    original <- do.call(L0glm.fit, c(args, list(y = Y)))
    order <- c(3, 1, 2)
    shuffled <- do.call(L0glm.fit, c(args, list(y = Y[, order])))
    expect_equal(as.matrix(original$coefficients)[, order], as.matrix(shuffled$coefficients), tolerance = 1e-9)
    expect_true(all(Matrix::rowSums(original$coefficients != 0)[1:5] > 0))
  }
})
