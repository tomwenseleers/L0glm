test_that("pilot row augmentation preserves the exact ridge system", {
  set.seed(1711)
  A <- matrix(0, 1100, 1000)
  A[cbind(seq_len(1000), seq_len(1000))] <- runif(1000, 1, 3)
  A[cbind(2:1001, seq_len(1000))] <- runif(1000)
  penalty <- c(0, rep(1e-6, 999))
  original <- A
  K <- L0glm_pilot_gram(A, penalty)
  expect_s4_class(K, "sparseMatrix")
  expect_identical(A, original)
  expect_equal(as.matrix(K), crossprod(A) + diag(penalty), tolerance = 1e-12)
  rhs <- matrix(rnorm(2000), 1000, 2)
  expect_equal(as.matrix(Matrix::solve(K, rhs)),
               solve(crossprod(A) + diag(penalty), rhs), tolerance = 1e-10)
  small <- A[1:20, 1:10]
  expect_identical(L0glm_pilot_gram(small, rep(1, 10)),
                   Matrix::crossprod(small) + Matrix::Diagonal(10))
  dense <- matrix(rnorm(1100 * 1000), 1100, 1000)
  expect_false(inherits(L0glm_pilot_gram(dense, rep(1, 1000)), "sparseMatrix"))
})
