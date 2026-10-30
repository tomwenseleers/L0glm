background_start_fixture <- function(sparse = TRUE) {
  set.seed(864)
  n <- 240L
  g <- convkernel("gaussian", 21L, list(sigma = 2), perc = 0)
  X <- convmatrix(g, n)
  A <- matrix(0, n, 3L)
  A[c(30, 85, 145, 205), ] <- matrix(seq(6000, 18000, length.out = 12L), 4L, 3L)
  B <- splines::bs(seq_len(n), df = 4L, intercept = TRUE)
  alpha <- matrix(rep(c(300, 900, 1800), each = 4L), 4L, 3L)
  mu <- as.matrix(X %*% A) + B %*% alpha
  Y <- matrix(rpois(length(mu), mu), n, 3L)
  Z <- cbind(X, Matrix::Matrix(B, sparse = TRUE))
  if (!sparse) Z <- as.matrix(Z)
  list(X = Z, y = Y, peaks = A, mandatory = n + 1:4, background = B %*% alpha)
}

test_that("automatic identity starts recover peaks above compulsory backgrounds", {
  for (sparse in c(FALSE, TRUE)) {
    d <- background_start_fixture(sparse)
    results <- lapply(c("none", "l2", "rms", "sd"), function(mode)
      L0glm.fit(d$X, d$y, family = poisson(identity), nonnegative = TRUE,
        always.include = d$mandatory, lambda = "gic", normalize = mode,
        parallelize = FALSE))
    for (a in results) {
      selected <- as.numeric(Matrix::rowSums(coef(a)[1:240, ])) > 0
      truth <- rowSums(d$peaks) > 0
      expect_true(all(selected[truth]))
      expect_lte(sum(selected & !truth), 2L)
      expect_equal(a$nuisance$initializer, "signal_backprojection")
      expect_identical(inherits(a$X, "sparseMatrix"), sparse)
      expect_equal(fitted(a), as.matrix(d$X %*% coef(a)), tolerance = 1e-9, ignore_attr = TRUE)
    }
    for (a in results[-1L]) expect_equal(as.numeric(coef(a)), as.numeric(coef(results[[1L]])), tolerance = 1e-6)
    native <- L0nnpois(d$X, d$y, always.include = d$mandatory, lambda = "gic", parallelize = FALSE)
    expect_equal(coef(native), coef(results[[1L]]), tolerance = 0)
  }
})

test_that("the signal pilot is equivariant to individual predictor scales", {
  d <- background_start_fixture()
  W <- matrix(rep(seq(.5, 2, length.out = 240L), 3L), 240L, 3L)
  original <- L0glm:::L0glm_signal_start(d$X, d$y, matrix(0, 240, 3), W, d$mandatory)
  factors <- c(seq(.5, 2, length.out = 240L), rep(1, 4))
  rescaled <- d$X %*% Matrix::Diagonal(x = factors)
  changed <- L0glm:::L0glm_signal_start(rescaled, d$y, matrix(0, 240, 3), W, d$mandatory)
  expect_equal(changed * factors, original, tolerance = 1e-10)
  expect_equal(as.matrix(rescaled %*% changed), as.matrix(d$X %*% original), tolerance = 1e-10)
  expect_identical(inherits(rescaled, "sparseMatrix"), TRUE)
})

test_that("known-offset automatic starts agree between native and general Poisson", {
  d <- background_start_fixture()
  X <- d$X[, 1:240, drop = FALSE]
  for (mode in c("none", "l2")) {
    a <- L0nnpois(X, d$y, offset = d$background, normalize = mode,
      lambda = "gic", parallelize = FALSE)
    b <- L0glm.fit(X, d$y, family = poisson(identity), nonnegative = TRUE,
      offset = d$background, normalize = mode, openmp = FALSE,
      lambda = "gic", parallelize = FALSE)
    selected <- as.numeric(Matrix::rowSums(coef(a))) > 0
    expect_true(all(selected[rowSums(d$peaks) > 0]))
    expect_equal(coef(a), coef(b), tolerance = 1e-6)
    expect_equal(fitted(a), as.matrix(X %*% coef(a)) + d$background, ignore_attr = TRUE)
  }
})

test_that("explicit starting zeros and zero-information channels retain their contracts", {
  d <- background_start_fixture()
  S <- matrix(1, ncol(d$X), 3L); S[30, ] <- 0
  a <- L0nnpois(d$X, d$y, start = S, always.include = d$mandatory, parallelize = FALSE)
  expect_true(all(coef(a)[30, ] == 0))
  expect_equal(a$nuisance$initializer, "joint")
  d$y[, 3L] <- 0
  b <- L0nnpois(d$X, d$y, always.include = d$mandatory, parallelize = FALSE)
  expect_true(all(coef(b)[, 3L] == 0))
  expect_equal(dim(coef(b)), c(244L, 3L))
})
