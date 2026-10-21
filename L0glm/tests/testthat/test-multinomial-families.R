test_that("coupled multinomial IRLS handles count matrices and agrees with multinom", {
  skip_if_not_installed("mgcv"); skip_if_not_installed("nnet")
  set.seed(1806)
  X <- matrix(rnorm(400 * 3), 400, 3)
  B <- matrix(c(.2, .6, -.4, .1, -.3, -.2, .3, .6), 4, 2)
  P <- L0glm_ext_softmax(cbind(0, cbind(1, X) %*% B))$p
  Y <- t(vapply(seq_len(400), function(i) as.numeric(rmultinom(1, 3, P[i, ])), numeric(3)))
  colnames(Y) <- c("reference", "A", "B")
  w <- sample(1:2, 400, TRUE)
  ref <- nnet::multinom(Y ~ X, weights = w, trace = FALSE, maxit = 300,
                        abstol = 1e-12, reltol = 1e-12)
  for (mode in c("irls", "bound")) {
    fit <- L0glm.fit(X, Y, family = mgcv::multinom(K = 2), lambda = 0,
      weights = w, relaxed = FALSE, maxit = 400, tol = 1e-9,
      control.family = list(multinomial = mode))
    ref_clr <- cbind(0, t(coef(ref))); ref_clr <- ref_clr - rowMeans(ref_clr)
    expect_equal(as.matrix(fit$coefficients), unname(ref_clr),
                 ignore_attr = TRUE, tolerance = 1e-5)
    expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)) +
      sum(w * (lgamma(rowSums(Y) + 1) - rowSums(lgamma(Y + 1)))), tolerance = 1e-6)
    expect_equal(rowSums(predict(fit)), rep(1, 400), tolerance = 1e-13)
    expect_equal(as.numeric(rowSums(fit$coefficients)), rep(0, 4), tolerance = 1e-13)
    expect_true(all(fit$convergence$converged))
  }
})

test_that("multinomial grouped penalties do not cancel opposite signs", {
  B <- rbind(c(1, -1), c(.5, .5), c(0, 0))
  D <- L0glm_ext_penalty(B, matrix(TRUE, 3, 2), integer(), 8, TRUE)
  expect_equal(D[1, ], c(1, 1))
  expect_equal(D[2, ], c(4, 4))
  expect_true(all(is.infinite(D[3, ])))
  threshold <- L0glm_ext_threshold(B, matrix(TRUE, 3, 2), integer(), .01, TRUE)
  expect_equal(threshold$B[1, ], c(1, -1))
  expect_true(all(!threshold$allowed[3, ]))
})

test_that("multinomial CLR normalization respects row masks and free category starts", {
  skip_if_not_installed("mgcv")
  set.seed(1808)
  X <- matrix(rnorm(300 * 4), 300, 4)
  P <- L0glm_ext_softmax(cbind(0, .5 + .6*X[, 1] - .5*X[, 2], -.3 - .5*X[, 1] + .7*X[, 2]))$p
  y <- factor(vapply(seq_len(300), function(i) sample.int(3, 1, prob = P[i, ]), integer(1)))
  a <- L0glm.fit(X, y, family = mgcv::multinom(K = 2), lambda = 0,
    relaxed = FALSE, normalize = "rms", maxit = 100, tol = 1e-9)
  scale <- c(.01, .5, 2, 50)
  b <- L0glm.fit(sweep(X, 2, scale, "*"), y, family = mgcv::multinom(K = 2),
    lambda = 0, relaxed = FALSE, normalize = "rms", maxit = 100, tol = 1e-9)
  expect_equal(predict(a), predict(b), tolerance = 1e-9)
  expect_equal(as.matrix(a$coefficients), as.matrix(b$coefficients) * c(1, scale), tolerance = 1e-8)
  initial <- matrix(1, 5, 3); initial[1, ] <- initial[3, ] <- initial[5, ] <- 0
  initial[2, 1] <- 0
  fit <- L0glm.fit(X, y, family = mgcv::multinom(K = 2), start = initial,
    always.include = 3L, lambda = 0, grouped = TRUE, relaxed = FALSE,
    maxit = 150, tol = 1e-8)
  expect_equal(dim(fit$coefficients), c(5L, 3L))
  expect_true(all(fit$coefficients[5, ] == 0))
  expect_true(any(fit$coefficients[2, ] != 0))
  expect_true(any(fit$coefficients[3, ] != 0))
  expect_true(any(fit$coefficients[1, ] != 0))
  expect_s4_class(fit$coefficients, "sparseMatrix")
  data <- L0glm_ext_prepare(X, y, mgcv::multinom(K = 2), NULL, NULL, TRUE,
    NULL, FALSE, "l2", 1, NULL, list())
  expect_equal(as.numeric(data$B), rep(0, 15), tolerance = 1e-12)
  expect_true(all(data$allowed))
  unit <- L0glm.fit(X, y, family = mgcv::multinom(K = 2), start = 1,
    lambda = 0, relaxed = FALSE, normalize = "rms", maxit = 100, tol = 1e-9)
  expect_equal(unit$coefficients, a$coefficients, tolerance = 1e-9)
  expect_warning(signed <- L0glm.fit(X, y, family = mgcv::multinom(K = 2),
    nonnegative = TRUE, lambda = 0, relaxed = FALSE, normalize = "rms", maxit = 100, tol = 1e-9), "ignored")
  expect_equal(signed$coefficients, a$coefficients, tolerance = 1e-9)
})

test_that("CLR quadratic solve obeys equality KKT conditions through OSQP", {
  skip_if_not_installed("mgcv"); skip_if_not_installed("osqp")
  set.seed(1809)
  X <- matrix(rnorm(120 * 4), 120, 4)
  y <- factor(sample(1:3, 120, TRUE, prob = c(.6, .3, .1)))
  data <- L0glm_ext_prepare(X, y, mgcv::multinom(K = 2), NULL, NULL,
    TRUE, NULL, FALSE, "rms", NULL, NULL, list())
  B <- matrix(0, data$p, data$m); D <- matrix(2, data$p, data$m); D[1, ] <- 0
  args <- list(solver = "osqp", method = 1L,
    control.solver = control.solver.gen("osqp", abstol = 1e-9, reltol = 1e-9))
  candidate <- L0glm_ext_multinomial_working(data, B, D, list(), args)
  H <- kronecker(diag(1/3, 3) - matrix(1/9, 3, 3), crossprod(data$Z))
  score <- as.numeric(crossprod(data$Z, data$y - 1/3))
  gradient <- as.numeric(H %*% as.numeric(candidate)) - score + as.numeric(D * candidate)
  gradient <- matrix(gradient, data$p, data$m)
  expect_lt(max(abs(gradient-rowMeans(gradient))), 1e-5)
  expect_lt(max(abs(rowSums(candidate))), 1e-12)
  expect_true(any(candidate < 0))
})
