test_that("NB coefficients match Fisher scoring at fixed or Pearson-estimated shape", {
  skip_if_not_installed("mgcv"); skip_if_not_installed("MASS")
  set.seed(1803)
  X <- matrix(rnorm(500 * 3), 500, 3)
  offset <- runif(500, -.2, .2); w <- sample(0:3, 500, TRUE)
  y <- rnbinom(500, mu = exp(.4 + X %*% c(.5, -.3, .2) + offset), size = 2)
  fixed <- l0glmFit(X, y, family = mgcv::nb(theta = 2), intercept = TRUE,
    weights = w, offset = offset, lambda = 0,  maxit = 150, tol = 1e-9)
  ref <- glm(y ~ X + offset(offset), family = MASS::negative.binomial(2), weights = w,
             control = glm.control(epsilon = 1e-12, maxit = 100))
  expect_equal(as.numeric(fixed$coefficients), unname(coef(ref)), tolerance = 1e-7)
  expect_equal(as.numeric(logLik(fixed)), as.numeric(logLik(ref)), tolerance = 1e-7)
  family <- mgcv::nb(theta = -1)
  fit <- l0glmFit(X, y, family = family, intercept = TRUE, weights = w, offset = offset,
    lambda = 0,  maxit = 200, tol = 1e-9)
  ref <- glm(y ~ X + offset(offset), family = MASS::negative.binomial(fit$theta), weights = w,
             control = glm.control(epsilon = 1e-12, maxit = 100))
  expect_equal(as.numeric(fit$coefficients), unname(coef(ref)), tolerance = 1e-6)
  expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)), tolerance = 1e-8)
  mu <- as.numeric(fitted(fit))
  pearson <- sum(w * (y-mu)^2 / (mu+mu^2/fit$theta))/mean(w[w>0])
  expect_equal(pearson, sum(w>0)-sum(coef(fit)!=0), tolerance = 1e-7)
  expect_equal(family$getTheta(TRUE), 1)
  expect_equal(fit$family$getTheta(TRUE), fit$theta)
  expect_equal(fit$metrics$nuisance, 1)
  small <- l0glmMetrics(fit, weights = rep(.001, 500), mu_floor = 0)
  expect_equal(small$mean_deviance, small$deviance/.5)
})

test_that("NB multichannel masks, theta and original-unit coefficients are retained", {
  skip_if_not_installed("mgcv")
  set.seed(1804)
  X <- matrix(rexp(200 * 6), 200, 6)
  truth <- matrix(c(2, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0), 6, 2)
  mu <- .2 + X %*% truth
  y <- matrix(rnbinom(length(mu), mu = mu, size = rep(c(2, 8), each = 200)), 200, 2)
  start <- matrix(1, 6, 2); start[5:6, ] <- 0; start[2, 1] <- 0
  fit <- l0glmFit(X, y, family = mgcv::nb(link = "identity"), offset = .2,
    nonnegative = TRUE, start = start, normalize = "rms", lambda = 2,
    grouped = TRUE,  maxit = 200, tol = 1e-8)
  expect_s4_class(fit$coefficients, "sparseMatrix")
  expect_equal(dim(fit$coefficients), c(6L, 2L))
  expect_true(all(fit$coefficients[c(5, 6), ] == 0))
  expect_equal(as.numeric(fit$coefficients[2, 1]), 0)
  expect_true(all(fit$coefficients >= 0))
  expect_length(fit$theta, 2)
  expect_true(all(is.finite(fit$theta) & fit$theta > 0))
  expect_equal(predict(fit), as.matrix(X %*% fit$coefficients) + .2, tolerance = 1e-12)
  expect_true(is.finite(fit$objective))
  expect_equal(fit$objective, -2 * fit$metrics$loglik + 2 * fit$metrics$selected)
})
