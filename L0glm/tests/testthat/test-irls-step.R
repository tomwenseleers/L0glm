test_that("shared backtracking checks valid means and a frozen target across families", {
  families <- list(gaussian(), gaussian("log"), poisson("log"), poisson("identity"),
    binomial("logit"), binomial("log"), Gamma("log"), Gamma("inverse"),
    inverse.gaussian("log"), inverse.gaussian("1/mu^2"))
  for (family in families) {
    y <- if (family$family == "binomial") .7 else 2
    previous <- family$linkfun(y)
    proposal <- if (family$link %in% c("inverse", "1/mu^2", "log") && family$family != "binomial")
      previous + 100 else previous - 100
    score <- function(b, e) L0glm_irls_deviance(y, e, 1, family) + .01*b^2
    step <- L0glm_step(proposal, previous, score, eta = proposal, eta_previous = previous)
    expect_true(step$accepted, info = paste(family$family, family$link))
    expect_gt(step$halvings, 0)
    expect_lte(score(step$beta, step$eta), score(previous, previous) + 1e-10)
    expect_equal(step$eta, step$beta)
  }
})

test_that("shared backtracking preserves convex bounds and absorbing zeroes", {
  old <- c(1, 0, .1); proposal <- c(10, 0, 1)
  step <- L0glm_step(proposal, old,
    function(b, e) sum(b^2), eta = proposal, eta_previous = old)
  expect_true(step$accepted)
  expect_equal(step$beta[2], 0)
  expect_true(all(step$beta >= 0))
  expect_lte(sum(step$beta^2), sum(old^2) + 1e-10)
  failed <- L0glm_step(proposal, old, function(b, e) Inf,
    eta = proposal, eta_previous = old, baseline = Inf, maxit = 3L)
  expect_false(failed$accepted)
  expect_equal(failed$beta, old)
  expect_equal(failed$eta, old)
})

test_that("first and discrete-support proposals still require finite valid means", {
  family <- poisson("identity")
  score <- function(b, e) L0glm_irls_deviance(2, e, 1, family)
  finite <- L0glm_step(4, 2, score, 4, 2, baseline = Inf)
  expect_true(finite$accepted)
  expect_equal(finite$halvings, 0L)
  guarded <- L0glm_step(-4, 2, score, -4, 2, baseline = Inf)
  expect_true(guarded$accepted)
  expect_gt(guarded$halvings, 0)
  expect_gte(guarded$eta, 0)
})

test_that("an infeasible starting direction can retain the valid null model", {
  X <- matrix(c(1, 0, -1, 1, 0), ncol = 1L)
  expect_warning(fit <- L0glm.fit.internal(X, c(1, 2, 3, 4, 5), start = 1,
    family = poisson("identity")), "retaining previous")
  expect_equal(fit$coefficients, 0)
  expect_true(all(is.finite(fit$weights)))
  expect_false(fit$convergence$converged)
})

test_that("special-family selection retains tiny terms needed for valid means", {
  skip_if_not_installed("mgcv")
  X <- cbind(c(1, 1, 0), c(0, 0, 1))
  fit <- L0glm.fit(X, c(10, 10, 1), family = mgcv::nb(theta = 2, link = "identity"),
    start = 1, normalize = "none", nonnegative = TRUE, lambda = .001, thresh = .5)
  expect_true(all(fitted(fit) > 0))
  expect_true(is.finite(as.numeric(logLik(fit))))
})

test_that("feasible first solves retain eligible terms from an invalid start", {
  X <- matrix(c(1, -2, 1, 1, 1, 2), ncol = 2L, byrow = TRUE)
  y <- c(3, 6, 7)
  args <- list(X = X, y = y, family = poisson("identity"), start = 1,
    normalize = "none", nonnegative = TRUE, lambda = 1e-4,
    relaxed = FALSE, parallelize = FALSE)
  a <- do.call(L0glm.fit, c(args, list(openmp = FALSE)))
  b <- do.call(L0glm.fit, c(args, list(openmp = TRUE)))
  expect_true(all(coef(a) > 0))
  expect_true(all(fitted(a) > 0))
  expect_equal(coef(a), coef(b), tolerance = 1e-7)
  expect_lt(max(abs(fitted(a) - y)), .01)
})

test_that("native infeasible starts retain a finite anchor with a warning", {
  X <- matrix(c(1, 0, -1, 1, 0), ncol = 1L)
  expect_warning(fit <- L0nnpois(X, c(1, 2, 3, 4, 5), start = 1,
    normalize = "none", lambda = .01, relaxed = FALSE, parallelize = FALSE))
  expect_equal(as.numeric(coef(fit)), 0)
  expect_false(fit$convergence$converged)
})
