validation_family_sample <- function(family, seed = 9201L) {
  set.seed(seed)
  X <- matrix(rnorm(120 * 6), 120, 6)
  V <- matrix(rnorm(60 * 6), 60, 6)
  a <- if (family$family %in% c("binomial", "quasibinomial"))
    if (family$link == "log") -1 else if (family$link == "identity") .5 else .1 else 1.5
  b <- if (family$family %in% c("binomial", "quasibinomial") && family$link == "identity") .03 else .1
  eta <- a + b * X[, 1]; ev <- a + b * V[, 1]
  draw <- function(eta) {
    mu <- family$linkinv(eta)
    switch(family$family,
      gaussian = mu + rnorm(length(mu), sd = .1),
      poisson = rpois(length(mu), mu), quasipoisson = rpois(length(mu), mu),
      binomial = rbinom(length(mu), 1, mu), quasibinomial = rbinom(length(mu), 1, mu),
      Gamma = rgamma(length(mu), 10, scale = mu / 10),
      inverse.gaussian = statmod::rinvgauss(length(mu), mean = mu, dispersion = .1))
  }
  list(X = X, V = V, y = draw(eta), v = draw(ev))
}

test_that("validation evaluates ordinary family deviances for every standard link", {
  families <- c(lapply(c("identity", "log", "inverse"), gaussian),
    lapply(c("identity", "log", "sqrt"), poisson),
    lapply(c("logit", "probit", "cloglog", "cauchit", "log"), binomial),
    lapply(c("identity", "log", "inverse"), Gamma),
    lapply(c("identity", "log", "inverse", "1/mu^2"), inverse.gaussian),
    list(quasipoisson(), quasibinomial()))
  for (family in families) {
    d <- validation_family_sample(family)
    fit <- L0glm.fit(d$X, d$y, family = family, intercept = TRUE, lambda = "gic",
      dispersion = 1, relaxed = TRUE, maxit = 100, parallelize = FALSE,
      validation = list(X = d$V, y = d$v))
    mu <- as.numeric(predict(fit, d$V))
    if (family$family %in% c("poisson", "quasipoisson", "Gamma", "inverse.gaussian"))
      mu <- pmax(mu, 1e-8)
    expected <- mean(family$dev.resids(d$v, mu, rep(1, length(mu))))
    info <- paste(family$family, family$link)
    expect_equal(fit$validation$deviance, expected, tolerance = 1e-8, info = info)
    expect_lte(expected, min(fit$validation$starts$validation_deviance) + 1e-10)
    expect_equal(dim(coef(fit)), c(7L, 1L), info = info)
  }
})

test_that("special-family validation retains training parameters and uses fixed references", {
  skip_if_not_installed("mgcv"); skip_if_not_installed("survival")
  d <- validation_family_sample(poisson())
  set.seed(9202)
  y <- rnbinom(120, mu = exp(.4 + d$X[, 1]), size = 2)
  v <- rnbinom(60, mu = exp(.4 + d$V[, 1]), size = 2)
  for (theta in c(-1, 2)) {
    family <- mgcv::nb(theta = theta)
    fit <- L0glm.fit(d$X, y, family = family, intercept = TRUE, lambda = "gic",
      relaxed = TRUE, maxit = 80, parallelize = FALSE, validation = list(X = d$V, y = v))
    mu <- as.numeric(predict(fit, d$V))
    saturated <- if (theta < 0) dpois(v, v, log = TRUE) else
      dnbinom(v, mu = v, size = theta, log = TRUE)
    expected <- 2 * mean(saturated - dnbinom(v, mu = pmax(mu, 1e-8), size = fit$theta, log = TRUE))
    expect_equal(fit$validation$deviance, expected, tolerance = 1e-8)
    expect_true(is.finite(fit$theta) && fit$theta > 0)
    for (stage in c("selection", "relaxed")) if (!is.null(fit$validation[[stage]]))
      expect_lte(fit$validation[[stage]]$deviance,
        min(fit$validation[[stage]]$history$deviance) + 1e-8)
  }
  labels <- function(n) factor(sample(letters[1:3], n, TRUE), levels = letters[1:3])
  ym <- labels(120); vm <- labels(60)
  mn <- L0glm.fit(d$X, ym, family = mgcv::multinom(2), lambda = "gic",
    relaxed = TRUE, maxit = 50, validation = list(X = d$V, y = vm))
  expect_equal(mn$validation$deviance,
    L0glm.metrics(mn, newdata = d$V, y = vm)$mean_deviance, tolerance = 1e-8)
  expect_equal(as.numeric(rowSums(as.matrix(coef(mn)))), rep(0, 7), tolerance = 1e-12)
  time <- rexp(120, exp(d$X[, 1])); tv <- rexp(60, exp(d$V[, 1]))
  yc <- survival::Surv(pmin(time, 2), time <= 2)
  vc <- survival::Surv(pmin(tv, 2), tv <= 2)
  cx <- L0glm.fit(d$X, yc, family = mgcv::cox.ph(), lambda = "gic",
    relaxed = TRUE, maxit = 50, validation = list(X = d$V, y = vc, rule = "1se"))
  expect_equal(cx$validation$deviance,
    L0glm.metrics(cx, newdata = d$V, y = vc)$mean_deviance, tolerance = 1e-8)
  expect_equal(dim(coef(cx)), c(6L, 1L))
})
