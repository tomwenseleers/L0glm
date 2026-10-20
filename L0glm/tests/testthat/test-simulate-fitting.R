test_that("simulated data and intercept truth can be fitted and scored directly", {
  families <- list(gaussian(), poisson(), binomial(), Gamma("log"), inverse.gaussian("log"))
  if (requireNamespace("mgcv", quietly=TRUE) && requireNamespace("survival", quietly=TRUE))
    families <- c(families, list(mgcv::nb(theta=3), mgcv::multinom(K=2), mgcv::cox.ph()))
  for (family in families) {
    kind <- L0glm_ext_kind(family)
    constant <- if (identical(kind, "multinomial")) c(.3, 0, -.3) else .5
    sim <- simul("normal", n=180, p=8, n_nonzero=2, family=family, intercept=constant,
      options=list(beta_range=c(-.6, .6)), seed=11)
    fit <- L0glm.fit(sim$X, sim$y, family=sim$family, weights=sim$weights,
      intercept=sim$intercept, nonnegative=FALSE, grouped=TRUE,
      maxit=100, tol=1e-8, parallelize=FALSE)
    expect_equal(dim(coef(fit)), dim(sim$beta_true))
    expect_true(all(is.finite(as.matrix(coef(fit)))))
    expect_true(is.finite(L0glm.metrics(fit, beta_true=sim$beta_true)$gic))
    expect_true(all(is.finite(predict(fit))))
  }
})

test_that("sparseglm honours simulation intercept metadata unless overridden", {
  for (intercept in list(FALSE, 0, 1.5)) {
    sim <- simul("normal", n=80, p=8, n_nonzero=2, intercept=intercept, seed=42)
    fit <- sparseglm(sim, method="L0glm", nlambda=2, lambda=c(10, 2))
    expect_identical(fit$data$intercept, sim$intercept)
    expect_equal(dim(coef(fit)$L0glm), dim(sim$beta_true))
    expect_identical(unname(fit$beta_true), unname(sim$beta_true))
    expect_identical(fit$timings$status, "ok")
  }
  sim <- simul("normal", n=80, p=8, n_nonzero=2, intercept=FALSE, seed=42)
  fit <- sparseglm(sim, method="L0glm", intercept=TRUE, nlambda=2, lambda=c(10, 2))
  expect_true(fit$data$intercept)
  expect_equal(nrow(coef(fit)$L0glm), 9)
})

test_that("special families retain valid shapes and explicit family constraints", {
  skip_if_not_installed("mgcv")
  sim <- simul("normal", n=40, p=4, family=mgcv::nb(theta=3))
  expect_equal(sim$theta, 3)
  expect_equal(simul("normal", n=40, p=4, family=mgcv::nb(theta=3),
    options=list(theta=3))$theta, 3)
  expect_error(simul("normal", n=40, p=4, family=mgcv::nb(theta=3),
    options=list(theta=2)), "conflicts")
  sim <- simulate_spike_train(n=40, p=40, k=5, family=mgcv::nb(theta=3, link="identity"))
  expect_equal(sim$theta, 3)
  expect_error(simul("normal", family="poisson", options=list(dispersion=2)), "fixed dispersion")
  sim <- simul("normal", n=30, p=3, beta=Matrix::Matrix(cbind(c(1, 0, 0), c(0, 1, 0)), sparse=TRUE),
    m=2, sparse=TRUE)
  expect_s4_class(sim$X, "sparseMatrix")
  expect_s4_class(sim$beta_true, "sparseMatrix")
})
