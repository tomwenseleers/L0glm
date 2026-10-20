test_that("inverse-Gaussian deviance remains stable at extreme finite means", {
  y <- c(.2, 1, 4)
  eta <- c(-2, .5, 3)
  w <- c(.5, 2, 0)
  deviance <- function(y, eta, w) L0glm_irls_deviance(y, eta, w, inverse.gaussian("log"))
  expect_equal(deviance(y, eta, w),
    sum(inverse.gaussian("log")$dev.resids(y, exp(eta), w)))
  expect_equal(deviance(c(2, 3), rep(700, 2), c(1, 2)), 1/2+2/3)
  expect_identical(deviance(y, c(-1000, 0, 0), w), Inf)
  expect_identical(deviance(y, c(1000, 0, 0), w), Inf)
})

test_that("inverse-Gaussian backtracking checks the frozen penalized objective", {
  step <- function(beta, previous, lambda) L0glm_step(beta, previous,
    eta=beta, eta_previous=previous,
    objective=function(b, e) L0glm_irls_deviance(2, e, 1, inverse.gaussian("log")) + lambda*b^2)
  damped <- step(15, 0, 1)
  expect_true(damped$accepted)
  expect_gt(damped$halvings, 0)
  expect_lte(L0glm_irls_deviance(2, damped$eta, 1, inverse.gaussian("log"))+damped$beta^2, .5)
  # Shrinkage can legitimately increase deviance while decreasing its sum
  # with the current ridge penalty. A deviance-only check would reject this.
  shrink <- step(0, log(2), 10)
  expect_true(shrink$accepted)
  expect_identical(shrink$halvings, 0L)
  expect_equal(shrink$beta, 0)
})

test_that("high-dispersion inverse-Gaussian tuning retains finite models", {
  sim <- simul("normal", n=180, p=10, n_nonzero=3,
    family=inverse.gaussian("log"), intercept=.5,
    options=list(beta_range=c(-.8, .8)), seed=21)
  expect_no_warning(tuned <- L0glm.tune(sim$X, sim$y, family=sim$family,
    weights=sim$weights, intercept=TRUE, nonnegative=FALSE, tune="gic",
    nlambda=3, maxit=200, tol=1e-8, parallelize=FALSE, plot=FALSE))
  expect_true(all(is.finite(tuned$fit_best$coefficients)))
  expect_true(all(is.finite(fitted(tuned$fit_best))))
  expect_true(all(fitted(tuned$fit_best) > 0))
  expect_true(all(tuned$fit_best$convergence$converged))
})

test_that("inverse-Gaussian backtracking does not repeat ridge solves", {
  sim <- simul("normal", n=180, p=10, n_nonzero=3,
    family=inverse.gaussian("log"), intercept=.5,
    options=list(beta_range=c(-.8, .8)), seed=21)
  solves <- halvings <- 0L
  original_solve <- lmridge
  original_step <- L0glm_step
  testthat::local_mocked_bindings(
    lmridge=function(...) {
      solves <<- solves+1L
      original_solve(...)
    },
    L0glm_step=function(...) {
      answer <- original_step(...)
      halvings <<- halvings+answer$halvings
      answer
    })
  fit <- L0glm.fit(sim$X, sim$y, family=sim$family, weights=sim$weights,
    intercept=TRUE, lambda="gic", maxit=200, tol=1e-8,
    relaxed=FALSE, parallelize=FALSE)
  expect_gt(halvings, 0L)
  expect_equal(solves, sum(fit$convergence$iterations))
})

test_that("unpenalized inverse-Gaussian IRLS agrees with likelihood optimization", {
  sim <- simul("normal", n=300, p=4, n_nonzero=4,
    family=inverse.gaussian("log"), intercept=.3,
    options=list(beta_range=c(-.5, .5), dispersion=2), seed=1)
  X <- cbind(1, sim$X)
  y <- as.numeric(sim$y)
  w <- seq(.5, 1.5, length.out=length(y)); w[1:6] <- 0
  offset <- seq(-.2, .2, length.out=length(y))
  objective <- function(beta) L0glm_irls_deviance(y, as.numeric(offset+X%*%beta), w, sim$family)
  gradient <- function(beta) {
    inv <- exp(-as.numeric(offset+X%*%beta))
    as.numeric(2*crossprod(X, w*(1-y*inv)*inv))
  }
  reference <- optim(rep(0, ncol(X)), objective, gradient, method="BFGS",
    control=list(reltol=1e-13, maxit=1000))
  expect_equal(reference$convergence, 0)
  for (sparse in c(FALSE, TRUE)) {
    design <- if (sparse) Matrix::Matrix(sim$X, sparse=TRUE) else sim$X
    expect_no_warning(fit <- L0glm.fit(design, y, family=sim$family,
      weights=w, offset=offset, intercept=TRUE, lambda=0, thresh=0,
      maxit=200, tol=1e-9, parallelize=FALSE))
    expect_true(all(fit$convergence$converged))
    expect_equal(as.numeric(fit$coefficients), reference$par, tolerance=5e-5)
    expect_equal(objective(as.numeric(fit$coefficients)), reference$value, tolerance=1e-8)
  }
})

test_that("the initial inverse-Gaussian working fit cannot falsely converge", {
  y <- exp(seq(-2, 2, length.out=80))
  w <- seq(.5, 1.5, length.out=80)
  fit <- L0glm.fit(matrix(numeric(), 80, 0), y, weights=w,
    family=inverse.gaussian("log"), intercept=TRUE, dispersion=1,
    lambda=0, maxit=100, tol=1e-9, parallelize=FALSE)
  expect_gt(fit$convergence$iterations, 1L)
  expect_true(fit$convergence$converged)
  expect_equal(as.numeric(fit$coefficients), log(weighted.mean(y, w)), tolerance=1e-8)
})

test_that("grouped inverse-Gaussian channels preserve masks on both execution paths", {
  sim <- simul("normal", n=120, p=5, n_nonzero=2, m=3,
    family=inverse.gaussian("log"), intercept=.4,
    options=list(beta_range=c(-.4, .4), dispersion=c(.1, 1, 2)), seed=9)
  start <- matrix(.1, 6, 3); start[6, ] <- 0; start[2, 1] <- 0
  args <- list(X=sim$X, y=sim$y, family=sim$family, intercept=TRUE,
    start=start, lambda="gic", grouped=TRUE, relaxed=TRUE,
    group.maxit=2, maxit=200, tol=1e-8)
  expect_no_warning(serial <- do.call(L0glm.fit, c(args, list(parallelize=FALSE))))
  expect_no_warning(parallel <- do.call(L0glm.fit, c(args, list(parallelize=TRUE, ncores=2))))
  expect_equal(parallel$coefficients, serial$coefficients)
  expect_equal(parallel$dispersion, serial$dispersion)
  expect_equal(dim(serial$coefficients), c(6L, 3L))
  expect_true(all(serial$coefficients[6, ] == 0))
  expect_equal(as.numeric(serial$coefficients[2, 1]), 0)
  expect_true(all(is.finite(fitted(serial))))
})

test_that("inverse-Gaussian damping respects exclusions, constraints and refits", {
  sim <- simul("normal", n=250, p=8, n_nonzero=3,
    family=inverse.gaussian("log"), intercept=.5,
    options=list(beta_range=c(.2, .6), dispersion=2), seed=7)
  start <- rep(.1, 9); start[c(3, 5)] <- 0
  for (normalize in c("none", "rms")) {
    expect_no_warning(fit <- L0glm.fit(sim$X, sim$y, family=sim$family,
      intercept=TRUE, start=start, prior.mean=rep(.01, 9),
      nonnegative=TRUE, normalize=normalize, lambda="gic", relaxed=TRUE,
      maxit=200, tol=1e-8, parallelize=FALSE))
    expect_s4_class(fit$coefficients, "sparseMatrix")
    expect_equal(dim(fit$coefficients), c(9L, 1L))
    expect_true(all(fit$coefficients >= 0))
    expect_true(all(fit$coefficients[c(3,5), ] == 0))
    expect_true(all(is.finite(fitted(fit))))
  }
})
