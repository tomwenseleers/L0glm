test_that("normal regression supports arbitrary dimensions, support and intercepts", {
  for (dimensions in list(c(20L, 70L), c(70L, 20L), c(1L, 1L))) {
    n <- dimensions[1L]; p <- dimensions[2L]; k <- min(4L, p)
    sim <- simul("normal", n=n, p=p, n_nonzero=k, m=2,
      intercept=c(0, -2), seed=35)
    expect_true(is.matrix(sim$X))
    expect_equal(dim(sim$X), c(n, p))
    expect_equal(dim(sim$beta_true), c(p+1L, 2L))
    expect_equal(sum(rowSums(abs(sim$beta_true[-1L, , drop=FALSE])) > 0), k)
    expect_equal(sim$beta_true[1L, ], c(y1=0, y2=-2))
    expect_equal(unname(sim$eta_true), unname(cbind(1, sim$X) %*% sim$beta_true))
    expect_equal(sim$y_true, sim$eta_true, ignore_attr=TRUE)
    expect_true(sim$intercept)
  }
  sim <- simul(5, n=20, p=25, seed=2, options=list(n_nonzero=3, beta_range=c(-2, -.5)))
  expect_equal(sum(sim$beta_true != 0), 3)
  expect_true(all(sim$beta_true[sim$beta_true != 0] >= -2 & sim$beta_true[sim$beta_true != 0] <= -.5))
  expect_false(sim$intercept)
  expect_equal(sim$y, simul(5, n=20, p=25, seed=2,
    options=list(n_nonzero=3, beta_range=c(-2, -.5)))$y)
})

test_that("simulation seeds are reproducible without disturbing caller randomness", {
  set.seed(381)
  before <- .Random.seed
  a <- simul("normal", n=15, p=40, seed=42)
  expect_identical(.Random.seed, before)
  b <- simul("normal", n=15, p=40, seed=42)
  expect_identical(a$X, b$X)
  expect_identical(a$beta_true, b$beta_true)
  expect_identical(a$y, b$y)
  expect_false(identical(a$y, simul("normal", n=15, p=40, seed=43)$y))
  expect_error(simul("normal", n=10, p=2, family=poisson("identity"), beta=c(1, 1)), "Invalid means")
  expect_identical(.Random.seed, before)
  invisible(simul("normal", n=15, p=40, seed=NULL))
  expect_false(identical(.Random.seed, before))
  rm(".Random.seed", envir=.GlobalEnv)
  on.exit(assign(".Random.seed", before, envir=.GlobalEnv), add=TRUE)
  invisible(simul("normal", n=15, p=40, seed=42))
  expect_false(exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE))
})

test_that("adding validation observations does not change the training problem", {
  for (method in c("normal", "abess", "bestsubset")) {
    training <- simul(method, n=50, p=15, n_nonzero=3, seed=16)
    validation <- simul(method, n=50, p=15, n_nonzero=3, seed=16, options=list(nval=20))
    for (name in c("X", "y", "beta_true", "dispersion"))
      expect_identical(training[[name]], validation[[name]])
  }
  sim <- simul("normal", n=20, p=3, family=binomial(), options=list(
    trials=rep(5, 20), nval=10, trialsval=rep(7, 10)))
  expect_true(all(sim$weights == 5))
  expect_true(all(sim$weightsval == 7))
})

test_that("standard GLM link functions preserve their specified means", {
  cases <- list(
    list(make=gaussian, links=c("identity", "log", "inverse"), eta=c(2, log(2), .5)),
    list(make=poisson, links=c("identity", "log", "sqrt"), eta=c(2, log(2), sqrt(2))),
    list(make=binomial, links=c("logit", "probit", "cauchit", "log", "cloglog"), eta=c(.3, .3, .3, -.3, .3)),
    list(make=Gamma, links=c("inverse", "identity", "log"), eta=c(.5, 2, log(2))),
    list(make=inverse.gaussian, links=c("1/mu^2", "inverse", "identity", "log"), eta=c(.25, .5, 2, log(2))))
  for (case in cases) for (j in seq_along(case$links)) {
    family <- case$make(case$links[j])
    sim <- simul("normal", n=20, p=2, beta=c(0, 0), intercept=case$eta[j], family=family)
    expect_equal(as.numeric(sim$y_true), rep(family$linkinv(case$eta[j]), 20))
    expect_true(all(is.finite(sim$y)))
  }
})

test_that("predictor covariance and population Gaussian SNR are correct", {
  b <- c(1, -2, 0, .5)
  for (structure in list(list(correlation="ar1", rho=.6),
      list(correlation="ar1", rho=-.6), list(correlation="equicorrelated", rho=.4),
      list(correlation="equicorrelated", rho=-.2), list(correlation="independent", rho=0))) {
    sim <- simul("normal", n=30000, p=4, beta=b, intercept=3,
      options=c(structure, list(snr=3, nval=200)), seed=18)
    Sigma <- if (structure$correlation == "ar1") structure$rho^abs(outer(1:4, 1:4, "-")) else {
      z <- matrix(structure$rho, 4, 4); diag(z) <- 1; z
    }
    expect_lt(max(abs(cov(sim$X)-Sigma)), .035)
    expect_equal(unname(sim$dispersion), as.numeric(crossprod(b, Sigma %*% b))/3)
    expect_equal(unname(sim$yval_true), unname(cbind(1, sim$Xval) %*% sim$beta_true))
    expect_equal(dim(sim$yval), c(200L, 1L))
    expect_lt(abs(var(as.numeric(sim$y-sim$y_true))/sim$dispersion-1), .04)
  }
})

test_that("all five best-subset patterns have the prescribed signal and SNR", {
  p <- 12L; k <- 4L; rho <- .4; snr <- 2
  patterns <- list(replace(numeric(p), round(seq(1, p, length.out=k)), 1),
    c(rep(1, k), rep(0, p-k)), c(seq(10, .5, length.out=k), rep(0, p-k)),
    c(-10, -6, -2, 2, 6, 10, rep(0, p-6)), c(rep(1, k), .5^seq_len(p-k)))
  Sigma <- rho^abs(outer(seq_len(p), seq_len(p), "-"))
  for (type in 1:5) {
    sim <- simul("bestsubset", n=80, p=p, n_nonzero=k, intercept=2,
      options=list(beta.type=type, rho=rho, snr=snr, nval=30))
    expect_equal(as.numeric(sim$beta_true), c(2, patterns[[type]]))
    expect_equal(unname(sim$dispersion), as.numeric(crossprod(patterns[[type]], Sigma %*% patterns[[type]]))/snr)
    expect_equal(unname(sim$y_true), unname(cbind(1, sim$X) %*% sim$beta_true))
  }
  expect_equal(sum(simul("bestsubset", p=5, n_nonzero=0, options=list(beta.type=2))$beta_true != 0), 0)
  expect_equal(sum(simul("bestsubset", p=5, n_nonzero=5, options=list(beta.type=5))$beta_true != 0), 5)
  expect_error(simul("bestsubset", p=5, options=list(beta.type=4)), "p>=6")
})

test_that("ABESS designs and coefficient magnitudes follow their specified setups", {
  n <- 200; p <- 30; k <- 4; a <- 5*sqrt(2*log(p)/n)
  sim <- simul("abess", n=n, p=p, n_nonzero=k, options=list(cortype=2, rho=.5))
  expect_equal(sim$covariance, list(type="ar1", rho=.5))
  expect_equal(sum(sim$beta_true != 0), k)
  expect_true(all(sim$beta_true[sim$beta_true != 0] >= a & sim$beta_true[sim$beta_true != 0] <= 100*a))
  for (name in c("binomial", "poisson")) {
    sim <- simul("abess", n=n, p=p, n_nonzero=k, family=name, intercept=1)
    b <- sim$beta_true[-1L, 1L]; b <- b[b != 0]
    expect_equal(length(b), k)
    expect_true(all(if (name == "binomial") b >= 2*a & b <= 10*a else abs(b) <= 2*a))
    expect_equal(unname(sim$y_true), unname(matrix(sim$family$linkinv(sim$eta_true), n)))
  }
  multi <- simul("abess", n=100, p=10, n_nonzero=3, family="mgaussian", options=list(y.dim=4))
  expect_equal(dim(multi$y), c(100L, 4L))
  expect_equal(length(multi$dispersion), 4L)
})

test_that("parametric and quasi-family draws have the correct conditional moments", {
  cases <- list(
    list(family=gaussian(), intercept=2, options=list(dispersion=.5), mu=2, variance=.5),
    list(family=poisson("identity"), intercept=3, options=list(), mu=3, variance=3),
    list(family=binomial(), intercept=qlogis(.3), options=list(trials=10), mu=.3, variance=.021),
    list(family=Gamma("log"), intercept=log(2), options=list(dispersion=.5), mu=2, variance=2),
    list(family=inverse.gaussian("log"), intercept=log(2), options=list(dispersion=.25), mu=2, variance=2),
    list(family=quasipoisson(), intercept=log(3), options=list(dispersion=2), mu=3, variance=6),
    list(family=quasibinomial(), intercept=qlogis(.3), options=list(trials=10, dispersion=2), mu=.3, variance=.042))
  if (requireNamespace("mgcv", quietly=TRUE)) cases <- c(cases,
    list(list(family=mgcv::nb(theta=2), intercept=log(3), options=list(), mu=3, variance=7.5)))
  for (case in cases) {
    sim <- simul("normal", n=60000, p=1, beta=0, family=case$family,
      intercept=case$intercept, options=case$options, seed=137)
    expect_equal(as.numeric(sim$y_true), rep(case$mu, 60000))
    expect_lt(abs(mean(sim$y)-case$mu), 6*sqrt(case$variance/60000))
    expect_lt(abs(var(as.numeric(sim$y))/case$variance-1), .08)
  }
})

test_that("multinomial simulations use CLR truth and valid count/factor responses", {
  skip_if_not_installed("mgcv")
  B <- rbind(c(1, -.5, 0), c(0, 0, 0), c(-1, .5, 0), c(0, 0, 0))
  sim <- simul("normal", n=150, p=4, beta=B, family=mgcv::multinom(K=2),
    intercept=c(1, 0, -.5), options=list(trials=3, nval=40), seed=11)
  expect_equal(dim(sim$beta_true), c(5L, 3L))
  expect_lt(max(abs(rowSums(sim$beta_true))), 1e-14)
  expect_lt(max(abs(rowSums(sim$y_true)-1)), 1e-14)
  expect_equal(rowSums(sim$y), rep(3, 150))
  expect_equal(rowSums(sim$yval), rep(3, 40))
  shifted <- simul("normal", n=150, p=4, beta=B+10, family=mgcv::multinom(K=2),
    intercept=c(1, 0, -.5)+7, options=list(trials=3, nval=40), seed=11)
  expect_equal(sim$y_true, shifted$y_true)
  factor_sim <- simul("normal", n=150, p=4, beta=B, family="multinomial", m=3,
    intercept=c(1, 0, -.5), options=list(response="factor"))
  expect_s3_class(factor_sim$y, "factor")
  expect_equal(levels(factor_sim$y), c("0", "1", "2"))
  null <- simul("normal", n=20, p=4, family="multinomial", beta=rep(1, 4), intercept=TRUE)
  expect_true(all(null$beta_true == 0))
  expect_equal(as.numeric(null$y_true), rep(1/3, 60))
})

test_that("Cox simulation handles a constant through the baseline hazard", {
  skip_if_not_installed("mgcv"); skip_if_not_installed("survival")
  args <- list(method="normal", n=5000, p=4, beta=c(.5, -.5, 0, 0), family="cox",
    options=list(weibull.shape=2, baseline.rate=.7, uniform.max=2), seed=45)
  plain <- do.call(simul, args)
  shifted <- do.call(simul, c(args, list(intercept=1.5)))
  expect_s3_class(shifted$y, "Surv")
  expect_false(shifted$intercept)
  expect_equal(dim(shifted$beta_true), c(4L, 1L))
  expect_equal(shifted$baseline.rate, .7*exp(1.5))
  expect_equal(shifted$event_time, plain$event_time*exp(-1.5/2))
  expect_equal(shifted$censor_time, plain$censor_time)
  expect_equal(as.numeric(shifted$y[, "status"]), as.numeric(shifted$event_time <= shifted$censor_time))
  exponential <- shifted$baseline.rate*as.numeric(shifted$y_true)*shifted$event_time^2
  expect_lt(abs(mean(exponential)-1), .04)
})

test_that("ordinal ABESS examples have ordered probabilities and thresholds", {
  sim <- simul("abess", n=100, p=8, beta=c(.5, rep(0, 7)), family="ordinal",
    m=4, intercept=.3, options=list(cutpoints=c(-2, 0, 2)))
  expect_true(is.ordered(sim$y))
  expect_equal(dim(sim$beta_true), c(9L, 1L))
  expect_equal(rowSums(sim$y_true), rep(1, 100))
  expect_equal(as.numeric(sim$y_true[, 1L]), plogis(-2-as.numeric(sim$eta_true)))
  expect_true(all(sim$y_true >= 0))
})

test_that("simulation validation rejects invalid domains and ambiguous options", {
  expect_error(simul(n=0), "n must")
  expect_error(simul("normal", p=4, n_nonzero=5), "n_nonzero")
  expect_error(simul("normal", options=list(typo=1)), "Unknown simulation option")
  expect_error(simul("normal", options=list(cortype=2, correlation="ar1")), "not both")
  expect_error(simul("normal", p=4, options=list(correlation="equicorrelated", rho=-.5)), "equicorrelation")
  expect_error(simul("normal", family=quasi()), "sampling distribution")
  expect_error(simul("normal", family=gaussian("log"), options=list(snr=3)), "identity link")
  expect_error(simul("normal", family=quasipoisson(), options=list(dispersion=.5)), "dispersion>=1")
  expect_error(simul("normal", family=quasibinomial(), options=list(dispersion=2)), "trial total")
  expect_error(simul("normal", options=list(beta_range=1)), "beta_range")
})
