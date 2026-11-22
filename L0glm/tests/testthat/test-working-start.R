test_that("the first Poisson step uses supplied full means on dense and sparse designs", {
  X <- cbind(c(1, 1, 2, 0, 1, 2), c(0, 1, 1, 2, 1, 0))
  B <- cbind(c(4, 3), c(5, 2))
  off <- seq(1, 2, length.out = nrow(X))
  mu <- X %*% B + off
  Y <- round(mu + c(-1, 2, 0, -1, 1, 0))
  baseline <- sapply(1:2, function(j) as.numeric(stats::lm.wfit(X, Y[,j] - off, 1 / mu[,j])$coef))
  for (design in list(X, methods::as(Matrix::Matrix(X, sparse = TRUE), "CsparseMatrix"))) {
    fit <- function(...) l0glmFit(design, Y, start=B, offset=off, family=poisson("identity"),
      nonnegative=TRUE, normalize="none", lambda=0, group.maxit=1L, maxit=1L,
      relaxed=FALSE, local.search=FALSE, screen="none", ncores=2L, ...)
    expect_equal(unname(as.matrix(coef(fit(mustart=mu)))), baseline, tolerance=1e-9)
    expect_equal(coef(fit(mustart=0, etastart=mu)), coef(fit(mustart=mu)))
    expect_equal(coef(fit(mustart=mu[,1])), coef(fit(mustart=cbind(mu[,1],mu[,1]))))
    expect_error(fit(mustart=0), "mean domain")
    expect_error(fit(etastart=c(2,3)), "length n")
    expect_error(fit(mustart=matrix(1,6,3)), "n x m")
  }
})

test_that("mean and predictor starts reproduce a glm.fit Fisher step with an offset", {
  X <- cbind(Intercept=1, x=seq(-1,1,length.out=20))
  off <- rep(c(0.2, 0.4), 10)
  eta <- c(X %*% c(1,0.3)) + off
  Y <- round(exp(eta) + rep(c(0.4,-1,2,0),5))
  reference <- suppressWarnings(stats::glm.fit(X, Y, offset=off, family=poisson(),
    etastart=eta, control=stats::glm.control(maxit=1L)))$coefficients
  for (native in c(TRUE,FALSE)) {
    run <- function(...) l0glmFit(X,Y, family=poisson(), intercept=FALSE, nonnegative=FALSE,
      start=c(1,0.3), offset=off, normalize="none", lambda=0, maxit=1,
      relaxed=FALSE, local.search=FALSE,ncores=1L,control.solver=list(.native=native), ...)
    expect_equal(as.numeric(coef(run(mustart=exp(eta)))), unname(reference), tolerance=1e-8)
    expect_equal(coef(run(etastart=eta)), coef(run(mustart=exp(eta))))
  }
})

test_that("grouped predictors and joint likelihoods accept corresponding observation starts", {
  skip_if_not_installed("mgcv")
  set.seed(314)
  X <- matrix(rnorm(180),60,3)
  counts <- rpois(60,3)
  eta <- seq(-0.4,0.5,length.out=60)
  run <- function(y, family, ...) l0glmFit(X,y, family=family, nonnegative=FALSE,
    lambda=0.1, start=1, maxit=2L, ncores=1L, relaxed=FALSE, local.search=FALSE, ...)
  expect_equal(coef(run(counts,poisson(),etastart=eta,group=1:3)),
    coef(run(counts,poisson(),mustart=exp(eta),group=1:3)))
  expect_equal(coef(run(counts,mgcv::nb(theta=2),etastart=eta)),
    coef(run(counts,mgcv::nb(theta=2),mustart=exp(eta))))
  status <- cbind(seq_len(60), rep(0:1,30))
  expect_equal(coef(run(status,mgcv::cox.ph(),etastart=eta)),
    coef(run(status,mgcv::cox.ph(),mustart=exp(eta))))
  z <- cbind(eta, -eta, 0)
  p <- exp(z)/rowSums(exp(z))
  y <- factor(rep(1:3,20))
  expect_equal(coef(run(y,mgcv::multinom(K=2),etastart=z)),
    coef(run(y,mgcv::multinom(K=2),mustart=p)),tolerance=1e-8)
  expect_error(run(y,mgcv::multinom(K=2),mustart=p/2), "summing to one")
})
