biglasso_fixture <- function() {
  set.seed(8501)
  X <- matrix(rnorm(120*40),120,40)
  list(X=X,y=2*X[,1]-2*X[,2]+rnorm(120),foldid=rep(1:3,40))
}

test_that("native biglasso paths retain the selected coefficients and predictions", {
  skip_if_not_installed("biglasso");skip_if_not_installed("bigmemory")
  d <- biglasso_fixture();M <- bigmemory::as.big.matrix(d$X)
  rows <- c(seq(2,120,2),seq(1,119,2))
  for(input in list(d$X[rows,],l0glmData(M)[rows,,drop=FALSE])) {
    fit <- sparseGLM(input,d$y[rows],method="biglasso",tune="bic",nlambda=7,
      intercept=TRUE,screen="none")
    expect_identical(fit$timings$status,"ok")
    grid <- fit$tuning[[1]]$candidates$parameter
    reference <- biglasso::biglasso(bigmemory::as.big.matrix(d$X[rows,]),d$y[rows],
      family="gaussian",lambda=grid,eps=1e-7,max.iter=10000)
    expected <- reference$beta[,fit$tuning[[1]]$index,drop=FALSE]
    expect_equal(as.numeric(coef(fit)[[1]]),as.numeric(expected),tolerance=1e-9)
    expect_equal(as.numeric(predict(fit,l0glmData(M)[1:5,,drop=FALSE])[[1]]),
      as.numeric(cbind(1,d$X[1:5,])%*%expected),tolerance=1e-9)
    expect_s4_class(coef(fit)[[1]],"sparseMatrix")
  }
})

test_that("file-backed folds agree with resident Gaussian, logistic and Cox fits", {
  skip_if_not_installed("biglasso");skip_if_not_installed("bigmemory")
  skip_if_not_installed("mgcv");skip_if_not_installed("survival")
  d <- biglasso_fixture();source <- l0glmData(bigmemory::as.big.matrix(d$X))
  for(family in list(gaussian(),binomial(),mgcv::cox.ph())) {
    y <- if(family$family=="gaussian")d$y else if(family$family=="binomial")
      rbinom(120,1,plogis(d$X[,1]))else {
        time <- rexp(120,exp(.5*d$X[,1]))
        survival::Surv(pmin(time,2),time<=2)
      }
    args <- list(y=y,family=family,method="biglasso",tune="cv",nlambda=5,
      foldid=d$foldid,screen="none")
    a <- do.call(sparseGLM,c(list(X=d$X),args))
    b <- do.call(sparseGLM,c(list(X=source),args))
    expect_identical(a$timings$status,"ok");expect_identical(b$timings$status,"ok")
    expect_equal(as.numeric(coef(a)[[1]]),as.numeric(coef(b)[[1]]),tolerance=1e-9)
    expect_equal(a$tuning,b$tuning,tolerance=1e-9)
  }
})

test_that("biglasso shares Gaussian grouping, offsets and zero-weight exclusions", {
  skip_if_not_installed("biglasso");skip_if_not_installed("bigmemory")
  d <- biglasso_fixture();Y <- cbind(d$y,-d$y+rnorm(120))
  off <- matrix(runif(240),120,2);w <- rep(1,120);w[seq(1,120,10)] <- 0
  for(grouped in c(TRUE,FALSE)) {
    args <- list(y=Y,method="biglasso",tune="cv",nlambda=5,foldid=d$foldid,
      weights=w,offset=off,grouped=grouped,screen="none")
    a <- do.call(sparseGLM,c(list(X=d$X),args))
    b <- do.call(sparseGLM,c(list(X=bigmemory::as.big.matrix(d$X)),args))
    expect_identical(a$timings$status,"ok");expect_identical(b$timings$status,"ok")
    expect_equal(coef(a),coef(b),tolerance=1e-9)
    expect_equal(a$tuning,b$tuning,tolerance=1e-9)
    expect_identical(a$timings$scope,if(grouped)"grouped"else"independent")
  }
})

test_that("native biglasso row views fit beyond the dense-read memory guard", {
  skip_if_not_installed("biglasso");skip_if_not_installed("bigmemory")
  set.seed(8502)
  X <- matrix(rnorm(200*2000),200,2000);y <- 3*X[,1]-2*X[,2]+rnorm(200)
  source <- l0glmData(bigmemory::as.big.matrix(X),control=list(max.cells=30000))
  expect_error(as.matrix(source),"max.cells")
  fit <- sparseGLM(source,y,method="biglasso",screen="none",tune="cv",
    nfolds=2,nlambda=4,options=list(lambda.min=.05,ncores=2L))
  expect_identical(fit$timings$status,"ok")
  expect_s3_class(fit$fits[[1]]$X,"L0glm_design")
  expect_equal(dim(coef(fit)[[1]]),c(2001L,1L))
  expect_length(predict(fit,source[1:3,,drop=FALSE])[[1]],3L)
  restored <- unserialize(serialize(source,NULL))
  expect_true(bigmemory::is.nil(restored$object@address))
  again <- sparseGLM(restored,y,method="biglasso",screen="none",tune="bic",
    lambda=fit$tuning[[1]]$candidates$parameter)
  expect_identical(again$timings$status,"ok")
  expect_equal(predict(again,restored[1:3,,drop=FALSE])[[1]],
    predict(again,source[1:3,,drop=FALSE])[[1]])
})

test_that("source imputation is learned separately within each biglasso fold", {
  skip_if_not_installed("biglasso");skip_if_not_installed("bigmemory")
  d <- biglasso_fixture();X <- d$X
  X[cbind(c(1,4,9),c(1,2,3))] <- NA_real_
  fit <- sparseGLM(l0glmData(bigmemory::as.big.matrix(X),missing="mean"),d$y,
    method="biglasso",screen="none",tune="cv",foldid=d$foldid,nlambda=4)
  expect_identical(fit$timings$status,"ok")
  grid <- fit$tuning[[1]]$candidates$parameter
  for(k in 1:3) {
    train <- d$foldid!=k;A <- X;means <- colMeans(X[train,],na.rm=TRUE)
    where <- which(is.na(A),arr.ind=TRUE);A[where] <- means[where[,2]]
    ref <- biglasso::biglasso(bigmemory::as.big.matrix(A[train,]),d$y[train],
      family="gaussian",lambda=grid,eps=1e-7,max.iter=10000)
    expected <- vapply(seq_along(grid),function(j)
      mean((d$y[!train]-drop(cbind(1,A[!train,])%*%ref$beta[,j]))^2),0.0)
    expect_equal(fit$tuning[[1]]$cv$deviance[,k],expected,tolerance=1e-9)
  }
})

test_that("biglasso capability checks preserve common constraints and options", {
  skip_if_not_installed("biglasso");skip_if_not_installed("bigmemory")
  d <- biglasso_fixture()
  for(extra in list(list(intercept=FALSE),list(nonnegative=TRUE),list(weights=2),
      list(family=poisson(identity)),list(family=binomial(),offset=.2))) {
    y <- if(identical(extra$family$family,"poisson"))rpois(120,2)else
      if(identical(extra$family$family,"binomial"))rbinom(120,1,.5)else d$y
    expect_warning(fit <- do.call(sparseGLM,c(list(X=d$X,y=y,method="biglasso",
      nlambda=3,screen="none"),extra)),"unavailable")
    expect_identical(fit$timings$status,"unsupported")
  }
  expect_warning(fit <- sparseGLM(d$X,d$y,method="biglasso",screen="none",
    options=list(row.idx=1:30)),"failed")
  expect_match(fit$timings$note,"row.idx",fixed=TRUE)
})
