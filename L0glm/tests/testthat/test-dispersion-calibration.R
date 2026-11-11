test_that("wide signed Gaussian fits use the row pilot and BIC only for an omitted penalty", {
  set.seed(1701)
  n <- 60L; p <- 120L; m <- 8L
  X <- matrix(rnorm(n*p),n,p)
  B <- matrix(0,p,m)
  B[1:8,] <- matrix(rnorm(8*m,sd=1.5),8,m)
  Y <- X%*%B + matrix(rnorm(n*m),n,m)
  args <- list(X=X,y=Y,intercept=TRUE,parallelize=FALSE,local.search=FALSE,screen="none")
  fit <- do.call(l0glmFit,args)
  score <- l0glmMetrics(fit,beta_true=rbind(0,B))
  expect_equal(fit$init.lambda,log(n))
  expect_identical(fit$dispersionCalibration$method,"adaptive ridge")
  expect_true(any(fit$dispersionCalibration$reliable))
  expect_gte(score$TP,7)
  expect_equal(score$df,sum(coef(fit)!=0)+sum(coef(fit)[fit$always.include,,drop=FALSE]==0)+m)
  expect_equal(coef(do.call(l0glmFit,c(args,list(lambda=log(n))))),coef(fit))
  expect_equal(do.call(l0glmFit,c(args,list(lambda="gic")))$init.lambda,presetLambda("gic",p,n))

  supplied <- do.call(l0glmFit,c(args,list(start=c(1,1,rep(0,p-1L)),dispersion=2)))
  expect_null(supplied$dispersionCalibration)
  expect_equal(supplied$dispersion,rep(2,m))
  expect_true(all(coef(supplied)[-c(1L,2L),]==0))
  disabled <- do.call(l0glmFit,c(args,list(control.family=list(dispersion.calibration=FALSE))))
  expect_null(disabled$dispersionCalibration)
  expect_error(do.call(l0glmFit,c(args,list(control.family=list(dispersion.calibration="never")))),
    "dispersion.calibration")
  screened <- do.call(l0glmFit,modifyList(args,list(screen="SIS",nsis=8L)))
  expect_equal(screened$init.lambda,presetLambda("gic",p,n))
  expect_null(screened$dispersionCalibration)
})

test_that("native and R group passes use the same final Pearson residual scale", {
  set.seed(1702)
  n <- 70L
  X <- matrix(rnorm(n*5L),n,5L)
  Y <- X[,1:2] %*% matrix(c(2,-1,-2,1.5),2L) + matrix(rnorm(n*2L),n,2L)
  w <- c(0,runif(n-1L,.5,2))
  args <- list(X=X,y=Y,weights=w,group=seq_len(ncol(X)),intercept=TRUE,lambda=0,
    grouped=TRUE,dispersion=c(-1,2),parallelize=FALSE,local.search=FALSE)
  native <- do.call(l0glmFit,args)
  reference <- do.call(l0glmFit,c(args,list(openmp=FALSE)))
  k <- sum(coef(native)[,1]!=0) + sum(coef(native)[native$always.include,1]==0)
  target <- sum(w*(Y[,1]-fitted(native)[,1])^2)/(sum(w>0)-k)
  expect_equal(native$dispersion,c(target,2))
  expect_equal(coef(native),coef(reference),tolerance=1e-8)
  expect_equal(native$dispersion,reference$dispersion,tolerance=1e-8)
})

test_that("the negative-binomial pilot starts shape estimation and honors known shapes", {
  skip_if_not_installed("mgcv")
  set.seed(1703)
  n <- 100L
  X <- matrix(rnorm(n*6L),n,6L)
  Y <- vapply(c(.45,.7),function(b) rnbinom(n,mu=exp(1+b*X[,1]),size=3),numeric(n))
  args <- list(X=X,y=Y,intercept=TRUE,lambda="gic",parallelize=FALSE,local.search=FALSE)
  free <- do.call(l0glmFit,c(args,list(family=mgcv::nb())))
  expect_identical(free$dispersionCalibration$method,"Poisson pilot")
  expect_equal(free$initialTheta,free$dispersionCalibration$calibrated)
  expect_true(all(free$theta.estimated))
  expect_true(any(abs(free$theta-free$initialTheta)>1e-4))
  mu <- fitted(free)
  expect_equal(as.numeric(logLik(free)),sum(dnbinom(Y,mu=mu,size=rep(free$theta,each=n),log=TRUE)))
  expect_equal(attr(logLik(free),"df"),
    sum(coef(free)!=0)+sum(coef(free)[free$always.include,,drop=FALSE]==0)+ncol(Y))
  fixed <- do.call(l0glmFit,c(args,list(family=mgcv::nb(theta=3))))
  expect_null(fixed$dispersionCalibration)
  expect_equal(fixed$theta,rep(3,2L))
  expect_false(any(fixed$theta.estimated))
  warm <- do.call(l0glmFit,c(args,list(family=mgcv::nb(theta=-5))))
  expect_null(warm$dispersionCalibration)
  expect_equal(warm$initialTheta,rep(5,2L))
})

test_that("the NB shape equation uses the residual count and respects its bounds", {
  set.seed(1704)
  n <- 120L; mu <- runif(n,2,10); y <- rnbinom(n,mu=mu,size=3)
  w <- runif(n,.5,2); w[1] <- 0
  control <- L0glm_ext_control(list())
  k <- 4L
  shape <- L0glm_nb_pearson_shape(y,mu,w,k,control)
  expect_equal(sum(w*(y-mu)^2/(mu+mu^2/shape))/mean(w[w>0]),sum(w>0)-k,tolerance=1e-6)
  expect_equal(L0glm_nb_pearson_shape(y,mu,7*w,k,control),shape)
  expect_equal(L0glm_nb_pearson_shape(mu,mu,w,k,control),control$theta.max)
  expect_true(is.na(L0glm_nb_pearson_shape(y,mu,w,sum(w>0),control)))
})

test_that("the constrained relaxed NB fallback respects bounds and eligibility", {
  skip_if_not_installed("mgcv")
  set.seed(1705)
  n <- 70L
  X <- cbind(1,matrix(runif(n*3L),n,3L))
  B <- matrix(c(1,.8,0,0),4,1)
  mu <- as.numeric(X%*%B)+.2
  data <- list(kind="nb",Z=Matrix::Matrix(X,sparse=TRUE),m=1L,
    y=matrix(rnbinom(n,mu=mu,size=3),n,1L),offset=matrix(.2,n,1L),
    weights=matrix(1,n,1L),family=mgcv::nb(link="identity"),
    lower=rep(0,4L),prior=matrix(0,4L,1L))
  D <- matrix(c(0,0,0,Inf),4L,1L)
  state <- list(theta=3)
  objective <- function(C) L0glm_ext_loss(data,C,state)
  bad <- B; bad[1L] <- 100
  expect_gt(objective(bad),objective(B))
  result <- L0glm_nb_refit_step(data,bad,B,D,state,objective)
  expect_true(result$accepted)
  expect_true(all(result$beta>=0))
  expect_equal(result$beta[4L],0)
  expect_lt(objective(result$beta),objective(B))
})
