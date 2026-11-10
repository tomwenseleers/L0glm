calibration_fixture <- function() {
  set.seed(718)
  n <- 100L; p <- 500L; m <- 40L; k <- 40L
  active <- sample.int(p,k)
  B <- matrix(0,p,m);time <- seq(0,2*pi,length.out=m)
  for(j in active)B[j,] <- sin((1+rnorm(1))*time+3*rnorm(1))
  X <- matrix(rnorm(n*p),n,p)
  y <- X%*%B+matrix(rnorm(n*m),n,m)
  V <- matrix(rnorm(300*p),300,p)
  list(X=X,y=y,B=B,V=V,mu=V%*%B)
}

test_that("channel splitting is reproducible, balanced and preserves the RNG", {
  set.seed(22);rng <- .Random.seed
  fold <- L0glm_calibration_fold(11L)
  expect_identical(.Random.seed,rng)
  expect_identical(fold,L0glm_calibration_fold(11L))
  expect_lte(diff(range(table(fold))),1)
  expect_setequal(unique(fold),1:2)
})

test_that("calibration improves the sparse multitask fit without changing its public specification", {
  a <- calibration_fixture()
  options <- list(X=a$X,y=a$y,lambda="gic",intercept=TRUE,
    parallelize=FALSE,nthreads=1L,local.search=FALSE)
  old <- do.call(l0glmFit,c(options,list(control.family=list(dispersion.calibration=FALSE))))
  rng <- .Random.seed
  fit <- do.call(l0glmFit,options)
  expect_identical(.Random.seed,rng)
  d <- fit$dispersionCalibration
  expect_true(d$applied)
  expect_identical(d$method,"cross-channel")
  expect_true(all(d$residualDf>=20L))
  expect_equal(fit$dispersion,d$calibrated)
  expect_true(all(fit$dispersion.estimated))
  expect_true(all(fit$dispersion>0 & is.finite(fit$dispersion)))
  expect_equal(dim(coef(fit)),c(501L,40L))
  expect_s4_class(coef(fit),"sparseMatrix")
  expect_lt(d$objective[["after"]],d$objective[["before"]])
  mse <- function(f)mean((as.matrix(predict(f,a$V))-a$mu)^2)
  expect_lt(mse(fit),.5*mse(old))
  truth <- rowSums(a$B!=0)>0
  recovered <- Matrix::rowSums(coef(fit)[-1,,drop=FALSE]!=0)>0
  expect_gte(sum(truth&recovered),30)
  expect_lte(sum(!truth&recovered),3)
  expect_null(fit$specification$start)
  expect_null(fit$specification$dispersion)
  expect_null(fit$specification$control.solver$.dispersion.fixed)
  metrics <- l0glmMetrics(fit)
  expect_equal(metrics$df,sum(coef(fit)!=0)+NCOL(a$y))
  expect_equal(as.numeric(GIC(fit)),metrics$gic)
  duplicate <- do.call(l0glmFit,c(list(X=a$X,y=a$y),fit$specification))
  expect_identical(coef(duplicate),coef(fit))
  expect_identical(duplicate$dispersion,fit$dispersion)
})

test_that("calibrated fits preserve response units, predictor units and sparse storage", {
  a <- calibration_fixture()
  fit <- function(X,y)l0glmFit(X,y,lambda="gic",intercept=TRUE,
    parallelize=FALSE,nthreads=1L,local.search=FALSE)
  base <- fit(a$X,a$y)
  units <- rep(c(.2,1,5),length.out=ncol(a$X))
  changed <- fit(t(t(a$X)*units),a$y*3)
  expect_true(changed$dispersionCalibration$applied)
  expect_equal(as.matrix(coef(changed))*c(1,units)/3,as.matrix(coef(base)),tolerance=1e-7)
  expect_equal(changed$dispersion/9,base$dispersion,tolerance=1e-7)
  sparse <- fit(Matrix::Matrix(a$X,sparse=TRUE),a$y)
  expect_s4_class(sparse$X,"sparseMatrix")
  expect_equal(as.matrix(coef(sparse)),as.matrix(coef(base)),tolerance=1e-7)
  expect_equal(sparse$dispersion,base$dispersion,tolerance=1e-7)
})

test_that("local refinement retains calibrated noise and its diagnostics", {
  a <- calibration_fixture()
  fit <- l0glmFit(a$X,a$y,lambda="gic",intercept=TRUE,parallelize=FALSE,nthreads=1L,
    local.search=TRUE,control.local=list(maxit=1L,max.evals=3L,top=1L,drop=1L,max.swap=1L))
  expect_true(fit$dispersionCalibration$applied)
  expect_true(fit$local.search$applied)
  expect_equal(fit$dispersion,fit$dispersionCalibration$calibrated)
  expect_true(all(fit$dispersion.estimated))
  expect_equal(fit$local.search$common.dispersion,fit$dispersion)
  expect_null(fit$specification$dispersion)
})

test_that("calibration preserves compact weights and offsets in serial and parallel fits", {
  a <- calibration_fixture();n <- nrow(a$X);m <- ncol(a$y)
  w <- seq(.5,2,length.out=n);w[1:5] <- 0;offset <- seq(-.2,.2,length.out=n)
  mean <- a$X%*%a$B
  y <- mean+(a$y-mean)/sqrt(pmax(w,.5))+offset
  args <- list(X=a$X,y=y,weights=w,offset=offset,lambda="gic",intercept=TRUE,
    parallelize=FALSE,nthreads=1L,local.search=FALSE)
  compact <- do.call(l0glmFit,args)
  expect_true(compact$dispersionCalibration$applied)
  expect_equal(dim(compact$prior.weights),c(n,1L))
  expect_equal(dim(compact$offset),c(n,1L))
  args$weights <- matrix(w,n,m);args$offset <- matrix(offset,n,m)
  expanded <- do.call(l0glmFit,args)
  expect_identical(coef(compact),coef(expanded))
  expect_identical(compact$dispersion,expanded$dispersion)
  args$parallelize <- TRUE;args$ncores <- 2L
  parallel <- do.call(l0glmFit,args)
  expect_equal(as.matrix(coef(parallel)),as.matrix(coef(compact)),tolerance=1e-9)
  expect_equal(parallel$dispersion,compact$dispersion,tolerance=1e-9)
})

test_that("single-channel fits retain the existing estimator across estimated-scale families", {
  set.seed(719);X <- matrix(rnorm(100*6),100,6)
  for(kind in c("gaussian","Gamma","inverse.gaussian","quasipoisson","quasibinomial","nb")) {
    w <- NULL
    if(kind=="gaussian") {family<-gaussian();y<-X[,1]+rnorm(100)}
    if(kind=="Gamma") {family<-Gamma(log);mu<-exp(.3+.2*X[,1]);y<-rgamma(100,5,scale=mu/5)}
    if(kind=="inverse.gaussian") {family<-inverse.gaussian(log);mu<-exp(.3+.2*X[,1]);y<-statmod::rinvgauss(100,mean=mu,dispersion=.1)}
    if(kind=="quasipoisson") {family<-quasipoisson(log);y<-3*rpois(100,exp(1+.2*X[,1]))}
    if(kind=="quasibinomial") {family<-quasibinomial();w<-rep(20,100);y<-rbinom(100,20,plogis(X[,1]))/20}
    if(kind=="nb") {skip_if_not_installed("mgcv");family<-mgcv::nb(theta=-1);y<-rnbinom(100,mu=exp(1+.2*X[,1]),size=3)}
    args <- list(X=X,y=y,weights=w,family=family,lambda="gic",intercept=TRUE,
      parallelize=FALSE,local.search=FALSE,maxit=60L)
    old <- do.call(l0glmFit,c(args,list(control.family=list(dispersion.calibration=FALSE))))
    fit <- do.call(l0glmFit,c(args,list(control.family=list(dispersion.calibration=TRUE))))
    expect_identical(coef(fit),coef(old),info=kind)
    expect_identical(fit$dispersion,old$dispersion,info=kind)
    expect_identical(fit$theta,old$theta,info=kind)
    expect_identical(fit$dispersionCalibration$reason,"single_channel",info=kind)
  }
})

test_that("explicit choices, independent channels and unsupported fits use the safe fallback", {
  set.seed(720);X <- matrix(rnorm(60*8),60,8);y <- cbind(X[,1]+rnorm(60),X[,1]+rnorm(60))
  fit <- l0glmFit(X,y,lambda="gic",parallelize=FALSE,local.search=FALSE,
    control.family=list(dispersion.calibration=FALSE))
  original <- fit
  scenarios <- list(
    list(start=1),list(prior.mean=0),list(dispersion=-1),list(validation=list()),
    list(control.solver=list(.nuisance.df=1)))
  for(extra in scenarios) {
    args <- c(list(control.family=list(dispersion.calibration=TRUE)),extra)
    answer <- L0glm_calibrated_fit(fit,args)
    expect_false(answer$dispersionCalibration$applied)
    expect_identical(coef(answer),coef(original))
    expect_identical(answer$dispersion,original$dispersion)
  }
  independent <- fit;independent$grouped <- FALSE
  expect_identical(L0glm_calibrated_fit(independent,list())$dispersionCalibration$reason,"independent_supports")
  local_mocked_bindings(L0glm_channel_calibration=function(...)stop("unavailable calibration"),.package="L0glm")
  answer <- L0glm_calibrated_fit(fit,list(control.family=list(dispersion.calibration=TRUE)))
  expect_identical(coef(answer),coef(fit))
  expect_identical(answer$dispersionCalibration$reason,"calibration_failed")
  expect_match(answer$dispersionCalibration$details,"unavailable")
})

test_that("a healthy wide model and fixed-noise families keep their fitted values", {
  set.seed(721);X <- matrix(rnorm(120*150),120,150)
  y <- outer(X[,1],seq(.5,2,length.out=8))+matrix(rnorm(120*8,sd=.2),120,8)
  args <- list(X=X,y=y,lambda="gic",parallelize=FALSE,local.search=FALSE)
  old <- do.call(l0glmFit,c(args,list(control.family=list(dispersion.calibration=FALSE))))
  fit <- do.call(l0glmFit,args)
  expect_identical(coef(fit),coef(old));expect_identical(fit$dispersion,old$dispersion)
  expect_false(isTRUE(fit$dispersionCalibration$applied))
  args$dispersion <- 1
  fixed <- do.call(l0glmFit,args)
  expect_null(fixed$dispersionCalibration)
  expect_equal(fixed$dispersion,rep(1,8))
  poisson <- l0glmFit(abs(X[,1:6]),matrix(rpois(120*2,2),120,2),family=poisson(identity),
    nonnegative=TRUE,parallelize=FALSE,local.search=FALSE)
  expect_null(poisson$dispersionCalibration)
})

test_that("substantial cross-channel residual dependence retains the original fit", {
  set.seed(725);X <- matrix(rnorm(60*120),60,120)
  y <- matrix(rnorm(60),60,12)+matrix(rnorm(60*12,sd=.001),60,12)
  args <- list(X=X,y=y,lambda="gic",intercept=TRUE,parallelize=FALSE,local.search=FALSE)
  old <- do.call(l0glmFit,c(args,list(control.family=list(dispersion.calibration=FALSE))))
  fit <- do.call(l0glmFit,c(args,list(control.family=list(dispersion.calibration=TRUE))))
  expect_false(fit$dispersionCalibration$applied)
  expect_identical(fit$dispersionCalibration$reason,"dependent_channels")
  expect_identical(coef(fit),coef(old));expect_identical(fit$dispersion,old$dispersion)
})

test_that("NB moment calibration respects bounds, weights and known shape", {
  skip_if_not_installed("mgcv")
  set.seed(722);n <- 10000L;X <- cbind(1,rnorm(n));mu <- exp(2+.2*X[,2])
  y <- rnbinom(n,mu=mu,size=4);family <- mgcv::nb(theta=-1)
  shape <- L0glm_calibration_theta(y,mu,rep(1,n),X,family,n-2L)
  expect_lt(abs(shape/4-1),.1)
  w <- rep(c(0,.5,1,2),length.out=n);df <- sum(w>0)-2L
  shape <- L0glm_calibration_theta(y,mu,w,X,family,df)
  expect_lt(abs(shape/4-1),.15)
  expect_equal(L0glm_calibration_theta(y,mu,w*7,X,family,df),shape,tolerance=1e-6)
  expect_equal(L0glm_calibration_theta(mu,mu,rep(1,n),X,family,n-2L,c(.01,100)),100)
  expect_equal(L0glm_calibration_theta(y,mu,rep(1,n),X,family,n-2L,c(.01,1)),1)
})

test_that("calibration runs independently inside CV training data", {
  set.seed(723);X <- matrix(rnorm(40*10),40,10)
  y <- outer(X[,1],1:3)+matrix(rnorm(120),40,3);seen <- list()
  local_mocked_bindings(L0glm_channel_calibration=function(base,...) {
    seen[[length(seen)+1L]] <<- dim(base$y)
    list(reason="consistent_scale")
  },.package="L0glm")
  result <- l0glmTune(X,y,family=gaussian(),tune="cv",nfolds=2L,foldid=rep(1:2,20),
    nlambda=3L,prop_range=.2,plot=FALSE,printbest=FALSE,parallelize=FALSE,
    local.search=FALSE,control.family=list(dispersion.calibration=TRUE))
  sizes <- do.call(rbind,seen)
  expect_true(all(c(20L,40L)%in%sizes[,1L]))
  expect_true(all(sizes[,2L]==3L))
  expect_true(is.finite(result$lambda_best))
  expect_null(result$fit_best$specification$start)
})

test_that("calibration controls are validated without forcing unused solver arguments", {
  X <- matrix(seq_len(30),10,3);y <- cbind(seq_len(10),seq_len(10)+1)
  for(mode in list(NA,1,"invalid",NULL))if(!is.null(mode))
    expect_error(l0glmFit(X,y,parallelize=FALSE,local.search=FALSE,
      control.family=list(dispersion.calibration=mode)),"dispersion.calibration")
  expect_no_error(l0glmFit(X,y[,1],parallelize=FALSE,ncores=stop("unused ncores"),
    local.search=FALSE,control.family=list(dispersion.calibration=FALSE)))
})

test_that("a calibrated scale stays frozen while remaining data-estimated", {
  set.seed(724);X <- matrix(rnorm(80*6),80,6)
  y <- outer(X[,1],1:3)+matrix(rnorm(240,sd=3),80,3)
  for(group in list(NULL,seq_len(6))) {
    fit <- l0glmFit(X,y,group=group,dispersion=-2,lambda="gic",relaxed=TRUE,
      control.solver=c(controlSolverGen("eigen",1L),list(.dispersion.fixed=TRUE)),
      parallelize=FALSE,local.search=FALSE)
    expect_identical(fit$dispersion,rep(2,3))
    expect_true(all(fit$dispersion.estimated))
    expect_equal(l0glmMetrics(fit)$df,sum(coef(fit)!=0)+3L)
  }
})

test_that("validation checkpoints do not overwrite the recorded initial scale", {
  set.seed(728);X <- matrix(rnorm(60*2),60,2);y <- 20*X[,1]+rnorm(60)
  validation <- L0glm_validation_context(list(X=X,y=y),X,y,gaussian(),FALSE)
  fit <- L0glm.fit.internal(X,y,family=gaussian(),dispersion=-2,start=c(1,1),
    maxit=1L,validation=validation)
  expect_identical(fit$initialDispersion,2)
})

test_that("NB calibration starts and frozen estimated shapes work for predictor groups", {
  skip_if_not_installed("mgcv")
  set.seed(726);X <- matrix(rnorm(60*4),60,4)
  y <- matrix(rnbinom(120,mu=4,size=3),60,2)
  control <- c(controlSolverGen("eigen",1L),list(.theta.start=c(2,3),.dispersion.fixed=TRUE))
  fit <- l0glmFit(X,y,family=mgcv::nb(theta=-1),group=seq_len(4),lambda="gic",
    maxit=1L,relaxed=FALSE,control.solver=control,parallelize=FALSE,local.search=FALSE)
  expect_equal(fit$theta,c(2,3));expect_true(all(fit$theta.estimated))
  expect_equal(fit$initialTheta,c(2,3))
  expected <- sum(vapply(1:2,function(j)sum(dnbinom(y[,j],mu=fitted(fit)[,j],size=c(2,3)[j],log=TRUE)),1))
  expect_equal(as.numeric(logLik(fit)),expected)
})

test_that("common-noise scoring respects dimension-based and specified group costs", {
  set.seed(727);X <- matrix(rnorm(60*4),60,4);y <- cbind(X[,1]+rnorm(60),X[,2]+rnorm(60))
  for(penalty in list("gic",2)) {
    fit <- l0glmFit(X,y,group=c(1,1,2,2),lambda=penalty,parallelize=FALSE,
      local.search=FALSE,control.family=list(dispersion.calibration=FALSE))
    fit$coefficients[1,1] <- 0
    common <- c(1,2);model <- fit;model$dispersion <- common;model$dispersion.estimated[] <- FALSE
    metrics <- l0glmMetrics(model,mu_floor=1e-8)
    expected <- if(is.character(penalty))-2*metrics$loglik+fit$init.lambda*metrics$df else metrics$group_objective
    expect_equal(L0glm_calibration_score(fit,common,FALSE),expected)
  }
})
