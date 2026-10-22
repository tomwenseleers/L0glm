test_that("screened fits match standalone SIS/ISIS with usable original dimensions", {
  set.seed(2610); n <- 150; X <- matrix(rnorm(n*120),n,120)
  colnames(X) <- paste0("x",1:120)
  y <- 1+X[,1]*2-X[,2]+rnorm(n,sd=.5)
  w <- sample(1:3,n,TRUE); off <- runif(n,-.2,.2)
  for (A in list(X,Matrix::Matrix(X,sparse=FALSE),Matrix::Matrix(X,sparse=TRUE))) {
    for (method in c("SIS","ISIS")) {
      rounds <- if(method=="SIS") 1L else 3L
      actual <- L0glm.fit(A,y,weights=w,offset=off,screen=method,nsis=25,
        screen.maxit=3,intercept=TRUE,always.include=3)
      expected <- sis_glm(A,y,weights=w,offset=off,nsis=25,iterative=TRUE,niter=rounds,
        return_stats=FALSE,intercept=TRUE,always.include=3)
      expect_equal(coef(actual),coef(expected$fit),tolerance=1e-10)
      expect_equal(actual$screening$history,expected$history,tolerance=1e-10)
      expect_identical(actual$screening$method,method)
      expect_false(actual$screening$automatic)
      expect_identical(dim(coef(actual)),c(121L,1L));expect_s4_class(coef(actual),"sparseMatrix")
      expect_equal(GIC(actual),min(actual$screening$history$gic),tolerance=1e-8)
      expect_equal(predict(actual,A,offset=off),fitted(actual),tolerance=1e-10)
      expect_equal(predict(actual,A[,120:1],offset=off),fitted(actual),tolerance=1e-10)
      expect_equal(inherits(actual$X,"sparseMatrix"),inherits(A,"sparseMatrix"))
      expect_output(print(actual),paste0("Screening: ",method))
      expect_identical(summary(actual)$screening,actual$screening)
    }
  }
})

test_that("screened defaults reproduce the small-batch GIC relaxed workflow", {
  set.seed(2611);n <- 80L;p <- 1000L;X <- matrix(rnorm(n*p),n,p)
  y <- X[,1]*2-X[,2]*2+rnorm(n)
  auto <- L0glm.fit(X,y)
  expect_identical(auto$screening$method,"ISIS");expect_true(auto$screening$automatic)
  expect_equal(auto$screening$screen_size,max(1L,round(n/(log(p)*log(log(n))))))
  expect_equal(auto$init.lambda,preset.lambda("gic",p,n))
  explicit <- L0glm.fit(X,y,screen="ISIS")
  expect_equal(coef(auto),coef(explicit),tolerance=1e-10)
  expect_equal(auto$screening$history,explicit$screening$history,tolerance=1e-10)
  expect_false(explicit$screening$automatic)
  expect_true(any(auto$convergence$stage=="relaxed"))
  sis <- L0glm.fit(X,y,screen="SIS")
  expect_equal(sis$screening$screen_size,2L*floor(n/log(n)))
  expect_equal(sis$screening$iterations,1L)
})

test_that("auto uses eligible rows and does not screen ordinary designs or masked small problems", {
  set.seed(2612);X <- matrix(rnorm(80*15),80,15);y <- rnorm(80)
  expect_identical(L0glm.fit(X,y),L0glm.fit(X,y,screen="none"))
  expect_identical(L0glm_screen_method("auto",Matrix::Matrix(0,20,9999,sparse=TRUE),NULL,FALSE,NULL),"ISIS")
  expect_identical(L0glm_screen_method("auto",Matrix::Matrix(0,20,200,sparse=TRUE),NULL,FALSE,NULL),"none")
  expect_identical(L0glm_screen_method("auto",Matrix::Matrix(0,5,51,sparse=TRUE),NULL,FALSE,NULL),"ISIS")
  A <- Matrix::Matrix(0,20,201,sparse=TRUE)
  expect_identical(L0glm_screen_method("auto",A,NULL,FALSE,NULL),"ISIS")
  expect_identical(L0glm_screen_method("auto",A,NULL,FALSE,1),"none")
  expect_identical(L0glm_screen_method("auto",A,0,FALSE,NULL),"none")
  start <- matrix(1,201,2);start[1,] <- 0
  expect_identical(L0glm_screen_method("auto",A,start,FALSE,NULL),"none")
  start[,2] <- -1;start[1,] <- c(1,-1)
  expect_identical(L0glm_screen_method("auto",A,start,FALSE,NULL),"ISIS")
  expect_identical(L0glm_screen_method("auto",A,Matrix::Matrix(start,sparse=TRUE),FALSE,NULL),"ISIS")
  expect_identical(L0glm_screen_method("auto",A,NULL,FALSE,1:200),"none")
  expect_identical(L0glm_screen_method("auto",Matrix::Matrix(0,1000,10000,sparse=TRUE),NULL,FALSE,NULL),"none")
  expect_identical(L0glm_screen_method("auto",Matrix::Matrix(0,1001,10000,sparse=TRUE),NULL,FALSE,NULL),"none")
  expect_identical(L0glm_screen_method("auto",Matrix::Matrix(0,10000,10000,sparse=TRUE),NULL,FALSE,NULL),"none")
  expect_identical(L0glm_screen_method("auto",Matrix::Matrix(0,500,1000000,sparse=TRUE),NULL,FALSE,NULL),"ISIS")
})

test_that("explicit fitting controls and all starting-zero exclusions survive screening", {
  set.seed(2613);n <- 110;X <- matrix(rnorm(n*40),n,40);y <- 1+X[,1]*2+rnorm(n)
  Y <- cbind(y,-y);start <- matrix(1,41,2);start[5,] <- 0;start[6,1] <- 0
  start[1,] <- 0;start[8,] <- 0
  prior <- matrix(0,41,2);prior[5,] <- 100
  fit <- L0glm.fit(X,Y,screen="ISIS",nsis=15,screen.maxit=2,intercept=TRUE,
    start=start,prior.mean=prior,always.include=7,grouped=FALSE,parallelize=FALSE,
    lambda=3,relaxed=FALSE,normalize="none",maxit=5,tol=1e-6)
  expected <- sis_glm(X,Y,nsis=15,iterative=TRUE,niter=2,return_stats=FALSE,intercept=TRUE,
    start=start,prior.mean=prior,always.include=7,grouped=FALSE,parallelize=FALSE,
    lambda=3,relaxed=FALSE,normalize="none",maxit=5,tol=1e-6)
  expect_equal(coef(fit),coef(expected$fit),tolerance=1e-10)
  expect_true(all(coef(fit)[5,]==0));expect_equal(coef(fit)[6,1],0)
  expect_false(4L%in%fit$screening$candidates);expect_true(7L%in%fit$screening$candidates)
  expect_equal(fit$always.include,c(1L,8L));expect_equal(fit$init.lambda,3)
  expect_false(any(fit$convergence$stage=="relaxed"));expect_equal(fit$preprocessing$mode,"none")
  zero <- L0glm.fit(X,y,screen="ISIS",start=0)
  expect_true(all(coef(zero)==0));expect_length(zero$screening$indices,0)
  expect_equal(zero$screening$iterations,1L)
})

test_that("screened fits support nonnegative grouped Poisson and all special families", {
  skip_if_not_installed("mgcv");skip_if_not_installed("survival")
  set.seed(2614);n <- 120;X <- Matrix::rsparsematrix(n,15,.2);P <- abs(X)
  Y <- matrix(rpois(n*2,2+as.numeric(P[,1])*2),n,2)
  for (openmp in c(TRUE,FALSE)) {
    fit <- L0glm.fit(P,Y,family=poisson(identity),nonnegative=TRUE,openmp=openmp,
      screen="ISIS",nsis=6,screen.maxit=2,parallelize=FALSE)
    expected <- sis_glm(P,Y,family=poisson(identity),nonnegative=TRUE,openmp=openmp,
      nsis=6,niter=2,iterative=TRUE,return_stats=FALSE)
    expect_equal(coef(fit),coef(expected$fit),tolerance=1e-8)
    expect_true(all(coef(fit)>=0));expect_identical(dim(coef(fit)),c(15L,2L))
  }
  off <- runif(n,-.1,.1);w <- sample(1:3,n,TRUE)
  counts <- rnbinom(n,mu=exp(.4+off),size=2)
  nb <- L0glm.fit(X,counts,family=mgcv::nb(),weights=w,offset=off,
    intercept=TRUE,screen="SIS",nsis=6)
  expect_true(all(is.finite(nb$theta)&nb$theta>0));expect_identical(dim(coef(nb)),c(16L,1L))
  labels <- factor(sample(letters[1:3],n,TRUE));start <- matrix(1,16,3);start[5,] <- 0
  for (grouped in c(TRUE,FALSE)) {
    mn <- L0glm.fit(X,labels,family=mgcv::multinom(2),screen="ISIS",nsis=6,
      screen.maxit=2,start=start,grouped=grouped)
    expect_identical(dim(coef(mn)),c(16L,3L));expect_true(all(coef(mn)[5,]==0))
    expect_equal(as.numeric(Matrix::rowSums(coef(mn))),rep(0,16),tolerance=1e-10)
    expect_equal(rowSums(predict(mn,X)),rep(1,n),tolerance=1e-10)
  }
  time <- rexp(n,exp(.5*as.numeric(X[,1])))
  cox <- L0glm.fit(X,survival::Surv(pmin(time,2),time<=2),family=mgcv::cox.ph(),
    screen="ISIS",nsis=6,screen.maxit=2,weights=w)
  expect_identical(dim(coef(cox)),c(15L,1L));expect_true(all(is.finite(predict(cox,X))))
  quasi <- L0glm.fit(X,rpois(n,2),family=quasipoisson(),intercept=TRUE,
    screen="ISIS",nsis=6,screen.maxit=2)
  expect_true(all(is.na(quasi$screening$history$gic)))
  expect_equal(quasi$screening$best_iteration,quasi$screening$iterations)
})

test_that("cross-validation screens its training rows rather than the full response", {
  set.seed(2615);X <- matrix(rnorm(90*12),90,12);y <- X[,1]+rnorm(90)
  sizes <- integer();original <- L0glm_screen_fit
  with_mocked_bindings({
    fit <- L0glm.tune(X,y,family=gaussian(),tune="cv",nlambda=3,
      foldid=rep(1:3,30),screen="ISIS",nsis=5,screen.maxit=1,plot=FALSE,printbest=FALSE)
    expect_true(is.finite(fit$lambda_best))
  },L0glm_screen_fit=function(X,...) {sizes <<- c(sizes,nrow(X));original(X,...)})
  expect_true(all(c(60L,90L)%in%sizes));expect_true(all(sizes%in%c(60L,90L)))
})

test_that("screen controls validate clearly and internal ridge pilots disable screening", {
  X <- matrix(1:150,30,5);y <- rnorm(30)
  for (choice in list(TRUE,1,NA,character(),c("SIS","none"),"invalid"))
    expect_error(L0glm.fit(X,y,screen=choice),"screen must")
  for (bad in list(0,-1,1.2,NA,Inf,"3",c(1,2))) {
    expect_error(L0glm.fit(X,y,screen="ISIS",nsis=bad),"nsis must")
    expect_error(L0glm.fit(X,y,screen="ISIS",screen.maxit=bad),"screen.maxit must")
  }
  expect_identical(L0glm_screen_method("sis",X,NULL,FALSE,NULL),"SIS")
  expect_identical(L0glm_screen_method("NONE",X,NULL,FALSE,NULL),"none")
  set.seed(2616);X <- matrix(rnorm(70*12),70,12);y <- X[,1]+rnorm(70)
  original <- L0glm.fit;calls <- list()
  with_mocked_bindings({
    fit <- sparseglm(X,y,method="adaptive_lasso",nlambda=3,
      options=list(L0glm=list(screen="ISIS",nsis=2)))
    expect_identical(fit$timings$status,"ok")
  },L0glm.fit=function(...) {calls[[length(calls)+1L]] <<- list(...);original(...)})
  expect_identical(calls[[1L]]$screen,"none")
  expect_equal(calls[[1L]]$maxit,1)
})
