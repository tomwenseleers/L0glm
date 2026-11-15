test_that('workspace controls reject unsafe values and guard direct allocations', {
  expect_error(controlSolverGen('eigen',1L,max.memory=-1),'max.memory')
  expect_error(controlSolverGen('eigen',1L,max.memory=NA_real_),'max.memory')
  expect_null(controlSolverGen('eigen',1L)$max.memory)
  set.seed(91);X <- matrix(rnorm(120*20),120,20);y <- rnorm(120)
  for(native in c(TRUE,FALSE)) {
    ctl <- controlSolverGen('eigen',1L,max.memory=0);ctl$.native <- native
    expect_error(l0glmFit(X,y,dispersion=1,intercept=FALSE,local.search=FALSE,
      start=1,maxit=1,control.solver=ctl),'workspace budget')
  }
})

test_that('matrix-free ridge matches LLT for primal and Woodbury systems', {
  set.seed(551)
  for(wide in c(FALSE,TRUE))for(sparse in c(FALSE,TRUE)) {
    n <- if(wide)35L else 90L;p <- if(wide)100L else 18L
    X <- matrix(rnorm(n*p),n,p);X[,2] <- .95*X[,1]+.05*X[,2]
    X[,1] <- 1
    if(sparse)X <- Matrix::Matrix(X*(abs(X)>.4),sparse=TRUE)
    y <- rnorm(n);d <- seq(.3,1.5,length.out=p);d[1] <- 0
    center <- rnorm(p)/10;start <- numeric(p)
    reference <- lmridge.eigen(X,y,d,start=start,prior.mean=center)
    for(method in 3:4) {
      ctl <- controlSolverGen('eigen',method,max.memory=0,tol=1e-10,maxit=2000L)
      answer <- lmridge.eigen(X,y,d,start=start,prior.mean=center,control.solver=ctl)
      expect_equal(answer,reference,tolerance=2e-8)
      gradient <- as.numeric(crossprod(X,X%*%answer-y))+d*(answer-center)
      rhs <- as.numeric(crossprod(X,y))+d*center
      expect_lte(sqrt(sum(gradient^2))/sqrt(sum(rhs^2)),1.1e-10)
    }
  }
})

test_that('matrix-free native fits and dispersion pilots retain the fitted solution', {
  set.seed(92);X <- matrix(rnorm(150*25),150,25)
  y <- X[,1]*2-X[,4]+rnorm(150)
  reference <- l0glmFit(X,y,local.search=FALSE,ncores=1)
  for(native in c(TRUE,FALSE)) {
    ctl <- controlSolverGen('eigen',3L,max.memory=0,tol=1e-10,maxit=2000L);ctl$.native <- native
    answer <- l0glmFit(X,y,local.search=FALSE,ncores=1,method=3L,control.solver=ctl)
    expect_identical(as.numeric(coef(answer))!=0,as.numeric(coef(reference))!=0)
    expect_equal(as.numeric(coef(answer)),as.numeric(coef(reference)),tolerance=1e-6)
    expect_equal(answer$dispersion,reference$dispersion,tolerance=1e-6)
  }
})
test_that('the adaptive Gaussian pilot respects the workspace budget', {
  set.seed(825);n <- 80L;p <- 180L
  data <- list(Z=matrix(rnorm(n*p),n,p),y=matrix(rnorm(n*2),n,2),
    weights=matrix(1,n,1),offset=matrix(0,n,1),included=integer())
  probes <- matrix(sample(c(-1,1),n*4,TRUE),n,4)
  reference <- L0glm_gaussian_pilot_scale(data,rep(.5,p),probes,1)
  actual <- L0glm_gaussian_pilot_scale(data,rep(.5,p),probes,1,
    controlSolverGen('eigen',3L,max.memory=0))
  expect_true(all(is.finite(actual)))
  expect_equal(actual,reference,tolerance=1e-7)
})

test_that('serial special-family adapters guard information before allocation', {
  set.seed(826);X <- matrix(rnorm(80*8),80,8)
  ctl <- controlSolverGen('eigen',1L,max.memory=0);ctl$.native <- FALSE
  args <- list(solver='eigen',method=1L,control.solver=ctl)
  data <- L0glm_ext_prepare(X,cbind(rexp(80),rep(c(1,0),40)),mgcv::cox.ph(),
    NULL,NULL,FALSE,NULL,FALSE,'none',1,NULL,list(cox='direct'))
  expect_error(with_mocked_bindings(L0glm_ext_cox_working(data,data$B,matrix(1,8,1),list(),args),
    L0glm_ext_cox_information=function(...)stop('information allocated'),.package='L0glm'),'workspace budget')
  data <- L0glm_ext_prepare(X,rep(0:2,length.out=80),mgcv::multinom(2),
    NULL,NULL,FALSE,NULL,FALSE,'none',NULL,NULL,list())
  expect_error(with_mocked_bindings(L0glm_ext_multinomial_working(data,data$B,matrix(1,8,3),list(),args),
    L0glm_ext_quadratic=function(...)stop('information allocated'),.package='L0glm'),'workspace budget')
})
