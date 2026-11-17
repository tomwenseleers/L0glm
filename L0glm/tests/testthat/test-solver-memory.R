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

test_that('sparse workspace refinement counts unique overlap and Cholesky fill', {
  n <- 600L
  X <- Matrix::bandSparse(n,k=-12:12,diagonals=lapply(-12:12,function(k)rep(1,n-abs(k))))
  X <- as(X,'dgCMatrix');budget <- .01*1024^3
  coarse <- solver_workspace_estimate_cpp(X,TRUE,FALSE)
  refined <- solver_sparse_workspace_cpp(X,FALSE,budget)
  expect_gt(coarse,budget)
  expect_lt(refined,budget)
  expect_true(L0glm_workspace_fits(X,controlSolverGen('eigen',1L,max.memory=.01)))
  expect_true(is.infinite(solver_sparse_workspace_cpp(X,FALSE,1024)))
  # One centre connected to every leaf has little Gram storage but dense
  # fill under the supported natural ordering; the guard must count that fill.
  star <- Matrix::sparseMatrix(i=rep(seq_len(299L),2L),
    j=c(rep(1L,299L),2:300),x=1,dims=c(299L,300L))
  expect_true(is.infinite(solver_sparse_workspace_cpp(star,FALSE,.5*1024^2)))
  Z <- X;Z@x[seq(1L,length(Z@x),by=4L)] <- 0
  expect_equal(solver_sparse_workspace_cpp(Z,FALSE,budget),
    solver_sparse_workspace_cpp(Matrix::drop0(Z),FALSE,budget))
})

test_that('the Poisson guard accounts for profiled nuisance columns', {
  n <- 600L
  X <- as(Matrix::bandSparse(n,k=-12:12,
    diagonals=lapply(-12:12,function(k)rep(1,n-abs(k)))),'dgCMatrix')
  B <- Matrix::Matrix(splines::bs(seq_len(n),df=4L,intercept=TRUE),sparse=TRUE)
  A <- as(cbind(X,B),'dgCMatrix');ctl <- controlSolverGen('eigen',1L,max.memory=.01)
  expect_error(L0glm_require_workspace(A,ctl,rows=TRUE),'workspace budget')
  expect_no_error(L0glm_require_workspace(A,ctl,rows=TRUE,
    unpenalized=ncol(X)+seq_len(ncol(B))))
})

test_that('a corrected concurrent Poisson guard retains the fitted coefficients', {
  set.seed(827);n <- 600L
  X <- as(Matrix::bandSparse(n,k=-12:12,
    diagonals=lapply(-12:12,function(k)rep(exp(-k*k/30),n-abs(k)))),'dgCMatrix')
  beta <- matrix(0,n,4L);beta[c(130L,280L,460L),] <- 30
  y <- matrix(rpois(n*4,as.numeric(X%*%beta)+.2),n,4L)
  args <- list(X=X,y=y,family=poisson('identity'),nonnegative=TRUE,
    grouped=TRUE,intercept=FALSE,offset=.2,lambda='gic',local.search=FALSE)
  reference <- do.call(l0glmFit,c(args,list(ncores=1L)))
  actual <- do.call(l0glmFit,c(args,list(ncores=4L,
    control.solver=controlSolverGen('eigen',1L,max.memory=.04))))
  expect_equal(as.matrix(coef(actual)),as.matrix(coef(reference)),tolerance=1e-10)
  expect_equal(actual$deviance,reference$deviance,tolerance=1e-10)
})
