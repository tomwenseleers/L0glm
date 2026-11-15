test_that('NB coordination uses native steps and the requested OpenMP channel budget', {
  set.seed(806);X <- matrix(rnorm(120*12),120,12)
  y <- matrix(rnbinom(120*4,mu=3,size=4),120,4)
  original <- glm_extended_cpp;used <- integer()
  wrapped <- function(...,nthreads) {used <<- c(used,nthreads);original(...,nthreads=nthreads)}
  fit <- with_mocked_bindings(l0glmFit(X,y,family=mgcv::nb(),start=1,lambda=1,
    maxit=4L,local.search=FALSE,ncores=4L,control.family=list(dispersion.calibration=FALSE)),
    glm_extended_cpp=wrapped,
    L0glm_ext_nb_working=function(...)stop('R NB working loop should not run'),.package='L0glm')
  expect_s3_class(fit,'L0glm_fit')
  expect_true(length(used)>0L)
  expect_true(all(used==4L))
})

test_that('native NB preserves weights, priors, start masks, links and grouped selection', {
  set.seed(807);n <- 140L;p <- 18L;m <- 3L
  X <- matrix(runif(n*p),n,p);beta <- matrix(0,p,m);beta[c(2,4,9),] <- .7
  w <- runif(n,.5,1.5);w[1:5] <- 0
  start <- matrix(1,p,m);start[6,] <- 0;start[8,2] <- 0
  prior <- matrix(.03,p,m)
  for(link in c('identity','log','sqrt'))for(sparse in c(FALSE,TRUE))for(grouped in c(FALSE,TRUE)) {
    family <- do.call(mgcv::nb,list(link=link))
    mu <- family$linkinv(X%*%beta+2)
    y <- matrix(rnbinom(n*m,mu=mu,size=4),n,m)
    design <- if(sparse)as(Matrix::Matrix(X,sparse=TRUE),'dgCMatrix')else X
    args <- list(X=design,y=y,family=family,lambda=1,weights=w,offset=rep(2,n),
      start=start,prior.mean=prior,nonnegative=TRUE,always.include=1L,
      grouped=grouped,intercept=FALSE,local.search=FALSE,ncores=2L,maxit=25L,
      control.family=list(dispersion.calibration=FALSE))
    ctl <- controlSolverGen('eigen',1L);ctl$.native <- FALSE
    reference <- do.call(l0glmFit,c(args,list(control.solver=ctl)))
    actual <- do.call(l0glmFit,args)
    expect_equal(as.matrix(coef(actual)),as.matrix(coef(reference)),tolerance=2e-6)
    expect_identical(as.numeric(coef(actual))!=0,as.numeric(coef(reference))!=0)
    expect_equal(actual$theta,reference$theta,tolerance=2e-6)
    expect_true(all(coef(actual)[6,]==0))
    expect_equal(as.numeric(coef(actual)[8,2]),0)
  }
})

test_that('joint multinomial coordination stays native with CLR masks and wide systems', {
  set.seed(820)
  for(wide in c(FALSE,TRUE))for(grouped in c(FALSE,TRUE)) {
    n <- if(wide)30L else 100L;p <- if(wide)55L else 12L;C <- 3L
    X <- matrix(rnorm(n*p),n,p);Y <- matrix(rmultinom(n,5,c(.3,.3,.4)),n,C,byrow=TRUE)
    start <- matrix(rnorm(p*C,0,.1),p,C);start[2,] <- 0
    args <- list(X=X,y=Y,family=mgcv::multinom(2),start=start,lambda=1,intercept=FALSE,
      maxit=20L,local.search=FALSE,grouped=grouped,ncores=3L)
    ctl <- controlSolverGen('eigen',1L);ctl$.native <- FALSE
    reference <- do.call(l0glmFit,c(args,list(control.solver=ctl)))
    actual <- with_mocked_bindings(do.call(l0glmFit,args),
      L0glm_ext_multinomial_working=function(...)stop('R joint working loop should not run'),.package='L0glm')
    expect_equal(as.matrix(coef(actual)),as.matrix(coef(reference)),tolerance=1e-7)
    expect_identical(as.numeric(coef(actual))!=0,as.numeric(coef(reference))!=0)
    expect_lt(max(abs(rowSums(coef(actual)))),1e-10)
    expect_true(all(coef(actual)[2,]==0))
  }
})
