guard_pair <- function(args) {
  selected <- suppressWarnings(do.call(l0glmFit,c(args,list(relaxed=FALSE))))
  refitted <- suppressWarnings(do.call(l0glmFit,c(args,list(relaxed=TRUE))))
  expect_true(all(refitted$relaxationSkipped))
  expect_identical(coef(refitted),coef(selected))
  expect_equal(refitted$dispersion,selected$dispersion)
  expect_s4_class(coef(refitted),'dgCMatrix')
  refitted
}

test_that('ordinary and specialized drivers retain unsupported selected fits', {
  set.seed(9001)
  n <- 24L;p <- 36L;X <- matrix(runif(n*p,0,.1),n,p)
  counts <- matrix(rpois(2*n,10),n,2)
  positive <- matrix(rgamma(2*n,3),n,2)
  common <- list(intercept=FALSE,start=1,lambda=.01,normalize='none',
    grouped=TRUE,maxit=1L,group.maxit=1L,local.search=FALSE,
    parallelize=FALSE,nthreads=1L,screen='none')
  for(sparse in c(FALSE,TRUE)) {
    common$X <- if(sparse)Matrix::Matrix(X,sparse=TRUE)else X
    for(f in list(gaussian(),poisson(identity),poisson(log),binomial(),
      Gamma(log),inverse.gaussian(log),quasipoisson(log),quasibinomial())) {
      y <- if(f$family=='gaussian')matrix(rnorm(2*n),n,2)else
        if(f$family%in%c('binomial','quasibinomial'))matrix(rbinom(2*n,1,.4),n,2)else
        if(f$family%in%c('Gamma','inverse.gaussian'))positive else counts
      args <- c(common,list(y=y,family=f))
      if(f$link=='identity'&&f$family=='poisson')args$offset <- rep(10,n)
      fit <- guard_pair(args)
      expect_equal(dim(coef(fit)),c(p,2L))
    }
    A <- matrix(.1,n,p)
    common$X <- if(sparse)Matrix::Matrix(A,sparse=TRUE)else A
    guard_pair(c(common,list(y=counts,family=poisson(identity),nonnegative=TRUE)))
  }
})

test_that('special-family refits check their joint information space', {
  skip_if_not_installed('mgcv');skip_if_not_installed('survival')
  set.seed(9004);n <- 24L;p <- 36L
  common <- list(intercept=FALSE,lambda=.01,normalize='none',start=1,
    maxit=1L,grouped=TRUE,group.maxit=1L,local.search=FALSE,parallelize=FALSE,
    screen='none',nthreads=1L)
  for(sparse in c(FALSE,TRUE)) {
    X <- matrix(rnorm(n*p),n,p)
    common$X <- if(sparse)Matrix::Matrix(X,sparse=TRUE)else X
    guard_pair(c(common,list(y=rpois(n,5),family=mgcv::nb(theta=2,link='log'))))
    guard_pair(c(common,list(y=rpois(n,5),family=mgcv::nb(theta=-2,link='log'))))
    args <- c(common,list(y=factor(rep(0:2,length.out=n)),family=mgcv::multinom(K=2)))
    args$start <- matrix(rnorm(p*3,sd=.1),p,3)
    fit <- guard_pair(args)
    expect_equal(as.numeric(Matrix::rowSums(coef(fit))),rep(0,p),tolerance=1e-12)
    guard_pair(c(common,list(y=survival::Surv(rexp(n),rep(c(1,1,0),length.out=n)),
      family=mgcv::cox.ph())))
  }
})

test_that('rank deficiency is checked even when selected size is below n', {
  set.seed(9002);A <- cbind(rnorm(60),rnorm(60));A <- cbind(A,A[,1])
  common <- list(X=A,intercept=TRUE,start=1,lambda=.01,normalize='rms',maxit=1L,
    local.search=FALSE,parallelize=FALSE,screen='none',nthreads=1L)
  for(f in list(gaussian(),poisson(log),binomial(),Gamma(log))) {
    y <- if(f$family=='gaussian')rnorm(60)else if(f$family=='binomial')rbinom(60,1,.5)else
      if(f$family=='Gamma')rgamma(60,3)else rpois(60,5)
    guard_pair(c(common,list(y=y,family=f)))
  }
})

test_that('channel masks preserve supported refits and original starting exclusions', {
  set.seed(9005);n <- 40L;p <- 12L;X <- matrix(rnorm(n*p),n,p)
  Y <- matrix(rnorm(2*n),n,2);weights <- matrix(1,n,2);weights[9:n,2] <- 0
  start <- matrix(1,p,2);start[p,] <- 0
  args <- list(X=X,y=Y,weights=weights,dispersion=1,start=start,lambda=.01,
    normalize='none',maxit=1L,group.maxit=1L,local.search=FALSE,
    parallelize=FALSE,screen='none',nthreads=1L)
  selected <- do.call(l0glmFit,c(args,list(relaxed=FALSE)))
  fit <- do.call(l0glmFit,c(args,list(relaxed=TRUE)))
  expect_identical(fit$relaxationSkipped,c(FALSE,TRUE))
  expect_identical(coef(fit)[,2,drop=FALSE],coef(selected)[,2,drop=FALSE])
  expect_equal(as.numeric(coef(fit)[p,]),c(0,0))
  expected <- solve(crossprod(X[,-p]),crossprod(X[,-p],Y[,1]))
  expect_equal(as.numeric(coef(fit)[-p,1]),as.numeric(expected),tolerance=1e-9)
})

test_that('independent channels refit their own selected supports', {
  set.seed(9006);n <- 30L;p <- 42L;X <- matrix(rnorm(n*p),n,p)
  start <- matrix(0,p,2);start[1:21,1] <- 1;start[22:42,2] <- 1
  Y <- matrix(rnorm(2*n),n,2)
  fit <- l0glmFit(X,Y,start=start,lambda=.01,dispersion=1,grouped=FALSE,
    relaxed=TRUE,maxit=1L,normalize='none',local.search=FALSE,
    parallelize=FALSE,screen='none')
  expect_false(any(fit$relaxationSkipped))
  for(k in 1:2) {
    use <- which(start[,k]!=0)
    expected <- solve(crossprod(X[,use]),crossprod(X[,use],Y[,k]))
    expect_equal(as.numeric(coef(fit)[use,k]),as.numeric(expected),tolerance=1e-8)
    expect_equal(as.numeric(coef(fit)[-use,k]),rep(0,p-length(use)))
  }
})

test_that('independent native Poisson channels do not pool their refit support', {
  X <- cbind(diag(30),diag(30));Y <- matrix(0,30,2);Y[1:20,] <- 10
  start <- matrix(0,60,2);start[1:20,1] <- 1;start[31:50,2] <- 1
  fit <- l0nnPois(X,Y,start=start,lambda=.01,grouped=FALSE,relaxed=TRUE,
    maxit=3L,normalize='none',local.search=FALSE,parallelize=FALSE)
  expect_false(any(fit$relaxationSkipped))
  expect_equal(as.matrix(coef(fit))[start!=0],rep(10,40),tolerance=1e-9)
  expect_equal(as.matrix(coef(fit))[start==0],rep(0,80))
})

test_that('predictor groups and the R fallback share the refit guard', {
  set.seed(9007);n <- 24L;p <- 36L;X <- matrix(rnorm(n*p),n,p)
  common <- list(X=X,y=rnorm(n),start=1,lambda=.01,maxit=1L,
    normalize='none',local.search=FALSE,parallelize=FALSE,screen='none')
  guard_pair(c(common,list(group=rep(1:12,each=3))))
  testthat::local_mocked_bindings(L0glm_native_family=function(...)NULL,.package='L0glm')
  guard_pair(common)
  common$y <- matrix(rnorm(n*2),n,2);common$group.maxit <- 1L
  guard_pair(common)
})

test_that('identifiable Gaussian refits retain ordinary least-squares accuracy', {
  set.seed(9003);X <- matrix(rnorm(150*8),150,8)
  Y <- X[,1:3]%*%matrix(rnorm(6),3,2)+matrix(rnorm(300),150,2)
  fit <- l0glmFit(X,Y,lambda='gic',intercept=TRUE,relaxed=TRUE,
    parallelize=FALSE,nthreads=1L,local.search=FALSE)
  active <- which(Matrix::rowSums(abs(coef(fit)))>0)
  A <- cbind(1,X)[,active,drop=FALSE]
  expect_false(any(fit$relaxationSkipped))
  expect_equal(as.matrix(coef(fit))[active,,drop=FALSE],
    solve(crossprod(A),crossprod(A,Y)),tolerance=1e-9,ignore_attr=TRUE)
})

test_that('fixed scalar signal families share the guard', {
  set.seed(9008);X <- matrix(runif(24*36,0,.1),24,36);Y <- matrix(rpois(48,10),24,2)
  for(f in list(nb1(3,log),gammaScale(2,log,censor=.5),gammaCensored(.3,log)))
    guard_pair(list(X=X,y=Y,family=f,start=1,lambda=.01,maxit=1L,group.maxit=1L,
      grouped=TRUE,local.search=FALSE,parallelize=FALSE,normalize='none',screen='none'))
})

test_that('an estimated scale requires residual degrees of freedom', {
  set.seed(9009);X <- matrix(rnorm(36),6,6);y <- rnorm(6)
  guard_pair(list(X=X,y=y,start=1,lambda=.001,maxit=1L,normalize='none',
    local.search=FALSE,parallelize=FALSE,screen='none'))
  X <- matrix(rnorm(24*23),24,23)
  expect_false(L0glm_relaxed_standard(X,NULL,matrix(TRUE,23,1),TRUE,nuisance_df=1L))
})

test_that('external native solvers preserve skipped channels', {
  skip_if_not_installed('glmnet')
  set.seed(9010);X <- matrix(rnorm(40*12),40,12);Y <- matrix(rnorm(80),40,2)
  weights <- matrix(1,40,2);weights[9:40,2] <- 0
  args <- list(X=X,y=Y,start=1,weights=weights,dispersion=1,lambda=.001,
    maxit=1L,group.maxit=1L,normalize='none',local.search=FALSE,
    solver='bigGlm',parallelize=FALSE,screen='none')
  selected <- suppressWarnings(do.call(l0glmFit,c(args,list(relaxed=FALSE))))
  fit <- suppressWarnings(do.call(l0glmFit,c(args,list(relaxed=TRUE))))
  expect_identical(fit$relaxationSkipped,c(FALSE,TRUE))
  expect_identical(coef(fit)[,2,drop=FALSE],coef(selected)[,2,drop=FALSE])
  expect_equal(summary(fit)$convergence$skipped,c(0L,1L))
  args$parallelize <- TRUE;args$ncores <- 2L;args$nthreads <- 1L
  parallel_fit <- suppressWarnings(do.call(l0glmFit,c(args,list(relaxed=TRUE))))
  expect_identical(parallel_fit$relaxationSkipped,fit$relaxationSkipped)
  expect_equal(coef(parallel_fit),coef(fit),tolerance=1e-12)
})

test_that('entrywise CLR information can be full rank with more rows than observations', {
  set.seed(9011)
  allowed <- rbind(c(TRUE,TRUE,FALSE),c(TRUE,FALSE,TRUE),
    c(FALSE,TRUE,TRUE),c(TRUE,TRUE,FALSE))
  data <- list(kind='multinomial',m=3L,Z=matrix(rnorm(12),3,4),
    y=matrix(1,3,3),weights=matrix(1,3,1))
  expect_true(L0glm_relaxed_check(data,matrix(0,4,3),allowed))
  data$Z[,4] <- data$Z[,1]
  expect_false(L0glm_relaxed_check(data,matrix(0,4,3),allowed))
})
