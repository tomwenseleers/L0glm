test_that('cached Gaussian paths retain full ISIS quality and original coordinates', {
  skip_if_not_installed('bigmemory')
  set.seed(811);n <- 160L;p <- 2600L
  A <- matrix(rnorm(n*p),n,p);for(j in 2:p)A[,j] <- .75*A[,j-1L]+sqrt(1-.75^2)*A[,j]
  colnames(A) <- paste0('x',1:p)
  w <- sample(0:2,n,TRUE);off <- runif(n,-.5,.5)
  y <- 2+3*A[,9]-2*A[,10]+A[,300]+off+rnorm(n,sd=.5)
  folder <- tempfile();dir.create(folder)
  M <- bigmemory::filebacked.big.matrix(n,p,type='double',backingfile='design.bin',
    descriptorfile='design.desc',backingpath=folder,dimnames=dimnames(A));M[,] <- A
  X <- l0glmData(M,control=list(max.cells=1e5,working.columns=500L))
  options <- list(y=y,weights=w,offset=off,intercept=TRUE,always.include=18L,
    screen='ISIS',nsis=12L,screen.maxit=6L,local.search=FALSE)
  for(dispersion in list(NULL,1,-1))for(start in list(NULL,1)) {
    args <- c(options,list(dispersion=dispersion,start=start))
    baseline <- X;baseline$control$working.columns <- 0L
    original <- do.call(l0glmFit,c(list(X=baseline),args))
    result <- do.call(l0glmFit,c(list(X=X),args))
    expect_type(result$screening$working.set,'list')
    expect_lte(GIC(result),GIC(original)+1e-7)
    expect_identical(dim(coef(result)),c(p+1L,1L))
    expect_s4_class(coef(result),'sparseMatrix')
    expect_length(result$preprocessing$scale,p+1L)
    expect_equal(result$p,p-1L)
    expect_equal(GIC(result),min(result$screening$history$gic),tolerance=1e-9)
    expect_equal(result$screening$best_iteration,which.min(result$screening$history$gic))
    expect_true(18L%in%result$screening$indices)
    expect_equal(predict(result,X,offset=off),fitted(result),tolerance=1e-10)
    expect_equal(predict(result,A[,p:1],offset=off),fitted(result),tolerance=1e-10)
    expect_lte(result$screening$working.set$global.scans,6L)
    expect_lte(result$screening$working.set$columns,floor(1e5/n))
    expect_lte(length(result$screening$candidates)*n,1e5)
  }
  automatic <- l0glmFit(X,y,weights=w,offset=off,intercept=TRUE,always.include=18L,
    screen='ISIS',nsis=12L,screen.maxit=6L)
  expect_true(is.finite(GIC(automatic)))
  expect_gte(automatic$screening$nthreads,1L)
  expect_true(isTRUE(automatic$local.search$applied))
  expect_equal(predict(automatic,X,offset=off),fitted(automatic),tolerance=1e-9)
})

test_that('working pools preserve unsupported-option and starting-mask paths', {
  skip_if_not_installed('bigmemory')
  set.seed(812);A <- matrix(rnorm(80*1200),80,1200);y <- A[,1]*2+rnorm(80)
  folder <- tempfile();dir.create(folder)
  M <- bigmemory::filebacked.big.matrix(80,1200,type='double',backingfile='design.bin',
    descriptorfile='design.desc',backingpath=folder);M[,] <- A
  X <- l0glmData(M,control=list(max.cells=10000,working.columns=100L))
  call <- function(...)l0glmFit(X,y,screen='ISIS',nsis=8L,screen.maxit=3L,
    local.search=FALSE,...)
  for(result in list(call(lambda='aic'),call(maxit=1L),
    call(validation=list(X=A,y=y)),call(prior.mean=0)))
    expect_null(result$screening$working.set)
  start <- rep(1,1201);start[2] <- 0;start[3] <- 0
  fit <- call(start=start,intercept=TRUE,always.include=2L)
  expect_null(fit$screening$working.set)
  expect_equal(as.numeric(coef(fit)[2,1]),0);expect_true(2L%in%fit$screening$indices)
  one <- l0glmFit(X,y,screen='SIS',nsis=8L,local.search=FALSE)
  expect_null(one$screening$working.set)
  two <- l0glmFit(X,y,screen='ISIS',nsis=8L,screen.maxit=2L,local.search=FALSE)
  expect_null(two$screening$working.set)
  for(bad in list(-1,.5,NA,Inf,NULL,'100'))
    expect_error(l0glmData(M,control=list(working.columns=bad)),'working.columns')
})

test_that('cached paths handle null models and no intercept', {
  skip_if_not_installed('bigmemory')
  set.seed(813);A <- matrix(rnorm(100*1600),100,1600)
  folder <- tempfile();dir.create(folder)
  M <- bigmemory::filebacked.big.matrix(100,1600,type='double',backingfile='design.bin',
    descriptorfile='design.desc',backingpath=folder);M[,] <- A
  X <- l0glmData(M,control=list(max.cells=20000,working.columns=150L))
  for(intercept in c(FALSE,TRUE)) {
    fit <- l0glmFit(X,rnorm(100),screen='ISIS',nsis=10L,screen.maxit=4L,
      intercept=intercept,local.search=FALSE)
    expect_identical(dim(coef(fit)),c(1600L+as.integer(intercept),1L))
    expect_true(is.finite(GIC(fit)))
    expect_true(all(is.finite(predict(fit,X))))
  }
})

test_that('tuning and CV preserve original dimensions with cached source controls', {
  skip_if_not_installed('bigmemory')
  set.seed(814);A <- matrix(rnorm(60*900),60,900);y <- 3*A[,1]-2*A[,2]+rnorm(60)
  folder <- tempfile();dir.create(folder)
  M <- bigmemory::filebacked.big.matrix(60,900,type='double',backingfile='design.bin',
    descriptorfile='design.desc',backingpath=folder);M[,] <- A
  X <- l0glmData(M,control=list(max.cells=6000,working.columns=80L))
  for(tune in c('gic','cv')) {
    result <- l0glmTune(X,y,family=gaussian(),tune=tune,nlambda=3L,prop_range=.1,
      nfolds=3L,intercept=TRUE,screen='ISIS',nsis=6L,screen.maxit=3L,
      local.search=FALSE,screen.path='preset',plot=FALSE,printbest=FALSE)
    expect_true(is.finite(result$lambda_best))
    expect_identical(dim(coef(result$fit_best)),c(901L,1L))
    expect_equal(result$fit_best$p,900L)
    expect_true(all(is.finite(predict(result$fit_best,X))))
  }
})
