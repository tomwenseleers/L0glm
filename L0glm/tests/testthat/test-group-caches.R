test_that('cached group solves agree with equivalent family callbacks', {
  set.seed(920);n <- 140L;p <- 24L
  X <- matrix(rnorm(n*p),n,p);X[,p] <- 0
  Y <- cbind(1+X[,1]-X[,2]+rnorm(n),-1+X[,1]+X[,2]+rnorm(n))
  w <- cbind(rep(c(0,1,2),length.out=n),rep(c(2,1,0),length.out=n))
  off <- cbind(runif(n),runif(n))
  start <- matrix(1,p+1L,2L);start[5,1] <- 0
  callback <- gaussian();callback$variance <- function(mu)rep(1,length(mu))
  for(sparse in c(FALSE,TRUE))for(grouped in c(FALSE,TRUE)) {
    A <- if(sparse)Matrix::Matrix(X,sparse=TRUE)else X
    args <- list(X=A,y=Y,weights=w,offset=off,group=rep(1:8,each=3),
      start=start,always.include=22L,intercept=TRUE,dispersion=c(.7,1.2),
      lambda='gic',relaxed=TRUE,grouped=grouped,prior.mean=.03,
      maxit=80L,local.search=FALSE,parallelize=FALSE)
    cached <- do.call(l0glmFit,c(args,list(family=gaussian())))
    reference <- do.call(l0glmFit,c(args,list(family=callback)))
    expect_equal(coef(cached),coef(reference),tolerance=1e-9)
    expect_equal(GIC(cached),GIC(reference),tolerance=1e-10)
    expect_equal(predict(cached),predict(reference),tolerance=1e-9)
    expect_identical(as.matrix(coef(cached)!=0),as.matrix(coef(reference)!=0))
    expect_equal(dim(coef(cached)),c(p+1L,2L))
    expect_equal(as.numeric(coef(cached)[5,1]),0)
  }
})

test_that('wide unpenalized group fits retain observation-space least squares', {
  set.seed(921);X <- matrix(rnorm(18*50),18,50);y <- rnorm(18)
  fit <- l0glmFit(X,y,group=rep(1:10,each=5),intercept=TRUE,lambda=0,
    dispersion=1,start=1,relaxed=TRUE,local.search=FALSE,parallelize=FALSE)
  expect_equal(as.numeric(fitted(fit)),y,tolerance=1e-8)
  expect_true(all(is.finite(coef(fit))))
})

test_that('group prediction caches are isolated between fits and offsets', {
  set.seed(922);X <- matrix(rnorm(150*12),150,12)
  y <- rbinom(150,1,plogis(.4+X[,1]-X[,2]));off <- seq(-.4,.4,length.out=150)
  call <- function(y,offset)l0glmFit(X,y,offset=offset,family=binomial(),
    group=rep(1:4,each=3),intercept=TRUE,lambda='gic',relaxed=TRUE,local.search=FALSE)
  first <- call(y,off);first_prediction <- fitted(first)
  second <- call(1-y,-off)
  repeated <- call(y,off)
  expect_equal(coef(first),coef(repeated),tolerance=1e-12)
  expect_equal(fitted(first),first_prediction,tolerance=1e-12)
  expect_equal(as.numeric(fitted(second)),1-as.numeric(first_prediction),tolerance=1e-7)
})
