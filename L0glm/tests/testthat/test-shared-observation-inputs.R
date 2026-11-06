test_that("shared row weights and offsets stay compact and match channel matrices", {
  set.seed(57);n<-70L;p<-8L;m<-3L
  X<-matrix(rexp(n*p),n,p);B<-matrix(0,p,m);B[c(2,6),]<-3
  off<-seq(.5,1.5,length.out=n);w<-runif(n,.7,1.3)
  y<-matrix(rpois(n*m,as.matrix(X%*%B)+off),n,m)
  for(openmp in c(FALSE,TRUE))for(grouped in c(FALSE,TRUE)) {
    args<-list(X=X,y=y,family=poisson("identity"),lambda="gic",nonnegative=TRUE,
      openmp=openmp,parallelize=FALSE,grouped=grouped,relaxed=TRUE)
    a<-do.call(l0glmFit,c(args,list(weights=w,offset=off)))
    b<-do.call(l0glmFit,c(args,list(weights=matrix(w,n,m),offset=matrix(off,n,m))))
    expect_equal(coef(a),coef(b),tolerance=1e-8)
    expect_equal(fitted(a),fitted(b),tolerance=1e-8)
    expect_equal(GIC(a),GIC(b),tolerance=1e-8)
    expect_equal(dim(a$prior.weights),c(n,1L))
    expect_equal(dim(a$offset),c(n,1L))
    expect_lt(as.numeric(object.size(a$offset)),as.numeric(object.size(b$offset)))
  }
  expect_error(l0glmFit(X,y,weights=numeric(2)),"weights")
  expect_error(l0nnPois(X,y,offset=numeric(2)),"offset")
  expect_error(l0nnPois(X,y,weights=-1),"nonnegative")
  groups<-Matrix::sparseMatrix(i=1:3,j=c(1L,1L,2L),x=1)
  expect_equal(L0glm_precluster_offset(off,n,m,groups),as.matrix(matrix(off,n,m)%*%groups))
})

test_that("tuning preserves compact shared inputs and equivalent CV scores", {
  set.seed(558); n<-48L; p<-4L; m<-3L
  X<-matrix(rexp(n*p),n,p); B<-matrix(0,p,m);B[2L,]<-c(2,3,4)
  off<-seq(.3,.8,length.out=n);w<-rep(c(0,.7,1.2,1.7),length.out=n)
  y<-matrix(rpois(n*m,X%*%B+off),n,m)
  args<-list(X=X,y=y,family=poisson("identity"),tune="cv",nlambda=3L,
    foldid=rep(1:3,each=n/3),nfolds=3L,nonnegative=TRUE,
    grouped=TRUE,openmp=TRUE,parallelize=FALSE,plot=FALSE,printbest=FALSE)
  a<-do.call(l0glmTune,c(args,list(weights=w,offset=off)))
  b<-do.call(l0glmTune,c(args,list(weights=matrix(w,n,m),offset=matrix(off,n,m))))
  expect_equal(a$cv$deviance,b$cv$deviance,tolerance=1e-10)
  expect_equal(a$cv$fold_weights,b$cv$fold_weights,tolerance=1e-12)
  expect_equal(a$lambda_best,b$lambda_best)
  expect_equal(a$lambda_1se,b$lambda_1se)
  expect_equal(coef(a$fit_best),coef(b$fit_best),tolerance=1e-10)
  expect_equal(dim(a$fit_args$weights),c(n,1L))
  expect_equal(dim(a$fit_args$offset),c(n,1L))
  expect_equal(predict(a$fit_best,X,offset=off),
    predict(a$fit_best,X,offset=matrix(off,n,m)),tolerance=1e-12)
})

test_that("negative-binomial CV accounts for every channel with shared weights", {
  skip_if_not_installed("mgcv")
  set.seed(559);n<-36L;m<-2L;X<-matrix(rnorm(n*3L),n,3L)
  w<-rep(c(0,.7,1.2),length.out=n);off<-seq(-.2,.2,length.out=n)
  y<-matrix(rnbinom(n*m,size=5,mu=exp(1+X[,2L]+off)),n,m)
  args<-list(X=X,y=y,family=mgcv::nb(theta=5,link="log"),tune="cv",
    nlambda=3L,foldid=rep(1:3,each=n/3),nfolds=3L,maxit=15L,
    parallelize=FALSE,plot=FALSE,printbest=FALSE)
  a<-do.call(l0glmTune,c(args,list(weights=w,offset=off)))
  b<-do.call(l0glmTune,c(args,list(weights=matrix(w,n,m),offset=matrix(off,n,m))))
  expect_equal(a$cv$deviance,b$cv$deviance,tolerance=1e-10)
  expect_equal(a$cv$fold_weights,b$cv$fold_weights,tolerance=1e-12)
  expect_equal(a$lambda_best,b$lambda_best)
  expect_equal(dim(a$fit_args$weights),c(n,1L))
  expect_equal(dim(a$fit_args$offset),c(n,1L))
})
