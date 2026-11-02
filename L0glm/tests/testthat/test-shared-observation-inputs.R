test_that("shared row weights and offsets stay compact and match channel matrices", {
  set.seed(57);n<-70L;p<-8L;m<-3L
  X<-matrix(rexp(n*p),n,p);B<-matrix(0,p,m);B[c(2,6),]<-3
  off<-seq(.5,1.5,length.out=n);w<-runif(n,.7,1.3)
  y<-matrix(rpois(n*m,as.matrix(X%*%B)+off),n,m)
  for(openmp in c(FALSE,TRUE))for(grouped in c(FALSE,TRUE)) {
    args<-list(X=X,y=y,family=poisson("identity"),lambda="gic",nonnegative=TRUE,
      openmp=openmp,parallelize=FALSE,grouped=grouped,relaxed=TRUE)
    a<-do.call(L0glm.fit,c(args,list(weights=w,offset=off)))
    b<-do.call(L0glm.fit,c(args,list(weights=matrix(w,n,m),offset=matrix(off,n,m))))
    expect_equal(coef(a),coef(b),tolerance=1e-8)
    expect_equal(fitted(a),fitted(b),tolerance=1e-8)
    expect_equal(GIC(a),GIC(b),tolerance=1e-8)
    expect_equal(dim(a$prior.weights),c(n,1L))
    expect_equal(dim(a$offset),c(n,1L))
    expect_lt(as.numeric(object.size(a$offset)),as.numeric(object.size(b$offset)))
  }
  expect_error(L0glm.fit(X,y,weights=numeric(2)),"weights")
  expect_error(L0nnpois(X,y,offset=numeric(2)),"offset")
  expect_error(L0nnpois(X,y,weights=-1),"nonnegative")
})
