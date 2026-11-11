test_that('native support-only Poisson refits match the guarded refit pass', {
  set.seed(7148);n<-100L;p<-20L;m<-4L
  A<-matrix(runif(n*p),n,p);B<-matrix(0,p,m);B[c(2,8),]<-c(10,6)
  Y<-matrix(rpois(n*m,1+A%*%B),n,m)
  for(sparse in c(FALSE,TRUE))for(normalize in c('none','rms')) {
    X<-Matrix::Matrix(A,sparse=sparse)
    args<-list(X=X,y=Y,offset=1,start=1,lambda='gic',grouped=TRUE,group.maxit=2L,
      normalize=normalize,maxit=25L,tol=1e-10,parallelize=FALSE,local.search=FALSE)
    selected<-do.call(L0glm_internal_pois,args)
    old<-do.call(L0glm_internal_pois,c(args,list(relaxed=TRUE)))
    rows<-which(Matrix::rowSums(coef(selected)!=0)>0)
    start<-as.matrix(coef(selected)[rows,,drop=FALSE])
    start[start==0]<-1e-6
    direct<-L0glm_internal_pois(X[,rows,drop=FALSE],Y,offset=1,start=start,lambda=0,
      normalize=normalize,grouped=TRUE,relaxed=TRUE,maxit=25L,tol=1e-10,
      parallelize=FALSE,local.search=FALSE,control.solver=list(.refit.only=TRUE))
    expect_equal(as.matrix(coef(direct)),as.matrix(coef(old)[rows,,drop=FALSE]),tolerance=1e-7)
    expect_true(all(is.finite(direct$weights)))
    expect_true(all(direct$convergence$stage=='relaxed'))
    expect_identical(inherits(direct$X,'sparseMatrix'),sparse)
    automatic<-do.call(l0nnPois,modifyList(args,list(local.search='auto',control.local=list(size.max=0L))))
    expect_lte(GIC(automatic),GIC(selected)+1e-7)
  }
})

test_that('native support-only controls reject unidentified supports', {
  X<-cbind(1:10,1:10);Y<-cbind(1:10,2:11)
  expect_error(L0glm_internal_pois(X,Y,start=matrix(1,2,2),lambda=0,relaxed=TRUE,
    local.search=FALSE,control.solver=list(.refit.only=TRUE)),'insufficient information')
  expect_error(L0glm_internal_pois(diag(3),1:3,lambda=0,relaxed=TRUE,
    local.search=FALSE,control.solver=list(.refit.only=TRUE)),'starting coefficients')
})
