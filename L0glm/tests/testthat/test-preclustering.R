precluster_fixture<-function(sparse=FALSE) {
  set.seed(921);n<-100L;p<-24L;m<-8L
  X<-matrix(rexp(n*p),n,p);X[X<1]<-0
  if(sparse) X<-Matrix::Matrix(X,sparse=TRUE)
  B<-matrix(0,p,m);B[c(3,9,17),]<-matrix(runif(3*m,.5,2),3,m)
  Y<-matrix(rpois(n*m,as.numeric(X%*%B)+.2),n,m)
  colnames(X)<-paste0("x",seq_len(p));colnames(Y)<-paste0("ion",seq_len(m))
  list(X=X,y=Y,truth=B,group=rep(1:4,each=2))
}

test_that("both public entry points reproduce the same original-channel workflow",{
  for(sparse in c(FALSE,TRUE)) {
    d<-precluster_fixture(sparse)
    controls<-list(clusterids=d$group)
    a<-L0nnpois(d$X,d$y,lambda="gic",offset=.2,precluster=TRUE,
      parallelize=FALSE,control.precluster=controls)
    b<-L0glm.fit(d$X,d$y,family=poisson(identity),nonnegative=TRUE,
      lambda="gic",offset=.2,precluster=TRUE,parallelize=FALSE,control.precluster=controls)
    expect_identical(coef(a),coef(b));expect_s4_class(coef(a),"sparseMatrix")
    expect_equal(dim(coef(a)),c(ncol(d$X),ncol(d$y)))
    expect_identical(colnames(coef(a)),colnames(d$y))
    expect_identical(inherits(a$X,"sparseMatrix"),sparse)
    expect_equal(fitted(a),as.matrix(d$X%*%coef(a))+.2,ignore_attr=TRUE)
    expect_equal(as.matrix(a$weights),1/pmax(fitted(a),.Machine$double.eps),ignore_attr=TRUE)
    expect_equal(a$init.lambda,L0glm:::preset.lambda("gic",ncol(d$X),nrow(d$X)))
    expect_equal(a$preclustering$lambda.screen,.1*a$init.lambda)
    expect_equal(a$preclustering$control$refit.passes,2)
    expect_true(all(coef(a)[setdiff(seq_len(ncol(d$X)),a$preclustering$candidates),,drop=FALSE]==0))
    expect_equal(predict(a,d$X[1:3,,drop=FALSE],offset=.2),fitted(a)[1:3,,drop=FALSE])
    expect_true(is.finite(L0glm.metrics(a)$gic))
  }
})

test_that("preclustering matches an explicit aggregated screen and thin refit",{
  d<-precluster_fixture(TRUE);X<-d$X;Y<-d$y;g<-d$group;m<-ncol(Y)
  lambda<-L0glm:::preset.lambda("gic",ncol(X),nrow(X))
  agg<-glm_aggregate(Y,family=poisson(),clusterids=g)
  reduced<-L0nnpois(X,agg$Y_agg,offset=.4,lambda=.1*lambda,parallelize=FALSE)
  rows<-which(Matrix::rowSums(abs(coef(reduced)))>0)
  start<-matrix(as.numeric(Matrix::rowSums(coef(reduced)[rows,,drop=FALSE]))/m,length(rows),m)
  manual<-L0nnpois(X[,rows,drop=FALSE],Y,offset=.2,lambda=lambda,start=start,
    group.maxit=2L,parallelize=FALSE)
  fit<-L0nnpois(X,Y,offset=.2,lambda=lambda,precluster=TRUE,parallelize=FALSE,
    control.precluster=list(clusterids=g))
  expect_equal(as.matrix(coef(fit)[rows,,drop=FALSE]),as.matrix(coef(manual)),tolerance=0)
  expect_identical(fit$preclustering$candidates,as.integer(rows))
  expect_equal(fit$preclustering$initialization$values,start[,1])
})

test_that("aggregate-only fits retain the target penalty and summed response units",{
  d<-precluster_fixture();n<-nrow(d$X);w<-runif(n,.5,2)
  fit<-L0nnpois(d$X,d$y,weights=w,offset=.2,lambda="gic",precluster=TRUE,
    cluster.refit=FALSE,parallelize=FALSE,control.precluster=list(clusterids=d$group))
  agg<-glm_aggregate(d$y,weights=w,family=poisson(),clusterids=d$group)
  manual<-L0nnpois(d$X,agg$Y_agg,weights=w,offset=.4,
    lambda=fit$preclustering$lambda.target,parallelize=FALSE)
  expect_identical(coef(fit),coef(manual));expect_equal(fit$y,agg$Y_agg)
  expect_equal(dim(coef(fit)),c(ncol(d$X),4L));expect_equal(fit$offset,matrix(.4,n,4))
  expect_equal(fit$prior.weights,matrix(w,n,4));expect_null(fit$preclustering$lambda.screen)
  expect_equal(dim(predict(fit,d$X[1:3,,drop=FALSE])),c(3L,4L))
})

test_that("starts, compulsory terms, priors and normalized scalar one keep their meaning",{
  d<-precluster_fixture(TRUE);p<-ncol(d$X);m<-ncol(d$y)
  S<-matrix(1,p,m);S[5,]<-0;S[3,1:2]<-0;S[1,]<-0
  fit<-L0nnpois(d$X,d$y,start=S,prior.mean=.1,intercept=TRUE,
    always.include=1L,offset=.2,precluster=TRUE,parallelize=FALSE,
    control.precluster=list(clusterids=d$group))
  expect_equal(dim(coef(fit)),c(p+1L,m));expect_true(fit$added_intercept)
  expect_equal(fit$always.include,c(1L,2L));expect_true(all(coef(fit)[6,]==0))
  expect_true(all(coef(fit)[4,1:2]==0));expect_true(2L%in%fit$preclustering$candidates)
  expect_true(all(coef(fit)[1:2,]>=0))
  for(normalize in c("l2","rms")) {
    one<-L0nnpois(d$X,d$y,start=1,normalize=normalize,offset=.2,
      precluster=TRUE,parallelize=FALSE,control.precluster=list(clusterids=d$group))
    rows<-one$preclustering$candidates
    scales<-L0glm:::L0glm_column_scale(d$X,normalize)
    expect_equal(as.numeric(one$preclustering$initialization$values),1/scales[rows])
  }
})

test_that("preclustering threads and disabled defaults preserve fit behaviour",{
  d<-precluster_fixture(TRUE)
  a<-L0nnpois(d$X,d$y,offset=.2,precluster=TRUE,parallelize=FALSE,
    control.precluster=list(k=4L))
  b<-L0nnpois(d$X,d$y,offset=.2,precluster=TRUE,ncores=4L,
    control.precluster=list(k=4L,nthreads=4L))
  expect_identical(coef(a),coef(b))
  base<-L0nnpois(d$X,d$y,offset=.2,parallelize=FALSE)
  off<-L0nnpois(d$X,d$y,offset=.2,parallelize=FALSE,precluster=FALSE)
  expect_identical(coef(base),coef(off));expect_null(off$preclustering)
  zero<-L0nnpois(d$X,d$y,offset=.2,lambda=0,parallelize=FALSE,precluster=TRUE)
  ordinary<-L0nnpois(d$X,d$y,offset=.2,lambda=0,parallelize=FALSE)
  expect_identical(coef(zero),coef(ordinary));expect_false(zero$preclustering$used)
  one<-L0nnpois(d$X,d$y[,1],offset=.2,precluster=TRUE,parallelize=FALSE)
  expect_false(one$preclustering$used)
  expect_equal(L0nnpois(d$X,d$y,offset=.2,precluster=TRUE,group.maxit=3L,
    parallelize=FALSE,control.precluster=list(clusterids=d$group))$preclustering$control$refit.passes,3L)
})

test_that("validation retains original channels or the requested aggregate scale",{
  d<-precluster_fixture();X<-d$X;Y<-d$y;training<-1:75;testing<-76:100
  for(refit in c(TRUE,FALSE)) {
    fit<-L0nnpois(X[training,],Y[training,],offset=.2,intercept=TRUE,
      precluster=TRUE,cluster.refit=refit,parallelize=FALSE,
      control.precluster=list(clusterids=d$group),
      validation=list(X=X[testing,ncol(X):1],y=Y[testing,],offset=.2,patience=2L))
    expect_true(is.finite(fit$validation$deviance))
    expect_equal(ncol(coef(fit)),if(refit) ncol(Y) else 4L)
  }
  colnames(X)<-NULL
  fit<-L0nnpois(X[training,],Y[training,],intercept=TRUE,offset=.2,
    precluster=TRUE,parallelize=FALSE,control.precluster=list(clusterids=d$group),
    validation=list(X=X[testing,],y=Y[testing,],offset=.2))
  expect_true(is.finite(fit$validation$deviance))
})

test_that("preclustering rejects unsupported specifications without changing other families",{
  d<-precluster_fixture()
  expect_error(L0glm.fit(d$X,d$y,precluster=TRUE),"Poisson")
  expect_error(L0glm.fit(d$X,d$y,family=poisson("log"),nonnegative=TRUE,precluster=TRUE),"identity")
  expect_error(L0nnpois(d$X,d$y,precluster="yes"),"TRUE or FALSE")
  expect_error(L0nnpois(d$X,d$y,precluster=TRUE,cluster.refit=NA),"TRUE or FALSE")
  expect_error(L0nnpois(d$X,d$y,precluster=TRUE,weights=matrix(runif(length(d$y)),nrow(d$y))),"shared")
  expect_error(L0nnpois(d$X,d$y,precluster=TRUE,control.precluster=list(unknown=1)),"supported")
  expect_error(L0nnpois(d$X,d$y,precluster=TRUE,control.precluster=list(screen.factor=2)),"screening")
  expect_error(L0nnpois(d$X,d$y,precluster=TRUE,control.precluster=list(clusterids=1)),"clusterids")
})


test_that("working-family controls preserve channel dimensions and dispersion",{
  d<-precluster_fixture(TRUE)
  for(refit in c(TRUE,FALSE)) {
    fit<-L0glm.fit(d$X,d$y,family=poisson(identity),nonnegative=TRUE,
      dispersion=rep(1,ncol(d$y)),lambda="gic",offset=.2,grouped=FALSE,
      openmp=FALSE,parallelize=FALSE,precluster=TRUE,cluster.refit=refit,
      control.precluster=list(clusterids=d$group))
    expect_equal(dim(coef(fit)),c(ncol(d$X),if(refit) ncol(d$y) else 4L))
    expect_equal(length(fit$dispersion),if(refit) ncol(d$y) else 4L)
    expect_true(all(fit$dispersion==1));expect_true(is.finite(L0glm.metrics(fit)$gic))
  }
})
