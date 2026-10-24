poisson_distance_reference<-function(Y,w=1,s=colSums(Y)) {
  w<-matrix(w,nrow(Y),ncol(Y));D<-matrix(0,ncol(Y),ncol(Y))
  if(ncol(Y)>1) for(j in 2:ncol(Y)) for(k in seq_len(j-1)) {
    den<-w[,j]*s[j]+w[,k]*s[k];rate<-ifelse(den>0,(w[,j]*Y[,j]+w[,k]*Y[,k])/den,0)
    for(z in c(j,k)) {
      mu<-s[z]*rate;good<-w[,z]>0 & mu>0
      ratio<-(w[good,z]*Y[good,z]+1)/(w[good,z]*mu[good]+1)
      D[j,k]<-D[j,k]+sum(w[good,z]*(mu[good]-mu[good]*ratio+Y[good,z]*log(ratio)))
    }
    D[k,j]<-D[j,k]
  }
  as.dist(D)
}

test_that("native Poisson distances preserve the smoothed likelihood gain",{
  set.seed(814);Y<-matrix(rpois(91*9,3),91,9);colnames(Y)<-letters[1:9]
  for(w in list(NULL,2,runif(91),matrix(runif(length(Y)),nrow(Y)))) {
    a<-glm_dist(Y,weights=w,method="poissondist",nthreads=1L)
    b<-poisson_distance_reference(Y,if(is.null(w)) 1 else w)
    expect_equal(as.numeric(a),as.numeric(b),tolerance=1e-12)
    expect_identical(attr(a,"Labels"),colnames(Y))
    expect_identical(a,glm_dist(Y,weights=w,method="poissondist",nthreads=4L))
    expect_equal(a,glm_dist(Matrix::Matrix(Y,sparse=TRUE),weights=w,method="poissondist",nthreads=1L))
  }
  s<-runif(9,.5,3)
  expect_equal(as.numeric(glm_dist(Y,method="poissondist",size_fac=s)),
    as.numeric(poisson_distance_reference(Y,s=s)),tolerance=1e-12)
  expect_equal(glm_dist(Y,method="poissondist",size_fac=s),
    glm_dist(Y,method="poissondist",size_fac=17*s),tolerance=1e-12)
})

test_that("packed distances handle zero and singleton channels",{
  Y<-cbind(a=rep(0,12),b=1:12,c=rep(0,12))
  d<-glm_dist(Y,method="poissondist")
  expect_length(d,3L);expect_true(all(is.finite(d)));expect_true(all(d>=0))
  expect_identical(as.matrix(d),t(as.matrix(d)))
  expect_length(glm_dist(Y[,1,drop=FALSE],method="poissondist"),0L)
  expect_equal(as.numeric(glm_dist(matrix(0,12,3),method="poissondist")),rep(0,3))
  expect_equal(as.numeric(glm_dist(Y,weights=0,method="poissondist")),rep(0,3))
  expect_equal(as.numeric(glm_dist(Y,method="poissondist",size_fac=rep(1,3))),
    as.numeric(poisson_distance_reference(Y,s=rep(1,3))))
})

test_that("clustering inputs fail clearly before native worker dispatch",{
  Y<-matrix(1,4,3)
  expect_error(glm_dist(Y,method="poissondist",size_fac=c(0,1,1)),"positive")
  expect_error(glm_dist(Y,method="poissondist",size_fac=1),"size_fac")
  expect_error(glm_dist(-Y,method="poissondist"),"nonnegative")
  expect_error(glm_dist(Y,weights=-1),"weights")
  expect_error(glm_dist(Y,weights=matrix(1,3,4)),"dimensions")
  expect_error(glm_dist(Y,weights=1:5),"weights")
  expect_error(glm_dist(Y,nthreads=1.5),"nthreads")
  expect_error(glm_dist(matrix(NA_real_,4,3)),"finite")
  expect_error(glm_dist(matrix(numeric(),0,3)),"nonempty")
})

test_that("pooled deviance retains its scalar-family definition",{
  Y<-matrix(c(1,2,4,3,0,1),3,2)
  expect_equal(as.numeric(glm_dist(Y,family=gaussian(),method="deviance")),sum((Y[,1]-Y[,2])^2)/2)
  expect_length(glm_dist(Y[,1,drop=FALSE],method="deviance"),0)
  expect_equal(as.numeric(glm_dist(Y,weights=0,method="deviance")),0)
})
