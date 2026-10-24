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


cluster_loss_reference<-function(Y,group,family,w=1,s=rep(1,ncol(Y))) {
  W<-matrix(w,nrow(Y),ncol(Y));loss<-0
  for(g in unique(group)) {
    js<-which(group==g);mass<-rowSums(sweep(W[,js,drop=FALSE],2,s[js],"*"))
    mu<-rowSums(W[,js,drop=FALSE]*Y[,js,drop=FALSE])/pmax(mass,1e-300)
    for(j in js) {good<-W[,j]>0;means<-s[j]*mu[good]
      terms<-family$dev.resids(Y[good,j],means,W[good,j]);terms[Y[good,j]==means]<-0
      loss<-loss+sum(terms)
    }
  }
  loss
}

test_that("native scalar unit deviances agree with the family functions",{
  set.seed(715)
  families<-list(gaussian(),poisson(),binomial(),Gamma(),inverse.gaussian(),MASS::negative.binomial(2))
  for(f in families) {
    Y<-switch(tolower(f$family),gaussian=matrix(rnorm(55),11),
      poisson=matrix(rpois(55,3),11),binomial=matrix(rbinom(55,1,.4),11),
      matrix(rexp(55)+.2,11))
    if(grepl("negative",f$family,ignore.case=TRUE)) Y<-matrix(rnbinom(55,mu=4,size=2),11)
    for(w in list(1,runif(11),matrix(runif(55),11))) {
      W<-matrix(w,11,5);D<-matrix(0,5,5)
      for(j in 2:5) for(i in seq_len(j-1L)) {
        mu<-(W[,i]*Y[,i]+W[,j]*Y[,j])/(W[,i]+W[,j])
        a<-f$dev.resids(Y[,i],mu,W[,i]);b<-f$dev.resids(Y[,j],mu,W[,j])
        a[Y[,i]==mu]<-0;b[Y[,j]==mu]<-0
        D[i,j]<-D[j,i]<-sum(a+b)
      }
      native<-glm_dist(Y,weights=w,family=f,method="deviance",nthreads=1)
      expect_equal(as.numeric(native),as.numeric(as.dist(D)),tolerance=1e-11)
      expect_identical(native,glm_dist(Y,weights=w,family=f,method="deviance",nthreads=4))
    }
  }
})

test_that("Bregman merges minimise the actual pooled-deviance increase",{
  set.seed(913)
  families<-list(gaussian(),poisson(),binomial(),Gamma(),inverse.gaussian(),MASS::negative.binomial(3))
  for(f in families) {
    Y<-matrix(rexp(36)+.2,6,6)
    if(f$family=="gaussian") Y<-matrix(rnorm(36),6,6)
    if(f$family=="poisson"||grepl("Negative",f$family)) Y<-matrix(rpois(36,4),6,6)
    if(f$family=="binomial") Y<-matrix(rbinom(36,1,.5),6,6)
    W<-matrix(runif(36,.3,2),6,6)
    sizes<-if(f$family=="poisson") runif(6,.4,2) else rep(1,6)
    method<-if(f$family=="poisson") "poissondist" else "deviance"
    out<-glm_cluster(Y,weights=W,family=f,method=method,size_fac=sizes,
      linkage="bregman",k=1L,nthreads=1L)
    groups<-setNames(lapply(1:6,identity),as.character(-(1:6)))
    for(step in 1:5) {
      before<-vapply(groups,function(js) cluster_loss_reference(Y[,js,drop=FALSE],rep(1,length(js)),f,W[,js,drop=FALSE],sizes[js]),numeric(1))
      pairs<-combn(seq_along(groups),2)
      costs<-apply(pairs,2,function(ab) {
        js<-unlist(groups[ab]);cluster_loss_reference(Y[,js,drop=FALSE],rep(1,length(js)),f,W[,js,drop=FALSE],sizes[js])-sum(before[ab])
      })
      merge<-as.character(out$hc$merge[step,]);js<-unlist(groups[merge])
      actual<-cluster_loss_reference(Y[,js,drop=FALSE],rep(1,length(js)),f,W[,js,drop=FALSE],sizes[js])-sum(before[match(merge,names(groups))])
      expect_equal(actual,min(costs),tolerance=1e-10)
      expect_equal(out$hc$height[step],actual,tolerance=1e-10)
      groups[merge]<-NULL;groups[[as.character(step)]]<-js
    }
    expect_equal(sort(out$hc$order),1:6)
    four<-glm_cluster(Y,weights=W,family=f,method=method,size_fac=sizes,
      linkage="bregman",k=1L,nthreads=4L)
    expect_identical(out$hc$merge,four$hc$merge);expect_identical(out$hc$height,four$hc$height)
  }
})

test_that("Gaussian Bregman hierarchy agrees with squared-Euclidean Ward",{
  set.seed(1);Y<-matrix(rnorm(63),9,7)
  out<-glm_cluster(Y,family=gaussian(),k=3L,nthreads=2L)
  ward<-hclust(dist(t(Y))^2/2,method="ward.D")
  expect_equal(out$hc$height,ward$height,tolerance=1e-12)
  for(k in 1:7) {
    a<-cutree(out$hc,k=k);b<-cutree(ward,k=k)
    expect_equal(outer(a,a,"=="),outer(b,b,"=="),ignore_attr=TRUE)
  }
})

test_that("Poisson cut criteria score a common original-channel likelihood",{
  set.seed(16);Y<-matrix(rpois(80*7,3),80,7);s<-runif(7,1,3)
  for(linkage in c("ward","bregman")) {
    fit<-glm_cluster(Y,size_fac=s,linkage=linkage,k.max=5L,nthreads=2L)
    for(j in seq_len(nrow(fit$history))) {
      k<-fit$history$k[j];g<-cutree(fit$hc,k=k);counts<-tabulate(g)
      loss<-cluster_loss_reference(Y,g,poisson(),s=s)
      score<-loss-2*sum(counts*log(counts/7))+log(80)*log(log(80))*(80*k+k-1L)
      expect_equal(fit$history$deviance[j],loss,tolerance=1e-10)
      expect_equal(fit$history$IC[j],score,tolerance=1e-10)
    }
    expect_equal(max(fit$clusterids),fit$history$k[which.min(fit$history$IC)])
  }
  binned<-glm_cluster(Y,size_fac=s,rowbins=4L,k=3L)
  expect_identical(binned$size_fac,s)
  expect_identical(binned$D,glm_dist(rowsum(Y,rep(1:20,each=4)),method="poissondist",size_fac=s))
  expect_equal(dim(binned$Y_agg),c(80L,3L))
})

test_that("aggregation retains correct Poisson weights and Gaussian precisions",{
  Y<-matrix(1:24,6,4);g<-c(10L,10L,20L,20L);w<-seq_len(6)
  a<-glm_aggregate(Y,weights=w,family=poisson(),clusterids=g)
  expect_equal(unname(a$Y_agg),cbind(Y[,1]+Y[,2],Y[,3]+Y[,4]))
  expect_equal(unname(a$weights_agg),unname(cbind(w,w)))
  expect_null(glm_aggregate(Y,family=poisson(),clusterids=g)$weights_agg)
  expect_error(glm_aggregate(Y,weights=matrix(1:24,6,4),family=poisson(),clusterids=g),"shared")
  b<-glm_aggregate(Y,family=gaussian(),clusterids=g)
  expect_equal(unname(b$Y_agg),cbind(rowMeans(Y[,1:2]),rowMeans(Y[,3:4])))
  expect_equal(unname(b$weights_agg),matrix(2,6,2))
  c<-glm_aggregate(Y%%2,family=binomial(),clusterids=g)
  expect_true(all(c$Y_agg>=0 & c$Y_agg<=1));expect_equal(unname(c$weights_agg),matrix(2,6,2))
  one<-glm_cluster(Y[,1,drop=FALSE],k.max=32L)
  expect_equal(one$clusterids,1L);expect_null(one$hc)
  expect_equal(unname(one$Y_agg),Y[,1,drop=FALSE])
  expect_error(glm_cluster(Y,rowbins=1.2),"integers")
  expect_error(glm_cluster(Y,family=gaussian(),rowbins=2),"Poisson")
  expect_error(glm_cluster(Y,dispersion=2),"dispersion")
  expect_error(glm_cluster(Y,k=5),"integers")
})


test_that("clustering remains available when weighted counts cannot be aggregated",{
  Y<-matrix(1:24,6,4);w<-matrix(seq_len(24),6,4)
  cl<-glm_cluster(Y,weights=w,linkage="bregman",k=2L)
  expect_length(cl$clusterids,4L);expect_null(cl$Y_agg)
  expect_match(cl$aggregation_note,"weights differ")
  normal<-glm_cluster(Y,family=gaussian(),dispersion=c(1,2,3,4),k=1)
  W<-matrix(rep(1/c(1,2,3,4),each=6),6,4)
  expect_equal(unname(normal$Y_agg),matrix(rowSums(Y*W)/rowSums(W),6,1))
  expect_equal(unname(normal$weights_agg),matrix(rowSums(W),6,1))
})


test_that("quasi and custom families retain working-deviance clustering",{
  set.seed(51);Y<-matrix(rpois(200,4),20,10)
  q<-glm_cluster(Y,family=quasipoisson())
  expect_equal(length(q$clusterids),10L);expect_identical(q$linkage,"ward")
  expect_equal(as.numeric(glm_dist(Y,family=quasipoisson())),as.numeric(glm_dist(Y,family=poisson())),tolerance=1e-12)
  expect_error(glm_cluster(Y,family=quasipoisson(),criterion="gic"),"likelihood IC")
  fixed<-glm_cluster(Y,family=quasipoisson(),linkage="bregman",k=2)
  expect_true(all(is.na(fixed$history$IC)));expect_equal(max(fixed$clusterids),2)
  custom<-gaussian();custom$family<-"custom scalar"
  z<-glm_cluster(Y,family=custom,k=3)
  expect_equal(max(z$clusterids),3L);expect_identical(z$linkage,"ward")
})

test_that("both Poisson linkages remove pure intensity differences using size factors", {
  Y <- outer(c(0,1,3,8,2),c(1,4,12))
  for (linkage in c("ward","bregman")) {
    cl <- glm_cluster(Y,method="poissondist",linkage=linkage,
      size_fac=colSums(Y),k=1L,nthreads=1L)
    raw <- glm_cluster(Y,method="deviance",linkage=linkage,k=1L,nthreads=1L)
    expect_equal(cl$history$deviance,0,tolerance=1e-10)
    expect_gt(raw$history$deviance,0)
    expect_equal(cl$Y_agg,matrix(rowSums(Y),ncol=1),ignore_attr=TRUE)
  }
})
