
nuisance_fixture <- function(sparse=FALSE) {
 set.seed(671);n<-110L;p<-14L;m<-3L
 X<-matrix(runif(n*p),n,p);B<-splines::bs(seq(0,1,length.out=n),df=4,intercept=TRUE)
 beta<-matrix(0,p,m);beta[c(2,8),]<-matrix(c(8,4,9),2,m,byrow=TRUE)
 Y<-matrix(rpois(n*m,X%*%beta+B%*%matrix(rep(c(.2,.6,.8,.3),m),4,m)+.1),n,m)
 Z<-cbind(X,B);colnames(Z)<-paste0("x",seq_len(ncol(Z)))
 if(sparse) Z<-Matrix::Matrix(Z,sparse=TRUE)
 list(X=Z,y=Y,p=p,mandatory=p+1:4)
}
test_that("both entry points use the same original-dimension nuisance workflow", {
 for(sparse in c(FALSE,TRUE)) for(cluster in c(FALSE,TRUE)) {
  d<-nuisance_fixture(sparse)
  opt<-list(X=d$X,y=d$y,always.include=d$mandatory,offset=.1,lambda="gic",
   parallelize=FALSE,precluster=cluster,control.precluster=list(k.max=2L))
  a<-do.call(L0nnpois,opt)
  b<-do.call(L0glm.fit,c(opt,list(family=poisson(identity),nonnegative=TRUE)))
  expect_equal(coef(a),coef(b),tolerance=0)
  expect_equal(dim(coef(a)),c(ncol(d$X),ncol(d$y)))
  expect_identical(inherits(a$X,"sparseMatrix"),sparse)
  expect_s4_class(coef(a),"sparseMatrix")
  expect_equal(fitted(a),as.matrix(d$X%*%coef(a))+.1,ignore_attr=TRUE)
  h<-subset(a$nuisance$history,accepted)
  expect_true(all(diff(h$score)<=1e-8*pmax(1,head(abs(h$score),-1))))
  expect_equal(tail(h$score,1),L0glm.metrics(a,mu_floor=0)$gic,tolerance=1e-9)
 }
})
test_that("starts, normalization, priors, and original offsets survive nuisance fitting", {
 d<-nuisance_fixture(TRUE);p<-ncol(d$X);m<-ncol(d$y)
 S<-matrix(1,p,m);S[4,]<-0;S[2,2]<-0;S[d$mandatory,]<-0
 fit<-L0nnpois(d$X,d$y,start=S,always.include=d$mandatory,lambda="gic",
   offset=matrix(seq(.1,.3,length.out=length(d$y)),nrow(d$y),m),
   prior.mean=matrix(.01,p,m),parallelize=FALSE)
 expect_true(all(coef(fit)[4,]==0));expect_equal(as.numeric(coef(fit)[2,2]),0)
 expect_true(any(coef(fit)[d$mandatory,]>0))
 for(mode in c("l2","rms")) {
  prep<-L0glm:::L0glm_preprocessing(d$X,mode,rep(0,p),d$mandatory)
  a<-L0nnpois(d$X,d$y,start=1,normalize=mode,always.include=d$mandatory,parallelize=FALSE)
  b<-L0nnpois(d$X,d$y,start=matrix(1/prep$scale,p,m),normalize=mode,
    always.include=d$mandatory,parallelize=FALSE)
  expect_equal(coef(a),coef(b),tolerance=0)
 }
})
test_that("single-iteration ridge pilots retain their primitive semantics", {
 d<-nuisance_fixture()
 opt<-list(X=d$X,y=d$y,lambda=2,start=1,always.include=d$mandatory,
  maxit=1L,parallelize=FALSE,relaxed=FALSE)
 a<-do.call(L0nnpois,opt)
 b<-do.call(L0nnpois,c(opt,list(control.solver=list(.nuisance=FALSE))))
 expect_equal(coef(a),coef(b),tolerance=0);expect_null(a$nuisance)
})
test_that("Gaussian compulsory solves retain known dispersions and signed binomial intercepts", {
 d<-nuisance_fixture();set.seed(82)
 Y<-as.matrix(d$X[,c(2,8)]%*%matrix(c(2,1,3),2,3,byrow=TRUE))+matrix(rnorm(length(d$y)),nrow(d$y),3)
 fit<-L0glm.fit(d$X,Y,family=gaussian(),nonnegative=TRUE,dispersion=c(1,2,3),
   always.include=d$mandatory,lambda="gic",parallelize=FALSE,relaxed=TRUE)
 expect_equal(fit$dispersion,c(1,2,3))
 expect_true(all(coef(fit)>=0))
 set.seed(83);X<-matrix(runif(400*5),400,5)
 y<-rbinom(400,1,plogis(-2+X[,1]+X[,3]))
 b<-L0glm.fit(X,y,family=binomial("logit"),intercept=TRUE,always.include=1L,
  nonnegative=TRUE,lambda="aic",parallelize=FALSE,relaxed=TRUE)
 expect_lt(coef(b)[1,1],0);expect_true(all(coef(b)[-1,]>=0))
})
test_that("zero compulsory columns and automatic intercepts keep valid dimensions", {
 set.seed(84);X<-cbind(matrix(runif(120*4),120,4),zero=0)
 y<-rpois(120,2*X[,1]+1)
 fit<-L0nnpois(X,y,always.include=5L,offset=.5,lambda=2,parallelize=FALSE)
 expect_equal(as.numeric(coef(fit)[5,1]),0);expect_true(all(is.finite(fitted(fit))))
 X<-cbind(1,X[,-5]);auto<-L0nnpois(X,y,parallelize=FALSE)
 expect_equal(auto$nuisance$terms,1L)
})
test_that("negative-binomial shape is optimized without changing the supplied family", {
 skip_if_not_installed("mgcv");set.seed(85)
 B<-splines::bs(seq(0,1,length.out=250),df=4,intercept=TRUE)
 y<-rnbinom(250,mu=as.numeric(B%*%c(2,8,10,3)),size=7)
 f<-mgcv::nb(link="identity");before<-f$getTheta(TRUE)
 fit<-L0glm.fit(B,y,family=f,nonnegative=TRUE,always.include=1:4,lambda="gic",parallelize=FALSE)
 expect_identical(f$getTheta(TRUE),before)
 expect_true(is.finite(fit$theta)&&fit$theta>0);expect_true(fit$theta.estimated)
 expect_true(all(coef(fit)>=0))
 fixed<-L0glm.fit(B,y,family=mgcv::nb(theta=7,link="identity"),nonnegative=TRUE,
   always.include=1:4,lambda="gic",parallelize=FALSE)
 expect_equal(fixed$theta,7);expect_false(fixed$theta.estimated)
})
test_that("Cox nuisance updates retain partial-likelihood predictions", {
 skip_if_not_installed("mgcv");skip_if_not_installed("survival");set.seed(86)
 X<-matrix(runif(140*7),140,7);event<-rexp(140,exp(X[,1]+X[,3]));censor<-rexp(140,.4)
 y<-survival::Surv(pmin(event,censor),event<censor)
 fit<-L0glm.fit(X,y,family=mgcv::cox.ph(),always.include=1L,nonnegative=TRUE,
  lambda="aic",parallelize=FALSE,relaxed=TRUE,maxit=40L)
 expect_equal(dim(coef(fit)),c(7L,1L));expect_true(all(coef(fit)>=0))
 expect_true(is.finite(L0glm.metrics(fit)$loglik));expect_true(all(is.finite(predict(fit,X,type="risk"))))
 expect_warning(with_intercept<-L0glm.fit(X,y,family=mgcv::cox.ph(),always.include=1L,
   nonnegative=TRUE,intercept=TRUE,lambda="aic",parallelize=FALSE,relaxed=TRUE,maxit=40L),
   "not identifiable")
 expect_equal(coef(with_intercept),coef(fit),tolerance=0)
})


test_that("conditional scoring supports nonidentity scalar-family links", {
  skip_if_not_installed("mgcv");skip_if_not_installed("statmod")
  set.seed(972);n<-150L;B<-matrix(runif(n*3,.2,1),n,3)
  offset<-matrix(.2,n,1);W<-matrix(runif(n,.5,2),n,1)
  fs<-list(gaussian("log"),poisson("log"),poisson("sqrt"),Gamma("log"),
    Gamma("inverse"),inverse.gaussian("log"),mgcv::nb(theta=10,link="log"))
  for(f in fs) {
    mu<-f$linkinv(as.numeric(offset+B%*%c(.2,.6,.9)))
    Y<-matrix(if(f$family=="gaussian") rnorm(n,mu,.3) else
      if(f$family=="poisson") rpois(n,mu) else
      if(f$family=="Gamma") rgamma(n,shape=12,scale=mu/12) else
      if(f$family=="inverse.gaussian") statmod::rinvgauss(n,mean=mu,shape=20) else
      rnbinom(n,mu=mu,size=10),n,1)
    storage.mode(Y)<-"double"
    fit<-L0glm:::L0glm_nuisance_conditional(B,Y,offset,W,matrix(.5,3,1),rep(0,3),list(f))
    loss<-function(a) L0glm:::L0glm_irls_deviance(Y,offset+B%*%a,W,f)
    oracle<-optim(as.numeric(fit$coefficients),loss,method="L-BFGS-B",lower=rep(0,3),
      control=list(factr=1,pgtol=1e-10))
    expect_lte(loss(fit$coefficients),oracle$value+1e-6*max(1,abs(oracle$value)))
    expect_true(all(fit$coefficients>=0));expect_lt(fit$kkt[1],1e-5)
  }
})

test_that("quasi families and data-frame designs retain valid fitted means", {
  set.seed(973);X<-as.data.frame(matrix(runif(200*6),200,6))
  y<-rpois(200,1+4*X[[2]])
  for(f in list(quasipoisson("identity"),quasipoisson("log"))) {
    fit<-L0glm.fit(X,y,family=f,nonnegative=TRUE,intercept=TRUE,
      always.include=1L,lambda=2,parallelize=FALSE,relaxed=TRUE)
    expect_true(all(is.finite(fitted(fit))));expect_true(all(coef(fit)>=0))
    expect_true(all(is.finite(fit$dispersion)))
  }
})
