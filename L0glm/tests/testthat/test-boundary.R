
test_that("boundary tests improve the conditional objective without size cutoffs", {
  X <- matrix(1,6,1); weights <- c(0,rep(1,5)); lower <- 0
  for(sparse in c(FALSE,TRUE)) for(family in list(gaussian(),poisson("identity"))) {
    Z <- Matrix(X,sparse=sparse)
    out <- L0glm_boundary(Z,rep(0,6),1e-18,rep(1e-18,6),weights,family,lower)
    expect_equal(out$beta,0);expect_equal(out$removed,1L)
  }
  for(scale in c(1,1e12,1e25)) {
    X <- matrix(scale,6,1)
    out <- L0glm_boundary(X,rep(1,6),1/scale,rep(1,6),rep(1,6),poisson("identity"),0)
    expect_identical(out$beta,1/scale)
  }
  # A negative design contribution can be removed when the endpoint is better.
  out <- L0glm_boundary(matrix(-1,6,1),rep(2,6),1,rep(1,6),rep(1,6),poisson("identity"),0)
  expect_equal(out$beta,0)
  signed <- L0glm_boundary(matrix(1,6,1),rep(0,6),1,rep(1,6),rep(1,6),gaussian(),-Inf)
  expect_equal(signed$beta,1)
})

test_that("family boundary passes preserve domains and reduce deviance", {
  set.seed(119)
  n<-80L;p<-6L;X<-matrix(runif(n*p),n,p);beta<-rep(.2,p);w<-runif(n);w[1:10]<-0
  families<-list(gaussian(),gaussian("log"),gaussian("inverse"),
    poisson("identity"),poisson("log"),poisson("sqrt"),
    Gamma("identity"),Gamma("log"),Gamma("inverse"),
    inverse.gaussian("identity"),inverse.gaussian("log"),inverse.gaussian("inverse"),
    inverse.gaussian("1/mu^2"),binomial("logit"),binomial("probit"),
    binomial("cloglog"),binomial("cauchit"),binomial("log"),binomial("identity"),
    quasipoisson("log"),quasibinomial("logit"))
  for(family in families) {
    off<-if(family$link=="log" && family$family=="binomial") -3 else
      if(family$link=="identity" && family$family=="binomial") .1 else 1
    if(family$family=="binomial" && family$link=="identity") b<-beta/10 else b<-beta
    eta<-as.numeric(X%*%b)+off
    mu<-family$linkinv(eta)
    y<-if(family$family %in% c("binomial","quasibinomial")) pmin(1,pmax(0,mu-.1)) else pmax(.1,mu-.5)
    for(sparse in c(FALSE,TRUE)) {
      result<-L0glm_boundary(Matrix(X,sparse=sparse),y,b,eta,w,family,rep(0,p))
      expect_lte(L0glm_irls_deviance(y,result$eta,w,family),
        L0glm_irls_deviance(y,eta,w,family)+1e-10)
      expect_equal(result$eta,as.numeric(X%*%result$beta)+off,tolerance=1e-12)
      expect_true(all(result$beta>=0))
    }
  }
  expect_equal(L0glm_boundary(matrix(1,10,1),rep(1,10),1,rep(1,10),rep(1,10),
    Gamma("identity"),0)$beta,1)
})

test_that("scaling, storage and zero starts are retained", {
  set.seed(120)
  X<-matrix(rexp(240),40,6);beta<-c(.1,0,.2,.3,1e-18,.1)
  eta<-as.numeric(.3+X%*%beta);y<-rpois(40,eta);w<-runif(40)
  for(family in list(gaussian(),poisson("identity"),poisson("log"),binomial())) {
    response<-if(family$family=="binomial") as.numeric(y>0) else y
    a<-L0glm_boundary(X,response,beta,eta,w,family,rep(0,6))
    d<-exp(seq(-10,10,length.out=6))
    b<-L0glm_boundary(sweep(X,2,d,"*"),response,beta/d,eta,w,family,rep(0,6))
    c<-L0glm_boundary(as(Matrix(X,sparse=TRUE),"dgCMatrix"),response,beta,eta,w,family,rep(0,6))
    expect_equal(a$beta,b$beta*d,tolerance=1e-9)
    expect_equal(a$beta,c$beta,tolerance=1e-12)
    expect_equal(a$beta[2],0)
  }
})

test_that("ordinary relaxed fits refresh their working quantities", {
  set.seed(121)
  X<-matrix(runif(200*8),200,8);beta<-c(2,0,1,rep(0,5));offset<-rep(.2,200)
  for(family in list(gaussian(),poisson("identity"),poisson("log"),binomial())) {
    eta<-as.numeric(X%*%beta)+offset;mu<-family$linkinv(eta)
    y<-switch(family$family,gaussian=rnorm(200,mu),poisson=rpois(200,mu),binomial=rbinom(200,1,mu))
    fit<-L0glm.fit(X,y,family=family,offset=offset,nonnegative=TRUE,intercept=TRUE,
      openmp=FALSE,relaxed=TRUE,parallelize=FALSE,lambda=2)
    expect_true(all(is.finite(coef(fit))))
    expect_s4_class(coef(fit),"sparseMatrix")
    expect_equal(predict(fit,type="link"),as.numeric(fit$X%*%coef(fit))+offset,ignore_attr=TRUE)
    if(family$family=="binomial") expect_true(is.finite(coef(fit)[1]))
  }
})

test_that("special-family boundary passes use the fitted nuisance parameters", {
  skip_if_not_installed("mgcv");skip_if_not_installed("survival")
  set.seed(122)
  X<-matrix(runif(120*5),120,5);b<-matrix(c(.2,1e-18,.1,0,.1),5,1)
  y<-rnbinom(120,mu=exp(X%*%b+.2),size=2)
  for(fixed in c(TRUE,FALSE)) {
    family<-mgcv::nb(theta=if(fixed) 2 else -2)
    data<-L0glm_ext_prepare(X,y,family,NULL,.2,FALSE,NULL,TRUE,"none",1,NULL,list())
    result<-list(B=b,state=list(theta=2,estimated=!fixed))
    out<-L0glm_ext_boundary(data,result)
    expect_gte(L0glm_ext_loglik(data,out$B,out$state),L0glm_ext_loglik(data,b,result$state)-1e-10)
    expect_equal(out$state,result$state)
    fit<-L0glm.fit(X,y,family=family,offset=.2,nonnegative=TRUE,relaxed=TRUE,
      parallelize=FALSE,lambda=2,maxit=40)
    expect_true(is.finite(as.numeric(logLik(fit))))
  }
  event<-rexp(120,exp(X%*%b));y<-survival::Surv(round(event,1),rep(c(1,1,0),40))
  data<-L0glm_ext_prepare(X,y,mgcv::cox.ph(),NULL,.2,FALSE,NULL,TRUE,"none",1,NULL,list())
  result<-list(B=b,state=list());out<-L0glm_ext_boundary(data,result)
  expect_gte(L0glm_ext_loglik(data,out$B,out$state),L0glm_ext_loglik(data,b,result$state)-1e-10)
  fit<-L0glm.fit(X,y,family=mgcv::cox.ph(),offset=.2,nonnegative=TRUE,relaxed=TRUE,
    parallelize=FALSE,lambda=2,maxit=40)
  expect_true(is.finite(as.numeric(logLik(fit))))
})
test_that("validation rejects a training-only boundary improvement", {
  X<-matrix(1,12,1);V<-matrix(1,8,1);family<-poisson("identity")
  context<-L0glm_validation_context(list(X=V,y=rep(1,8),offset=.1),X,matrix(0,12,1),family,FALSE)
  context$baseline<-context$refit_baseline<-matrix(1,1,1)
  r<-L0glm.fit.internal(X,rep(0,12),family=family,start=1,offset=rep(.1,12),
    lower=0,lambda=1e-30,hard_threshold=FALSE,validation=context)
  native<-L0nnpois_dense_cpp(X,matrix(0,12,1),matrix(1,12,1),matrix(.1,12,1),
    matrix(0,1,1),matrix(1,1,1),integer(),1e-30,0L,1,FALSE,TRUE,1L,10L,1e-5,1e-10,1L,FALSE,context)
  expect_equal(as.numeric(r$coefficients),1)
  expect_equal(as.numeric(native$coefficients),1)
})
test_that("NB endpoint tests retain fitted theta and remove redundant positive terms", {
  skip_if_not_installed("mgcv")
  X<-matrix(1,15,2);b<-matrix(c(1e-18,.5),2,1)
  for(link in c("log","identity")) {
    data<-L0glm_ext_prepare(X,rep(0,15),do.call(mgcv::nb,list(theta=-1,link=link)),NULL,.3,
      FALSE,NULL,TRUE,"none",1,NULL,list())
    state<-list(theta=7,estimated=TRUE)
    z<-L0glm_ext_boundary(data,list(B=b,state=state))
    expect_true(all(z$B==0));expect_identical(z$state,state)
  }
})
test_that("Poisson subtraction respects signed offsets and accumulated roundoff", {
  X<-matrix(1,6,200);beta<-rep(1e-18,200);eta<-as.numeric(X%*%beta)
  for(sparse in c(FALSE,TRUE)) {
    result<-L0glm_boundary(Matrix::Matrix(X,sparse=sparse),rep(0,6),beta,eta,
      rep(1,6),poisson("identity"),rep(0,200),rep(0,6))
    expect_true(all(result$beta==0))
  }
  result<-L0glm_boundary(matrix(1,6,1),rep(0,6),1.1,rep(1,6),rep(1,6),
    poisson("identity"),0,rep(-.1,6))
  expect_equal(result$beta,1.1)
})


test_that("native boundary weights agree with the final reconstructed signal", {
  n <- 12L; p <- 200L; m <- 2L
  for (sparse in c(FALSE,TRUE)) {
    X <- if (sparse) as(Matrix::Matrix(matrix(1,n,p),sparse=TRUE),"dgCMatrix") else matrix(1,n,p)
    native <- if (sparse) L0nnpois_sparse_cpp else L0nnpois_dense_cpp
    fit <- withCallingHandlers(native(X,matrix(0,n,m),matrix(1,n,m),matrix(0,n,m),
      matrix(0,p,m),matrix(1e-18,p,m),integer(),1e-30,0L,1,FALSE,TRUE,
      2L,2L,1e-5,1e-10,1L,FALSE,NULL),warning=function(w) {
        expect_match(conditionMessage(w),"No decreasing finite Poisson step",fixed=TRUE)
        invokeRestart("muffleWarning")
      })
    mu <- as.matrix(X%*%fit$coefficients)
    expect_equal(as.matrix(fit$weights),1/pmax(mu,.Machine$double.eps))
    expect_true(all(fit$coefficients==0))
  }
})
