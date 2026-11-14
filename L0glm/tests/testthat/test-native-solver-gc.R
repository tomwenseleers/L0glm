test_that('native solver callback arguments survive forced garbage collection', {
  set.seed(1217)
  X<-matrix(runif(40L*6L,.1,1),40L,6L)
  y<-1+rowSums(X);prior<-seq_len(6L)/10
  received<-new.env(parent=emptyenv())
  factory<-function(solver,control) {
    # Collection during argument wrapping exposes unprotected temporary SEXPs.
    gctorture2(1L)
    function(A,response,penalty,center,guess,lower,upper,iteration) {
      gctorture2(0L)
      received$arguments<-list(A=A,response=response,penalty=penalty,center=center,
        guess=guess,lower=lower,upper=upper,iteration=iteration)
      as.numeric(guess)
    }
  }
  local_mocked_bindings(L0glm_native_solver_bridge=factory)
  on.exit(gctorture2(0L),add=TRUE)
  for(sparse in c(FALSE,TRUE)) {
    A<-if(sparse)Matrix::Matrix(X,sparse=TRUE)else X
    fit<-l0glmFit(A,y,family=gaussian(),dispersion=1,start=1,offset=1,
      prior.mean=prior,lambda=1,nonnegative=TRUE,normalize='none',
      solver='bigGlm',maxit=1L,openmp=TRUE,parallelize=FALSE,nthreads=1L,
      local.search=FALSE,screen='none',intercept=FALSE)
    gctorture2(0L)
    args<-received$arguments
    expect_equal(as.matrix(args$A),X)
    expect_equal(as.numeric(args$response),as.numeric(y-1))
    expect_equal(as.numeric(args$penalty),rep(.25,6))
    expect_equal(as.numeric(args$center),prior)
    expect_equal(as.numeric(args$guess),rep(1,6))
    expect_equal(as.numeric(args$lower),rep(0,6))
    expect_equal(as.numeric(args$upper),rep(.Machine$double.xmax,6))
    expect_equal(as.numeric(coef(fit)),rep(1,6))
  }
})

test_that('native bigGlm matches R fits for correlated dense and sparse designs', {
  skip_if_not_installed('glmnet')
  for(family in list(poisson(identity),gaussian(identity))) {
    poisson_case<-family$family=='poisson'
    d<-simul(method=if(poisson_case)2L else 1L,n=1000L,p=1000L,
      sparse=TRUE,seed=1L,options=list(sd_noise=1))
    for(sparse in c(TRUE,FALSE)) {
      X<-if(sparse)d$X else as.matrix(d$X)
      args<-list(X=X,y=d$y,family=family,lambda='gic',start=NULL,
        nonnegative=TRUE,local.search=FALSE,solver='bigGlm',parallelize=FALSE,nthreads=1L)
      r<-do.call(l0glmFit,c(args,list(openmp=FALSE)))
      native<-do.call(l0glmFit,c(args,list(openmp=TRUE)))
      expect_equal(as.numeric(coef(native)),as.numeric(coef(r)),tolerance=1e-8)
      expect_identical(as.numeric(coef(native))!=0,as.numeric(coef(r))!=0)
      met<-l0glmMetrics(native,beta_true=d$beta_true)
      expect_gt(met$R2_coefs,if(poisson_case).95 else .99999)
    }
  }
})
