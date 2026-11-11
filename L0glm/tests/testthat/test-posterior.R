test_that('Gaussian posterior moments match the frozen-prior closed form', {
  set.seed(7140);X<-matrix(rnorm(140*5),140,5);w<-runif(140,.5,2);off<-rnorm(140,sd=.1)
  y<-1+2*X[,1]-X[,2]+off+rnorm(140,sd=sqrt(2))
  for(normalize in c('none','rms','l2')) {
    fit<-l0glmFit(X,y,weights=w,offset=off,dispersion=2,intercept=TRUE,normalize=normalize,
      prior.mean=rep(.1,6),lambda='gic',local.search=FALSE)
    original<-fit
    p<-summary(fit,inference='posterior',control.inference=list(nsim=500L,seed=1))$inference
    rows<-p$coefficients$row;A<-as.matrix(fit$X[,rows,drop=FALSE]);D<-diag(p$prior.precision)
    Q<-crossprod(A,A*w)/2+D
    expected<-as.numeric(solve(Q,crossprod(A,w*(y-off))/2+D%*%p$prior.mean))
    expect_equal(p$coefficients$estimate,expected,tolerance=1e-8)
    expect_equal(unname(p$covariance),unname(solve(Q)),tolerance=1e-8)
    expect_true(p$converged)
    expect_equal(fit,original)
    expect_identical(attr(p$draws,'coefficient.dim'),dim(coef(fit)))
    expect_equal(dim(p$draws),c(6L,500L))
    expect_true(all(p$prior.precision[p$coefficients$row==1L]==0))
  }
})

test_that('nonnegative draws sample a truncated Gaussian at a boundary', {
  X<-matrix(1,25,1,dimnames=list(NULL,'Intercept'))
  fit<-l0glmFit(X,rep(-.5,25),nonnegative=TRUE,lambda=1,dispersion=1,local.search=FALSE)
  p<-summary(fit,inference='posterior',control.inference=list(nsim=10000L,seed=12))$inference
  expected<--.5+.2*dnorm(2.5)/pnorm(2.5,lower.tail=FALSE)
  expect_lt(abs(p$coefficients$estimate-expected),.003)
  expect_equal(p$sampling$quadratic$mean,-.5,tolerance=1e-10)
  expect_equal(p$coefficients$mode,0)
  expect_true(all(p$draws>0))
  expect_true(p$coefficients$cred.low>0)
  expect_identical(p$sampler,'truncated-Gaussian Gibbs')
})

test_that('truncated Gibbs preserves correlation and tail stability', {
  V<-matrix(c(1,.7,.7,1),2);Q<-solve(V);mu<-c(-.3,.4)
  set.seed(7141)
  raw<-t(chol(V))%*%matrix(rnorm(8e5),2)+mu
  reference<-raw[,colSums(raw<0)==0,drop=FALSE]
  draws<-posterior_truncated_cpp(mu,Q,c(TRUE,TRUE),c(0,0),30000L,1000L,1L)
  expect_equal(rowMeans(draws),rowMeans(reference),tolerance=.025)
  expect_equal(cov(t(draws)),cov(t(reference)),tolerance=.035)
  extreme<-posterior_truncated_cpp(-100,matrix(1,1),TRUE,0,100L,10L,1L)
  expect_true(all(is.finite(extreme)&extreme>=0))
})

test_that('posterior sampling respects support masks, RNG state and cached models', {
  set.seed(7142);X<-matrix(rnorm(100*5),100,5);y<-2*X[,1]+rnorm(100)
  fit<-l0glmFit(X,cbind(y,-y),lambda=4,start=rbind(1,1,1,1,0),local.search=FALSE)
  state<-.Random.seed
  p<-summary(fit,inference='posterior',control.inference=list(nsim=300L,seed=5))$inference
  a<-simulate(fit,nsim=300L,seed=5)
  b<-simulate(p,nsim=300L,seed=5)
  expect_identical(.Random.seed,state)
  expect_identical(a,b)
  expect_equal(as.matrix(a),as.matrix(p$draws))
  expect_true(all(a[c(5,10),]==0))
  expect_s4_class(a,'sparseMatrix')
  expect_error(simulate(fit,nsim=0),'positive integer')
  expect_error(simulate(fit,seed=NA),'seed')
  expect_error(simulate(p,unknown=1),'Unused')
  expect_error(summary(fit,inference='posterior',control.inference=list(nsim=1)),'nsim')
  expect_error(summary(fit,inference='posterior',control.inference=list(max.parameters=1)),'limit|parameters')
  expect_error(simulate(fit,nsim=1e6,control.inference=list(max.cells=10)),'memory')
  zero<-l0glmFit(X,y,start=0,local.search=FALSE)
  empty<-summary(zero,inference='posterior',control.inference=list(nsim=10L))$inference
  expect_equal(nrow(empty$coefficients),0L)
  expect_true(all(empty$draws==0))
  expect_equal(dim(simulate(zero)),c(5L,1L))
})

test_that('accepted full support refits use a flat prior', {
  set.seed(7143);X<-matrix(rnorm(120*6),120,6);y<-2*X[,1]+rnorm(120)
  fit<-l0glmFit(X,y,lambda='gic',parallelize=FALSE)
  expect_true(fit$support.refit)
  p<-summary(fit,inference='posterior',control.inference=list(nsim=100L))$inference
  expect_identical(p$prior.source,'flat prior on the refitted support')
  expect_true(all(p$prior.precision==0))
  A<-as.matrix(fit$X[,p$coefficients$row,drop=FALSE])
  expected<-lm.fit(A,y)$coefficients
  expect_equal(p$coefficients$estimate,as.numeric(expected),tolerance=1e-7)
})

test_that('automatic multichannel support-only refits use a flat prior', {
  set.seed(7145);n<-200L
  X<-qr.Q(qr(matrix(rnorm(n*30),n,30)))*sqrt(n)
  B<-matrix(2,30,3)
  Y<-X%*%B+matrix(rnorm(n*3),n,3)
  fit<-l0glmFit(X,Y,dispersion=1,parallelize=FALSE)
  expect_identical(fit$local.search$stopping,'support_refit_only')
  expect_true(fit$support.refit)
  p<-summary(fit,inference='posterior',control.inference=list(nsim=100L,seed=8))$inference
  expect_true(all(p$prior.precision==0))
  expect_identical(p$prior.source,'flat prior on the refitted support')
  expect_equal(p$coefficients$estimate,p$coefficients$selected.estimate,tolerance=1e-7)
})

test_that('ordinary and special posterior families converge on their valid spaces', {
  skip_if_not_installed('mgcv');skip_if_not_installed('survival')
  set.seed(7144);n<-120L;X<-matrix(rnorm(n*4),n,4);mu<-exp(.3+.5*X[,1])
  families<-list(poisson(),binomial('probit'),Gamma('log'),quasipoisson(),
    inverse.gaussian('log'),mgcv::nb(),mgcv::multinom(2),mgcv::cox.ph(),
    nb1(2,'log'),gammaScale(.3,'log'),gammaCensored(.1,'log'))
  outcomes<-list(rpois(n,mu),rbinom(n,1,pnorm(.2+X[,1])),rgamma(n,4,scale=mu/4),
    rpois(n,mu),statmod::rinvgauss(n,mu,dispersion=1),rnbinom(n,mu=mu,size=3),
    factor(sample(1:3,n,TRUE)),survival::Surv(rexp(n,exp(X[,1]/2)),rep(c(1,1,0),40)),
    rpois(n,mu),rgamma(n,4,scale=mu/4),{v<-rgamma(n,4,scale=mu/4);ifelse(v<.5,0,v)})
  for(i in seq_along(families)) {
    fit<-l0glmFit(X,outcomes[[i]],family=families[[i]],lambda='aic',
      intercept=i!=8L,local.search=FALSE,parallelize=FALSE)
    p<-summary(fit,inference='posterior',control.inference=list(nsim=100L,seed=6))$inference
    expect_true(p$converged,info=families[[i]]$family)
    expect_true(all(is.finite(p$draws)),info=families[[i]]$family)
    if(i==7L) {
      B<-matrix(p$draws[,1L],nrow(coef(fit)),ncol(coef(fit)))
      expect_equal(rowSums(B),rep(0,nrow(B)),tolerance=1e-10)
    }
  }
})
