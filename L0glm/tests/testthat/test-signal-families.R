test_that("signal-family densities and scores match the specified distributions", {
  y <- c(0, 1, 3, 7, 40); mu <- c(.1, .5, 4, 9, 50)
  families <- list(nb1(3), Gamma_scale(2, censor=.5), Gamma_censored(.3))
  reference <- list(dnbinom(y, size=mu/2, mu=mu, log=TRUE),
    ifelse(y==0, pgamma(.5, shape=mu/2, scale=2, log.p=TRUE), dgamma(y, shape=mu/2, scale=2, log=TRUE)),
    ifelse(y==0, pgamma(.5, shape=1/.3, scale=mu*.3, log.p=TRUE), dgamma(y, shape=1/.3, scale=mu*.3, log=TRUE)))
  for (j in seq_along(families)) {
    f <- families[[j]]
    expect_equal(f$log_density(y,mu),reference[[j]],tolerance=1e-10)
    eps <- mu*1e-5
    numerical <- -(f$log_density(y,mu+eps)-f$log_density(y,mu-eps))/(2*eps)
    expect_equal(f$score(y,mu)[,1], numerical, tolerance=2e-6)
    expect_true(all(f$dev.resids(y,mu,rep(1,5))>=0))
  }
  expect_equal(nb1(1)$log_density(y,mu),dpois(y,mu,log=TRUE))
  expect_equal(nb1(1.000001)$log_density(y,mu),dnbinom(y,size=mu/1e-6,mu=mu,log=TRUE),tolerance=1e-8)
  expect_error(nb1(.5),"noise parameter")
  expect_error(Gamma_scale(-1),"noise parameter")
  expect_identical(nb1(3,identity)$link,"identity")
  expect_identical(Gamma_scale(2,log)$link,"log")
})

test_that("the exact family scores fit, screen, tune and validate", {
  set.seed(831)
  X <- matrix(rnorm(180*5),180,5); mu <- exp(.8+.25*X[,1]-.2*X[,2])
  ys <- list(rnbinom(180, size=mu/2,mu=mu),rgamma(180,shape=mu/2,scale=2),
    rgamma(180,shape=4,scale=mu/4))
  fs <- list(nb1(3), Gamma_scale(2,censor=.5),Gamma_censored(.25))
  for (j in seq_along(fs)) {
    f <- fs[[j]]; y <- ys[[j]]
    if(j>1) y[y<.5] <- 0
    fit <- L0glm.fit(X,y,family=f,intercept=TRUE,lambda=0,maxit=100,tol=1e-9)
    objective <- function(b) -sum(f$log_density(y,exp(cbind(1,X)%*%b)))
    optimal <- optim(as.numeric(coef(fit)),objective,method="BFGS",control=list(reltol=1e-11))
    expect_equal(-as.numeric(logLik(fit)),optimal$value,tolerance=2e-5)
    expect_equal(as.numeric(logLik(fit)),sum(f$log_density(y,fitted(fit))),tolerance=1e-9)
    screen <- sis_glm(X,y,family=f,nsis=3L)
    expect_true(length(screen)>0)
    for(tune in c("gic","cv")) {
      z <- L0glm.tune(X,y,family=f,intercept=TRUE,nlambda=3L,tune=tune,
        nfolds=2L,plot=FALSE,printbest=FALSE,maxit=60)
      expect_true(is.finite(as.numeric(logLik(z$fit_best))))
    }
    a <- L0glm.fit(X,y,family=f,intercept=TRUE,lambda=2,
      validation=list(X=X,y=y),maxit=40)
    expect_true(is.finite(as.numeric(logLik(a))))
  }
})

test_that("signal clustering uses fitted likelihood ratios and is deterministic", {
  set.seed(955)
  Y <- matrix(rpois(60,5),10,6);W <- matrix(runif(60),10,6);W[1,] <- 0
  for(f in list(nb1(3),Gamma_scale(2,censor=.5),Gamma_censored(.3))) {
    a <- glm_dist(Y,weights=W,family=f,nthreads=1L)
    b <- glm_dist(Y,weights=W,family=f,nthreads=2L)
    expect_identical(a,b)
    expect_true(all(is.finite(a)&a>=0))
    for(i in 1:5) for(j in (i+1):6) {
      separate <- shared <- 0
      for(r in 2:10) {
        ll <- function(mu,y,w) sum(w*f$log_density(y,rep(mu,length(y))))
        sep <- vapply(c(i,j),function(k) optimize(function(mu) -ll(mu,Y[r,k],W[r,k]),c(1e-8,100))$objective,0)
        joint <- optimize(function(mu) -ll(mu,Y[r,c(i,j)],W[r,c(i,j)]),c(1e-8,100))$objective
        separate <- separate+sum(sep);shared <- shared+joint
      }
      expect_equal(as.matrix(a)[i,j],2*(shared-separate),tolerance=2e-6)
    }
    z <- glm_cluster(Y,family=f,k.max=3L,weights=rep(1,10),nthreads=2L)
    expect_length(z$clusterids,6)
    expect_true(all(is.finite(z$history$IC)))
  }
})

test_that("quasipoisson scores return scaled likelihood kernels", {
  X <- matrix(1,12,1);y <- rep(c(0,1,3,6),3)
  z <- L0glm.fit(X,y,family=quasipoisson("identity"),dispersion=3,
    lambda=0,nonnegative=TRUE,normalize="none")
  expect_equal(as.numeric(logLik(z)),sum(dpois(y,as.numeric(fitted(z)),log=TRUE))/3,tolerance=1e-8)
  expect_true(attr(logLik(z),"quasi"))
  expect_true(is.finite(GIC(z)))
})

test_that("signal families simulate both regression and blurred multichannel data", {
  for(f in list(nb1(3),Gamma_scale(2,censor=.5),Gamma_censored(.3))) {
    a <- simul("normal",n=60,p=5,n_nonzero=2,family=f,intercept=TRUE,m=2L)
    expect_equal(dim(a$y),c(60L,2L));expect_true(f$validate_response(a$y))
    f <- switch(f$scalar$kind-3L,Gamma_censored(.3,"identity"),nb1(3,"identity"),
      Gamma_scale(2,"identity",censor=.5))
    b <- simul("deconvolution",n=100,p=100,n_nonzero=5,m=3L,family=f)
    expect_equal(dim(b$y),c(100L,3L));expect_true(f$validate_response(b$y))
  }
})
