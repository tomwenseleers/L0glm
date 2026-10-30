slcn_test_data <- function() {
  set.seed(12);n <- 180L
  X <- convmatrix(convkernel("gaussian",15L,list(sigma=1.5)),n)
  H <- cbind(X[,c(25,70,135)],pmax(0,1-abs(seq_len(n)-100)/15))
  H <- H%*%Matrix::Diagonal(x=1/Matrix::colSums(H))
  S <- matrix(rexp(24,rate=1/3000),4,6)
  B <- splines::bs(seq_len(n),df=4,intercept=TRUE)
  mu <- as.matrix(H%*%S)+B%*%matrix(20,4,6)
  list(X=X,Y=matrix(rpois(length(mu),mu),n,6),B=B)
}

slcn_check_model <- function(fit) {
  expect_s3_class(fit,"splogconNMF");expect_s4_class(fit$H,"dgCMatrix")
  expect_equal(unname(fitted(fit)),unname(as.matrix(fit$H%*%fit$S)+fit$background%*%fit$alpha+fit$offset),tolerance=1e-8)
  expect_true(all(fit$S>=0)&&all(fit$alpha>=0))
  expect_equal(as.numeric(Matrix::colSums(fit$H)),rep(1,ncol(fit$H)),tolerance=1e-10)
  for(j in seq_len(ncol(fit$H))) {
    h <- as.numeric(fit$H[,j]);rows <- which(h>0)
    expect_true(all(diff(rows)==1L))
    expect_true(all(diff(log(h[rows]),differences=2L)<=1e-6))
  }
  expect_true(is.finite(fit$loglik))
  expect_true(all(fit$history$delta[fit$history$accepted]<=1e-6))
  expect_equal(fit$objective, -2*fit$loglik+fit$lambda*fit$rank,tolerance=1e-8)
  expect_lt(abs(fit$ledger_error),1e-3+1e-8*abs(fit$complete_delta))
  expect_equal(as.numeric(logLik(fit)),fit$loglik)
  expect_equal(AIC(fit),-2*fit$loglik+2*fit$spec$information_factor*attr(logLik(fit),"df"))
}

test_that("all noise models preserve original dimensions, exact scores and constraints", {
  d <- slcn_test_data();Y <- cbind(d$Y,silent=0);W <- matrix(1,nrow(Y),ncol(Y));W[1:5,] <- 0
  controls <- list(broad_repeats=0L,full_sweeps=0L,final_polishes=0L)
  for(family in c("gaussian","poisson","quasipoisson","nb","Gamma","nb1","Gamma_scale")) {
    phi <- switch(family,gaussian=16,quasipoisson=5,Gamma=.3,nb1=5,Gamma_scale=5,NULL)
    cfg <- list(X=d$X,y=Y,family=family,background=d$B,weights=W,dispersion=phi,
      theta=if(family=="nb") 10 else NULL,offset=.1,correlation="none",control=controls)
    a <- do.call(splogconNMF,c(cfg,list(ncores=1L)))
    b <- do.call(splogconNMF,c(cfg,list(ncores=2L)))
    slcn_check_model(a);slcn_check_model(b)
    expect_equal(dim(a$mu),dim(Y));expect_equal(ncol(a$S),7L)
    expect_true(all(a$S[,7L]==0))
    expect_identical(a$S,b$S);expect_identical(a$H,b$H)
    expect_identical(a$mu,b$mu);expect_identical(a$alpha,b$alpha)
    reference <- switch(family,
      gaussian=dnorm(Y,a$mu,4,log=TRUE),poisson=dpois(Y,a$mu,log=TRUE),
      quasipoisson=dpois(Y,a$mu,log=TRUE)/5,nb=dnbinom(Y,mu=a$mu,size=10,log=TRUE),
      nb1=dnbinom(Y,mu=a$mu,size=a$mu/4,log=TRUE),
      Gamma=ifelse(Y==0,pgamma(.5,shape=1/.3,scale=a$mu*.3,log.p=TRUE),dgamma(Y,shape=1/.3,scale=a$mu*.3,log=TRUE)),
      Gamma_scale=ifelse(Y==0,pgamma(.5,shape=a$mu/5,scale=5,log.p=TRUE),dgamma(Y,shape=a$mu/5,scale=5,log=TRUE)))
    expect_equal(a$loglik,sum(W[W>0]*reference[W>0]),tolerance=1e-8)
    cfg$start <- list(H=a$H*2,S=a$S/2,alpha=a$alpha)
    z <- do.call(splogconNMF,c(cfg,list(ncores=1L)))
    slcn_check_model(z);expect_lte(z$objective,a$objective+1e-4)
  }
})

test_that("preclustering refits every original channel and supports both linkages", {
  d <- slcn_test_data()
  for(linkage in c("ward","bregman")) {
    z <- splogconNMF(d$X,d$Y,background=d$B,precluster=3L,correlation="none",ncores=2L,
      control=list(linkage=linkage))
    slcn_check_model(z)
    expect_equal(z$aggregation$channels,3L);expect_equal(ncol(z$S),6L)
    expect_length(z$aggregation$group,6L)
  }
})

test_that("empty models, numeric penalties and family constructors are supported", {
  X <- convmatrix(convkernel("gaussian",9L),25L)
  z <- splogconNMF(X,matrix(0,25,2),ncores=1L,correlation="none")
  slcn_check_model(z)
  expect_equal(dim(z$H),c(25L,0L));expect_equal(dim(z$S),c(0L,2L))
  expect_true(all(z$mu==0));expect_output(print(z),"Sparse log-concave NMF")
  d <- slcn_test_data()
  a <- splogconNMF(d$X,d$Y,background=d$B,correlation="none",ncores=1L)
  b <- splogconNMF(as.matrix(d$X),d$Y,background=d$B,lambda=a$lambda,correlation="none",ncores=1L)
  expect_identical(a$S,b$S);expect_identical(a$mu,b$mu)
  for(family in list(nb1(4,link="identity"),Gamma_scale(3,link="identity",censor=.5))) {
    fit <- splogconNMF(d$X,d$Y,background=d$B,family=family,correlation="none",ncores=1L,
      control=list(broad_repeats=0L,full_sweeps=0L,final_polishes=0L))
    expect_equal(fit$spec$dispersion,family$scalar$parameter)
  }
  expect_error(splogconNMF(X,matrix(.5,25,2)),"counts")
  expect_error(splogconNMF(X,matrix(0,24,2)),"nonnegative columns")
  expect_error(splogconNMF(X,rep(0,25),background=rep(-1,25)),"background")
  expect_error(splogconNMF(X,rep(0,25),control=list(unknown=1)),"Unknown")
})

test_that("temporal noise calibration is parallel and recovers marginal variance", {
  set.seed(943);Y <- matrix(rnorm(18000,sd=3),3000,6)
  a <- temporal_dispersion(Y,gaussian(),method="moment",ncores=1L)
  b <- temporal_dispersion(Y,gaussian(),method="moment",ncores=2L)
  expect_identical(a,b)
  expect_equal(unname(a$pooled[["dispersion"]]),9,tolerance=.05)
})

test_that("native log-concavity projection agrees with the concave envelope", {
  set.seed(717)
  project <- getFromNamespace("spolo_lc","L0glm")
  for(n in c(2L,5L,31L,100L)) for(i in 1:10) {
    h <- rexp(n);z <- project(h)
    expect_equal(sum(z$h),1,tolerance=1e-12)
    expect_true(all(diff(log(z$h),differences=2L)<=1e-12))
    normalized <- h/sum(h)
    envelope <- log(z$h*(1+z$gap))
    expect_true(all(envelope>=log(normalized)-1e-12))
  }
})
