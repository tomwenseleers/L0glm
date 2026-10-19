sis_reference <- function(X,y,f,w,off,N,one_sided=FALSE) {
  fit <- glm.fit(N,y,weights=w,offset=off,family=f,control=glm.control(epsilon=1e-10,maxit=100))
  mu <- fit$fitted.values;d <- f$mu.eta(fit$linear.predictors);v <- f$variance(mu)
  phi <- if(f$family%in%c("poisson","binomial")||grepl("negative binomial",f$family,ignore.case=TRUE)) 1 else
    sum(w*(y-mu)^2/v)/max(1,sum(w>0)-fit$rank)
  s <- w*(y-mu)*d/v/phi;h <- w*d^2/v/phi
  G <- solve(crossprod(N,h*N));U0 <- crossprod(N,s)
  vapply(seq_len(ncol(X)),function(j) {
    x <- X[,j];a <- crossprod(N,h*x)
    U <- sum(x*s)-sum(a*(G%*%U0));I <- sum(h*x^2)-sum(a*(G%*%a))
    if(I<=sum(h*x^2)*1e-12) 0 else (if(one_sided) max(U,0) else U)^2/I
  },numeric(1))
}

test_that("native SIS matches efficient scores across ordinary families and links", {
  set.seed(2501);n <- 180;p <- 12
  X <- cbind(matrix(rnorm(n*(p-1)),n,p-1),1)
  N <- cbind(1,runif(n,-.2,.2));w <- sample(0:3,n,TRUE);off <- runif(n,-.02,.02)
  families <- list(gaussian(),gaussian("log"),gaussian("inverse"),binomial(),
    binomial("probit"),binomial("cloglog"),binomial("cauchit"),binomial("identity"),
    binomial("log"),poisson(),poisson("identity"),poisson("sqrt"),Gamma("log"),
    Gamma("identity"),Gamma("inverse"),inverse.gaussian("log"),inverse.gaussian("identity"),
    inverse.gaussian("inverse"),inverse.gaussian("1/mu^2"),quasibinomial(),quasipoisson(),
    MASS::negative.binomial(2))
  for(f in families) {
    eta <- if(f$family%in%c("binomial","quasibinomial")) switch(f$link,identity=.4+off,log=-1+off,off) else
      if(f$link%in%c("identity","inverse","sqrt","1/mu^2")) 2+off else .7+off
    mu <- f$linkinv(eta)
    y <- if(f$family%in%c("binomial","quasibinomial")) rbinom(n,1,mu) else
      if(f$family%in%c("poisson","quasipoisson")) rpois(n,mu) else
        if(grepl("Negative Binomial",f$family)) rnbinom(n,mu=mu,size=2) else
          if(f$family=="gaussian") mu+rnorm(n,0,.1) else rgamma(n,40,scale=mu/40)
    expected <- sis_reference(X,y,f,w,off,N)
    actual <- sis_glm(X,y,f,weights=w,offset=off,nuisance=N,nsis=5)
    expect_equal(actual$stats$stat,expected,tolerance=1e-7,info=paste(f$family,f$link))
    expect_identical(actual$indices,head(order(expected,decreasing=TRUE),5L))
    sparse <- sis_glm(Matrix::Matrix(X,sparse=TRUE),y,f,weights=w,offset=off,nuisance=N,nsis=5)
    expect_equal(sparse$stats$stat,expected,tolerance=1e-7)
    expect_equal(actual$stats$stat[p],0)
  }
})

test_that("SIS weights, signs, nuisance projections and column scales are respected", {
  set.seed(2502);n <- 150;X <- matrix(rnorm(n*8),n,8);N <- cbind(1,runif(n))
  off <- sin(seq_len(n))*.3;w <- sample(0:3,n,TRUE);f <- binomial("probit")
  y <- rbinom(n,1,pnorm(off+.3*N[,2]))
  score <- sis_glm(X,y,f,weights=w,offset=off,nuisance=N,nsis=4)
  expect_equal(score$stats$stat,sis_reference(X,y,f,w,off,N),tolerance=1e-8)
  one <- sis_glm(X,y,f,weights=w,offset=off,nuisance=N,nsis=4,nonnegative=TRUE)
  expect_equal(one$stats$stat,sis_reference(X,y,f,w,off,N,TRUE),tolerance=1e-8)
  scaled <- sweep(X,2,10^seq(-12,12,length.out=8),"*")
  expect_equal(sis_glm(scaled,y,f,weights=w,offset=off,nuisance=N)$stats$stat,
    score$stats$stat,tolerance=1e-7)
  expect_equal(sis_glm(X+10,y,f,weights=w,offset=off,nuisance=N)$stats$stat,
    score$stats$stat,tolerance=1e-7)
  # Scaling or duplicating nuisance columns leaves their projection unchanged.
  NN <- cbind(N[,1],N[,2]*1e-10,N[,2]*1e10)
  expect_equal(sis_glm(X,y,f,weights=w,offset=off,nuisance=NN)$stats$stat,
    score$stats$stat,tolerance=1e-6)
  fit <- glm.fit(N,y,weights=w,offset=off,family=f)
  context <- L0glm_sis_context(y,f,n,w,off,N)
  ll <- function(b) sum(w*dbinom(y,1,pnorm(fit$linear.predictors+X[,1]*b),log=TRUE))
  expect_equal(sum(X[,1]*context$score),(ll(1e-5)-ll(-1e-5))/2e-5,tolerance=1e-7)
})

test_that("independent channels avoid cancellation and quadratic channel storage", {
  set.seed(2503);n <- 120;X <- Matrix::rsparsematrix(n,30,.08);y <- rnorm(n)
  a <- sis_glm(X,y,nsis=8);b <- sis_glm(X,cbind(y,-y),nsis=8)
  expect_equal(b$stats$stat,2*a$stats$stat,tolerance=1e-9)
  expect_identical(a$indices,b$indices);expect_true(all(is.na(b$stats$sign)))
  Y <- matrix(rnorm(n*80),n,80);w <- matrix(runif(n*80,.2,2),n,80)
  context <- L0glm_sis_context(Y,gaussian(),n,w,NULL,NULL)
  expect_identical(dim(context$weight),c(as.integer(n),80L))
  expect_length(context$inverse,80L)
  expect_lt(as.numeric(object.size(context)),1e6)
  actual <- sis_glm(X,Y,weights=w,nsis=8,return_stats=FALSE)
  expected <- Reduce(`+`,lapply(1:80,function(j) sis_glm(X,Y[,j],weights=w[,j])$stats$stat))
  expect_identical(actual$indices,head(order(expected,decreasing=TRUE),8L))
})

test_that("SIS handles fixed and estimated NB theta, including multichannel offsets", {
  skip_if_not_installed("mgcv")
  set.seed(2504);n <- 220;X <- matrix(rnorm(n*12),n,12);N <- cbind(1,runif(n,-.2,.2))
  off <- cbind(runif(n,-.1,.1),runif(n,-.2,.2));w <- cbind(sample(0:3,n,TRUE),sample(0:3,n,TRUE))
  y <- cbind(rnbinom(n,mu=exp(.7+off[,1]),size=2),rnbinom(n,mu=exp(.4+off[,2]),size=4))
  for(f in list(mgcv::nb(theta=2),mgcv::nb())) {
    score <- sis_glm(X,y,f,weights=w,offset=off,nuisance=N,nsis=5)
    expect_length(score$theta,2);expect_true(all(is.finite(score$theta)&score$theta>0))
    expected <- Reduce(`+`,lapply(1:2,function(j)
      sis_glm(X,y[,j],f,weights=w[,j],offset=off[,j],nuisance=N)$stats$stat))
    expect_equal(score$stats$stat,expected,tolerance=1e-6)
    expect_equal(sis_glm(Matrix::Matrix(X,sparse=TRUE),y,f,weights=w,offset=off,
      nuisance=N)$stats$stat,score$stats$stat,tolerance=1e-7)
  }
  # No nuisance terms: theta is still estimated with fixed known means.
  free <- sis_glm(X,y,mgcv::nb(),weights=w,offset=off,nuisance=matrix(0,n,0))
  expect_true(all(is.finite(free$theta)&free$theta>0))
})

test_that("joint multinomial SIS matches explicit coupled information and class permutations", {
  skip_if_not_installed("mgcv")
  set.seed(2505);n <- 100;X <- matrix(rnorm(n*6),n,6);N <- cbind(1,runif(n,-.2,.2))
  y <- matrix(rpois(n*3,2),n,3);w <- sample(0:3,n,TRUE);off <- matrix(rnorm(n*3,0,.1),n,3)
  f <- mgcv::multinom(2);context <- L0glm_sis_context(y,f,n,w,off,N)
  P <- context$probability;Q <- context$basis;q <- ncol(N);r <- ncol(Q)
  expected <- vapply(1:6,function(j) {
    Z <- cbind(N,X[,j]);H <- matrix(0,(q+1)*r,(q+1)*r)
    for(i in 1:n) H <- H+kronecker(crossprod(Q,(diag(P[i,])-tcrossprod(P[i,]))%*%Q)*w[i]*sum(y[i,]),tcrossprod(Z[i,]))
    U <- as.numeric(crossprod(Z,context$score));candidate <- seq_len(r)*(q+1)
    nuisance <- setdiff(seq_along(U),candidate);G <- solve(H[nuisance,nuisance])
    u <- U[candidate]-H[candidate,nuisance]%*%G%*%U[nuisance]
    I <- H[candidate,candidate]-H[candidate,nuisance]%*%G%*%H[nuisance,candidate]
    as.numeric(crossprod(u,solve(I,u)))
  },numeric(1))
  actual <- sis_glm(X,y,f,weights=w,offset=off,nuisance=N)
  expect_equal(actual$stats$stat,expected,tolerance=1e-7)
  expect_equal(sis_glm(Matrix::Matrix(X,sparse=TRUE),y,f,weights=w,offset=off,nuisance=N)$stats$stat,expected,tolerance=1e-7)
  expect_equal(sis_glm(X,y[,3:1],f,weights=w,offset=off[,3:1],nuisance=N)$stats$stat,expected,tolerance=1e-7)
  classes <- factor(sample(letters[1:3],n,TRUE))
  a <- sis_glm(X,classes,f);b <- sis_glm(X,factor(classes,levels=rev(levels(classes))),f)
  expect_equal(a$stats$stat,b$stats$stat,tolerance=1e-7)
  expect_equal(sis_glm(X,as.integer(classes)-1L,f)$stats$stat,a$stats$stat,tolerance=1e-7)
  expect_equal(sis_glm(X,classes,f,nonnegative=TRUE)$stats$stat,a$stats$stat)
  expect_lt(as.numeric(object.size(context)),1e5)
})

test_that("Cox SIS matches literal weighted tied risk-set calculations", {
  skip_if_not_installed("mgcv");skip_if_not_installed("survival")
  set.seed(2506);n <- 130;X <- cbind(matrix(rnorm(n*6),n,6),1)
  time <- sample(1:35,n,TRUE);event <- rbinom(n,1,.65);y <- survival::Surv(time,event)
  w <- sample(0:3,n,TRUE);off <- rnorm(n,0,.1);f <- mgcv::cox.ph()
  for(N in list(matrix(0,n,0),matrix(rnorm(n),n,1))) {
    eta <- off
    if(ncol(N)) {
      fit <- L0glm.fit(N,y,family=f,weights=w,offset=off,lambda=0,maxit=100,tol=1e-10)
      eta <- as.numeric(N%*%coef(fit))+off
    }
    expected <- vapply(1:6,function(j) {
      Z <- cbind(N,X[,j]);U <- colSums(w*event*Z);H <- matrix(0,ncol(Z),ncol(Z))
      for(t in unique(time[event==1])) {
        risk <- time>=t;e <- w[risk]*exp(eta[risk]);d <- sum(w[time==t]*event[time==t])
        M <- colSums(e*Z[risk,,drop=FALSE])/sum(e)
        U <- U-d*M;H <- H+d*(crossprod(Z[risk,,drop=FALSE],e*Z[risk,,drop=FALSE])/sum(e)-tcrossprod(M))
      }
      q <- ncol(N)
      if(q) {G <- solve(H[1:q,1:q,drop=FALSE]);u <- U[q+1]-H[q+1,1:q,drop=FALSE]%*%G%*%U[1:q]
        I <- H[q+1,q+1]-H[q+1,1:q,drop=FALSE]%*%G%*%H[1:q,q+1,drop=FALSE]
      } else {u <- U[1];I <- H[1,1]}
      as.numeric(u^2/I)
    },numeric(1))
    actual <- sis_glm(X,y,f,weights=w,offset=off,nuisance=N)
    expect_equal(actual$stats$stat[1:6],expected,tolerance=1e-6)
    expect_equal(actual$stats$stat[7],0)
    expect_equal(sis_glm(Matrix::Matrix(X,sparse=TRUE),y,f,weights=w,offset=off,nuisance=N)$stats$stat,
      actual$stats$stat,tolerance=1e-7)
  }
})

test_that("SIS returns stable top-k output, supports matrix types and validates inputs", {
  set.seed(2507);X <- matrix(rnorm(200),40,5);X <- cbind(X,X[,1],0,1);y <- rnorm(40)
  full <- sis_glm(X,y,nsis=100);small <- sis_glm(X,y,nsis=3,return_stats=FALSE)
  expect_identical(small$indices,head(full$ranking,3L));expect_null(small$stats);expect_null(small$ranking)
  expect_equal(small$top_statistics,full$stats$stat[small$indices])
  expect_lt(match(1L,full$ranking),match(6L,full$ranking))
  expect_length(full$indices,ncol(X));expect_length(sis_glm(X,y)$indices,8)
  expect_equal(sis_glm(Matrix::Matrix(X,sparse=FALSE),y)$stats$stat,full$stats$stat)
  expect_length(sis_glm(Matrix::Diagonal(40),y,return_stats=FALSE)$indices,floor(40/log(40)))
  expect_length(sis_glm(matrix(1:200,40),y)$indices,5)
  for(k in list(0,-1,1.5,NA,Inf,c(1,2),"3")) expect_error(sis_glm(X,y,nsis=k),"nsis")
  bad <- X;bad[1,1] <- NA_real_;expect_error(sis_glm(bad,y),"finite")
  expect_error(sis_glm(X,y[-1]),"one row")
  expect_error(sis_glm(X,y,weights=rep(0,40)),"weights")
  expect_error(sis_glm(X,y,offset=NA_real_),"offset")
  expect_error(sis_glm(X,y,nuisance=matrix(1,5,2)),"nuisance")
  expect_error(sis_glm(X,y,return_stats=NA),"TRUE or FALSE")
  expect_error(sis_glm(X,y,nonnegative=1),"TRUE or FALSE")
  expect_error(sis_glm(matrix("a",40,2),y),"numeric")
})
