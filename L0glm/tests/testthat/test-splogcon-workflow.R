test_that("aggregation keeps the original-channel rank penalty and avoids repeated selection", {
  set.seed(916);n <- 70L;m <- 8L
  X <- convmatrix(convkernel("gaussian",11L,list(sigma=1.5)),n)
  H <- X[,c(20L,47L)];H <- H %*% Matrix::Diagonal(x=1/Matrix::colSums(H))
  S <- matrix(rexp(2L*m,rate=1/2500),2L,m)
  Y <- matrix(rpois(n*m,as.matrix(H%*%S)+2),n,m)
  f <- splogconNMF(X,Y,background=matrix(1,n,1L),precluster=2L,ncores=1L,
    correlation="none",control=list(profile_cycles=0L,profile_fusion=FALSE))
  a <- f$aggregation$model
  expect_equal(a$objective,-2*a$loglik+f$lambda*ncol(a$H),tolerance=1e-10)
  expect_equal(f$lambda,m*preset.lambda("gic",p=n,n=n))
  expect_identical(f$history$stage,c("initial","family_refit"))
  expect_equal(ncol(f$S),m)
  expect_lt(abs(f$ledger_error),1e-4)
  expect_error(splogconNMF(X,Y,control=list(profile_cycles=-1)),"profile_cycles")
  expect_error(splogconNMF(X,Y,control=list(fusion_passes=.5)),"fusion_passes")
})

test_that("native profile updates respect hard concavity and the original likelihood", {
  set.seed(123);n <- 41L;h <- dnorm(seq(-3,3,length.out=n));s <- c(50,80,10)
  mu <- 2+h%o%s;Y <- matrix(rpois(length(mu),mu),n,3L)
  for(kind in 0:6) for(mode in 0:2) {
    y <- if(kind%in%c(4L,6L))pmax(Y,.5) else Y
    a <- slcn_log_profile(y,mu,h,s,kind,if(kind==2L)3 else .5,5,1,n,
      maxit=4L,q=if(mode==0L)12L else 32L,nthreads=1L,mode=mode)
    b <- slcn_log_profile(y,mu,h,s,kind,if(kind==2L)3 else .5,5,1,n,
      maxit=4L,q=if(mode==0L)12L else 32L,nthreads=2L,mode=mode)
    expect_identical(a,b)
    expect_true(all(is.finite(a$profile)&a$profile>0))
    expect_lt(max(diff(log(a$profile),differences=2L)),1e-6)
    expect_true(all(a$history<=1e-7))
  }
})
