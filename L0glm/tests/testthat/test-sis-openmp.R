test_that("OpenMP SIS is bitwise identical across all supported family/link contexts", {
  for (case in sis_thread_cases()) {
    context <- L0glm_sis_context(case$y,case$family,nrow(case$X),case$weights,case$offset,case$nuisance)
    for (X in list(case$X,Matrix::Matrix(case$X,sparse=TRUE))) {
      for (full in c(TRUE,FALSE)) {
        reference <- sis_score_cpp(X,context,7L,full,1L)
        for (threads in c(2L,4L,8L)) {
          actual <- sis_score_cpp(X,context,7L,full,threads)
          expect_identical(actual[names(actual)!="nthreads"],reference[names(reference)!="nthreads"])
          expect_gte(actual$nthreads,1L);expect_lte(actual$nthreads,threads)
        }
      }
    }
  }
})

test_that("top-k heap merges preserve ties across workers and interrupt batches", {
  set.seed(2711);column <- rnorm(30);X <- matrix(rep(column,20000),30,20000)
  context <- L0glm_sis_context(rnorm(30),gaussian(),30,NULL,NULL,NULL)
  for (A in list(X,Matrix::Matrix(X,sparse=TRUE))) for (full in c(TRUE,FALSE))
    for (threads in c(1L,2L,4L,8L)) {
      fit <- sis_score_cpp(A,context,11L,full,threads)
      expect_identical(fit$indices,seq_len(11L))
      expect_true(all(fit$top_statistics==fit$top_statistics[1L]))
    }
})

test_that("public SIS/ISIS threading preserves results and does not alter BLAS settings", {
  skip_if_not_installed("mgcv");skip_if_not_installed("survival")
  set.seed(2712);X <- matrix(rnorm(100*160),100,160);y <- X[,1]*2+rnorm(100)
  blas <- RhpcBLASctl::blas_get_num_procs();omp <- RhpcBLASctl::omp_get_max_threads()
  first <- sis_glm(X,y,nthreads=1L)
  for (threads in c(2L,4L,8L)) {
    fit <- sis_glm(X,y,nthreads=threads)
    expect_identical(fit[names(fit)!="nthreads"],first[names(first)!="nthreads"])
    expect_equal(RhpcBLASctl::blas_get_num_procs(),blas)
    expect_equal(RhpcBLASctl::omp_get_max_threads(),omp)
  }
  for (f in list(gaussian(),poisson(identity),binomial(),mgcv::nb(),mgcv::multinom(2),mgcv::cox.ph())) {
    kind <- L0glm_ext_kind(f)
    yy <- if (identical(kind,"multinomial")) factor(sample(letters[1:3],100,TRUE)) else
      if (identical(kind,"cox")) survival::Surv(rexp(100),rbinom(100,1,.7)) else
        switch(f$family,gaussian=y,poisson=rpois(100,2),binomial=rbinom(100,1,.5),rnbinom(100,mu=2,size=2))
    A <- if (f$family=="poisson") abs(X) else X
    included_intercept <- !identical(kind,"cox")
    reference <- L0glm.fit(A,yy,family=f,screen="ISIS",nsis=10,screen.maxit=2,screen.nthreads=1,
      intercept=included_intercept,nonnegative=f$family=="poisson")
    for (threads in c(2L,4L,8L)) {
      actual <- L0glm.fit(A,yy,family=f,screen="ISIS",nsis=10,screen.maxit=2,screen.nthreads=threads,
        intercept=included_intercept,nonnegative=f$family=="poisson")
      expect_identical(coef(actual),coef(reference));expect_identical(actual$screening$history,reference$screening$history)
      expect_identical(fitted(actual),fitted(reference))
      expect_gte(actual$screening$nthreads,1L);expect_lte(actual$screening$nthreads,threads)
    }
  }
  standalone <- sis_glm(X,y,nsis=10,iterative=TRUE,niter=2,nthreads=4)
  main <- L0glm.fit(X,y,screen="ISIS",nsis=10,screen.maxit=2,screen.nthreads=4)
  expect_identical(coef(standalone$fit),coef(main));expect_identical(standalone$history,main$screening$history)
  cv <- L0glm.tune(X,y,family=gaussian(),tune="cv",nlambda=3,nfolds=3,screen="SIS",nsis=10,
    screen.nthreads=4,plot=FALSE,printbest=FALSE)
  expect_true(is.finite(cv$lambda_best));expect_gte(cv$fit_best$screening$nthreads,1L)
})

test_that("threaded screening reports errors on the main thread and validates counts", {
  set.seed(2713);X <- matrix(rnorm(40*20),40,20);y <- rnorm(40)
  for (bad in list(0,-1,1.2,NA,Inf,"4",c(1,2),NULL,TRUE,1e15)) {
    expect_error(sis_glm(X,y,nthreads=bad),"nthreads must be a positive integer")
    expect_error(L0glm.fit(X,y,screen="SIS",screen.nthreads=bad),"nthreads must be a positive integer")
  }
  context <- L0glm_sis_context(y,gaussian(),40,NULL,NULL,NULL)
  for (value in c(NA_real_,NaN,Inf,-Inf)) {
    bad <- X;bad[2,19] <- value
    for (A in list(bad,Matrix::Matrix(bad,sparse=TRUE)))
      expect_error(sis_score_cpp(A,context,5L,FALSE,8L),"finite")
  }
  huge <- context;huge$score[] <- 1e300
  expect_error(sis_score_cpp(X,huge,5L,FALSE,8L),"Undefined screening statistic")
  expect_error(sis_score_cpp(X,context,5L,FALSE,0L),"nthreads must be positive")
  small <- sis_score_cpp(X[,1:2],context,2L,FALSE,8L)
  expect_lte(small$nthreads,2L)
})
