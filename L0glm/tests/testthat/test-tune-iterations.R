test_that("fitters and tuning accept one positive iteration limit", {
  X <- matrix(seq_len(30)/30,10L)
  y <- rep(c(1,2),5)
  for (bad in list(0,-1,NA_real_,Inf,2.5,numeric(),"x",c(1,2),c(4,4),2^32)) {
    expect_error(l0glmFit(X,y,maxit=bad),"single positive integer")
    expect_error(l0nnPois(X,y,maxit=bad),"single positive integer")
    expect_error(l0glmTune(X,y,family=gaussian(),maxit=bad,plot=FALSE),"single positive integer")
    expect_error(l0glmTune(X,y,family=gaussian(),maxit=bad,tune="cv",plot=FALSE),"single positive integer")
  }
})

test_that("scalar iteration limits retain ordinary lambda CV and replay", {
  set.seed(7101)
  X <- matrix(rnorm(48L*5L),48L)
  y <- 2*X[,1]-X[,3]+rnorm(48L)
  args <- list(X=X,y=y,family=gaussian(),tune="cv",foldid=rep(1:3,16L),maxit=4L,
    nlambda=3L,prop_range=.1,local.search=FALSE,parallelize=FALSE,
    plot=FALSE,printbest=FALSE)
  fit <- do.call(l0glmTune,args)
  expect_length(fit$lambdas,3L)
  expect_equal(dim(fit$cv$deviance),c(3L,3L))
  expect_false(any(c("maxit_best","maxit_1se","maxits") %in% names(fit)))
  expect_null(fit$cv$maxit)
  expect_null(fit$metrics$maxit)
  expect_identical(fit$fit_args$maxit,4L)
  for (which in c("best","1se")) {
    replay <- do.call(l0glmFit,c(list(X=X,y=y,
      lambda=fit[[paste0("lambda_",which)]]),fit$fit_args))
    expect_equal(coef(replay),coef(fit[[paste0("fit_",which)]]))
  }
  file <- tempfile(fileext=".pdf")
  grDevices::pdf(file)
  expect_silent(L0glm_plot_tuning(fit))
  grDevices::dev.off()
  expect_gt(file.info(file)$size,1000)
})

test_that("scalar limits retain special-family lambda CV", {
  skip_if_not_installed("mgcv")
  skip_if_not_installed("survival")
  set.seed(7102)
  X <- matrix(rnorm(48L*3L),48L)
  cases <- list(
    list(family=mgcv::nb(),y=rnbinom(48L,mu=exp(.2+.5*X[,1]),size=3)),
    list(family=mgcv::multinom(K=2),y=factor(rep(1:3,16))),
    list(family=mgcv::cox.ph(),y=survival::Surv(rexp(48L,exp(.5*X[,1])),rep(c(1,1,0),16)))
  )
  for (case in cases) {
    tuned <- do.call(l0glmTune,c(list(X=X,tune="cv",foldid=rep(1:3,each=16),
      maxit=5L,nlambda=2L,prop_range=.1,relaxed=FALSE,local.search=FALSE,
      plot=FALSE,printbest=FALSE),case))
    expect_equal(dim(tuned$cv$deviance),c(2L,3L))
    expect_true(all(is.finite(tuned$cv$mean)))
    expect_identical(tuned$fit_args$maxit,5L)
    replay <- do.call(l0glmFit,c(list(X=X,y=case$y,lambda=tuned$lambda_1se),tuned$fit_args))
    expect_equal(coef(replay),coef(tuned$fit_1se))
  }
})
