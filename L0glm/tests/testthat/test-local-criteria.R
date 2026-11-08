test_that('special-family paths score the same IC as their public methods', {
  skip_if_not_installed('mgcv')
  set.seed(8810)
  X <- matrix(rnorm(75*4),75,4)
  Ys <- list(matrix(rnbinom(150,mu=exp(.5*X[,1]),size=3),75,2),
    factor(rep(1:3,25)),survival::Surv(rexp(75,exp(.7*X[,1])),rep(1,75)))
  families <- list(mgcv::nb(theta=3),mgcv::multinom(K=2),mgcv::cox.ph())
  for (j in seq_along(families)) for (grouped in c(FALSE,TRUE)) {
    for (criterion in c('gic','bic','aic','ebic','mbic')) {
      tuned <- l0glmTune(X,Ys[[j]],family=families[[j]],tune=criterion,
        nlambda=3,prop_range=.1,grouped=grouped,parallelize=FALSE,
        plot=FALSE,printbest=FALSE,local.search='auto',
        control.local=list(maxit=1L,max.evals=8L,ncores=1L))
      expect_equal(tuned$metrics$tuning_score,tuned$metrics[[criterion]])
      expect_equal(tuned$lambda_best,tuned$metrics$lambda[which.min(tuned$metrics[[criterion]])])
      expect_equal(l0glmMetrics(tuned$fit_best,mu_floor=tuned$mu_floor)[[criterion]],
        min(tuned$metrics[[criterion]]),tolerance=1e-8)
    }
  }
})

test_that('parameter and row objectives remain explicit for multichannel search', {
  set.seed(8811)
  X <- matrix(rexp(120*7),120,7)
  B <- matrix(0,7,3);B[1:2,] <- c(3,1);B[2,2] <- 0
  y <- matrix(rpois(360,as.numeric(X%*%B)),120,3)
  a <- list(X=X,y=y,family=poisson('identity'),nonnegative=TRUE,parallelize=FALSE,
    lambda='gic',grouped=TRUE,local.search=TRUE)
  for (target in c('IC','L0')) {
    fit <- do.call(l0glmFit,c(a,list(control.local=list(target=target,maxit=2L,max.evals=18L))))
    met <- l0glmMetrics(fit,mu_floor=1e-4)
    df <- if(target=='IC')met$df else met$support_size_estimated
    expect_equal(fit$local.search$final$value,-2*met$loglik+fit$init.lambda*df,tolerance=1e-8)
  }
})

test_that('local-search CV scores can be reproduced by training-only fits', {
  set.seed(8812)
  X <- matrix(rnorm(54*5),54,5);y <- 2*X[,1]+rnorm(54)
  foldid <- rep(1:3,18)
  control <- list(maxit=1L,max.evals=8L,ncores=1L)
  tuned <- l0glmTune(X,y,family=gaussian(),tune='cv',nlambda=3,prop_range=.1,
    foldid=foldid,parallelize=FALSE,plot=FALSE,printbest=FALSE,
    local.search=TRUE,control.local=control)
  keep <- foldid!=1
  expected <- vapply(seq_along(tuned$lambdas),function(i) {
    fit <- l0glmFit(X[keep,],y[keep],lambda=tuned$cv$lambda_train[i,1],
      parallelize=FALSE,local.search=TRUE,
      control.local=c(control,list(mu.floor=tuned$mu_floor)))
    mean((y[!keep]-as.numeric(predict(fit,X[!keep,])))^2)
  },numeric(1))
  expect_equal(tuned$cv$deviance[,1],expected,tolerance=1e-9,ignore_attr=TRUE)
})
