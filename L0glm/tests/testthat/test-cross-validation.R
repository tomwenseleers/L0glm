test_that("Gaussian CV matches independent weighted held-out predictions", {
  set.seed(501)
  n <- 36L
  X <- matrix(rnorm(n*4),n,4)
  y <- cbind(1+2*X[,1]+rnorm(n),-1-X[,2]+rnorm(n,sd=2))
  w <- cbind(seq(.5,2,length.out=n),seq(2,.5,length.out=n))
  w[1,1] <- w[2,2] <- 0
  off <- cbind(seq(0,.5,length.out=n),seq(.5,0,length.out=n))
  foldid <- rep(1:3,length.out=n)
  args <- list(X=X,y=y,weights=w,offset=off,family=gaussian(),
    dispersion=c(2,-3),intercept=TRUE,grouped=FALSE,relaxed=TRUE,
    parallelize=FALSE,maxit=40)
  tuned <- do.call(L0glm.tune,c(args,list(tune="cv",foldid=foldid,
    nlambda=3,prop_range=.1,plot=FALSE,printbest=FALSE)))
  expected <- matrix(NA_real_,3,3)
  total <- sum(w)
  fold_weights <- vapply(1:3,function(k)sum(w[foldid==k,,drop=FALSE]),numeric(1))
  for(k in 1:3) {
    valid <- foldid==k
    a <- args
    a$X <- X[!valid,,drop=FALSE]; a$y <- y[!valid,,drop=FALSE]
    a$weights <- w[!valid,,drop=FALSE]; a$offset <- off[!valid,,drop=FALSE]
    for(i in seq_along(tuned$lambdas)) {
      fit <- do.call(L0glm.fit,c(a,list(lambda=tuned$lambdas[i]*sum(w[!valid,,drop=FALSE])/total)))
      prediction <- as.matrix(cbind(1,X[valid,,drop=FALSE]) %*% fit$coefficients)+off[valid,,drop=FALSE]
      expected[i,k] <- sum(w[valid,,drop=FALSE]*(y[valid,,drop=FALSE]-prediction)^2)/fold_weights[k]
    }
  }
  proportions <- fold_weights/sum(fold_weights)
  mean <- as.numeric(expected %*% proportions)
  se <- vapply(seq_len(3),function(i)
    sqrt(sum(proportions*(expected[i,]-mean[i])^2)/2),numeric(1))
  expect_equal(unname(tuned$cv$deviance),expected)
  expect_equal(tuned$cv$mean,mean)
  expect_equal(tuned$cv$se,se)
  expect_equal(tuned$metrics$cv_deviance,mean)
  expect_equal(tuned$metrics$cv_se,se)
  expect_equal(tuned$cv$lambda_train,outer(tuned$lambdas,1-fold_weights/total))
  expect_identical(tuned$cv$nfolds,3L)
  expect_equal(tuned$lambda_best,tuned$lambdas[which.min(mean)])
  expect_equal(tuned$lambda_preset,preset.lambda("gic",ncol(X),n))
  full <- do.call(L0glm.fit,c(args,list(lambda=tuned$lambda_best)))
  expect_equal(tuned$fit_best$coefficients,full$coefficients)
  expect_equal(tuned$fit_best$dispersion,full$dispersion)
  full_1se <- do.call(L0glm.fit,c(args,list(lambda=tuned$lambda_1se)))
  expect_equal(tuned$fit_1se,full_1se)
  expect_equal(tuned$fit_1se$coefficients,
               tuned$coefficients[[match(tuned$lambda_1se,tuned$lambdas)]])
})

test_that("CV rows and supplied free dispersion starts stay inside each training fold", {
  set.seed(502)
  n <- 24L
  X <- matrix(rnorm(n*3),n,3,dimnames=list(as.character(seq_len(n)),NULL))
  y <- X[,1]+rnorm(n)
  folds <- rep(1:3,length.out=n)
  seen <- list()
  original <- L0glm.fit
  testthat::local_mocked_bindings(L0glm.fit=function(X,y,...) {
    args <- list(...)
    if(nrow(X)<n) seen[[length(seen)+1L]] <<- list(rows=as.integer(rownames(X)),
      y=y,dispersion=args$dispersion,start=args$start)
    original(X,y,...)
  })
  tuned <- L0glm.tune(X,y,family=gaussian(),dispersion=-2,tune="cv",
    foldid=folds,nlambda=3,prop_range=.1,parallelize=FALSE,plot=FALSE,printbest=FALSE)
  expect_length(seen,9)
  for(k in 1:3) for(i in 1:3) {
    call <- seen[[(k-1)*3+i]]
    expect_identical(call$rows,which(folds!=k))
    expect_equal(as.numeric(call$y),as.numeric(y[folds!=k]))
    expect_identical(call$dispersion,-2)
    expect_null(call$start)
  }
  expect_true(all(is.finite(tuned$cv$mean)))
})

test_that("Poisson CV uses all held-out channels and the numerical scoring floor", {
  set.seed(503)
  X <- rbind(c(0,0,0),matrix(rexp(19*3),19,3))
  beta <- matrix(c(4,0,1,0,3,1),3,2)
  y <- matrix(rpois(40,as.numeric(X %*% beta)),20,2)
  y[1,] <- c(1,2)
  folds <- rep(1:5,4)
  args <- list(X=X,y=y,family=poisson(identity),nonnegative=TRUE,
               grouped=TRUE,group.maxit=2,relaxed=TRUE,parallelize=FALSE)
  tuned <- do.call(L0glm.tune,c(args,list(tune="CV",foldid=folds,nlambda=3,
    prop_range=.1,plot=FALSE,printbest=FALSE)))
  expect_identical(tuned$tune,"cv")
  expect_identical(tuned$cv$foldid,folds)
  expect_true(all(is.finite(tuned$cv$deviance)))
  for(k in 1:5) {
    valid <- folds==k
    fit <- L0glm.fit(X[!valid,,drop=FALSE],y[!valid,,drop=FALSE],
      family=poisson(identity),nonnegative=TRUE,grouped=TRUE,group.maxit=2,
      relaxed=TRUE,parallelize=FALSE,lambda=tuned$cv$lambda_train[1,k])
    mu <- pmax(as.matrix(X[valid,,drop=FALSE] %*% fit$coefficients),1e-8)
    expected <- sum(poisson()$dev.resids(as.numeric(y[valid,,drop=FALSE]),as.numeric(mu),1))/sum(valid)/2
    expect_equal(unname(tuned$cv$deviance[1,k]),expected)
  }
})

test_that("default CV folds are balanced, reproducible, and have the usual equal-weight SE", {
  set.seed(504)
  X <- matrix(rnorm(30*3),30,3)
  y <- 2*X[,1]+rnorm(30)
  args <- list(X=X,y=y,family=gaussian(),dispersion=1,tune="cv",
               nlambda=2,prop_range=.1,parallelize=FALSE,plot=FALSE,printbest=FALSE)
  set.seed(77)
  a <- do.call(L0glm.tune,args)
  set.seed(77)
  b <- do.call(L0glm.tune,args)
  expect_identical(a,b)
  expect_identical(a$cv$nfolds,5L)
  expect_equal(tabulate(a$cv$foldid),rep(6L,5))
  expect_equal(a$cv$mean,rowMeans(a$cv$deviance))
  expect_equal(a$cv$se,apply(a$cv$deviance,1,sd)/sqrt(5))
  c <- do.call(L0glm.tune,c(args,list(foldid=a$cv$foldid)))
  expect_identical(c$cv,a$cv)
})

test_that("CV supports binomial and quasi families without relying on IC scores", {
  set.seed(505)
  X <- matrix(rnorm(30*3),30,3)
  w <- rep(20,30)
  y <- round(w*plogis(-.7+X[,1]))/w
  bin <- L0glm.tune(X,y,weights=w,family=binomial(),tune="cv",nfolds=3,
    nlambda=2,prop_range=.1,parallelize=FALSE,plot=FALSE,printbest=FALSE)
  expect_true(all(is.finite(bin$cv$mean)))
  expect_equal(nrow(bin$fit_best$coefficients),ncol(X)+1)
  X <- cbind(1,abs(X))
  y <- rpois(30,as.numeric(X %*% c(2,1,0,0)))
  quasi <- L0glm.tune(X,y,family=quasipoisson(identity),dispersion=2,
    nonnegative=TRUE,tune="cv",nfolds=3,nlambda=2,prop_range=.1,
    parallelize=FALSE,plot=FALSE,printbest=FALSE)
  expect_true(all(is.na(quasi$metrics$gic)))
  expect_true(all(is.finite(quasi$cv$mean)))
  expect_false(is.null(quasi$fit_best))
})

test_that("CV validates fold controls and rejects the removed IC argument clearly", {
  X <- matrix(seq_len(30),10,3)
  args <- list(X=X,y=seq_len(10),family=gaussian(),dispersion=1,tune="cv",
               nlambda=2,parallelize=FALSE,plot=FALSE,printbest=FALSE)
  for(bad in list(1,1.5,11,NA_real_,Inf,"x",numeric())) {
    a <- args; a$nfolds <- bad
    expect_error(do.call(L0glm.tune,a),"nfolds")
  }
  for(bad in list(rep(NA,10),rep(0,10),rep(c(1,3),5),1:3)) {
    a <- args; a$foldid <- bad
    expect_error(do.call(L0glm.tune,a),"foldid")
  }
  a <- args; a$foldid <- rep(1:2,5); a$nfolds <- 3
  expect_error(do.call(L0glm.tune,a),"match")
  a <- args; a$foldid <- rep(1:2,5); a$weights <- rep(c(0,1),5)
  expect_error(do.call(L0glm.tune,a),"Every fold")
  a <- args; a$weights <- rep(-1,10)
  expect_error(do.call(L0glm.tune,a),"nonnegative")
  a <- args; a$offset <- matrix(0,3,3)
  expect_error(do.call(L0glm.tune,a),"offset")
  a <- args; a$tune <- "gic"; a$foldid <- rep(1:2,5)
  expect_error(do.call(L0glm.tune,a),"only used")
  a <- args; a$tune <- "unknown"
  expect_error(do.call(L0glm.tune,a),"tune must be one of")
  a <- args; a$IC <- "gic"
  expect_error(do.call(L0glm.tune,a),"IC was renamed to tune")
})

test_that("a one-variable CV model retains zero GIC preset with a positive grid reference", {
  X <- matrix(seq(-1,1,length.out=15),ncol=1)
  y <- 2*X[,1]+sin(seq_len(15))
  tuned <- L0glm.tune(X,y,family=gaussian(),tune="cv",nfolds=3,
    dispersion=1,nlambda=2,prop_range=.1,parallelize=FALSE,plot=FALSE,printbest=FALSE)
  expect_identical(tuned$lambda_preset,0)
  expect_identical(tuned$lambda_reference,2)
  expect_true(all(is.finite(tuned$cv$mean)))
})

test_that("the one-SE choice uses the minimum's SE and the largest eligible lambda", {
  means <- c(1.1, 1, 1.6, 1.2, 1.4)
  ses <- c(100, .25, 100, 100, 100)
  testthat::local_mocked_bindings(L0glm_cross_validate = function(...) {
    list(mean = means, se = ses)
  })
  args <- list(X = diag(6), y = c(4, 1, .2, -.1, .1, .2),
    family = gaussian(), dispersion = 1, nlambda = 5, nfolds = 2,
    parallelize = FALSE, plot = FALSE, printbest = FALSE)
  a <- do.call(L0glm.tune, c(args, list(tune = "cv")))
  # The cutoff is 1 + .25. The last candidate fails despite its large SE,
  # and the fourth qualifies despite the intervening third candidate failing.
  expect_identical(a$lambda_best, a$lambdas[2])
  expect_identical(a$lambda_1se, a$lambdas[4])
  expect_identical(a$fit_1se$coefficients, a$coefficients[[4]])

  means <- c(1.3, 1, 1.2, 1, Inf)
  ses <- c(100, 0, 100, 0, 0)
  b <- do.call(L0glm.tune, c(args, list(tune = "cv")))
  # With zero SE only tied minima qualify; choose the largest of those lambdas.
  expect_identical(b$lambda_best, b$lambdas[2])
  expect_identical(b$lambda_1se, b$lambdas[4])

  means <- c(1.3, 1, 1.2, 1.1, Inf)
  tied <- do.call(L0glm.tune, c(args, list(tune = "cv")))
  expect_identical(tied$lambda_1se, tied$lambda_best)
  expect_identical(tied$fit_1se, tied$fit_best)

  ic <- do.call(L0glm.tune, c(args, list(tune = "gic")))
  expect_null(ic$lambda_1se)
  expect_null(ic$fit_1se)
})
