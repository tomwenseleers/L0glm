test_that("independent tuning can recover features excluded at the preset", {
  X <- diag(6)
  y <- c(4, 1, .2, -.1, .1, .2)
  args <- list(X=X, y=y, family=gaussian(), dispersion=1,
               parallelize=FALSE, grouped=FALSE, relaxed=FALSE, maxit=100)
  tuned <- do.call(L0glm.tune, c(args, list(IC="aic", nlambda=5, plot=FALSE,
                                         printbest=FALSE)))
  preset <- do.call(L0glm.fit, c(args, list(lambda=tuned$lambda_preset)))
  expect_equal(as.numeric(preset$coefficients)[2], 0)
  expect_gt(abs(as.numeric(tuned$coefficients[[1]])[2]), 0)
  for (i in seq_along(tuned$lambdas)) {
    direct <- do.call(L0glm.fit, c(args, list(lambda=tuned$lambdas[i])))
    expect_equal(tuned$coefficients[[i]], direct$coefficients)
  }
  expect_true(tuned$endpoint_zero)
  expect_equal(sum(tail(tuned$coefficients,1)[[1]] != 0), 0)
  expect_equal(tail(tuned$lambdas,1), tuned$lambda_max)
  expect_equal(tuned$fit_best$coefficients,
               tuned$coefficients[[which.min(tuned$metrics$aic)]])
})

test_that("automatic scale does not cancel a large fitted penalty", {
  set.seed(41)
  X <- matrix(rnorm(100*8),100)
  y <- 3*X[,1] + rnorm(100)
  a <- L0glm.fit(X,y,family=gaussian(),lambda="gic",parallelize=FALSE)
  b <- L0glm.fit(X,y,family=gaussian(),lambda=1e10,parallelize=FALSE)
  expect_gt(sum(a$coefficients!=0),0)
  expect_equal(sum(b$coefficients!=0),0)
})

test_that("Poisson null likelihood is infinite rather than floored during tuning", {
  tuned <- L0glm.tune(diag(6), c(30,10,2,0,1,0), family=poisson(identity),
    nonnegative=TRUE, nlambda=4, maxit=100, plot=FALSE, printbest=FALSE,
    parallelize=FALSE, mu_floor=0)
  expect_true(tuned$endpoint_zero)
  expect_identical(tail(tuned$metrics$gic,1), Inf)
  expect_true(is.finite(tuned$lambda_best))
})

test_that("default Poisson tuning scores zero predictions in every channel", {
  X <- rbind(diag(6),rep(0,6))
  y <- cbind(c(30,10,2,0,1,0,1),c(0,8,0,2,0,1,2))
  tuned <- L0glm.tune(X,y,family=poisson(identity),nonnegative=TRUE,
    grouped=FALSE,parallelize=FALSE,nlambda=4,plot=FALSE,printbest=FALSE)
  expect_true(all(is.finite(tuned$metrics$gic)))
  expect_true(tuned$lambda_best %in% tuned$lambdas)
  expect_true(tuned$endpoint_zero)
  b <- tuned$fit_best$coefficients
  expected <- sum(dpois(y,lambda=pmax(as.matrix(X %*% b),1e-8),log=TRUE))
  expect_equal(tuned$metrics$loglik[which.min(tuned$metrics$gic)],expected)
  expect_equal(tuned$mu_floor,1e-8)
  expect_error(L0glm.tune(X,y,family=poisson(identity),nonnegative=TRUE,
    grouped=FALSE,parallelize=FALSE,nlambda=2,plot=FALSE,printbest=FALSE,
    mu_floor=0),"positive mu_floor")
})

test_that("endpoint refinement removes a conservative working-model bound", {
  testthat::local_mocked_bindings(maxlambda=function(...) 1e12)
  tuned <- L0glm.tune(diag(6),c(4,1,.2,-.1,.1,.2),family=gaussian(),
    dispersion=1,IC="aic",nlambda=3,maxit=100,parallelize=FALSE,
    plot=FALSE,printbest=FALSE)
  expect_lt(tuned$lambda_max,100)
  expect_true(tuned$endpoint_zero)
  expect_lte(tuned$lambda_max/tuned$lambda_nonzero,1.05)
  nonzero <- L0glm.fit(diag(6),c(4,1,.2,-.1,.1,.2),family=gaussian(),
    dispersion=1,lambda=tuned$lambda_nonzero,maxit=100,parallelize=FALSE)
  expect_gt(sum(nonzero$coefficients!=0),0)
  expect_equal(sum(tail(tuned$coefficients,1)[[1]]!=0),0)
})

test_that("weighted offset tuning scores the fitted Gaussian model", {
  set.seed(42)
  X <- matrix(rnorm(80*3),80)
  w <- seq(.5,2,length.out=80); off <- seq(3,6,length.out=80)
  y <- off + X[,1]*2 + rnorm(80)/sqrt(w)
  tuned <- L0glm.tune(X,y,family=gaussian(),dispersion=2,weights=w,offset=off,
    nlambda=3,prop_range=.2,parallelize=FALSE,plot=FALSE,printbest=FALSE)
  for(i in seq_along(tuned$lambdas)) {
    b <- tuned$coefficients[[i]]
    ll <- sum(dnorm(y,as.numeric(X %*% b)+off,sd=sqrt(2/w),log=TRUE))
    expect_equal(tuned$metrics$loglik[i],ll)
    expect_equal(tuned$metrics$gic[i],-2*ll+preset.lambda("gic",3,80)*sum(b!=0))
  }
})

test_that("Gaussian likelihood profiles scale once per channel and counts p correctly", {
  X <- cbind(1,seq(-1,1,length.out=30))
  y <- cbind(rep(c(-1,1),15),rep(c(-2,2),15))
  b <- matrix(0,2,2)
  met <- get.metrics(b,X,y,family=gaussian())
  expected <- sum(dnorm(y[,1],0,1,log=TRUE)) + sum(dnorm(y[,2],0,2,log=TRUE))
  expect_equal(met$p,2)
  expect_equal(met$loglik,expected)
  expect_equal(met$gic,-2*expected+preset.lambda("gic",2,30)*2)
})

test_that("scoring retains channel dimensions through nonidentity inverse links", {
  X <- cbind(1,seq(-1,1,length.out=12))
  b <- matrix(c(-.5,1,.3,-.7),2,2)
  prob <- plogis(X %*% b)
  w <- matrix(rep(c(10,20),each=12),12,2)
  y <- round(w*prob)/w
  met <- get.metrics(b,X,y,weights=w,family=binomial())
  expect_equal(met$loglik,sum(dbinom(round(w*y),size=w,prob=prob,log=TRUE)))
  expect_equal(met$gic,-2*met$loglik+preset.lambda("gic",2,12)*sum(b!=0))
  yy <- exp(X %*% b) + .2
  phi <- c(.5,.8)
  met <- get.metrics(b,X,yy,family=Gamma(log),dispersion=phi)
  expected <- sum(dgamma(yy[,1],shape=1/phi[1],scale=exp(X %*% b)[,1]*phi[1],log=TRUE)) +
    sum(dgamma(yy[,2],shape=1/phi[2],scale=exp(X %*% b)[,2]*phi[2],log=TRUE))
  expect_equal(met$loglik,expected)
})

test_that("the working-model threshold checks every channel at its own scale", {
  X <- diag(3)
  z <- cbind(c(8,0,0),c(0,2,0))
  phi <- c(100,1)
  bound <- maxlambda(X,matrix(1,3,2),z,dispersion=phi,coefs=z,family=gaussian())
  # Enumerate all subsets: the zero model wins at this threshold on an
  # orthogonal, known-scale least-squares problem, including the second channel.
  expect_equal(bound,4)
  masks <- as.matrix(expand.grid(rep(list(c(FALSE,TRUE)),6)))
  scores <- apply(masks,1,function(mask) sum(sweep((z-z*mask)^2,2,phi,"/"))+bound*sum(mask))
  expect_equal(scores[1],min(scores))
  expect_equal(maxlambda(cbind(X,0),matrix(1,3,2),z,dispersion=phi,
                         coefs=rbind(z,0),family=gaussian()),bound)
})

test_that("tuning handles an empty preset, intercepts, and invalid controls", {
  args <- list(X=matrix(seq(-1,1,length.out=30),ncol=1), y=rep(1,30),
               family=gaussian(), dispersion=1, parallelize=FALSE,
               plot=FALSE,printbest=FALSE,IC="aic",nlambda=2)
  a <- do.call(L0glm.tune,args)
  expect_true(a$endpoint_zero)
  b <- do.call(L0glm.tune,c(args,list(intercept=TRUE)))
  expect_equal(nrow(b$coefficients[[1]]),2)
  expect_equal(b$metrics$p,rep(2,2))
  expect_error(do.call(L0glm.tune,c(args,list(start=c(1)))),"independently")
  args$nlambda <- 1
  expect_error(do.call(L0glm.tune,args),"at least 2")
  args$nlambda <- 3; args$prop_range <- 0
  expect_error(do.call(L0glm.tune,args),"prop_range")
  args$prop_range <- 1; args$mu_floor <- -1
  expect_error(do.call(L0glm.tune,args),"mu_floor")
})
