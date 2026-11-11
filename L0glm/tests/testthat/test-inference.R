test_that('selected-support Gaussian estimates and errors agree with weighted lm', {
  set.seed(87);n <- 160L;X <- matrix(rnorm(n*6),n,6)
  Y <- cbind(1+2*X[,1]-2*X[,2]+rnorm(n),2-X[,1]+X[,2]+rnorm(n))
  w <- cbind(runif(n,.5,1.5),runif(n,.5,1.5));offset <- matrix(rnorm(n*2,0,.2),n,2)
  start <- matrix(1,7,2);start[2,1] <- 0
  fit <- l0glmFit(X,Y,weights=w,offset=offset,intercept=TRUE,lambda='gic',
    start=start,local.search=FALSE)
  point <- coef(fit)
  inference <- summary(fit,inference='conditional')$inference
  for(j in 1:2) {
    table <- inference$coefficients[inference$coefficients$column==j,]
    A <- as.matrix(fit$X[,table$row,drop=FALSE])
    reference <- lm.wfit(A,Y[,j]-offset[,j],w[,j])
    V <- solve(crossprod(A,w[,j]*A))*sum(w[,j]*reference$residuals^2)/reference$df.residual
    expect_equal(table$estimate,as.numeric(reference$coefficients),tolerance=1e-7)
    expect_equal(table$std.error,as.numeric(sqrt(diag(V))),tolerance=1e-7)
  }
  expect_identical(coef(fit),point)
  expect_true(all(nzchar(inference$coefficients$term)))
  repeated <- do.call(l0glmFit,c(list(X=X,y=Y),fit$specification))
  expect_identical(coef(repeated),point)
})

test_that('multinomial CLR and Cox errors agree with reference likelihood fits', {
  skip_if_not_installed('mgcv');skip_if_not_installed('nnet');skip_if_not_installed('survival')
  set.seed(325);n <- 300L;X <- matrix(rnorm(n*4),n,4)
  Y <- t(vapply(seq_len(n),function(i)as.numeric(rmultinom(1,10,
    exp(c(.7*X[i,1],-.7*X[i,1],.7*X[i,2])))),numeric(3)))
  fit <- l0glmFit(X,Y,family=mgcv::multinom(2),lambda='gic',intercept=TRUE,
    local.search=FALSE,maxit=100L)
  details <- summary(fit,inference='conditional')$inference
  A <- fit$X[,details$covariance.rows,drop=FALSE];p <- ncol(A);C <- ncol(Y)
  reference <- nnet::multinom(Y~A-1,trace=FALSE,Hess=TRUE,maxit=500,abstol=1e-10,reltol=1e-12)
  beta <- t(rbind(0,coef(reference)));beta <- beta-rowMeans(beta)
  T <- kronecker(rbind(rep(0,C-1L),diag(C-1L))-matrix(1/C,C,C-1L),diag(p))
  expect_equal(as.numeric(details$refit[details$covariance.rows,]),as.numeric(beta),tolerance=1e-5)
  expect_equal(details$coefficients$std.error,sqrt(diag(T%*%vcov(reference)%*%t(T))),tolerance=1e-5)
  y <- survival::Surv(rexp(n,exp(.7*X[,1]-.7*X[,2])),rbinom(n,1,.8))
  w <- runif(n,.8,1.2);off <- rnorm(n,0,.1)
  fit <- l0glmFit(X,y,family=mgcv::cox.ph(),lambda='gic',weights=w,offset=off,
    local.search=FALSE,maxit=100L)
  details <- summary(fit,inference='conditional')$inference
  A <- fit$X[,details$covariance.rows,drop=FALSE]
  reference <- survival::coxph(y~A+offset(off),weights=w,ties='breslow',robust=FALSE,
    control=survival::coxph.control(eps=1e-10,iter.max=100))
  expect_equal(details$coefficients$estimate,as.numeric(coef(reference)),tolerance=1e-5)
  expect_equal(details$coefficients$std.error,as.numeric(sqrt(diag(vcov(reference)))),tolerance=1e-5)
})

test_that('bootstrap reselection includes zeros and restores a supplied seed', {
  set.seed(18);n <- 120L;X <- matrix(rnorm(n*6),n,6)
  y <- 1+2*X[,1]-2*X[,2]+X[,1]*X[,2]+rnorm(n)
  fit <- l0glmFit(X,y,intercept=TRUE,interactions=1L,lambda='gic',
    local.search=FALSE)
  seed <- .Random.seed
  control <- list(nboot=8L,seed=27)
  details <- summary(fit,inference='bootstrap',control.inference=control)$inference
  expect_identical(.Random.seed,seed)
  expect_equal(details$successful,8L)
  expect_s4_class(details$draws,'sparseMatrix')
  expect_true(all(details$coefficients$selection.frequency>=0 & details$coefficients$selection.frequency<=1))
  repeated <- summary(fit,inference='bootstrap',control.inference=control)$inference
  expect_equal(repeated$draws,details$draws)
  expect_true(all(details$coefficients$conf.low<=details$coefficients$conf.high))
  for(i in seq_len(nrow(details$coefficients))) {
    row <- details$coefficients$row[i]+(details$coefficients$column[i]-1L)*nrow(details$selection.frequency)
    expect_equal(details$coefficients$selection.frequency[i],mean(details$draws[row,]!=0))
  }
})

test_that('bootstrap failures, empty supports and memory guards are explicit', {
  set.seed(70);X <- matrix(rnorm(80*4),80,4);y <- 1+2*X[,1]+rnorm(80)
  fit <- l0glmFit(X,y,lambda='gic',local.search=FALSE)
  b <- 0L
  refit <- function(args) {b <<- b+1L;if(b==2L)stop('failed replicate');do.call(l0glmFit,args)}
  expect_warning(details <- summary(fit,inference='bootstrap',
    control.inference=list(nboot=4L,seed=7,refit=refit))$inference,'1 of 4')
  expect_equal(details$successful,3L)
  expect_equal(details$diagnostics$status[2],'failed replicate')
  expect_equal(ncol(details$draws),3L)
  expect_error(summary(fit,inference='bootstrap',control.inference=list(max.cells=1)),'max.cells')
  expect_error(summary(fit,inference='bootstrap',control.inference=list(seed=NA)),'seed')
  empty <- l0glmFit(X,y,lambda=1e8,local.search=FALSE)
  expect_equal(nrow(summary(empty,inference='conditional')$inference$coefficients),0L)
  boot <- summary(empty,inference='bootstrap',control.inference=list(nboot=3L,seed=1))$inference
  expect_equal(nrow(boot$coefficients),0L)
  expect_true(all(boot$selection.frequency==0))
})
