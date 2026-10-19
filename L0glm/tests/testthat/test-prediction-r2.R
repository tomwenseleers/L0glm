test_that("Gaussian weighted prediction R-squared agrees with weighted lm", {
  set.seed(2420)
  x <- matrix(rnorm(160),80,2)
  y <- 2 + x %*% c(1,-.6) + rnorm(80,sd=.5)
  w <- seq(.1,2,length.out=80);w[1:3] <- 0
  for (intercept in c(FALSE,TRUE)) {
    reference <- if (intercept) lm(y~x,weights=w) else lm(y~0+x,weights=w)
    expected <- summary(reference)
    X <- model.matrix(reference)
    for (sparse in c(FALSE,TRUE)) {
      scores <- get.metrics(coef(reference),Matrix::Matrix(X,sparse=sparse),y,weights=w,dispersion=1)
      expect_equal(scores$R2_weighted,expected$r.squared)
      expect_equal(scores$adj_R2_weighted,expected$adj.r.squared)
    }
  }
})

test_that("prediction scores use fitted GLM working weights rather than initial weights", {
  set.seed(2421)
  x <- matrix(runif(120,-1,1),60,2)
  X <- cbind(Intercept=1,x)
  prior <- rep(c(1,2,3),20)
  for (f in list(poisson(identity),binomial("probit"),Gamma("log"),inverse.gaussian("log"))) {
    beta <- matrix(c(2,.4,-.3),3,1)
    mu <- f$linkinv(as.numeric(X %*% beta))
    y <- if (f$family == "binomial") rbinom(60,1,mu) else if (f$family == "poisson")
      rpois(60,mu) else mu*runif(60,.7,1.3)
    w <- prior*f$mu.eta(f$linkfun(mu))^2/f$variance(mu)
    rss <- sum(w*(y-mu)^2);tss <- sum(w*(y-weighted.mean(y,w))^2)
    scores <- get.metrics(beta,X,y,weights=prior,family=f,dispersion=1)
    expect_equal(scores$R2_weighted,1-rss/tss)
    expect_equal(scores$adj_R2_weighted,1-(rss/tss)*59/57)
  }
  # Poisson zero means use the fitting guard, independent of deviance floors.
  X <- matrix(c(0,0,1,2),4,1);y <- c(0,0,2,2)
  exact <- get.metrics(1,X,y,family=poisson(identity),mu_floor=0)
  floored <- get.metrics(1,X,y,family=poisson(identity),mu_floor=.1)
  expect_equal(exact$R2_weighted,5/6)
  expect_equal(exact[c("R2_weighted","adj_R2_weighted")],floored[c("R2_weighted","adj_R2_weighted")])
})

test_that("multichannel prediction scores pool within-channel sums and active parameters", {
  set.seed(2422)
  X <- cbind(Intercept=1,matrix(rnorm(120),60,2))
  beta <- rbind(c(2,-3),c(1,0),c(0,.5))
  mu <- X %*% beta
  y <- mu + cbind(rnorm(60,.1,.3),rnorm(60,-.2,1))
  w <- cbind(seq(.2,1,length.out=60),seq(2,.5,length.out=60));w[1:5,2] <- 0
  rss <- sum(w*(y-mu)^2)
  tss <- sum(vapply(1:2,function(j)sum(w[,j]*(y[,j]-weighted.mean(y[,j],w[,j]))^2),numeric(1)))
  scores <- get.metrics(beta,X,y,weights=w,dispersion=1)
  expect_equal(scores$R2_weighted,1-rss/tss)
  expect_equal(scores$adj_R2_weighted,1-(rss/tss)*(115-2)/(115-4))
  # A mandatory term is still a fitted parameter when its coefficient is zero.
  beta[1,2] <- 0
  scores <- get.metrics(beta,X,y,weights=w,dispersion=1)
  expect_equal((1-scores$adj_R2_weighted)/(1-scores$R2_weighted),(115-2)/(115-4))
  # A common rescaling of precision weights cannot change either score.
  scaled <- get.metrics(beta,X,y,weights=w*1e7,dispersion=1)
  expect_equal(scores[c("R2_weighted","adj_R2_weighted")],scaled[c("R2_weighted","adj_R2_weighted")])
})

test_that("prediction scores preserve negatives and handle undefined adjustments", {
  X <- cbind(Intercept=1,x=1:6)
  bad <- get.metrics(c(100,1),X,1:6,dispersion=1)
  expect_lt(bad$R2_weighted,0)
  expect_lt(bad$adj_R2_weighted,bad$R2_weighted)
  constant <- get.metrics(1,matrix(1,6,1),rep(1,6),dispersion=1)
  expect_true(is.na(constant$R2_weighted))
  saturated <- get.metrics(rep(1,4),diag(4),1:4,dispersion=1)
  expect_true(is.finite(saturated$R2_weighted))
  expect_true(is.na(saturated$adj_R2_weighted))
  empty <- get.metrics(c(1,1),X,1:6,weights=rep(0,6),dispersion=1)
  expect_true(is.na(empty$R2_weighted))
  base <- get.metrics(c(1,1),X,1:6,weights=c(0,rep(1,5)),dispersion=1)
  ignored <- get.metrics(c(1,1),X,c(1e200,2:6),weights=c(0,rep(1,5)),dispersion=1)
  expect_equal(base[c("R2_weighted","adj_R2_weighted")],ignored[c("R2_weighted","adj_R2_weighted")])
})

test_that("S3 prediction metrics use fitted NB theta and score new observations", {
  skip_if_not_installed("mgcv")
  set.seed(2423)
  X <- matrix(rnorm(160),80,2)
  y <- matrix(rnbinom(160,size=2,mu=exp(1 + X[,1]*.4)),80,2)
  fit <- L0glm.fit(X,y,family=mgcv::nb(theta=2),intercept=TRUE,
    lambda=2,grouped=FALSE,relaxed=FALSE)
  mu <- fitted(fit)
  w <- mu/(1+mu/2)
  rss <- sum(w*(y-mu)^2)
  tss <- sum(vapply(1:2,function(j)sum(w[,j]*(y[,j]-weighted.mean(y[,j],w[,j]))^2),numeric(1)))
  scores <- L0glm.metrics(fit)
  expect_equal(scores$R2_weighted,1-rss/tss)
  k <- sum(coef(fit) != 0) + sum(coef(fit)[1,,drop=FALSE] == 0)
  expect_equal(scores$adj_R2_weighted,1-(rss/tss)*158/(160-k))
  future <- L0glm.metrics(fit,newdata=X[1:12,],y=y[1:12,])
  expect_true(is.finite(future$adj_R2_weighted))
  expect_equal(future$R2_weighted,L0glm_prediction_metrics(y[1:12,],mu[1:12,],matrix(1,12,2),
    fit$families,coef(fit),1L,TRUE)$R2_weighted)
})
