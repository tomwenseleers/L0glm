test_that("public fitters retain penalized estimates and omit relaxed controls", {
  expect_false("relaxed"%in%names(formals(l0glmFit)))
  expect_false("relaxed"%in%names(formals(l0nnPois)))
  set.seed(7110)
  X <- matrix(runif(90L*6L),90L)
  Y <- matrix(rpois(270,2+outer(X[,1],c(8,4,2))),90L)
  for(openmp in c(FALSE,TRUE)) {
    f <- l0glmFit(X,Y,family=poisson(identity),offset=2,nonnegative=TRUE,
      lambda="gic",openmp=openmp,local.search=FALSE,parallelize=FALSE)
    expect_false(any(f$convergence$stage=="relaxed"))
    expect_true(all(coef(f)>=0))
    expect_false("relaxed"%in%names(f$specification))
  }
  expect_error(l0glmFit(X,Y,relaxed=TRUE),"unused argument")
  expect_error(l0nnPois(X,Y,relaxed=TRUE),"unused argument")
  expect_error(l0glmTune(X,Y,family=poisson(),relaxed=TRUE),"not a fitting argument")
})

test_that("screening, special families and predictor groups retain penalties", {
  skip_if_not_installed("mgcv")
  skip_if_not_installed("survival")
  set.seed(7111)
  X <- matrix(rnorm(90L*6L),90L)
  cases <- list(
    list(y=cbind(X[,1]+rnorm(90),-X[,2]+rnorm(90)),family=gaussian()),
    list(y=cbind(rnbinom(90,mu=exp(.5*X[,1]),size=3),rnbinom(90,mu=exp(.5*X[,2]),size=3)),family=mgcv::nb()),
    list(y=factor(rep(1:3,30)),family=mgcv::multinom(2)),
    list(y=survival::Surv(rexp(90,exp(.5*X[,1])),rep(c(1,1,0),30)),family=mgcv::cox.ph()))
  for(case in cases)for(group in list(NULL,c(1,1,2,2,3,3))) {
    f <- do.call(l0glmFit,c(list(X=X,group=group,lambda="gic",local.search=FALSE,
      parallelize=FALSE,screen="none"),case))
    expect_false(any(f$convergence$stage=="relaxed"))
    expect_true(any(as.numeric(f$lambdas)>0)|is.null(f$lambdas))
  }
  wide <- matrix(rnorm(60L*800L),60L)
  y <- 2*wide[,1]-wide[,2]+rnorm(60)
  for(screen in c("SIS","ISIS")) {
    f <- l0glmFit(wide,y,screen=screen,nsis=8L,screen.maxit=2L,local.search=FALSE)
    expect_false(any(f$convergence$stage=="relaxed"))
    expect_equal(dim(coef(f)),c(800L,1L))
  }
})

test_that("explicitly unpenalized fits use a single ordinary fitting pass", {
  set.seed(7112)
  X <- matrix(rnorm(80L*4L),80L)
  y <- 1+2*X[,1]+rnorm(80)
  f <- l0glmFit(X,y,intercept=TRUE,lambda=0,local.search=FALSE)
  expected <- stats::lm.fit(cbind(1,X),y)$coefficients
  expect_equal(as.numeric(coef(f)),as.numeric(expected),tolerance=1e-8)
  expect_false(any(f$convergence$stage=="relaxed"))
})

test_that('large multichannel models get one criterion-checked support refit', {
  skip_if_not_installed('mgcv')
  set.seed(7113)
  X <- matrix(runif(100*6),100,6)
  B <- matrix(0,6,3);B[c(1,3),] <- c(4,7)
  for(family in list(poisson(identity),gaussian(),mgcv::nb(link='identity'))) {
    mu <- 2+X%*%B
    Y <- if(family$family=='gaussian')mu+matrix(rnorm(300),100)else
      if(family$family=='poisson')matrix(rpois(300,mu),100)else matrix(rnbinom(300,mu=mu,size=5),100)
    start <- matrix(1,6,3);start[6,] <- 0
    fit <- l0glmFit(X,Y,family=family,offset=2,start=start,lambda='gic',nonnegative=TRUE,
      parallelize=FALSE,control.local=list(size.max=0L))
    search <- fit$local.search
    expect_identical(search$stopping,'support_refit_only')
    expect_equal(search$evaluations,1L)
    expect_lte(search$final$value,search$initial$value)
    expect_lte(nrow(search$history),2L)
    expect_true(all(coef(fit)[6,]==0))
    expect_equal(search$final$GIC,GIC(fit,mu_floor=search$scoring.floor),tolerance=1e-8)
  }
})
