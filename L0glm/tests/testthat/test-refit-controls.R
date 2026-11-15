test_that("public fits retain penalized estimates by default and allow relaxed refits", {
  expect_identical(formals(l0glmFit)$relaxed,FALSE)
  set.seed(7110)
  X <- matrix(runif(90L*6L),90L)
  Y <- matrix(rpois(270,2+outer(X[,1],c(8,4,2))),90L)
  for(openmp in c(FALSE,TRUE)) {
    f <- l0glmFit(X,Y,family=poisson(identity),offset=2,nonnegative=TRUE,
      lambda="gic",control.solver=utils::modifyList(list(),list(.native=openmp)),local.search=FALSE,ncores=1L)
    expect_false(any(f$convergence$stage=="relaxed"))
    expect_true(all(coef(f)>=0))
    expect_false("relaxed"%in%names(f$specification))
  }
  selected_fit <- l0glmFit(X,Y,family=poisson(identity),offset=2,
    nonnegative=TRUE,lambda="gic",local.search=FALSE,ncores=1L)
  selected <- rowSums(as.matrix(coef(selected_fit)) != 0) > 0
  relaxed <- l0glmFit(X,Y,family=poisson(identity),offset=2,
    nonnegative=TRUE,lambda="gic",relaxed=TRUE,local.search=FALSE,ncores=1L)
  relaxed_selected <- rowSums(as.matrix(coef(relaxed)) != 0) > 0
  expect_true(any(relaxed$convergence$stage=="relaxed"))
  expect_true(all(which(relaxed_selected)%in%which(selected)))
  expect_true(all(coef(relaxed)>=0))
  expect_true(relaxed$specification$relaxed)
  expect_error(l0glmFit(X,Y,relaxed=NA),"relaxed must be TRUE or FALSE")
  expect_error(l0glmTune(X,Y,family=poisson(),relaxed=TRUE),
    "not supported by l0glmTune")
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
    f <- do.call(l0glmFit,c(list(X=X,group=group,lambda="gic",local.search=FALSE,ncores=1L,screen="none"),case))
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
  refitted <- l0glmFit(wide,y,screen="SIS",nsis=8L,relaxed=TRUE,
    local.search=FALSE)
  expect_true(any(refitted$convergence$stage=="relaxed"))
  expect_true(refitted$specification$relaxed)
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
    fit <- l0glmFit(X,Y,family=family,offset=2,start=start,lambda='gic',nonnegative=TRUE,ncores=1L,control.local=list(size.max=0L))
    search <- fit$local.search
    expect_identical(search$stopping,'support_refit_only')
    expect_equal(search$evaluations,1L)
    expect_lte(search$final$value,search$initial$value)
    expect_lte(nrow(search$history),2L)
    expect_true(all(coef(fit)[6,]==0))
    expect_equal(search$final$GIC,GIC(fit,mu_floor=search$scoring.floor),tolerance=1e-8)
  }
})
