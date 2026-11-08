test_that("screened local search preserves its original-coordinate search", {
  set.seed(9711)
  X <- matrix(rnorm(160*100),160,100)
  colnames(X) <- paste0("x",seq_len(ncol(X)))
  y <- 1+2*X[,2]-X[,7]+rnorm(160,sd=.5)
  args <- list(X=X,y=y,intercept=TRUE,always.include=4L,screen="SIS",nsis=12L,
    screen.nthreads=1L,parallelize=FALSE,local.search=TRUE,
    control.local=list(maxit=3L,max.evals=45L))
  args$start <- c(0,rep(1,100));args$start[c(5,100)] <- 0
  fast <- do.call(l0glmFit,args)
  slow <- local({
    local_mocked_bindings(.package="L0glm",L0glm_local_workspace=function(...)NULL)
    do.call(l0glmFit,args)
  })
  expect_equal(coef(fast),coef(slow),tolerance=1e-12)
  expect_equal(fast$local.search$initial,slow$local.search$initial)
  expect_equal(fast$local.search$final,slow$local.search$final)
  expect_equal(fast$local.search$history,slow$local.search$history)
  expect_equal(fast$local.search$proposals,slow$local.search$proposals)
  expect_equal(fast$local.search$pool,slow$local.search$pool)
  expect_equal(fast$lambdas,slow$lambdas)
  expect_equal(fast$nuisance,slow$nuisance)
  expect_identical(fast$p,slow$p)
  expect_identical(fast$always.include,slow$always.include)
  expect_equal(GIC(fast),GIC(slow),tolerance=1e-10)
  expect_equal(predict(fast,X),predict(slow,X),tolerance=1e-12)
  expect_identical(dim(coef(fast)),c(101L,1L))
  expect_s4_class(coef(fast),"sparseMatrix")
  expect_lt(fast$local.search$workspace[["working.rows"]],nrow(coef(fast)))
  expect_equal(as.numeric(coef(fast)[100,]),0)
  expect_true(5L%in%fast$always.include)
  all <- do.call(l0glmFit,modifyList(args,list(control.local=list(maxit=1L,max.evals=8L,scope="all"))))
  expect_null(all$local.search$workspace)
})

test_that("compact workspaces preserve sparse channel masks and validation", {
  set.seed(9712)
  X <- Matrix::Matrix(matrix(rnorm(140*80),140,80),sparse=TRUE)
  colnames(X) <- paste0("x",1:80)
  B <- matrix(0,80,2);B[2:3,] <- c(2,-1)
  y <- as.matrix(X%*%B)+matrix(rnorm(280,sd=.4),140,2)
  start <- matrix(1,80,2);start[79,] <- 0;start[2,1] <- 0
  args <- list(X=X,y=y,start=start,screen="ISIS",nsis=12,screen.maxit=2,
    screen.nthreads=1,grouped=TRUE,parallelize=FALSE,local.search=TRUE,
    control.local=list(maxit=2L,max.evals=20L),
    validation=list(X=X[1:30,80:1],y=y[1:30,],starts="default"))
  fast <- do.call(l0glmFit,args)
  slow <- local({
    local_mocked_bindings(.package="L0glm",L0glm_local_workspace=function(...)NULL)
    do.call(l0glmFit,args)
  })
  expect_equal(coef(fast),coef(slow),tolerance=1e-10)
  expect_equal(fast$local.search$history,slow$local.search$history,tolerance=1e-10)
  expect_equal(fast$validation$deviance,slow$validation$deviance,tolerance=1e-10)
  expect_s4_class(fast$X,"sparseMatrix")
  expect_identical(dim(coef(fast)),c(80L,2L))
  expect_equal(as.numeric(coef(fast)[79,]),c(0,0))
  expect_equal(as.numeric(coef(fast)[2,1]),0)
})

test_that("special-family compact searches retain original IC dimensions", {
  skip_if_not_installed("mgcv")
  set.seed(9713)
  X <- matrix(rnorm(130*70),130,70)
  families <- list(mgcv::nb(theta=3),mgcv::multinom(K=2),mgcv::cox.ph())
  for (family in families) {
    kind <- L0glm_ext_kind(family)
    y <- switch(kind,nb=rnbinom(130,mu=exp(.6*X[,1]),size=3),
      multinomial=factor(rep(1:3,length.out=130)),
      cox=survival::Surv(rexp(130,exp(.7*X[,1])),rep(1,130)))
    args <- list(X=X,y=y,family=family,screen="SIS",nsis=10L,
      screen.nthreads=1L,parallelize=FALSE,local.search=TRUE,
      control.local=list(maxit=2L,max.evals=20L))
    fast <- do.call(l0glmFit,args)
    slow <- local({
      local_mocked_bindings(.package="L0glm",L0glm_local_workspace=function(...)NULL)
      do.call(l0glmFit,args)
    })
    expect_equal(coef(fast),coef(slow),tolerance=1e-10)
    expect_equal(fast$local.search$history,slow$local.search$history,tolerance=1e-10)
    expect_equal(fast$local.search$pool,slow$local.search$pool)
    expect_equal(GIC(fast,mu_floor=1e-4),GIC(slow,mu_floor=1e-4),tolerance=1e-9)
    expect_equal(l0glmMetrics(fast)$p,l0glmMetrics(slow)$p)
    expect_identical(dim(coef(fast)),dim(coef(slow)))
  }
})

test_that("an unchanged compact model retains its initial nuisance state", {
  set.seed(9714)
  X <- matrix(rexp(140*70),140,70)
  y <- rpois(140,3+2*X[,2])
  args <- list(X=X,y=y,family=poisson("identity"),nonnegative=TRUE,
    intercept=TRUE,always.include=4L,lambda=10000,screen="SIS",nsis=8L,
    screen.nthreads=1L,parallelize=FALSE,local.search=TRUE,
    control.local=list(maxit=1L,max.evals=8L,tol=1e-4))
  fast <- do.call(l0glmFit,args)
  slow <- local({
    local_mocked_bindings(.package="L0glm",L0glm_local_workspace=function(...)NULL)
    do.call(l0glmFit,args)
  })
  expect_equal(nrow(fast$local.search$history),1L)
  expect_equal(coef(fast),coef(slow))
  expect_equal(fast$lambdas,slow$lambdas)
  keys <- setdiff(names(slow$nuisance),"seconds")
  expect_equal(fast$nuisance[keys],slow$nuisance[keys])
})
