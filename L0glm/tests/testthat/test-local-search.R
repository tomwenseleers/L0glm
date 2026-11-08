local_budget <- list(maxit=2L,max.evals=20L,ncores=1L,nthreads=1L)

test_that('local support refinement is opt-in, monotone and reproducible', {
  set.seed(8801)
  X <- matrix(rnorm(180*16),180,16)
  y <- 2*X[,1]-1.5*X[,2]+rnorm(180)
  a <- list(X=X,y=y,lambda='gic',intercept=TRUE,relaxed=TRUE,parallelize=FALSE)
  plain <- do.call(l0glmFit,c(a,list(local.search=FALSE)))
  off <- do.call(l0glmFit,c(a,list(local.search=FALSE)))
  fit <- do.call(l0glmFit,c(a,list(local.search=TRUE,control.local=local_budget)))
  expect_identical(coef(off),coef(plain))
  expect_null(off$local.search)
  expect_identical(fit$local.search$initial.fit$coefficients,coef(plain))
  expect_true(fit$local.search$applied)
  expect_lte(fit$local.search$final$value,fit$local.search$initial$value)
  expect_true(all(diff(fit$local.search$history$value)<0))
  expect_equal(GIC(fit,mu_floor=1e-4),fit$local.search$final$GIC,tolerance=1e-8)
  expect_identical(dim(coef(fit)),c(17L,1L))
  expect_s4_class(coef(fit),'sparseMatrix')
  expect_equal(as.numeric(predict(fit,X)),as.numeric(fitted(fit)),tolerance=1e-10)
  expect_error(l0glmFit(X,y,local.search=NA),'local.search')
  expect_error(l0glmFit(X,y,local.search=TRUE,control.local=list(unknown=1)),'control.local')
})

test_that('starting exclusions, mandatory terms and sparse designs survive search', {
  set.seed(8802)
  X <- Matrix::rsparsematrix(160,18,.35)
  B <- matrix(0,18,3); B[1:3,] <- c(2,-1.5,1)
  y <- as.matrix(X%*%B)+matrix(rnorm(480,sd=.4),160,3)
  start <- matrix(1,18,3);start[3,1] <- start[18,] <- start[1,] <- 0
  fit <- l0glmFit(X,y,start=start,always.include=1,grouped=TRUE,
    parallelize=FALSE,local.search=TRUE,control.local=local_budget)
  expect_equal(as.numeric(coef(fit)[18,]),rep(0,3))
  expect_equal(as.numeric(coef(fit)[3,1]),0)
  expect_true(any(coef(fit)[1,]!=0))
  expect_identical(dim(coef(fit)),c(18L,3L))
  expect_s4_class(fit$X,'sparseMatrix')
  expect_equal(predict(fit,X),fitted(fit),tolerance=1e-8)
})

test_that('search runs after screening and original-channel preclustering refit', {
  set.seed(8803)
  X <- matrix(rnorm(100*70),100,70)
  y <- 2*X[,1]-X[,2]+rnorm(100)
  for (screen in c('SIS','ISIS')) {
    fit <- l0glmFit(X,y,screen=screen,nsis=12,screen.maxit=2,parallelize=FALSE,
      screen.nthreads=1,local.search=TRUE,control.local=local_budget)
    expect_true(fit$local.search$after.screening)
    expect_equal(fit$screening$method,screen)
    expect_identical(dim(coef(fit)),c(70L,1L))
    expect_equal(fit$local.search$scope,'screened')
  }
  X <- Matrix::Diagonal(70)
  y <- matrix(0,70,6);y[c(10,30,50),] <- rep(c(20,30,40),6)
  fit <- l0nnPois(X,y,parallelize=FALSE,precluster=TRUE,
    control.precluster=list(k.max=2),local.search=TRUE,control.local=local_budget)
  expect_true(fit$preclustering$refitted)
  expect_equal(ncol(coef(fit)),6)
  expect_lte(fit$local.search$final$value,fit$local.search$initial$value)
})

test_that('Poisson entry points share refinement and validation selection', {
  set.seed(8804)
  X <- matrix(rexp(120*9),120,9)
  B <- matrix(0,9,2);B[1:2,] <- c(4,2)
  y <- matrix(rpois(240,as.numeric(X%*%B)),120,2)
  a <- list(X=X,y=y,family=poisson('identity'),parallelize=FALSE,
    local.search=TRUE,control.local=local_budget)
  direct <- do.call(l0nnPois,a)
  general <- do.call(l0glmFit,c(a,list(nonnegative=TRUE)))
  expect_equal(coef(direct),coef(general),tolerance=1e-8)
  validation <- list(X=X[1:40,],y=y[1:40,])
  val <- do.call(l0nnPois,c(a,list(validation=validation)))
  expect_identical(val$local.search$target,'validation')
  expect_lte(val$validation$local.search[['final']],val$validation$local.search[['initial']])
})

test_that('ordinary and special-family tuning refines each training model', {
  skip_if_not_installed('mgcv')
  set.seed(8805)
  X <- matrix(rnorm(80*5),80,5)
  for (family in list(gaussian(),mgcv::nb(theta=3),mgcv::multinom(K=2),mgcv::cox.ph())) {
    kind <- L0glm_ext_kind(family)
    y <- if(is.null(kind))2*X[,1]+rnorm(80) else if(kind=='nb')
      rnbinom(80,mu=exp(.4*X[,1]),size=3) else if(kind=='multinomial')
      factor(rep(1:3,length.out=80)) else survival::Surv(rexp(80,exp(.6*X[,1])),rep(1,80))
    fit <- l0glmTune(X,y,family=family,nlambda=3,prop_range=.15,
      tune='cv',nfolds=2,foldid=rep(1:2,40),plot=FALSE,printbest=FALSE,
      parallelize=FALSE,local.search=TRUE,control.local=local_budget)
    expect_true(fit$fit_best$local.search$applied)
    expect_true(fit$fit_1se$local.search$applied)
    expect_equal(fit$fit_best$local.search$scoring.floor,fit$mu_floor)
    expect_true(all(is.finite(fit$cv$deviance)))
    replay <- do.call(l0glmFit,c(list(X=X,y=y,lambda=fit$lambda_best),fit$fit_args))
    expect_equal(coef(replay),coef(fit$fit_best),tolerance=1e-8)
  }
})

test_that('local search is disabled for pilot-only initialization fits', {
  set.seed(8806)
  X <- matrix(rnorm(60*4),60,4);y <- X[,1]+rnorm(60)
  calls <- 0L
  local_mocked_bindings(.package='L0glm',L0glm_local_refine=function(fit,args,control) {
    calls <<- calls+1L
    stop('unexpected pilot refinement')
  })
  expect_silent(pilot <- L0glm_pilot_start(list(X=X,y=y,parallelize=FALSE,
    local.search=TRUE,control.local=local_budget),c(2,5)))
  expect_equal(calls,0L)
  expect_equal(nrow(pilot$pilot),10)
})

test_that('automatic search uses independent coefficient dimensions', {
  for (fun in list(l0glmFit,l0nnPois,l0glmTune))
    expect_identical(formals(fun)$local.search,'auto')
  set.seed(8807)
  X <- matrix(rnorm(120*8),120,8);y <- 2*X[,1]+rnorm(120)
  a <- list(X=X,y=y,parallelize=FALSE)
  plain <- do.call(l0glmFit,c(a,list(local.search=FALSE)))
  automatic <- do.call(l0glmFit,c(a,list(local.search='auto',control.local=local_budget)))
  on <- do.call(l0glmFit,c(a,list(local.search=TRUE,control.local=local_budget)))
  expect_true(automatic$local.search$applied)
  expect_identical(coef(automatic),coef(on))
  off <- do.call(l0glmFit,c(a,list(local.search='auto',
    control.local=c(local_budget,list(size.max=0L)))))
  expect_false(off$local.search$applied)
  expect_identical(coef(off),coef(plain))
  expect_identical(off$local.search$stopping,'size_limit')
  expect_output(print(off),'skipped')
  fake <- plain;fake$coefficients <- Matrix::Matrix(matrix(c(1,0,1,0,0,1),2,3),sparse=TRUE)
  fake$always.include <- integer();fake$grouped <- TRUE
  expect_equal(L0glm_local_size(fake),6)
  fake$kind <- 'multinomial'
  expect_equal(L0glm_local_size(fake),4)
  fake$grouped <- FALSE
  expect_equal(L0glm_local_size(fake),1)
  expect_error(l0glmFit(X,y,local.search='auto',control.local=list(size.max=-1)),'size.max')
  default <- do.call(l0glmFit,c(a,list(control.local=local_budget)))
  expect_true(default$local.search$applied)
  expect_identical(coef(default),coef(automatic))
  pilot <- l0glmFit(X,y,maxit=1L,parallelize=FALSE)
  expect_identical(pilot$local.search$stopping,'initialization_only')
  unpenalized <- l0glmFit(X,y,lambda=0,parallelize=FALSE)
  expect_identical(unpenalized$local.search$stopping,'unpenalized')
  tune <- l0glmTune(X,y,family=gaussian(),nlambda=3,prop_range=.1,
    local.search='auto',control.local=local_budget,parallelize=FALSE,plot=FALSE,printbest=FALSE)
  expect_true(tune$fit_best$local.search$applied)
  disabled <- l0glmTune(X,y,family=gaussian(),nlambda=3,prop_range=.1,
    local.search=FALSE,parallelize=FALSE,plot=FALSE,printbest=FALSE)
  expect_null(disabled$fit_best$local.search)
})

test_that('default refinement executes once after each complete orchestrated fit', {
  set.seed(8813)
  original <- L0glm_local_refine
  calls <- list()
  local_mocked_bindings(.package='L0glm',L0glm_local_refine=function(fit,args,control) {
    calls[[length(calls)+1L]] <<- list(screening=fit$screening,clustering=fit$preclustering)
    original(fit,args,control)
  })
  X <- matrix(rnorm(80*40),80,40);y <- 2*X[,1]+rnorm(80)
  screened <- l0glmFit(X,y,screen='ISIS',nsis=8,screen.maxit=2,
    parallelize=FALSE,control.local=local_budget)
  expect_length(calls,1L)
  expect_equal(calls[[1]]$screening$method,'ISIS')
  X <- Matrix::Diagonal(60);y <- matrix(0,60,4);y[c(10,30,50),] <- 20
  clustered <- l0glmFit(X,y,family=poisson('identity'),nonnegative=TRUE,
    precluster=TRUE,control.precluster=list(k.max=2),parallelize=FALSE,
    control.local=local_budget)
  expect_length(calls,2L)
  expect_true(calls[[2]]$clustering$refitted)
  expect_equal(ncol(coef(clustered)),4)
})

test_that('an unbounded initial likelihood cannot be replaced by a finite score', {
  fit <- l0glmFit(diag(3),cbind(c(1,0,0),c(-1,0,0)),lambda=.01,
    grouped=TRUE,relaxed=TRUE,normalize=FALSE,parallelize=FALSE)
  expect_identical(fit$local.search$stopping,'unbounded_criterion')
  expect_equal(nrow(fit$local.search$history),1L)
  expect_identical(fit$local.search$final$value,-Inf)
  expect_false(L0glm_local_better(0,-Inf,1e-8))
  expect_true(L0glm_local_better(0,Inf,1e-8))
})
