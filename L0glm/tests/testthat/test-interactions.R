test_that('interaction settings add pairwise then three-way products', {
  set.seed(171)
  X <- matrix(rnorm(300*12),300,12,dimnames=list(NULL,letters[1:12]))
  y <- 1+2*X[,1]-2*X[,2]+2*X[,3]+2*X[,1]*X[,2]+
    2*X[,1]*X[,2]*X[,3]+rnorm(300)
  args <- list(X=X,y=y,lambda='gic',dispersion=1,intercept=TRUE,
    local.search=FALSE,maxit=70L,parallelize=FALSE)
  main <- do.call(l0glmFit,args)
  for(order in 0L) {
    fit <- do.call(l0glmFit,c(args,list(interactions=order)))
    expect_identical(coef(fit),coef(main));expect_identical(GIC(fit),GIC(main))
    expect_null(fit$interactions)
  }
  fit <- do.call(l0glmFit,c(args,list(interactions=2L)))
  expect_equal(fit$interactions$history$order,0:2)
  expect_equal(fit$interactions$selected_order,2L)
  pair <- do.call(l0glmFit,c(args,list(interactions=1L)))
  expect_true(all(lengths(pair$interactions$terms)==2L))
  expect_equal(pair$interactions$history$order,0:1)
  degree <- lengths(fit$interactions$terms)
  expect_true(all(degree%in%2:3));expect_true(any(degree==3L))
  expect_true(all(vapply(fit$interactions$terms,function(j)!anyDuplicated(j),TRUE)))
  expect_lt(mean((y-as.numeric(predict(fit)))^2),1.3)
  expect_equal(predict(fit,X),predict(fit),tolerance=1e-10)
  expect_equal(predict(fit,X[,12:1]),predict(fit),tolerance=1e-10)
  expect_equal(nrow(coef(fit)),13L+length(degree))
  expect_s4_class(coef(fit),'CsparseMatrix')
  expect_s3_class(fit$X,'L0glm_design')
  replay <- do.call(l0glmFit,c(list(X=X,y=y),fit$interactions$specification))
  expect_equal(coef(replay),coef(fit),tolerance=1e-10)
  expect_error(predict(fit,X[,1:3]),'original main-effect')
  for(order in list(-1,NA,Inf,1.5,1e30,c(1,2),TRUE))
    expect_error(l0glmFit(X,y,interactions=order),'nonnegative integer')
  expect_error(do.call(l0glmFit,c(args,list(interactions=1L,
    control.interactions=list(max.terms=1)))),'max.terms')
  expect_error(do.call(l0glmFit,c(args,list(interactions=1L,
    control.interactions=list(max.cells=1)))),'max.cells')
})

test_that('interaction views preserve storage, indexing and original predictor units', {
  set.seed(32);X <- matrix(rnorm(140*6),140,6,dimnames=list(NULL,letters[1:6]))
  terms <- list(c(1L,2L),c(1L,2L,3L))
  explicit <- cbind(X,`a:b`=X[,1]*X[,2],`a:b:c`=X[,1]*X[,2]*X[,3])
  for(sparse in c(FALSE,TRUE))for(intercept in c(FALSE,TRUE)) {
    view <- L0glm_design_view(if(sparse)as(X,'dgCMatrix')else X,terms,intercept)
    expected <- if(intercept)cbind(Intercept=1,explicit)else explicit
    expect_equal(as.matrix(view[,,drop=FALSE]),expected)
    expect_equal(as.matrix(view[c(20,1,3),c(2,ncol(view),1),drop=FALSE]),
      expected[c(20,1,3),c(2,ncol(view),1),drop=FALSE])
    expect_identical(inherits(view[,,drop=FALSE],'sparseMatrix'),sparse)
  }
  y <- 1+2*X[,1]-2*X[,2]+2*X[,1]*X[,2]+rnorm(140,sd=.4)
  args <- list(X=X,y=y,interactions=1L,lambda='gic',intercept=TRUE,
    dispersion=1,local.search=FALSE,maxit=60L,parallelize=FALSE)
  fit <- do.call(l0glmFit,args)
  for(A in list(as(X,'dgCMatrix'),X*rep(10^seq(-2,2,length.out=6),each=140))) {
    alternative <- do.call(l0glmFit,modifyList(args,list(X=A)))
    expect_equal(predict(alternative),predict(fit),tolerance=1e-8)
  }
})

test_that('interaction starts preserve channel masks and compulsory main effects', {
  set.seed(193);X <- matrix(rnorm(160*6),160,6)
  Y <- cbind(1+2*X[,1]-2*X[,2]+2*X[,1]*X[,2],
    2-X[,1]+2*X[,2]-X[,1]*X[,2])+matrix(rnorm(320,sd=.2),160)
  start <- matrix(1,7,2);start[2,1] <- 0;start[4,] <- 0;start[7,] <- 0
  fit <- l0glmFit(X,Y,interactions=1L,intercept=TRUE,always.include=6L,
    start=start,prior.mean=matrix(.1,7,2),lambda='gic',
    dispersion=1,local.search=FALSE,maxit=60L,parallelize=FALSE)
  expect_equal(as.numeric(coef(fit)[2,1]),0)
  expect_true(all(coef(fit)[4,]==0))
  expect_true(7L%in%fit$always.include)
  expect_false(any(vapply(fit$interactions$terms,function(j)any(j%in%c(3L,6L)),TRUE)))
  for(i in seq_along(fit$interactions$terms))if(1L%in%fit$interactions$terms[[i]])
    expect_equal(as.numeric(coef(fit)[7L+i,1]),0)
})

test_that('generic quasi families can select interaction stages', {
  set.seed(92);X <- matrix(rnorm(150*5),150,5)
  y <- 1+2*X[,1]-2*X[,2]+X[,1]*X[,2]+rnorm(150)
  fit <- l0glmFit(X,y,family=quasi('identity','constant'),lambda='gic',
    interactions=1L,intercept=TRUE,local.search=FALSE)
  expect_true(all(is.finite(fit$interactions$history$score)))
  expect_true(all(is.finite(predict(fit,X))))
})

test_that('native Poisson interactions accept compact weights, offsets and preclustering', {
  set.seed(144);n <- 160L;X <- matrix(runif(n*5),n,5)
  mu <- 1+5*X[,1]+5*X[,2]+20*X[,1]*X[,2]
  Y <- matrix(rpois(n*4,rep(mu,4)),n,4)
  w <- runif(n,.7,1.3);offset <- rep(1,n)
  args <- list(X=X,y=Y,interactions=1L,lambda='gic',weights=w,offset=offset,
    local.search=FALSE,maxit=40L,parallelize=FALSE)
  fast <- do.call(l0nnPois,args)
  general <- do.call(l0glmFit,c(args,list(family=poisson(identity),nonnegative=TRUE)))
  expanded <- do.call(l0nnPois,modifyList(args,list(weights=matrix(w,n,4),offset=matrix(offset,n,4))))
  expect_equal(coef(fast),coef(general),tolerance=1e-10)
  expect_equal(coef(fast),coef(expanded),tolerance=1e-10)
  expect_true(all(coef(fast)>=0));expect_true(is.finite(GIC(fast)))
  expect_error(l0nnPois(X,Y,family=gaussian(),interactions=1L),'poisson')
  expect_error(do.call(l0nnPois,c(args,list(precluster=TRUE,cluster.refit=FALSE))),
    'cluster.refit=TRUE')
  clustered <- do.call(l0nnPois,c(args,list(precluster=TRUE,
    control.precluster=list(k.max=2L))))
  expect_equal(ncol(coef(clustered)),4L)
  expect_true(all(is.finite(predict(clustered,X))))
})

test_that('group costs and fold-specific interaction generation remain reproducible', {
  set.seed(321);X <- matrix(rnorm(140*6),140,6,dimnames=list(NULL,letters[1:6]))
  y <- 1+2*X[,1]-2*X[,2]+X[,1]*X[,2]+rnorm(140)
  fold <- rep(1:2,70)
  cv <- l0glmTune(X,y,family=gaussian(),interactions=1L,intercept=TRUE,
    nlambda=3L,tune='cv',nfolds=2L,foldid=fold,plot=FALSE,
    printbest=FALSE,local.search=FALSE,maxit=60L)
  train <- fold!=1L
  manual <- l0glmFit(X[train,],y[train],interactions=1L,intercept=TRUE,
    lambda=cv$cv$lambda_train[1,1],local.search=FALSE,maxit=60L)
  expect_equal(unname(cv$cv$deviance[1,1]),
    mean((y[!train]-as.numeric(predict(manual,X[!train,])))^2),tolerance=1e-10)
  expect_identical(cv$cv$interaction_terms[[1]][[1]],manual$interactions$terms)
  i <- which.min(cv$metrics$cv_deviance)
  expect_equal(cv$coefficients[[i]],coef(cv$fit_best))
  expect_equal(GIC(cv$fit_best),cv$metrics$gic[i])
  expect_true(all(is.finite(predict(cv$fit_1se,X))))
  expect_error(l0glmTune(X,y,family=gaussian(),interactions=1L,start='pilot',
    plot=FALSE,printbest=FALSE),'pilot')
  grouped <- l0glmTune(X,y,family=gaussian(),interactions=1L,group=rep(1:3,each=2),
    intercept=TRUE,tune='gic',nlambda=3L,plot=FALSE,
    printbest=FALSE,local.search='auto',maxit=60L)
  expect_equal(grouped$fit_best$group.cost.type,'dimensions')
  expect_true(isTRUE(attr(grouped$fit_best$group.multiplier,'dimension_based')))
  expect_true(is.finite(GIC(grouped$fit_best)))
  truth <- matrix(c(1,2,-2,0,0,0,0,1),ncol=1,
    dimnames=list(c('Intercept',letters[1:6],'a:b'),NULL))
  path <- l0glmTune(X,y,family=gaussian(),interactions=1L,intercept=TRUE,
    nlambda=3L,tune='gic',plot=FALSE,printbest=FALSE,
    local.search=FALSE,maxit=60L,coefs.true=truth)
  expect_true(all(is.finite(path$metrics$F1_score)))
  expect_equal(path$metrics$tuning_score,path$metrics$gic)
  args <- L0glm_group_tune_cost(X,y,gaussian(),list(group=c(1,1,2,2,0,0),
    group.multiplier=c('1'=2,'2'=2)))
  expect_equal(args$always.include,c(5L,6L))
})
