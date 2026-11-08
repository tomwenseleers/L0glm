test_that('orthogonal group selection preserves opposing signs and matches block L0', {
  set.seed(7);n <- 200L
  X <- qr.Q(qr(matrix(rnorm(n*8L),n,8L)))*sqrt(n)
  truth <- c(.5,-.5,.04,-.04,0,0,.8,-.6)
  fit <- l0glmFit(X,as.numeric(X%*%truth),group=rep(1:4,each=2),
    family=gaussian(),dispersion=1,normalize='none',lambda=4,relaxed=TRUE,
    maxit=100L,local.search=FALSE,parallelize=FALSE)
  truth[3:6] <- 0
  expect_equal(as.numeric(coef(fit)),truth,tolerance=1e-9)
  expect_s4_class(coef(fit),'CsparseMatrix')
  expect_equal(l0glmMetrics(fit)$group_support_size,2L)
  expect_equal(l0glmMetrics(fit)$group_df,4L)
})

test_that('group presets, sparse storage and predictor units give reproducible fits', {
  set.seed(104);X <- matrix(rnorm(180*12),180,12)
  beta <- c(1.5,-1.5,0,0,.8,-.8,rep(0,6))
  y <- .5+as.numeric(X%*%beta)+rnorm(180,sd=.4)
  group <- rep(letters[1:6],each=2)
  args <- list(X=X,y=y,group=group,lambda='gic',intercept=TRUE,dispersion=1,
    relaxed=TRUE,local.search=FALSE,maxit=100L,parallelize=FALSE)
  fit <- do.call(l0glmFit,args)
  sparse <- do.call(l0glmFit,modifyList(args,list(X=as(X,'dgCMatrix'))))
  expect_equal(as.matrix(coef(fit)),as.matrix(coef(sparse)),tolerance=1e-8)
  numeric <- do.call(l0glmFit,modifyList(args,list(lambda=fit$init.lambda,
    group.multiplier=fit$group.multiplier)))
  expect_equal(coef(fit),coef(numeric),tolerance=1e-12)
  scaled <- do.call(l0glmFit,modifyList(args,list(X=sweep(X,2,10^seq(-3,3,length.out=12),'*'))))
  expect_equal(as.numeric(predict(fit)),as.numeric(predict(scaled)),tolerance=1e-8)
  expect_equal(GIC(fit),GIC(numeric),tolerance=1e-12)
  label <- factor(group,levels=rev(unique(group)))
  factored <- do.call(l0glmFit,modifyList(args,list(group=label)))
  expect_equal(coef(fit),coef(factored),tolerance=1e-12)
})

test_that('group zero, compulsory terms and starting-entry masks retain their meanings', {
  set.seed(13);X <- matrix(rnorm(120*8),120,8)
  Y <- cbind(1+X[,1]-X[,2]+rnorm(120),1-X[,1]+X[,2]+rnorm(120))
  group <- c(rep(1:3,each=2),0,0)
  start <- matrix(1,9,2);start[2,1] <- 0;start[5,] <- 0
  for(grouped in c(FALSE,TRUE)) {
    fit <- l0glmFit(X,Y,group=group,intercept=TRUE,start=start,
      lambda='gic',dispersion=1,relaxed=TRUE,grouped=grouped,
      maxit=70L,local.search=FALSE,parallelize=FALSE)
    expect_equal(fit$always.include,c(1L,8L,9L))
    expect_equal(dim(coef(fit)),c(9L,2L))
    expect_equal(as.numeric(coef(fit)[2,1]),0)
    expect_true(all(coef(fit)[5,]==0))
    expect_true(all(fit$group.multiplier==if(grouped)4 else 2))
  }
  fit <- l0glmFit(X,Y[,1],group=group,intercept=TRUE,start=0,
    lambda='gic',dispersion=1,local.search=FALSE)
  expect_true(all(coef(fit)[2:7,]==0))
})

test_that('group tuning and validation retain correct full-data models', {
  set.seed(20);X <- matrix(rnorm(100*6),100,6)
  y <- 1+X[,1]-X[,2]+rnorm(100,sd=.3)
  group <- rep(1:3,each=2)
  cv <- l0glmTune(X,y,family=gaussian(),group=group,intercept=TRUE,
    tune='cv',nfolds=3L,nlambda=3L,dispersion=1,relaxed=TRUE,
    maxit=60L,plot=FALSE,printbest=FALSE,local.search=FALSE)
  fit <- l0glmFit(X,y,family=gaussian(),group=group,intercept=TRUE,
    lambda=cv$lambda_best,group.multiplier=cv$fit_best$group.multiplier,
    dispersion=1,relaxed=TRUE,maxit=60L,local.search=FALSE)
  expect_equal(coef(cv$fit_best),coef(fit),tolerance=1e-9)
  expect_true(all(is.finite(predict(cv$fit_1se,X))))
  validated <- l0glmFit(X,y,group=group,intercept=TRUE,lambda='gic',
    dispersion=1,relaxed=TRUE,local.search=FALSE,
    validation=list(X=X,y=y,starts='default'))
  expect_true(is.finite(validated$validation$deviance))
})

test_that('invalid group specifications and unsupported combinations are explicit', {
  set.seed(9);X <- matrix(rnorm(120),30,4);y <- rnorm(30)
  expect_error(l0glmFit(X,y,group=c(1,NA,2,2)), 'nonmissing')
  expect_error(l0glmFit(X,y,group=1:3), 'one nonmissing label')
  expect_error(l0glmFit(X,y,group=rep(1:2,each=2),group.multiplier=c(1,-1)), 'positive')
  expect_error(l0glmFit(X,y,group=rep(1:2,each=2),local.search='invalid'), 'local.search')
  expect_error(l0nnPois(X,abs(y),family=gaussian(),group=c(1,1,2,2)), 'only defined')
})

test_that('grouped GIC uses the standard n, p and active parameter conventions', {
  set.seed(414);X <- matrix(rnorm(100*8),100,8)
  colnames(X) <- letters[1:8]
  y <- X[,1]-X[,2]+rnorm(100);w <- rep(c(0,1),50)
  fit <- l0glmFit(X,y,group=rep(1:4,each=2),always.include=8L,
    family=gaussian,weights=w,intercept=TRUE,lambda='gic',dispersion=1,
    relaxed=TRUE,local.search=FALSE)
  expect_equal(fit$init.lambda,presetLambda('gic',7,100))
  B <- coef(fit);k <- sum(B!=0)+sum(B[fit$always.include,,drop=FALSE]==0)
  scores <- l0glmMetrics(fit,mu_floor=0)
  expect_equal(scores$df,k)
  expect_equal(scores$gic,-2*scores$loglik+fit$init.lambda*k)
})

test_that('local refinement moves complete predictor groups and respects exclusions', {
  set.seed(171);X <- matrix(rnorm(180*12),180,12)
  group <- rep(letters[1:6],each=2);start <- rep(1,12);start[11] <- 0
  y <- 1+X[,1]-X[,2]+.5*X[,3]-.5*X[,4]+rnorm(180)
  for(joint in c(FALSE,TRUE)) {
    Y <- if(joint)cbind(y,-y)else y
    fit <- l0glmFit(X,Y,group=group,lambda='gic',start=c(0,start),dispersion=1,
      grouped=joint,intercept=TRUE,relaxed=TRUE,local.search=TRUE)
    expect_true(fit$local.search$applied)
    expect_true(all(diff(fit$local.search$history$value)<0))
    expect_lte(fit$local.search$final$value,fit$local.search$initial$value)
    expect_equal(as.numeric(coef(fit)[12,]),rep(0,NCOL(Y)))
    expect_true(is.finite(GIC(fit)))
    expect_equal(l0glmMetrics(fit)$df,sum(coef(fit)!=0)+
      sum(coef(fit)[fit$always.include,,drop=FALSE]==0))
    expect_equal(fit$local.search$final$objective,GIC(fit),tolerance=1e-7)
  }
})

test_that('nonnegative Poisson and scalar noise families support predictor groups', {
  set.seed(81);X <- matrix(runif(140*6),140,6)
  mu <- as.numeric(2+X%*%c(1,1,0,0,1,1))
  group <- rep(1:3,each=2)
  Y <- cbind(rpois(140,mu),rpois(140,mu))
  fit <- l0nnPois(X,Y,group=group,offset=2,lambda='gic',relaxed=TRUE,
    local.search=FALSE,maxit=60L,parallelize=FALSE)
  expect_true(all(coef(fit)>=0))
  expect_true(is.finite(GIC(fit)))
  expect_equal(dim(coef(fit)),c(6L,2L))
  noise <- list(nb1=nb1(1,'identity'),scale=gammaScale(1,'identity'),
    censored=gammaCensored(.3,'identity',censor=.5))
  response <- list(nb1=rnbinom(140,mu=mu,size=mu),
    scale=rgamma(140,shape=mu,scale=1),censored=rgamma(140,shape=1/.3,scale=mu*.3))
  response$censored[response$censored<.5] <- 0
  for(name in names(noise)) {
    fit <- l0glmFit(X,response[[name]],group=group,offset=2,lambda='gic',
      family=noise[[name]],nonnegative=TRUE,relaxed=TRUE,local.search=FALSE,maxit=70L)
    expect_true(all(coef(fit)>=0))
    expect_true(is.finite(GIC(fit)))
    expect_true(all(is.finite(predict(fit))))
  }
})
