test_that("one-round ISIS matches the direct candidate fit and original-p GIC", {
  set.seed(2511);n <- 180;X <- matrix(rnorm(n*80),n,80);colnames(X) <- paste0("x",1:80)
  off <- runif(n,-.3,.3);w <- sample(1:3,n,TRUE)
  y <- 2+X[,1]*1.4-X[,2]*1.3+off+rnorm(n,sd=.5)
  result <- sis_glm(X,y,weights=w,offset=off,nsis=25,iterative=TRUE,niter=1,
    intercept=TRUE,always.include=3,return_stats=TRUE)
  selected <- result$candidates
  direct <- L0glm.fit(X[,selected],y,weights=w,offset=off,lambda=result$lambda,
    intercept=TRUE,always.include=match(3,selected),relaxed=TRUE,parallelize=FALSE)
  rows <- c(1,selected+1)
  expect_equal(as.matrix(coef(result$fit)[rows,,drop=FALSE]),as.matrix(coef(direct)),tolerance=1e-8)
  expect_true(all(coef(result$fit)[setdiff(1:81,rows),]==0))
  expect_s4_class(coef(result$fit),"sparseMatrix")
  expect_identical(dim(coef(result$fit)),c(81L,1L))
  expect_equal(result$lambda,preset.lambda("gic",79,n))
  expect_equal(result$fit$p,79)
  expect_equal(GIC(result$fit),result$history$gic[1],tolerance=1e-8)
  expect_equal(result$history$gic[1],-2*as.numeric(logLik(direct))+log(79)*log(log(n))*attr(logLik(direct),"df"))
  expect_equal(predict(result$fit,X,offset=off),predict(direct,X[,selected],offset=off),tolerance=1e-8)
  expect_equal(fitted(result$fit),predict(result$fit,X,offset=off))
  expect_equal(predict(result$fit,X[,80:1],offset=off),fitted(result$fit))
  expect_equal(summary(result$fit)$metrics$gic,result$history$gic[1],tolerance=1e-8)
  expect_true(3L%in%result$indices)
  expect_identical(result$stats$keep,result$stats$feature%in%result$indices)
  expect_length(result$ranking,80)
  expect_length(result$fit$preprocessing$scale,81)
  expect_length(result$fit$lambdas,81)
})

test_that("ISIS keeps the best round, starts afresh and stops on a repeated candidate set", {
  set.seed(2512);n <- 180;X <- matrix(rnorm(n*300),n,300)
  y <- as.numeric(X[,1:8]%*%rep(c(1,-1),4))+rnorm(n)
  result <- sis_glm(X,y,nsis=35,iterative=TRUE,niter=10,return_stats=FALSE)
  expect_equal(result$fit$metrics$gic,min(result$history$gic),tolerance=1e-8)
  expect_equal(result$best_iteration,which.min(result$history$gic))
  expect_true(result$converged);expect_lt(result$iterations,10)
  expect_null(result$stats);expect_null(result$ranking)
  expect_equal(result$nsis,length(result$indices));expect_equal(result$screen_size,35)
  direct <- L0glm.fit(X[,result$candidates],y,lambda=result$lambda,relaxed=TRUE,parallelize=FALSE)
  expect_equal(as.numeric(coef(result$fit)[result$candidates,]),as.numeric(coef(direct)),tolerance=1e-7)
  expect_equal(as.numeric(coef(result$fit)[-result$candidates,]),rep(0,300-length(result$candidates)))
  again <- sis_glm(X,y,nsis=35,iterative=TRUE,niter=10,return_stats=FALSE)
  expect_equal(coef(result$fit),coef(again$fit));expect_equal(result$history,again$history)
})

test_that("ISIS preserves scalar starts, initial zero masks, prior centers and mandatory rows", {
  set.seed(2513);n <- 140;X <- matrix(rnorm(n*30),n,30)
  y <- 1+X[,1]*1.5+rnorm(n)
  for(norm in c("none","rms","l2")) {
    result <- sis_glm(X,y,nsis=10,iterative=TRUE,niter=1,start=1,normalize=norm,intercept=TRUE)
    direct <- L0glm.fit(X[,result$candidates],y,start=1,normalize=norm,intercept=TRUE,
      lambda=result$lambda,relaxed=TRUE,parallelize=FALSE)
    expect_equal(as.numeric(coef(result$fit)[c(1,result$candidates+1),]),as.numeric(coef(direct)),tolerance=1e-8)
  }
  Y <- cbind(y,-y);start <- matrix(1,31,2)
  start[5,] <- 0;start[6,1] <- 0;start[8,] <- 0
  prior <- matrix(0,31,2);prior[5,] <- 100
  result <- sis_glm(X,Y,nsis=12,iterative=TRUE,niter=3,intercept=TRUE,
    start=start,prior.mean=prior,always.include=7,grouped=FALSE)
  expect_true(all(coef(result$fit)[5,]==0));expect_equal(coef(result$fit)[6,1],0)
  expect_false(4%in%result$candidates);expect_true(7%in%result$candidates)
  expect_true(7%in%result$indices);expect_identical(dim(coef(result$fit)),c(31L,2L))
  zero <- sis_glm(X,y,iterative=TRUE,start=0,return_stats=FALSE)
  expect_true(all(coef(zero$fit)==0));expect_length(zero$indices,0)
  expect_equal(zero$iterations,1);expect_true(zero$converged)
})

test_that("ISIS fits ordinary families and noncanonical links, with quasi fallback", {
  set.seed(2514);n <- 160;X <- matrix(rnorm(n*16),n,16);off <- runif(n,-.02,.02)
  cases <- list(gaussian(),gaussian("log"),binomial(),binomial("probit"),poisson(),poisson("identity"),
    Gamma("log"),Gamma("inverse"),inverse.gaussian("log"),inverse.gaussian("inverse"),quasibinomial(),quasipoisson())
  for(f in cases) {
    eta <- if(f$link%in%c("identity","inverse")) 2+off else .5+off
    mu <- f$linkinv(eta)
    y <- if(f$family%in%c("binomial","quasibinomial")) rbinom(n,1,mu) else
      if(f$family%in%c("poisson","quasipoisson")) rpois(n,mu) else
        if(f$family=="gaussian") mu+rnorm(n,0,.15) else rgamma(n,shape=30,scale=mu/30)
    result <- sis_glm(X,y,f,offset=off,nsis=6,iterative=TRUE,niter=2,intercept=TRUE,
      nonnegative=f$family=="poisson" && f$link=="identity",return_stats=FALSE)
    expect_s3_class(result$fit,"L0glm_fit");expect_true(all(is.finite(coef(result$fit))))
    expect_identical(dim(coef(result$fit)),c(17L,1L))
    expect_true(all(is.finite(predict(result$fit,X,offset=off))))
    if(f$family=="quasibinomial") {
      expect_true(all(is.na(result$history$gic)))
      expect_equal(result$best_iteration,result$iterations)
    }
    if(f$family=="quasipoisson") {
      expect_true(all(is.finite(result$history$gic)))
      expect_equal(result$best_iteration, which.min(result$history$gic))
    }
  }
})

test_that("ISIS supports grouped Poisson and special-family models with usable full-size fits", {
  skip_if_not_installed("mgcv");skip_if_not_installed("survival")
  set.seed(2515);n <- 180;X <- Matrix::rsparsematrix(n,20,.2)
  P <- abs(X);Y <- matrix(rpois(n*3,2+as.numeric(P[,1])*2),n,3)
  pois <- sis_glm(P,Y,poisson(identity),nonnegative=TRUE,nsis=8,iterative=TRUE,niter=2)
  expect_s4_class(pois$fit$X,"sparseMatrix")
  expect_true(all(coef(pois$fit)>=0));expect_identical(dim(coef(pois$fit)),c(20L,3L))
  expect_true(all(is.finite(L0glm.metrics(pois$fit)$gic)))
  off <- matrix(runif(n*2,-.1,.1),n,2);w <- sample(1:3,n,TRUE)
  counts <- matrix(rnbinom(n*2,mu=exp(.5+off),size=2),n,2)
  nb <- sis_glm(X,counts,mgcv::nb(),weights=w,offset=off,nsis=8,iterative=TRUE,niter=2,intercept=TRUE)
  expect_identical(dim(coef(nb$fit)),c(21L,2L));expect_length(nb$theta,2)
  expect_true(all(is.finite(nb$theta)&nb$theta>0))
  expect_equal(nb$fit$metrics$gic,min(nb$history$gic),tolerance=1e-7)
  labels <- factor(sample(letters[1:3],n,TRUE))
  start <- matrix(1,21,3);start[5,] <- 0
  mn <- sis_glm(X,labels,mgcv::multinom(2),nsis=8,iterative=TRUE,niter=2,start=start)
  expect_identical(dim(coef(mn$fit)),c(21L,3L));expect_true(all(coef(mn$fit)[5,]==0))
  expect_equal(unname(Matrix::rowSums(coef(mn$fit))),rep(0,21),tolerance=1e-10)
  prob <- predict(mn$fit,X,type="probabilities")
  expect_equal(rowSums(prob),rep(1,n),tolerance=1e-10)
  event <- rexp(n,exp(as.numeric(X[,1])*.5));outcome <- survival::Surv(pmin(event,2),event<=2)
  cox <- sis_glm(X,outcome,mgcv::cox.ph(),weights=w,nsis=8,iterative=TRUE,niter=2)
  expect_identical(dim(coef(cox$fit)),c(20L,1L))
  expect_true(all(is.finite(predict(cox$fit,X,type="risk"))))
  expect_true(all(is.finite(predict(cox$fit,X[1:3,],type="survival",times=c(0,1)))))
  expect_equal(cox$lambda,preset.lambda("gic",20,sum(w*outcome[,2])))
})

test_that("ISIS selects joint multinomial effects with grouped and CLR-entry penalties", {
  skip_if_not_installed("mgcv")
  set.seed(2517);n <- 220;X <- matrix(rnorm(n*15),n,15)
  logits <- cbind(2*X[,1],-2*X[,1],0)
  prob <- exp(logits-apply(logits,1,max));prob <- prob/rowSums(prob)
  y <- factor(vapply(1:n,function(i) sample(letters[1:3],1,prob=prob[i,]),character(1)),levels=letters[1:3])
  for(grouped in c(TRUE,FALSE)) {
    fit <- sis_glm(X,y,mgcv::multinom(2),nsis=5,iterative=TRUE,niter=2,grouped=grouped,
      return_stats=FALSE)
    expect_true(1L%in%fit$indices)
    expect_equal(unname(Matrix::rowSums(coef(fit$fit))),rep(0,16),tolerance=1e-10)
    expect_equal(rowSums(predict(fit$fit,X)),rep(1,n),tolerance=1e-10)
    expect_equal(fit$fit$metrics$gic,min(fit$history$gic),tolerance=1e-7)
  }
})

test_that("ISIS validates outer controls and forwards model-specific controls", {
  set.seed(2516);X <- matrix(rnorm(300),60,5);y <- rnorm(60)
  for(k in list(0,-1,1.2,NA,Inf,c(1,2),"2")) expect_error(sis_glm(X,y,iterative=TRUE,niter=k),"niter")
  expect_error(sis_glm(X,y,iterative=NA),"TRUE or FALSE")
  expect_error(sis_glm(X,y,iterative=TRUE,nuisance=matrix(1,60,1)),"always.include")
  expect_error(sis_glm(X,y,lambda="gic"),"iterative=TRUE")
  expect_error(sis_glm(X,y,iterative=TRUE,invalid=3),"named L0glm.fit")
  expect_error(sis_glm(X,y,iterative=TRUE,start=rep(1,4)),"dimensions")
  fit <- sis_glm(X,y,iterative=TRUE,niter=1,lambda=2,maxit=3,tol=1e-6,normalize="none")
  expect_equal(fit$lambda,2);expect_equal(fit$fit$preprocessing$mode,"none")
})

test_that("large identity-Poisson offsets use a feasible null-fit retry", {
  set.seed(2519);n <- 80;X <- matrix(rnorm(n*6),n,6)
  off <- rep(c(0,100),each=n/2);y <- rep(1,n)
  expect_no_warning(result <- sis_glm(X,y,poisson(identity),offset=off,nsis=3))
  b <- uniroot(function(b) sum(y/(off+b)-1),c(1e-8,2),tol=1e-12)$root
  context <- L0glm_sis_context(y,poisson(identity),n,NULL,off,NULL)
  expect_equal(context$score[,1],(y-(off+b))/(off+b),tolerance=1e-7)
  expect_true(all(is.finite(result$stats$stat)))
  expect_error(sis_glm(X,-y,poisson(identity),offset=off),"negative")
})
