test_that("cached Cox products, information and ridge steps equal exact derivatives", {
  skip_if_not_installed("mgcv")
  set.seed(1860)
  for(storage in c("dense","sparse")) for(p in c(6L,60L)) {
    n <- 40L
    X <- matrix(rnorm(n*p),n,p); X[abs(X)<.6] <- 0
    if(storage=="sparse") X <- Matrix::Matrix(X,sparse=TRUE)
    y <- cbind(time=round(rexp(n),1),status=rep(c(1,1,0),length.out=n))
    weights <- sample(0:3,n,TRUE); offset <- runif(n,-.2,.2)
    data <- L0glm_ext_prepare(X,y,mgcv::cox.ph(),weights,offset,FALSE,NULL,FALSE,
      "none",NULL,NULL,list(cox="pcg"))
    expect_identical(inherits(data$Z,"sparseMatrix"),storage=="sparse")
    expect_equal(as.matrix(data$Z),as.matrix(X)[data$cox.cache$order,],ignore_attr=TRUE)
    B <- matrix(rnorm(p,sd=.1),p,1); D <- matrix(runif(p,.5,2),p,1)
    data$prior <- matrix(rnorm(p,sd=.2),p,1)
    s <- L0glm_ext_cox_stats(data,B); A <- data$Z
    a <- rev(cumsum(rev(ifelse(s$S0>0,s$d/s$S0^2,0))))
    W <- diag(s$q)-tcrossprod(s$e)*matrix(a[outer(s$group,s$group,pmax)],n,n)
    H <- as.matrix(Matrix::crossprod(A,W%*%A))
    expect_equal(L0glm_ext_cox_information(A,s),H,tolerance=1e-10)
    v <- rnorm(p)/sqrt(p)
    expect_equal(L0glm_ext_cox_hvp(v,A,s),as.numeric(H%*%v),tolerance=1e-10)
    h <- 1e-5
    gradient <- function(b)as.numeric(Matrix::crossprod(A,L0glm_ext_cox_stats(data,b)$score))
    numerical <- (gradient(B-h*v)-gradient(B+h*v))/(2*h)
    expect_equal(L0glm_ext_cox_hvp(v,A,s),numerical,tolerance=1e-6)
    pcg <- L0glm_ext_cox_pcg_step(data,B,D,s)
    score <- as.numeric(Matrix::crossprod(A,s$score))
    expected <- solve(H+diag(as.numeric(D)),H%*%B+score+as.numeric(D*data$prior))
    expect_equal(as.numeric(pcg),as.numeric(expected),tolerance=1e-8)
    expect_true(attr(pcg,"cox.solve")$converged)
    expect_lte(attr(pcg,"cox.solve")$residual,1e-10)
  }
})

test_that("Cox cache is built once and public observation order is preserved", {
  skip_if_not_installed("mgcv")
  set.seed(1861)
  X <- matrix(rnorm(180*5),180,5); y <- cbind(rexp(180),rep(c(1,0,1),60))
  w <- sample(0:2,180,TRUE); off <- rnorm(180,sd=.1)
  original <- L0glm_ext_cox_cache; builds <- 0L
  fit <- testthat::with_mocked_bindings({
    ans <- L0glm.fit(X,y,family=mgcv::cox.ph(),weights=w,offset=off,lambda=2,
      relaxed=TRUE,maxit=100L,tol=1e-8)
    logLik(ans); L0glm.metrics(ans); L0glm_ext_cox_baseline(ans)
    ans
  },L0glm_ext_cox_cache=function(...) {builds<<-builds+1L;original(...)},.package="L0glm")
  expect_equal(builds,1L)
  expect_equal(fit$X,X)
  expect_equal(as.numeric(fit$weights),as.numeric(w))
  expect_equal(as.numeric(fit$offset),off)
  expect_equal(unname(fit$y),unname(y))
  expect_identical(fit$initialization$backend,"pcg")
  expect_true(fit$initialization$converged)
  expect_lte(fit$initialization$residual,1e-10)
  expect_true(all(vapply(fit$history,function(h)all(h$backend%in%c("llt","pcg","empty")),logical(1))))
  changed <- L0glm.metrics(fit,weights=w+1)
  eta <- as.numeric(predict(fit,type="link"))
  expected <- L0glm_ext_cox_eta(eta,fit$response,w+1)$loglik
  expect_equal(changed$loglik,expected,tolerance=1e-10)
})

test_that("Cox initializer is matrix free and a safeguarded ridge step", {
  skip_if_not_installed("mgcv")
  set.seed(1862)
  X <- Matrix::rsparsematrix(50,250,.15)
  y <- cbind(round(rexp(50),1),rep(c(1,1,0,1,0),10))
  data <- L0glm_ext_prepare(X,y,mgcv::cox.ph(),NULL,NULL,FALSE,NULL,FALSE,
    "rms",NULL,NULL,list(max.working=2L))
  initial <- testthat::with_mocked_bindings(L0glm_ext_cox_initial(data),
    L0glm_ext_cox_direct=function(...)stop("dense solver used"),
    lmridge=function(...)stop("dense ridge initializer used"),
    L0glm_ext_cox_information=function(...)stop("dense information used"),.package="L0glm")
  expect_true(all(is.finite(initial)))
  expect_true(attr(initial,"cox.solve")$converged)
  f <- function(b)-2*L0glm_ext_cox_stats(data,b)$loglik+1e-3*sum(b^2)
  expect_lte(f(initial),f(matrix(0,250,1)))
  expect_true(inherits(data$Z,"sparseMatrix"))
  s <- L0glm_ext_cox_stats(data,matrix(0,250,1)); good <- s$q>0
  ref <- lmridge(data$Z[good,,drop=FALSE],s$score[good]/s$q[good],weights=s$q[good],
    lambdas=rep(1e-3,250),lower=data$lower)
  expect_equal(as.numeric(initial)*2^attr(initial,"cox.solve")$step_halvings,as.numeric(ref),tolerance=1e-7)
})

test_that("Cox direct, PCG and automatic fits match coxph with weights and offsets", {
  skip_if_not_installed("mgcv"); skip_if_not_installed("survival")
  set.seed(1863)
  n <- 350L; p <- 5L
  X <- matrix(rnorm(n*p),n,p); offset <- rnorm(n,sd=.2); w <- sample(0:3,n,TRUE)
  event <- rexp(n,exp(.6*X[,1]-.4*X[,2]+offset)); censor <- rexp(n,.4)
  y <- survival::Surv(round(pmin(event,censor),2),as.numeric(event<=censor))
  ref <- survival::coxph(y~X+offset(offset),weights=replace(w,w==0,1),subset=w>0,
    ties="breslow",control=survival::coxph.control(eps=1e-11))
  for(mode in c("auto","direct","pcg")) {
    fit <- L0glm.fit(X,y,family=mgcv::cox.ph(),weights=w,offset=offset,
      lambda=0,relaxed=FALSE,maxit=100L,tol=1e-9,control.family=list(cox=mode))
    expect_equal(as.numeric(fit$coefficients),unname(coef(ref)),tolerance=1e-7)
    expect_equal(as.numeric(logLik(fit)),as.numeric(logLik(ref)),tolerance=1e-7)
    expect_true(all(fit$convergence$converged))
  }
  shifted <- L0glm.fit(X,y,family=mgcv::cox.ph(),weights=w,offset=offset+1000,
    lambda=0,relaxed=FALSE,maxit=100L,tol=1e-9,control.family=list(cox="pcg"))
  ordinary <- L0glm.fit(X,y,family=mgcv::cox.ph(),weights=w,offset=offset,
    lambda=0,relaxed=FALSE,maxit=100L,tol=1e-9)
  expect_equal(shifted$coefficients,ordinary$coefficients,tolerance=1e-7)
  expect_equal(predict(shifted,type="survival",times=c(0,.1,1,10)),
    predict(ordinary,type="survival",times=c(0,.1,1,10)),tolerance=1e-7)
})

test_that("Cox masks, solver routing, guards and failure diagnostics are explicit", {
  skip_if_not_installed("mgcv")
  set.seed(1864)
  X <- matrix(rnorm(120*8),120,8); y <- cbind(rexp(120),rep(c(1,0,1),40))
  fit <- L0glm.fit(X,y,family=mgcv::cox.ph(),start=c(1,1,rep(0,6)),
    always.include=3L,lambda=2,relaxed=TRUE,maxit=100L,tol=1e-8,
    control.family=list(cox="pcg"))
  expect_null(fit$initialization)
  expect_true(all(fit$coefficients[4:8,]==0))
  expect_true(fit$coefficients[3,]!=0)
  expect_error(L0glm.fit(X,y,family=mgcv::cox.ph(),nonnegative=TRUE,
    control.family=list(cox="pcg")),"unrestricted")
  data <- L0glm_ext_prepare(X,y,mgcv::cox.ph(),NULL,NULL,FALSE,NULL,FALSE,
    "none",1,NULL,list(cox="pcg",cox.cg.maxit=1L))
  B <- data$B; D <- matrix(1,8,1)
  args <- list(solver="eigen",method=1L,control.solver=control.solver.gen("eigen",1L))
  expect_error(L0glm_ext_cox_working(data,B,D,list(),args),"did not converge")
  data$control$cox <- "auto"; data$control$cox.llt.max <- 0L
  expect_warning(other <- L0glm_ext_cox_working(data,B,D,list(),args),"falling back")
  expect_identical(attr(other,"cox.solve")$backend,"direct_fallback")
  data$control$cox <- "direct"; data$control$max.working <- 2L
  expect_error(L0glm_ext_cox_working(data,B,D,list(),args),"max.working")
  data$control$cox <- "pcg"; data$control$cox.cg.maxit <- 1000L
  expect_silent(L0glm_ext_cox_working(data,B,D,list(),args))
  for(control in list(list(cox="bad"),list(cox.cg.tol=0),list(cox.cg.tol=1),
    list(cox.cg.maxit=0),list(cox.llt.max=-1),list(cox.cg.maxit=c(1,2)),list(cox.llt.max=c(1,2))))
    expect_error(L0glm_ext_control(control))
})

test_that("a complete forced-PCG Cox fit never calls a dense or ridge backend", {
  skip_if_not_installed("mgcv")
  set.seed(1865)
  X <- Matrix::rsparsematrix(100,220,.2)
  y <- cbind(rexp(100),rep(c(1,1,0,1),25))
  fit <- testthat::with_mocked_bindings(
    L0glm.fit(X,y,family=mgcv::cox.ph(),lambda=2,relaxed=TRUE,maxit=100L,tol=1e-8,
      control.family=list(cox="pcg",max.working=2L)),
    L0glm_ext_cox_direct=function(...)stop("dense backend called"),
    L0glm_ext_cox_information=function(...)stop("dense Hessian called"),
    L0glm_ext_whiten=function(...)stop("dense whitening called"),
    lmridge=function(...)stop("ridge backend called"),.package="L0glm")
  expect_true(all(is.finite(fit$coefficients)))
  expect_true(inherits(fit$X,"sparseMatrix"))
  expect_equal(fit$initialization$space,"observations")
  expect_lte(fit$initialization$dimension,100L)
})

test_that("Cox automatic routing respects storage and direct constraint solvers", {
  skip_if_not_installed("mgcv"); skip_if_not_installed("osqp")
  set.seed(1866)
  X <- matrix(rnorm(200*4),200,4)
  y <- cbind(round(rexp(200,exp(.7*X[,1]-.6*X[,2])),2),rep(c(1,1,0,1),50))
  dense <- L0glm.fit(X,y,family=mgcv::cox.ph(),lambda=0,relaxed=FALSE,maxit=100L,tol=1e-9)
  sparse <- L0glm.fit(Matrix::Matrix(X,sparse=TRUE),y,family=mgcv::cox.ph(),lambda=0,
    relaxed=FALSE,maxit=100L,tol=1e-9)
  expect_true(all(dense$history$selection$backend=="llt"))
  expect_true(all(sparse$history$selection$backend=="pcg"))
  expect_equal(sparse$coefficients,dense$coefficients,tolerance=1e-7)
  constrained <- L0glm.fit(X,y,family=mgcv::cox.ph(),lambda=0,nonnegative=TRUE,
    solver="osqp",relaxed=FALSE,maxit=100L,tol=1e-8,
    control.solver=control.solver.gen("osqp",abstol=1e-9,reltol=1e-9))
  expect_true(all(constrained$coefficients>=-1e-8))
  expect_true(all(constrained$history$selection$backend=="osqp"))
  s <- L0glm_ext_cox_eta(as.numeric(predict(constrained,type="link")),constrained$response,
    constrained$weights[,1L],constrained$cox.cache)
  score <- as.numeric(crossprod(X[s$order,],s$score)); beta <- as.numeric(constrained$coefficients)
  expect_lt(max(abs(score[beta>1e-6])),1e-4)
  expect_true(all(score[beta<=1e-6]<1e-4))
})

test_that("Cox CV and replay are equivalent across hybrid and matrix-free paths", {
  skip_if_not_installed("mgcv")
  set.seed(1867)
  X <- matrix(rnorm(90*4),90,4)
  y <- cbind(round(rexp(90,exp(.5*X[,1]-.4*X[,2])),2),rep(c(1,1,0),30))
  args <- list(X=X,y=y,family=mgcv::cox.ph(),tune="cv",foldid=rep(1:3,each=30),
    nlambda=3,plot=FALSE,printbest=FALSE,maxit=100L,tol=1e-8)
  a <- do.call(L0glm.tune,c(args,list(control.family=list(cox="auto"))))
  b <- do.call(L0glm.tune,c(args,list(control.family=list(cox="pcg"))))
  expect_equal(a$lambdas,b$lambdas,tolerance=1e-7)
  expect_equal(a$cv$deviance,b$cv$deviance,tolerance=1e-6)
  expect_equal(a$fit_best$coefficients,b$fit_best$coefficients,tolerance=1e-7)
  replay <- do.call(L0glm.fit,c(list(X=X,y=y,lambda=b$lambda_1se),b$fit_args))
  expect_equal(replay$coefficients,b$fit_1se$coefficients,tolerance=1e-10)
})

test_that("Cox CG verifies scaled residuals across large diagonal ranges", {
  diagonal <- 10^seq(-20,20,length.out=20)
  A <- diag(diagonal)+matrix(1e-22,20,20)
  truth <- seq(.1,2,length.out=20)/sqrt(diagonal); rhs <- as.numeric(A%*%truth)
  answer <- L0glm_ext_cox_cg(function(v)as.numeric(A%*%v),rhs,diag(A),1e-10,100L)
  expect_true(answer$converged)
  expect_lte(answer$residual,1e-10)
  expect_equal(answer$value*sqrt(diagonal),truth*sqrt(diagonal),tolerance=1e-8)
  zero <- L0glm_ext_cox_cg(function(v)v,numeric(3),rep(1,3),1e-10,100L)
  expect_equal(zero$iterations,0L)
})
