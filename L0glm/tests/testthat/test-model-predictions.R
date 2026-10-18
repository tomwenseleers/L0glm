test_that("ordinary predictions match glm on training and original-scale new data", {
  set.seed(2601)
  n <- 120L
  X <- sweep(matrix(rnorm(n*3),n,3),2,c(.3,2,5),"*")
  X <- sweep(X,2,c(2,-1,3),"+")
  colnames(X) <- c("a","b","c"); rownames(X) <- paste0("row",seq_len(n))
  off <- runif(n,-.2,.2); w <- runif(n,.5,2)
  eta <- as.numeric(.2 + X %*% c(.4,-.2,.1) + off)
  newX <- X[1:5,,drop=FALSE] + .1; newoff <- off[1:5] + .05
  for (family in list(gaussian(),poisson(),binomial(),binomial("probit"),quasibinomial())) {
    mu <- family$linkinv(eta)
    y <- switch(family$family,gaussian=rnorm(n,mu,sd=.5),
      poisson=rpois(n,mu),binomial=rbinom(n,1,mu),quasibinomial=rbinom(n,1,mu))
    weights <- if (family$family %in% c("binomial","quasibinomial")) rep(1:3,length.out=n) else w
    frame <- data.frame(y=y,X,off=off,w=weights)
    ref <- glm(y ~ a+b+c+offset(off),data=frame,weights=w,family=family,
      control=glm.control(epsilon=1e-12,maxit=100))
    for (normalization in c("none","rms")) {
      fit <- L0glm.fit(X,y,family=family,weights=weights,offset=off,intercept=TRUE,
        normalize=normalization,lambda=0,thresh=0,relaxed=FALSE,
        parallelize=FALSE,maxit=100,tol=1e-14)
      expect_equal(as.numeric(fitted(fit)),as.numeric(fitted(ref)),tolerance=1e-6)
      expect_equal(fitted.values(fit),fitted(fit))
      expect_equal(predict(fit),fitted(fit))
      for (type in c("response","link")) {
        expected <- as.numeric(predict(ref,newdata=data.frame(newX,off=newoff),type=type))
        for (design in list(newX,newX[,3:1],as.data.frame(newX[,3:1]),
          Matrix::Matrix(newX,sparse=TRUE),Matrix::Matrix(newX,sparse=FALSE))) {
          value <- predict(fit,design,type=type,offset=newoff)
          expect_equal(dim(value),c(5L,1L))
          expect_equal(rownames(value),rownames(newX))
          expect_equal(as.numeric(value),expected,tolerance=1e-6)
        }
      }
      expect_equal(as.numeric(predict(fit,newX,type="link")),
        as.numeric(cbind(1,newX) %*% coef(fit)))
      expect_equal(predict(fit,newX[1,,drop=FALSE]),predict(fit,newX[1,]),ignore_attr=TRUE)
      expect_equal(dim(predict(fit,newX[FALSE,,drop=FALSE])),c(0L,1L))
      if (family$family %in% c("binomial","quasibinomial")) expect_equal(predict(fit,type="probabilities"),fitted(fit))
      else expect_error(predict(fit,type="probabilities"),"Prediction type")
    }
  }
})

test_that("native multichannel predictions retain offsets, intercepts and response dimensions", {
  set.seed(2602)
  X <- matrix(runif(180),60,3,dimnames=list(NULL,c("a","b","c")))
  B <- matrix(c(2,1,0,0,2,1),3,2)
  off <- cbind(runif(60,.5,1),runif(60,1,2))
  Y <- matrix(rpois(120,as.numeric(X %*% B + off)),60,2,
    dimnames=list(NULL,c("mass1","mass2")))
  for (sparse in c(FALSE,TRUE)) {
    fit <- L0nnpois(Matrix::Matrix(X,sparse=sparse),Y,offset=off,
      intercept=TRUE,parallelize=FALSE,lambda=2)
    expected <- as.matrix(cbind(1,X) %*% coef(fit)) + off
    expect_equal(fitted(fit),expected,ignore_attr=TRUE)
    expect_equal(colnames(fitted(fit)),colnames(Y))
    expect_equal(predict(fit,X[1:4,],offset=off[1:4,]),expected[1:4,,drop=FALSE],ignore_attr=TRUE)
    expect_equal(predict(fit,X[1:4,]),expected[1:4,,drop=FALSE]-off[1:4,],ignore_attr=TRUE)
    expect_equal(fitted(fit,offset=0),expected-off,ignore_attr=TRUE)
    zero <- fit; zero$coefficients[,] <- 0
    expect_equal(as.numeric(fitted(zero,offset=0)),numeric(length(Y)))
  }
  # Manual intercept columns are retained when no intercept was added by the API.
  manual <- L0glm.fit(cbind(Intercept=1,X),Y,offset=off,family=poisson(identity),
    nonnegative=TRUE,intercept=FALSE,lambda=2,parallelize=FALSE)
  expect_equal(predict(manual,cbind(Intercept=1,X)[1:4,],offset=off[1:4,]),
    fitted(manual)[1:4,,drop=FALSE])
})

test_that("prediction validates dimensions, names, finite values and zero-start columns", {
  set.seed(2603)
  X <- matrix(rnorm(240),80,3,dimnames=list(NULL,c("a","b","c")))
  y <- X[,1]+rnorm(80)
  fit <- L0glm.fit(X,y,intercept=TRUE,start=c(0,1,0,0),lambda=2,parallelize=FALSE)
  expect_equal(as.numeric(coef(fit)[3:4,]),c(0,0))
  expect_equal(as.numeric(predict(fit,X[1:4,])),as.numeric(cbind(1,X[1:4,]) %*% coef(fit)))
  expect_error(predict(fit,X[,1:2]),"prediction columns")
  wrong <- X; colnames(wrong) <- c("a","b","unknown")
  expect_error(predict(fit,wrong),"column names")
  wrong <- X; wrong[1,1] <- NA_real_
  expect_error(predict(fit,wrong),"finite numeric")
  expect_error(predict(fit,data.frame(a=letters[1:3],b=1:3,c=1:3)),"numeric predictors")
  expect_error(predict(fit,X[1:3,],offset=matrix(0,1,3)),"offset")
  expect_error(predict(fit,X[1:3,],offset=Inf),"finite")
  expect_error(predict(fit,X[1:3,],offset="bad"),"numeric")
  expect_error(predict(fit,type="unknown"),"Prediction type")
  expect_error(fitted(fit,newdata=X),"Use predict")
  single <- L0glm.fit(X[,1,drop=FALSE],y,lambda=2,parallelize=FALSE)
  expect_equal(predict(single,X[1:4,1]),predict(single,X[1:4,1,drop=FALSE]))
  unnamed <- L0glm.fit(unname(X),y,lambda=2,intercept=TRUE,parallelize=FALSE)
  expect_equal(predict(unnamed,as.data.frame(X[1:4,])),predict(unnamed,unname(X[1:4,])))
})

test_that("multinomial predictions are stable probabilities and deterministic classes", {
  skip_if_not_installed("mgcv"); skip_if_not_installed("nnet")
  set.seed(2604)
  X <- matrix(rnorm(360),120,3,dimnames=list(NULL,c("a","b","c")))
  Y <- factor(rep(c("A","B","C"),40))
  off <- matrix(rnorm(360,sd=.1),120,3)
  fit <- L0glm.fit(X,Y,family=mgcv::multinom(K=2),offset=off,
    intercept=TRUE,lambda=0,relaxed=FALSE,maxit=100,tol=1e-9)
  expected_eta <- as.matrix(cbind(1,X) %*% coef(fit)) + off
  expected_eta <- expected_eta-rowMeans(expected_eta)
  expect_equal(predict(fit,type="link"),expected_eta,ignore_attr=TRUE)
  probability <- fitted(fit)
  expect_equal(rowSums(probability),rep(1,120),tolerance=1e-13)
  expect_equal(colnames(probability),levels(Y))
  expect_equal(probability,predict(fit,type="probabilities"))
  expect_equal(probability,predict(fit,offset=off+1000),tolerance=1e-12)
  expect_equal(predict(fit,X[1:5,],offset=off[1:5,]),probability[1:5,,drop=FALSE])
  # Compare against an independently fitted multinomial likelihood without offsets.
  plain <- L0glm.fit(X,Y,family=mgcv::multinom(K=2),intercept=TRUE,
    lambda=0,relaxed=FALSE,maxit=100,tol=1e-9)
  ref <- nnet::multinom(Y ~ a+b+c,data=data.frame(Y=Y,X),trace=FALSE,
    maxit=300,abstol=1e-12,reltol=1e-12)
  expect_equal(predict(plain,X[1:5,]),predict(ref,newdata=as.data.frame(X[1:5,]),type="probs"),ignore_attr=TRUE,tolerance=1e-5)
  fit$coefficients[,] <- 0
  set.seed(9); before <- .Random.seed
  expect_equal(as.character(predict(fit,X[1:4,],type="class")),rep("A",4))
  expect_identical(.Random.seed,before)
  expect_equal(dim(predict(fit,X[FALSE,,drop=FALSE])),c(0L,3L))
  expect_equal(length(predict(fit,X[FALSE,,drop=FALSE],type="class")),0)
})

test_that("NB means and Cox risks use the original linear predictor", {
  skip_if_not_installed("mgcv"); skip_if_not_installed("survival")
  set.seed(2605)
  X <- matrix(rnorm(300),100,3,dimnames=list(NULL,c("a","b","c")))
  off <- runif(100,-.3,.3)
  nb <- L0glm.fit(X,rnbinom(100,mu=exp(.2+.4*X[,1]+off),size=2),
    family=mgcv::nb(theta=2),offset=off,lambda=2,intercept=TRUE,relaxed=FALSE)
  expect_equal(fitted(nb),exp(as.matrix(cbind(1,X) %*% coef(nb)) + off),ignore_attr=TRUE)
  expect_equal(predict(nb,X[1:4,],offset=off[1:4]),fitted(nb)[1:4,,drop=FALSE])
  event <- rexp(100,exp(.5*X[,1]+off)); censor <- rexp(100,.3)
  Y <- survival::Surv(pmin(event,censor),event <= censor)
  fit <- L0glm.fit(X,Y,family=mgcv::cox.ph(),offset=off,lambda=0,
    relaxed=FALSE,maxit=100,tol=1e-9)
  ref <- survival::coxph(Y ~ X+offset(off),ties="breslow")
  expected <- exp(as.numeric(X %*% coef(fit)) + off)
  expect_equal(fitted(fit),expected)
  expect_equal(fitted(fit,type="risk"),expected)
  expect_equal(predict(fit,X[1:5,],offset=off[1:5],type="risk"),expected[1:5])
  expect_equal(fitted(fit,type="link"),predict(fit,type="link"))
  expect_equal(as.numeric(coef(fit)),unname(coef(ref)),tolerance=1e-6)
  # predict.coxph centers offsets even with reference="zero"; remove that shift.
  expect_equal(fitted(fit),as.numeric(predict(ref,type="risk",reference="zero"))*exp(mean(off)),tolerance=1e-6)
  expect_error(predict(fit,type="survival"),"Supply times")
  expect_error(predict(fit,type="survival",times=-1),"nonnegative")
  expect_equal(dim(predict(fit,X[FALSE,,drop=FALSE],type="risk")),NULL)
  survival <- fitted(fit,type="survival",times=c(0,.1,1,10))
  expect_true(all(survival >= 0 & survival <= 1))
  expect_true(all(apply(survival,1,function(x) all(diff(x)<=0))))
  fresh <- fit; fresh$baseline <- NULL
  expect_equal(predict(fresh,type="survival",times=c(0,.1,1,10)),survival)
})

test_that("collection predictions retain method order and unavailable entries without refitting", {
  set.seed(2606)
  X <- matrix(rnorm(240),80,3)
  y <- X[,1]+rnorm(80)
  expect_warning(result <- sparseglm(X,y,offset=1,method=c("L0glm","lasso"),nlambda=3,
    options=list(L0glm=list(solver="invalid"))),"failed")
  predictions <- predict(result,X[1:5,],offset=2)
  expect_identical(names(predictions),names(result$fits))
  expect_null(predictions$L0glm)
  expect_equal(predictions$lasso,predict(result$fits$lasso,X[1:5,],offset=2))
  expect_equal(fitted(result),predict(result))
  expect_equal(fitted(result)$lasso,fitted(result$fits$lasso))
  expect_error(fitted(result,newdata=X),"Use predict")
  local_mocked_bindings(L0glm.fit=function(...) stop("Must not refit"),
    L0nnpois=function(...) stop("Must not refit"),L0glm.tune=function(...) stop("Must not tune"))
  expect_equal(predict(result,X[1:5,],offset=2),predictions)
})

test_that("special-family comparisons provide probabilities, risks and lazy survival predictions", {
  skip_if_not_installed("mgcv"); skip_if_not_installed("survival")
  set.seed(2607)
  X <- matrix(rnorm(240),80,3)
  Y <- factor(rep(1:3,length.out=80))
  multi <- sparseglm(X,Y,family="multinomial",method=c("abess","lasso"),nlambda=3,support.size=0:2)
  for (value in predict(multi,X[1:4,],type="probabilities")) {
    expect_equal(dim(value),c(4L,3L))
    expect_equal(rowSums(value),rep(1,4),tolerance=1e-13)
  }
  event <- rexp(80,exp(.5*X[,1])); censor <- rexp(80,.3)
  Y <- survival::Surv(pmin(event,censor),event<=censor)
  cox <- sparseglm(X,Y,family="cox",method=c("abess","lasso"),nlambda=3,support.size=0:2)
  risks <- predict(cox,X[1:4,],type="risk")
  expect_true(all(vapply(risks,function(x) length(x)==4 && all(x>0),logical(1))))
  for (method in names(cox$fits)) {
    expect_null(cox$fits[[method]]$baseline)
    result <- predict(cox$fits[[method]],X[1:4,],type="survival",times=c(0,1,10))
    expect_equal(dim(result),c(4L,3L))
    expect_true(all(result>=0 & result<=1))
  }
})
