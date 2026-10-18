test_that("multinomial adapters reconstruct reference coefficients and score counts in CLR", {
  skip_if_not_installed("mgcv"); skip_if_not_installed("nnet")
  set.seed(2410)
  X <- matrix(rnorm(300),100,3)
  B <- rbind(c(.5,-.5,0),c(0,.3,-.3),c(0,0,0))
  P <- L0glm_ext_softmax(X %*% B)$p
  Y <- t(vapply(seq_len(100),function(i) rmultinom(1,3,P[i,]),numeric(3)))
  colnames(Y) <- c("A","B","C")
  result <- sparseglm(X,Y,family="multinomial",method=c("abess","lasso"),
    support.size=3,lambda=list(lasso=0),nlambda=3,beta_true=B,
    options=list(abess=list(max.newton.iter=100,newton.thresh=1e-10),lasso=list(thresh=1e-12)))
  expect_true(all(result$timings$status == "ok"))
  ref <- nnet::multinom(Y ~ X,trace=FALSE,maxit=300,abstol=1e-12,reltol=1e-12)
  constants <- sum(lgamma(rowSums(Y)+1)-rowSums(lgamma(Y+1)))
  for (fit in result$fits) {
    expect_s4_class(coef(fit),"sparseMatrix")
    expect_equal(dim(coef(fit)),c(4L,3L))
    expect_equal(as.numeric(Matrix::rowSums(coef(fit))),rep(0,4),tolerance=1e-13)
    expect_equal(as.numeric(logLik(fit)),as.numeric(logLik(ref))+constants,tolerance=1e-6)
    expect_equal(predict(fit),fitted(ref),ignore_attr=TRUE,tolerance=1e-5)
  }
  expect_true(all(c("MCC","within_row_MCC") %in% names(L0glm.metrics(result))))
  label <- factor(apply(Y,1,which.max),levels=1:3)
  fit <- sparseglm(X,label,family=mgcv::multinom(K=2),method="adaptive_lasso",
    tune="cv",nlambda=4,nfolds=2,grouped=FALSE)
  expect_equal(fit$timings$status,"ok")
  expect_equal(fit$timings$scope,"ungrouped")
  expect_equal(length(fit$tuning$adaptive_lasso$adaptive$penalty.factor),3)
  for (k in 1:2) {
    train <- fit$foldid != k
    pilot <- L0glm.fit(X[train,],label[train],family=mgcv::multinom(K=2),
      intercept=TRUE,lambda="gic",start=NULL,maxit=1,grouped=FALSE,relaxed=FALSE,
      thresh=0,parallelize=FALSE)
    expect_equal(fit$tuning$adaptive_lasso$cv$adaptive[[k]]$coefficients,coef(pilot))
  }
  expect_warning(unsupported <- sparseglm(X,label,family="multinomial",
    method=c("MCP","SCAD","L0Learn"),nlambda=3),"unavailable")
  expect_true(all(L0glm.metrics(unsupported)$status == "unsupported"))
})

test_that("Cox adapters use common Breslow likelihood and full-risk-set CV subtraction", {
  skip_if_not_installed("mgcv"); skip_if_not_installed("survival")
  set.seed(2411)
  X <- matrix(rnorm(300),100,3)
  event <- rexp(100,exp(.6*X[,1])); censor <- rexp(100,.3)
  Y <- survival::Surv(pmin(event,censor),event <= censor)
  result <- sparseglm(X,Y,family="cox",method=c("abess","lasso"),
    support.size=3,lambda=list(lasso=0),nlambda=3,
    options=list(abess=list(max.newton.iter=100,newton.thresh=1e-10),lasso=list(thresh=1e-12)))
  ref <- survival::coxph(Y ~ X,ties="breslow",control=survival::coxph.control(eps=1e-11))
  expect_true(all(result$timings$status == "ok"))
  for (fit in result$fits) {
    expect_equal(dim(coef(fit)),c(3L,1L))
    expect_equal(as.numeric(logLik(fit)),as.numeric(logLik(ref)),tolerance=1e-7)
    expect_equal(as.numeric(coef(fit)),unname(coef(ref)),tolerance=2e-5)
    expect_equal(nobs(fit),sum(Y[,2]))
  }
  weights <- rep(c(.5,1,2,3),25)
  expect_warning(weighted <- sparseglm(X,Y,family="cox",method=c("abess","lasso"),weights=weights,
    support.size=3,nlambda=3,lambda=list(lasso=0),
    options=list(lasso=list(thresh=1e-12))),"unavailable")
  expect_equal(weighted$timings$status,c("unsupported","ok"))
  ref_weighted <- survival::coxph(Y ~ X,weights=weights,ties="breslow",
    control=survival::coxph.control(eps=1e-11))
  expect_equal(as.numeric(coef(weighted$fits$lasso)),unname(coef(ref_weighted)),tolerance=2e-5)
  expect_equal(as.numeric(logLik(weighted$fits$lasso)),as.numeric(logLik(ref_weighted)),tolerance=1e-7)
  tied <- Y; tied[,1] <- round(tied[,1],1)
  expect_warning(no <- sparseglm(X,tied,family="cox",method="abess",nlambda=3),"unavailable")
  expect_equal(no$timings$status,"unsupported")
  with_zero <- sparseglm(X,tied,family="cox",method="lasso",lambda=0,nlambda=3,
    options=list(lasso=list(thresh=1e-12)))
  ref_tied <- survival::coxph(tied ~ X,ties="breslow",
    control=survival::coxph.control(eps=1e-11))
  expect_equal(as.numeric(coef(with_zero$fits$lasso)),unname(coef(ref_tied)),tolerance=2e-5)
  expect_equal(as.numeric(logLik(with_zero$fits$lasso)),as.numeric(logLik(ref_tied)),tolerance=1e-7)
  folds <- rep(1:2,50)
  cvfit <- sparseglm(X,Y,family=mgcv::cox.ph(),method=c("L0glm","lasso","adaptive_relaxed_lasso"),
    tune="cv",foldid=folds,nlambda=4,lambda=list(L0glm=c(20,6,2)))
  expect_true(all(cvfit$timings$status == "ok"))
  cv <- cvfit$tuning$lasso$cv
  train <- folds != 1
  model <- glmnet::glmnet(X[train,],Y[train,],family="cox",cox.ties="breslow",
    lambda=cvfit$tuning$lasso$candidates$parameter)
  beta <- as.numeric(coef(model)[,2])
  partial <- function(X,y) {
    eta <- as.numeric(X %*% beta)
    event <- which(y[,2] == 1)
    sum(vapply(event,function(i) eta[i] - log(sum(exp(eta[y[,1] >= y[i,1]]))),numeric(1)))
  }
  expected <- -2*(partial(X,Y)-partial(X[train,],Y[train,]))/sum(Y[!train,2])
  expect_equal(cv$deviance[2,1],expected,tolerance=1e-10)
  cv <- cvfit$tuning$L0glm$cv
  expect_equal(cv$parameter_train,outer(cvfit$tuning$L0glm$candidates$parameter,
    1-cv$fold_weights/sum(cv$fold_weights)))
  expect_error(sparseglm(X,Y,family="cox",intercept=TRUE),"not identifiable")
  expect_warning(unsupported <- sparseglm(X,Y,family="cox",method="L0Learn"),"unavailable")
  expect_equal(L0glm.metrics(unsupported)$status,"unsupported")
  if (requireNamespace("ncvreg",quietly=TRUE)) {
    fits <- sparseglm(X,Y,family="cox",method=c("MCP","SCAD"),nlambda=4)
    expect_true(all(fits$timings$status == "ok"))
    expect_true(all(is.finite(L0glm.metrics(fits)$loglik)))
    tied <- Y; tied[,1] <- round(tied[,1],1)
    expect_warning(no <- sparseglm(X,tied,family="cox",method="MCP",nlambda=3),"unavailable")
    expect_equal(L0glm.metrics(no)$status,"unsupported")
  }
})
