screening_fixture <- function(m=1L, sparse=FALSE) {
  set.seed(8401)
  X <- matrix(rnorm(80*24),80,24,dimnames=list(NULL,paste0("x",1:24)))
  B <- matrix(0,24,m); B[1:3,] <- c(3,-2,1.5)
  Y <- X%*%B+matrix(rnorm(80*m,sd=.4),80,m)
  list(X=if (sparse) Matrix::Matrix(X,sparse=TRUE) else X,y=Y,beta_true=B)
}

test_that("every installed adapter supports shared SIS and ISIS screening", {
  sim <- screening_fixture()
  for (mode in c("SIS","ISIS")) {
    result <- sparseglm(sim,method="all",screen=mode,nsis=8,screen.maxit=3,
      nlambda=4,support.size=0:5,lambda=list(L0glm=c(10,5,2)))
    expect_true(all(result$timings$status=="ok"))
    expect_true(all(result$timings$screen==mode))
    for (name in names(result$fits)) {
      B <- coef(result$fits[[name]]); diagnostics <- result$tuning[[name]]$screening
      expect_s4_class(B,"sparseMatrix")
      expect_equal(dim(B),c(25L,1L))
      expect_identical(rownames(B),c("Intercept",colnames(sim$X)))
      expect_true(all(B[-c(1L,diagnostics$candidates+1L),,drop=FALSE]==0))
      expect_equal(result$fits[[name]]$p,24)
      expect_equal(min(result$tuning[[name]]$candidates$score),L0glm.metrics(result$fits[[name]])$gic)
      expect_equal(diagnostics$history$score[diagnostics$best_iteration],
        min(diagnostics$history$score))
    }
  }
  one <- sparseglm(sim,method="all",screen="SIS",nsis=1,nlambda=3,
    support.size=0:1,lambda=list(L0glm=5))
  expect_true(all(one$timings$status=="ok"))
  expect_true(all(vapply(coef(one),function(B) sum(B[-1L,]!=0)<=1L,logical(1))))
})

test_that("CV screens using only training rows and restores full coefficients", {
  sim <- screening_fixture()
  folds <- rep(1:2,40)
  for (mode in c("SIS","ISIS")) {
    result <- sparseglm(sim,method=c("lasso","adaptive_lasso"),screen=mode,
      nsis=6,screen.maxit=3,tune="cv",foldid=folds,nlambda=4)
    expect_true(all(result$timings$status=="ok"))
    for (name in names(result$fits)) for (k in 1:2) {
      train <- folds!=k
      reference <- sparseglm(sim$X[train,],sim$y[train,,drop=FALSE],method=name,
        screen=mode,nsis=6,screen.maxit=3,nlambda=4,
        lambda=result$tuning[[name]]$candidates$parameter)
      actual <- result$tuning[[name]]$cv$screening[[k]]
      expect_equal(actual$candidates,reference$tuning[[name]]$screening$candidates)
      expect_equal(actual$history,reference$tuning[[name]]$screening$history)
      d <- result$data
      d$X <- sim$X[train,,drop=FALSE]; d$design <- cbind(Intercept=1,d$X)
      for (field in c("y","weights","offset")) d[[field]] <- d[[field]][train,,drop=FALSE]
      path <- L0glm_benchmark_path(name,d,result$tuning[[name]]$candidates$parameter)
      for (j in seq_along(path$grid)) {
        fit <- L0glm_benchmark_fit(path$coefficients[[j]],d)
        loss <- L0glm.metrics(fit,newdata=sim$X[!train,,drop=FALSE],
          y=sim$y[!train,,drop=FALSE])$deviance/sum(!train)
        expect_equal(result$tuning[[name]]$cv$deviance[j,k],loss)
      }
    }
  }
})

test_that("multichannel screening keeps grouping, sparse designs and channel dimensions", {
  sim <- screening_fixture(3L,TRUE)
  for (grouped in c(TRUE,FALSE)) {
    result <- sparseglm(sim,method=c("L0glm","lasso","adaptive_relaxed_lasso"),
      grouped=grouped,screen="ISIS",nsis=6,screen.maxit=3,nlambda=3,
      lambda=list(L0glm=c(10,5,2)))
    expect_true(all(result$timings$status=="ok"))
    expect_true(all(result$timings$scope==if (grouped) "grouped" else "independent"))
    expect_true(all(vapply(coef(result),function(B) identical(dim(B),c(25L,3L)),logical(1))))
    expect_s4_class(result$data$X,"sparseMatrix")
    expect_equal(dim(result$tuning$adaptive_relaxed_lasso$adaptive$coefficients),c(25L,3L))
  }
  sim <- simul(method=4,n=80,p=80,s=.1,seed=2)
  channels <- order(colSums(sim$y),decreasing=TRUE)[1:3]
  for (name in c("y","beta_true")) sim[[name]] <- sim[[name]][,channels,drop=FALSE]
  result <- sparseglm(sim,family=poisson(identity),method=c("L0glm","abess","lasso"),
    nonnegative=TRUE,intercept=FALSE,screen="ISIS",nsis=12,screen.maxit=2,
    approximation="wls_grouped",nlambda=3,support.size=0:10,lambda=list(L0glm=c(10,5,2)))
  expect_true(all(result$timings$status=="ok"))
  expect_true(all(result$timings$scope=="grouped"))
  expect_true(all(vapply(coef(result),function(B) identical(dim(B),c(80L,3L)),logical(1))))
})

test_that("screening handles joint multinomial and Cox models", {
  skip_if_not_installed("mgcv"); skip_if_not_installed("survival")
  set.seed(8402)
  X <- matrix(rnorm(100*12),100,12)
  P <- L0glm_ext_softmax(cbind(X[,1],-X[,1],0))$p
  Y <- t(vapply(1:100,function(i) rmultinom(1,3,P[i,]),numeric(3)))
  event <- rexp(100,exp(X[,1])); censor <- rexp(100,.3)
  for (mode in c("SIS","ISIS")) {
    for (kind in c("multinomial","cox")) {
      y <- if (kind=="multinomial") Y else survival::Surv(pmin(event,censor),event<=censor)
      result <- sparseglm(X,y,family=kind,method="all",screen=mode,nsis=5,
        screen.maxit=2,tune="cv",nfolds=2,nlambda=3,support.size=0:3,
        lambda=list(L0glm=c(10,5,2)))
      expect_true(all(result$timings$status=="ok"))
      for (name in names(result$fits)) {
        expect_true(all(is.finite(result$tuning[[name]]$cv$mean)))
        expect_equal(nrow(coef(result$fits[[name]])),12L+(kind=="multinomial"))
        if (kind=="multinomial") expect_equal(as.numeric(Matrix::rowSums(coef(result$fits[[name]]))),rep(0,13),tolerance=1e-12)
      }
    }
  }
})

test_that("automatic screening and zero starts respect eligible and mandatory predictors", {
  set.seed(8403)
  X <- matrix(rnorm(30*301),30,301); y <- X[,1]*3+rnorm(30)
  result <- sparseglm(X,y,method="lasso",screen="auto",nsis=6,nlambda=3,screen.maxit=2)
  expect_identical(result$timings$screen,"ISIS")
  none <- sparseglm(X[,1:300],y,method="lasso",screen="auto",nlambda=3)
  expect_identical(none$timings$screen,"none")
  cv <- sparseglm(X[,1:300],y,method="lasso",screen="auto",nlambda=3,
    tune="cv",nfolds=2,nsis=5,screen.maxit=2)
  expect_identical(cv$timings$screen,"none")
  expect_true(all(vapply(cv$tuning$lasso$cv$screening,function(z)
    identical(z$method,"ISIS"),logical(1))))
  start <- c(1,rep(0,23))
  sim <- screening_fixture()
  fit <- sparseglm(sim,method="L0glm",screen="SIS",nsis=3,lambda=5,
    options=list(L0glm=list(start=start,always.include=2L)))
  expect_equal(fit$timings$status,"ok")
  expect_true(all(coef(fit)$L0glm[4:25,]==0))
  expect_true(2L %in% fit$tuning$L0glm$screening$candidates)
  expect_error(sparseglm(sim,screen="invalid"),"screen")
  expect_error(sparseglm(sim,nsis=0),"nsis")
  expect_error(sparseglm(sim,screen.maxit=0),"screen.maxit")
})
