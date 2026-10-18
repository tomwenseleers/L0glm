sparseglm_fixture <- function(m=1L, sparse=FALSE) {
  set.seed(2401)
  X <- matrix(rnorm(576), 72, 8)
  B <- matrix(0, 8, m)
  B[1:2, ] <- rep(c(2,-1), m)
  if (m > 1L) B[2, 2] <- 0
  y <- X %*% B + matrix(rnorm(72*m, sd=.5), 72, m)
  list(X=if (sparse) Matrix::Matrix(X, sparse=TRUE) else X, y=y, beta_true=B)
}

test_that("all requested single-response methods expose selected models, metrics and timings", {
  sim <- sparseglm_fixture()
  methods <- c("L0glm","abess","L0Learn","lasso","relaxed_lasso","adaptive_lasso","adaptive_relaxed_lasso")
  if (requireNamespace("ordinis",quietly=TRUE)) methods <- append(methods,c("MCP","SCAD"),after=3)
  result <- sparseglm(sim, method=methods, nlambda=5, support.size=0:4,
    lambda=list(L0glm=c(20,6,2)))
  expect_s3_class(result,"sparseglm")
  expect_identical(names(result$fits),methods)
  tab <- L0glm.metrics(result)
  expect_equal(nrow(tab),length(methods))
  expect_true(all(tab$status == "ok"))
  expect_true(all(tab$seconds >= 0 & tab$path_seconds >= 0 & tab$cv_seconds == 0))
  expect_true(all(tab$seconds + 1e-8 >= tab$path_seconds + tab$cv_seconds))
  expect_true(all(c("MCC","within_row_MCC","gic","loglik","note") %in% names(tab)))
  for (method in methods) {
    fit <- result$fits[[method]]
    expect_s4_class(coef(fit),"sparseMatrix")
    expect_equal(dim(coef(fit)),c(9L,1L))
    scores <- result$tuning[[method]]$candidates
    expect_equal(result$tuning[[method]]$index,which.min(scores$score))
    expect_equal(scores$score[result$tuning[[method]]$index],L0glm.metrics(fit,mu_floor=result$mu_floor)$gic)
    expect_equal(tab$loglik[tab$method == method],as.numeric(logLik(fit)))
    expect_equal(coef(result)[[method]],coef(fit))
  }
  expect_equal(sparseglm(sim,method="abess",support.size=c(0,2,4),nlambda=3)$tuning$abess$candidates$parameter,c(0,2,4))
})

test_that("every supported information criterion selects the common candidate score", {
  sim <- sparseglm_fixture()
  for (ic in c("aic","bic","gic","ebic","mbic")) {
    result <- sparseglm(sim,method="lasso",tune=toupper(ic),nlambda=8)
    path <- result$tuning$lasso
    expect_equal(path$best,path$candidates$parameter[which.min(path$candidates$score)])
    expect_equal(min(path$candidates$score),L0glm.metrics(result$fits$lasso)[[ic]])
    expect_equal(L0glm.metrics(result)[[ic]],min(path$candidates$score))
  }
})

test_that("CV uses shared row folds and exact fold-local adaptive pilots", {
  sim <- sparseglm_fixture(2L)
  folds <- rep(1:3,length.out=72)
  result <- sparseglm(sim,method=c("L0glm","lasso","adaptive_lasso","adaptive_relaxed_lasso","L0Learn"),
    grouped=FALSE,tune="cv",foldid=folds,nlambda=5,support.size=0:8,
    lambda=list(L0glm=c(20,6,2)))
  expect_true(all(result$timings$status == "ok"))
  expect_identical(result$foldid,folds)
  for (name in names(result$tuning)) {
    cv <- result$tuning[[name]]$cv
    expect_identical(cv$foldid,folds)
    expect_equal(cv$nfolds,3)
    expect_equal(cv$mean,as.numeric(cv$deviance %*% (cv$fold_weights/sum(cv$fold_weights))))
    if (name == "L0Learn") expect_true(all(is.finite(cv$mean)))
  }
  for (k in 1:3) {
    train <- folds != k
    pilot <- L0glm.fit(sim$X[train,],sim$y[train,],start=NULL,lambda="gic",maxit=1,
      grouped=FALSE,relaxed=FALSE,thresh=0,intercept=TRUE,parallelize=FALSE)
    expected <- 1/pmax(abs(as.matrix(coef(pilot))[-1,,drop=FALSE]),1e-8)
    for (method in c("adaptive_lasso","adaptive_relaxed_lasso")) {
      saved <- result$tuning[[method]]$cv$adaptive[[k]]
      expect_equal(saved$coefficients,coef(pilot))
      expect_equal(saved$penalty.factor,expected)
    }
  }
  cv <- result$tuning$L0glm$cv
  expect_equal(cv$parameter_train,outer(result$tuning$L0glm$candidates$parameter,
    1-cv$fold_weights/sum(cv$fold_weights)))
  # Direct held-out likelihood check for a nonadaptive LASSO fold.
  cv <- result$tuning$lasso$cv
  keep <- folds != 1L
  predictions <- vapply(1:2,function(k) {
    ref <- glmnet::glmnet(sim$X[keep,],sim$y[keep,k],family="gaussian",
      lambda=result$tuning$lasso$candidates$parameter,intercept=TRUE)
    as.numeric(cbind(1,sim$X[!keep,]) %*% coef(ref)[,2])
  },numeric(sum(!keep)))
  expected <- mean((sim$y[!keep,,drop=FALSE]-predictions)^2)
  expect_equal(cv$deviance[2,1],expected)
})

test_that("multichannel grouping preserves matrix weights and signed pilot rows", {
  sim <- sparseglm_fixture(3L,TRUE)
  sim$beta_true[1,] <- c(2,-2,0)
  sim$y <- as.matrix(sim$X %*% sim$beta_true) + matrix(rnorm(216,sd=.4),72,3)
  grouped <- sparseglm(sim,method=c("lasso","adaptive_lasso"),nlambda=5)
  expect_true(all(grouped$timings$scope == "grouped"))
  pilot <- grouped$tuning$adaptive_lasso$adaptive
  B <- as.matrix(pilot$coefficients)[-1,,drop=FALSE]
  expect_equal(pilot$penalty.factor,1/pmax(sqrt(rowMeans(B^2)),1e-8))
  expect_true(any(L0glm.metrics(grouped)$within_row_FP > 0))
  abess_grouped <- sparseglm(sim,method="abess",tune="cv",nfolds=3,
    nlambda=3,support.size=c(0,1,3))
  expect_equal(abess_grouped$timings$status,"ok")
  expect_equal(abess_grouped$timings$scope,"grouped")
  expect_true(all(is.finite(abess_grouped$tuning$abess$cv$mean)))
  expect_null(colnames(abess_grouped$data$y))
  W <- matrix(runif(216,.5,2),72,3)
  independent <- sparseglm(sim,method="lasso",weights=W,nlambda=5)
  expect_equal(independent$timings$scope,"independent")
  expect_equal(independent$fits$lasso$prior.weights,W)
  expected <- L0glm:::get.metrics(coef(independent$fits$lasso),cbind(1,sim$X),sim$y,
    weights=W,family=gaussian(),always.include=1L)$loglik
  expect_equal(L0glm.metrics(independent)$loglik,expected)
  # Support-size paths remain integer grids when ABESS fits channels separately.
  abess <- sparseglm(sim,method="abess",weights=W,nlambda=5,support.size=c(0,1,3,5))
  expect_equal(abess$timings$scope,"independent")
  expect_equal(abess$tuning$abess$candidates$parameter,c(0,1,3,5))
  expect_true(all(is.finite(abess$tuning$abess$candidates$score)))
  expect_warning(unsupported <- sparseglm(sim,method="lasso",weights=W,grouped=TRUE,nlambda=3),"unavailable")
  tab <- L0glm.metrics(unsupported)
  expect_equal(tab$status,"unsupported")
  expect_true(is.na(tab$gic))
  expect_match(tab$note,"grouped")
})

test_that("identity Poisson approximations are explicit and scored as Poisson", {
  set.seed(2402)
  X <- matrix(runif(480),60,8); B <- c(3,2,rep(0,6))
  y <- rpois(60,1+X %*% B)
  expect_warning(no <- sparseglm(X,y,family=poisson(identity),method="abess",
    offset=1,nonnegative=TRUE,support.size=0:3),"unavailable")
  expect_equal(L0glm.metrics(no)$status,"unsupported")
  methods <- c("L0Learn","abess","lasso")
  if (requireNamespace("ordinis",quietly=TRUE)) methods <- c(methods,"MCP","SCAD")
  result <- sparseglm(X,y,family=poisson(identity),method=methods,offset=1,
    nonnegative=TRUE,approximation="wls",support.size=0:3,nlambda=5,beta_true=B)
  expect_true(all(result$timings$status == "ok"))
  for (name in methods) {
    B <- as.numeric(coef(result$fits[[name]]))
    expect_true(all(B >= 0))
    expected <- sum(dpois(y,1+as.numeric(X %*% B),log=TRUE))
    expect_equal(L0glm.metrics(result)$loglik[result$timings$method == name],expected)
    expect_equal(result$timings$approximation[result$timings$method == name],
      if (name == "lasso") "none" else "identity_poisson_wls")
  }
})

test_that("input validation, RNG preservation and backend failures are explicit", {
  sim <- sparseglm_fixture()
  set.seed(18); before <- .Random.seed
  result <- sparseglm(sim,method="lasso",tune="cv",nlambda=3,nfolds=3,seed=5)
  expect_identical(.Random.seed,before)
  other <- sparseglm(sim$X,sim$y,method="lasso",tune="cv",nlambda=3,nfolds=3,seed=5)
  expect_identical(result$foldid,other$foldid)
  expect_identical(coef(result),coef(other))
  expect_error(sparseglm(sim,method="wrong"),"Unknown method")
  expect_error(sparseglm(sim,method=1),"method")
  expect_error(sparseglm(sim,nfolds=1,tune="cv"),"nfolds")
  expect_error(sparseglm(sim,foldid=rep(1:2,36)),"only used")
  expect_error(sparseglm(sim,lambda=c(1,NA)),"lambda")
  expect_error(sparseglm(sim,support.size=9),"support.size")
  expect_error(sparseglm(sim,weights=rep(-1,72)),"weights")
  expect_error(sparseglm(sim$X[-1,],sim$y),"aligned")
  expect_error(sparseglm(cbind(1,sim$X),sim$y),"intercept column")
  expect_warning(failed <- sparseglm(sim,method=c("L0glm","lasso"),nlambda=3,
    options=list(L0glm=list(solver="invalid"))),"failed")
  expect_equal(failed$timings$status,c("error","ok"))
  expect_true(is.na(L0glm.metrics(failed)$gic[1]))
  expect_true(is.finite(L0glm.metrics(failed)$gic[2]))
})
