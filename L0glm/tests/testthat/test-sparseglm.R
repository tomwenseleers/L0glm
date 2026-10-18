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
      "identity_poisson_wls")
  }
})

test_that("Poisson deconvolution evaluates ten finite glmnet candidates without an intercept", {
  sim <- simul(method=2,n=120,p=120,s=.1,seed=2)
  expect_true(any(sim$y == 0))
  methods <- c("lasso","relaxed_lasso","adaptive_lasso","adaptive_relaxed_lasso")
  result <- sparseglm(sim,family=poisson(identity),method=methods,
    approximation="wls",nonnegative=TRUE,intercept=FALSE,nlambda=10)
  expect_true(all(result$timings$status == "ok"))
  expect_true(all(result$timings$approximation == "identity_poisson_wls"))
  for (name in methods) {
    candidates <- result$tuning[[name]]$candidates
    expect_equal(nrow(candidates),10L)
    expect_equal(candidates$support_size[1],0)
    expect_true(all(is.finite(candidates$score)))
    expect_true(all(is.finite(candidates$parameter) & candidates$parameter > 0))
    expect_true(all(diff(candidates$parameter) < 0))
    expect_equal(diff(log(candidates$parameter)),rep(log(1e-4)/9,9))
    expect_equal(dim(coef(result)[[name]]),c(120L,1L))
    expect_true(all(coef(result)[[name]] >= 0))
  }
  # An explicit user grid, including zero, must be used without modification.
  specified <- sparseglm(sim,family=poisson(identity),method="lasso",
    approximation="wls",nonnegative=TRUE,intercept=FALSE,nlambda=10,lambda=c(1,.1,0))
  expect_equal(specified$tuning$lasso$candidates$parameter,c(1,.1,0))
  # ABESS must search beyond its own conservative default support ceiling.
  abess <- sparseglm(sim,family=poisson(identity),method="abess",approximation="wls",
    nonnegative=TRUE,intercept=FALSE,support.size=c(0,12,24),nlambda=10)
  expect_equal(abess$timings$status,"ok")
  expect_equal(abess$tuning$abess$candidates$parameter,c(0,12,24))
  expect_true(all(is.finite(abess$tuning$abess$candidates$score)))
})

test_that("all expands only to installed methods supporting the supplied arguments", {
  sim <- sparseglm_fixture()
  methods <- L0glm_benchmark_methods("all")
  for (nonnegative in c(FALSE,TRUE)) {
    result <- sparseglm(sim,method="all",nonnegative=nonnegative,intercept=TRUE,
      nlambda=3,support.size=0:2,lambda=list(L0glm=6))
    packages <- c(L0glm="L0glm",abess="abess",L0Learn="L0Learn",MCP="ordinis",SCAD="ordinis",
      lasso="glmnet",relaxed_lasso="glmnet",adaptive_lasso="glmnet",adaptive_relaxed_lasso="glmnet")
    expected <- methods[vapply(packages[methods],requireNamespace,logical(1),quietly=TRUE)]
    if (nonnegative) expected <- intersect(expected,c("L0glm","L0Learn"))
    expect_identical(names(result$fits),expected)
    expect_true(all(result$timings$status == "ok"))
  }
  sim <- simul(method=2,n=120,p=120,s=.1,seed=2)
  exact <- sparseglm(sim,family=poisson(identity),method="all",nonnegative=TRUE,
    intercept=FALSE,nlambda=3,lambda=list(L0glm=6))
  expect_identical(names(exact$fits),"L0glm")
  expect_identical(exact$timings$status,"ok")
  multi <- sparseglm_fixture(2)
  group <- sparseglm(multi,method="all",grouped=TRUE,weights=matrix(seq(.5,2,length.out=144),72,2),
    nlambda=3,support.size=0:2,lambda=list(L0glm=6))
  expect_identical(names(group$fits),"L0glm")
  expect_identical(group$timings$status,"ok")
  binary <- sparseglm_fixture()
  binary$y <- as.numeric(rep(c(1,0,0),24)); binary$beta_true <- NULL
  probit <- sparseglm(binary,family=binomial("probit"),method="all",intercept=TRUE,
    nonnegative=TRUE,nlambda=3,lambda=list(L0glm=6))
  expect_identical(names(probit$fits),c("L0glm","lasso","relaxed_lasso","adaptive_lasso","adaptive_relaxed_lasso"))
  expect_true(all(probit$timings$status == "ok"))
  testthat::local_mocked_bindings(L0glm_benchmark_require=function(package) {
    L0glm_benchmark_unavailable("Optional test backend is not installed.","missing_package")
  }, .package="L0glm")
  installed <- sparseglm(sim,family=poisson(identity),method="all",nonnegative=TRUE,
    intercept=FALSE,approximation="wls",nlambda=3,lambda=list(L0glm=6))
  expect_identical(names(installed$fits),"L0glm")
  # Explicit requests retain unsupported rows and their diagnostic.
  expect_warning(explicit <- sparseglm(sim,family=poisson(identity),method="lasso",
    nonnegative=TRUE,intercept=FALSE,nlambda=3),"unavailable")
  expect_identical(explicit$timings$status,"unsupported")
  expect_match(explicit$timings$note,"null mean")
})

test_that("multichannel Poisson WLS paths retain channel weights and valid CV scores", {
  sim <- simul(method=4,n=100,p=100,s=.1,seed=2)
  channels <- order(colSums(sim$y),decreasing=TRUE)[1:2]
  for (name in c("y","y_true","beta_true")) sim[[name]] <- sim[[name]][,channels,drop=FALSE]
  methods <- c("lasso","relaxed_lasso","adaptive_lasso","adaptive_relaxed_lasso")
  folds <- rep(1:2,50)
  result <- sparseglm(sim,family=poisson(identity),method=methods,
    nonnegative=TRUE,intercept=FALSE,approximation="wls",nlambda=10,
    tune="cv",foldid=folds)
  expect_true(all(result$timings$status == "ok"))
  expect_true(all(result$timings$scope == "independent"))
  for (name in methods) {
    expect_equal(dim(coef(result)[[name]]),c(100L,2L))
    expect_equal(nrow(result$tuning[[name]]$candidates),10L)
    expect_true(all(is.finite(result$tuning[[name]]$cv$mean)))
    expect_equal(result$fits[[name]]$prior.weights,matrix(1,100,2))
  }
  keep <- folds != 1L
  grid <- result$tuning$lasso$candidates$parameter
  predictions <- vapply(1:2,function(k) {
    ref <- glmnet::glmnet(sim$X[keep,],sim$y[keep,k],family="gaussian",
      weights=1/(sim$y[keep,k]+.1),lambda=grid,intercept=FALSE,lower.limits=0)
    as.numeric(sim$X[!keep,] %*% as.numeric(coef(ref)[-1,5]))
  },numeric(sum(!keep)))
  expected <- sum(poisson()$dev.resids(sim$y[!keep,,drop=FALSE],pmax(predictions,1e-8),1))/length(predictions)
  expect_equal(result$tuning$lasso$cv$deviance[5,1],expected)
})

test_that("shared-row Poisson WLS enables native grouped ABESS and glmnet comparisons", {
  sim <- simul(method=4,n=100,p=100,s=.1,seed=2)
  channels <- order(colSums(sim$y),decreasing=TRUE)[1:2]
  for (name in c("y","y_true","beta_true")) sim[[name]] <- sim[[name]][,channels,drop=FALSE]
  glmnet_methods <- c("lasso","relaxed_lasso","adaptive_lasso","adaptive_relaxed_lasso")
  result <- sparseglm(sim,family=poisson(identity),method=c("abess",glmnet_methods),
    nonnegative=TRUE,intercept=FALSE,approximation="wls_grouped",nlambda=10,support.size=0:15)
  expect_true(all(result$timings$status == "ok"))
  expect_true(all(result$timings$scope == "grouped"))
  expect_true(all(result$timings$approximation == "identity_poisson_row_wls"))
  row_w <- 1/(rowMeans(sim$y)+.1)
  setup <- L0glm_benchmark_setup("lasso",result$data)
  expect_equal(setup$data$weights,matrix(row_w,100,2))
  expect_equal(result$fits$lasso$prior.weights,matrix(1,100,2))
  for (method in glmnet_methods) expect_equal(nrow(result$tuning[[method]]$candidates),10L)
  reference <- glmnet::glmnet(sim$X,sim$y,family="mgaussian",weights=row_w,
    lambda=result$tuning$lasso$candidates$parameter,intercept=FALSE,lower.limits=0,
    standardize.response=FALSE)
  index <- result$tuning$lasso$index
  expected <- do.call(cbind,lapply(reference$beta,function(B)B[,index,drop=FALSE]))
  expect_equal(unname(as.matrix(coef(result)$lasso)),unname(as.matrix(expected)))
  expect_warning(unavailable <- sparseglm(sim,family=poisson(identity),method="abess",
    nonnegative=TRUE,intercept=FALSE,approximation="wls",support.size=0:2),"unavailable")
  expect_equal(unavailable$timings$status,"unsupported")
  expect_match(unavailable$timings$note,"channel-specific weights")
  auto <- sparseglm(sim,family=poisson(identity),method="all",approximation="wls",
    nonnegative=TRUE,intercept=FALSE,nlambda=3,support.size=0:2,lambda=list(L0glm=6))
  expect_false("abess" %in% names(auto$fits))
  expect_true(all(auto$timings$status == "ok"))
  independent <- sparseglm(sim,family=poisson(identity),method="abess",approximation="wls",
    grouped=FALSE,nonnegative=TRUE,intercept=FALSE,support.size=0:2)
  expect_equal(independent$timings$status,"ok")
  expect_equal(independent$timings$scope,"independent")
  # Fold-local shared row weights must still be scored as held-out Poisson.
  folds <- rep(1:2,50)
  cv <- sparseglm(sim,family=poisson(identity),method="lasso",approximation="wls_grouped",
    nonnegative=TRUE,intercept=FALSE,nlambda=4,tune="cv",foldid=folds)
  expect_equal(cv$timings$scope,"grouped")
  train <- folds != 1L
  reference <- glmnet::glmnet(sim$X[train,],sim$y[train,,drop=FALSE],family="mgaussian",
    weights=row_w[train],lambda=cv$tuning$lasso$candidates$parameter,
    intercept=FALSE,lower.limits=0,standardize.response=FALSE)
  B <- do.call(cbind,lapply(reference$beta,function(B)B[,2,drop=FALSE]))
  mu <- as.matrix(sim$X[!train,] %*% B)
  expected <- sum(poisson()$dev.resids(sim$y[!train,,drop=FALSE],pmax(mu,1e-8),1))/length(mu)
  expect_equal(cv$tuning$lasso$cv$deviance[2,1],expected)
})

test_that("ABESS Gaussian row weighting agrees with weighted least squares on a fixed support", {
  sim <- sparseglm_fixture()
  W <- seq(.2,2,length.out=72)
  for (intercept in c(FALSE,TRUE)) {
    fit <- sparseglm(sim,method="abess",weights=W,intercept=intercept,support.size=8)
    X <- if (intercept) cbind(1,sim$X) else sim$X
    expected <- stats::lm.wfit(X,as.numeric(sim$y),w=W)$coefficients
    expect_equal(fit$timings$status,"ok")
    expect_equal(as.numeric(coef(fit)$abess),unname(expected),tolerance=1e-8)
    expect_equal(fit$tuning$abess$candidates$parameter,8)
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
