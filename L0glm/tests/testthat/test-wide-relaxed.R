test_that("wide Gaussian relaxation agrees with native glmnet on each support", {
  skip_if_not_installed("glmnet")
  set.seed(610)
  X <- matrix(rnorm(60*110), 60, 110)
  weights <- runif(60, .5, 2)
  for (m in 1:2) for (sparse in c(FALSE, TRUE)) for (nonnegative in c(FALSE, TRUE)) {
    A <- if (sparse) Matrix::Matrix(X, sparse=TRUE) else X
    B <- matrix(0,110,m); B[1:3,] <- c(1.5,.7,.3)
    offset <- matrix(runif(60*m,-.2,.2),60,m)
    y <- A %*% B + offset + matrix(rnorm(60*m,sd=.3),60,m)+.4
    grid <- c(3,1,.3,.1)
    intercept <- !nonnegative
    native <- do.call(glmnet::glmnet,list(x=A,y=if(m==1L) as.numeric(y-offset) else y-offset,
      family=if(m==1L)"gaussian" else "mgaussian",weights=weights,
      lambda=grid,intercept=intercept,relax=TRUE,standardize.response=FALSE,
      lower.limits=if(nonnegative)0 else -Inf))$relaxed
    result <- sparseglm(A,y,weights=weights,offset=offset,method="relaxed_lasso",
      grouped=m>1L,intercept=intercept,nonnegative=nonnegative,lambda=grid,tune="bic")
    expect_equal(result$timings$status,"ok")
    j <- result$tuning$relaxed_lasso$index
    expected <- if(m==1L) native$beta[,j,drop=FALSE] else
      do.call(cbind,lapply(native$beta,function(z)z[,j]))
    if(intercept) expected <- rbind(if(m==1L) native$a0[j] else native$a0[,j],expected)
    expect_equal(as.matrix(coef(result$fits$relaxed_lasso)),as.matrix(expected),
      tolerance=1e-5,ignore_attr=TRUE)
    expect_equal(result$tuning$relaxed_lasso$candidates$parameter,grid)
  }
})

test_that("wide refits use only selected columns and preserve adaptive penalties", {
  skip_if_not_installed("glmnet")
  set.seed(611)
  X <- matrix(rnorm(60*140),60,140)
  y <- 2*X[,1]+rnorm(60,sd=.2)
  native <- glmnet::glmnet
  refit_widths <- integer()
  testthat::local_mocked_bindings(glmnet=function(x,y,...) {
    args <- list(...)
    if(identical(args$lambda,0)) refit_widths <<- c(refit_widths,ncol(x))
    native(x,y,...)
  },.package="glmnet")
  fit <- sparseglm(X,y,method="adaptive_relaxed_lasso",tune="cv",nfolds=3L,
    foldid=rep(1:3,20),nlambda=5L,intercept=FALSE)
  expect_equal(fit$timings$status,"ok")
  expect_gt(length(refit_widths),0L)
  expect_true(all(refit_widths<ncol(X)))
  expect_equal(dim(coef(fit$fits$adaptive_relaxed_lasso)),c(140L,1L))
  expect_true(any(is.finite(fit$tuning$adaptive_relaxed_lasso$cv$mean)))
})

test_that("shared Gaussian CV matrices agree with explicit training subsets", {
  skip_if_not_installed("glmnet")
  set.seed(612)
  X <- matrix(rnorm(60*90),60,90)
  folds <- rep(1:3,20)
  weights <- runif(60,.5,2)
  offset <- runif(60,-.3,.3)
  y <- .4+1.5*X[,1]+rnorm(60,sd=.5)+offset
  grid <- c(3,1,.3)
  fit <- sparseglm(X,y,method=c("lasso","relaxed_lasso"),tune="cv",
    lambda=grid,foldid=folds,weights=weights,offset=offset,intercept=TRUE)
  expect_true(all(fit$timings$status=="ok"))
  for(method in names(fit$fits)) for(k in 1:3) {
    train <- folds!=k
    native <- do.call(glmnet::glmnet,list(x=X[train,,drop=FALSE],
      y=(y-offset)[train],weights=weights[train],intercept=TRUE,
      lambda=grid,relax=method=="relaxed_lasso"))
    if(method=="relaxed_lasso") native <- native$relaxed
    pred <- as.matrix(cbind(1,X[!train,,drop=FALSE]) %*% rbind(native$a0,native$beta))+offset[!train]
    expected <- colSums(weights[!train]*(pred-y[!train])^2)/sum(weights[!train])
    expect_equal(fit$tuning[[method]]$cv$deviance[,k],expected,tolerance=1e-5,
      ignore_attr=TRUE)
  }
})
