model_summary_fixture <- function(m=1L) {
  set.seed(2501)
  X <- matrix(rnorm(320),64,5,dimnames=list(NULL,letters[1:5]))
  B <- matrix(rep(c(2,-1,0,0,0),m),5,m)
  y <- X %*% B + matrix(rnorm(64*m,sd=.5),64,m)
  fit <- L0glm.fit(X,y,lambda=0,thresh=0,intercept=TRUE,
    parallelize=FALSE,relaxed=FALSE)
  list(fit=fit,X=X,y=y,truth=rbind(0,B))
}

test_that("fit summaries retain sparse original-scale coefficients and common scores", {
  data <- model_summary_fixture(3L)
  fit <- data$fit
  details <- summary(fit,beta_true=data$truth)
  expect_s3_class(details,"summary.L0glm_fit")
  expect_identical(details$coefficients,coef(fit))
  expect_identical(coef(details),coef(fit))
  expect_s4_class(details$coefficients,"sparseMatrix")
  expect_equal(details$metrics,L0glm.metrics(fit,mu_floor=0,beta_true=data$truth))
  expect_equal(c(details$n,details$p,details$m),c(64,5,3))
  expect_equal(details$nobs,nobs(fit))
  expect_equal(details$selected_rows,setdiff(which(Matrix::rowSums(coef(fit)!=0)>0),fit$always.include))
  expect_equal(details$nonzero,sum(coef(fit)[-fit$always.include,,drop=FALSE]!=0))
  expect_identical(details$dispersion,fit$dispersion)
  expect_identical(details$dispersion.estimated,fit$dispersion.estimated)
  output <- capture.output(visible <- withVisible(print(details,max.rows=2,max.cols=1)))
  expect_false(visible$visible)
  expect_identical(visible$value,details)
  expect_true(any(grepl("Showing 2 of",output,fixed=TRUE)))
  expect_true(any(grepl("Recovery metrics",output,fixed=TRUE)))
  expect_equal(dim(details$coefficients),c(6L,3L))
  expect_error(print(details,max.rows=0),"positive integers")
  expect_error(print(details,max.cols=1.5),"positive integers")
  expect_output(print(details,max.rows=Inf,max.cols=Inf),"Coefficient estimates")
})

test_that("printing fits is compact and never evaluates likelihoods", {
  fit <- model_summary_fixture()$fit
  local_mocked_bindings(L0glm.metrics=function(...) stop("Scoring must not run"))
  output <- capture.output(visible <- withVisible(print(fit,digits=3)))
  expect_false(visible$visible)
  expect_identical(visible$value,fit)
  expect_lt(length(output),16)
  expect_true(any(grepl("Selected rows",output,fixed=TRUE)))
  expect_false(any(grepl("prior.weights|^\\$X|^\\$y",output)))
})

test_that("convergence displays use the final pass for each stage and channel", {
  fit <- model_summary_fixture(2L)$fit
  fit$convergence <- data.frame(stage=rep(c("selection","relaxed"),c(4,2)),
    pass=c(1,1,2,2,1,1),channel=rep(1:2,3),iterations=c(25,25,7,8,5,6),
    converged=c(FALSE,FALSE,TRUE,TRUE,TRUE,FALSE))
  convergence <- summary(fit)$convergence
  expect_equal(convergence$stage,c("selection","relaxed"))
  expect_equal(convergence$channels,c(2,2))
  expect_equal(convergence$converged,c(2,1))
  expect_equal(convergence$not_converged,c(0,1))
  expect_equal(convergence$max_iterations,c(8,6))
  fit$convergence <- NULL
  expect_null(summary(fit)$convergence)
  expect_output(print(fit),"Convergence diagnostics: unavailable")
  fit$coefficients[,] <- 0
  fit$always.include <- integer()
  expect_output(print(summary(fit)),"No selected coefficients")
})

test_that("native summaries distinguish exact scores from optional positive-mean floors", {
  set.seed(2502)
  X <- matrix(runif(120),40,3)
  y <- rpois(40,1+X[,1])
  fit <- L0nnpois(X,y,lambda=2,parallelize=FALSE)
  fit$coefficients[,] <- 0
  expect_s3_class(summary(fit),"summary.L0glm_fit")
  expect_equal(summary(fit)$metrics$loglik,-Inf)
  floored <- summary(fit,mu_floor=1e-4)
  expect_equal(floored$metrics,L0glm.metrics(fit,mu_floor=1e-4))
  expect_true(is.finite(floored$metrics$loglik))
  expect_output(print(floored),"Scoring mean floor")
  expect_output(print(fit),"Family: poisson")
})

test_that("special families inherit summary methods without losing their coordinates", {
  skip_if_not_installed("mgcv"); skip_if_not_installed("survival")
  set.seed(2503)
  X <- matrix(rnorm(180),60,3)
  responses <- list(
    nb=rnbinom(60,mu=exp(.2+.3*X[,1]),size=2),
    multinomial=factor(rep(1:3,20)),
    cox=survival::Surv(rexp(60),rep(c(1,1,0),20)))
  families <- list(nb=mgcv::nb(theta=2),multinomial=mgcv::multinom(K=2),cox=mgcv::cox.ph())
  for (kind in names(families)) {
    fit <- L0glm.fit(X,responses[[kind]],family=families[[kind]],
      lambda=2,relaxed=FALSE,maxit=30)
    details <- summary(fit)
    expect_s3_class(details,"summary.L0glm_fit")
    expect_identical(details$coefficients,coef(fit))
    expect_equal(details$metrics,L0glm.metrics(fit,mu_floor=0))
    expect_output(print(details),"Fit metrics")
    if (kind == "nb") {
      expect_identical(details$theta,fit$theta)
      expect_identical(details$theta.estimated,fit$theta.estimated)
      expect_output(print(details),"NB theta")
    }
    if (kind == "multinomial") expect_equal(as.numeric(Matrix::rowSums(details$coefficients)),rep(0,nrow(coef(fit))),tolerance=1e-12)
    if (kind == "cox") expect_output(print(details),"Criterion sample size \\(Cox\\)")
  }
})

test_that("comparison summaries retain metrics, CV diagnostics and method notes", {
  data <- model_summary_fixture()
  sim <- list(X=data$X,y=data$y,beta_true=data$truth[-1,,drop=FALSE])
  comparison <- sparseglm(sim,method=c("L0glm","lasso"),tune="cv",nfolds=2,
    nlambda=3,lambda=list(L0glm=c(10,4,2)))
  details <- summary(comparison)
  expect_s3_class(details,"summary.sparseglm")
  expect_identical(details$metrics,L0glm.metrics(comparison))
  expect_identical(details$call,comparison$call)
  expect_equal(details$nfolds,2)
  expect_equal(details$tune,"cv")
  expect_output(print(details),"cv_mean")
  expect_output(print(details),"F1_score")
  expect_output(print(comparison),"CV \\(2 folds\\)")
  capture.output(visible <- withVisible(print(details)))
  expect_false(visible$visible)
  expect_identical(visible$value,details)
  expect_true(all(is.na(summary(comparison,beta_true=NULL)$metrics$MCC)))
  expect_warning(failed <- sparseglm(sim,method=c("L0glm","lasso"),nlambda=3,
    options=list(L0glm=list(solver="invalid"))),"failed")
  expect_equal(summary(failed)$metrics$status,c("error","ok"))
  expect_output(print(summary(failed)),"Method notes")
  failed$fits <- failed$fits[1L]; failed$timings <- failed$timings[1L,,drop=FALSE]
  expect_output(print(summary(failed)),"error")
})

test_that("summary and print methods never refit models", {
  data <- model_summary_fixture()
  comparison <- sparseglm(data$X,data$y,method="lasso",nlambda=3)
  local_mocked_bindings(L0glm.fit=function(...) stop("Must not refit"),
    L0nnpois=function(...) stop("Must not refit"),L0glm.tune=function(...) stop("Must not tune"))
  expect_s3_class(summary(data$fit),"summary.L0glm_fit")
  expect_s3_class(summary(comparison),"summary.sparseglm")
  local_mocked_bindings(L0glm.metrics=function(...) stop("Must not score"))
  expect_output(print(comparison),"Sparse GLM comparison")
  expect_output(print(comparison$fits$lasso),"L0glm fit")
})
