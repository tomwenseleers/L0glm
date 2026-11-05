# Compare equal endpoints: screening ranks first, complete fitted models next.
library(L0glm)
stopifnot(requireNamespace("SIS",quietly=TRUE),requireNamespace("glmnet",quietly=TRUE))
RhpcBLASctl::blas_set_num_threads(1L)
n <- 500L
p <- 1000000L
nsis <- 100L
repetitions <- 3L
root <- if (file.exists("L0glm/DESCRIPTION")) "." else ".."
output <- file.path(root,"benchmarks","3_benchmarks_sis_vs_SIS_results")
dir.create(output,recursive=TRUE,showWarnings=FALSE)

make_data <- function(p) {
  set.seed(123,kind="Mersenne-Twister")
  X <- matrix(rnorm(as.double(n)*p),n,p)
  active <- sample.int(p,60L)
  beta <- numeric(p)
  beta[active] <- runif(60L,-.8,.8)
  set.seed(1)
  list(X=X,y=as.numeric(X[,active,drop=FALSE]%*%beta[active])+rnorm(n),
    beta_true=beta,intercept=FALSE)
}
sim <- make_data(p)

# SIS 1.5 exposes screening through its internal obtain.ix0 helper. Calling
# just this helper excludes its subsequent LASSO fit. Complete ranks are kept
# by all three calls. Gaussian ranks are invariant to column standardization.
sis_screen <- getFromNamespace("obtain.ix0","SIS")
screening <- list(); reference_rank <- NULL
for (repetition in seq_len(repetitions)) for (method in c("sis_glm_1","sis_glm_auto","SIS_screen")) {
  gc()
  seconds <- system.time(fit <- switch(method,
    sis_glm_1=sis_glm(sim$X,sim$y,nsis=nsis,nthreads=1L),
    sis_glm_auto=sis_glm(sim$X,sim$y,nsis=nsis),
    SIS_screen=sis_screen(sim$X,sim$y,s1=NULL,s2=NULL,family="gaussian",
      nsis=nsis,iter=FALSE,varISIS="vanilla",perm=FALSE,q=1,
      greedy=FALSE,greedy.size=1L,iterind=0L)))[["elapsed"]]
  ranking <- as.integer(if (method=="SIS_screen") fit$ix0all else fit$ranking)
  if (is.null(reference_rank)) reference_rank <- ranking
  row <- data.frame(method,repetition,seconds,
    threads=if (method=="SIS_screen") NA_integer_ else fit$nthreads,
    identical_top_100=identical(ranking[seq_len(nsis)],reference_rank[seq_len(nsis)]),
    identical_full_rank=identical(ranking,reference_rank),
    rank_correlation=cor(order(ranking),order(reference_rank)))
  stopifnot(row$identical_full_rank)
  screening[[length(screening)+1L]] <- row
  print(row)
}
screening <- do.call(rbind,screening)
write.csv(screening,file.path(output,"timings.csv"),row.names=FALSE)

wrap_sis <- function(fit,sim) {
  B <- Matrix::sparseMatrix(i=c(1L,fit$ix+1L),j=rep(1L,length(fit$ix)+1L),x=as.numeric(fit$coef.est),
    dims=c(ncol(sim$X)+1L,1L))
  structure(list(coefficients=B,X=cbind(Intercept=1,sim$X),y=matrix(sim$y),
    family=gaussian(),prior.weights=matrix(1,n,1),offset=matrix(0,n,1),
    dispersion=1,dispersion.estimated=TRUE,p=ncol(sim$X),always.include=1L),class="L0glm_fit")
}

# Every row includes screening AND fitting. L0glm uses its default GIC preset;
# LASSO variants use BIC because GIC often selects a null LASSO model here.
# The million-predictor ISIS L0glm call is exactly the homepage default fit.
run_models <- function(sim,screen,repetitions,
  methods=c("SIS_package","L0glm.fit","sparseglm_lasso","sparseglm_relaxed_lasso"),
  batch=nsis,rounds=3L) {
  p <- ncol(sim$X); result <- list()
  for (repetition in seq_len(repetitions)) {
    for (method in methods) {
      gc()
      seconds <- system.time(fit <- switch(method,
        SIS_package=SIS::SIS(sim$X,sim$y,family="gaussian",nsis=nsis,
          iter=screen=="ISIS",iter.max=rounds,penalty="lasso",tune="bic",
          standardize=FALSE,parallel=FALSE,seed=123L),
        L0glm.fit=if (p==1000000L && screen=="ISIS" && is.null(batch) && rounds==10L)
          L0glm.fit(sim$X,sim$y,family=gaussian()) else
          L0glm.fit(sim$X,sim$y,family=gaussian(),screen=screen,nsis=batch,screen.maxit=rounds),
        sparseglm_lasso=sparseglm(sim,method="lasso",screen=screen,nsis=batch,
          screen.maxit=rounds,tune="bic",nlambda=10L),
        sparseglm_relaxed_lasso=sparseglm(sim,method="relaxed_lasso",screen=screen,
          nsis=batch,screen.maxit=rounds,tune="bic",nlambda=10L)))[["elapsed"]]
      model <- if (method=="SIS_package") wrap_sis(fit,sim) else
        if (inherits(fit,"sparseglm")) fit$fits[[1L]] else fit
      truth <- if (nrow(coef(model))==p+1L) c(0,sim$beta_true) else sim$beta_true
      scores <- L0glm.metrics(model,beta_true=truth)
      row <- cbind(data.frame(n,p,screen,method,repetition,seconds),
        scores[c("gic","bic","TP","FP","FN","F1_score","MCC")])
      result[[length(result)+1L]] <- row
      print(row)
      rm(fit,model)
    }
  }
  do.call(rbind,result)
}
models <- list(run_models(sim,"SIS",repetitions))
models[[2L]] <- run_models(sim,"ISIS",repetitions,
  methods=c("L0glm.fit","sparseglm_lasso","sparseglm_relaxed_lasso"),batch=NULL,rounds=10L)
rm(sim);gc()
# Conditional per-predictor fits in SIS::SIS(iter=TRUE) are much more costly
# than its first correlation screen, so test both complete ISIS pipelines on
# a separate 5,000-predictor example.
sim <- make_data(5000L)
models[[3L]] <- run_models(sim,"ISIS",repetitions)
models <- do.call(rbind,models)
write.csv(models,file.path(output,"models.csv"),row.names=FALSE)
print(aggregate(seconds~method,screening,median))
print(aggregate(seconds~p+screen+method,models,median))
writeLines(trimws(capture.output(sessionInfo()),which="right"),file.path(output,"sessionInfo.txt"))
