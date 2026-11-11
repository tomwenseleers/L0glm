# A nested-design extension of the Gaussian/logistic support benchmarks.
# Outcomes and signal columns remain fixed when null predictors are added.
benchmark_support_difficulty <- function(output="benchmarks/plots/website",
    repetitions=3L, predictors=c(1000L,20000L),
    methods=c("abess","L0Learn","lasso","relaxed_lasso","adaptive_lasso",
              "adaptive_relaxed_lasso","MCP","SCAD")) {
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  rng <- RNGkind();on.exit(do.call(RNGkind,as.list(rng)),add=TRUE)
  RNGkind("Mersenne-Twister","Inversion","Rejection")
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
  ar_design <- function(n,p,rho=.3) {
    X <- matrix(rnorm(n*p),n,p)
    if(p>1L) for(j in 2:p) X[,j] <- rho*X[,j-1L]+sqrt(1-rho^2)*X[,j]
    X
  }
  auc <- function(y,mu) {
    yes <- y==1
    (sum(rank(mu)[yes])-sum(yes)*(sum(yes)+1)/2)/(sum(yes)*sum(!yes))
  }
  records <- list()
  for(model in c("gaussian","binomial")) for(seed in seq_len(repetitions)) {
    n <- if(model=="gaussian")500L else 1000L
    k <- if(model=="gaussian")40L else 20L
    family <- if(model=="gaussian")gaussian()else binomial()
    set.seed(seed+if(model=="gaussian")1000L else 2000L)
    active <- sample.int(min(predictors),k)
    beta <- numeric(max(predictors));beta[active] <- sample(c(-1,1),k,TRUE)*runif(k,.35,.8)
    train <- ar_design(n,max(predictors));test <- ar_design(n,max(predictors))
    eta <- as.numeric(train[,active]%*%beta[active]);eta.test <- as.numeric(test[,active]%*%beta[active])
    y <- if(model=="gaussian")eta+rnorm(n)else rbinom(n,1,plogis(eta))
    y.test <- if(model=="gaussian")eta.test+rnorm(n)else rbinom(n,1,plogis(eta.test))
    foldid <- sample(rep(1:5,length.out=n))
    outcome_hash <- digest::digest(list(y,y.test,active,beta[active]),algo="xxhash64")
    for(p in predictors) {
      X <- train[,seq_len(p),drop=FALSE];V <- test[,seq_len(p),drop=FALSE]
      data_hash <- digest::digest(list(X,y,V,y.test,beta[seq_len(p)]),algo="xxhash64")
      collect <- function(fit,method,tune,seconds,status="ok",note="") {
        if(status!="ok") return(data.frame(model,n,p,k,seed,data_hash,outcome_hash,method,tune,
          seconds,status,note,GIC=NA_real_,TP=NA,FP=NA,FN=NA,TN=NA,F1=NA,MCC=NA,prediction=NA))
        b <- as.numeric(coef(fit))[-1L];found <- abs(b)>1e-5;truth <- beta[seq_len(p)]!=0
        stopifnot(length(b)==p,all(is.finite(b)))
        tp <- sum(found&truth);fp <- sum(found&!truth);fn <- sum(!found&truth);tn <- sum(!found&!truth)
        denom <- sqrt(as.double(tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
        mu <- as.numeric(predict(fit,newdata=V,type="response"))
        stopifnot(length(mu)==n,all(is.finite(mu)))
        data.frame(model,n,p,k,seed,data_hash,outcome_hash,method,tune,seconds,status,note,
          GIC=as.numeric(L0glm::GIC(fit)),TP=tp,FP=fp,FN=fn,TN=tn,F1=2*tp/(2*tp+fp+fn),
          MCC=if(denom)(as.double(tp)*tn-fp*fn)/denom else 0,
          prediction=if(model=="gaussian")sqrt(mean((mu-y.test)^2))else auc(y.test,mu))
      }
      elapsed <- system.time(fit <- L0glm::l0glmFit(X,y,family=family,intercept=TRUE,
        lambda="gic",screen.nthreads=1L,parallelize=FALSE))[["elapsed"]]
      records[[length(records)+1L]] <- collect(fit,"L0glm","gic",elapsed)
      for(tuning in c("gic","cv")) {
        selected <- if(tuning=="gic")methods else c("lasso","relaxed_lasso")
        fits <- L0glm::sparseGLM(X,y,family=family,method=selected,tune=tuning,
          nlambda=40L,support.size=0:100,intercept=TRUE,screen="none",
          nfolds=5L,foldid=if(tuning=="cv")foldid else NULL,options=list(abess=list(thread=1L)))
        metrics <- L0glm::l0glmMetrics(fits)
        for(i in seq_along(fits$fits)) {
          records[[length(records)+1L]] <- collect(fits$fits[[i]],metrics$method[i],tuning,
            metrics$seconds[i],metrics$status[i],metrics$note[i])
        }
      }
      write.csv(do.call(rbind,records),file.path(output,"support-difficulty-runs.csv"),row.names=FALSE)
      message(model,", p=",p,", seed=",seed," completed")
      rm(fit,fits);gc()
    }
  }
  runs <- do.call(rbind,records)
  ok <- runs[runs$status=="ok",]
  summary <- merge(aggregate(cbind(GIC,TP,FP,F1,MCC,prediction)~model+n+p+k+method+tune,ok,mean),
    aggregate(seconds~model+n+p+k+method+tune,ok,median),by=c("model","n","p","k","method","tune"))
  write.csv(summary,file.path(output,"support-difficulty.csv"),row.names=FALSE)
  writeLines(trimws(capture.output(sessionInfo()),which="right"),
    file.path(output,"support-difficulty-session.txt"))
  invisible(summary)
}
