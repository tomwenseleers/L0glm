library(L0glm)
library(grpreg)
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)
root <- if(file.exists('L0glm/DESCRIPTION')) '.' else '..'
out <- file.path(root,'benchmarks/plots/predictor-groups')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
replications <- 3L
nlambda <- 40L
cases <- list(
  gaussian_equal=list(family='gaussian',n=400L,size=rep(4L,50L),active=6L,rho=.3,strength=1),
  gaussian_wide=list(family='gaussian',n=400L,size=rep(c(2L,4L,8L,10L),25L),active=8L,rho=.5,strength=1),
  logistic_equal=list(family='binomial',n=600L,size=rep(4L,60L),active=8L,rho=.3,strength=1.8),
  logistic_unequal=list(family='binomial',n=1000L,size=rep(c(2L,4L,8L,10L),25L),active=10L,rho=.6,strength=2))
design <- function(n,p,rho) {
  X <- matrix(rnorm(n*p),n,p)
  for(j in seq.int(2L,p))X[,j] <- rho*X[,j-1L]+sqrt(1-rho^2)*X[,j]
  X
}
auc <- function(y,p) {
  a <- sum(y==1);b <- sum(y==0)
  (sum(rank(p)[y==1])-a*(a+1)/2)/(a*b)
}
results <- list();counter <- 0L
for(case in names(cases))for(replication in seq_len(replications)) {
  cfg <- cases[[case]]
  set.seed(700L+replication+match(case,names(cases))*100L)
  group <- c(0L,0L,rep(seq_along(cfg$size),cfg$size))
  p <- length(group);X <- design(cfg$n,p,cfg$rho);V <- design(2000L,p,cfg$rho)
  active <- sort(sample.int(length(cfg$size),cfg$active));beta <- numeric(p);beta[1:2] <- c(.3,-.3)
  for(g in active) {
    j <- which(group==g)
    beta[j] <- sample(c(-1,1),length(j),TRUE)*runif(length(j),.8,1.2)*cfg$strength/sqrt(length(j))
  }
  eta <- .3+as.numeric(X%*%beta);eta_test <- .3+as.numeric(V%*%beta)
  y <- if(cfg$family=='gaussian')eta+rnorm(cfg$n)else rbinom(cfg$n,1,plogis(eta))
  v <- if(cfg$family=='gaussian')eta_test+rnorm(2000L)else rbinom(2000L,1,plogis(eta_test))
  fold <- integer(cfg$n)
  for(k in if(cfg$family=='gaussian')list(seq_len(cfg$n))else split(seq_len(cfg$n),y))
    fold[k] <- sample(rep(1:5,length.out=length(k)))
  penalty <- presetLambda('gic',p=sum(group!=0),n=cfg$n)
  group_df <- function(b) {
    chosen <- unique(group[which(b[-1L]!=0)]);chosen <- chosen[chosen!=0]
    c(groups=length(chosen),df=3+sum(group!=0 & b[-1L]!=0))
  }
  common <- function(b) {
    predicted <- as.numeric(b[1L]+X%*%b[-1L])
    ll <- if(cfg$family=='gaussian')sum(dnorm(y,predicted,1,log=TRUE))else
      sum(dbinom(y,1,plogis(predicted),log=TRUE))
    -2*ll+penalty*group_df(b)[['df']]
  }
  record <- function(method,tuning,b,seconds,converged=NA,edge=FALSE,note='') {
    selected <- setdiff(unique(group[which(b[-1L]!=0)]),0L)
    tp <- sum(selected%in%active);fp <- sum(!selected%in%active)
    fn <- length(active)-tp;tn <- length(cfg$size)-tp-fp-fn
    den <- sqrt(as.double(tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    estimate <- as.numeric(b[1L]+V%*%b[-1L])
    counter <<- counter+1L
    results[[counter]] <<- data.frame(case,replication,family=cfg$family,n=cfg$n,p,
      groups=length(cfg$size),true_groups=length(active),method,tune=tuning,seconds,GIC=common(b),
      TP=tp,FP=fp,FN=fn,TN=tn,F1=2*tp/(2*tp+fp+fn),
      MCC=if(den>0)(as.double(tp)*tn-fp*fn)/den else NA_real_,
      test_RMSE=if(cfg$family=='gaussian')sqrt(mean((v-estimate)^2))else NA_real_,
      test_AUC=if(cfg$family=='binomial')auc(v,estimate)else NA_real_,
      test_logloss=if(cfg$family=='binomial')-mean(dbinom(v,1,plogis(estimate),log=TRUE))else NA_real_,
      converged,grid_edge=edge,note)
    cat(case,replication,method,tuning,seconds,'TP',tp,'FP',fp,'\n')
    write.csv(do.call(rbind,results),file.path(out,'runs.csv'),row.names=FALSE)
  }
  f <- if(cfg$family=='gaussian')gaussian()else binomial()
  args <- list(X=X,y=y,family=f,group=group,lambda='gic',intercept=TRUE,
    dispersion=1,local.search='auto',parallelize=FALSE,nthreads=1L)
  seconds <- system.time(fit <- do.call(l0glmFit,args),gcFirst=FALSE)[['elapsed']]
  record('L0glm','GIC preset',as.numeric(coef(fit)),seconds,all(fit$convergence$converged),
    note='Known Gaussian variance 1; binomial dispersion 1')
  if(cfg$family=='gaussian') {
    args$dispersion <- NULL
    seconds <- system.time(estimated <- do.call(l0glmFit,args),gcFirst=FALSE)[['elapsed']]
    record('L0glm estimated scale','GIC preset',as.numeric(coef(estimated)),seconds,
      all(estimated$convergence$converged),note='Estimated fitting scale; common scores use true variance 1')
    args$dispersion <- 1
  }
  args$lambda <- NULL
  seconds <- system.time(cv <- do.call(l0glmTune,c(args,list(tune='cv',nfolds=5L,foldid=fold,
    nlambda=10L,plot=FALSE,printbest=FALSE))),gcFirst=FALSE)[['elapsed']]
  record('L0glm','5-fold CV',as.numeric(coef(cv$fit_best)),seconds,
    all(cv$fit_best$convergence$converged),cv$lambda_best%in%range(cv$lambdas))
  for(penalty_name in c('grLasso','grMCP','grSCAD')) {
    messages <- character()
    seconds <- system.time(path <- withCallingHandlers(grpreg(X,y,group=group,
      family=cfg$family,penalty=penalty_name,nlambda=nlambda),warning=function(w) {
        messages <<- c(messages,conditionMessage(w));invokeRestart('muffleWarning')
      }),gcFirst=FALSE)[['elapsed']]
    selection_seconds <- system.time({
      coefficient <- as.matrix(path$beta)
      stopifnot(nrow(coefficient)==p+1L)
      score <- apply(coefficient,2L,common);best <- which.min(score)
    },gcFirst=FALSE)[['elapsed']]
    seconds <- seconds+selection_seconds
    record(penalty_name,'GIC path',coefficient[,best],seconds,edge=best%in%c(1,ncol(coefficient)),
      note=paste(unique(messages),collapse='; '))
    messages <- character()
    seconds <- system.time(cv <- withCallingHandlers(cv.grpreg(X,y,group=group,
      family=cfg$family,penalty=penalty_name,nlambda=nlambda,nfolds=5L,fold=fold),warning=function(w) {
        messages <<- c(messages,conditionMessage(w));invokeRestart('muffleWarning')
      }),gcFirst=FALSE)[['elapsed']]
    record(penalty_name,'5-fold CV',as.numeric(cv$fit$beta[,cv$min]),seconds,
      edge=cv$min%in%c(1,ncol(cv$fit$beta)),note=paste(unique(messages),collapse='; '))
  }
}
writeLines(capture.output(sessionInfo()),file.path(out,'session.txt'))
runs <- do.call(rbind,results)
key <- interaction(runs$case,runs$method,runs$tune,drop=TRUE)
summary <- do.call(rbind,lapply(split(runs,key),function(d) {
  data.frame(d[1,c('case','family','n','p','groups','true_groups','method','tune')],
    replications=nrow(d),seconds=median(d$seconds),GIC=mean(d$GIC),
    TP=mean(d$TP),FP=mean(d$FP),F1=mean(d$F1),MCC=mean(d$MCC),
    test_RMSE=mean(d$test_RMSE),test_AUC=mean(d$test_AUC),
    test_logloss=mean(d$test_logloss),converged=mean(d$converged))
}))
rownames(summary) <- NULL
write.csv(summary,file.path(out,'summary.csv'),row.names=FALSE)
print(summary,row.names=FALSE)
