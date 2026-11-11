# Run from the checkout root or its L0glm package directory.
# Set run_full=TRUE for the published 1.73-million-predictor throughput case.
library(L0glm)
run_full <- FALSE
root <- if(file.exists('L0glm/DESCRIPTION'))getwd()else dirname(getwd())
cfg <- if(run_full)list(n=4000L,p=1733480L,train=3200L,nsis=1000L)else
  list(n=400L,p=2000L,train=320L,nsis=100L)
bed_file <- file.path(tempdir(),paste0('genotypes-',cfg$n,'-',cfg$p,'.bed'))
out <- file.path(root,'benchmarks/plots/file-backed')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
if(!file.exists(bed_file)) {
  Rcpp::sourceCpp(file.path(root,'benchmarks/helpers/packed_genotypes.cpp'))
  benchmarkBedData(bed_file,cfg$n,cfg$p,101L)
}
X <- l0glmData(bed_file,n=cfg$n,p=cfg$p,
  control=list(block.columns=512L,block.cells=4e6,max.cells=5e6))
ids <- if(run_full)as.integer(c(1,2300,20000,1000000,cfg$p))else c(1L,3L,23L,1154L,2000L)
beta <- c(.9,-.9,.8,-.8,.9)
set.seed(101)
signal <- drop(as.matrix(X[,ids,drop=FALSE])%*%beta)
offset <- rnorm(cfg$n,0,.1)
y <- rbinom(cfg$n,1,plogis(-.5+signal+offset))
weights <- sample(1:3,cfg$n,TRUE)
training <- seq_len(cfg$train);test <- setdiff(seq_len(cfg$n),training)
truth <- Matrix::sparseMatrix(i=ids+1L,j=rep(1L,length(ids)),x=beta,dims=c(cfg$p+1L,1L))

run <- function() {
  old <- RhpcBLASctl::blas_get_num_procs()
  on.exit(RhpcBLASctl::blas_set_num_threads(old))
  RhpcBLASctl::blas_set_num_threads(1L)
  results <- list();reference <- NULL
  for(threads in c(1L,2L,4L,8L)) {
    elapsed <- system.time(fit <- l0glmFit(X[training,,drop=FALSE],y[training],
      family=binomial(),weights=weights[training],offset=offset[training],
      intercept=TRUE,screen='SIS',nsis=cfg$nsis,screen.nthreads=threads,
      lambda='gic',local.search=FALSE,nthreads=1L))[[3L]]
    if(is.null(reference))reference <- coef(fit)else stopifnot(identical(reference,coef(fit)))
    score <- l0glmMetrics(fit,beta_true=truth)[c('gic','TP','FP','FN','F1_score','MCC')]
    results[[length(results)+1L]] <- cbind(method='SIS',threads,seconds=elapsed,score,
      logloss=-mean(dbinom(y[test],1,predict(fit,X[test,,drop=FALSE],offset=offset[test]),log=TRUE)))
  }
  elapsed <- system.time(fit <- l0glmFit(X[training,,drop=FALSE],y[training],
    family=binomial(),weights=weights[training],offset=offset[training],
    intercept=TRUE,screen='ISIS',screen.maxit=3L,screen.nthreads=8L,
    lambda='gic',local.search='auto',nthreads=1L))[[3L]]
  score <- l0glmMetrics(fit,beta_true=truth)[c('gic','TP','FP','FN','F1_score','MCC')]
  results[[length(results)+1L]] <- cbind(method='ISIS',threads=8L,seconds=elapsed,score,
    logloss=-mean(dbinom(y[test],1,predict(fit,X[test,,drop=FALSE],offset=offset[test]),log=TRUE)))
  do.call(rbind,results)
}
results <- run()
print(results,row.names=FALSE)
write.csv(results,file.path(out,if(run_full)'results.csv'else 'small-example.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'session.txt'))
