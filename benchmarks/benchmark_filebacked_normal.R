# Reproduce a Gaussian ISIS fit without loading the complete design into RAM.
# Run from the checkout root. The full double-precision file occupies 80 GB.
benchmark_filebacked_normal <- function(n=10000L,p=1000000L,n_nonzero=60L,
  train=floor(.8*n),block.columns=256L,fit.block.columns=4096L,intercept=TRUE,
  screen.nthreads=8L,
  cache='codex_sandbox/filebacked_normal',output='benchmarks/plots/file-backed') {
  stopifnot(requireNamespace('bigmemory',quietly=TRUE),n>train,train>2L,
    p>=n_nonzero,n_nonzero>0L,block.columns>0L,fit.block.columns>0L)
  library(L0glm)
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
  dir.create(cache,recursive=TRUE,showWarnings=FALSE)
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  cache <- normalizePath(cache,winslash='/',mustWork=TRUE)
  config <- list(n=as.integer(n),p=as.integer(p),n_nonzero=as.integer(n_nonzero),seed=123L)
  stem <- paste0('normal-',n,'-',p)
  descriptor <- file.path(cache,paste0(stem,'.desc'))
  metadata <- file.path(cache,paste0(stem,'.rds'))
  if(file.exists(descriptor)||file.exists(metadata)) {
    stopifnot(file.exists(descriptor),file.exists(metadata))
    saved <- readRDS(metadata)
    stopifnot(identical(saved$config,config))
    backing <- bigmemory::attach.big.matrix(descriptor)
  } else {
    generation <- system.time({
      backing <- bigmemory::filebacked.big.matrix(nrow=n,ncol=p,type='double',
        backingfile=paste0(stem,'.bin'),descriptorfile=paste0(stem,'.desc'),
        backingpath=cache)
      set.seed(123,kind='Mersenne-Twister',normal.kind='Inversion',sample.kind='Rejection')
      for(first in seq.int(1L,p,by=block.columns)) {
        columns <- seq.int(first,min(p,first+block.columns-1L))
        backing[,columns] <- matrix(rnorm(as.double(n)*length(columns)),n)
        if(first%%(block.columns*100L)==1L)message('Generated columns ',first,' / ',p)
      }
      active <- sample.int(p,n_nonzero)
      beta <- runif(n_nonzero,-.8,.8)
      set.seed(1)
      y <- drop(backing[,active,drop=FALSE]%*%beta)+rnorm(n)
    })[['elapsed']]
    saved <- list(config=config,active=active,beta=beta,y=y,generation_seconds=generation)
    saveRDS(saved,metadata)
  }
  X <- l0glmData(backing,control=list(block.columns=fit.block.columns,
    block.cells=as.double(train)*fit.block.columns,max.cells=1e8))
  training <- seq_len(train);test <- seq.int(train+1L,n)
  truth <- Matrix::sparseMatrix(i=saved$active+as.integer(intercept),j=rep(1L,n_nonzero),
    x=saved$beta,dims=c(p+as.integer(intercept),1L))
  gc()
  seconds <- system.time(fit <- l0glmFit(X[training,,drop=FALSE],saved$y[training],
    family=gaussian(),intercept=intercept,screen='ISIS',screen.nthreads=screen.nthreads))[['elapsed']]
  metrics <- l0glmMetrics(fit,beta_true=truth)
  metrics$seconds <- seconds
  metrics$test_RMSE <- sqrt(mean((predict(fit,X[test,,drop=FALSE])-saved$y[test])^2))
  metrics$observations <- n;metrics$training <- train;metrics$test <- length(test)
  metrics$file_GB <- file.info(file.path(cache,paste0(stem,'.bin')))$size/1e9
  metrics$generation_seconds <- saved$generation_seconds
  metrics$screening_threads <- fit$screening$nthreads
  metrics$block_columns <- fit.block.columns;metrics$intercept <- intercept
  stopifnot(identical(dim(coef(fit)),c(as.integer(p+intercept),1L)),inherits(coef(fit),'sparseMatrix'))
  write.csv(metrics,file.path(output,paste0(stem,'-metrics.csv')),row.names=FALSE)
  write.csv(fit$screening$history,file.path(output,paste0(stem,'-history.csv')),row.names=FALSE)
  saveRDS(list(coefficients=coef(fit),screening=fit$screening,metrics=metrics),
    file.path(cache,paste0(stem,'-fit.rds')))
  print(metrics[c('observations','training','p','seconds','file_GB','gic','TP','FP','FN',
    'F1_score','MCC','test_RMSE')],row.names=FALSE)
  invisible(list(metrics=metrics,fit=fit,source=X,data=saved))
}

if(sys.nframe()==0L)benchmark_filebacked_normal()
