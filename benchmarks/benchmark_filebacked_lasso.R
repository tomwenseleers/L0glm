# Native file-backed LASSO comparison on the reusable 80 GB Gaussian design.
# Run benchmark_filebacked_normal.R first to create the descriptor and response.
benchmark_filebacked_lasso <- function(n=10000L,p=1000000L,train=floor(.8*n),
  nlambda=10L,ncores=8L,cache='codex_sandbox/filebacked_normal',
  output='benchmarks/plots/file-backed') {
  stopifnot(requireNamespace('bigmemory',quietly=TRUE),requireNamespace('biglasso',quietly=TRUE))
  library(L0glm)
  RhpcBLASctl::blas_set_num_threads(1L);RhpcBLASctl::omp_set_num_threads(1L)
  prefix <- file.path(cache,paste0('normal-',n,'-',p))
  saved <- readRDS(paste0(prefix,'.rds'))
  stopifnot(saved$config$n==n,saved$config$p==p,train>2L,train<n)
  M <- bigmemory::attach.big.matrix(paste0(prefix,'.desc'))
  X <- l0glmData(M,control=list(block.columns=4096L,block.cells=as.double(train)*4096L,
    max.cells=1e8))
  training <- seq_len(train);test <- seq.int(train+1L,n)
  truth <- Matrix::sparseMatrix(i=saved$active+1L,j=rep(1L,length(saved$active)),
    x=saved$beta,dims=c(p+1L,1L))
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  path <- file.path(output,paste0('normal-',n,'-',p,'-lasso.csv'))
  results <- list()
  for(method in c('L0glm','biglasso')) {
    message('START ',method,' ',Sys.time())
    gc()
    seconds <- system.time({
      if(method=='L0glm')fit <- l0glmFit(X[training,,drop=FALSE],saved$y[training],
        family=gaussian(),intercept=TRUE,screen='ISIS',screen.nthreads=ncores) else {
        comparison <- sparseGLM(X[training,,drop=FALSE],saved$y[training],
          family=gaussian(),intercept=TRUE,method='biglasso',screen='none',tune='bic',
          nlambda=nlambda,options=list(ncores=ncores,lambda.min=.05,verbose=TRUE))
        stopifnot(comparison$timings$status=='ok')
        fit <- comparison$fits[[1L]]
      }
    })[['elapsed']]
    metrics <- l0glmMetrics(fit,beta_true=truth)
    metrics$method <- method;metrics$seconds <- seconds
    metrics$tune <- if(method=='L0glm')'GIC preset'else'BIC path'
    metrics$screen <- if(method=='L0glm')'ISIS'else'none'
    metrics$nlambda <- if(method=='L0glm')1L else nrow(comparison$tuning[[1]]$candidates)
    metrics$grid_interior <- if(method=='L0glm')NA else
      comparison$tuning[[1]]$index>1L && comparison$tuning[[1]]$index<metrics$nlambda
    metrics$test_RMSE <- sqrt(mean((predict(fit,X[test,,drop=FALSE])-saved$y[test])^2))
    metrics$observations <- n;metrics$training <- train;metrics$test <- length(test)
    metrics$file_GB <- file.info(paste0(prefix,'.bin'))$size/1e9
    metrics$ncores <- ncores;metrics$intercept <- TRUE
    results[[method]] <- metrics
    write.csv(do.call(rbind,results),path,row.names=FALSE)
    checkpoint <- list(coefficients=coef(fit),metrics=metrics,
      screening=fit$screening,lambda=if(method=='L0glm')fit$init.lambda else comparison$tuning[[1]]$best,
      tuning=if(method=='biglasso')comparison$tuning[[1]]else NULL)
    saveRDS(checkpoint,paste0(prefix,'-',method,'-comparison.rds'))
    print(metrics[c('method','seconds','gic','TP','FP','FN','F1_score','MCC','test_RMSE')])
    rm(fit);if(exists('comparison',inherits=FALSE))rm(comparison)
  }
  writeLines(capture.output(sessionInfo()),file.path(output,'filebacked-lasso-sessionInfo.txt'))
  invisible(do.call(rbind,results))
}

if(sys.nframe()==0L)benchmark_filebacked_lasso()
