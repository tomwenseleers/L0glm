# Reproduce the BIC and five-fold-CV LASSO comparisons in the README and website.
# Run from the repository root. The 1,000 x 1,000,000 design occupies 8 GB.
# A fresh R worker per workflow bounds memory retained by native allocations.
million_lasso_worker <- function(input,method,tune,n,p,n_nonzero,nlambda,folds,libs) {
  .libPaths(libs)
  library(L0glm)
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
  data <- readRDS(input)
  stopifnot(identical(dim(data$X),as.integer(c(n,p))),sum(data$beta!=0)==n_nonzero)
  options <- list(type.gaussian='naive')
  if('control'%in%names(formals(glmnet::glmnet)))options$control <- list(pmax=2*n+20L)else
    options$pmax <- 2*n+20L
  gc()
  result <- sparseGLM(data$X,data$y,family=gaussian(),method=method,
    tune=tune,nlambda=nlambda,nfolds=length(unique(folds)),
    foldid=if(tune=='cv')folds else NULL,intercept=FALSE,beta_true=data$beta,
    screen='none',options=setNames(list(options),method))
  scores <- l0glmMetrics(result)
  if(any(scores$status!='ok'))stop(scores$note,call.=FALSE)
  scores$method <- paste0(method,'_',tune)
  scores$tune <- if(tune=='cv')paste0(length(unique(folds)),'-fold CV')else 'BIC path'
  list(scores=scores,session=capture.output(sessionInfo()))
}

benchmark_million_lasso <- function(n=1000L,p=1000000L,n_nonzero=60L,
  nlambda=40L,nfolds=5L,output='benchmarks/plots/website',input=NULL,
  cache='codex_sandbox/website_large_examples') {
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  dir.create(cache,recursive=TRUE,showWarnings=FALSE)
  if(is.null(input)) {
    set.seed(123,kind='Mersenne-Twister',normal.kind='Inversion',sample.kind='Rejection')
    X <- matrix(rnorm(as.double(n)*p),n,p)
    active <- sample.int(p,n_nonzero)
    beta <- numeric(p);beta[active] <- runif(n_nonzero,-.8,.8)
    set.seed(1)
    y <- drop(X[,active,drop=FALSE]%*%beta[active])+rnorm(n)
    input <- file.path(cache,paste0('million-lasso-',n,'-',p,'.rds'))
    saveRDS(list(X=X,y=y,beta=beta),input,compress=FALSE)
    rm(X,y,beta);gc()
  }
  input <- normalizePath(input,winslash='/',mustWork=TRUE)
  set.seed(42)
  folds <- sample(rep(seq_len(nfolds),length.out=n))
  records <- list()
  for(tune in c('bic','cv'))for(method in c('lasso','relaxed_lasso')) {
    message('Starting ',method,' / ',tune,' at ',Sys.time())
    cluster <- parallel::makePSOCKcluster(1L,outfile=file.path(cache,'million-lasso-worker.log'))
    result <- tryCatch(parallel::clusterCall(cluster,million_lasso_worker,input,
      method,tune,n,p,n_nonzero,nlambda,folds,.libPaths())[[1L]],
      finally=parallel::stopCluster(cluster))
    records[[paste(method,tune)]] <- result$scores
    scores <- do.call(rbind,records)
    write.csv(scores,file.path(output,'million-lasso-metrics.csv'),row.names=FALSE)
    writeLines(sub('[[:blank:]]+$','',result$session),file.path(output,'million-lasso-sessionInfo.txt'))
    print(result$scores[c('method','seconds','gic','TP','FP','F1_score','MCC')],row.names=FALSE)
    gc()
  }
  invisible(scores)
}

if(sys.nframe()==0L)benchmark_million_lasso()
