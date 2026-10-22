# Parallel cv.glmnet for standalone benchmark scripts (including Windows).
# Reuse an existing parallel backend; otherwise manage a temporary PSOCK cluster.
cv_glmnet_parallel <- function(x, y, ..., parallel=TRUE, ncores=NULL) {
  dots <- list(...)
  if (!parallel || foreach::getDoParWorkers() > 1L)
    return(do.call(glmnet::cv.glmnet,c(list(x=x,y=y,parallel=parallel),dots)))
  if (!requireNamespace("doParallel",quietly=TRUE))
    stop("Install doParallel to run parallel glmnet cross-validation.",call.=FALSE)
  folds <- if (!is.null(dots$foldid)) length(unique(dots$foldid)) else
    if (!is.null(dots$nfolds)) dots$nfolds else 10L
  if (is.null(ncores)) {
    ncores <- suppressWarnings(parallel::detectCores(logical=FALSE))
    if (is.na(ncores)) ncores <- 1L
    # PSOCK workers copy dense designs: bound concurrent copies for large inputs.
    if (as.numeric(object.size(x)) > 512*1024^2) ncores <- min(ncores,2L)
  }
  ncores <- max(1L,min(as.integer(ncores),folds))
  if (ncores == 1L) return(do.call(glmnet::cv.glmnet,c(list(x=x,y=y,parallel=FALSE),dots)))
  cluster <- parallel::makePSOCKcluster(ncores)
  on.exit({parallel::stopCluster(cluster);foreach::registerDoSEQ()},add=TRUE)
  parallel::clusterCall(cluster,function() {
    if (requireNamespace("RhpcBLASctl",quietly=TRUE)) {
      RhpcBLASctl::blas_set_num_threads(1L)
      RhpcBLASctl::omp_set_num_threads(1L)
    }
    NULL
  })
  doParallel::registerDoParallel(cluster)
  do.call(glmnet::cv.glmnet,c(list(x=x,y=y,parallel=TRUE),dots))
}
