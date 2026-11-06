# All nine adapters on the homepage data: one timed run per method and screen.
# L0glm uses a GIC preset; competing methods search BIC-selected native paths.
# These are per-method fitting/scoring times; shared input preparation is separate.
library(L0glm)
RhpcBLASctl::blas_set_num_threads(1)
root <- if (file.exists("L0glm/DESCRIPTION")) "." else ".."
output <- file.path(root,"benchmarks/3_benchmarks_sis_vs_SIS_results")
dir.create(output,recursive=TRUE,showWarnings=FALSE)
set.seed(123, kind="Mersenne-Twister", normal.kind="Inversion", sample.kind="Rejection")
n <- 500L;p <- 1000000L
X <- matrix(rnorm(n*p),n,p)
active <- sample.int(p,60L);beta <- numeric(p);beta[active] <- runif(60L,-.8,.8)
set.seed(1);y <- as.numeric(X[,active]%*%beta[active])+rnorm(n)
others <- c('abess','L0Learn','MCP','SCAD','lasso','relaxed_lasso',
  'adaptive_lasso','adaptive_relaxed_lasso')
results <- list()
for (mode in c('SIS','ISIS')) {
  size <- if (mode=='SIS') 100L else NULL
  l0 <- sparseGLM(X,y,method='L0glm',screen=mode,nsis=size,tune='gic',
    intercept=FALSE,nlambda=10,beta_true=beta,
    lambda=list(L0glm=presetLambda('gic',p,n)),options=list(L0glm=list(relaxed=TRUE)))
  competitors <- sparseGLM(X,y,method=others,screen=mode,nsis=size,tune='bic',
    intercept=FALSE,nlambda=10,support.size=0:100,beta_true=beta)
  rows <- rbind(l0glmMetrics(l0),l0glmMetrics(competitors))
  stopifnot(all(rows$status=='ok'))
  print(rows[,c('method','screen','tune','seconds','gic','TP','FP','FN','F1_score','MCC')])
  results[[mode]] <- rows
}
write.csv(do.call(rbind,results),file.path(output,'all_methods.csv'),row.names=FALSE)
