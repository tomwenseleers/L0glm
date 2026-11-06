# Reproduce the BIC and five-fold-CV LASSO comparisons in the README and website.
# Run from the repository root. Defaults allocate a 500 x 1,000,000 dense design.
benchmark_million_lasso <- function(n=500L, p=1000000L, n_nonzero=60L,
  nlambda=40L, nfolds=5L, output="benchmarks/plots/website") {
  library(L0glm)
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
  set.seed(123, kind="Mersenne-Twister", normal.kind="Inversion", sample.kind="Rejection")
  X <- matrix(rnorm(n*p), n, p)
  active <- sample.int(p, n_nonzero)
  beta <- numeric(p)
  beta[active] <- runif(n_nonzero, -.8, .8)
  set.seed(1)
  y <- as.numeric(X[, active, drop=FALSE] %*% beta[active]) + rnorm(n)
  methods <- c("lasso", "relaxed_lasso")

  comparison <- sparseGLM(X, y, family=gaussian(), method=methods,
    tune="bic", nlambda=nlambda, intercept=FALSE, beta_true=beta, screen="none")
  bic <- l0glmMetrics(comparison)
  stopifnot(all(bic$status == "ok"))
  bic$method <- paste0(methods, "_bic")
  bic$tune <- "BIC path"
  rm(comparison); gc()

  set.seed(42)
  folds <- sample(rep(seq_len(nfolds), length.out=n))
  path_options <- list(control=list(pmax=2*n+20L))
  comparison <- sparseGLM(X, y, family=gaussian(), method=methods,
    tune="cv", nfolds=nfolds, foldid=folds, nlambda=nlambda,
    intercept=FALSE, beta_true=beta, screen="none",
    options=list(lasso=path_options, relaxed_lasso=path_options))
  cv <- l0glmMetrics(comparison)
  stopifnot(all(cv$status == "ok"))
  cv$method <- paste0(methods, "_cv")
  cv$tune <- paste0(nfolds, "-fold CV")
  scores <- rbind(bic, cv)
  dir.create(output, recursive=TRUE, showWarnings=FALSE)
  write.csv(scores, file.path(output, "million-lasso-metrics.csv"), row.names=FALSE)
  fields <- c("method", "tune", "seconds", "gic", "TP", "FP", "FN", "TN", "F1_score", "MCC")
  print(scores[fields], row.names=FALSE)
  invisible(scores)
}

if (sys.nframe() == 0L) benchmark_million_lasso()
