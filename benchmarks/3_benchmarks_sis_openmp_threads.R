# OpenMP SIS and end-to-end ISIS on the original million-predictor example.
# Install the current package first. X alone occupies about 4 GB.
library(L0glm)
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)

# Editable settings; no environment-variable configuration.
n <- 500L
p <- 1000000L
n_nonzero <- 60L
thread_counts <- c(1L,2L,4L,8L)
sis_repetitions <- 5L
isis_repetitions <- 3L
nsis <- 2L*floor(n/log(n))
root_dir <- if (file.exists("L0glm/DESCRIPTION")) "." else ".."
out_dir <- file.path(root_dir,"benchmarks/3_benchmarks_sis_openmp_results")
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)

set.seed(123)
X <- matrix(rnorm(n*p),n,p)
active <- sample.int(p,n_nonzero)
beta_true <- numeric(p);beta_true[active] <- runif(n_nonzero,-.8,.8)
y_true <- as.numeric(X[,active,drop=FALSE]%*%beta_true[active])
set.seed(1);y <- y_true+rnorm(n)

# Check every statistic, sign and rank before timing the top-only production path.
full_reference <- sis_glm(X,y,nsis=nsis,return_stats=TRUE,nthreads=1L)
for (threads in thread_counts) {
  actual <- sis_glm(X,y,nsis=nsis,return_stats=TRUE,nthreads=threads)
  stopifnot(actual$nthreads==threads,
    identical(actual[names(actual)!="nthreads"],full_reference[names(full_reference)!="nthreads"]))
}
rm(actual,full_reference)
reference <- sis_glm(X,y,nsis=nsis,return_stats=FALSE,nthreads=1L)
records <- list()
record <- function(stage,repeat_id,threads,used,seconds,metrics=NULL) {
  row <- data.frame(stage,repeat_id,requested_threads=threads,actual_threads=used,seconds,
    gic=NA_real_,F1_score=NA_real_,MCC=NA_real_,selected=NA_integer_,
    max_coefficient_difference=NA_real_)
  if (!is.null(metrics)) {
    row[c("gic","F1_score","MCC")] <- metrics[c("gic","F1_score","MCC")]
    row$selected <- metrics$support_size_estimated
    row$max_coefficient_difference <- 0
  }
  records[[length(records)+1L]] <<- row
  write.csv(do.call(rbind,records),file.path(out_dir,"timings.csv"),row.names=FALSE)
  print(row,row.names=FALSE,digits=6)
  flush.console()
}

# Rotate execution order across repetitions to reduce ordering/temperature bias.
for (repeat_id in seq_len(sis_repetitions)) {
  run_order <- thread_counts[1L+(seq_along(thread_counts)+repeat_id-2L)%%length(thread_counts)]
  for (threads in run_order) {
    invisible(gc())
    seconds <- system.time(fit <- sis_glm(X,y,nsis=nsis,return_stats=FALSE,nthreads=threads))[["elapsed"]]
    stopifnot(identical(fit[names(fit)!="nthreads"],reference[names(reference)!="nthreads"]))
    record("SIS",repeat_id,threads,fit$nthreads,seconds)
  }
}

reference_fit <- NULL
for (repeat_id in seq_len(isis_repetitions)) {
  run_order <- if (repeat_id%%2L) thread_counts else rev(thread_counts)
  for (threads in run_order) {
    invisible(gc())
    seconds <- system.time(fit <- L0glm.fit(X,y,screen.nthreads=threads,nthreads=1L))[["elapsed"]]
    if (is.null(reference_fit)) reference_fit <- fit
    stopifnot(identical(coef(fit),coef(reference_fit)),
      identical(fit$screening$history,reference_fit$screening$history),
      fit$screening$nthreads==threads,inherits(coef(fit),"sparseMatrix"),
      identical(dim(coef(fit)),c(p,1L)))
    metrics <- L0glm.metrics(fit,beta_true=beta_true)
    record("ISIS + L0glm",repeat_id,threads,fit$screening$nthreads,seconds,metrics)
  }
}
table <- do.call(rbind,records)
summary <- do.call(rbind,lapply(split(table,list(table$stage,table$requested_threads),drop=TRUE),function(x)
  data.frame(stage=x$stage[1L],threads=x$requested_threads[1L],repetitions=nrow(x),
    median_seconds=median(x$seconds),min_seconds=min(x$seconds),max_seconds=max(x$seconds))))
summary <- summary[order(summary$stage,summary$threads),]
summary$speedup <- NA_real_
for (stage in unique(summary$stage)) {
  keep <- summary$stage==stage
  baseline <- summary$median_seconds[keep&summary$threads==1L]
  summary$speedup[keep] <- baseline/summary$median_seconds[keep]
}
write.csv(summary,file.path(out_dir,"summary.csv"),row.names=FALSE)
write.csv(reference_fit$screening$history,file.path(out_dir,"isis_history.csv"),row.names=FALSE)
writeLines(sub("[[:blank:]]+$","",capture.output(sessionInfo())),file.path(out_dir,"sessionInfo.txt"))
print(summary,row.names=FALSE,digits=5)
cat("ALL STATISTICS, RANKINGS, COEFFICIENTS AND ISIS HISTORIES MATCH EXACTLY ACROSS THREAD COUNTS\n")
