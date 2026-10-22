# Direct GIC-preset fits and native glmnet tuning for the website examples.
# CV uses glmnet's weighted-Gaussian prediction error, with common random folds.
website_deconvolution <- function(case=c("single", "multi"), n=10000L, p=n,
  n_nonzero=1000L, seed=2L, ncores=4L, cache="codex_sandbox/website_large_examples",
  output="benchmarks/plots/website", reuse=TRUE) {
  case <- match.arg(case)
  pkgload::load_all("L0glm", recompile=FALSE, quiet=TRUE)
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
  source("benchmarks/cv_glmnet_parallel.R", local=TRUE)
  dir.create(cache, recursive=TRUE, showWarnings=FALSE)
  dir.create(output, recursive=TRUE, showWarnings=FALSE)
  sim <- simul(method=if (case=="single") 2L else 4L,
    n=n, p=p, n_nonzero=n_nonzero, seed=seed)
  set.seed(42L)
  folds <- sample(rep(seq_len(5L), length.out=n))
  message(case, ": ", n, " rows, ", p, " predictors, ", NCOL(sim$y), " channels")
  checkpoint <- function(name, code) {
    path <- file.path(cache, paste0(case, "-", name, ".rds"))
    if (reuse && file.exists(path)) return(readRDS(path))
    message("Fitting ", case, "/", name, " at ", Sys.time())
    value <- force(code)
    saveRDS(value, path)
    value
  }
  preset <- checkpoint("preset-four-core", {
    times <- numeric(3L)
    for (run in seq_along(times)) {
      gc()
      times[run] <- system.time(fit <- L0glm.fit(sim$X, sim$y,
        family=poisson("identity"), lambda="gic", nonnegative=TRUE,
        intercept=FALSE, relaxed=TRUE, grouped=TRUE,
        parallelize=TRUE, ncores=ncores))[["elapsed"]]
      if (run>1L) stopifnot(identical(coef(fit), previous))
      previous <- coef(fit)
    }
    list(fit=fit, seconds=median(times), repeats=times)
  })
  gic <- lapply(c("lasso", "relaxed_lasso"), function(method) checkpoint(method,
    sparseglm(sim, family=poisson("identity"), method=method,
      tune="gic", nlambda=10L, nonnegative=TRUE, intercept=FALSE,
      approximation=if (case=="single") "wls" else "wls_grouped")))
  stopifnot(all(vapply(gic, function(z) all(z$timings$status=="ok"), logical(1))))
  stopifnot(all(vapply(gic, function(z)
    identical(unname(z$beta_true), unname(as.matrix(sim$beta_true))), logical(1))))
  grid <- gic[[1L]]$tuning[[1L]]$candidates$parameter
  stopifnot(length(grid)==10L, all(is.finite(grid)), all(grid>0),
    all(diff(grid)<0), identical(grid, gic[[2L]]$tuning[[1L]]$candidates$parameter))
  row_weights <- if (case=="single") 1/(as.numeric(sim$y)+.1) else
    1/(rowMeans(sim$y)+.1)
  cv <- lapply(c(FALSE, TRUE), function(relaxed) checkpoint(
    if (relaxed) "relaxed-native-cv" else "lasso-native-cv", {
      gc()
      elapsed <- system.time(fit <- cv_glmnet_parallel(sim$X,
        if (case=="single") as.numeric(sim$y) else as.matrix(sim$y),
        family=if (case=="single") "gaussian" else "mgaussian",
        weights=row_weights, lambda=grid, nfolds=5L, foldid=folds,
        intercept=FALSE, lower.limits=0, standardize.response=FALSE,
        relax=relaxed, gamma=0, type.measure="mse",
        parallel=TRUE, ncores=ncores))[["elapsed"]]
      list(fit=fit, seconds=elapsed, foldid=folds)
    }))
  result <- gic[[1L]]
  keys <- c("L0glm", "lasso_gic", "relaxed_lasso_gic", "lasso_cv", "relaxed_lasso_cv")
  result$fits <- setNames(vector("list", length(keys)), keys)
  result$tuning <- setNames(vector("list", length(keys)), keys)
  timing <- gic[[1L]]$timings[rep(1L, length(keys)), , drop=FALSE]
  timing$method <- keys
  timing$seconds <- c(preset$seconds, vapply(gic, function(z) z$timings$seconds, numeric(1)),
    vapply(cv, `[[`, numeric(1), "seconds"))
  timing$tune <- c("GIC preset", "GIC path", "GIC path", "5-fold CV", "5-fold CV")
  timing$cv_mean <- timing$cv_se <- NA_real_
  timing$path_seconds <- timing$cv_seconds <- NA_real_
  timing$note <- ""
  timing$approximation[1L] <- "none"
  timing$selected_parameter[1L] <- preset$fit$init.lambda
  result$fits[[1L]] <- preset$fit
  result$tuning[[1L]] <- list(parameter="lambda", best=preset$fit$init.lambda,
    scope=if (case=="single") "single" else "grouped", approximation="none")
  coefficients <- function(path, lambda) {
    B <- predict(path, s=lambda, type="coefficients")
    if (is.list(B)) B <- do.call(cbind, B)
    B[-1L, , drop=FALSE]
  }
  for (j in 1:2) {
    result$fits[[j+1L]] <- gic[[j]]$fits[[1L]]
    result$tuning[[j+1L]] <- gic[[j]]$tuning[[1L]]
    timing$selected_parameter[j+1L] <- gic[[j]]$timings$selected_parameter
    timing$tuning_score[j+1L] <- gic[[j]]$timings$tuning_score
    native <- cv[[j]]$fit
    selected <- if (j==1L) native$lambda.min else native$relaxed$lambda.min
    path <- if (j==1L) native$glmnet.fit else native$glmnet.fit$relaxed
    B <- coefficients(path, selected)
    # The CV full-data path must reproduce the GIC fit at its selected lambda.
    stopifnot(isTRUE(all.equal(as.numeric(coefficients(path,
      gic[[j]]$timings$selected_parameter)), as.numeric(coef(gic[[j]]$fits[[1L]])),
      tolerance=1e-10)))
    stopifnot(all(is.finite(B)), min(B)>=0,
      identical(dim(B), dim(as.matrix(sim$beta_true))),
      identical(cv[[j]]$foldid, folds))
    result$fits[[j+3L]] <- L0glm:::L0glm_benchmark_fit(B, result$data)
    curve <- if (j==1L) native else native$relaxed$statlist[["g:0"]]
    index <- which.min(abs(curve$lambda-selected))
    stopifnot(length(curve$lambda)==10L, is.finite(curve$cvm[index]),
      abs(curve$lambda[index]-selected)<1e-12*max(1, selected),
      curve$cvm[index]==min(curve$cvm, na.rm=TRUE))
    result$tuning[[j+3L]] <- list(parameter="lambda", best=selected,
      scope=if (case=="single") "single" else "grouped",
      approximation=timing$approximation[j+3L],
      candidates=data.frame(parameter=curve$lambda, score=curve$cvm, se=curve$cvsd),
      cv=list(foldid=folds, nfolds=5L, loss="weighted Gaussian MSE"))
    timing$selected_parameter[j+3L] <- selected
    timing$tuning_score[j+3L] <- timing$cv_mean[j+3L] <- curve$cvm[index]
    timing$cv_se[j+3L] <- curve$cvsd[index]
  }
  result$timings <- timing
  result$foldid <- folds
  scores <- L0glm.metrics(result)
  timing$tuning_score[1L] <- scores$gic[1L]
  result$timings <- timing
  scores$tuning_score[1L] <- timing$tuning_score[1L]
  stopifnot(all(scores$status=="ok"))
  write.csv(scores, file.path(output, paste0(case,"-metrics.csv")), row.names=FALSE)
  write.csv(data.frame(run=seq_along(preset$repeats), seconds=preset$repeats),
    file.path(output, paste0(case,"-preset-timings.csv")), row.names=FALSE)
  write.csv(data.frame(observation=seq_len(n), fold=folds),
    file.path(output, paste0(case,"-folds.csv")), row.names=FALSE)
  candidates <- do.call(rbind, lapply(2:5, function(j) {
    values <- result$tuning[[j]]$candidates
    values$method <- keys[j]
    values$selected <- values$parameter==result$tuning[[j]]$best
    if (!"se" %in% names(values)) values$se <- NA_real_
    values[c("method", "parameter", "score", "se", "selected")]
  }))
  write.csv(candidates, file.path(output, paste0(case,"-tuning.csv")), row.names=FALSE)
  writeLines(sub("[[:blank:]]+$", "", capture.output(sessionInfo())),
    file.path(output, "sessionInfo.txt"))
  saveRDS(result, file.path(cache, paste0(case,"-comparison.rds")))
  print(scores[c("method", "tune", "seconds", "gic", "TP", "FP", "FN", "F1_score", "MCC")])
  invisible(result)
}
