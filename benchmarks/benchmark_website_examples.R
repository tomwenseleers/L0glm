# Reproduce the large examples displayed on the pkgdown website.
# Run from the repository root with R 4.5.0 and one BLAS thread.
# Intermediate fitted objects stay in the ignored codex_sandbox directory.
website_fit_examples <- function(cases=c("million", "single", "multi", "selection")) {
  pkgload::load_all("L0glm", recompile=FALSE, quiet=TRUE)
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
  cache <- "codex_sandbox/website_large_examples"
  output <- "benchmarks/plots/website"
  dir.create(cache, recursive=TRUE, showWarnings=FALSE)
  dir.create(output, recursive=TRUE, showWarnings=FALSE)
  message("Package ", packageVersion("L0glm"), "; BLAS threads: 1")
  for (case in cases) {
    message("Starting ", case, " at ", Sys.time())
    if (case == "million") {
      set.seed(123)
      n <- 500L; p <- 1000000L
      X <- matrix(rnorm(n*p), n, p)
      active <- sample.int(p, 60L)
      beta <- numeric(p); beta[active] <- runif(60L, -.8, .8)
      set.seed(1)
      y <- as.numeric(X[, active] %*% beta[active]) + rnorm(n)
      timings <- numeric(3L)
      for (run in seq_along(timings)) {
        gc()
        timings[run] <- system.time(fit <- L0glm.fit(X, y, family=gaussian()))[["elapsed"]]
        if (run > 1L) stopifnot(identical(coef(fit), B))
        B <- coef(fit)
      }
      elapsed <- median(timings)
      scores <- L0glm.metrics(fit, beta_true=beta)
      scores$seconds <- elapsed
      scores$screening_threads <- fit$screening$nthreads
      write.csv(scores, file.path(output, "million-metrics.csv"), row.names=FALSE)
      write.csv(data.frame(run=seq_along(timings), seconds=timings),
        file.path(output, "million-timings.csv"), row.names=FALSE)
      message("Million-predictor preset: ", elapsed, " s; TP=", scores$TP,
        ", FP=", scores$FP, ", screen=", fit$screening$method)
      saveRDS(list(coefficients=B, beta_true=beta, metrics=scores,
        screening=fit$screening), file.path(cache, "million-summary.rds"))
      elapsed_path <- system.time(path <- L0glm.tune(X, y, family=gaussian(),
        tune="gic", nlambda=10L, prop_range=.25, plot=FALSE,
        printbest=FALSE, coefs.true=beta))[["elapsed"]]
      path$fit_best$X <- NULL
      if (!is.null(path$fit_1se)) path$fit_1se$X <- NULL
      saveRDS(path, file.path(cache, "million-path.rds"))
      write.csv(data.frame(seconds=elapsed_path, lambda_best=path$lambda_best,
        lambda_preset=path$lambda_preset), file.path(output, "million-tuning.csv"),
        row.names=FALSE)
      rm(X, y, fit, path, B, beta); gc()
    } else if (case %in% c("single", "multi")) {
      sim <- simul(method=if (case == "single") 2L else 4L,
        n=10000L, p=10000L, n_nonzero=1000L, seed=2L)
      message(case, ": ", nrow(sim$X), " rows, ", ncol(sim$X),
        " predictors, ", ncol(as.matrix(sim$y)), " channels")
      pieces <- lapply(c("L0glm", "lasso", "relaxed_lasso"), function(method) {
        destination <- file.path(cache, paste0(case, "-", method, ".rds"))
        if (file.exists(destination)) return(readRDS(destination))
        message("Fitting ", case, "/", method, " at ", Sys.time())
        result <- sparseglm(sim, family=poisson("identity"), method=method,
          tune="gic", nlambda=10L, nonnegative=TRUE, intercept=FALSE,
          approximation=if (case == "single") "wls" else "wls_grouped",
          lambda=list(L0glm=log(ncol(sim$X))*log(log(nrow(sim$X)))),
          options=list(L0glm=list(relaxed=TRUE, ncores=8L)))
        stopifnot(all(result$timings$status == "ok"))
        saveRDS(result, destination)
        message("Finished ", case, "/", method, ": ", result$timings$seconds, " s")
        result
      })
      result <- pieces[[1L]]
      result$fits <- do.call(c, lapply(pieces, `[[`, "fits"))
      result$tuning <- do.call(c, lapply(pieces, `[[`, "tuning"))
      result$timings <- do.call(rbind, lapply(pieces, `[[`, "timings"))
      scores <- L0glm.metrics(result)
      write.csv(scores, file.path(output, paste0(case, "-metrics.csv")), row.names=FALSE)
      saveRDS(result, file.path(cache, paste0(case, "-comparison.rds")))
      print(scores[c("method", "seconds", "TP", "FP", "FN", "F1_score", "MCC")])
      rm(sim, result, scores, pieces); gc()
    } else if (case == "selection") {
      sim <- simul("normal", n=1000L, p=5000L, n_nonzero=100L, seed=42L,
        options=list(beta_range=c(-2, 2)))
      elapsed <- system.time(fit <- L0glm.fit(sim$X, sim$y,
        family=gaussian(), lambda="gic", relaxed=TRUE))[["elapsed"]]
      scores <- L0glm.metrics(fit, beta_true=sim$beta_true); scores$seconds <- elapsed
      write.csv(scores, file.path(output, "selection-metrics.csv"), row.names=FALSE)
      gic <- L0glm.tune(sim$X, sim$y, family=gaussian(), tune="gic",
        nlambda=10L, prop_range=.25, relaxed=TRUE, plot=FALSE, printbest=FALSE,
        coefs.true=sim$beta_true)
      set.seed(10)
      folds <- sample(rep(1:5, length.out=nrow(sim$X)))
      cv <- L0glm.tune(sim$X, sim$y, family=gaussian(), tune="cv",
        nlambda=10L, prop_range=.25, nfolds=5L, foldid=folds, relaxed=TRUE,
        plot=FALSE, printbest=FALSE, coefs.true=sim$beta_true)
      saveRDS(list(gic=gic, cv=cv, beta_true=sim$beta_true),
        file.path(cache, "selection-paths.rds"))
      write.csv(data.frame(lambda_best=gic$lambda_best, lambda_preset=gic$lambda_preset,
        cv_best=cv$lambda_best, cv_1se=cv$lambda_1se),
        file.path(output, "selection-tuning.csv"), row.names=FALSE)
      rm(sim, fit, gic, cv, scores); gc()
    }
    message("Completed ", case, " at ", Sys.time())
  }
  writeLines(sub("[[:blank:]]+$", "", capture.output(sessionInfo())),
    file.path(output, "sessionInfo.txt"))
  invisible(NULL)
}

if (sys.nframe() == 0L) website_fit_examples()
