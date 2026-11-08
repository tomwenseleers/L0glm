# Paired local-search comparisons on the abess Gaussian/logistic designs.
# Run after installing the checkout; all variants use exactly the same data.
benchmark_local_search <- function(output = "benchmarks/plots/website",
                                   repetitions = 10L,
                                   controls = list(auto = list())) {
  stopifnot(requireNamespace("abess", quietly = TRUE))
  rng <- RNGkind()
  on.exit(do.call(RNGkind, as.list(rng)), add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  auc <- function(y, mu) {
    positive <- y == 1
    (sum(rank(mu)[positive]) - sum(positive) * (sum(positive) + 1) / 2) /
      (sum(positive) * sum(!positive))
  }
  records <- list()
  for (model in c("gaussian", "binomial")) for (rho in c(.1, .7)) {
    n <- if (model == "gaussian") 500L else 1000L
    generator <- list(n = n, p = 1000L, rho = rho, cortype = if (rho == .1) 2L else 3L,
      family = model, sigma = if (model == "gaussian") 1 else 0)
    for (seed in seq_len(repetitions)) {
      train <- do.call(abess::generate.data, c(generator, list(support.size = 10L, seed = seed)))
      test <- do.call(abess::generate.data, c(generator, list(beta = train$beta, seed = seed + 100L)))
      data_hash <- digest::digest(list(train$x, train$y, train$beta, test$x, test$y),
        algo = "xxhash64")
      for (variant in c("disabled", names(controls))) {
        seconds <- system.time(fit <- L0glm::l0glmFit(train$x, train$y,
          family = if (model == "gaussian") gaussian() else binomial(), lambda = "gic",
          intercept = TRUE, relaxed = TRUE, parallelize = FALSE, ncores = 1L, screen = "none",
          local.search = if (variant == "disabled") FALSE else "auto",
          control.local = if (variant == "disabled") list() else controls[[variant]]))[["elapsed"]]
        metrics <- L0glm::l0glmMetrics(fit, beta_true = matrix(c(0, train$beta), ncol = 1L),
          mu_floor = 1e-4)
        mu <- as.numeric(predict(fit, test$x))
        records[[length(records) + 1L]] <- data.frame(model, rho, seed, data_hash, variant, seconds,
          GIC = metrics$gic, TP = metrics$TP, FP = metrics$FP, F1 = metrics$F1_score, MCC = metrics$MCC,
          prediction = if (model == "gaussian") sqrt(mean((test$y - mu)^2)) else auc(test$y, mu))
        write.csv(do.call(rbind, records), file.path(output, "local-search-runs.csv"), row.names = FALSE)
      }
    }
  }
  runs <- do.call(rbind, records)
  summary <- merge(aggregate(cbind(GIC, TP, FP, F1, MCC, prediction) ~ model + rho + variant, runs, mean),
    aggregate(seconds ~ model + rho + variant, runs, median), by = c("model", "rho", "variant"))
  write.csv(summary, file.path(output, "local-search.csv"), row.names = FALSE)
  invisible(summary)
}

# Confirm that automatic size limits preserve large deconvolution fits.
benchmark_local_deconvolution <- function(output = "benchmarks/plots/website",
                                          sizes = c(1000L, 10000L)) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  records <- list()
  for (n in sizes) for (multichannel in c(FALSE, TRUE)) {
    sim <- L0glm::simulateSpikeTrain(n = n, p = n, k = n %/% 10L,
      m = if (n == 1000L && multichannel) 8L else NULL,
      multichannel = multichannel, seed = 123L, Plot = FALSE)
    for (mode in list(FALSE, "auto", TRUE)) {
      seconds <- system.time(fit <- L0glm::l0glmFit(sim$X, sim$y,
        family = poisson("identity"), lambda = "gic", nonnegative = TRUE,
        grouped = multichannel, relaxed = TRUE, parallelize = multichannel,
        ncores = if (multichannel) 8L else 1L, screen = "none", local.search = mode,
        control.local = list(maxit = 3L, max.evals = 36L, nthreads = 8L)))[["elapsed"]]
      if (identical(mode, FALSE)) plain <- fit
      if (identical(mode, "auto"))
        stopifnot(!fit$local.search$applied, identical(coef(fit), coef(plain)))
      metrics <- L0glm::l0glmMetrics(fit, beta_true = sim$beta_true, mu_floor = 1e-4)
      stopifnot(metrics$gic <= L0glm::GIC(plain, mu_floor = 1e-4) + 1e-5)
      records[[length(records) + 1L]] <- data.frame(n, p = n, m = NCOL(sim$y),
        mode = as.character(mode), seconds,
        metrics[c("gic", "TP", "FP", "F1_score", "MCC", "within_row_F1_score", "within_row_MCC")])
      write.csv(do.call(rbind, records), file.path(output, "local-search-deconvolution.csv"),
        row.names = FALSE)
    }
  }
  invisible(do.call(rbind, records))
}

# Compare the same search in original and retained-screening coordinates.
# testthat supplies a temporary binding change for the full-coordinate control.
benchmark_local_million <- function(output = "benchmarks/plots/website", repetitions = 3L) {
  stopifnot(requireNamespace("testthat", quietly = TRUE))
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
  set.seed(123, kind="Mersenne-Twister", normal.kind="Inversion", sample.kind="Rejection")
  n <- 500L; p <- 1000000L
  X <- matrix(rnorm(n*p), n, p)
  active <- sample.int(p, 60L)
  beta <- numeric(p); beta[active] <- runif(60L, -.8, .8)
  set.seed(1)
  y <- as.numeric(X[,active] %*% beta[active]) + rnorm(n)
  records <- list(); reference <- NULL
  for (repetition in seq_len(repetitions)) {
    order <- if (repetition %% 2L) c("original", "screened") else c("screened", "original")
    for (coordinates in order) {
      gc()
      seconds <- system.time(fit <- if (coordinates == "original")
        testthat::with_mocked_bindings(L0glm::l0glmFit(X,y),
          L0glm_local_workspace=function(...)NULL, .package="L0glm") else
        L0glm::l0glmFit(X,y))[["elapsed"]]
      search <- fit$local.search
      answer <- list(coefficients=coef(fit), history=search$history, proposals=search$proposals)
      if (is.null(reference)) reference <- answer else
        stopifnot(isTRUE(all.equal(answer, reference, tolerance=1e-10)))
      metrics <- L0glm::l0glmMetrics(fit, beta_true=beta)
      records[[length(records)+1L]] <- data.frame(repetition, coordinates, seconds,
        local_seconds=search$seconds, screen_seconds=search$screen.seconds,
        refit_seconds=search$refit.seconds, evaluations=search$evaluations,
        metrics[c("gic","TP","FP","F1_score","MCC")])
      write.csv(do.call(rbind, records), file.path(output,"local-search-million-runs.csv"),
        row.names=FALSE)
      rm(fit, search, answer)
    }
  }
  invisible(do.call(rbind, records))
}
