# Reproduce the Gaussian and logistic data designs in the abess power vignette.
# Use one common GIC for support selection and an independently drawn test set.
benchmark_abess_simulations <- function(output = "benchmarks/plots/website",
                                        repetitions = 10L) {
  stopifnot(requireNamespace("abess", quietly = TRUE))
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
  methods <- c("L0glm", "abess", "L0Learn", "lasso", "relaxed_lasso",
    "adaptive_lasso", "adaptive_relaxed_lasso", "MCP", "SCAD")
  auc <- function(y, mu) {
    positive <- y == 1
    m <- sum(positive)
    (sum(rank(mu)[positive]) - m * (m + 1) / 2) / (m * sum(!positive))
  }
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  runs <- list()
  for (model in c("gaussian", "binomial")) {
    gaussian <- model == "gaussian"
    family <- if (gaussian) gaussian() else binomial()
    n <- if (gaussian) 500L else 1000L
    p <- 1000L
    for (rho in c(.1, .7)) {
      generator <- list(n = n, p = p, rho = rho, cortype = if (rho == .1) 2L else 3L,
        family = model, sigma = if (gaussian) 1 else 0)
      for (seed in seq_len(repetitions)) {
        train <- do.call(abess::generate.data, c(generator, list(support.size = 10L, seed = seed)))
        test <- do.call(abess::generate.data, c(generator, list(beta = train$beta, seed = seed + 100L)))
        fit <- L0glm::sparseGLM(train$x, train$y, family = family, method = methods,
          tune = "gic", nlambda = 40L, support.size = 0:60, intercept = TRUE,
          lambda = list(L0glm = L0glm::presetLambda("gic", p, n)),
          options = list(L0glm = list(relaxed = TRUE, ncores = 1L),
            abess = list(thread = 1L)))
        metrics <- L0glm::l0glmMetrics(fit)
        if (any(metrics$status != "ok")) {
          print(metrics[metrics$status != "ok", c("method", "status", "note")])
          stop("A comparator failed for ", model, ", rho=", rho, ", seed=", seed)
        }
        active <- abs(train$beta) > 1e-5
        results <- lapply(seq_along(fit$fits), function(i) {
          f <- fit$fits[[i]]
          beta <- as.numeric(coef(f))[-1L]
          mu <- as.numeric(predict(f, newdata = test$x, type = "response"))
          stopifnot(length(beta) == p, all(is.finite(mu)))
          found <- abs(beta) > 1e-5
          tp <- sum(active & found)
          fp <- sum(!active & found)
          fn <- sum(active & !found)
          tn <- sum(!active & !found)
          mcc_denom <- sqrt(as.double(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
          data.frame(model, rho, seed, method = names(fit$fits)[i],
            seconds = metrics$seconds[i], GIC = metrics$gic[i], TP = tp, FP = fp,
            TPR = tp / sum(active), FPR = fp / sum(!active),
            MCC = if (mcc_denom > 0) (tp * tn - fp * fn) / mcc_denom else 0,
            prediction = if (gaussian) sqrt(mean((test$y - mu)^2)) else auc(test$y, mu))
        })
        runs[[length(runs) + 1L]] <- do.call(rbind, results)
        message(model, " rho=", rho, " replicate ", seed)
      }
    }
  }
  runs <- do.call(rbind, runs)
  write.csv(runs, file.path(output, "abess-simulations-runs.csv"), row.names = FALSE)
  summary <- aggregate(cbind(seconds, GIC, TPR, FPR, MCC, prediction) ~ model + rho + method,
    runs, mean)
  summary$seconds <- aggregate(seconds ~ model + rho + method, runs, median)$seconds
  write.csv(summary, file.path(output, "abess-simulations.csv"), row.names = FALSE)
  writeLines(trimws(capture.output(sessionInfo()), which = "right"),
    file.path(output, "abess-simulations-session.txt"))
  invisible(summary)
}
