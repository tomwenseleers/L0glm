## Compare native Poisson validation with the general driver and check thread invariance.
library(L0glm)
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)

n <- p <- 1000L
m <- 16L
n_nonzero <- 60L
repetitions <- 5L
sim <- simul(n = n, p = p, m = m, n_nonzero = n_nonzero, seed = 201028L)
set.seed(201029)
y.validation <- matrix(rpois(length(sim$y_true), sim$y_true), n, m)
validation <- list(X = sim$X, y = y.validation)
settings <- data.frame(backend = c("General driver", rep("l0nnPois", 4)),
  ncores = c(1L, 1L, 2L, 4L, 8L))

fits <- vector("list", nrow(settings))
results <- lapply(seq_len(nrow(settings)), function(i) {
  args <- list(X = sim$X, y = sim$y, family = poisson("identity"),
    nonnegative = TRUE, lambda = "gic", relaxed = TRUE, group.maxit = 3L,
    openmp = settings$backend[i] == "l0nnPois", parallelize = settings$ncores[i] > 1L,
    ncores = settings$ncores[i], validation = validation)
  fit <- do.call(l0glmFit, args)
  seconds <- replicate(repetitions, system.time(fit <- do.call(l0glmFit, args))[["elapsed"]])
  fits[[i]] <<- fit
  metrics <- l0glmMetrics(fit, beta_true = sim$beta_true)
  data.frame(settings[i, ], seconds = median(seconds),
    validation_deviance = fit$validation$deviance,
    metrics[, c("gic", "F1_score", "MCC")])
})
results <- do.call(rbind, results)
for (fit in fits[-1L]) stopifnot(isTRUE(all.equal(coef(fit), coef(fits[[1L]]), tolerance = 1e-7)))
print(results, row.names = FALSE)
