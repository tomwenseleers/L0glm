# Extended sparse-truth grid. Settings are ordinary R values.
# Source from the repository root or the L0glm package directory.
repo <- if (file.exists("L0glm/DESCRIPTION")) "." else ".."
stopifnot(file.exists(file.path(repo, "L0glm/DESCRIPTION")))
args <- character()
variant <- "current"
root <- file.path(repo, "benchmarks/3_benchmark_fastglm_sparse_and_documented")
out <- file.path(root, "sparse_truth")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
pkgload::load_all(file.path(repo, "L0glm"), compile = FALSE, quiet = TRUE)
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)
families <- c("gaussian_identity", "binomial_logit", "binomial_probit", "binomial_cloglog",
              "poisson_log", "Gamma_log", "inverse_gaussian_log", "nb_log")
shapes <- data.frame(shape = c("tiny_tall", "tall", "large_tall", "broad_tall", "square", "wide", "large_wide"),
  n = c(1000L, 10000L, 50000L, 10000L, 400L, 150L, 400L),
  p = c(25L, 100L, 300L, 800L, 400L, 600L, 2400L), k = c(4L, 8L, 12L, 40L, 12L, 8L, 16L))
cases <- do.call(rbind, lapply(c("dense", "sparse"), function(storage) transform(shapes, storage = storage)))
cases$rho <- 0
cases <- rbind(cases,
  data.frame(shape = c("correlated_tall", "correlated_wide", "weak_tall"),
    n = c(4000L, 200L, 10000L), p = c(200L, 1000L, 800L), k = c(10L, 10L, 80L),
    storage = "dense", rho = c(.9, .9, 0)))
cases$case <- paste(cases$shape, cases$storage, sep = "_")
write.csv(cases, file.path(out, "cases.csv"), row.names = FALSE)
family_for <- function(name) switch(name,
  gaussian_identity = gaussian(), binomial_logit = binomial(), binomial_probit = binomial("probit"),
  binomial_cloglog = binomial("cloglog"), poisson_log = poisson(), Gamma_log = Gamma("log"),
  inverse_gaussian_log = inverse.gaussian("log"), nb_log = MASS::negative.binomial(2))
design <- function(n, p, storage, rho) {
  if (storage == "sparse") return(Matrix::rsparsematrix(n, p, density = .05)/sqrt(.05))
  X <- matrix(rnorm(n*p), n, p)
  if (rho) for (j in 2:p) X[, j] <- rho*X[, j-1L] + sqrt(1-rho^2)*X[, j]
  X
}
generate <- function(i, name) {
  cfg <- cases[i, ]
  set.seed(6100L + i)
  X <- design(cfg$n, cfg$p, cfg$storage, cfg$rho)
  test <- design(1000L, cfg$p, cfg$storage, cfg$rho)
  supported <- sample(which(colSums(as.matrix(X)^2) > 0), cfg$k)
  b <- numeric(cfg$p)
  b[supported] <- runif(cfg$k, if (cfg$shape == "weak_tall") .025 else .2,
    if (cfg$shape == "weak_tall") .07 else .5) * rep(c(-1, 1), length.out = cfg$k)
  intercept <- if (startsWith(name, "binomial")) -.3 else .2
  family <- family_for(name)
  mu <- family$linkinv(as.numeric(intercept + X %*% b))
  mu_test <- family$linkinv(as.numeric(intercept + test %*% b))
  set.seed(16100L + 100L*i + match(name, families))
  y <- switch(name, gaussian_identity = rnorm(cfg$n, mu, .8),
    binomial_logit = rbinom(cfg$n, 1, mu), binomial_probit = rbinom(cfg$n, 1, mu),
    binomial_cloglog = rbinom(cfg$n, 1, mu), poisson_log = rpois(cfg$n, mu),
    Gamma_log = rgamma(cfg$n, shape = 4, rate = 4/mu),
    inverse_gaussian_log = statmod::rinvgauss(cfg$n, mean = mu, shape = 5),
    nb_log = rnbinom(cfg$n, mu = mu, size = 2))
  list(X = X, test = test, y = y, mu = mu, mu_test = mu_test,
       beta = c(intercept, b), family = family)
}
bound_mu <- function(mu, name) {
  if (startsWith(name, "binomial")) return(pmin(pmax(mu, 1e-8), 1-1e-8))
  if (name != "gaussian_identity") pmax(mu, 1e-8) else mu
}
score <- function(beta, sim, name, identifiable = TRUE) {
  estimated <- beta[-1L] != 0; truth <- sim$beta[-1L] != 0
  TP <- sum(estimated & truth); FP <- sum(estimated & !truth)
  FN <- sum(!estimated & truth); TN <- sum(!estimated & !truth)
  denom <- sqrt(as.double(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  mu <- bound_mu(sim$family$linkinv(as.numeric(beta[1L] + sim$X %*% beta[-1L])), name)
  test <- bound_mu(sim$family$linkinv(as.numeric(beta[1L] + sim$test %*% beta[-1L])), name)
  dev <- sum(sim$family$dev.resids(sim$y, mu, rep(1, length(mu))))
  dev_true <- sum(sim$family$dev.resids(sim$mu, mu, rep(1, length(mu))))
  if (name == "nb_log") {
    LL <- sum(dnbinom(sim$y, size = 2, mu = mu, log = TRUE))
    GIC <- -2*LL + preset.lambda("gic", ncol(sim$X), nrow(sim$X))*(sum(estimated)+1L)
  } else {
    met <- get.metrics(beta, cbind(1, sim$X), matrix(sim$y), family = sim$family,
      beta_true = sim$beta, always.include = 1L, mu_floor = 1e-8)
    LL <- met$loglik; GIC <- met$gic
  }
  if (!identifiable) LL <- GIC <- NA_real_
  data.frame(selected = sum(estimated), TP, FP, FN, TN,
    F1 = 2*TP/(2*TP+FP+FN), MCC = if (denom > 0) (TP*TN-FP*FN)/denom else NA_real_,
    coefficient_RMSE = sqrt(mean((beta-sim$beta)^2)),
    prediction_RMSE_train = sqrt(mean((mu-sim$mu)^2)),
    prediction_RMSE_test = sqrt(mean((test-sim$mu_test)^2)),
    deviance = dev, deviance_true = dev_true, loglik = LL, GIC = GIC)
}
rows <- fits <- list()
for (i in seq_len(nrow(cases))) for (name in if (length(args) > 1L) args[-1L] else families) {
  cfg <- cases[i, ]
  sim <- generate(i, name)
  wide <- cfg$p + 1L >= cfg$n
  methods <- c("fastglmPure", if (!wide) "L0_unpenalized", "L0_gic")
  for (method in methods) {
    id <- paste(cfg$case, name, method, sep = "__")
    cat(format(Sys.time(), "%H:%M:%S"), variant, id, "\n"); flush.console()
    run <- if (method == "fastglmPure") function() fastglm::fastglmPure(
      cbind(1, sim$X), sim$y, family = sim$family,
      method = if (wide && cfg$storage == "dense") 0L else 2L, tol = 1e-10, maxit = 25L) else
      function() L0glm.fit(sim$X, sim$y, family = if (name == "nb_log") mgcv::nb(theta = 2) else sim$family,
        intercept = TRUE, nonnegative = FALSE, lambda = if (method == "L0_gic") "gic" else 0,
        thresh = if (method == "L0_gic") 1e-5 else 0,
        relaxed = FALSE, maxit = 25L, tol = 1e-10, nthreads = 1L, parallelize = FALSE)
    fit <- NULL; error <- NA_character_; warnings <- character()
    reps <- if (cfg$n*cfg$p <= 1e6 && !wide) 3L else 1L
    gc()
    timed <- withCallingHandlers(tryCatch(microbenchmark::microbenchmark(fit <- run(), times = reps),
      error = function(e) { error <<- conditionMessage(e); NULL }),
      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
    beta <- if (is.null(fit)) NULL else as.numeric(fit$coefficients)
    row <- data.frame(variant, case = cfg$case, family = name, method,
      n = cfg$n, p = cfg$p, true_support = cfg$k, storage = cfg$storage,
      seconds = if (is.null(timed)) NA_real_ else median(timed$time)/1e9, repetitions = reps,
      full_design_tall = !wide, iterations = if (is.null(fit)) NA_integer_ else
        if (is.null(fit$convergence)) fit$iter else max(fit$convergence$iterations),
      converged = if (is.null(fit)) FALSE else if (is.null(fit$convergence)) isTRUE(fit$converged) else
        all(fit$convergence$converged), rank = if (is.null(fit$rank)) NA_integer_ else fit$rank,
      warnings = paste(unique(warnings), collapse = " | "), error)
    if (!is.null(beta) && all(is.finite(beta))) {
      metrics <- tryCatch(score(beta, sim, name, !wide || method != "fastglmPure"),
        error = function(e) { row$error <<- paste("scoring:", conditionMessage(e)); NULL })
      if (!is.null(metrics)) row <- cbind(row, metrics)
    }
    rows[[id]] <- row
    fits[[id]] <- if (is.null(beta)) list(error = error) else
      list(beta = beta, convergence = fit$convergence, dispersion = fit$dispersion)
    columns <- unique(unlist(lapply(rows, names)))
    aligned <- lapply(rows, function(x) { x[setdiff(columns, names(x))] <- NA; x[columns] })
    write.csv(do.call(rbind, aligned), file.path(out, "timings_metrics.csv"), row.names = FALSE)
    saveRDS(fits, file.path(out, "fits.rds"))
    print(row[, intersect(c("method", "seconds", "selected", "F1", "converged", "error"), names(row))])
  }
  rm(sim); gc()
}
cat("GRID COMPLETE", variant, "\n")
