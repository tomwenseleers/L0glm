# Reproduce the GLM data examples in fastglm README and vignettes.
# Source URLs and adaptations are recorded in each result row.
# Source from the repository root or the L0glm package directory.
repo <- if (file.exists("L0glm/DESCRIPTION")) "." else ".."
stopifnot(file.exists(file.path(repo, "L0glm/DESCRIPTION")))
args <- character()
variant <- "current"
root <- file.path(repo, "benchmarks/3_benchmark_fastglm_sparse_and_documented")
out <- file.path(root, "documented_examples")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
pkgload::load_all(file.path(repo, "L0glm"), compile = FALSE, quiet = TRUE)
RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)
source_url <- "https://github.com/jaredhuling/fastglm/blob/master/"
examples <- list()
make <- function(X, y, family, beta = NULL, mu = NULL, mandatory = integer(),
                 source = "README.Rmd", note = "", weights = rep(1, nrow(X)),
                 offset = rep(0, nrow(X)), coefficient_truth = TRUE) {
  list(X = X, y = as.numeric(y), family = family, beta = beta, mu = mu,
    mandatory = mandatory, source = paste0(source_url, source), note = note,
    weights = weights, offset = offset, coefficient_truth = coefficient_truth)
}
for (size in c(5000L, 50000L)) local({
  n <- size
  examples[[paste0("readme_dense_logit_", n)]] <<- function() {
    set.seed(123); p <- if (n == 5000L) 30L else 300L
    X <- matrix(rnorm(n*p), n, p); b <- rnorm(p)*.05
    mu <- plogis(X %*% b); y <- rbinom(n, 1, mu)
    make(X, y, binomial(), b, mu, note = if (n == 5000L) "Exact README example" else
      "README generator scaled to the user's n=50000,p=300 example")
  }
})
examples$readme_speed <- function() {
  set.seed(123); X <- matrix(rnorm(10000*100, sd = 3), 10000, 100)
  b <- numeric(100); b[1:25] <- runif(25, -.1, .1)
  eta <- as.numeric(X %*% b); y <- as.numeric(eta > rnorm(10000))
  make(X, y, binomial(), b, pnorm(eta), coefficient_truth = FALSE,
    note = "Exact README speed DGP: probit latent response, fitted with logit")
}
for (size in c(1000L, 10000L, 100000L, 500000L)) for (kind in c("logit", "poisson")) local({
  n <- size; type <- kind
  examples[[paste0("vignette_", type, "_", n)]] <<- function() {
    set.seed(1); X <- matrix(rnorm(n*30), n, 30); b <- rnorm(30)*.05
    family <- if (type == "logit") binomial() else poisson()
    mu <- if (type == "logit") plogis(X %*% b) else exp(X %*% b + 1)
    y <- if (type == "logit") rbinom(n, 1, mu) else rpois(n, mu)
    make(X, y, family, b, mu, coefficient_truth = type == "logit", source = "vignettes/benchmarks-fastglm.Rmd",
      note = paste("Documented generator and n; seed reset to 1 per case.",
        if (type == "poisson") "As written, fitting omits the true +1 intercept." else ""))
  }
})
for (correct_intercept in c(FALSE, TRUE)) local({
  corrected <- correct_intercept
  examples[[paste0("overview_poisson_", if (corrected) "with_intercept" else "as_written")]] <<- function() {
    set.seed(7); n <- 50000L; X <- matrix(rnorm(n*30), n, 30); b <- rnorm(30)*.05
    mu <- exp(X %*% b + 1); y <- rpois(n, mu)
    if (corrected) { X <- cbind(1, X); b <- c(1, b) }
    make(X, y, poisson(), b, mu, coefficient_truth = corrected, mandatory = if (corrected) 1L else integer(),
      source = "vignettes/fastglm-overview.Rmd", note = if (corrected)
        "Same documented data, with the omitted +1 intercept included in the fitted design" else
        "Exact documented generator and fitted design; true +1 intercept omitted")
  }
})
for (sparse in c(FALSE, TRUE)) local({
  sp <- sparse
  examples[[paste0("vignette_factor_", if (sp) "sparse" else "dense")]] <<- function() {
    set.seed(2); n <- 50000L
    g <- factor(sample.int(100, n, replace = TRUE)); xn <- rnorm(n)
    X <- model.matrix(~ g*xn); b <- rnorm(ncol(X))/4
    mu <- plogis(X %*% b); y <- rbinom(n, 1, mu)
    if (sp) X <- as(X, "CsparseMatrix")
    make(X, y, binomial(), b, mu, mandatory = 1L,
      source = "vignettes/benchmarks-fastglm.Rmd", note = "Exact one-hot interaction example; all coefficients nonzero")
  }
})
for (sparse in c(FALSE, TRUE)) local({
  sp <- sparse
  examples[[paste0("large_data_factor_", if (sp) "sparse" else "dense")]] <<- function() {
    set.seed(1); n <- 5000L
    g <- factor(sample.int(50, n, replace = TRUE)); xnum <- rnorm(n)
    X <- model.matrix(~ g*xnum); b <- rnorm(ncol(X))/4
    mu <- plogis(X %*% b); y <- rbinom(n, 1, mu)
    if (sp) X <- as(X, "CsparseMatrix")
    make(X, y, binomial(), b, mu, mandatory = 1L,
      source = "vignettes/large-data-fastglm.Rmd", note = "Exact 50-level factor generator; seed=1")
  }
})
examples$large_data_dense <- function() {
  set.seed(1); n <- 20000L; X <- matrix(rnorm(n*20), n, 20); b <- rnorm(20)*.05
  mu <- plogis(X %*% b); y <- rbinom(n, 1, mu)
  make(X, y, binomial(), b, mu, source = "vignettes/large-data-fastglm.Rmd",
    note = "Documented big.matrix data generator; ordinary dense fit, seed reset to 1")
}
examples$large_data_streaming <- function() {
  set.seed(1); n <- 10000L; X <- cbind(1, matrix(rnorm(n*4), n, 4))
  b <- c(.2,.5,-.3,.4,-.2); mu <- plogis(X %*% b); y <- rbinom(n, 1, mu)
  make(X, y, binomial(), b, mu, mandatory = 1L, source = "vignettes/large-data-fastglm.Rmd",
    note = "Documented streaming example's data; full-matrix fit, seed reset to 1")
}
examples$large_data_weights_offset <- function() {
  set.seed(1); n <- 8000L; X <- cbind(1, matrix(rnorm(n*3), n, 3))
  offset <- runif(n, -.1, .1); weights <- runif(n, .5, 1.5); b <- c(.2,.4,-.2,.3)
  mu <- exp(X %*% b + offset); y <- rpois(n, mu)
  make(X, y, poisson(), b, mu, mandatory = 1L, offset = offset, weights = weights,
    source = "vignettes/large-data-fastglm.Rmd", note = "Documented weighted/offset Poisson data; seed reset to 1")
}
examples$vignette_gamma_sqrt <- function() {
  set.seed(1); X <- matrix(rnorm(10000*100), ncol = 100)
  b <- numeric(101); b[c(2,4,5,6)] <- c(.25,-.25,.5,-.5)
  eta <- as.numeric(cbind(1, X) %*% b)
  y <- exp(eta + rnorm(10000)) + .1
  make(cbind(1, X), y, Gamma("sqrt"), b, exp(eta+.5)+.1, mandatory = 1L,
    source = "vignettes/fastglm.Rmd", coefficient_truth = FALSE,
    note = "Exact stability example: lognormal response and misspecified Gamma sqrt link")
}
examples$overview_esoph <- function() {
  data(esoph, package = "datasets")
  X <- model.matrix(~ agegp + unclass(tobgp) + unclass(alcgp), data = esoph)
  totals <- esoph$ncases + esoph$ncontrols
  make(X, esoph$ncases/totals, binomial("cloglog"), mandatory = 1L, weights = totals,
    source = "vignettes/fastglm-overview.Rmd", coefficient_truth = FALSE,
    note = "Documented esoph data; equivalent proportions plus trial weights for the L0glm interface")
}
examples$vignette_nb_challenging <- function() {
  set.seed(1); n <- 200L; X <- cbind(1, matrix(rnorm(n*8), n, 8))
  b <- c(4, 1.5,-1.2,1.5,-1.2,1.5,-1.2,1.5,-1.2)
  mu <- exp(X %*% b); y <- MASS::rnegbin(n, mu = mu, theta = 1)
  make(X, y, MASS::negative.binomial(1), b, mu, mandatory = 1L,
    source = "vignettes/nb-convergence-fastglm.Rmd", note = "First documented replicate; theta=1 fixed for matched Pure/L0 comparison")
}
for (seed in 421:423) local({
  simulation_seed <- seed
  examples[[paste0("additional_dense_gaussian_", seed)]] <<- function() {
    set.seed(simulation_seed); X <- matrix(rnorm(400*100), 400, 100)
    b <- runif(100, .3, .6) * rep(c(-1, 1), 50)
    mu <- as.numeric(X %*% b); y <- rnorm(400, mu, .3)
    d <- make(X, y, gaussian(), b, mu,
      note = "Additional dense-truth control, not a fastglm example: n=400,p=100, every slope nonzero, noise SD=0.3")
    d$source <- NA_character_
    d
  }
})
metrics <- function(beta, d) {
  mu <- d$family$linkinv(as.numeric(d$X %*% beta) + d$offset)
  if (d$family$family == "binomial") mu <- pmin(pmax(mu, 1e-8), 1-1e-8) else
    if (d$family$family != "gaussian") mu <- pmax(mu, 1e-8)
  slopes <- setdiff(seq_along(beta), d$mandatory)
  selected <- beta[slopes] != 0
  if (startsWith(d$family$family, "Negative Binomial")) {
    LL <- sum(dnbinom(d$y, mu = mu, size = 1, log = TRUE))
    GIC <- -2*LL + log(length(slopes))*log(log(nrow(d$X)))*sum(beta != 0)
  } else {
    met <- get.metrics(beta, d$X, matrix(d$y), family = d$family,
      weights = matrix(d$weights), offset = matrix(d$offset), always.include = d$mandatory,
      p = max(1L, length(slopes)), mu_floor = 1e-8)
    LL <- met$loglik; GIC <- met$gic
  }
  truth <- if (is.null(d$beta)) rep(NA, length(slopes)) else d$beta[slopes] != 0
  TP <- sum(selected & truth); FP <- sum(selected & !truth); FN <- sum(!selected & truth); TN <- sum(!selected & !truth)
  denom <- sqrt(as.double(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  data.frame(selected = sum(selected), true_support = sum(truth),
    support_fraction = mean(selected), F1 = if (is.na(TP)) NA_real_ else 2*TP/(2*TP+FP+FN),
    MCC = if (!is.na(denom) && denom > 0) (TP*TN-FP*FN)/denom else NA_real_,
    coefficient_RMSE = if (!d$coefficient_truth || is.null(d$beta)) NA_real_ else sqrt(mean((beta-d$beta)^2)),
    prediction_RMSE_true = if (is.null(d$mu)) NA_real_ else sqrt(mean((mu-d$mu)^2)),
    deviance = sum(d$family$dev.resids(d$y, mu, d$weights)), loglik = LL, GIC = GIC)
}
rows <- fits <- list()
resume <- FALSE # Set TRUE only to resume the same package revision and settings.
if (resume && file.exists(file.path(out, "timings_metrics.csv"))) {
  saved <- read.csv(file.path(out, "timings_metrics.csv"), stringsAsFactors = FALSE)
  saved <- saved[is.na(saved$error) | saved$error != "Coefficient dimension mismatch", , drop = FALSE]
  saved <- saved[is.na(saved$seconds) | saved$seconds > 0, , drop = FALSE]
  rows <- split(saved, paste(saved$case, saved$method, sep = "__"))
  fits <- readRDS(file.path(out, "fits.rds"))
  fits <- fits[intersect(names(rows), names(fits))]
}
for (id in names(examples)) {
  d <- examples[[id]]()
  methods <- c(c("fastglmPure_QR", "fastglmPure_LLT", "abess_default"),
               "L0_unpenalized", "L0_default", "L0_gic",
               if (id %in% c("vignette_gamma_sqrt", "vignette_nb_challenging")) "L0_gic_100")
  for (method in methods) {
    key <- paste(id, method, sep = "__")
    if (key %in% names(rows)) next
    cat(format(Sys.time(), "%H:%M:%S"), variant, key, "\n"); flush.console()
    search_max <- NA_integer_
    run <- if (startsWith(method, "fastglmPure")) function() {
      fastglm::fastglmPure(d$X, d$y, family = d$family, weights = d$weights, offset = d$offset,
        method = if (method == "fastglmPure_QR") 0L else 2L, tol = 1e-10, maxit = 25L)
    } else if (method == "abess_default") function() {
      if (id == "readme_dense_logit_50000") stop("Not run: expanded user example, beyond the documented n=5000,p=30 case")
      if (nrow(d$X) >= 100000L && id != "vignette_logit_100000")
        stop("Not run: large abess comparison reserved for the n=100000 logit example")
      if (!d$family$family %in% c("gaussian", "poisson", "binomial") ||
          d$family$family == "binomial" && d$family$link != "logit") stop("abess does not fit this family/link")
      if (any(d$offset != 0) || any(d$weights != 1)) stop("Weighted/offset case excluded from the abess default-interface comparison")
      free <- setdiff(seq_len(ncol(d$X)), d$mandatory)
      result <- callr::r(function(X, y, family, intercept) {
        RhpcBLASctl::blas_set_num_threads(1L)
        RhpcBLASctl::omp_set_num_threads(1L)
        timed <- microbenchmark::microbenchmark(fit <- abess::abess(X, y,
          family = family, fit.intercept = intercept, num.threads = 1L), times = 1L)
        b <- as.numeric(coef(fit, support.size = fit$best.size))
        list(coefficients = if (intercept) b else b[-1L],
          max_support = max(fit$support.size), seconds = median(timed$time)/1e9)
      }, args = list(X = as.matrix(d$X[, free, drop = FALSE]), y = d$y,
        family = d$family$family, intercept = length(d$mandatory) > 0L), timeout = 120)
      search_max <<- result$max_support
      result
    } else function() L0glm.fit(d$X, d$y, family = if (id == "vignette_nb_challenging") mgcv::nb(theta = 1) else d$family,
      weights = d$weights, offset = d$offset, intercept = FALSE, always.include = d$mandatory,
      nonnegative = FALSE, relaxed = FALSE, parallelize = FALSE, nthreads = 1L,
      lambda = switch(method, L0_unpenalized = 0, L0_default = 1, "gic"),
      thresh = if (method == "L0_unpenalized") 0 else 1e-5,
      maxit = if (method == "L0_gic_100") 100L else 25L, tol = 1e-10)
    warnings <- character(); error <- NA_character_; fit <- NULL
    repetitions <- if (method != "abess_default" && nrow(d$X)*(ncol(d$X)-length(d$mandatory)) <= 1e6) 3L else 1L
    gc()
    timing <- withCallingHandlers(tryCatch(microbenchmark::microbenchmark(fit <- run(), times = repetitions),
      error = function(e) { error <<- conditionMessage(e); NULL }),
      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
    row <- data.frame(variant, case = id, method, n = nrow(d$X), p = ncol(d$X)-length(d$mandatory),
      family = d$family$family, link = d$family$link, sparse = inherits(d$X, "sparseMatrix"),
      seconds = if (is.null(timing)) NA_real_ else if (!is.null(fit$seconds)) fit$seconds else median(timing$time)/1e9,
      repetitions = repetitions,
      converged = if (is.null(fit$convergence)) if (is.null(fit$converged)) NA else fit$converged else all(fit$convergence$converged),
      iterations = if (is.null(fit$convergence)) if (is.null(fit$iter)) NA_integer_ else fit$iter else max(fit$convergence$iterations),
      abess_max_support = search_max, warnings = paste(unique(warnings), collapse = " | "), error,
      source = d$source, note = d$note)
    beta <- if (is.null(fit)) NULL else as.numeric(fit$coefficients)
    if (!is.null(beta) && length(beta) != ncol(d$X)) row$error <- "Coefficient dimension mismatch"
    if (!is.null(beta) && length(beta) == ncol(d$X) && all(is.finite(beta))) {
      met <- tryCatch(metrics(beta, d), error = function(e) { row$error <<- paste("scoring:", conditionMessage(e)); NULL })
      if (!is.null(met)) row <- cbind(row, met)
    }
    rows[[key]] <- row; fits[[key]] <- list(beta = beta, error = error)
    columns <- unique(unlist(lapply(rows, names)))
    aligned <- lapply(rows, function(x) { x[setdiff(columns, names(x))] <- NA; x[columns] })
    write.csv(do.call(rbind, aligned), file.path(out, "timings_metrics.csv"), row.names = FALSE)
    saveRDS(fits, file.path(out, "fits.rds"))
    print(row[, intersect(c("method", "seconds", "selected", "true_support", "abess_max_support", "error"), names(row))])
  }
  rm(d); gc()
}
cat("DOCUMENTED EXAMPLES COMPLETE", variant, "\n")
