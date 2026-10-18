# Run from the repository root or the L0glm package directory.
package_dir <- if (file.exists("L0glm/DESCRIPTION")) "L0glm" else "."
pkgload::load_all(package_dir, quiet = TRUE)
root_dir <- if (package_dir == "L0glm") "." else ".."
output_dir <- file.path(root_dir, "codex_sandbox", "benchmark_methods", "examples")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Editable settings; timings include path fitting, scoring and CV when requested.
nlambda <- 12L
nfolds <- 3L
set.seed(123)
n <- 150L
p <- 12L
X <- matrix(rnorm(n * p), n, p)
beta <- c(2, -1.5, 1, rep(0, p - 3L))
cases <- list(
  gaussian = list(X = X, y = as.numeric(X %*% beta + rnorm(n)),
                  beta_true = beta, family = gaussian()),
  binomial = list(X = X, y = rbinom(n, 1, plogis(X %*% beta / 2)),
                  beta_true = beta / 2, family = binomial())
)

# Sparse multichannel design with both row sparsity and within-row sparsity.
A <- Matrix::rsparsematrix(n, p, density = .2)
B <- matrix(0, p, 4L)
B[1L, ] <- c(2, -2, 0, 1)
B[2L, ] <- c(0, 1, -1, 0)
B[3L, ] <- c(-2, 0, 1, 0)
cases$multichannel <- list(X = A, y = as.matrix(A %*% B) +
  matrix(rnorm(n * 4L, sd = .5), n, 4L), beta_true = B, family = gaussian())

# Unsupported identity-Poisson backends use explicitly requested WLS; all
# candidates are nevertheless scored using the original Poisson likelihood.
P <- abs(A)
positive_beta <- abs(beta)
cases$identity_poisson <- list(X = P,
  y = rpois(n, 1 + as.numeric(P %*% positive_beta)), beta_true = positive_beta,
  family = poisson(identity), offset = 1, nonnegative = TRUE)

if (requireNamespace("mgcv", quietly = TRUE) &&
    requireNamespace("survival", quietly = TRUE)) {
  event <- rexp(n, exp(as.numeric(X %*% beta / 2)))
  censor <- rexp(n, .35)
  cases$cox <- list(X = X,
    y = survival::Surv(pmin(event, censor), event <= censor),
    beta_true = beta / 2, family = mgcv::cox.ph())
}
if (requireNamespace("mgcv", quietly = TRUE)) {
  M <- matrix(0, p, 3L)
  M[1L, ] <- c(1, -1, 0)
  M[2L, ] <- c(0, .8, -.8)
  eta <- X %*% M
  probability <- exp(eta - apply(eta, 1L, max))
  probability <- probability / rowSums(probability)
  counts <- t(vapply(seq_len(n), function(i)
    as.vector(rmultinom(1L, 4L, probability[i, ])), integer(3L)))
  cases$multinomial <- list(X = X, y = counts, beta_true = M,
                            family = mgcv::multinom(K = 2L))
}

scores <- paths <- list()
for (name in names(cases)) for (criterion in c("gic", "cv")) {
  message(name, ": ", criterion)
  data <- cases[[name]]
  nonnegative <- isTRUE(data$nonnegative)
  result <- sparseglm(data, method = "all", tune = criterion,
    nlambda = nlambda, nfolds = nfolds, support.size = 0:8,
    nonnegative = nonnegative, approximation = if (nonnegative) "wls" else "none")
  table <- L0glm.metrics(result)
  stopifnot(!any(table$status == "error"))
  key <- paste(name, criterion, sep = "_")
  scores[[key]] <- cbind(case = name, table)
  for (method in names(result$tuning)) {
    path <- result$tuning[[method]]
    if (is.null(path)) next
    candidates <- path$candidates
    if (!"se" %in% names(candidates)) candidates$se <- NA_real_
    paths[[paste(key, method, sep = "_")]] <- cbind(case = name,
      method = method, tune = criterion, parameter_type = path$parameter, candidates)
  }
  saveRDS(result, file.path(output_dir, paste0(key, ".rds")))
  if (name == "gaussian" && criterion == "gic") {
    png(file.path(output_dir, "gaussian_coefficients.png"),
        width = 1800, height = 2200, res = 160)
    plot(result)
    dev.off()
  }
}
# Special-family scores add model-specific columns; align the tables by name.
columns <- unique(unlist(lapply(scores, names)))
scores <- lapply(scores, function(x) {
  for (column in setdiff(columns, names(x))) x[[column]] <- NA
  x[columns]
})
write.csv(do.call(rbind, scores), file.path(output_dir, "metrics.csv"), row.names = FALSE)
write.csv(do.call(rbind, paths), file.path(output_dir, "paths.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
