# These component checks isolate IRLS/tuning from final support refinement.
test_that("scalar preclustering aggregates offsets and refits original channels", {
  skip_if_not_installed("mgcv")
  set.seed(86)
  n <- 80L
  X <- matrix(runif(n * 8L), n, 8L)
  groups <- c(1L, 1L, 2L, 2L)
  w <- seq(0.8, 1.2, length.out=n)
  O <- matrix(rep(c(2, 3, 4, 5), each=n), n)
  mu <- O + outer(X[, 1], c(18, 21, 25, 30)) + outer(X[, 3], c(10, 11, 12, 14))
  counts <- function(size) matrix(rnbinom(length(mu), mu=mu, size=size), n)
  positive <- matrix(rgamma(length(mu), shape=20, scale=mu/20), n)
  censored <- positive
  censored[seq_len(10L), ] <- 0
  models <- list(
    quasi=list(quasipoisson("identity"), counts(mu/2), 3),
    nb1=list(nb1(3, "identity"), counts(mu/2), NULL),
    nb2=list(mgcv::nb(theta=20, link="identity"), counts(20), NULL),
    gamma=list(Gamma("identity"), positive, 0.05),
    gamma_scale=list(gammaScale(0.5, "identity"), positive, NULL),
    gamma_censored=list(gammaCensored(0.05, "identity"), censored, NULL),
    inverse_gaussian=list(inverse.gaussian("identity"), positive, 0.01)
  )
  for (x in models) {
    family <- x[[1L]]
    Y <- x[[2L]]
    args <- list(X=X, y=Y, family=family, lambda=8, offset=O, weights=w,
      dispersion=x[[3L]], nonnegative=TRUE, parallelize=FALSE, precluster=TRUE,
      control.precluster=list(clusterids=groups))
    reduced <- do.call(l0glmFit, utils::modifyList(c(args, list(cluster.refit=FALSE)), list(local.search=FALSE)))
    full <- do.call(l0glmFit, utils::modifyList(args, list(local.search=FALSE)))
    expect_identical(dim(coef(reduced)), c(8L, 2L))
    expect_identical(dim(coef(full)), c(8L, 4L))
    expect_equal(unname(as.matrix(reduced$offset)), cbind(rowSums(O[, 1:2]), rowSums(O[, 3:4])))
    expect_equal(unname(as.matrix(reduced$y)), cbind(rowSums(Y[, 1:2]), rowSums(Y[, 3:4])))
    expect_equal(unname(as.matrix(full$offset)), O)
    expect_true(all(is.finite(fitted(full))))
    expect_true(is.finite(l0glmMetrics(full)$gic))
    expect_true(all(which(Matrix::rowSums(abs(coef(full))) > 0) %in% full$preclustering$candidates))
  }
})

test_that("Gaussian preclustering retains signed starts and scales case weights", {
  set.seed(13)
  X <- matrix(rnorm(120 * 6L), 120, 6L)
  offset <- cbind(1, 2, 3) + matrix(0, 1L, 3L)
  offset <- offset[rep(1L, nrow(X)), ]
  y <- outer(X[, 2], c(-3, -4, -5)) + offset + matrix(rnorm(360, sd=.2), 120, 3)
  start <- rep(0, ncol(X)); start[2L] <- -1
  w <- seq(0.8, 1.2, length.out=nrow(X))
  args <- list(X=X, y=y, offset=offset, weights=w, family=gaussian(),
    start=start, dispersion=c(1, 2, 4), precluster=TRUE, lambda=6,
    parallelize=FALSE, control.precluster=list(clusterids=c(1L, 1L, 2L)))
  reduced <- do.call(l0glmFit, utils::modifyList(c(args, list(cluster.refit=FALSE)), list(local.search=FALSE)))
  full <- do.call(l0glmFit, utils::modifyList(args, list(local.search=FALSE)))
  expect_equal(unname(as.matrix(reduced$y)), cbind(rowMeans(y[, 1:2]), y[, 3]))
  expect_equal(as.numeric(reduced$prior.weights[, 1]), w * 2)
  expect_equal(as.numeric(reduced$dispersion), c(1.5, 4))
  expect_true(all(coef(full)[2L, ] < 0))
  expect_true(all(coef(full)[-2L, ] == 0))
  expect_error(l0glmFit(local.search=FALSE, X, y, family=gaussian("log"), precluster=TRUE), "identity-link")
})

test_that("identity-link NB2 tuning learns partitions in its training folds", {
  skip_if_not_installed("mgcv")
  set.seed(73)
  n <- 54L
  X <- matrix(runif(n * 6L), n)
  mu <- 2 + outer(X[, 2L], c(18, 16, 2, 1)) + outer(X[, 4L], c(1, 2, 12, 14))
  Y <- matrix(rnbinom(length(mu), mu=mu, size=30), n)
  folds <- rep(1:3, length.out=n)
  fit <- l0glmTune(local.search=FALSE, X, Y, family=mgcv::nb(theta=30, link="identity"), tune="cv",
    offset=2, nonnegative=TRUE, parallelize=FALSE, precluster=TRUE, foldid=folds,
    nlambda=3L, prop_range=.2, control.precluster=list(k=2L),
    plot=FALSE, printbest=FALSE)
  expect_true(all(is.finite(fit$cv$mean)))
  expect_identical(dim(coef(fit$fit_best)), c(6L, 4L))
  expect_length(fit$cv$preclustering, 3L)
  for (i in 1:3) {
    train <- folds != i
    labels <- glmCluster(Y[train, , drop=FALSE], family=poisson(),
      k=2L, linkage="ward", nthreads=1L)$clusterids
    expect_identical(fit$cv$preclustering[[i]]$clusterids, labels)
  }
})
