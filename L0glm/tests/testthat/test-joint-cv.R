test_that("joint CV reproduces independent weighted multichannel training fits", {
  set.seed(7091)
  n <- 48L
  X <- matrix(rnorm(n * 7L), n)
  y <- cbind(1 + 2 * X[, 1] + rnorm(n), -X[, 3] + rnorm(n))
  w <- cbind(seq(.5, 2, length.out = n), seq(2, .5, length.out = n))
  w[1, 1] <- w[2, 2] <- 0
  off <- seq(0, .2, length.out = n)
  folds <- rep(1:3, length.out = n)
  args <- list(family = gaussian(), weights = w, offset = off,
    intercept = TRUE, dispersion = c(1, -2), start = 1,
    relaxed = FALSE, grouped = FALSE, local.search = FALSE, parallelize = FALSE,
    screen = "SIS", nsis = 4L, screen.nthreads = 1L)
  out <- do.call(l0glmTune, c(list(X = X, y = y, tune = "cv", foldid = folds,
    maxit = c(1L, 4L), nlambda = 3L, prop_range = .1,
    printbest = FALSE, plot = FALSE), args))
  expect_identical(out$maxits, c(1L, 4L))
  expect_equal(nrow(out$metrics), 6L)
  expect_equal(out$cv$maxit, out$metrics$maxit)
  expected <- out$cv$deviance
  for (k in 1:3) {
    train <- folds != k
    a <- args
    a$weights <- w[train, , drop = FALSE]
    a$offset <- off[train]
    for (i in seq_along(out$lambdas)) {
      fit <- do.call(l0glmFit, c(list(X = X[train, ], y = y[train, ],
        lambda = out$cv$lambda_train[i, k], maxit = out$cv$maxit[i]), a))
      pred <- predict(fit, X[!train, ], offset = off[!train])
      expected[i, k] <- sum(w[!train, ] * (y[!train, ] - pred)^2) / sum(w[!train, ])
    }
  }
  expect_equal(out$cv$deviance, expected, tolerance = 1e-10)
  proportions <- out$cv$fold_weights / sum(out$cv$fold_weights)
  expect_equal(out$cv$mean, drop(expected %*% proportions), tolerance = 1e-10)
  expect_equal(out$cv$se, sqrt(drop((expected - out$cv$mean)^2 %*% proportions) / 2),
    tolerance = 1e-10)
  best <- which.min(out$cv$mean)
  expect_equal(out$lambda_best, out$lambdas[best])
  expect_equal(out$maxit_best, out$cv$maxit[best])
  eligible <- which(out$cv$mean <= out$cv$mean[best] + out$cv$se[best])
  one <- eligible[order(-out$lambdas[eligible], -out$cv$maxit[eligible])][1L]
  expect_equal(out$lambda_1se, out$lambdas[one])
  expect_equal(out$maxit_1se, out$cv$maxit[one])
  replay <- do.call(l0glmFit, c(list(X = X, y = y, lambda = out$lambda_best), out$fit_args))
  expect_equal(coef(replay), coef(out$fit_best))
  replay_args <- out$fit_args
  replay_args$maxit <- out$maxit_1se
  replay <- do.call(l0glmFit, c(list(X = X, y = y, lambda = out$lambda_1se), replay_args))
  expect_equal(coef(replay), coef(out$fit_1se))
})

test_that("iteration tuning supports native Poisson, quasi, Gamma and binomial fits", {
  set.seed(7092)
  X <- matrix(runif(48L * 4L), 48L)
  mu <- 2 + outer(X[, 1], c(6, 3))
  cases <- list(
    list(family = poisson("identity"), y = matrix(rpois(length(mu), mu), 48L),
      nonnegative = TRUE, offset = 2, grouped = TRUE, openmp = TRUE),
    list(family = quasipoisson("identity"), y = matrix(rpois(length(mu), mu), 48L),
      nonnegative = TRUE, offset = 2, grouped = TRUE),
    list(family = Gamma("log"), y = rgamma(48L, 4, scale = exp(X[, 1])/4)),
    list(family = inverse.gaussian("log"), y = exp(X[, 1] + rnorm(48L, sd = .1)),
      dispersion = .01),
    list(family = nb1(2, "log"), y = rpois(48L, exp(X[, 1]))),
    list(family = gammaScale(.5, "log"), y = rgamma(48L, 4, scale = exp(X[, 1])/4)),
    list(family = gammaCensored(.1, "log"), y = {
      value <- rgamma(48L, 4, scale = exp(X[, 1])/4)
      ifelse(value < .5, 0, value)
    }),
    list(family = binomial("probit"), y = rbinom(48L, 1, pnorm(2 * X[, 1] - 1)))
  )
  for (case in cases) {
    out <- do.call(l0glmTune, c(list(X = X, tune = "cv", maxit = c(1L, 5L),
      foldid = rep(1:3, 16), nlambda = 2, prop_range = .1,
      relaxed = FALSE, local.search = FALSE, parallelize = FALSE,
      printbest = FALSE, plot = FALSE), case))
    expect_true(all(is.finite(out$cv$deviance)))
    expect_equal(dim(out$cv$deviance), c(4L, 3L))
    expect_true(out$maxit_best %in% c(1L, 5L))
    expect_equal(coef(out$fit_1se), out$coefficients[[which(out$lambdas == out$lambda_1se &
      out$cv$maxit == out$maxit_1se)]])
  }
})

test_that("NB, multinomial and Cox iteration tuning score the corresponding training fits", {
  skip_if_not_installed("mgcv")
  skip_if_not_installed("survival")
  set.seed(7093)
  n <- 60L
  X <- matrix(rnorm(n * 4L), n)
  cases <- list(
    list(family = mgcv::nb(), y = rnbinom(n, mu = exp(.3 + .5 * X[, 1]), size = 3)),
    list(family = mgcv::multinom(K = 2), y = factor(sample(1:3, n, TRUE))),
    list(family = mgcv::cox.ph(),
      y = survival::Surv(rexp(n, exp(.5 * X[, 1])), rep(c(1, 1, 0, 1, 0), n/5)))
  )
  folds <- rep(1:3, length.out = n)
  for (case in cases) {
    out <- do.call(l0glmTune, c(list(X = X, tune = "cv", foldid = folds,
      maxit = c(2L, 6L), nlambda = 2, prop_range = .1, relaxed = FALSE,
      local.search = FALSE, printbest = FALSE, plot = FALSE), case))
    expect_true(all(is.finite(out$cv$deviance)))
    expect_equal(dim(out$cv$deviance), c(4L, 3L))
    train <- folds != 1L
    for (i in c(1L, 3L)) {
      a <- out$fit_args
      a$maxit <- out$cv$maxit[i]
      a$weights <- a$weights[train, , drop = FALSE]
      a$offset <- a$offset[train, , drop = FALSE]
      yy <- if (is.null(dim(case$y))) case$y[train] else case$y[train, , drop = FALSE]
      fit <- do.call(l0glmFit, c(list(X = X[train, ], y = yy,
        lambda = out$cv$lambda_train[i, 1]), a))
      if (out$kind == "cox") {
        loss <- -2 * (l0glmMetrics(fit, X, case$y)$loglik - as.numeric(logLik(fit))) /
          sum(case$y[!train, 2L])
      } else {
        yy <- if (is.null(dim(case$y))) case$y[!train] else case$y[!train, , drop = FALSE]
        loss <- l0glmMetrics(fit, X[!train, ], yy)$mean_deviance
        if (out$kind == "nb") expect_equal(out$cv$theta[[1L]][[i]], fit$theta)
      }
      expect_equal(unname(out$cv$deviance[i, 1L]), loss, tolerance = 1e-9)
    }
    replay <- do.call(l0glmFit, c(list(X = X, y = case$y, lambda = out$lambda_best), out$fit_args))
    expect_equal(coef(replay), coef(out$fit_best), tolerance = 1e-10)
    if (out$kind == "multinomial") expect_equal(rowSums(predict(out$fit_best)), rep(1, n))
  }
})

test_that("joint CV retains train-only imputation and selected-parent interaction recipes", {
  set.seed(7094)
  n <- 45L
  X <- matrix(rnorm(n * 9L), n)
  y <- 3 * X[, 1] - 2 * X[, 2] + 2 * X[, 1] * X[, 2] + rnorm(n)
  folds <- rep(1:3, 15)
  out <- l0glmTune(X, y, family = gaussian(), tune = "cv", foldid = folds,
    maxit = c(1L, 4L), nlambda = 2, prop_range = .1, interactions = 1L,
    relaxed = FALSE, local.search = FALSE, parallelize = FALSE, screen = "none",
    plot = FALSE, printbest = FALSE)
  expect_true(all(is.finite(out$cv$mean)))
  expect_length(out$coefficients, 4L)
  expect_equal(out$metrics$maxit, out$cv$maxit)
  a <- out$fit_args; a$maxit <- out$maxit_1se
  replay <- do.call(l0glmFit, c(list(X = X, y = y, lambda = out$lambda_1se), a))
  expect_equal(predict(replay, X), predict(out$fit_1se, X), tolerance = 1e-10)
  X[2L, 3L] <- X[7L, 2L] <- NA_real_
  source <- l0glmData(X, missing = "mean", control = list(max.cells = 200, block.columns = 3L))
  out <- l0glmTune(source, y, family = gaussian(), tune = "cv", foldid = folds,
    maxit = c(1L, 4L), nlambda = 2, prop_range = .1, screen = "SIS", nsis = 3L,
    screen.nthreads = 1L, relaxed = FALSE, local.search = FALSE, intercept = TRUE,
    parallelize = FALSE, plot = FALSE, printbest = FALSE)
  for (k in 1:3) {
    train <- folds != k
    means <- colMeans(X[train, ], na.rm = TRUE)
    impute <- function(A) {
      na <- which(is.na(A), arr.ind = TRUE)
      A[na] <- means[na[, 2L]]
      A
    }
    for (i in seq_along(out$lambdas)) {
      fit <- l0glmFit(impute(X[train, ]), y[train], family = gaussian(),
        lambda = out$cv$lambda_train[i, k], maxit = out$cv$maxit[i],
        screen = "SIS", nsis = 3L, screen.nthreads = 1L, relaxed = FALSE,
        local.search = FALSE, intercept = TRUE, parallelize = FALSE)
      loss <- mean((y[!train] - predict(fit, impute(X[!train, ])))^2)
      expect_equal(unname(out$cv$deviance[i, k]), loss, tolerance = 1e-9)
    }
  }
})

test_that("joint CV relearns channel partitions and supports predictor groups and pilot starts", {
  set.seed(7096)
  n <- 42L
  X <- matrix(runif(n * 6L), n)
  mu <- 2 + outer(X[, 1L], c(15, 12, 2, 1)) + outer(X[, 3L], c(1, 2, 12, 15))
  Y <- matrix(rpois(length(mu), mu), n)
  folds <- rep(1:3, length.out = n)
  out <- l0glmTune(X, Y, family = poisson(identity), offset = 2, nonnegative = TRUE,
    tune = "cv", foldid = folds, maxit = c(1L, 4L), precluster = TRUE,
    control.precluster = list(k = 2L), nlambda = 2L, prop_range = .1,
    relaxed = FALSE, local.search = FALSE, parallelize = FALSE, plot = FALSE, printbest = FALSE)
  expect_identical(dim(coef(out$fit_best)), c(6L, 4L))
  for (k in 1:3) {
    labels <- glmCluster(Y[folds != k, ], family = poisson(), k = 2L,
      linkage = "ward", nthreads = 1L)$clusterids
    expect_identical(out$cv$preclustering[[k]]$clusterids, labels)
  }
  y <- 2 * X[, 1] - X[, 3] + rnorm(n)
  out <- l0glmTune(X, y, family = gaussian(), tune = "cv", foldid = folds,
    maxit = c(1L, 4L), group = c(1, 1, 2, 2, 3, 3), start = "pilot",
    nlambda = 2L, prop_range = .1, relaxed = FALSE, local.search = FALSE,
    parallelize = FALSE, plot = FALSE, printbest = FALSE)
  expect_true(all(is.finite(out$cv$mean)))
  expect_length(out$cv$pilots, 3L)
  expect_true(all(vapply(out$cv$pilots, function(x) nrow(x) >= 1L, logical(1))))
  replay <- do.call(l0glmFit, c(list(X = X, y = y, lambda = out$lambda_best), out$fit_args))
  expect_equal(coef(replay), coef(out$fit_best))
})

test_that("iteration grids validate inputs, preserve scalar fits and plot without refitting", {
  set.seed(7095)
  X <- matrix(rnorm(30L * 3L), 30L)
  y <- 2 * X[, 1] + rnorm(30)
  args <- list(X = X, y = y, family = gaussian(), tune = "cv", foldid = rep(1:3, 10),
    nlambda = 2L, prop_range = .1, plot = FALSE, printbest = FALSE, local.search = FALSE)
  for (bad in list(0, -1, NA_real_, Inf, 2.5, numeric(), "x", c(1, NA), 2^32))
    expect_error(do.call(l0glmTune, c(args, list(maxit = bad))), "maxit")
  a <- args; a$tune <- "gic"; a$foldid <- NULL; a$maxit <- c(1, 3)
  expect_error(do.call(l0glmTune, a), "requires tune")
  a <- do.call(l0glmTune, c(args, list(maxit = c(4, 4))))
  b <- do.call(l0glmTune, c(args, list(maxit = 4L)))
  expect_identical(a, b)
  joint <- do.call(l0glmTune, c(args, list(maxit = c(4, 1))))
  file <- tempfile(fileext = ".pdf")
  grDevices::pdf(file)
  expect_silent(L0glm_plot_tuning(joint))
  grDevices::dev.off()
  expect_gt(file.info(file)$size, 1000)
  choice <- L0glm_cv_choice(list(mean = c(1, 1.05, 1.02, 1.1), se = c(.2, 0, 0, 0)),
    c(2, 4, 2, 4), c(1L, 1L, 5L, 5L))
  expect_identical(choice$best, 1L)
  expect_identical(choice$one_se, 4L)
})

test_that("shared ISIS pools use the largest iteration limit within each training fold", {
  set.seed(7099)
  n <- 36L
  X <- matrix(rnorm(n * 400L), n)
  y <- 3 * X[, 2] - 2 * X[, 9] + rnorm(n)
  folds <- rep(1:3, length.out = n)
  args <- list(family = gaussian(), dispersion = 1, intercept = TRUE,
    screen = "ISIS", nsis = 5L, screen.maxit = 2L, screen.nthreads = 1L,
    relaxed = FALSE, local.search = FALSE, parallelize = FALSE)
  out <- do.call(l0glmTune, c(list(X = X, y = y, tune = "cv", foldid = folds,
    maxit = c(1L, 3L), screen.path = "preset", nlambda = 3L, prop_range = .1,
    plot = FALSE, printbest = FALSE), args))
  anchor <- which(out$lambdas == out$lambda_reference & out$cv$maxit == 3L)
  for (k in 1:3) {
    train <- folds != k
    fit <- do.call(l0glmFit, c(list(X = X[train, ], y = y[train], maxit = 3L,
      lambda = out$cv$lambda_train[anchor, k]), args))
    expect_identical(out$cv$screen.candidates[[k]], fit$screening$candidates)
  }
})
