test_that("criterion presets match the documented penalties in any case", {
  for (n in c(40, 400)) {
    for (p in c(3, 30)) {
      expected <- c(aic = 2, bic = log(n), gic = log(p) * log(log(n)),
                    mebic = 2 * log(p) * log(log(n)),
                    ebic = log(n) + 2 * max(0, 1 - log(n)/(2 * log(p))) * log(p),
                    mbic = log(n * p^2 / 16))
      for (ic in names(expected)) {
        mixed <- paste0(toupper(substr(ic, 1, 1)), substring(ic, 2))
        for (name in c(ic, toupper(ic), mixed)) {
          expect_equal(L0glm:::preset.lambda(name, p = p, n = n),
                       unname(expected[ic]), info = name)
        }
      }
      scores <- L0glm:::get.IC(loglik = c(-12, -18), n = n, p = p,
                              edf = c(1, 3))
      for (ic in names(expected)) {
        expect_equal(scores[[ic]], c(24, 36) + expected[[ic]] * c(1, 3),
                     info = ic)
      }
    }
  }
  expect_equal(L0glm:::preset.lambda("EBIC", p = 20, n = 100, gamma_ebic = 1),
               log(100) + 2 * log(20))
})

test_that("both fitters resolve every named preset before fitting", {
  set.seed(71)
  X <- matrix(runif(160), 40, 4)
  y <- rpois(40, 2 + 4 * X[, 1])
  args <- list(X = X, y = y, family = poisson(identity), intercept = TRUE,
               parallelize = FALSE, grouped = FALSE, relaxed = FALSE,
               maxit = 5)

  for (fitter in list(L0glm.fit, L0nnpois)) {
    if ("nonnegative" %in% names(formals(fitter))) args$nonnegative <- TRUE
    else args$nonnegative <- NULL
    for (ic in c("aic", "bic", "gic", "mebic", "ebic", "mbic")) {
      penalty <- L0glm:::preset.lambda(ic, p = 5, n = 40)
      numeric_fit <- do.call(fitter, c(args, list(lambda = penalty)))
      mixed <- paste0(toupper(substr(ic, 1, 1)), substring(ic, 2))
      for (name in c(ic, toupper(ic), mixed)) {
        named_fit <- do.call(fitter, c(args, list(lambda = name)))
        expect_equal(named_fit$init.lambda, penalty, info = name)
        expect_identical(named_fit$coefficients, numeric_fit$coefficients,
                         info = name)
        expect_identical(named_fit$lambdas, numeric_fit$lambdas, info = name)
      }
    }
  }
})

test_that("MEBIC works for sparse grouped and relaxed fits", {
  set.seed(72)
  X <- Matrix::rsparsematrix(40, 8, density = 0.5)
  X@x <- abs(X@x)
  beta <- matrix(0, 8, 2)
  beta[c(2, 6), ] <- c(3, 2)
  y <- matrix(rpois(80, as.vector(1 + X %*% beta)), 40, 2)
  penalty <- 2 * log(8) * log(log(40))
  args <- list(X = X, y = y, family = poisson(identity),
               offset = rep(1, 40), intercept = FALSE,
               grouped = TRUE, parallelize = FALSE, group.maxit = 2,
               maxit = 8)
  for (fitter in list(L0glm.fit, L0nnpois)) {
    if ("nonnegative" %in% names(formals(fitter))) args$nonnegative <- TRUE
    else args$nonnegative <- NULL
    for (relaxed in c(FALSE, TRUE)) {
      args$relaxed <- relaxed
      numeric_fit <- do.call(fitter, c(args, list(lambda = penalty)))
      named_fit <- do.call(fitter, c(args, list(lambda = "MeBiC")))
      expect_equal(named_fit$init.lambda, penalty)
      expect_identical(named_fit$coefficients, numeric_fit$coefficients)
      expect_identical(named_fit$lambdas, numeric_fit$lambdas)
    }
  }
})

test_that("MEBIC metrics are available without known true coefficients", {
  X <- cbind(1, seq(0, 1, length.out = 20))
  beta <- c(1, 2)
  metrics <- L0glm:::get.metrics(beta, X, y = matrix(rep(1:4, 5)),
                                family = poisson(identity))
  expect_equal(metrics$mebic + 2 * metrics$loglik,
               2 * (metrics$gic + 2 * metrics$loglik))
  expect_false("hbic" %in% names(metrics))
  expect_true(is.na(metrics$MD_selcoefs))
})

test_that("tuning selects MEBIC scores with case-insensitive names", {
  set.seed(73)
  X <- matrix(rnorm(240), 60, 4)
  y <- 8 * X[, 1] + 2 * X[, 2] + rnorm(60, sd = 0.2)
  args <- list(X = X, y = y, family = gaussian(), dispersion = 1,
               nlambda = 3, prop_range = 0.2, printbest = FALSE, plot = FALSE,
               grouped = FALSE, parallelize = FALSE, maxit = 20)
  lower <- do.call(L0glm.tune, c(args, list(IC = "mebic")))
  upper <- do.call(L0glm.tune, c(args, list(IC = "MEBIC")))
  expect_identical(upper, lower)
  expect_identical(upper$IC, "mebic")
  expect_equal(upper$lambda_preset, 2 * log(4) * log(log(60)))
  expect_equal(upper$lambda_best,
               upper$metrics$lambda[which.min(upper$metrics$mebic)])
  default <- do.call(L0glm.tune, args)
  expect_identical(default$IC, "gic")
  expect_error(do.call(L0glm.tune, c(args, list(IC = "HBIC"))),
               "Information criterion must be one of")
})

test_that("invalid criteria and lambda values fail clearly in both fitters", {
  X <- cbind(1, seq(0, 1, length.out = 12))
  for (fitter in list(L0glm.fit, L0nnpois)) {
    args <- list(X = X, y = rep(2, 12), parallelize = FALSE)
    for (bad in list("unknown", "hbic", "HBIC", "", NA_character_, c("aic", "bic"))) {
      expect_error(do.call(fitter, c(args, list(lambda = bad))),
                   "Information criterion must be one of.*mebic")
    }
    for (bad in list(NA_real_, Inf, c(1, 2), numeric(0))) {
      expect_error(do.call(fitter, c(args, list(lambda = bad))),
                   "lambda must be a finite numeric scalar")
    }
  }
})
