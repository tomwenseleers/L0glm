validation_data <- function() {
  set.seed(301)
  X <- matrix(rnorm(150 * 20), 150, 20)
  V <- matrix(rnorm(70 * 20), 70, 20)
  beta <- c(1, -.7, rep(0, 18))
  list(X = X, V = V, y = as.numeric(X %*% beta) + rnorm(150),
       v = as.numeric(V %*% beta) + rnorm(70))
}

test_that("validation returns the best checkpoint and records both fitting stages", {
  d <- validation_data()
  fit <- L0glm.fit(d$X, d$y, lambda = "gic", relaxed = TRUE, maxit = 100,
    validation = list(X = d$V, y = d$v), parallelize = FALSE)
  loss <- mean((d$v - as.numeric(predict(fit, d$V)))^2)
  expect_equal(loss, fit$validation$deviance, tolerance = 1e-12)
  expect_lte(loss, mean(d$v^2))
  for (stage in c("selection", "relaxed")) {
    info <- fit$validation[[stage]]
    if (is.null(info)) next
    expect_lte(info$deviance, min(info$history$deviance) + 1e-8)
    expect_equal(info$deviance,
      info$history$deviance[info$history$iteration == info$best_iteration])
    if (stage == "selection") {
      selected <- info$history$selected[info$history$iteration > 0]
      expect_true(all(diff(selected) <= 0))
    }
  }
  expect_s4_class(coef(fit), "sparseMatrix")
  expect_equal(dim(coef(fit)), c(20L, 1L))
  # A validation target opposed to the training signal favors the null model.
  null <- L0glm.fit(d$X, d$y, lambda = "gic", relaxed = TRUE, maxit = 100,
    validation = list(X = d$V, y = numeric(70), patience = 2L, minit = 3L),
    parallelize = FALSE)
  expect_true(all(coef(null) == 0))
  expect_equal(null$validation$selection$best_iteration, 0L)
  first <- L0glm.fit(d$X, d$y, start = 1, dispersion = 1, lambda = "gic",
    maxit = 1, parallelize = FALSE)
  stopped <- L0glm.fit(d$X, d$y, start = 1, dispersion = 1, lambda = "gic", maxit = 100,
    validation = list(X = d$V, y = as.numeric(predict(first, d$V)),
                      patience = 2L, minit = 3L), parallelize = FALSE)
  expect_equal(coef(stopped), coef(first), tolerance = 1e-12)
  expect_equal(stopped$validation$selection$best_iteration, 1L)
  expect_true(stopped$validation$selection$stopped)
  expect_equal(stopped$validation$selection$iterations, 3L)
})

test_that("weighted validation respects offsets, predictor units and column names", {
  d <- validation_data()
  colnames(d$X) <- colnames(d$V) <- paste0("x", 1:20)
  w <- seq(.4, 2, length.out = 150); wv <- seq(.2, 1, length.out = 70)
  off <- seq(-1, 1, length.out = 150); ov <- seq(-.5, .5, length.out = 70)
  yy <- 2 + d$y + off; vv <- 2 + d$v + ov
  run <- function(X, V) L0glm.fit(X, yy, weights = w, offset = off,
    intercept = TRUE, lambda = "gic", maxit = 100, relaxed = TRUE,
    validation = list(X = V, y = vv, weights = wv, offset = ov), parallelize = FALSE)
  a <- run(d$X, d$V[, 20:1])
  expected <- weighted.mean((vv - as.numeric(predict(a, d$V, offset = ov)))^2, wv)
  expect_equal(expected, a$validation$deviance, tolerance = 1e-10)
  scale <- seq(.1, 5, length.out = 20); shift <- seq(-2, 2, length.out = 20)
  XT <- sweep(sweep(d$X, 2, scale, "*"), 2, shift, "+")
  VT <- sweep(sweep(d$V, 2, scale, "*"), 2, shift, "+")
  b <- run(XT, VT)
  expect_equal(as.numeric(predict(a, d$V, offset = ov)),
               as.numeric(predict(b, VT, offset = ov)), tolerance = 1e-7)
  context <- L0glm_validation_context(list(X = d$V, y = d$v,
    weights = c(0, rep(1, 69))), d$X, d$y, gaussian(), FALSE)
  context$y[1] <- 1e300
  expect_equal(L0glm_validation_loss(context, numeric(20)), mean(d$v[-1]^2))
})

test_that("checkpointing preserves start masks, constraints and fixed dispersion", {
  d <- validation_data()
  start <- c(0, 1, rep(0, 19))
  fit <- L0glm.fit(d$X, d$y, intercept = TRUE, start = start, dispersion = 2,
    lambda = "gic", relaxed = TRUE, maxit = 100,
    validation = list(X = d$V, y = d$v), parallelize = FALSE)
  expect_equal(dim(coef(fit)), c(21L, 1L))
  expect_true(all(coef(fit)[3:21, 1] == 0))
  expect_false(any(fit$validation$starts$start == "marginal"))
  expect_identical(fit$dispersion, 2)
  scaled <- L0glm.fit(d$X, 10 * d$y, intercept = TRUE, start = start * 10,
    dispersion = 200, lambda = "gic", relaxed = TRUE, maxit = 100,
    validation = list(X = d$V, y = 10 * d$v), parallelize = FALSE)
  expect_equal(coef(fit), coef(scaled) / 10, tolerance = 1e-7)
  positive <- L0glm.fit(abs(d$X), abs(d$y), intercept = TRUE, nonnegative = TRUE,
    always.include = 3, lambda = "gic", relaxed = TRUE,
    validation = list(X = abs(d$V), y = abs(d$v)), parallelize = FALSE)
  expect_true(all(coef(positive) >= 0))
  expect_equal(positive$always.include, c(1L, 4L))
})

test_that("validation selects batches and starts by minimum error or the one-SE rule", {
  d <- validation_data()
  run <- function(rule) L0glm.fit(d$X, d$y, screen = "isis", nsis = c(1:5, 10),
    screen.maxit = 4, lambda = "gic", relaxed = TRUE, maxit = 100,
    validation = list(X = d$V, y = d$v, rule = rule), parallelize = FALSE)
  for (rule in c("min", "1se")) {
    fit <- run(rule)
    table <- fit$validation$starts
    expect_true(all(c("default", "marginal", "null") %in% table$start))
    expect_equal(nrow(table), 13L)
    expect_equal(fit$validation$minimum, min(table$validation_deviance))
    limit <- fit$validation$minimum + if (rule == "1se") fit$validation$se else 0
    expect_lte(fit$validation$deviance, limit + 1e-12)
    eligible <- table[table$validation_deviance <= limit, ]
    expect_equal(sum(coef(fit) != 0), min(eligible$selected))
    expect_true(is.na(fit$validation$selected_nsis) || fit$validation$selected_nsis %in% c(1:5, 10))
    expect_equal(dim(coef(fit)), c(20L, 1L))
    expect_equal(as.numeric(predict(fit, d$V)), as.numeric(d$V %*% coef(fit)), tolerance = 1e-12)
  }
  fixed <- L0glm.fit(d$X, d$y, screen = "ISIS", lambda = "gic", relaxed = TRUE,
    validation = list(X = d$V, y = d$v, starts = "default"), parallelize = FALSE)
  expected <- min(ncol(d$X), max(1, round(nrow(d$X) /
    (log(ncol(d$X)) * log(log(nrow(d$X)))))))
  expect_equal(fixed$validation$starts$nsis[fixed$validation$starts$start == "default"], expected)
  expect_equal(fixed$init.lambda, preset.lambda("gic", ncol(d$X), nrow(d$X)))
  # With one equally weighted validation observation the empirical SE is zero.
  one <- L0glm.fit(d$X, d$y, lambda = "gic", validation = list(
    X = d$V[1, , drop = FALSE], y = d$v[1], rule = "1se"), parallelize = FALSE)
  expect_equal(one$validation$se, 0)
})

test_that("validation batch grids and selection controls validate clearly", {
  d <- validation_data()
  validation <- list(X = d$V, y = d$v)
  for (bad in list(numeric(), 0, 1.5, NA_real_, c(1, Inf), "2"))
    expect_error(L0glm.fit(d$X, d$y, screen = "ISIS", nsis = bad,
      validation = validation), "positive integers")
  expect_error(L0glm.fit(d$X, d$y, nsis = 1:3, validation = validation), "batch grid requires")
  expect_error(L0glm.fit(d$X, d$y,
    validation = c(validation, list(starts = "supplied"))), "Validation starts")
  expect_error(L0glm.fit(d$X, d$y,
    validation = c(validation, list(rule = "unknown"))), "Validation rule")
})

test_that("a one-predictor marginal start retains original coefficient units", {
  X <- matrix(seq(-2, 2, length.out = 40), ncol = 1)
  start <- L0glm_marginal_start(X, X[, 1], NULL, NULL, FALSE, NULL)
  expect_equal(dim(start), c(1L, 1L))
  expect_equal(as.numeric(start), 1)
  fit <- L0glm.fit(X, X[, 1], start = start, lambda = 0, dispersion = 1,
    normalize = "rms", maxit = 1, parallelize = FALSE)
  expect_equal(as.numeric(coef(fit)), 1, tolerance = 1e-12)
})

test_that("screening uses validation checkpoints with original-dimension coefficients", {
  d <- validation_data()
  for (mode in c("SIS", "ISIS")) {
    fit <- L0glm.fit(Matrix::Matrix(d$X, sparse = TRUE), d$y, lambda = "gic",
      screen = mode, nsis = 3L, screen.maxit = 4L, relaxed = TRUE, maxit = 100,
      validation = list(X = Matrix::Matrix(d$V, sparse = TRUE), y = d$v), parallelize = FALSE)
    expect_equal(dim(coef(fit)), c(20L, 1L))
    expect_equal(ncol(fit$X), 20L)
    expect_equal(mean((d$v - as.numeric(predict(fit, d$V)))^2),
                 min(fit$screening$history$validation_deviance), tolerance = 1e-10)
    expect_equal(fit$validation$deviance, min(fit$screening$history$validation_deviance))
  }
})

test_that("validation rejects malformed inputs and unsupported families", {
  d <- validation_data()
  run <- function(validation, ...) L0glm.fit(d$X, d$y, validation = validation, ...)
  expect_error(run(list(y = d$v)), "contain X and y")
  expect_error(run(list(X = d$V[, 1:3], y = d$v)), "dimensions")
  expect_error(run(list(X = d$V, y = letters[1:70])), "numeric dimensions")
  expect_error(run(list(X = d$V, y = d$v, weights = c(1, 2))), "incompatible dimensions")
  expect_error(run(list(X = d$V, y = d$v, weights = 0)), "Invalid validation")
  expect_error(run(list(X = d$V, y = d$v, patience = 0)), "positive integers")
  expect_error(run(list(X = d$V, y = d$v, tol = "a")), "tol must")
  expect_error(run(list(X = d$V, y = d$v, unknown = 1)), "named validation")
  expect_error(L0glm.fit(d$X, abs(d$y), family = poisson(identity),
    validation = list(X = d$V, y = d$v)), "nonnegative")
  expect_error(L0glm.fit(d$X, cbind(d$y, d$y),
    validation = list(X = d$V, y = d$v)), "dimensions")
})
