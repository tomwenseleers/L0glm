native_validation_data <- function() {
  set.seed(201026)
  X <- matrix(runif(100 * 14), 100, 14)
  V <- matrix(runif(60 * 14), 60, 14)
  beta <- matrix(0, 14, 3); beta[c(2, 7, 11), ] <- c(15, 9, 20)
  list(X = X, V = V, y = matrix(rpois(300, X %*% beta + .3), 100, 3),
    v = matrix(rpois(180, V %*% beta + .2), 60, 3))
}

test_that("native validation agrees with R for storage, groups, weights and offsets", {
  d <- native_validation_data()
  for (sparse in c(FALSE, TRUE)) for (vsparse in c(FALSE, TRUE))
    for (grouped in c(FALSE, TRUE)) for (relaxed in c(FALSE, TRUE)) {
      args <- list(X = Matrix::Matrix(d$X, sparse = sparse), y = d$y,
        weights = matrix(seq(.5, 1.5, length.out = 300), 100, 3),
        offset = .3, family = poisson("identity"), lambda = 2,
        start = 1, intercept = TRUE, normalize = "rms", nonnegative = TRUE,
        parallelize = FALSE, grouped = grouped, relaxed = relaxed,
        maxit = 20L, group.maxit = 3L,
        validation = list(X = Matrix::Matrix(d$V, sparse = vsparse), y = d$v,
          weights = matrix(seq(.6, 1.3, length.out = 180), 60, 3), offset = .2))
      # The public ordinary-GLM weight/offset format includes a shared vector.
      args$offset <- rep(.3, 100)
      r <- do.call(L0glm.fit, c(args, list(openmp = FALSE)))
      native <- do.call(L0glm.fit, c(args, list(openmp = TRUE)))
      expect_equal(coef(native), coef(r), tolerance = 1e-7)
      expect_equal(native$validation$deviance, r$validation$deviance, tolerance = 1e-9)
      expect_equal(as.matrix(native$weights), as.matrix(r$weights), tolerance = 1e-7)
      expect_equal(as.matrix(native$lambdas), as.matrix(r$lambdas), tolerance = 1e-7)
      expect_identical(inherits(native$X, "sparseMatrix"), sparse)
      expect_equal(native$validation$selection$best_iteration, r$validation$selection$best_iteration)
    }
})

test_that("native stopping restores iteration one and its working weights", {
  d <- native_validation_data(); y <- d$y[, 1]
  for (sparse in c(FALSE, TRUE)) {
    args <- list(X = Matrix::Matrix(d$X, sparse = sparse), y = y, start = 1,
      lambda = 2, relaxed = FALSE, parallelize = FALSE)
    first <- do.call(L0nnpois, c(args, list(maxit = 1L)))
    validation <- list(X = d$V, y = as.numeric(predict(first, d$V)),
      patience = 2L, minit = 3L, tol = 1e-10)
    stopped <- do.call(L0nnpois, c(args, list(maxit = 30L, validation = validation)))
    expect_equal(coef(stopped), coef(first), tolerance = 1e-10)
    expect_equal(stopped$validation$selection$best_iteration, 1L)
    expect_equal(stopped$validation$selection$iterations, 3L)
    expect_true(stopped$validation$selection$stopped)
    expect_equal(stopped$weights, 1 / as.numeric(fitted(stopped)), tolerance = 1e-10)
    expect_equal(stopped$validation$deviance, 0, tolerance = 1e-12)
  }
})

test_that("validation stays in the native backend with default initialization", {
  local_mocked_bindings(L0glm.fit.internal = function(...) stop("R IRLS fallback"))
  d <- native_validation_data()
  fit <- L0nnpois(d$X, d$y, lambda = 2, grouped = TRUE, relaxed = TRUE,
    ncores = 2L, validation = list(X = d$V, y = d$v))
  expect_true(is.finite(fit$validation$deviance))
  expect_true(any(coef(fit) > 0))
  scalar <- L0nnpois(d$X, d$y, weights = 2, offset = .3, lambda = 2,
    ncores = 2L, validation = list(X = d$V, y = d$v, weights = 2, offset = .3))
  expanded <- L0nnpois(d$X, d$y, weights = rep(2, 100), offset = rep(.3, 100), lambda = 2,
    ncores = 2L, validation = list(X = d$V, y = d$v, weights = rep(2, 60), offset = rep(.3, 60)))
  expect_equal(coef(scalar), coef(expanded), tolerance = 1e-12)
  expect_equal(scalar$validation$deviance, expanded$validation$deviance, tolerance = 1e-12)
})

test_that("OpenMP native validation preserves explicit masks and checkpoint histories", {
  d <- native_validation_data()
  start <- matrix(1, 14, 3); start[c(1, 4), ] <- 0; start[5, 2] <- 0
  for (sparse in c(FALSE, TRUE)) for (grouped in c(FALSE, TRUE)) {
    args <- list(X = Matrix::Matrix(d$X, sparse = sparse), y = d$y, start = start,
      always.include = 1L, lambda = 2, maxit = 15L, group.maxit = 3L,
      grouped = grouped, relaxed = TRUE, validation = list(X = d$V, y = d$v))
    reference <- do.call(L0nnpois, c(args, list(parallelize = FALSE)))
    for (threads in c(1L, 2L, 4L, 8L)) {
      fit <- do.call(L0nnpois, c(args, list(ncores = threads, parallelize = TRUE)))
      # Candidate timing is the only intentionally variable diagnostic.
      fit$validation$starts$seconds <- reference$validation$starts$seconds
      expect_equal(fit, reference, tolerance = 1e-12)
      expect_equal(dim(coef(fit)), c(14L, 3L))
      expect_s4_class(coef(fit), "sparseMatrix")
      expect_equal(as.numeric(coef(fit)[4, ]), numeric(3))
      expect_equal(as.numeric(coef(fit)[5, 2]), 0)
      expect_true(any(coef(fit)[1, ] != 0))
    }
  }
})

test_that("native validation handles wide systems, zero weights and mandatory null fits", {
  set.seed(201027)
  X <- matrix(runif(20 * 50), 20, 50); V <- X[1:12, ]
  beta <- c(12, rep(0, 49)); y <- rpois(20, X %*% beta + .1)
  for (sparse in c(FALSE, TRUE)) {
    args <- list(X = Matrix::Matrix(X, sparse = sparse), y = y, start = 1,
      family = poisson("identity"), lambda = "gic", maxit = 12L, relaxed = TRUE, intercept = TRUE,
      validation = list(X = V, y = rpois(12, V %*% beta + .2),
        weights = c(0, rep(1, 11)), offset = rep(.2, 12)), parallelize = FALSE)
    fit <- do.call(L0nnpois, args)
    r <- do.call(L0glm.fit, c(args, list(nonnegative = TRUE, openmp = FALSE)))
    expect_equal(coef(fit), coef(r), tolerance = 1e-7)
    expect_equal(fit$validation$deviance, r$validation$deviance, tolerance = 1e-8)
    expect_true(all(is.finite(fitted(fit))))
    args$start <- 0
    null <- do.call(L0nnpois, args)
    expect_equal(as.numeric(coef(null)[-1, ]), numeric(50))
    expect_true(coef(null)[1, ] > 0)
  }
})
