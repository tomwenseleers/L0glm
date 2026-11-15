test_that("batched bounded solves match independent support enumeration", {
  for (seed in 1:6) {
    set.seed(7200 + seed)
    p <- 6L
    X <- matrix(rnorm(80 * p), 80, p)
    X[, 2] <- X[, 1] + .03 * X[, 2]
    y <- as.numeric(X %*% rnorm(p) * 20) + rnorm(80)
    d <- runif(p, .01, 2); center <- rnorm(p)
    G <- crossprod(X) + diag(d)
    h <- as.numeric(crossprod(X, y)) + d * center
    best <- numeric(p); value <- Inf
    objective <- function(b) sum((y - X %*% b)^2) + sum(d * (b - center)^2)
    for (mask in 0:(2^p - 1L)) {
      free <- which(as.logical(intToBits(mask)[seq_len(p)]))
      b <- numeric(p)
      if (length(free)) b[free] <- solve(G[free, free, drop = FALSE], h[free])
      gradient <- as.numeric(G %*% b) - h
      if (any(b < -1e-9) || any(gradient[setdiff(seq_len(p), free)] < -1e-8)) next
      if (objective(b) < value) { best <- b; value <- objective(b) }
    }
    expect_true(is.finite(value))
    start <- pmax(0, as.numeric(solve(G, h)))
    for (sparse in c(FALSE, TRUE)) {
      A <- if (sparse) as(as(Matrix::Matrix(G, sparse = TRUE), "generalMatrix"), "CsparseMatrix") else G
      fit <- solve_bounded_gram(A, h, start, rep(0, p), rep(Inf, p), batched = TRUE)
      expect_true(fit$success)
      expect_lte(fit$error, 1e-10)
      expect_equal(as.numeric(fit$coefficients), best, tolerance = 1e-8)
      expect_lte(objective(fit$coefficients), objective(start) + 1e-8)
    }
  }
  empty <- solve_bounded_gram(matrix(0, 0, 0), numeric(), numeric(), numeric(), numeric(), TRUE)
  expect_true(empty$success)
  expect_error(solve_bounded_gram(diag(2), 1, 1, 0, Inf, TRUE), "dimensions")
  # Fixed and unidentified coordinates remain feasible; singular free faces fail safely.
  zero <- solve_bounded_gram(diag(c(2, 0)), c(4, 0), c(1, 3), c(0, 0), c(Inf, Inf), TRUE)
  expect_true(zero$success)
  expect_equal(as.numeric(zero$coefficients), c(2, 3))
  singular <- solve_bounded_gram(matrix(1, 2, 2), c(1, 1), c(1, 1), c(0, 0), c(Inf, Inf), TRUE)
  expect_false(singular$success)
})

test_that("all existing Gram solve paths can polish the same constrained ridge", {
  set.seed(16021)
  X <- matrix(rnorm(80 * 8), 80, 8)
  X[, 2] <- X[, 1] + .05 * X[, 2]
  y <- as.numeric(X %*% c(3, -2, 2, -1, 0, 0, 1, 0)) + rnorm(80)
  w <- runif(80, .2, 2); d <- exp(seq(-2, 2, length.out = 8)); center <- seq(-.3, .4, .1)
  G <- crossprod(X, w * X) + diag(d)
  h <- as.numeric(crossprod(X, w * y)) + d * center
  lower <- c(-Inf, rep(0, 7)); upper <- c(Inf, rep(2, 7))
  reference <- solve_bounded_gram(G, h, pmax(lower, pmin(upper, solve(G, h))), lower, upper)$coefficients
  for (sparse in c(FALSE, TRUE)) for (method in 1:6) {
    A <- if (sparse) Matrix::Matrix(X, sparse = TRUE) else X
    control <- controlSolverGen("eigen", method, rowaug = FALSE, nonnegative.method="constrained")
    fit <- lmRidge(A, y, weights = w, lambdas = d, lower = lower, upper = upper,
      prior.mean = center, method = method, control.solver = control)
    expect_equal(fit, as.numeric(reference), tolerance = 1e-7, info = paste(sparse, method))
    projected <- as.numeric(G %*% fit) - h
    projected[(fit <= lower & projected > 0) | (fit >= upper & projected < 0)] <- 0
    expect_lt(max(abs(projected)), 1e-7)
  }
})

test_that("polishing reuses the cached weighted Gram without accessing raw data", {
  set.seed(16022)
  X <- matrix(rnorm(100 * 8), 100, 8); y <- rnorm(100); w <- runif(100, .1, 3)
  control <- controlSolverGen("eigen", 2L, nonnegative.method="constrained")
  expected <- lmRidge(X, y, weights = w, lambdas = rep(.3, 8), nonnegative = TRUE,
    method = 2L, control.solver = control)
  control$.crossprod <- list(xx = crossprod(X, w * X), xy = as.numeric(crossprod(X, w * y)))
  actual <- lmridge.eigen(matrix(NA_real_, 100, 8), rep(NA_real_, 100), lambdas = rep(.3, 8),
    lower = rep(0, 8), control.solver = control)
  expect_equal(actual, expected, tolerance = 1e-9)
})

test_that("default, signed, rectangular and Woodbury paths bypass the polish", {
  set.seed(16023)
  X <- matrix(rnorm(40 * 12), 40, 12); y <- rnorm(40)
  solve_original <- solve_bounded_gram
  calls <- 0L
  local_mocked_bindings(solve_bounded_gram = function(...) {
    calls <<- calls + 1L
    solve_original(...)
  })
  for (sparse in c(FALSE, TRUE)) for (method in 1:6) {
    A <- if (sparse) Matrix::Matrix(X, sparse = TRUE) else X
    ordinary <- lmRidge(A, y, lambdas = rep(1, 12), nonnegative = TRUE, method = method)
    explicit <- lmRidge(A, y, lambdas = rep(1, 12), nonnegative = TRUE, method = method,
      control.solver = controlSolverGen("eigen", method, nonnegative.method="clamp"))
    expect_identical(explicit, ordinary)
    signed <- lmRidge(A, y, lambdas = rep(1, 12), method = method)
    enabled <- lmRidge(A, y, lambdas = rep(1, 12), method = method,
      control.solver = controlSolverGen("eigen", method, nonnegative.method="constrained"))
    expect_identical(enabled, signed)
    wide <- t(A)
    off <- lmRidge(wide, y[1:12], lambdas = rep(1, 40), nonnegative = TRUE, method = method)
    on <- lmRidge(wide, y[1:12], lambdas = rep(1, 40), nonnegative = TRUE, method = method,
      control.solver = controlSolverGen("eigen", method, nonnegative.method="constrained"))
    expect_identical(on, off)
    if (method %in% c(4L, 5L, 6L)) {
      control <- controlSolverGen("eigen", method, rowaug = TRUE, nonnegative.method="constrained")
      on <- lmRidge(A, y, lambdas = rep(1, 12), nonnegative = TRUE, method = method, control.solver = control)
      control$nonnegative.method <- "clamp"
      off <- lmRidge(A, y, lambdas = rep(1, 12), nonnegative = TRUE, method = method, control.solver = control)
      expect_identical(on, off)
    }
  }
  expect_identical(calls, 0L)
})

test_that("native polish respects channel masks, mandatory terms and threading", {
  set.seed(612)
  X <- matrix(rexp(120 * 12), 120, 12)
  Y <- sapply(1:4, function(j) rpois(120, 3 + X[, 1] * 2 + X[, 5] + j / 3))
  start <- matrix(1, 12, 4); start[2, ] <- 0; start[6, 1] <- 0; start[12, ] <- 0
  prior <- matrix(.1, 12, 4); prior[2, ] <- 2
  for (sparse in c(FALSE, TRUE)) for (method in c(1L, 2L, 3L, 4L)) {
    A <- if (sparse) Matrix::Matrix(X, sparse = TRUE) else X
    fits <- lapply(c(1L, 4L), function(threads) l0glmFit(A, Y, family = poisson(identity),
      start = start, prior.mean = prior, always.include = 12L, nonnegative = TRUE,
      weights = seq(.5, 2, length.out = 120), offset = rep(3, 120),
      local.search = FALSE, ncores = threads, nthreads = 1L,
      method = method, control.solver = controlSolverGen("eigen", method,
        rowaug = FALSE, nonnegative.method="constrained")))
    for (fit in fits) {
      b <- as.matrix(coef(fit))
      expect_equal(dim(b), c(12L, 4L))
      expect_true(all(is.finite(b) & b >= 0))
      expect_true(all(b[2, ] == 0))
      expect_equal(b[6, 1], 0)
      expect_true(all(as.matrix(fit$lambdas)[12, ] <= 1e-20))
    }
    expect_equal(as.matrix(coef(fits[[1]])), as.matrix(coef(fits[[2]])), tolerance = 1e-10)
  }
})

test_that("LLT and LDLT polish recover the correlated Poisson peak pair", {
  d <- simul(method = 2L, n = 1000L, p = 1000L, sparse = TRUE, seed = 1L,
    options = list(sd_noise = 1))
  for (native in c(FALSE, TRUE)) for (method in c(1L, 2L)) {
    fit <- l0glmFit(d$X, d$y, family = poisson(identity), nonnegative = TRUE,
      start = NULL, lambda = "gic", local.search = FALSE,ncores=1L,
      method = method,control.solver=utils::modifyList(controlSolverGen("eigen", method, nonnegative.method="constrained"),list(.native=native)))
    b <- as.numeric(coef(fit))
    expect_true(all(b[c(836, 838)] > 0))
    expect_equal(b[c(837, 839)], c(0, 0))
    expect_gt(l0glmMetrics(fit, beta_true = d$beta_true)$R2_coefs, .95)
  }
  direct <- l0nnPois(d$X, d$y, start = NULL, lambda = "gic", local.search = FALSE,ncores=1L,
    control.solver = controlSolverGen("eigen", 1L, nonnegative.method="constrained"))
  expect_gt(l0glmMetrics(direct, beta_true = d$beta_true)$R2_coefs, .95)
})

test_that("native Gaussian caches and predictor groups accept the polish", {
  set.seed(16024)
  X <- matrix(rnorm(100 * 12), 100, 12)
  X[, 2] <- X[, 1] + .03 * X[, 2]
  Y <- cbind(2 + 3 * X[, 1] + X[, 5], 1 + 2 * X[, 1] + 2 * X[, 5]) + rnorm(200)
  for (group in list(NULL, rep(1:6, each = 2L))) {
    fits <- lapply(c(FALSE, TRUE), function(native) l0glmFit(X, Y, family = gaussian(),
      nonnegative = TRUE, intercept = TRUE, group = group, local.search = FALSE,ncores=1L,control.solver=utils::modifyList(controlSolverGen("eigen", 1L, nonnegative.method="constrained"),list(.native=native))))
    for (fit in fits) expect_true(all(is.finite(coef(fit)) & coef(fit) >= 0))
    expect_equal(as.matrix(coef(fits[[1]])), as.matrix(coef(fits[[2]])), tolerance = 2e-5)
  }
})

test_that("a binomial intercept stays unrestricted with nonnegative polishing", {
  set.seed(16025)
  X <- matrix(rexp(180 * 8), 180, 8)
  y <- rbinom(180, 1, plogis(-4 + X[, 1] + .7 * X[, 3]))
  for (native in c(FALSE, TRUE)) for (method in c(1L, 2L)) {
    fit <- l0glmFit(X, y, family = binomial(), intercept = TRUE, nonnegative = TRUE,
      local.search = FALSE,ncores=1L, method = method,control.solver=utils::modifyList(controlSolverGen("eigen", method, nonnegative.method="constrained"),list(.native=native)))
    b <- as.numeric(coef(fit))
    expect_lt(b[1], 0)
    expect_true(all(b[-1] >= 0))
  }
})

test_that("tuning forwards the constrained-polish control to candidate fits", {
  set.seed(16026)
  X <- matrix(rexp(80 * 8), 80, 8)
  y <- rpois(80, 2 * X[, 1] + X[, 3])
  control <- controlSolverGen("eigen", 2L, nonnegative.method="constrained")
  tuned <- l0glmTune(X, y, family = poisson(identity), nonnegative = TRUE,
    tune = "gic", nlambda = 3L, plot = FALSE,ncores=1L,
    method = 2L, control.solver = control, local.search = FALSE)
  reference <- l0glmFit(X, y, family = poisson(identity), nonnegative = TRUE,
    lambda = tuned$lambda_best,ncores=1L,
    method = 2L, control.solver = control, local.search = FALSE)
  expect_equal(as.numeric(coef(tuned$fit_best)), as.numeric(coef(reference)), tolerance = 1e-8)
})
