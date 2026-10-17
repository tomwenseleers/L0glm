test_that("wide ridge profiling agrees with the full normal equations", {
  set.seed(1660)
  X <- cbind(1, seq(-1, 1, length.out = 18), matrix(rnorm(18 * 60), 18, 60) + 2)
  y <- rnorm(18); w <- runif(18, 0.2, 2)
  d <- c(0, 0, exp(seq(-1, 2, length.out = 60)))
  center <- rnorm(ncol(X)) / 10
  G <- crossprod(X, w * X) + diag(d)
  rhs <- crossprod(X, w * y) + d * center
  reference <- as.numeric(solve(G, rhs))
  for (sparse in c(FALSE, TRUE)) for (method in 1:6) {
    beta <- lmridge(Matrix::Matrix(X, sparse = sparse), y, weights = w,
                    lambdas = d, prior.mean = center, method = method,
                    control.solver = control.solver.gen("eigen", method, tol = 1e-12))
    expect_equal(beta, reference, tolerance = 1e-7)
    expect_lt(max(abs(G %*% beta - rhs)), 1e-7)
  }
})

test_that("profiling preserves dense or sparse observation-space systems", {
  X <- cbind(1, matrix(seq_len(8 * 31) %% 17, 8, 31))
  original <- eigen_ridge_solve
  seen <- list()
  local_mocked_bindings(eigen_ridge_solve = function(A, b, start, control, iteration) {
    seen[[length(seen) + 1L]] <<- list(dim = dim(A), sparse = issparse(A), rhs = dim(b))
    original(A, b, start, control, iteration)
  })
  for (sparse in c(FALSE, TRUE)) {
    lmridge(Matrix::Matrix(X, sparse = sparse), 1:8, lambdas = c(0, rep(1, 31)))
    call <- seen[[length(seen)]]
    expect_equal(call$dim, c(8L, 8L))
    expect_identical(call$sparse, sparse)
    expect_equal(call$rhs, c(8L, 2L))
  }
})

test_that("native profiling handles wide models with several mandatory terms", {
  set.seed(1661)
  n <- 18L; p <- 60L
  X <- cbind(Intercept = 1, matrix(runif(n * (p - 1)), n, p - 1))
  y <- rpois(n, 3); w <- 1 / (y + 0.1)
  start <- rep(0.3, p); start[1:2] <- 0
  d <- rep(2 / 4 / 0.3^2, p); d[1:2] <- 0
  G <- crossprod(X, w * X) + diag(d)
  reference <- pmax(0, as.numeric(solve(G, crossprod(X, w * y))))
  for (sparse in c(FALSE, TRUE)) for (natural in c(FALSE, TRUE)) {
    fit <- L0nnpois(Matrix::Matrix(X, sparse = sparse), y,
                    always.include = 2L, lambda = 2, start = start,
                    maxit = 1L, relaxed = FALSE, thresh = 0, parallelize = FALSE,
                    control.solver = control.solver.gen("eigen", 1L,
                                                        NaturalOrdering = natural))
    expect_equal(as.numeric(fit$coefficients), reference, tolerance = 1e-9)
    expect_equal(fit$always.include, 1:2)
  }
})

test_that("GIC counts the original selectable pool despite intercepts and masks", {
  set.seed(1662)
  X <- matrix(runif(40 * 7), 40, 7)
  y <- rpois(40, 3)
  for (fitter in list(L0glm.fit, L0nnpois)) for (intercept in c(FALSE, TRUE)) {
    args <- list(X = X, y = y, family = poisson(identity), lambda = "gic",
                 intercept = intercept, always.include = 2L,
                 start = c(1, 0, rep(0, 5)), maxit = 1L, parallelize = FALSE)
    if (intercept) args$start <- c(0, args$start)
    if (identical(fitter, L0glm.fit)) args$nonnegative <- TRUE
    fit <- do.call(fitter, args)
    expect_equal(fit$p, 6L)
    expect_equal(fit$init.lambda, preset.lambda("gic", 6, 40))
  }
  XX <- cbind(Intercept = 1, X)
  met <- get.metrics(matrix(0, 8, 1), XX, matrix(y), family = poisson(identity),
                      always.include = 3L)
  expect_equal(met$p, 6L)
  expect_equal(met$support_size_estimated, 0)
  expect_equal(met$gic, -2 * met$loglik + preset.lambda("gic", 6, 40) * 2)
})
