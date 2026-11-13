test_that("all Eigen methods solve tall and wide weighted ridge objectives", {
  set.seed(1601)
  for (shape in list(c(45L, 12L), c(12L, 45L))) {
    X <- matrix(rnorm(prod(shape)), shape[1L], shape[2L])
    y <- rnorm(nrow(X))
    w <- runif(nrow(X), 0.2, 2)
    d <- exp(seq(-2, 2, length.out = ncol(X)))
    center <- rnorm(ncol(X)) / 10
    G <- crossprod(X, w * X) + diag(d)
    rhs <- crossprod(X, w * y) + d * center
    reference <- as.numeric(solve(G, rhs))
    for (sparse in c(FALSE, TRUE)) for (method in 1:6) {
      for (rowaug in c(FALSE, TRUE)) {
        control <- controlSolverGen("eigen", method, rowaug = rowaug,
                                      tol = 1e-11, maxit = 3000L)
        beta <- lmRidge(Matrix::Matrix(X, sparse = sparse), y, weights = w,
                        lambdas = d, prior.mean = center, method = method,
                        control.solver = control)
        expect_equal(beta, reference, tolerance = 1e-6,
                     info = paste(shape, sparse, method, rowaug))
        expect_lt(max(abs(G %*% beta - rhs)) / max(abs(rhs)), 1e-6)
      }
    }
  }
})

test_that("wide systems use the requested kernel in observation coordinates", {
  set.seed(1602)
  X <- matrix(rnorm(8 * 31), 8, 31)
  original <- eigen_ridge_solve
  calls <- list()
  local_mocked_bindings(eigen_ridge_solve = function(A, b, start, control, iteration, ...) {
    calls[[length(calls) + 1L]] <<- list(dim = dim(A), sparse = issparse(A),
                                       method = control$method, start = length(start))
    original(A, b, start, control, iteration, ...)
  })
  for (sparse in c(FALSE, TRUE)) for (method in 1:6) {
    beta <- lmRidge(Matrix::Matrix(X, sparse = sparse), rnorm(8),
                    lambdas = rep(2, 31), method = method)
    call <- calls[[length(calls)]]
    expect_equal(call$dim, c(8L, 8L))
    expect_identical(call$sparse, sparse)
    expect_identical(call$method, controlSolverGen("eigen", method)$method)
    expect_equal(call$start, 8L)
    expect_length(beta, 31L)
  }
})

test_that("unpenalized wide solves return the minimum norm solution", {
  set.seed(1603)
  X <- matrix(rnorm(8 * 31), 8, 31)
  y <- rnorm(8)
  reference <- as.numeric(t(X) %*% solve(tcrossprod(X), y))
  for (sparse in c(FALSE, TRUE)) for (method in 1:6) {
    beta <- lmRidge(Matrix::Matrix(X, sparse = sparse), y, method = method)
    expect_equal(beta, reference, tolerance = 1e-6)
  }
})

test_that("dense and sparse CG honor the selected preconditioner", {
  set.seed(1605)
  X <- matrix(rnorm(80 * 25), 80, 25)
  A <- crossprod(X) + diag(25)
  b <- matrix(rnorm(25))
  sparse <- as(as(Matrix::Matrix(A, sparse = TRUE), "generalMatrix"), "CsparseMatrix")
  original <- A + 0
  for (precond in 0:2) for (limit in c(1L, 100L)) {
    dense <- solve_dense_CG(A, b, matrix(0, 25), tol = 1e-12,
                            maxit = limit, precond = precond)
    sparse_fit <- solve_sparse_CG(sparse, b, matrix(0, 25), tol = 1e-12,
                                  maxit = limit, precond = precond)
    expect_identical(dense$status, sparse_fit$status)
    expect_equal(dense$coefficients, sparse_fit$coefficients, tolerance = 1e-10)
    if (limit == 100L) {
      expect_true(dense$status)
      expect_lt(sqrt(sum((A %*% dense$coefficients - b)^2) / sum(b^2)), 1e-12)
    }
  }
  diagonal <- solve_dense_CG(A, b, matrix(0, 25), tol = 1e-12, maxit = 1L, precond = 1L)
  ic <- solve_dense_CG(A, b, matrix(0, 25), tol = 1e-12, maxit = 1L, precond = 2L)
  expect_false(diagonal$status)
  expect_true(ic$status)
  expect_equal(ic$coefficients, solve(A, b), tolerance = 1e-10)
  expect_equal(solve_dense_CG(A, b, matrix(0, 25), tol = 1e-12, maxit = 1L), ic)
  expect_equal(solve_sparse_CG(sparse, b, matrix(0, 25), tol = 1e-12, maxit = 1L), ic,
               tolerance = 1e-10)
  expect_identical(A, original)
})

test_that("default CG accuracy preserves LLT quality in the Gaussian PSF fit", {
  coefficients <- list()
  for (sparse in c(FALSE, TRUE)) for (native in c(FALSE, TRUE)) {
    sim <- simul(method = 1L, n = 500L, p = 500L, sparse = sparse,
                 seed = 1L, options = list(sd_noise = 1))
    lambda <- presetLambda("gic", 500, 500)
    start <- init_coefs(sim$X, sim$y, mu = sim$y, family = gaussian(), lambda = lambda)
    fits <- lapply(c(1L, 3L), function(method)
      l0glmFit(sim$X, sim$y, family = gaussian(), start = start, lambda = lambda,
        nonnegative = TRUE,  parallelize = FALSE,
        openmp = native, ncores = 1L, method = method))
    reference <- get.metrics(fits[[1L]]$coefficients, sim$X, sim$y,
                             family = gaussian(), beta_true = sim$beta_true)
    actual <- get.metrics(fits[[2L]]$coefficients, sim$X, sim$y,
                          family = gaussian(), beta_true = sim$beta_true)
    expect_identical(as.numeric(fits[[1L]]$coefficients) != 0,
                     as.numeric(fits[[2L]]$coefficients) != 0)
    expect_lt(abs(actual$deviance_true - reference$deviance_true) / 500, 5e-5)
    expect_equal(actual$R2_coefs, reference$R2_coefs, tolerance = 1e-5)
    coefficients[[length(coefficients) + 1L]] <- as.numeric(fits[[2L]]$coefficients)
  }
  for (beta in coefficients[-1L]) {
    expect_identical(beta != 0, coefficients[[1L]] != 0)
    expect_equal(beta, coefficients[[1L]], tolerance = 1e-5)
  }
})
