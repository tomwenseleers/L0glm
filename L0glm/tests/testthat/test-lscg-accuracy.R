test_that("LSCG diagnostics check the actual normal residual for every right-hand side", {
  set.seed(1861)
  A <- matrix(rnorm(60 * 12), 60, 12)
  b <- matrix(rnorm(60 * 3), 60, 3)
  start <- matrix(0, 12, 3)
  for (sparse in c(FALSE, TRUE)) {
    design <- if (sparse) as(as(Matrix::Matrix(A, sparse=TRUE), "generalMatrix"), "CsparseMatrix") else A
    solver <- if (sparse) solve_sparse_LSCG else solve_dense_LSCG
    fit <- solver(design, b, start, tol=1e-10, maxit=300L)
    residual <- crossprod(A, b - A %*% fit$coefficients)
    rhs <- crossprod(A, b)
    error <- max(sqrt(colSums(residual^2) / colSums(rhs^2)))
    expect_equal(fit$error, error, tolerance=1e-12)
    expect_identical(fit$status, is.finite(error) && error <= 1e-10)
    expect_true(fit$status)
    expect_equal(fit$coefficients, solve(crossprod(A), rhs), tolerance=1e-8)
  }
})

test_that("inexact LSCG refits the correct working system in R and native drivers", {
  set.seed(1862)
  X <- matrix(rnorm(60 * 12), 60, 12)
  y <- rnorm(60)
  args <- list(y=y, family=gaussian(), lambda=2, start=rep(.2, 12),
    dispersion=1, weights=seq(.7, 1.3, length.out=60), offset=seq(0, .2, length.out=60),
    prior.mean=rep(.03, 12), always.include=1L, normalize="none",
    maxit=1L, parallelize=FALSE, ncores=1L, local.search=FALSE)
  for (sparse in c(FALSE, TRUE)) for (native in c(FALSE, TRUE)) {
    common <- c(list(X=Matrix::Matrix(X, sparse=sparse), openmp=native), args)
    reference <- do.call(l0glmFit, c(common, list(method=1L)))
    for (rowaug in c(FALSE, TRUE)) {
      expect_warning(actual <- do.call(l0glmFit, c(common, list(method=4L,
        control.solver=controlSolverGen("eigen", 4L, tol=1e-14, maxit=1L, rowaug=rowaug)))),
        "Cholesky LLT")
      expect_equal(actual$coefficients, reference$coefficients, tolerance=1e-10)
      expect_equal(predict(actual), predict(reference), tolerance=1e-10)
    }
  }
})

test_that("default LSCG preserves LLT Gaussian PSF quality for both storage types and drivers", {
  for (sparse in c(FALSE, TRUE)) {
    sim <- simul(method=1L, n=500L, p=500L, sparse=sparse, seed=1L,
                 options=list(sd_noise=1))
    lambda <- presetLambda("gic", 500L, 500L)
    start <- init_coefs(sim$X, sim$y, mu=sim$y, family=gaussian(), lambda=lambda)
    for (native in c(FALSE, TRUE)) {
      args <- list(X=sim$X, y=sim$y, family=gaussian(), start=start, lambda=lambda,
        nonnegative=TRUE, parallelize=FALSE, openmp=native, ncores=1L)
      reference <- do.call(l0glmFit, c(args, list(method=1L)))
      warnings <- character()
      actual <- withCallingHandlers(do.call(l0glmFit, c(args, list(method=4L))),
        warning=function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")})
      expect_true(any(grepl("Cholesky LLT", warnings)))
      expect_identical(as.numeric(coef(actual)) != 0, as.numeric(coef(reference)) != 0)
      expect_equal(coef(actual), coef(reference), tolerance=1e-6)
      expect_equal(GIC(actual), GIC(reference), tolerance=1e-8)
      estimate <- get.metrics(coef(actual), sim$X, sim$y, family=gaussian(), beta_true=sim$beta_true)
      truth <- get.metrics(coef(reference), sim$X, sim$y, family=gaussian(), beta_true=sim$beta_true)
      expect_lt(abs(estimate$deviance_true - truth$deviance_true)/500, 5e-5)
      expect_equal(estimate$R2_coefs, truth$R2_coefs, tolerance=1e-7)
    }
  }
})

test_that("default LSCG preserves grouped Poisson deconvolution and threaded fits", {
  for (sparse in c(FALSE, TRUE)) {
    sim <- simulateSpikeTrain(n=200L, k=20L, m=4L, mean_beta=2000,
      family=poisson("identity"), sparse=sparse, seed=2L)
    lambda <- presetLambda("gic", 200L, 200L)
    Y <- as.matrix(sim$y)
    start <- vapply(seq_len(ncol(Y)), function(j)
      init_coefs(sim$X, Y[, j], mu=pmax(Y[, j], .1), family=sim$family,
        lambda=lambda), numeric(200L))
    for (native in c(FALSE, TRUE)) {
      args <- list(X=sim$X, y=sim$y, family=sim$family, start=start, lambda=lambda,
        nonnegative=TRUE, grouped=TRUE, parallelize=native, openmp=native,
        ncores=2L, local.search=FALSE)
      reference <- suppressWarnings(do.call(l0glmFit, c(args, list(method=1L))))
      actual <- suppressWarnings(do.call(l0glmFit, c(args, list(method=4L))))
      expect_identical(as.numeric(coef(actual)) != 0, as.numeric(coef(reference)) != 0)
      expect_equal(coef(actual), coef(reference), tolerance=1e-6)
      expect_equal(GIC(actual), GIC(reference), tolerance=1e-8)
    }
  }
})

test_that("invalid LSCG coefficients use the central fallback for quadratic family solvers", {
  set.seed(1863)
  X <- matrix(rnorm(60 * 12), 60, 12)
  A <- crossprod(X) + diag(12)
  b <- matrix(rnorm(12))
  control <- controlSolverGen("eigen", 4L)
  reference <- as.numeric(solve(A, b))
  for (coefficients in list(NULL, rep(NA_real_, 12), c(Inf, rep(1, 11)))) {
    solver <- function(...) list(status=TRUE, coefficients=coefficients)
    expect_warning(actual <- with_mocked_bindings(
      eigen_ridge_solve(A, b, rep(0, 12), control, 1L),
      solve_dense_LSCG=solver, .package="L0glm"), "Cholesky LLT")
    expect_equal(actual, reference, tolerance=1e-10)
  }
})
