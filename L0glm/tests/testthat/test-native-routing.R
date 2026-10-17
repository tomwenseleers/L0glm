test_that("native wide ridge updates match the explicit weighted objective", {
  set.seed(1610)
  n <- 18L; p <- 65L
  X <- matrix(runif(n * p), n, p)
  y <- rpois(n, 3)
  start <- seq(0.1, 1, length.out = p)
  target <- rep(0.02, p)
  d <- (3 / 4) / start^2
  w <- 1 / (y + 0.1)
  reference <- pmax(0, as.numeric(solve(crossprod(X, w * X) + diag(d),
                                          crossprod(X, w * y) + d * target)))
  for (sparse in c(FALSE, TRUE)) for (natural in c(FALSE, TRUE)) {
    fit <- L0nnpois(Matrix::Matrix(X, sparse = sparse), y,
                    start = start, prior.mean = target, lambda = 3,
                    maxit = 1L, thresh = 0, parallelize = FALSE,
                    control.solver = control.solver.gen("eigen", 1L,
                                                        NaturalOrdering = natural))
    expect_equal(as.numeric(fit$coefficients), reference, tolerance = 1e-10)
    expect_identical(inherits(fit$X, "sparseMatrix"), sparse)
  }
})

test_that("native and R wide multichannel fits agree with storage preserved", {
  set.seed(1611)
  X <- matrix(runif(18 * 55), 18, 55)
  Y <- matrix(rpois(54, 3), 18, 3)
  for (sparse in c(FALSE, TRUE)) for (grouped in c(FALSE, TRUE)) {
    args <- list(X = Matrix::Matrix(X, sparse = sparse), y = Y,
                 family = poisson(identity), lambda = 2, start = 0.5,
                 grouped = grouped, relaxed = TRUE, maxit = 8L,
                 group.maxit = 2L, parallelize = FALSE)
    native <- do.call(L0nnpois, args)
    general <- do.call(L0glm.fit, c(args, list(nonnegative = TRUE, openmp = FALSE)))
    expect_equal(as.matrix(native$coefficients), as.matrix(general$coefficients),
                 tolerance = 1e-7)
    expect_identical(inherits(native$X, "sparseMatrix"), sparse)
  }
})
