test_that("LSCG solves finite ridge systems with nonzero prior means", {
  set.seed(715)
  X <- matrix(rnorm(80 * 20), 80, 20)
  y <- rnorm(80)
  lambdas <- 10 ^ seq(-6, 8, length.out = ncol(X))
  target <- c(0.1, -0.2, 0.3, rep(0, 17))
  start <- 1 / sqrt(pmax(lambdas, 1))
  reference <- lmridge(X, y, lambdas = lambdas, prior.mean = target)
  control <- control.solver.gen("eigen", method = 4L, tol = 1e-11, maxit = 300L)
  for (sparse in c(FALSE, TRUE)) {
    estimate <- lmridge(Matrix::Matrix(X, sparse = sparse), y,
                        lambdas = lambdas, prior.mean = target, start = start,
                        method = 4L, control.solver = control)
    expect_equal(estimate, reference, tolerance = 1e-8)
  }
})

test_that("LSCG iteration limits retain finite iterates and report nonconvergence", {
  set.seed(716)
  X <- matrix(rnorm(80 * 20), 80, 20)
  y <- rnorm(80)
  expect_warning(
    estimate <- lmridge(X, y, lambdas = rep(1, 20), method = 4L,
                        control.solver = control.solver.gen("eigen", method = 4L,
                                                            tol = 1e-14, maxit = 1L)),
    "LSCG did not converge.*finite approximate solution"
  )
  native <- solve_dense_LSCG(A=rbind(X,diag(20)),b=matrix(c(y,rep(0,20))),
                             x0=matrix(0,20),tol=1e-14,maxit=1L,precond=TRUE)
  expect_equal(estimate,as.numeric(native$coefficients),tolerance=1e-12)
  expect_false(isTRUE(native$status))
  expect_null(lscg_coefficients(list(status=FALSE,coefficients=NULL),1))
  expect_null(lscg_coefficients(list(status=FALSE,coefficients=c(1,Inf)),1))
})

test_that("OSQP preserves bounds and nonzero ridge targets in original coordinates", {
  skip_if_not_installed("osqp")
  column_scale <- c(0.5, 1, 4, 1)
  X <- diag(column_scale)
  y <- c(1, 2, 3, 4)
  lambdas <- c(0.1, 2, 1, 20)
  target <- c(0.1, -0.2, 0.3, 0)
  lower <- c(0, -0.3, 0, 0)
  upper <- c(0.8, 0.2, 1, Inf)
  reference <- pmin(upper, pmax(lower,
    (column_scale * y + lambdas * target) / (column_scale ^ 2 + lambdas)))
  for (sparse in c(FALSE, TRUE)) {
    for (rowaug in c(FALSE, TRUE)) {
      estimate <- lmridge(Matrix::Matrix(X, sparse = sparse), y,
                          lambdas = lambdas, prior.mean = target,
                          lower = lower, upper = upper, start = c(0.5, 0.1, 0.1, 0.1),
                          solver = "osqp",
                          control.solver = control.solver.gen("osqp", rowaug = rowaug,
                                                              abstol = 1e-9, reltol = 1e-9))
      expect_equal(estimate, reference, tolerance = 1e-6)
      expect_true(all(estimate >= lower & estimate <= upper))
    }
  }
})

test_that("OSQP iteration limits are reported instead of silently accepted", {
  skip_if_not_installed("osqp")
  set.seed(717)
  X <- matrix(rnorm(80 * 20), 80, 20)
  expect_warning(
    lmridge(X, rnorm(80), lambdas = rep(1, 20), nonnegative = TRUE, solver = "osqp",
            control.solver = control.solver.gen("osqp", maxit = 1L)),
    "OSQP did not converge.*maximum iterations reached"
  )
})

test_that("OSQP failure statuses retain bounded starts rather than placeholders", {
  start <- c(-0.1,0.4,2)
  lower <- c(0,0,0); upper <- c(1,1,1)
  for (status in c("problem non convex","primal infeasible","dual infeasible",
                    "primal infeasible inaccurate","unsolved")) {
    # Infeasibility/nonconvexity may come with a finite OSQP sentinel vector.
    result <- list(info=list(status=status),x=rep(2.14e9,3))
    expect_warning(beta <- osqp_coefficients(result,start,lower,upper,3),
                   "OSQP failed in IRLS iteration 3.*retaining previous coefficients")
    expect_equal(beta,c(0,0.4,1))
  }
  for (x in list(NULL,c(1,NA,2),c(1,Inf,2),numeric(),1)) {
    expect_warning(beta <- osqp_coefficients(list(info=list(status="solved"),x=x),
                                             start,lower,upper,3),
                   "invalid coefficient vector.*retaining previous coefficients")
    expect_equal(beta,c(0,0.4,1))
  }
  expect_equal(osqp_coefficients(list(info=list(status="solved"),x=c(.1,.2,.3)),
                                start,lower,upper,3),c(.1,.2,.3))
})

test_that("OSQP setup exceptions are warnings and preserve starting coefficients", {
  skip_if_not_installed("osqp")
  testthat::local_mocked_bindings(osqp=function(...) stop("nonconvex setup"),.package="osqp")
  expect_warning(beta <- lmridge(diag(3),c(1,2,3),solver="osqp",
                                 nonnegative=TRUE,start=c(.1,.2,.3)),
                 "OSQP failed.*nonconvex setup.*retaining previous coefficients")
  expect_equal(beta,c(.1,.2,.3))
})

test_that("dense storage of an exactly sparse design preserves native fits", {
  set.seed(718)
  X <- Matrix::rsparsematrix(80, 50, density = 0.04)
  X@x <- abs(X@x)
  X[1, 1] <- 1e-30
  y <- matrix(rpois(160, 2), 80, 2)
  args <- list(y = y, lambda = 0.5, weights = c(0, rep(1, 79)),
               offset = rep(0.2, 80), prior.mean = rep(0.01, 50),
               normalize = TRUE, grouped = TRUE, relaxed = TRUE,
               parallelize = FALSE, maxit = 10L, group.maxit = 2L)
  sparse <- do.call(L0nnpois, c(list(X = X), args))
  dense <- do.call(L0nnpois, c(list(X = as.matrix(X)), args))
  expect_equal(dense, sparse, tolerance = 1e-12)
})
