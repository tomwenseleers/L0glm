# Testing suite for control.solver.gen

test_that("control.solver.gen warns about unsupported standardization for non-bigGlm solvers", {
  output <- capture.output({
    params <- control.solver.gen(solver = "osqp", standardize = TRUE, verbose = TRUE)
  })
  expect_match(output, "Selected solver does not support standardization. To enable standardization, please set the solver to bigGlm\\.")
})


test_that("control.solver.gen handles invalid solvers", {
  expect_error(control.solver.gen(solver = "invalid_solver"), 
               "Chosen solver is not available.")
})

test_that("control.solver.gen handles missing method for eigen solver", {
  output <- capture.output({
    params <- control.solver.gen(solver = "eigen", verbose = TRUE)
  })
  expect_match(output, "No method for Eigen solver was provided. Method is set to 7.")
  expect_equal(params$method, "chol")
})

test_that("control.solver.gen handles osqp solver defaults and custom parameters", {
  # Default parameters
  params <- control.solver.gen(solver = "osqp")
  expect_equal(params$max_iter, 1e5)
  expect_equal(params$eps_abs, 1e-7)
  expect_equal(params$eps_rel, 1e-7)
  expect_true(is.null(params$rowaug) || params$rowaug)
  
  # Custom parameters
  params_custom <- control.solver.gen(solver = "osqp", maxit = 500, reltol = 1e-5, abstol = 1e-6, rowaug = FALSE)
  expect_equal(params_custom$max_iter, 500)
  expect_equal(params_custom$eps_abs, 1e-6)
  expect_equal(params_custom$eps_rel, 1e-5)
  expect_false(params_custom$rowaug)
})

test_that("control.solver.gen handles bigGlm solver defaults, aliases, and custom parameters", {
  # Default parameters
  params <- control.solver.gen(solver = "bigGlm")
  expect_equal(params$maxit, 1e6)
  expect_equal(params$thresh, 1e-10)
  expect_true(params$standardize)
  expect_true(is.null(params$rowaug) || params$rowaug)
  expect_equal(control.solver.gen(solver = "glmnet"), params)
  
  # Custom parameters
  params_custom <- control.solver.gen(solver = "bigGlm", maxit = 1000, abstol = 1e-5, standardize = FALSE, rowaug = TRUE)
  expect_equal(params_custom$maxit, 1000)
  expect_equal(params_custom$thresh, 1e-5)
  expect_false(params_custom$standardize)
  expect_true(is.null(params_custom$rowaug) || params_custom$rowaug)
})

test_that("control.solver.gen does not warn if rowaug is FALSE for bigGlm", {
  params <- control.solver.gen(solver = "bigGlm", rowaug = FALSE)
  expect_true(is.null(params$rowaug) || !params$rowaug)
})

test_that("control.solver.gen handles nnls solver", {
  params <- control.solver.gen(solver = "nnls")
  expect_equal(params, list())
})

test_that("control.solver.gen sets default method for eigen solver", {
  params <- control.solver.gen(solver = "eigen")
  expect_equal(params$method, "chol")
  expect_true(is.null(params$rowaug) || params$rowaug)
})

test_that("control.solver.gen handles all methods for eigen solver", {
  for (method in 1:11) {
    params <- control.solver.gen(solver = "eigen", method = method)
    
    if (method == 1 || method == 2) {
      expect_equal(params$method, "QR")
      expect_false(params$rowaug)
      expect_equal(params$pivot, ifelse(method == 1, 0L, 1L))  # Explicit pivot check
      expect_equal(params$ordering, 1L)
    } else if (method == 3 || method == 4) {
      expect_equal(params$method, "LSCG")
      expect_equal(params$tol, 1E-7)
      expect_true(is.null(params$rowaug) || params$rowaug)
      expect_equal(params$precond, ifelse(method == 3, 0L, 1L))  # Preconditioner check
    } else if (method %in% 5:6 || method == 9) {
      expect_equal(params$method, "CG")
      expect_equal(params$tol, 1E-7)
      expect_true(is.null(params$rowaug) || params$rowaug)
      expect_equal(params$precond, ifelse(method == 5, 0L, ifelse(method == 6, 1L, 2L)))  # Preconditioner check
    } else if (method == 7 || method == 8) {
      expect_equal(params$method, "chol")
      expect_equal(params$pivot, ifelse(method == 7, 0L, 1L))  # Explicit pivot check
      expect_true(is.null(params$rowaug) || params$rowaug)
    } else if (method == 10 || method == 11) {
      expect_equal(params$method, "SVD")
      expect_equal(params$svd, ifelse(method == 10, 0L, 1L))
      expect_equal(params$precond, 0L)
      expect_equal(params$threshold, 0)
      expect_true(is.null(params$rowaug) || params$rowaug)
    }
  }
})


test_that("control.solver.gen handles invalid methods for eigen solver", {
  expect_error(control.solver.gen(solver = "eigen", method = 12L), 
               "Chosen method for Eigen solvers is unavailable. Please choose a method \\(1-11\\).")
})


test_that("control.solver.gen handles negative maxit values", {
  # osqp
  output_osqp <- capture.output({
    params <- control.solver.gen(solver = "osqp", maxit = -1, verbose = TRUE)
  })
  expect_match(output_osqp, "Control.solver\\$maxit needs to be a positive integer fot the chosen solver\\. Value was set to solver default \\(10000\\)\\.")
  
  # bigGlm
  output_glmnet <- capture.output({
    params <- control.solver.gen(solver = "bigGlm", maxit = -1, verbose = TRUE)
  })
  expect_match(output_glmnet, "Control.solver\\$maxit needs to be a positive integer fot the chosen solver\\. Value was set to solver default \\(1e5\\)\\.")
})


test_that("control.solver.gen handles missing parameters gracefully", {
  # osqp
  params_osqp <- control.solver.gen(solver = "osqp")
  expect_equal(params_osqp$max_iter, 1e5)
  expect_equal(params_osqp$eps_abs, 1e-7)
  expect_equal(params_osqp$eps_rel, 1e-7)
  
  # bigGlm
  params_glmnet <- control.solver.gen(solver = "bigGlm")
  expect_equal(params_glmnet$maxit, 1e6)
  expect_equal(params_glmnet$thresh, 1e-10)
})

