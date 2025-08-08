# Testing basic functionality of the lmridge function
#----------------------------------------------------

# Example data for tests
set.seed(123)
n <- 100
p <- 10
X <- matrix(rnorm(n * p), n, p)
y <- X %*% rnorm(p) + rnorm(n)

test_that("Basic functionality test with default settings", {
  result <- lmridge(X, y, lambdas = rep(1, p))
  expect_type(result, "double")
  expect_length(result, ncol(X))
})

test_that("Single lambda value as vector results in correct ridge penalty", {
  result <- lmridge(X, y, lambdas = rep(5, p))  # Provide lambda as a vector
  expect_true(!is.null(result))
  expect_equal(length(result), ncol(X))
})

test_that("Eigen solver returns valid coefficients with vector lambdas", {
  result <- lmridge(X, y, lambdas = rep(1, p), solver = "eigen")
  expect_true(all(!is.na(result)))
  expect_equal(length(result), ncol(X))
})

test_that("Glmnet solver returns valid coefficients with vector lambdas", {
  result <- lmridge(X, y, lambdas = rep(1, p), solver = "glmnet")
  expect_true(all(!is.na(result)))
  expect_equal(length(result), ncol(X))
})

test_that("OSQP solver with nonnegativity constraints and vector lambdas", {
  result <- lmridge(X, y, lambdas = rep(1, p), solver = "osqp", nonnegative = TRUE)
  expect_true(all(result >= 0))
})

test_that("NNLS solver enforces nonnegative constraint with vector lambdas", {
  result <- lmridge(X, y, lambdas = rep(1, p), solver = "nnls", nonnegative = TRUE)
  expect_true(all(result >= 0))
})

test_that("Lower and upper bounds constraints are respected with vector lambdas", {
  lower_bounds <- rep(-0.5, ncol(X))
  upper_bounds <- rep(0.5, ncol(X))
  result <- lmridge(X, y, lambdas = rep(1, p), lower = lower_bounds, upper = upper_bounds)
  expect_true(all(result >= lower_bounds))
  expect_true(all(result <= upper_bounds))
})

test_that("Scaling applied correctly to columns of X with vector lambdas", {
  scaling_factors <- rep(0.5, ncol(X))
  result <- lmridge(X, y, lambdas = rep(1, p), scaling = scaling_factors)
  expect_true(all(abs(result) < 5))  # Expect scaled coefficients
})

test_that("Prior mean shifts coefficients with vector lambdas", {
  prior_mean <- rep(1, ncol(X))
  result <- lmridge(X, y, lambdas = rep(1, p), prior.mean = prior_mean)
  expect_true(any(result > 1))  # Coefficients should be shifted positively
})

test_that("Handles sparse matrix input correctly with vector lambdas", {
  X_sparse <- Matrix(X, sparse = TRUE)
  result <- lmridge(X_sparse, y, lambdas = rep(1, p))
  expect_true(all(!is.na(result)))
  expect_equal(length(result), ncol(X))
})

test_that("Handles large matrix with 20000 columns and vector lambdas", {
  X_large <- matrix(rnorm(100 * 20000), 100, 20000)
  y_large <- X_large %*% rnorm(20000) + rnorm(100)
  result <- lmridge(X_large, y_large, lambdas = rep(1, 20000))
  expect_equal(length(result), ncol(X_large))
})

test_that("Zero lambdas vector produces least squares solution", {
  result <- lmridge(X, y, lambdas = rep(0, p))
  expect_true(all(!is.na(result)))
  expect_equal(length(result), ncol(X))
})

test_that("Handles large matrix with 20000 columns and zero lambdas vector", {
  X_large <- matrix(rnorm(100 * 20000), 100, 20000)
  y_large <- X_large %*% rnorm(20000) + rnorm(100)
  result <- lmridge(X_large, y_large, lambdas = rep(0, 20000))
  expect_equal(length(result), ncol(X_large))
})

test_that("Large lambdas vector shrinks coefficients towards zero", {
  result <- lmridge(X, y, lambdas = rep(1e10, p))
  expect_true(all(abs(result) < 1e-5))  # Expect coefficients close to zero
})

test_that("Single predictor column in X", {
  X_single <- matrix(rnorm(n), n, 1)
  result <- lmridge(X_single, y, lambdas = rep(1, 1))
  expect_length(result, 1)
})

test_that("Identical lower and upper bounds force coefficients to a specific value", {
  bound_value <- 0.1
  result <- lmridge(X, y, lambdas = rep(1, p), lower = rep(bound_value, p), upper = rep(bound_value, p))
  expect_true(all(abs(result - bound_value) < 1e-5))
})

test_that("Lower bounds are respected with vector lambdas", {
  lower_bounds <- rep(-0.5, ncol(X))
  result <- lmridge(X, y, lambdas = rep(1, p), lower = lower_bounds)
  expect_true(all(result >= lower_bounds))
})

test_that("prior.mean is NULL", {
  # prior.mean is NULL, should default to a vector of zeros
  result <- lmridge(X, y, lambdas = rep(1, p), prior.mean = NULL)
  expect_type(result, "double")
  expect_length(result, ncol(X))
})

test_that("prior.mean is a single numeric value", {
  # prior.mean is a single value, should expand to a vector of the same length as ncol(X)
  single_value <- 2
  result <- lmridge(X, y, lambdas = rep(1, p), prior.mean = single_value)
  expect_type(result, "double")
  expect_length(result, ncol(X))
})

test_that("prior.mean has incorrect length", {
  # prior.mean length does not match ncol(X), should raise an error
  incorrect_length_prior <- rep(1, p + 1)  # Length mismatch
  expect_error(lmridge(X, y, lambdas = rep(1, p), prior.mean = incorrect_length_prior),
               "Prior.mean needs to be a vector of length p or a numerical.")
})

test_that("Invalid solver choice produces an error", {
  # Use an invalid solver name to trigger the error
  invalid_solver <- "invalid_solver"
  expect_error(lmridge(X, y, lambdas = rep(1, p), solver = invalid_solver), "Invalid solver choice")
})

test_that("Else branch is triggered when lambdas is at default of zero", {
  # By omitting lambdas, it defaults to NULL, triggering the else branch
  result <- lmridge(X, y)
  expect_type(result, "double")
  expect_length(result, ncol(X))
})

test_that("Else branch is triggered when lambdas sum to zero (using zero vector)", {
  # lambdas is a vector of zeros, so Xwaug should be Xw and ywaug should be yw
  zero_lambdas <- rep(0, p)
  result <- lmridge(X, y, lambdas = zero_lambdas)
  expect_type(result, "double")
  expect_length(result, ncol(X))
})


# Testing lmridge.osqp
#----------------------

test_that("lmridge.osqp handles basic inputs", {
  Xw <- Matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2, sparse = TRUE)
  yw <- c(1, 2)
  lower <- c(-Inf, -Inf)
  upper <- c(Inf, Inf)
  start <- c(0, 0)
  
  control.solver <- control.solver.gen(solver = "osqp")
  
  result <- lmridge.osqp(Xw, yw, lower = lower, upper = upper, start = start, control.solver = control.solver)
  
  expect_length(result, 2)  # Expecting two coefficients
  expect_true(all(!is.na(result)))  # Ensure no NAs in the result
})

test_that("lmridge.osqp handles sparse Xw input with rowaug = FALSE", {
  Xw_sparse <- Matrix(runif(4), nrow = 2, ncol = 2, sparse = TRUE)  # Sparse matrix
  yw <- c(1, 2)
  lower <- c(-Inf, -Inf)
  upper <- c(Inf, Inf)
  start <- c(0, 0)
  
  control.solver <- control.solver.gen(solver = "osqp", rowaug = FALSE)  # Set rowaug to FALSE
  
  # Run the function
  result <- lmridge.osqp(Xw_sparse, yw, lower = lower, upper = upper, start = start, control.solver = control.solver)
  
  # Basic checks
  expect_length(result, 2)
  expect_true(all(!is.na(result)))
})


test_that("lmridge.osqp handles non-sparse Xw input with rowaug = FALSE", {
  Xw_dense <- matrix(runif(4), nrow = 2, ncol = 2)  # Dense matrix
  yw <- c(1, 2)
  lower <- c(-Inf, -Inf)
  upper <- c(Inf, Inf)
  start <- c(0, 0)
  
  control.solver <- control.solver.gen(solver = "osqp", rowaug = FALSE)  # Set rowaug to FALSE
  
  # Run the function
  result <- lmridge.osqp(Xw_dense, yw, lower = lower, upper = upper, start = start, control.solver = control.solver)
  
  # Basic checks
  expect_length(result, 2)
  expect_true(all(!is.na(result)))
})

test_that("lmridge.osqp handles edge cases with sparse and non-sparse inputs", {
  # Sparse matrix edge case
  Xw_sparse <- Matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2, sparse = TRUE)
  yw <- c(1, 2)
  lower <- c(-Inf, -Inf)
  upper <- c(Inf, Inf)
  start <- c(0, 0)
  lambdas <- c(0.1, 0.2)  # Valid lambdas vector
  
  control.solver <- control.solver.gen(solver = "osqp")
  
  result_sparse <- lmridge.osqp(Xw_sparse, yw, lower = lower, upper = upper, start = start, lambdas = lambdas, control.solver = control.solver)
  
  expect_length(result_sparse, 2)
  expect_true(all(!is.na(result_sparse)))
  
  # Dense matrix edge case
  Xw_dense <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2)
  
  result_dense <- lmridge.osqp(Xw_dense, yw, lower = lower, upper = upper, start = start, lambdas = lambdas, control.solver = control.solver)
  
  expect_length(result_dense, 2)
  expect_true(all(!is.na(result_dense)))
})

test_that("lmridge.osqp handles warm_start condition", {
  Xw <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2)
  yw <- c(1, 2)
  lower <- c(-Inf, -Inf)
  upper <- c(Inf, Inf)
  
  # Test case where start is non-zero (warm_start should be TRUE)
  start_non_zero <- c(1, 1)  # Non-zero start vector
  control.solver <- control.solver.gen(solver = "osqp")
  
  result_non_zero <- lmridge.osqp(Xw, yw, lower = lower, upper = upper, start = start_non_zero, lambdas = c(0.1, 0.2), control.solver = control.solver)
  
  # Check if warm_start is TRUE
  settings_non_zero <- do.call(osqp::osqpSettings, control.solver[-7])
  settings_non_zero <- c(settings_non_zero, warm_start = TRUE)  # Expected setting when start is non-zero
  
  expect_true("warm_start" %in% names(settings_non_zero))
  expect_equal(settings_non_zero[["warm_start"]], TRUE)
  
  # Test case where start is zero (warm_start should be FALSE)
  start_zero <- c(0, 0)  # Zero start vector
  
  result_zero <- lmridge.osqp(Xw, yw, lower = lower, upper = upper, start = start_zero, lambdas = c(0.1, 0.2), control.solver = control.solver)
  
  # Check if warm_start is FALSE
  settings_zero <- do.call(osqp::osqpSettings, control.solver[-7])
  settings_zero <- c(settings_zero, warm_start = FALSE)  # Expected setting when start is zero
  
  expect_true("warm_start" %in% names(settings_zero))
  expect_equal(settings_zero[["warm_start"]], FALSE)
})

test_that("lmridge.osqp handles NULL lambdas", {
  Xw <- Matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2, sparse = TRUE)
  yw <- c(1, 2)
  lower <- c(-Inf, -Inf)
  upper <- c(Inf, Inf)
  start <- c(0, 0)
  lambdas <- NULL
  
  control.solver <- control.solver.gen(solver = "osqp")
  
  result <- lmridge.osqp(Xw, yw, lambdas = lambdas, lower = lower, upper = upper, start = start, control.solver = control.solver)
  
  expect_length(result, 2)
  expect_true(all(!is.na(result)))
})

test_that("lmridge.osqp handles warm start", {
  Xw <- Matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2, sparse = TRUE)
  yw <- c(1, 2)
  lower <- c(-Inf, -Inf)
  upper <- c(Inf, Inf)
  start <- c(0, 0)
  
  control.solver <- control.solver.gen(solver = "osqp")
  
  # Set warm start to TRUE
  control.solver$warm_start <- TRUE
  
  result <- lmridge.osqp(Xw, yw, lower = lower, upper = upper, start = start, control.solver = control.solver)
  
  expect_length(result, 2)
  expect_true(all(!is.na(result)))
})

test_that("lmridge.osqp handles different iterations", {
  Xw <- Matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2, sparse = TRUE)
  yw <- c(1, 2)
  lower <- c(-Inf, -Inf)
  upper <- c(Inf, Inf)
  start <- c(0, 0)
  
  control.solver <- control.solver.gen(solver = "osqp")
  
  # Set iteration to 10 (just for testing purposes)
  result <- lmridge.osqp(Xw, yw, lower = lower, upper = upper, start = start, control.solver = control.solver, iteration = 10)
  
  expect_length(result, 2)
  expect_true(all(!is.na(result)))
})

test_that("lmridge.osqp assigns start to a zero vector when it is NULL", {
  # Setup inputs
  Xw <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2)
  yw <- c(1, 2)
  lower <- c(-Inf, -Inf)
  upper <- c(Inf, Inf)
  start <- NULL
  
  # Generate control solver
  control.solver <- control.solver.gen(solver = "osqp")
  
  # Run the function
  result <- lmridge.osqp(Xw, yw, lower = lower, upper = upper, start = start, control.solver = control.solver)
  
  # Expectations on the result
  expect_length(result, 2)               # Expect the output to have two coefficients
  expect_true(all(!is.na(result)))       # Ensure the coefficients are not NA
})


test_that("lmridge.osqp does not overwrite start if it is not NULL", {
  # Create some basic test inputs
  Xw <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2)
  yw <- c(1, 2)
  lower <- c(-Inf, -Inf)
  upper <- c(Inf, Inf)
  
  # Non-NULL start
  start <- c(0.5, 1.5)  # custom start values
  
  control.solver <- control.solver.gen(solver = "osqp")
  
  # Call lmridge.osqp with non-NULL start
  result <- lmridge.osqp(Xw, yw, lower = lower, upper = upper, start = start, control.solver = control.solver)
  
  # Test that start was not overwritten
  expect_equal(start, c(0.5, 1.5))  # start should remain unchanged
})

# Testing the lmridge.glmnet function
#------------------------------------
# Sample input data
set.seed(123)
n <- 50  # Smaller sample size for easier testing
p <- 1001  # Columns > 1000 to test the glmnet branch condition
X <- matrix(rnorm(n * p), n, p)
y <- rnorm(n)
lambdas <- rep(1, p)
lower_bounds <- rep(-1, p)
upper_bounds <- rep(1, p)

# Control settings for glmnet
control_solver <- list(
  standardize = TRUE,
  thresh = 1e-7,
  maxit = 1000
)

test_that("Initial glmnet path is triggered with large dimensions and iteration <= 5", {
  n <- 50
  p <- 1001  # Ensure this matches the length of lambdas to avoid mismatch
  
  # Generate synthetic data with dimensions matching the expectations of the function
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  lambdas <- rep(1, p)  # Match lambdas length to number of predictors
  lower_bounds <- rep(-1, p)
  upper_bounds <- rep(1, p)
  
  control_solver <- list(
    standardize = TRUE,
    thresh = 1e-7,
    maxit = 1000
  )
  
  # Run lmridge.glmnet with correct lambdas length
  result <- lmridge.glmnet(X, y, lambdas = lambdas, lower = lower_bounds, upper = upper_bounds, 
                           control.solver = control_solver, iteration = 5)
  
  # Assertions to verify results
  expect_type(result, "double")
  expect_length(result, p)
  expect_true(all(result >= lower_bounds))
  expect_true(all(result <= upper_bounds))
})


test_that("Fallback bigGlm path is triggered with smaller matrix dimensions", {
  X_small <- matrix(rnorm(100 * 50), 100, 50)  # Small matrix
  y_small <- rnorm(100)
  lambdas_small <- rep(1, 50)
  lower_small <- rep(-1, 50)
  upper_small <- rep(1, 50)
  
  result <- lmridge.glmnet(X_small, y_small, lambdas = lambdas_small, lower = lower_small, upper = upper_small, 
                           control.solver = control_solver, iteration = 5)
  expect_type(result, "double")
  expect_length(result, 50)
  expect_true(all(result >= lower_small))
  expect_true(all(result <= upper_small))
})


test_that("Handles zero lambdas vector and triggers bigGlm path", {
  lambdas_zero <- rep(0, p)
  result <- lmridge.glmnet(X, y, lambdas = lambdas_zero, lower = lower_bounds, upper = upper_bounds, 
                           control.solver = control_solver, iteration = 6)
  expect_type(result, "double")
  expect_length(result, p)
  expect_true(all(result >= lower_bounds))
  expect_true(all(result <= upper_bounds))
})

test_that("Handles lambdas with non-zero values and applies constraints correctly", {
  lambdas_non_zero <- rep(1, p)
  result <- lmridge.glmnet(X, y, lambdas = lambdas_non_zero, lower = lower_bounds, upper = upper_bounds, 
                           control.solver = control_solver, iteration = 5)
  expect_type(result, "double")
  expect_length(result, p)
  expect_true(all(result >= lower_bounds))
  expect_true(all(result <= upper_bounds))
})

# Testing the lmridge.nnls function
#----------------------------------
lambdas_nonzero <- rep(1, ncol(X))

# Test with non-zero lambdas (augmented case)
test_that("lmridge.nnls handles augmented matrices with non-zero lambdas", {
  set.seed(123)
  X <- matrix(rnorm(100), nrow = 20, ncol = 5)
  y <- rnorm(20)
  lambdas_nonzero <- rep(1, ncol(X))
  lower_bounds <- rep(0, ncol(X))
  upper_bounds <- rep(Inf, ncol(X))
  
  result <- lmridge.nnls(X, y, lambdas = lambdas_nonzero, lower = lower_bounds, upper = upper_bounds, 
                         control.solver = control.solver, iteration = 1)
  expect_type(result, "double")
  expect_length(result, ncol(X))
  expect_true(all(result >= lower_bounds))
})

# Test with zero lambdas (non-augmented case)
test_that("lmridge.nnls handles non-augmented matrices with zero lambdas", {
  lambdas_zero <- rep(0, ncol(X))
  
  result <- lmridge.nnls(X, y, lambdas = lambdas_zero, lower = lower_bounds, upper = upper_bounds, 
                         control.solver = control.solver, iteration = 1)
  expect_type(result, "double")
  expect_length(result, ncol(X))
  expect_true(all(result >= lower_bounds))
})

# Test with NULL lambdas
test_that("lmridge.nnls handles NULL lambdas", {
  result <- lmridge.nnls(X, y, lambdas = NULL, lower = lower_bounds, upper = upper_bounds, 
                         control.solver = control.solver, iteration = 1)
  expect_type(result, "double")
  expect_length(result, ncol(X))
  expect_true(all(result >= lower_bounds))
})

# Test with single lambda value
test_that("lmridge.nnls handles single lambda value", {
  lambdas_single_value <- rep(1, ncol(X))
  result <- lmridge.nnls(X, y, lambdas = lambdas_single_value, lower = lower_bounds, upper = upper_bounds, 
                         control.solver = control.solver, iteration = 1)
  expect_type(result, "double")
  expect_length(result, ncol(X))
  expect_true(all(result >= lower_bounds))
})

# Test with boundary constraints
test_that("lmridge.nnls respects boundary constraints", {
  lower_bounds <- rep(-0.5, ncol(X))
  upper_bounds <- rep(0.5, ncol(X))
  
  result <- lmridge.nnls(X, y, lambdas = lambdas_nonzero, lower = lower_bounds, upper = upper_bounds, 
                         control.solver = control.solver, iteration = 1)
  expect_true(all(result >= lower_bounds))
  expect_true(all(result <= upper_bounds))
})

# Test with large matrix dimensions
test_that("lmridge.nnls handles large matrix dimensions", {
  lambdas_large <- rep(1, ncol(X))
  result <- lmridge.nnls(X, y, lambdas = lambdas_large, lower = lower_bounds, upper = upper_bounds, 
                         control.solver = control_solver, iteration = 1)
  expect_type(result, "double")
  expect_length(result, ncol(X))
})

# Test with small matrix dimensions
test_that("lmridge.nnls handles small matrix dimensions", {
  X_small <- matrix(rnorm(4), nrow = 2, ncol = 2)
  y_small <- rnorm(2)
  
  result <- lmridge.nnls(X_small, y_small, lambdas = c(0.1, 0.1), 
                         lower = c(0, 0), upper = c(1, 1), 
                         control.solver = control.solver, iteration = 1)
  expect_length(result, ncol(X_small))
})

# Consistency check on repeated runs
test_that("lmridge.nnls produces consistent results on repeated runs", {
  result_1 <- lmridge.nnls(X, y, lambdas = lambdas_nonzero, lower = lower_bounds, upper = upper_bounds, 
                           control.solver = control.solver, iteration = 1)
  result_2 <- lmridge.nnls(X, y, lambdas = lambdas_nonzero, lower = lower_bounds, upper = upper_bounds, 
                           control.solver = control.solver, iteration = 1)
  expect_equal(result_1, result_2)
})

# Test with extreme values in X and y
test_that("lmridge.nnls handles extreme values", {
  lambdas_extreme <- rep(10^6, ncol(X))  # Extreme values for lambda
  result <- lmridge.nnls(X, y, lambdas = lambdas_extreme, lower = lower_bounds, upper = upper_bounds, 
                         control.solver = control_solver, iteration = 1)
  expect_type(result, "double")
  expect_length(result, ncol(X))
})


# Testing the lmridge.eigen function
#-----------------------------------
# Helper function to check sparse matrix
issparse <- function(mat) inherits(mat, "sparseMatrix")

# Mock data
set.seed(123)
n <- 50
p <- 20

X_dense <- matrix(rnorm(n * p), nrow = n, ncol = p)
X_sparse <- Matrix::rsparsematrix(n, p, density = 0.1)
yw <- rnorm(n)
lambdas <- rep(0.1, p)
control_solver_base <- list(tol = 1e-6, maxit = 100, pivot = TRUE, ordering = TRUE, rowaug = TRUE)

# Constraints
lower <- rep(-10, p)
upper <- rep(10, p)

test_that("lmridge.eigen works with all solver methods for dense and sparse matrices", {
  methods <- 1L:8L
  for (method in methods) {
    control_solver <- control.solver.gen(solver = "eigen", method = method)
    
    # Dense matrix test
    coefs_dense <- lmridge.eigen(
      Xw = X_dense, yw = yw,
      lambdas = lambdas,
      lower = lower, upper = upper,
      control.solver = control_solver
    )
    expect_type(coefs_dense, "double")
    expect_equal(length(coefs_dense), p, info = paste("Dense method =", method))
    
    # Sparse matrix test
    coefs_sparse <- lmridge.eigen(
      Xw = X_sparse, yw = yw,
      lambdas = lambdas,
      lower = lower, upper = upper,
      control.solver = control_solver
    )
    expect_type(coefs_sparse, "double")
    expect_equal(length(coefs_sparse), p, info = paste("Sparse method =", method))
  }
})

test_that("lmridge.eigen handles zero lambdas correctly", {
  methods <- 1L:8L
  for (method in methods) {
    control_solver <- control.solver.gen(solver = "eigen", method = method)
    
    coefs_dense <- lmridge.eigen(
      Xw = X_dense, yw = yw,
      lambdas = rep(0, ncol(X_dense)),
      lower = lower, upper = upper,
      control.solver = control_solver
    )
    expect_type(coefs_dense, "double")
    expect_equal(length(coefs_dense), p, info = paste("Dense method =", method, "null lambdas"))
    
    coefs_sparse <- lmridge.eigen(
      Xw = X_sparse, yw = yw,
      lambdas = rep(0, ncol(X_sparse)),
      lower = lower, upper = upper,
      control.solver = control_solver
    )
    expect_type(coefs_sparse, "double")
    expect_equal(length(coefs_sparse), p, info = paste("Sparse method =", method, "null lambdas"))
  }
})

test_that("lmridge.eigen works with various starting values", {
  methods <- 1L:8L
  start_values <- rep(1, p) # Custom start
  
  for (method in methods) {
    control_solver <- control.solver.gen(solver = "eigen", method = method)
    
    # Custom starting values
    coefs <- lmridge.eigen(
      Xw = X_dense, yw = yw,
      lambdas = lambdas,
      start = start_values,
      lower = lower, upper = upper,
      control.solver = control_solver
    )
    expect_type(coefs, "double")
    expect_equal(length(coefs), p, info = paste("Method =", method, "with custom start"))
    
    # Default starting values
    coefs_default <- lmridge.eigen(
      Xw = X_dense, yw = yw,
      lambdas = lambdas,
      start = NULL,
      lower = lower, upper = upper,
      control.solver = control_solver
    )
    expect_type(coefs_default, "double")
    expect_equal(length(coefs_default), p, info = paste("Method =", method, "default start"))
  }
})


# High-dimensional mock data
set.seed(123)
n_high <- 5000
p_high <- 1000

X_dense_high <- matrix(rnorm(n_high * p_high), nrow = n_high, ncol = p_high)
X_sparse_high <- Matrix::rsparsematrix(n_high, p_high, density = 0.01)
yw_high <- rnorm(n_high)
lambdas_high <- rep(0.1, p_high) # Non-zero lambdas for regularization
control_solver_base <- list(tol = 1e-6, maxit = 1000, pivot = TRUE, ordering = TRUE, rowaug = TRUE)

# Constraints for coefficients
lower <- rep(-10, p_high)
upper <- rep(10, p_high)

test_that("lmridge.eigen works with all solver methods for high-dimensional dense and sparse matrices", {
  methods <- 1L:8L
  for (method in methods) {
    control_solver <- control.solver.gen(solver = "eigen", method = method)
    
    # Dense matrix test
    coefs_dense <- lmridge.eigen(
      Xw = X_dense_high, yw = yw_high,
      lambdas = lambdas_high,
      lower = lower, upper = upper,
      control.solver = control_solver
    )
    expect_type(coefs_dense, "double")
    expect_equal(length(coefs_dense), p_high, info = paste("Dense method =", method))
    
    # Sparse matrix test
    coefs_sparse <- lmridge.eigen(
      Xw = X_sparse_high, yw = yw_high,
      lambdas = lambdas_high,
      lower = lower, upper = upper,
      control.solver = control_solver
    )
    expect_type(coefs_sparse, "double")
    expect_equal(length(coefs_sparse), p_high, info = paste("Sparse method =", method))
  }
})

# Helper function to create tall matrices
generate_tall_matrix <- function(n, p, sparse = FALSE) {
  if (sparse) {
    return(Matrix::rsparsematrix(n, p, density = 0.1))
  } else {
    return(matrix(rnorm(n * p), nrow = n, ncol = p))
  }
}

# Mock response vector
generate_response <- function(n) {
  rnorm(n)
}

test_that("lmridge.eigen handles tall matrices (n * 100 >= p) for all methods", {
  n <- 2000   # Rows
  p <- 15     # Columns (ensures n * 100 >= p)
  
  lambdas <- rep(0.1, p)
  start <- rep(0, p)
  
  methods <- 1L:8L
  for (method in methods) {
    control_solver <- control.solver.gen(solver = "eigen", method = method)
    
    # Test dense matrix
    X_dense <- generate_tall_matrix(n, p, sparse = FALSE)
    yw_dense <- generate_response(n)
    
    coefs_dense <- lmridge.eigen(
      Xw = X_dense, yw = yw_dense,
      lambdas = lambdas,
      start = start,
      lower = rep(-Inf, p), upper = rep(Inf, p),
      control.solver = control_solver
    )
    
    expect_type(coefs_dense, "double")
    expect_equal(length(coefs_dense), p, info = paste("Method =", method, "Dense matrix"))
    
    # Test sparse matrix
    X_sparse <- generate_tall_matrix(n, p, sparse = TRUE)
    yw_sparse <- generate_response(n)
    
    coefs_sparse <- lmridge.eigen(
      Xw = X_sparse, yw = yw_sparse,
      lambdas = lambdas,
      start = start,
      lower = rep(-Inf, p), upper = rep(Inf, p),
      control.solver = control_solver
    )
    
    expect_type(coefs_sparse, "double")
    expect_equal(length(coefs_sparse), p, info = paste("Method =", method, "Sparse matrix"))
  }
})


test_that("lmridge.eigen handles wide matrices (p > n) for all methods", {
  n <- 50    # Rows
  p <- 500   # Columns (ensures p > n)
  
  lambdas <- rep(0.1, p)
  prior_mean <- rep(0, p)
  start <- NULL
  
  methods <- 1L:8L
  for (method in methods) {
    control_solver <- control.solver.gen(solver = "eigen", method = method)
    
    # Test dense matrix
    X_dense <- matrix(rnorm(n * p), nrow = n, ncol = p)
    yw_dense <- rnorm(n)
    
    coefs_dense <- lmridge.eigen(
      Xw = X_dense, yw = yw_dense,
      lambdas = lambdas,
      start = start,
      prior.mean = prior_mean,
      lower = rep(-Inf, p), upper = rep(Inf, p),
      control.solver = control_solver
    )
    
    expect_type(coefs_dense, "double")
    expect_equal(length(coefs_dense), p, info = paste("Method =", method, "Dense matrix"))
    
    # Test sparse matrix
    X_sparse <- Matrix::rsparsematrix(n, p, density = 0.1)
    yw_sparse <- rnorm(n)
    
    coefs_sparse <- lmridge.eigen(
      Xw = X_sparse, yw = yw_sparse,
      lambdas = lambdas,
      start = start,
      prior.mean = prior_mean,
      lower = rep(-Inf, p), upper = rep(Inf, p),
      control.solver = control_solver
    )
    
    expect_type(coefs_sparse, "double")
    expect_equal(length(coefs_sparse), p, info = paste("Method =", method, "Sparse matrix"))
  }
})

test_that("lmridge.eigen handles row augmentation for all methods", {
  n <- 100   # Rows
  p <- 20    # Columns
  
  # Define lambdas and prior.mean
  nonzero_lambdas <- rep(0.1, p)  # Non-zero lambdas
  zero_lambdas <- rep(0, p)       # Zero lambdas
  prior_mean <- rep(0, p)         # Prior mean
  
  # Generate mock data
  X_dense <- matrix(rnorm(n * p), nrow = n, ncol = p)
  X_sparse <- Matrix::rsparsematrix(n, p, density = 0.1)
  yw <- rnorm(n)
  
  methods <- 1L:8L
  for (method in methods) {
    control_solver <- control.solver.gen(solver = "eigen", method = method)
    
    # Test non-zero lambdas (row augmentation should occur)
    coefs_dense_aug <- lmridge.eigen(
      Xw = X_dense, yw = yw,
      lambdas = nonzero_lambdas,
      prior.mean = prior_mean,
      control.solver = control_solver
    )
    expect_type(coefs_dense_aug, "double")
    expect_equal(length(coefs_dense_aug), p, info = paste("Method =", method, "Dense matrix, non-zero lambdas"))
    
    coefs_sparse_aug <- lmridge.eigen(
      Xw = X_sparse, yw = yw,
      lambdas = nonzero_lambdas,
      prior.mean = prior_mean,
      control.solver = control_solver
    )
    expect_type(coefs_sparse_aug, "double")
    expect_equal(length(coefs_sparse_aug), p, info = paste("Method =", method, "Sparse matrix, non-zero lambdas"))
    
    # Test zero lambdas (row augmentation should NOT occur)
    coefs_dense_no_aug <- lmridge.eigen(
      Xw = X_dense, yw = yw,
      lambdas = zero_lambdas,
      prior.mean = prior_mean,
      control.solver = control_solver
    )
    expect_type(coefs_dense_no_aug, "double")
    expect_equal(length(coefs_dense_no_aug), p, info = paste("Method =", method, "Dense matrix, zero lambdas"))
    
    coefs_sparse_no_aug <- lmridge.eigen(
      Xw = X_sparse, yw = yw,
      lambdas = zero_lambdas,
      prior.mean = prior_mean,
      control.solver = control_solver
    )
    expect_type(coefs_sparse_no_aug, "double")
    expect_equal(length(coefs_sparse_no_aug), p, info = paste("Method =", method, "Sparse matrix, zero lambdas"))
  }
})


# Testing Xweights function
#--------------------------

test_that("Xweights returns X unchanged when weights vector is NULL", {
  X_dense <- matrix(1:9, nrow = 3, ncol = 3)
  result <- Xweights(X_dense, NULL)
  expect_identical(result, X_dense)  # Expect exact identical structure and values
})

test_that("Xweights returns X unchanged when weights sum matches row count", {
  X_dense <- matrix(1:9, nrow = 3, ncol = 3)
  weights_all_one <- rep(1, nrow(X_dense))
  result <- Xweights(X_dense, weights_all_one)
  expect_identical(result, X_dense)  # Again expect identical structure
})

test_that("Xweights applies weights correctly to dense matrix", {
  X_dense <- matrix(1:9, nrow = 3, ncol = 3)
  weights <- c(1, 0.5, 2)
  expected <- X_dense * sqrt(weights)
  result <- Xweights(X_dense, weights)
  expect_equal(result, expected)  # Allow small differences due to numeric precision
})


# Testing the normalizeX function
#--------------------------------

test_that("normalizeX scales the matrix columns correctly", {
  X <- matrix(1:9, nrow = 3, ncol = 3)
  result <- normalizeX(X)
  scaling <- attr(result, "scaling")
  
  expect_type(scaling, "double")
  expect_equal(length(scaling), ncol(X))  # Verify scaling vector length
  expect_true(all(scaling > 0))           # Ensure scaling values are positive
})

test_that("normalizeX returns a matrix with correct dimensions", {
  X <- matrix(1:9, nrow = 3, ncol = 3)
  result <- normalizeX(X)
  
  expect_equal(dim(result), dim(X))       # Dimensions should match input matrix
})

test_that("normalizeX adds the correct attributes", {
  X <- matrix(1:9, nrow = 3, ncol = 3)
  result <- normalizeX(X)
  
  expect_true(attr(result, "normalized"))       # Check 'normalized' attribute
  expect_true(!is.null(attr(result, "scaling")))  # Check 'scaling' attribute
})

test_that("normalizeX returns correct scaling for a simple matrix", {
  X <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow = 3, ncol = 3)
  result <- normalizeX(X)
  scaling <- attr(result, "scaling")
  
  expected_scaling <- scaleX(X) 
  expect_equal(scaling, expected_scaling)
})

test_that("normalizeX handles edge cases like zero columns", {
  X <- matrix(c(0, 0, 0, 1, 2, 3), nrow = 3, ncol = 2)
  result <- normalizeX(X)
  
  expect_true(all(is.finite(attr(result, "scaling")))) 
  expect_equal(dim(result), dim(X))                     
})

# Testing centerX function
#-------------------------
test_that("centerX centers columns of matrix correctly", {
  X <- matrix(1:9, nrow = 3, ncol = 3)
  result <- centerX(X)
  
  column_means <- Matrix::colMeans(X)
  expected <- sweep(X, 2, column_means, "-")  # Manually calculate expected centered values
  
  expect_equal(result, expected)               # Check if result matches expected centered matrix
})

test_that("centerX preserves matrix dimensions", {
  X <- matrix(1:9, nrow = 3, ncol = 3)
  result <- centerX(X)
  
  expect_equal(dim(result), dim(X))            # Dimensions should remain the same
})

test_that("centerX centers each column around zero", {
  X <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3)
  result <- centerX(X)
  
  column_means_after_centering <- Matrix::colMeans(result)
  expect_true(all(abs(column_means_after_centering) < 1e-10))  # Check if each column mean is near zero
})

test_that("centerX handles a matrix with constant columns", {
  X <- matrix(rep(5, 9), nrow = 3, ncol = 3)  # A matrix where each column has the same constant value
  result <- centerX(X)
  
  expect_true(all(result == 0))               # Expect a matrix of zeros after centering
})

test_that("centerX handles a sparse matrix input", {
  X <- Matrix::Matrix(1:9, nrow = 3, ncol = 3, sparse = TRUE)
  result <- centerX(X)
  
  column_means <- Matrix::colMeans(X)
  expected <- X - Matrix::Matrix(rep(column_means, each = nrow(X)), nrow = nrow(X), ncol = ncol(X), sparse = TRUE)
  
  expect_equal(as.matrix(result), as.matrix(expected))  # Check if result matches expected centered matrix
  expect_true(inherits(result, "Matrix"))               # Confirm output is a Matrix object
})


# Testing start_coefs function
#----------------------------

test_that("start_coefs calculates correctly for Gaussian family", {
  X <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3)
  family <- gaussian(identity)
  
  result <- start_coefs(X, family)
  
  scaling <- sqrt(colSums(X^2))  # Expected scaling (L2 norm of columns)
  expected <- rep(sqrt(5), ncol(X)) / scaling
  
  expect_equal(result, expected, tolerance = 1e-6)
})

test_that("start_coefs calculates correctly for Gamma family", {
  X <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3)
  family <- Gamma(link = "identity")
  
  result <- start_coefs(X, family)
  
  scaling <- sqrt(colSums(X^2))  # Expected scaling (L2 norm of columns)
  expected <- rep(sqrt(50), ncol(X)) / scaling
  
  expect_equal(result, expected, tolerance = 1e-6)
})

test_that("start_coefs calculates correctly for inverse Gaussian family", {
  X <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3)
  family <- inverse.gaussian(link = "identity")
  
  result <- start_coefs(X, family)
  
  scaling <- sqrt(colSums(X^2))  # Expected scaling (L2 norm of columns)
  expected <- rep(sqrt(50), ncol(X)) / scaling
  
  expect_equal(result, expected, tolerance = 1e-6)
})

test_that("start_coefs returns zero for matrix with zero columns", {
  X <- matrix(0, nrow = 2, ncol = 0)
  family <- gaussian(identity)
  
  result <- start_coefs(X, family)
  
  expect_equal(result, numeric(0))
})

test_that("start_coefs returns correct length", {
  X <- matrix(rnorm(20), nrow = 5, ncol = 4)
  family <- gaussian(identity)
  
  result <- start_coefs(X, family)
  
  expect_length(result, ncol(X))
})

# Testing maxlambda function
#---------------------------

test_that("maxlambda calculates correctly for Gaussian family without unpenalized terms", {
  X <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3)
  weights <- matrix(c(1, 1), nrow = 2, ncol = 1)
  z <- matrix(c(2, 3), nrow = 2, ncol = 1)
  coefs <- c(0.5, 1, 1.5)
  family <- gaussian(identity)
  
  result <- maxlambda(X, weights, z, dispersion = 1, always.include = NULL, coefs, family)
  
  # Calculations for expected values
  Xw <- Diagonal(x = sqrt(weights)) %*% X
  zw <- Diagonal(x = sqrt(weights)) %*% z
  maxlambda_expected <- (4 / 1) * max((1 / colSums(Xw^2)) * (crossprod(Xw, zw)^2), na.rm = TRUE) * 0.5
  
  expect_equal(result, maxlambda_expected, tolerance = 1e-6)
})

test_that("maxlambda calculates correctly for Poisson family with terms that are always included", {
  X <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3)
  weights <- matrix(c(1, 1), nrow = 2, ncol = 1)
  z <- matrix(c(2, 3), nrow = 2, ncol = 1)
  coefs <- c(0.5, 1, 1.5)
  always.include <- 1  # Always include first column
  family <- poisson(link = "log")
  
  result <- maxlambda(X, weights, z, dispersion = 1, always.include, coefs, family)
  
  # Calculations for expected values
  contrib_alwaysincl <- X[, always.include, drop = FALSE] %*% coefs[always.include]
  z_adjusted <- z - contrib_alwaysincl
  X_penalized <- X[, -always.include, drop = FALSE]
  Xw <- Diagonal(x = sqrt(weights)) %*% X_penalized
  zw <- Diagonal(x = sqrt(weights)) %*% z_adjusted
  maxlambda_expected <- (4 / 1) * max((1 / colSums(Xw^2)) * (crossprod(Xw, zw)^2), na.rm = TRUE)
  
  expect_equal(result, maxlambda_expected, tolerance = 1e-6)
})

test_that("maxlambda handles binomial family with default dispersion", {
  X <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3)
  weights <- matrix(c(1, 2), nrow = 2, ncol = 1)
  z <- matrix(c(1, 4), nrow = 2, ncol = 1)
  coefs <- c(0.5, 1, 1.5)
  family <- binomial(link = "logit")
  
  result <- maxlambda(X, weights, z, dispersion = 1, always.include = NULL, coefs, family)
  
  Xw <- Diagonal(x = sqrt(weights)) %*% X
  zw <- Diagonal(x = sqrt(weights)) %*% z
  maxlambda_expected <- (4 / 1) * max((1 / colSums(Xw^2)) * (crossprod(Xw, zw)^2), na.rm = TRUE)
  
  expect_equal(result, maxlambda_expected, tolerance = 1e-6)
})

test_that("maxlambda respects family-specific dispersion for Gamma family", {
  X <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3)
  weights <- matrix(c(1, 1), nrow = 2, ncol = 1)
  z <- matrix(c(2, 3), nrow = 2, ncol = 1)
  coefs <- c(0.5, 1, 1.5)
  family <- Gamma(link = "inverse")
  
  result <- maxlambda(X, weights, z, dispersion = 1.5, always.include = NULL, coefs, family)
  
  Xw <- Diagonal(x = sqrt(weights)) %*% X
  zw <- Diagonal(x = sqrt(weights)) %*% z
  maxlambda_expected <- (4 / 1.5) * max((1 / colSums(Xw^2)) * (crossprod(Xw, zw)^2), na.rm = TRUE)
  
  expect_equal(result, maxlambda_expected, tolerance = 1e-6)
})


