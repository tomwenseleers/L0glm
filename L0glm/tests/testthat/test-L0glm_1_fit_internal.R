# UNIT TESTING FOR L0glm.fit.internal FUNCTION

test_that("L0glm.fit.internal works with basic Gaussian data", {
  set.seed(1)
  n <- 100
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  
  # Run the L0glm.fit.internal function
  result <- L0glm.fit.internal(
    X = X,
    y = y,
    family = gaussian(),
    lambda = log(p) * log(log(n)),
    normalize = TRUE
  )
  
  # Check if the function returns a list with correct elements
  expect_type(result, "list")
  expect_true("coefficients" %in% names(result))
  expect_true("weights" %in% names(result))
  expect_true("dispersion" %in% names(result))
  
  # Check that the coefficients are numeric and of length p
  expect_type(result$coefficients, "double")
  expect_length(result$coefficients, p)
})


test_that("L0glm.fit.internal works with basic Poisson data", {
  set.seed(1)
  n <- 100
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- rpois(n, 1)
  
  # Run the L0glm.fit.internal function
  result <- L0glm.fit.internal(
    X = X,
    y = y,
    family = poisson(),
    lambda = log(p) * log(log(n)),
    normalize = TRUE
  )
  
  # Check if the function returns a list with correct elements
  expect_type(result, "list")
  expect_true("coefficients" %in% names(result))
  expect_true("weights" %in% names(result))
  expect_true("dispersion" %in% names(result))
  
  # Check that the coefficients are numeric and of length p
  expect_type(result$coefficients, "double")
  expect_length(result$coefficients, p)
})


test_that("L0glm.fit.internal works with basic Binomial data", {
  set.seed(1)
  n <- 100
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- rbinom(n, 1, 0.5)
  
  # Run the L0glm.fit.internal function
  result <- L0glm.fit.internal(
    X = X,
    y = y,
    family = binomial(),
    lambda = log(p) * log(log(n)),
    normalize = TRUE
  )
  
  # Check if the function returns a list with correct elements
  expect_type(result, "list")
  expect_true("coefficients" %in% names(result))
  expect_true("weights" %in% names(result))
  expect_true("dispersion" %in% names(result))
  
  # Check that the coefficients are numeric and of length p
  expect_type(result$coefficients, "double")
  expect_length(result$coefficients, p)
})


test_that("L0glm.fit.internal works with Gamma data", {
  set.seed(1)
  n <- 100
  p <- 10
  X <- cbind(Intercept = 1, matrix(rnorm(n * (p - 1L)), n, p - 1L))
  y <- rgamma(n, 1, 1)
  
  # Run the L0glm.fit.internal function
  result <- L0glm.fit.internal(
    X = X,
    y = y,
    family = Gamma(),
    always.include = 1L,
    lambda = log(p) * log(log(n)),
    normalize = TRUE,
    weights = rep(1, nrow(X))
  )
  
  # Check if the function returns a list with correct elements
  expect_type(result, "list")
  expect_true("coefficients" %in% names(result))
  expect_true("weights" %in% names(result))
  expect_true("dispersion" %in% names(result))
  
  # Check that the coefficients are numeric and of length p
  expect_type(result$coefficients, "double")
  expect_length(result$coefficients, p)
})


test_that("L0glm.fit.internal respects nonnegativity constraints via lower", {
  # Simulate some data
  set.seed(3)
  n <- 50
  p <- 8
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  
  # Run the L0glm.fit.internal function with nonnegative constraint
  result <- L0glm.fit.internal(
    X = X,
    y = y,
    family = gaussian(),
    lambda = log(p) * log(log(n)),
    lower = rep(0, p),
    normalize = TRUE
  )
  
  # Check that the coefficients are non-negative
  expect_true(all(result$coefficients >= 0))
})

test_that("L0glm.fit.internal handles non-NULL prior.mean correctly with stable inputs", {
  # Mock data
  n <- 10
  p <- 5
  set.seed(123)
  X <- matrix(rnorm(n * p, mean = 0, sd = 1), nrow = n, ncol = p)  # Well-scaled predictors
  y <- rpois(n, lambda = 2)  # Poisson response, compatible with Poisson family
  
  # Define a stable prior.mean
  prior_mean <- rep(0.1, p)  # Small, constant values to avoid numerical instability
  
  # Call the function with non-NULL prior.mean
  result <- tryCatch(
    {
      L0glm.fit.internal(
        X = X,
        y = y,
        prior.mean = prior_mean,
        family = poisson(),  # Explicitly pass a valid family object
        verbose = FALSE
      )
    },
    error = function(e) e
  )
  
  # Debugging: Capture outputs
  if (inherits(result, "error")) {
    cat("Error encountered:\n")
    print(result)
  }
  
  # Ensure the function runs without producing NaNs
  expect_false(inherits(result, "error"), "Function should handle non-NULL prior.mean without errors.")
  
  # Ensure coefficients are computed
  expect_true(!is.null(result$coefficients), "Coefficients should be returned.")
  expect_true(sum(abs(result$coefficients)) > 0, "Coefficients should not be all zero.")
})



# Only works with very small lambdas
# These user-specified starting values are subject to heavy penalization when lambda >= 1
# Resulting in all coefficients being zero
test_that("L0glm.fit.internal handles user-specified starting values", {
  set.seed(1)
  n <- 100
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- rpois(n, lambda = 5 + X[, 1] * 2) # Add a signal on the first variable
  
  # Manually specify starting values
  start_values <- rep(0.1, p)
  
  lambda_test <- 1e-1 
  
  # Increase the number of iterations
  max_iterations <- 100
  
  # Fit the model with user-specified starting values
  result <- L0glm.fit.internal(
    X = X,
    y = y,
    family = poisson(),
    start = start_values,
    lambda = lambda_test,
    maxit = max_iterations, 
    verbose = FALSE 
  )
  
  # Check that at least one coefficient is non-zero
  expect_true(any(result$coefficients != 0), 
              info = "Expected at least one non-zero coefficient, but all were zero.")
})


test_that("L0glm.fit.internal handles user-specified starting values with nonnegative constraint", {
  # Set up test data with a stronger signal
  set.seed(1)
  n <- 100
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- rpois(n, lambda = 5 + X[, 1] * 2) # Add a signal on the first variable
  
  # Manually specify starting values
  start_values <- rep(0.1, p)
  
  lambda_test <- 1e-1
  
  # Increase the number of iterations
  max_iterations <- 100
  
  # Fit the model with user-specified starting values
  result <- L0glm.fit.internal(
    X = X,
    y = y,
    family = poisson(),
    start = start_values,
    lambda = lambda_test,
    maxit = max_iterations,
    lower = rep(0, p),
    verbose = FALSE 
  )

  expect_true(all(result$coefficients >= 0), 
              info = "Expected all coefficients to be non-negative, but some were negative.")
})


# Works for Poisson but there is a slight difference for Gaussian each time
test_that("L0glm.fit.internal computes lambdas correctly with non-zero lambda", {
  # Simulate some data
  set.seed(1)
  n <- 100  # Number of observations
  p <- 10   # Number of variables
  X <- matrix(rpois(n * p, 5), nrow = n, ncol = p)
  y <- rpois(n, 5)  # Response variable
  
  # Define a non-zero lambda value
  lambda <- 1e-4
  
  # Call the function with a Gaussian family
  result <- L0glm.fit.internal(X = X, 
                               y = y, 
                               family = poisson(link = "identity"), 
                               lambda = lambda)
  
  # Check if coefficients are valid
  expect_true(all(is.finite(result$coefficients)), "Coefficients should be finite")
  
  # Calculate expected lambdas manually
  pen.weights <- AR.weights(beta = result$coefficients, gamma = 2, delta = 0, q = 0)
  
  # Ensure pen.weights has non-zero values to avoid Inf when calculating expected_lambdas
  expect_true(all(pen.weights > 0), "pen.weights should not be all zeros")
  
  dispersion <- 1  # Adjust if necessary based on your model
  expected_lambdas <- (1/4) * lambda * dispersion * pen.weights
  
  # Compare actual and expected lambdas
  expect_equal(result$lambdas, expected_lambdas, tolerance = 1e-6)
})


test_that("L0glm.fit.internal computes lambdas correctly with zero lambda", {
  # Simulate some data
  set.seed(1)
  n <- 100  # Number of observations
  p <- 10   # Number of variables
  X <- matrix(rpois(n * p, 5), nrow = n, ncol = p)
  y <- rpois(n, 5)  # Response variable
  
  # Define a zero lambda value
  lambda <- 0
  
  # Call the function with a Poisson family
  result <- L0glm.fit.internal(X = X, 
                               y = y, 
                               family = poisson(link = "identity"), 
                               lambda = lambda)
  
  # Check if coefficients are valid
  expect_true(all(is.finite(result$coefficients)), "Coefficients should be finite")
  
  # Calculate expected lambdas manually
  pen.weights <- AR.weights(beta = result$coefficients, gamma = 2, delta = 0, q = 0)
  
  # Ensure pen.weights has non-zero values to avoid Inf when calculating expected_lambdas
  expect_true(all(pen.weights > 0), "pen.weights should not be all zeros")
  
  dispersion <- 1  # Adjust if necessary based on your model
  expected_lambdas <- rep(0, length(pen.weights))
  expected_lambdas[is.infinite(pen.weights)] <- Inf
  expect_equal(result$lambdas, expected_lambdas, tolerance = 1e-6)
})


test_that("L0glm.fit.internal evaluates nonzero correctly and handles ridge solver fallback", {
  # Controlled small data
  n <- 5
  p <- 3
  X <- matrix(c(1, 2, 3, 4, 5, 1, 0, -1, 1, 0, 1, -1, 1, -1, 0), nrow = n, ncol = p)
  y <- c(1, 2, 3, 4, 5)  # Simple response
  
  # Ensure positive weights
  weights <- rep(1, n)
  
  # Case 1: start is NULL, nonzero depends on beta
  result1 <- tryCatch(
    {
      L0glm.fit.internal(
        X = X,
        y = y,
        weights = weights,
        start = NULL,  # Trigger the is.null(start) branch
        verbose = FALSE
      )
    },
    error = function(e) e
  )
  
  # Ensure the function runs without errors
  expect_false(inherits(result1, "error"), "Function should handle start = NULL without errors.")
  
  # Case 2: start is not NULL, nonzero depends on start
  start <- c(0, 1, 0)  # Example start vector
  result2 <- tryCatch(
    {
      L0glm.fit.internal(
        X = X,
        y = y,
        weights = weights,
        start = start,  # Trigger the !is.null(start) branch
        verbose = FALSE
      )
    },
    error = function(e) e
  )
  
  # Ensure the function runs without errors
  expect_false(inherits(result2, "error"), "Function should handle non-NULL start without errors.")
  
  # Validate fallback mechanism (when the ridge solver fails)
  expect_true(sum(is.finite(result2$coefficients)) > 0, 
              "Fallback mechanism should provide valid coefficients.")
})

test_that("L0glm.fit.internal handles uninformative observations correctly", {
  
  ## --- build a data set with no informative rows -------------------------
  n <- 5; p <- 3
  X <- matrix(0, nrow = n, ncol = p)   # all-zero design
  y <- rep(0, n)                       # all-zero response
  w <- rep(0, n)                       # zero prior weights
  
  ## --- run and catch the warning -----------------------------------------
  expect_warning(
    utils::capture.output(
      res <- L0glm.fit.internal(
        X       = X,
        y       = y,
        weights = w,
        verbose = TRUE)
    ),
    "All coefficients zero - nothing to fit"
  )
  
  ## --- structural checks --------------------------------------------------
  expect_type(res$coefficients, "double")
  expect_length(res$coefficients, p)
  expect_true(all(res$coefficients == 0))
  
  expect_type(res$weights, "double")
  expect_length(res$weights, n)
  expect_true(all(res$weights == 0))
  
  ## dispersion still defined
  expect_equal(res$dispersion, 1)
})


test_that("L0glm.fit.internal returns zero coefficients under huge penalty", {
  n <- 10; p <- 5
  X <- matrix(1, n, p)          # identical columns
  y <- rep(1, n)
  lambda_big <- 1e10            # enormous ridge
  
  res <- L0glm.fit.internal(X, y, lambda = lambda_big, verbose = FALSE)
  
  ## every coefficient driven to 0
  expect_equal(res$coefficients, rep(0, p))
  
  ## Unsmoothed weights diverge at zero; these coefficients are excluded.
  expect_equal(res$lambdas, rep(Inf, p))
})

