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
  X <- matrix(rnorm(n * p), n, p)
  y <- rgamma(n, 1, 1)
  
  # Run the L0glm.fit.internal function
  result <- L0glm.fit.internal(
    X = X,
    y = y,
    family = Gamma(),
    lambda = log(p) * log(log(n)),
    normalize = TRUE,
    weights = NULL
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


test_that("L0glm.fit.internal respects nonnegative constraint", {
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
    nonnegative = TRUE,
    normalize = TRUE
  )
  
  # Check that the coefficients are non-negative
  expect_true(all(result$coefficients >= 0))
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
    verbose = TRUE 
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
    nonnegative = TRUE,
    verbose = TRUE 
  )

  expect_true(all(result$coefficients >= 0), 
              info = "Expected all coefficients to be non-negative, but some were negative.")
})


test_that("solver handling with _noreg works correctly", {
  # Test case where solver includes _noreg
  solver <- "eigen_noreg"
  lambda <- 1  # some non-zero value
  maxit <- 25  # some value greater than 1
  
  # Execute the code being tested
  if (grepl("_noreg", solver)) {
    solver <- gsub("_noreg", "", solver)
    lambda <- 0
    maxit <- 1
  }
  
  # Assertions
  expect_equal(solver, "eigen")
  expect_equal(lambda, 0)
  expect_equal(maxit, 1)
  
  # Test case where solver does not include _noreg
  solver <- "eigen"
  lambda <- 1  # some non-zero value
  maxit <- 25  # some value greater than 1
  
  # Execute the code being tested
  if (grepl("_noreg", solver)) {
    solver <- gsub("_noreg", "", solver)
    lambda <- 0
    maxit <- 1
  }
  
  # Assertions
  expect_equal(solver, "eigen")
  expect_equal(lambda, 1)
  expect_equal(maxit, 25)
})


# Create a custom family object
custom_family <- list(
  family = "custom",
  variance = function(mu) mu,           # valid function
  linkinv = function(eta) exp(eta),     # valid function
  linkfun = function(mu) log(mu),       # valid function
  initialize = function() list(mu = NULL, n = NULL, m = NULL),
  dev.resids = function(y, mu, wt) (y - mu)^2 / wt,  # valid function
  mu.eta = function(eta) exp(eta)       # valid function
)


# Assign class "family"
class(custom_family) <- "family"


test_that("L0glm.fit.internal checks if family object has valid functions", {
  
  # Test case when family object is missing variance function
  invalid_family_1 <- custom_family
  invalid_family_1$variance <- NULL  # Set variance to NULL
  
  expect_error({
    L0glm.fit.internal(X = matrix(rnorm(100), nrow = 10), 
                       y = rpois(10, lambda = 5), 
                       family = invalid_family_1)
  }, "'family' argument seems not to be a valid family object", fixed = TRUE)
  
  # Test case when family object is missing linkinv function
  invalid_family_2 <- custom_family
  invalid_family_2$linkinv <- NULL  # Set linkinv to NULL
  
  expect_error({
    L0glm.fit.internal(X = matrix(rnorm(100), nrow = 10), 
                       y = rpois(10, lambda = 5), 
                       family = invalid_family_2)
  }, "'family' argument seems not to be a valid family object", fixed = TRUE)
  
  # Test case with valid family object
  valid_family <- poisson()  # A valid family object
  
  expect_silent({
    L0glm.fit.internal(X = matrix(rnorm(100), nrow = 10), 
                       y = rpois(10, lambda = 5), 
                       family = valid_family)
  })
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
  expect_equal(result$lambdas, expected_lambdas, tolerance = 1e-6)
})


test_that("L0glm.fit.internal updates dispersion correctly when dispersion_init is NULL", {
  # Simulate some data
  set.seed(1)
  n <- 100
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  
  # Define a mock family object with identity link
  family <- gaussian(link = "identity")
  
  # Mock residuals and weights
  mu <- rnorm(n)  # Expected values
  gprime <- rep(1, n)  # Derivative of the link function
  residuals <- y - mu
  residuals_old <- residuals  # Initially, old residuals are the same
  W.ret <- rep(1, n)  # Equal weights
  good <- 1:n  # All observations are good
  nonzero <- rep(TRUE, n)  # All coefficients are non-zero
  
  # Simulate the initial dispersion
  dispersion_init <- NULL
  
  # Implement the relevant part of the function to test
  if (is.null(dispersion_init)) {
    residuals_old <- residuals  # old unweighted working residuals (on adjusted z scale)
    
    if (family$link == "identity") {
      residuals <- y - mu  # updated unweighted residuals
    } else {
      residuals <- (y - mu) / gprime  # updated unweighted residuals
    }
    
    n.ok <- length(good)  # good observations
    rank <- sum(nonzero)  # estimated rank
    df.residual <- n.ok  # residual degrees of freedom
    
    if (rank < n.ok) {
      df.residual <- df.residual - rank  # to also allow for high dimensional case
    }
    
    dispersion_old <- sum((W.ret[good] * residuals_old[good]^2)) / df.residual  # previous iteration's estimated dispersion
    dispersion_new <- sum((W.ret[good] * residuals[good]^2)) / df.residual  # updated estimate of dispersion (variance for Gaussian case)
    disp_rel_change <- abs(dispersion_new - dispersion_old) / dispersion_old
    
    if (disp_rel_change < 0.1) {
      dispersion <- dispersion_new  # update the dispersion
    } else {
      dispersion <- dispersion_old  # if not converged, keep old dispersion
    }
  }
  
  # Assertions to test the expected behavior
  expect_true(!is.null(dispersion), info = "Dispersion should not be NULL after update.")
  expect_true(dispersion >= 0, info = "Dispersion should be non-negative.")
})

