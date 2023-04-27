
# TODO documentation
group.cpglm <- function(X, Y,
                            Weights = matrix(1, nrow = nrow(Y), ncol = ncol(Y)), # TODO for efficiency maybe use a single vector
                            family = poisson(identity),
                            intercept = FALSE,
                            Start = NULL,
                            constr = "nonneg",
                            lambdas = 0,
                            group.pen.fun = function(mat) sqrt(rowSums(mat^2)), # default = L2 norm of the estimated coefficients per component
                            no.pen = 0,
                            level = 0.95,
                            maxit = 25,
                            rel.tol = 1E-4,
                            thresh = 1E-6,
                            delta = 1E-5,
                            cv = "none",
                            verbose = TRUE){
  n <- nrow(Y) # nr of obersvations
  m <- ncol(Y) # nr of correlated channels
  p <- ncol(X) # nr of covariates

  # Check grouping penalty function
  test.fun <- group.pen.fun(matrix(rep(c(-1, 0, 1), nrow = 3, ncol = 5)))
  if(any(is.na(test.fun) | is.infinite(test.fun)) || !is.vector(test.fun) || length(test.fun) != 3)
    stop("Invalid function supplied. \nThe function should take a matrix and return a vector of length \'ncol(X)\' with finite and non-NA elements")

  # Check family argument
  if (is.character(family)) family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) family <- family()
  if (family$link!="identity" | !(family$family %in% c("gaussian","poisson","Gamma","inverse.gaussian","quasipoisson")) )
    stop("family(link) must be one of: gaussian(identity), poisson(identity), Gamma(identity), inverse.gaussian(identity), quasipoisson(identity)")
  # TODO: support any link function? also support binomial & quasibinomial

  # Add an intercept if required
  if(intercept){
    X <- cbind(Intercept = 1, X)
    no.pen <- unique(c(1, no.pen + 1)) # No penalization of the intercept
  }

  # Get penalization factors lambdas
  if(identical(lambdas, 0)) lambdas <- rep(0, m)
  if(any(lambdas < 0)) stop("Negative lambdas not allowed")
  if(cv == "none"){
    if(length(lambdas) != m) stop("You must provide a lambda for every channel, ie length(lambdas) must equal ncol(Y), or choose a cross-validation method")
  } else {
    # TODO
    stop("Not implemented")
  }
  lambdas <- rep(lambdas, p + intercept)
  lambdas[no.pen] <- 0
  # Fit the L0 penalized nonnegative GLM
  fit <- group.cpglm.fit(X = X, Y = Y, Weights = Weights, family = family,
                             lambdas = lambdas, group.pen.fun = group.pen.fun,
                             Start = Start, maxit = maxit, rel.tol = rel.tol,
                             constr = constr, thresh = thresh, delta = delta)

  # TODO compute glm slots: residuals, deviance, etc.
  # TODO compute effetive degrees of freedom
  # TODO perform inference on estimated coefficients

  names(fit)[names(fit) == "mu"] <- "fitted.values"
  return(structure(fit, class = c("cpglm"))) # TODO glm and lm compatibility? Maybe not possible because of multivariate Y

}



group.cpglm.fit <- function(X = X, Y = Y,
                                Weights = Weights,
                                family = family,
                                lambdas = rep(0, ncol(X)),
                                group.pen.fun = function(mat) sqrt(Matrix::rowSums(mat^2)), # default = L2 norm of the estimated coefficients per component
                                Start = NULL,
                                constr = "nonneg",
                                maxit = 25,
                                rel.tol = 1E-4, # Needed ?
                                thresh = 1E-6,
                                gamma = 2,
                                delta = 1E-5){
  n <- nrow(Y) # nr of obersvations
  m <- ncol(Y) # nr of correlated channels
  p <- ncol(X) # nr of covariates


  # The adaptive ridge loop
  Beta <- if(is.null(Start)) matrix(1, nrow = p, ncol = m) else Start
  Betas <- list() # list with coefficient estimates over successive iterations
  for (i in 1:maxit) {

    # Beta0 <- Beta # coefs from previous iteration
    Beta.group <- group.pen.fun(Beta)
    lambda_adap <- lambdas %o% AR.weights(beta = Beta.group, gamma = gamma, delta = delta) # ridge penalty matrix with adaptive weights

    # Fit the model with current adaptive weights
    Beta <- do.call(cbind, lapply(1:m, function(channel){
      # beta <- cpglm.iwls(X = X, y = Y[,channel], weights = Weights[,channel],
      #                  family = family, start = Start[channel,], lambda = lambda_adap[channel,],
      # beta <- block.fit(X = rbind(X*sqrt(Weights[,channel]), diag(sqrt(lambda_adap[channel,]))),
      #                   y = c(Y[,channel]*sqrt(Weights[,channel]), rep(0, p)),
      #                   control.fit = control.fit) # TODO put block size as an argument
      # microbenchmark(
      #   beta <- nnls(A = rbind(X*sqrt(Weights[,channel]), diag(sqrt(lambda_adap[channel,]))),
      #                b = c(Y[,channel]*sqrt(Weights[,channel]), rep(0, p)))$x
      # )
      nnls(A = rbind(X * sqrt(Weights[,channel]), diag(sqrt(lambda_adap[channel,]))),
           b = c(Y[,channel]*sqrt(Weights[,channel]), rep(0, p)))$x
      # TODO Tom's code is slightly faster (for computing the row augmented data), adapt to increase speed
    }))
    Start <- Beta

    # TODO needed ?
    # Betas <- c(Betas, list(Beta)) #this can become memory inefficient ! # store the new coefficients

    # Control iteration
    # break when the non zero subset doesn't change between consecutive iterations
    # stop.loop <- any(Matrix::rowSums(Beta) == 0) & all(!xor(Matrix::rowSums(Beta) == 0, Matrix::rowSums(Beta0) == 0)) # TODO improve this ?
    # if(stop.loop) break



  }
  return(list(coefficients = Beta,
              mu = X %*% Beta,
              prior.weights = Weights,
              # TODO weights = fit$weights,
              # TODO deviance = fit$deviance,
              family = family,
              # TODO lambda_adap = lambda_adap,
              # TODO converged = stop.loop
              niter = i))
}