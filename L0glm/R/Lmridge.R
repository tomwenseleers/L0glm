#' Model solvers
#'
#' @import nnls
#' @import Matrix
#' @import glmnet
#' @importFrom osqp osqpSettings osqp solve_osqp
#' @importFrom stats lm.fit
#' @importFrom utils flush.console
#'
#' @description





#' @describeIn L0glm
#'
#' Implementation of different solvers
#' It is called by the upper code \code{L0glm.iwls}.
#'
#' @export
lmridge <- function(X, y, lambda,
                    weights = NULL,
                    lower = rep(-.Machine$double.xmax,ncol(X)),
                    upper = rep(.Machine$double.xmax,ncol(X)),
                    nonnegative = FALSE,
                    solver = "osqp",
                    start = rep(0, ncol(X)),
                    control.solver = NULL,
                    verbose){
  
  
  # If nonnegative is true and no lower bound is provided, or the lower bound is larger than 0, it is clipped to 0
  if (nonnegative && all(lower < 0)) lower <- rep(0, ncol(X))
  
  # Initialize control.solver variable
  if(is.null(control.solver)){
    control.solver <- list(solver = solver)
  }
  
  if (solver == "osqp") {
    
    control.solver <- do.call("control.solver.gen", control.solver)
    coefs <- lmridge_osqp(X = X, y = y, weights = weights, lambdas = lambda,
                          lower = lower, upper = upper, control.solver = control.solver)
    
  } else if (solver == "glmnet") {
    control.solver <- do.call("control.solver.gen", control.solver)
    coefs <- lmridge_glmnet(X = X, y = y, weights = weights, lambdas = lambda,
                            lower = lower, upper = upper, control.solver = control.solver)
    
    
  } else if (solver == "LSCG") {
    if (issparse) { 
      coefs <- solve_sparse_LSCG(A = Xaug, b = yaug, x = start)$coefficients 
    } else { 
      coefs <- solve_dense_LSCG(A = Xaug, b = yaug, x = start)$coefficients }
    coefs <- pmax(lower, pmin(coefs, upper)) # clip coefficients within allowed range
  } else if (solver == "QR") {
    if (issparse) { coefs <- solve_sparse_QR(A = Xaug, b = yaug)$coefficients } else
    { coefs <- solve_dense_QR(A = Xaug, b = yaug)$coefficients }
    coefs <- pmax(lower, pmin(coefs, upper)) # clip coefficients within allowed range
  } else if (solver == "nnls") {
    coefs <- nnls(A=Xaug, b=yaug)$x
  } else {
    stop("Invalid solver choice") }
  
  ###Old solvers###
  # Weighted LS update, with optional nonnegativity constraints
  # if(is.null(control.fit$block.size)) control.fit$block.size <- p
  # block.ids <- ceiling((1:p)/control.fit$block.size)
  # beta[!converged_set] <- block.fit(X = Xaug[good,!converged_set,drop=F]*sqrt(W[good]),
  #                                   y = zaug[good]*sqrt(W[good]),
  #                                   block.id = block.ids[!converged_set],
  #                                   nonnegative = nonnegative,
  #                                   control = control.fit)
  
  
  # TODO clean this
  # if(F){
  #   if(nonnegative){
  #     beta[!converged_set] <- nnls(A = Xaug[good,!converged_set,drop=F]*sqrt(W[good]),
  #                                  b = zaug[good]*sqrt(W[good]))$x
  #   } else {
  #     beta[!converged_set] <- lm.fit(x = as.matrix(Xaug[good,!converged_set,drop=F]*sqrt(W[good])),
  #                                    y = zaug[good]*sqrt(W[good]),
  #                                    tol = 1E-07)$coefficients
  #     beta[is.na(beta) | is.infinite(beta)] <- 0
  #   }
  # }
  return(coefs)
}


lmridge_osqp <- function(X, y, weights, lambdas, lower, upper, control.solver){
  
  Xw <- Xweights(X, weights)
  yw <- Yweights(y, weights)
  
  # XtX = crossprod(Xw) + diaglambdas(lambdas, issparse(Xw))
  # Xty = crossprod(Xw, yw)
  
  Xwaug <- row_augX(Xw, lambda)
  ywaug <- row_augY(yw, lambda)

  XtX = crossprod(Xwaug)
  Xty = crossprod(Xwaug, ywaug)
  
  #Obtain settings  and create model
  settings <- do.call(osqpSettings, control.solver)
  if (sum(start)!=0) settings <- c(settings, warm_start = TRUE) else settings <- c(settings, warm_start = FALSE)
  
  model <- osqp(P = XtX, q = -Xty, A = .sparseDiagonal(p), l = lower, u = upper, pars = settings)
  
  if (!is.null(start)) model$WarmStart(x=start)
  
  #Solve model
  coefs <- model$Solve()$x
  coefs <- pmax(lower, pmin(coefs, upper)) # fitted coefficients sometimes go very slightly outside constraint zone due to numerical inaccuracies in solver - this is fixed here via clipping
}

lmridge_glmnet <- function(X, y, weights, lambdas, lower, upper, control.solver){
  # Set path depending on dimensions of X
  if(sum(lambda) == 0 && ncol(X) >= (nrow(X)/2)){ # avoid numerical instabilities if lambda is 0 and n is not much larger than p
    path = TRUE
  } else {
    path = FALSE
  }
  
  # Obtain Xwaug and ywaug
  Xw <- Xweights(X, weights)
  yw <- Yweights(y, weights)
  
  Xwaug <- row_augX(Xw, lambda)
  ywaug <- row_augY(yw, lambda)
  
  model <- bigGlm(x = Xwaug, y = ywaug,
                  intercept = FALSE, standardize = FALSE, standardize.response = FALSE,
                  lower.limits = c(-Inf, lower), upper.limits = c(Inf, upper),
                  family = "gaussian", path = FALSE, thresh = control.solver$thresh,
                  maxit = control.solver$maxit)
  coefs <- coef(model)[-1]
}

lmridge_eigen <- function(X, y, weights, lambdas, lower, upper, control.solver){
  
}


Xweights <- function(X, weights){
  if(!is.null(weights)){
    Xw <- X * sqrt(weights)
  } else {
    Xw <- X
  }
  
  return(Xw)
}

Yweights <- function(y, weights){
  if (!is.null(weights)) {
    yw <- y * sqrt(weights)
  } else {
    yw <- y
  }
  return(yw)
}

# Row augmentation of the covariate matrix if lambda!=0 to add ridge penalty
# Ref :
# https://stats.stackexchange.com/a/164546/27340
# https://stats.stackexchange.com/a/203699/27340
row_augX <- function(Xw, lambda){
  p <- ncol(Xw)
  lambda <- sqrt(lambda)
  if (sum(lambda) != 0) {
    Xaug <- rbind(Xw, diaglambdas(lambda, issparse(Xw)))
  } else {
    Xaug <- Xw
  }
  return(Xaug)
}

row_augY <- function(yw, lambda){
  p <- ncol(Xw)
  
  if (sum(lambda) != 0) {
    yaug <- c(yw, rep(0, p))
  } else {
    yaug <- yw
  }
  return(yaug)
}

diaglambdas <- function(lambdas, sparse){
  if(sparse) {
    dlambdas <- .sparseDiagonal(p, lambdas)
  } else {
    dlambdas <- diag(lambdas, p)
  }
  return (dlambdas)
}
issparse <- function(X){
  return(any(class(X) %in% c("dgCMatrix", "dsCMatrix")))
}
