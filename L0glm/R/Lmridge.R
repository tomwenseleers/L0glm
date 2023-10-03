#' Fit Constrained Ridge Regularised Weighted Least Squares Regression using chosen Solver
#'
#' @name lmridge
#' @aliases lmridge
#' @title Fits a ridge regularised weighted least squares regression using a
#'  chosen solver subject to optional nonnegativity or box constraints on the coefficients.
#'
#' @import Matrix
#' @importFrom nnls nnls
#' @import glmnet
#' @importFrom osqp osqpSettings osqp solve_osqp
#' @importFrom qlcMatrix colMax colMin
#'
#' @description
#' Implementation of different solvers
#' It is called by the upper \code{L0glm.fit.internal} code.
#' Note that
#' \code{lambda} should be a vector of length \code{p}, specifying
#' the ridge penalties on each of the coefficients.
#' 
#' @param X 
#' A matrix of predictor variables.
#' @param y 
#' A vector or matrix of response variables.
#' @param lambdas XXX
#' @param weights XXX
#' @param nonnegative XXX
#' @param lower XXX
#' @param upper XXX
#' @param scaling XXX
#' @param prior.mean XXX
#' @param solver XXX
#' @param method XXX
#' @param start XXX
#' @param control.solver XXX
#' @param verbose XXX
#' 
#' @return XXX
#' 
#' @export

# TODO finish documentation for this function

lmridge <- function(X, y, 
                    lambdas = NULL,
                    weights = NULL,
                    lower = rep(-.Machine$double.xmax,ncol(X)),
                    upper = rep(.Machine$double.xmax,ncol(X)),
                    nonnegative = FALSE,
                    scaling = rep(1, ncol(X)), # to allow scaling of columns of X, eg by L2 norm
                    prior.mean = NULL,
                    solver = "eigen",
                    method = ifelse(ncol(X)>=20000, 3L, 7L),
                    start = rep(1, ncol(X)),
                    control.solver = do.call("control.solver.gen", list(solver = solver, method = method)),
                    verbose = FALSE){
  # scale columns of X if desired
  if (sum(scaling)!=length(scaling)) { X <- Matrix(X %*% Matrix::Diagonal(x = 1 / scaling), sparse=TRUE)
  start <- start * scaling }
  
  # If nonnegative is true and no lower bound is provided, or the lower bound is larger than 0, it is clipped to 0
  if(nonnegative) lower[lower < 0] <- 0
  Xw <- Xweights(X, weights)
  yw <- Yweights(y, weights)
  if (!issparse(Xw)) { Xw = as.matrix(Xw) }
  yw[is.nan(yw)] <- 0
  
  if(is.null(prior.mean)){
    prior.mean = rep(0, ncol(X))
  } else if (length(as.vector(prior.mean)) == 1){
    prior.mean = rep(prior.mean, ncol(X))
  } else if (length(as.vector(prior.mean)) != ncol(X)) {
    stop("Prior.mean needs to be a vector of length p or a numerical.")
  }
  
  if (solver == "osqp") {
    coefs <- lmridge.osqp(Xw = Xw, yw = yw, lambdas = lambdas, start = start, 
                          prior.mean = prior.mean, 
                          lower = lower, upper = upper, 
                          control.solver = control.solver)
    
  } else if (solver == "glmnet") {
    coefs <- lmridge.glmnet(Xw = Xw, yw = yw, lambdas = lambdas, prior.mean = prior.mean,
                            lower = lower, upper = upper, control.solver = control.solver)
    
    
  } else if (solver == "eigen") {
    coefs <- lmridge.eigen(Xw = Xw, yw = yw, 
                           lambdas = lambdas, 
                           start = start, prior.mean = prior.mean,
                           lower = lower, upper = upper, control.solver = control.solver,
                           verbose = verbose)
    
  } else if (solver == "nnls") {
    coefs <- lmridge.nnls(Xw = Xw, yw = yw, lambdas = lambdas, prior.mean = prior.mean,
                          lower = lower, upper = upper, control.solver = control.solver)
  } else {
    stop("Invalid solver choice") }
  
  # rescale coefficients to original scale
  if (sum(scaling)!=length(scaling)) { coefs = coefs / scaling }
  
  if (sum(is.na(coefs))==length(coefs)) stop("Ridge solver failed.")
  
  return(coefs)
}



lmridge.osqp <- function(Xw, 
                         yw, 
                         lambdas = NULL, 
                         lower, upper, start, 
                         control.solver = do.call("control.solver.gen", list(solver = "eigen", method = 7L)), 
                         prior.mean = NULL){
  Xw = as(Matrix(Xw, sparse=T), "dgCMatrix")
  if(!is.null(lambdas)&(sum(lambdas)!=0)) {
    Xwaug <- row_augX(Xw, lambdas)
    ywaug <- augY(yw, ncol(Xwaug), lambdas = lambdas, prior.mean = prior.mean) 
  } else {
    Xwaug <- Xw
    ywaug <- yw
  }
  
  if(control.solver$rowaug == TRUE){
    XtX <- crossprod(Xwaug)
    Xty <- crossprod(Xwaug, ywaug)
  } else {
    if (issparse(Xw)){
      XtX <- as(Matrix(crossprod(Xwaug), sparse=T),"dgCMatrix")
    } else {
      XtX <- crossprod(Xwaug)
    }
    Xty <- as.matrix(crossprod(Xwaug, 
                               augY(yw, ncol(Xw), 
                                    lambdas = lambdas, 
                                    prior.mean = prior.mean)))
  }
  
  # Obtain settings  and create model
  settings <- do.call(osqp::osqpSettings, control.solver[-7])
  if (sum(start)!=0) settings <- c(settings, warm_start = TRUE) else 
    settings <- c(settings, warm_start = FALSE)
  
  model <- osqp::osqp(P = XtX, q = -Xty, 
                A = .sparseDiagonal(ncol(Xw)), # Diagonal(n = ncol(Xw)),
                l = lower, u = upper, pars = settings)
  if(is.null(start)) start <- rep(1, ncol(Xw))
  model$WarmStart(x=start)
  
  #Solve model
  coefs <- model$Solve()$x
  coefs <- clip(coefs, lower, upper, reflect_lower=F, reflect_upper=F) # fitted coefficients sometimes go very slightly outside constraint zone due to numerical inaccuracies in solver - this is fixed here via clipping
  return(coefs)
}


lmridge.glmnet <- function(Xw, yw, lambdas = NULL, 
                           lower, upper, 
                           control.solver, 
                           prior.mean = NULL){
  # Set path depending on dimensions of Xw
  if(sum(lambdas) == 0 && ncol(Xw) >= nrow(Xw)){ # avoid numerical instabilities if lambdas is 0 and p>=n
    path = TRUE
  } else {
    path = FALSE
  }
  
  # Obtain Xwaug and ywaug
  if(!is.null(lambdas)&(sum(lambdas)!=0)) {
    Xwaug <- row_augX(Xw, lambdas)
    ywaug <- augY(yw, ncol(Xwaug), 
                  lambdas = lambdas, 
                  prior.mean = prior.mean)
  } else {
    Xwaug <- Xw
    ywaug <- yw
  }
  
  model <- bigGlm(x = Xwaug, 
                  y = ywaug,
                  intercept = FALSE, 
                  standardize = control.solver$standardize, 
                  standardize.response = FALSE,
                  lower.limits = c(-Inf, lower), 
                  upper.limits = c(Inf, upper),
                  family = "gaussian", 
                  path = FALSE, 
                  thresh = control.solver$thresh,
                  maxit = control.solver$maxit)
  coefs <- coef(model)[-1]
  coefs <- clip(coefs, lower, upper, reflect_lower=F, reflect_upper=F) # check
  return(coefs)
}



lmridge.nnls <- function(Xw, yw, lambdas = NULL, lower, upper, 
                         control.solver, 
                         prior.mean = NULL){
  # Obtain Xwaug and ywaug
  if(!is.null(lambdas)&(sum(lambdas)!=0)) {
    Xwaug <- row_augX(Xw, lambdas)
    ywaug <- augY(yw, ncol(Xwaug), lambdas = lambdas, prior.mean = prior.mean)
  } else {
    Xwaug <- Xw
    ywaug <- yw
  }
  
  coefs <- nnls::nnls(A=as.matrix(Xwaug), b=ywaug)$x
  return(coefs)
}

lmridge.eigen <- function(Xw, yw, 
                          lambdas = rep(1, ncol(Xw)), 
                          lower = rep(-.Machine$double.xmax,ncol(Xw)), 
                          upper = rep(.Machine$double.xmax,ncol(Xw)), 
                          start = rep(0, ncol(Xw)), 
                          control.solver = do.call("control.solver.gen", 
                                                   list(solver = "eigen", method = 7L)), 
                          prior.mean = rep(0, ncol(Xw)),
                          verbose = FALSE){ 
  lambda = mean(lambdas)
  pen.weights = lambdas / lambda 
  if (!issparse(Xw)) { Xw = as.matrix(Xw) } else { Xw = as(Xw, "dgCMatrix") }
  # Xw = as(Matrix(Xw, sparse=T), "dgCMatrix") # test - always coerce X to sparse matrix
  
  if((ncol(Xw) > nrow(Xw))&control.solver$rowaug) { # &(control.solver$method != "LSCG")
    control.solver$rowaug <- FALSE
    if(verbose&control.solver$method != "LSCG") print("Row augmentation is not recommended for p > n with this solver. Using Woodbury identity optimized normal equations appropriate for p > n case.")
  }
  
  if (control.solver$method == "QR") {
    if (control.solver$rowaug) {
      # Obtain Xwaug and ywaug
      if (!is.null(lambdas)&(sum(lambdas)!=0)) {
        Xwaug <- row_augX(Xw, lambdas)
        ywaug <- augY(yw, ncol(Xwaug), lambdas = lambdas, prior.mean = prior.mean)
      } else {
        Xwaug <- Xw
        ywaug <- yw
      }
      
      if (issparse(Xwaug)){
        coefs <- solve_sparse_QR(A = Xwaug,
                                 b = ywaug)$coefficients
      } else {
        coefs <- solve_dense_QR(A = Xwaug,
                                b = ywaug, 
                                pivot = control.solver$pivot)$coefficients
      }
    } else {
      D <- Diagonal(x=1/pen.weights)
      if(nrow(Xw) >= ncol(Xw)){
        if(issparse(Xw)){
          coefs <- as.vector(solve_sparse_QR(A = D %*% Matrix::crossprod(Xw) + Diagonal(x=rep(lambda, ncol(Xw))),
                                             b = as.matrix(D %*% Matrix::crossprod(Xw,yw) + lambda*Matrix::crossprod(D, prior.mean)))$coefficients)
        } else { # if it's dense
          coefs <- as.vector(solve_dense_QR(A = as.matrix(D %*% Matrix::crossprod(Xw) + Diagonal(x=rep(lambda, ncol(Xw)))),
                                            b = as.matrix(D %*% Matrix::crossprod(Xw,yw) + lambda*Matrix::crossprod(D, prior.mean)),
                                            pivot = control.solver$pivot)$coefficients)
        }
      } else { # if p > n
        if(issparse(Xw)){
          # Precompute some values
          DXw <- Matrix::tcrossprod(D, Xw)
          y_star <- as.vector(yw - Xw %*% prior.mean)
          
          # Calculate β using Woodbury identity
          coefs <- as.vector(prior.mean + DXw %*% solve_sparse_QR(A = Xw %*% DXw + Diagonal(x=rep(lambda, nrow(Xw))),
                                                                  b = y_star)$coefficients) 
        } else { # if it's dense
          # Precompute some values
          DXw <- Matrix::tcrossprod(D, Xw)
          y_star <- as.vector(yw - Xw %*% prior.mean)
          
          # Calculate β using Woodbury identity
          coefs <- as.vector(prior.mean + DXw %*% solve_dense_QR(A = as.matrix(Xw %*% DXw + Diagonal(x=rep(lambda, nrow(Xw)))),
                                                                 b = y_star,
                                                                 pivot = control.solver$pivot)$coefficients)
        }
      }
    }
  } else if (control.solver$method == "LSCG") {
    if (control.solver$rowaug) {
      # Obtain Xwaug and ywaug
      if(!is.null(lambdas)) {
        Xwaug <- row_augX(Xw, lambdas)
        ywaug <- augY(yw, ncol(Xwaug), lambdas = lambdas, prior.mean = prior.mean)
      } else {
        Xwaug <- Xw
        ywaug <- yw
      }
      if(is.null(start)) start <- rep(0, ncol(Xw))
      # X.n <- sqrt(colSums(Xw^2)) 
      if(issparse(Xwaug)){
        coefs <- solve_sparse_LSCG(A = Xwaug,
                                   b = as.matrix(ywaug),
                                   x0 = as.matrix(start),
                                   tol = control.solver$tol,
                                   maxit = control.solver$maxit,
                                   precond = control.solver$precond)$coefficients 
      } else {
        coefs <- solve_dense_LSCG(A = Xwaug,
                                  b = as.matrix(ywaug),
                                  x0 = as.matrix(start),
                                  tol = control.solver$tol,
                                  maxit = control.solver$maxit,
                                  precond = control.solver$precond)$coefficients
      }
    } else { # Normal equations
      D <- Diagonal(x=1/pen.weights)
      
      if(nrow(Xw) >= ncol(Xw)){
        if(issparse(Xw)){
          coefs <- as.vector(solve_sparse_LSCG(A = D %*% Matrix::crossprod(Xw) + Diagonal(x=rep(lambda, ncol(Xw))),
                                               b = as.matrix(D %*% Matrix::crossprod(Xw,yw) + lambda*Matrix::crossprod(D, prior.mean)),
                                               x0 = as.matrix(start),
                                               tol = control.solver$tol,
                                               precond = control.solver$precond)$coefficients)
        } else { # if it's dense
          coefs <- as.vector(solve_dense_LSCG(A = as.matrix(D %*% Matrix::crossprod(Xw) + Diagonal(x=rep(lambda, ncol(Xw)))),
                                              b = as.matrix(D %*% Matrix::crossprod(Xw,yw) + lambda*Matrix::crossprod(D, prior.mean)),
                                              x0 = as.matrix(start),
                                              tol = control.solver$tol,
                                              maxit = control.solver$maxit,
                                              precond = control.solver$precond)$coefficients)
        }
      } else { # if p > n
        DXw <- Matrix::tcrossprod(D, Xw)
        if(issparse(Xw)){
          # Precompute some values
          y_star <- as.vector(yw - Xw %*% prior.mean)

          # Calculate β using Woodbury identity
          coefs <- as.vector(prior.mean + DXw %*% Matrix::solve(Xw %*% DXw + Diagonal(x=rep(lambda, nrow(Xw))),
                                                                y_star,
                                                                tol = control.solver$tol))
          # TODO; replace Matrix::solve with solve_sparse_LSCG
        } else { # if it's dense
          # Precompute some values
          y_star <- as.vector(yw - Xw %*% prior.mean)
          
          # Calculate β using Woodbury identity
          coefs <- as.vector(prior.mean + DXw %*% Matrix::solve(Xw %*% DXw + Diagonal(x=rep(lambda, nrow(Xw))),
                                                                y_star,
                                                                tol = control.solver$tol))
          # TODO: replace Matrix::solve with solve_dense_LSCG
        }
      }
    }
  } else if (control.solver$method == "CG") {
    if (control.solver$rowaug) {
      # Obtain Xwaug and ywaug
      if(!is.null(lambdas)&(sum(lambdas)!=0)) {
        Xwaug <- row_augX(Xw, lambdas)
        ywaug <- augY(yw, ncol(Xwaug), lambdas = lambdas, prior.mean = prior.mean)
      } else {
        Xwaug <- Xw
        ywaug <- yw
      }
      if(is.null(start)) start <- rep(1, ncol(Xw))
      if(issparse(Xwaug)){
        coefs <- solve_sparse_CG(A = as(crossprod(Xwaug), "dgCMatrix"), 
                                 b = as.matrix(crossprod(Xwaug, ywaug)),
                                 x0 = as.matrix(start),
                                 tol = control.solver$tol,
                                 maxit = control.solver$maxit,
                                 precond = control.solver$precond)$coefficients
        
      } else {
        coefs <- solve_dense_CG(A = crossprod(Xwaug),
                                b = as.matrix(crossprod(Xwaug, ywaug)),
                                x0 = as.matrix(start),
                                tol = control.solver$tol,
                                maxit = control.solver$maxit,
                                precond = control.solver$precond)$coefficients
      }
      
    } else { # normal equations
      D <- Diagonal(x=1/pen.weights)
      
      if(nrow(Xw) >= ncol(Xw)){
        if(issparse(Xw)){
          coefs <- as.vector(solve_sparse_CG(A = D %*% Matrix::crossprod(Xw) + Diagonal(x=rep(lambda, ncol(Xw))),
                                             b = as.matrix(D %*% Matrix::crossprod(Xw,yw) + lambda*Matrix::crossprod(D, prior.mean)),
                                             x0 = as.matrix(start),
                                             tol = control.solver$tol,
                                             maxit = control.solver$maxit,
                                             precond = control.solver$precond)$coefficients)
        } else { # if it's dense
          coefs <- as.vector(solve_dense_CG(A = as.matrix(D %*% Matrix::crossprod(Xw) + Diagonal(x=rep(lambda, ncol(Xw)))),
                                            b = as.matrix(D %*% Matrix::crossprod(Xw,yw) + lambda*Matrix::crossprod(D, prior.mean)),
                                            x0 = as.matrix(start),
                                            tol = control.solver$tol,
                                            maxit = control.solver$maxit,
                                            precond = control.solver$precond)$coefficients)
        }
      } else { # if p > n
        if(issparse(Xw)){
          # Precompute some values
          DXw <- Matrix::tcrossprod(D, Xw)
          y_star <- as.vector(yw - Xw %*% prior.mean)
          
          # Calculate β using Woodbury identity
          coefs <- as.vector(prior.mean + DXw %*% Matrix::solve(Xw %*% DXw + Diagonal(x=rep(lambda, nrow(Xw))),
                                                                y_star,
                                                                tol = control.solver$tol))
          # TODO replace Matrix::solve with solve_sparse_CG
        } else { # if it's dense
          # Precompute some values
          DXw <- Matrix::tcrossprod(D, Xw)
          y_star <- as.vector(yw - Xw %*% prior.mean)
          
          # Calculate β using Woodbury identity
          coefs <- as.vector(prior.mean + DXw %*% Matrix::solve(Xw %*% DXw + Diagonal(x=rep(lambda, nrow(Xw))),
                                                                y_star,
                                                                tol = control.solver$tol))
          coefs_solve_largep_dense <- as.vector(prior.mean + DXw %*% Matrix::solve(
            as.matrix(Xw %*% DXw + Diagonal(x=rep(lambda, nrow(Xw)))),
            y_star,
            tol=2*.Machine$double.eps))
          # TODO replace Matrix::solve with solve_dense_CG
        }
      }
    }
  } else if (control.solver$method == "chol"){
    if (control.solver$rowaug) {
      # Obtain Xwaug and ywaug
      if(!is.null(lambdas)&(sum(lambdas)!=0)) {
        Xwaug <- row_augX(Xw, lambdas)
        ywaug <- augY(yw, ncol(Xwaug), lambdas = lambdas, prior.mean = prior.mean)
      } else {
        Xwaug <- Xw
        ywaug <- yw
      }
      if(issparse(Xwaug)){
        coefs <- solve_sparse_chol(A = as(crossprod(Xwaug),"dgCMatrix"), 
                                   b = as.matrix(crossprod(Xwaug, ywaug)), 
                                   pivot = control.solver$pivot,
                                   ordering = control.solver$ordering)$coefficients
      } else {
        Xwaug = row_augX(Xw, lambdas)
        coefs <- solve_dense_chol(A = crossprod(Xwaug),
                                  b = as.matrix(crossprod(Xwaug, ywaug)),
                                  pivot = control.solver$pivot)$coefficients
      }
    } else { # normal equations
      D <- Diagonal(x=1/pen.weights)
      
      if(nrow(Xw) >= ncol(Xw)){
        if(issparse(Xw)){
          coefs <- as.vector(solve_sparse_chol(A = D %*% Matrix::crossprod(Xw) + Diagonal(x=rep(lambda, ncol(Xw))),
                                               b = as.matrix(D %*% Matrix::crossprod(Xw,yw) + lambda*Matrix::crossprod(D, prior.mean)),
                                               pivot = control.solver$pivot,
                                               ordering = control.solver$ordering)$coefficients)
        } else { # if it's dense
          coefs <- as.vector(solve_dense_chol(A = as.matrix(D %*% Matrix::crossprod(Xw) + Diagonal(x=rep(lambda, ncol(Xw)))),
                                              b = as.matrix(D %*% Matrix::crossprod(Xw,yw) + lambda*Matrix::crossprod(D, prior.mean)),
                                              pivot = control.solver$pivot)$coefficients)
        }
      } else { # if p > n
        if(issparse(Xw)){
          # Precompute some values
          DXw <- Matrix::tcrossprod(D, Xw)
          y_star <- as.vector(yw - Xw %*% prior.mean)
          # Calculate β using Woodbury identity
          coefs <- as.vector(prior.mean + DXw %*% solve_sparse_chol(A = Xw %*% DXw + Diagonal(x=rep(lambda, nrow(Xw))),
                                                                    b = y_star,
                                                                    pivot = control.solver$pivot,
                                                                    ordering = control.solver$ordering)$coefficients) 
        } else { # if it's dense
          # Precompute some values
          DXw <- Matrix::tcrossprod(D, Xw)
          y_star <- as.vector(yw - Xw %*% prior.mean)
          
          # Calculate β using Woodbury identity
          coefs <- as.vector(prior.mean + DXw %*% solve_dense_chol(A = as.matrix(Xw %*% DXw + Diagonal(x=rep(lambda, nrow(Xw)))),
                                                                   b = y_star,
                                                                   pivot = control.solver$pivot)$coefficients)
        }
      }
    }
  } else {
    stop("Invalid method.")
  }
  coefs <- clip(coefs, lower, upper, reflect_lower = F, reflect_upper = F) # clip coefficients within allowed range - TODO: check if reflect_lower=T might work better to deal with nonnegativity constraints
  ifelse(!is.null(coefs), return(coefs), return(start)) #TODO remove?
}


Xweights <- function(X, weights){
  if(!is.null(weights)){
    return(Diagonal(x=sqrt(weights)) %*% X) # Faster than X*sqrt(weights)
  } else {
    return(X)
  }
}

Yweights <- function(y, weights){
  if (!is.null(weights)) {
    return(y * sqrt(weights))
  } else {
    return(y)
  }
}




# UTILITY FUNCTIONS ####

# Row augmentation of the covariate matrix if lambda!=0 to add ridge penalty
# Ref :
# https://stats.stackexchange.com/a/164546/27340
# https://stats.stackexchange.com/a/203699/27340
row_augX <- function(Xw, lambdas = NULL){
  if (!is.null(lambdas)&(sum(lambdas)!=0)) {
    out <- rbind(Xw, Diagonal(x=sqrt(lambdas))) } else
    { out <- Xw }
  
  if ("dgeMatrix" %in% class(out)) { out = as(out, "dgCMatrix") } # now always returning sparse matrix, as that's most efficient when using row augmentation
  
  return(out)
}

# augment y
# NOTE: for zero-centered ridge regression augmentation is with zeros
# (in a Bayesian framework this would correspond to the MAP with a zero-centered Gaussian prior)
# a nonzero centered ridge regression, centered at the coefficient estimates of the
# previous iteration might also be worth considering (in a Bayesian sense corresponding to the
# MAP with a prior centered at the coefficient estimates of the previous iteration)
# in the latter case, we would have to augment not with zeros but with prior.mean * sqrt(lambdas)
# where prior.mean is a vector with the values at which the ridge penalty should be centered
# see e.g. https://www.sciencedirect.com/science/article/abs/pii/S0047259X21001603
# In an optimisation sense, this would add "momentum" to the optimisation process.
# Not sure how this would perform - might be worth trying quickly.
augY <- function(yw, p, prior.mean = NULL, lambdas = NULL){
  if (!is.null(lambdas)&(sum(lambdas)!=0)) {
    if(is.null(prior.mean) | is.null(lambdas)){
      return(c(yw, rep(0, p))) # zero centered row augmentation
    } else {
      return(c(yw, prior.mean * sqrt(lambdas)))  # nonzero centered row augmentation
    } } else { return(yw) }
}


# diaglambdas <- function(lambdas){ # this one always returns a sparse matrix
#   # if(sparse) {
#   #  dlambdas <- Diagonal(x=lambdas) # faster than .sparseDiagonal(p, lambdas)
#   # } else {
#   #  dlambdas <- Diagonal(x=lambdas) # as.matrix(Diagonal(x=lambdas)) # faster than diag(lambdas, p)
#   # }
#   return (Diagonal(x=lambdas))
# }


issparse <- function(X){
  return(any(class(X) %in% c("dgCMatrix", "dsCMatrix")))
}


# Calculate scaling of columns of matrix X in terms of L2 norm
# INPUT X matrix
scaleX = function (X) return(sqrt(Matrix::colSums(X^2)))

# L2 norm normalize columns of matrix X
# INPUT
#   X: input matrix
# OUTPUT
#   X: matrix with L2 norm of columns normalized to 1
#   attr(X, "scaling"): normalisation factors
#   attr(X, "normalized"): TRUE
normalizeX = function (X) {  
  scaling <-  scaleX(X)
  X <- as.matrix(X %*% Matrix::Diagonal(x = 1 / scaling))
  attr(X, "scaling") <- scaling
  attr(X, "normalized") <- TRUE
  return(X)
}

# center columns of matrix X to zero
centerX <- function(X) {
  Xcenter = Matrix::colMeans(X)
  return(X - rep(Xcenter, rep.int(nrow(X), ncol(X))))
}

# is X normalized?
isnormalized = function (X) { out = !is.null(attr(X, "normalized"))
if (out) out = attr(X, "normalized")
return(out) }

# function to calculate good coefficient initialisations to use as starting values in 
# L0glm.fit which can be used as start - uses an orthogonal approximation
# to a ridge regression
# INPUT:
# X : covariate matrix
# y: observations
# mu: expected value on response scale
# weights: prior observation weights of GLM
# mu: expected value on response scale
# family: GLM family
# lambda: individual ridge penalty on individual coefficients
init_coefs = function( X,
                       y,
                       mu,
                       weights = rep(1, nrow(X)),
                       offset = rep(0, nrow(X)),
                       family = poisson(identity),
                       lambda = 1) { 
  # mu = as.matrix(mu)
  # offset = as.matrix(offset)
  eta = as.matrix(family$linkfun(mu))
  varg = family$variance(mu) # variance
  if (length(lambda) == 1) lambda = rep(lambda, ncol(X))
  
  if (family$link == "identity") {
    # simplified formula and faster than formula in else statement
    W <- W.ret <- as.vector(weights / varg) # 1/variance (* weights) with identity link, where variance = 1 with Gaussian noise, variance = mu with Poisson noise, variance = mu^2 with Gamma noise
    good = which(!is.na(W) & !(W == 0) & !is.infinite(W)) # valid rows # (gprime != 0)
    z.ret <- (y - offset) # working (adjusted) response (returned non-subsetted version)
    z <- z.ret[good] # subset used for fit
  } else {
    gprime = family$mu.eta(eta) # derivative of the expected values
    good = which((weights != 0) &
                   (gprime != 0) &
                   !is.infinite(weights) & !is.na(weights)) # valid rows
    z.ret <- (eta - offset) + (y - mu) / gprime # working (adjusted) response
    z <- z.ret[good]
    W.ret <- weights * as.vector(gprime ^ 2 / varg) # updated weights
    W <- W.ret[good]
  }
  
  XW = Matrix(Diagonal(x = sqrt(W)) %*% X[good, , drop = F], sparse = issparse(X))
  XWaug = Matrix(rbind(XW, Diagonal(x = sqrt(lambda))), sparse = issparse(X))
  zwaug = as.vector(c(z * sqrt(W), rep(0, ncol(X))))
  # orthogonal initialisation (assuming variables are approx orthogonal)
  beta_init = as.vector(crossprod(XWaug, zwaug) / colSums(XWaug ^ 2))
  # colSums(XWaug^2) = diag(crossprod(XWaug)))
  # PS original variables in columns of X can be very far from orthogonal,
  # but the augmentation with a diagonal matrix with sqrt(lambdas) along the diagonal
  # makes the assumption hold better
  
  # possible alternative initialisation
  # diagonal initialisation (assuming system is dominated by diagonal of X)
  # for signal deconvolution problem with n=p this gives the best initial solution
  # only works if X is square though
  # beta_init_diagonal <- z*sqrt(W) / diag(XW)
  
  # scale.factor = 0.1 # (1+lambda[[1]]) # multiply coefs by some scale factor to compensate ridge bias a bit
  # beta_init = beta_init*scale.factor
  return(beta_init)
}


# alternative simpler constant coefficient initialisation
# with sensible, well-performing values for different GLM families
# INPUT
# X design matrix
# family (to allow different starting values for different families, 
# in a Bayesian sense corresponding to different SDs of Gaussian prior, 
# since lambda=1/start^2 and lambda=1/SD prior)
start_coefs = function (X, family=gaussian(identity)) { 
  p = ncol(X)
  scaling = scaleX(X) # L2 norm of columns of X
  if (family$family=="Gamma"|family$family=="inverse.gaussian") { sd = 50 } else { sd = 5}
  startvals = rep(sqrt(sd), p) / scaling
  # sd is corresponding in a Bayesian sense to SD of Gaussian prior with gamma=2 since
  # lambda=1/start^2 and lambda=1/SD prior
  return(startvals)
}


# maximum lambda that would be expected to penalize all penalized coefficients to zero
# this takes into account some terms that are unpenalized (e.g. intercept) & 
# prior coefficient estimates
# see https://stats.stackexchange.com/questions/416144/minimum-and-maximum-regularization-in-l0-pseudonorm-penalized-regression
maxlambda = function (X, weights, z, dispersion=1, no.pen=NULL, coefs, family) {
  z <- as.matrix(z)
  weights <- as.matrix(weights)
  sel_task <- which.max(colMeans(abs(z))) # for multitask select column with max mean abs value?
  z <- z[, sel_task, drop=F]
  weights <- weights[, sel_task, drop=F]
  if (family$family %in% c("binomial", "poisson")) { dispersion = 1 }
  
  # remove contribution of unpenalized terms
  if (!is.null(no.pen)) { X_penalized <- X[, -no.pen, drop=F]
  # Calculate the contribution to z of the unpenalized terms
  unpen_contrib = X[, no.pen, drop=F] %*% coefs[no.pen]
  # Adjust z by subtracting the unpenalized contribution
  z <- z - unpen_contrib
  X <- X_penalized
  }
  
  Xw <- Diagonal(x=sqrt(weights)) %*% X
  zw <- Diagonal(x=sqrt(weights)) %*% z
  
  maxlambda <- (4/dispersion) * max((1/colSums(Xw^2)) * (crossprod(Xw, zw)^2), na.rm=T) 
  # the factor of 4 is to put lambda on the lambda scale instead of the 
  # lambda tilde scale, which is how things are now done in the L0glm function 
  # (cf Frommlet & Nuel paper)
  if (family$family=="gaussian") maxlambda = maxlambda*(1/2) # TODO check this
  
  return(maxlambda)
}

# example
# ml = maxlambda(X=normalizeX(data_gaus$X), weights=L0glm_gaus_presetlambda$weights,
#                z=L0glm_gaus_presetlambda$z, no.pen=1, coefs=L0glm_gaus_presetlambda$coefficients)
# ml
# 2346.189
# abline(v=log10(ml), col="cyan")