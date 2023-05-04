#' Implementation of solvers
#'
#' @import nnls
#' @import Matrix
#' @import glmnet
#' @importFrom osqp osqpSettings osqp solve_osqp
#' @importFrom stats lm.fit
#' @importFrom utils flush.console
#' @importFrom qlcMatrix colMax
#'
#' @description
#' Fits L0 penalized model using specified solver. This code contains minimum 





#' @describeIn Lmridge
#'
#' Implementation of different solvers
#' It is called by the upper code \code{L0glm.iwls}.
#'
#' @export
lmridge <- function(X, y, lambdas,
                    pen.weights,
                    weights = NULL,
                    lower = rep(-.Machine$double.xmax,ncol(X)),
                    upper = rep(.Machine$double.xmax,ncol(X)),
                    nonnegative = FALSE,
                    block.size = NULL,
                    solver = "osqp",
                    method = NULL,
                    start = rep(0, ncol(X)),
                    control.solver = NULL,
                    verbose){

  # If nonnegative is true and no lower bound is provided, or the lower bound is larger than 0, it is clipped to 0
  if (nonnegative && all(lower < 0)) lower <- rep(0, ncol(X))
  # Initialize control.solver variable
  if(is.null(control.solver)) call.solver <- list(solver = solver, method = method)
  
  if(is.null(block.size) || block.size == ncol(X)) {
    Xw <- Xweights(X, weights)
    yw <- Yweights(y, weights)
  } else { # If block.size < ncol(X), Xw and yw are calculated in block.lmridge function
    Xw <- X
    yw <- y
  }
  
  
  if (solver == "osqp") {
    if(is.null(control.solver)) control.solver <- do.call("control.solver.gen", call.solver)
    coefs <- lmridge.osqp(Xw = Xw, yw = yw, weights = weights, lambdas = lambdas, start = start,
                          lower = lower, upper = upper, control.solver = control.solver)
    
  } else if (solver == "glmnet") {
    if(is.null(control.solver)) control.solver <- do.call("control.solver.gen", call.solver)
    coefs <- lmridge.glmnet(Xw = Xw, yw = yw, weights = weights, lambdas = lambdas,
                            lower = lower, upper = upper, control.solver = control.solver)
    
    
  } else if (solver == "eigen") {
    if(is.null(control.solver)) control.solver <- do.call("control.solver.gen", call.solver)
    coefs <- lmridge.eigen(Xw = Xw, yw = yw, weights = weights, pen.weights = pen.weights,
                           lambdas = lambdas, start = start,
                           lower = lower, upper = upper, control.solver = control.solver)
    
  } else if (solver == "nnls") {
    if(is.null(control.solver)) control.solver <- do.call("control.solver.gen", call.solver)
    coefs <- lmridge.nnls(Xw = Xw, yw = yw, weights = weights, lambdas = lambdas,
                          lower = lower, upper = upper, control.solver = control.solver)
  } else {
    stop("Invalid solver choice") }
  return(coefs)
}


lmridge.osqp <- function(Xw, yw, weights, lambdas, lower, upper, start, control.solver){
  Xwaug <- row_augX(Xw, lambdas)
  if(control.solver$rowaug == TRUE){
    XtX <- crossprod(Xwaug)
    Xty <- crossprod(Xwaug, augY(yw, lambdas, ncol(Xwaug)))
  } else {
    if (issparse(Xw)){
      XtX <- as(crossprod(Xwaug),"dgCMatrix")
    } else {
      XtX <- crossprod(Xwaug)
    }
    Xty <- as.matrix(crossprod(Xwaug, augY(yw, lambdas, ncol(Xw))))
  }
  
  #Obtain settings  and create model
  settings <- do.call(osqpSettings, control.solver[-7])
  if (sum(start)!=0) settings <- c(settings, warm_start = TRUE) else settings <- c(settings, warm_start = FALSE)
  
  model <- osqp(P = XtX, q = -Xty, A = .sparseDiagonal(ncol(Xw)), l = lower, u = upper, pars = settings)
  if (!is.null(start)) model$WarmStart(x=start)
  
  #Solve model
  coefs <- model$Solve()$x
  coefs <- pmax(lower, pmin(coefs, upper)) # fitted coefficients sometimes go very slightly outside constraint zone due to numerical inaccuracies in solver - this is fixed here via clipping
  return(coefs)
}

lmridge.glmnet <- function(Xw, yw, weights, lambdas, lower, upper, control.solver){
  # Set path depending on dimensions of Xw
  if(sum(lambdas) == 0 && ncol(Xw) >= nrow(Xw)){ # avoid numerical instabilities if lambdas is 0 and p>=n
    path = TRUE
  } else {
    path = FALSE
  }
  
  # Obtain Xwaug and ywaug
  Xwaug <- row_augX(Xw, lambdas)
  ywaug <- augY(yw, lambdas, ncol(Xw))
  model <- bigGlm(x = Xwaug, y = ywaug,
                  intercept = FALSE, standardize = FALSE, standardize.response = FALSE,
                  lower.limits = c(-Inf, lower), upper.limits = c(Inf, upper),
                  family = "gaussian", path = FALSE, thresh = control.solver$thresh,
                  maxit = control.solver$maxit)
  coefs <- coef(model)[-1]
  coefs <- pmax(lower, pmin(coefs, upper)) #TODO ask TOM
  return(coefs)
}

lmridge.nnls <- function(Xw, yw, weights, lambdas, lower, upper, control.solver){
  # Obtain Xwaug and ywaug
  Xwaug <- row_augX(Xw, lambdas)
  ywaug <- augY(yw, lambdas, ncol(Xw))
  
  coefs <- nnls(A=Xwaug, b=ywaug)$x
  return(coefs)
}

lmridge.eigen <- function(Xw, yw, weights, pen.weights, lambdas, lower, upper, start, control.solver = control.solver){
  if (is.null(start)) { start = rep(0,ncol(X))  }
  
  if(control.solver$method == "QR") {
    if(control.solver$rowaug == TRUE){
      if(issparse(Xw)){
        coefs <- solve_sparse_QR(A = row_augX(Xw, lambdas),
                                 b = augY(yw, lambdas, ncol(Xw)))$coefficients
      } else {
        coefs <- solve_dense_QR(A = row_augX(Xw, lambdas),
                                b = augY(yw, lambdas, ncol(Xw)), pivot = control.solver$pivot)$coefficients
      }
    } else {
      # Xwaug = row_augX(Xw, lambdas)
      # if(issparse(Xw)){
      #   coefs <- solve_sparse_QR(A = as(crossprod(Xwaug),"dgCMatrix"),
      #                            b = as.matrix(crossprod(Xwaug, augY(yw, lambdas, ncol(Xw)))))$coefficients
      # 
      #   
      #   } else {
      #   coefs <- solve_dense_QR(A = crossprod(Xwaug),
      #                           b = as.matrix(crossprod(Xwaug, augY(yw, lambdas, ncol(Xw)))),
      #                           pivot = control.solver$pivot)$coefficients}
      D <- diag(1/pen.weights) # cf notation Liu & Li 2016
      if(issparse(Xw)){
        Dsp <- Matrix(D, sparse=T)
        coefs <- solve_sparse_QR(A = Dsp %*% crossprod(Xw) + diaglambdas(lambdas, issparse(Xw), ncol(Xw)),
                                 b = as.vector(Dsp %*% crossprod(Xw,yw)))$coefficients # TODO add ordering options to control gen for QR?
      } else {
        coefs <- solve_dense_QR(A = D %*% crossprod(Xw) + diaglambdas(lambdas, issparse(Xw), ncol(Xw)),
                                b = as.vector(D %*% crossprod(Xw,yw)),
                                pivot = control.solver$pivot)$coefficients
      }
    }
    
  } else if (control.solver$method == "LSCG") {
    #browser()
    Xwaug = row_augX(Xw, lambdas)
    if(issparse(Xw)){
      coefs <- solve_sparse_LSCG(A = as(crossprod(Xwaug),"dgCMatrix"),
                                 b = as.matrix(crossprod(Xwaug, augY(yw, lambdas, ncol(Xw)))),
                                 x0 = start,
                                 tol = control.solver$tol,
                                 maxit = control.solver$maxit,
                                 precond = control.solver$precond)$coefficients
    } else {
      coefs <- solve_sparse_LSCG(A = crossprod(Xwaug),
                                 b = as.matrix(crossprod(Xwaug, augY(yw, lambdas, ncol(Xw)))),
                                 x0 = start,
                                 tol = control.solver$tol,
                                 maxit = control.solver$maxit,
                                 precond = control.solver$precond)$coefficients
    }
  } else if (control.solver$method == "CG") {
    Xwaug <- row_augX(Xw, lambdas)
    if(issparse(Xw)){
      coefs <- solve_sparse_CG(A = as(crossprod(Xwaug), "dgCMatrix"),
                               b = as.vector(crossprod(Xwaug, augY(yw, lambdas, ncol(Xw)))),
                               x0 = start,
                               tol = control.solver$tol,
                               maxit = control.solver$maxit,
                               precond = control.solver$precond)$coefficients
      
    } else {
      coefs <- solve_dense_CG(A = as(crossprod(Xwaug), "dgCMatrix"),
                              b = as.vector(crossprod(Xwaug, augY(yw, lambdas, ncol(Xw)))),
                              x0 = start,
                              tol = control.solver$tol,
                              maxit = control.solver$maxit,
                              precond = control.solver$precond)$coefficients
    }
  
    
  } else if (control.solver$method == "chol"){
    if(issparse(Xw)){
      Xwaug = row_augX(Xw, lambdas)
      coefs <- solve_sparse_chol(A = as(crossprod(Xwaug),"dgCMatrix"),
                                 b = as.matrix(crossprod(Xwaug, augY(yw, lambdas, ncol(Xw)))), 
                                 pivot = control.solver$pivot,
                                 ordering = control.solver$ordering)$coefficients
    } else {
      Xwaug = row_augX(Xw, lambdas)
      coefs <- solve_dense_chol(A = crossprod(Xwaug),
                                b = as.matrix(crossprod(Xwaug, augY(yw, lambdas, ncol(Xw)))),
                                pivot = control.solver$pivot)$coefficients
    }
  } else {
    stop("Invalid method.")
  }
  coefs <- pmax(lower, pmin(coefs, upper)) # clip coefficients within allowed range
  ifelse(!is.null(coefs),return(coefs), return(start))
}


Xweights <- function(X, weights){
  if(!is.null(weights)){
    return(X * sqrt(weights))
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

# Row augmentation of the covariate matrix if lambda!=0 to add ridge penalty
# Ref :
# https://stats.stackexchange.com/a/164546/27340
# https://stats.stackexchange.com/a/203699/27340
row_augX <- function(Xw, lambdas){
  lambdas <- sqrt(lambdas)
  return(rbind(Xw, diaglambdas(lambdas, issparse(Xw), ncol(Xw))))
}

augY <- function(yw, lambdas, p){
  #return(if (sum(lambdas)>0) c(yw, rep(0, p)) else yw) TODO remove?
  return(c(yw, rep(0, p)))
}

diaglambdas <- function(lambdas, sparse, p){
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
