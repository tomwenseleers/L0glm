#' Model solvers
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





#' @describeIn L0glm
#'
#' Implementation of different solvers
#' It is called by the upper code \code{L0glm.iwls}.
#'
#' @export
lmridge <- function(X, y, lambdas,
                    weights = NULL,
                    lower = rep(-.Machine$double.xmax,ncol(X)),
                    upper = rep(.Machine$double.xmax,ncol(X)),
                    nonnegative = FALSE,
                    solver = "osqp",
                    method = NULL,
                    start = rep(0, ncol(X)),
                    control.solver = NULL,
                    verbose){


  # If nonnegative is true and no lower bound is provided, or the lower bound is larger than 0, it is clipped to 0
  if (nonnegative && all(lower < 0)) lower <- rep(0, ncol(X))

  # Initialize control.solver variable
  if(is.null(control.solver)) call.solver <- list(solver = solver, method = method)


  if (solver == "osqp") {
    if(is.null(control.solver)) control.solver <- do.call("control.solver.gen", call.solver)
    coefs <- lmridge_osqp(X = X, y = y, weights = weights, lambdas = lambdas, start = start,
                          lower = lower, upper = upper, control.solver = control.solver)

  } else if (solver == "glmnet") {
    if(is.null(control.solver)) control.solver <- do.call("control.solver.gen", call.solver)
    coefs <- lmridge_glmnet(X = X, y = y, weights = weights, lambdas = lambdas,
                            lower = lower, upper = upper, control.solver = control.solver)


  } else if (solver == "eigen") {
    if(is.null(control.solver)) control.solver <- do.call("control.solver.gen", call.solver)
    coefs <- lmridge_eigen(X = X, y = y, weights = weights, pen_weights = pen_weights, lambdas = lambdas, start = start,
                           lower = lower, upper = upper, control.solver = control.solver)

  } else if (solver == "nnls") {
    if(is.null(control.solver)) control.solver <- do.call("control.solver.gen", call.solver)
    coefs <- lmridge_nnls(X = X, y = y, weights = weights, lambdas = lambdas,
                          lower = lower, upper = upper, control.solver = control.solver)
  } else {
    stop("Invalid solver choice") }

  return(coefs)
}


lmridge_osqp <- function(X, y, weights, lambdas, lower, upper, start, control.solver){
  Xw <- Xweights(X, weights)
  yw <- Yweights(y, weights)

  if(control.solver$rowaug == TRUE){
    Xwaug <- row_augX(Xw, lambdas)
    XtX <- crossprod(Xwaug)
    Xty <- crossprod(Xwaug, augY(yw, lambdas, ncol(Xw)))
  } else {
    D <- diag(1/pen_weights)
    if (issparse(Xw)){
      Dsp <- Matrix(D, sparse=T)
      XtX <- Dsp %*% crossprod(Xw) + diaglambdas(lambdas, issparse(Xw), ncol(Xw))
      Xty <- Dsp %*% crossprod(Xw, yw)
    } else {
      XtX <- D %*% crossprod(Xw) + diaglambdas(lambdas, issparse(Xw), ncol(Xw))
      Xty <- D %*% crossprod(Xw, yw)
    }
  }

  #Obtain settings  and create model
  settings <- do.call(osqpSettings, control.solver[-7])
  if (sum(start)!=0) settings <- c(settings, warm_start = TRUE) else settings <- c(settings, warm_start = FALSE)

  model <- osqp(P = XtX, q = -Xty, A = .sparseDiagonal(ncol(X)), l = lower, u = upper, pars = settings)

  if (!is.null(start)) model$WarmStart(x=start)

  #Solve model
  coefs <- model$Solve()$x
  coefs <- pmax(lower, pmin(coefs, upper)) # fitted coefficients sometimes go very slightly outside constraint zone due to numerical inaccuracies in solver - this is fixed here via clipping
  return(coefs)
}

lmridge_glmnet <- function(X, y, weights, lambdas, lower, upper, control.solver){
  # Set path depending on dimensions of X
  if(sum(lambdas) == 0 && ncol(X) >= nrow(X)){ # avoid numerical instabilities if lambdas is 0 and p>=n
    path = TRUE
  } else {
    path = FALSE
  }

  # Obtain Xwaug and ywaug
  Xw <- Xweights(X, weights)
  yw <- Yweights(y, weights)

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

lmridge_nnls <- function(X, y, weights, lambdas, lower, upper, control.solver){
  # Obtain Xwaug and ywaug
  Xw <- Xweights(X, weights)
  yw <- Yweights(y, weights)

  Xwaug <- row_augX(Xw, lambdas)
  ywaug <- augY(yw, lambdas, ncol(Xw))

  coefs <- nnls(A=Xwaug, b=ywaug)$x
  return(coefs)
}

lmridge_eigen <- function(X, y, weights, pen_weights, lambdas, lower, upper, start, control.solver = control.solver){
  Xw <- Xweights(X, weights)
  yw <- Yweights(y, weights)

  if (is.null(start)) { start = rep(0,ncol(X))  }

  if(control.solver$method == "LSCG") {
    D <- diag(1/pen_weights)
    if(issparse(Xw)){
      Dsp <- Matrix(D, sparse=T)
      coefs <- solve_sparse_LSCG(A = Dsp %*% crossprod(Xw) + diaglambdas(lambdas, issparse(Xw), ncol(Xw)),
                                 b = as.vector(Dsp %*% crossprod(Xw,yw)),
                                 x0 = start,
                                 tol = control.solver$tol,
                                 maxit = control.solver$maxit,
                                 precond = control.solver$precond)$coefficients
    } else {
      coefs <- solve_sparse_LSCG(A = D %*% crossprod(Xw) + diaglambdas(lambdas, issparse(Xw), ncol(Xw)),
                                 b = as.vector(D %*% crossprod(Xw,yw)),
                                 x0 = start,
                                 tol = control.solver$tol,
                                 maxit = control.solver$maxit,
                                 precond = control.solver$precond)$coefficients
    }

  } else if (control.solver$method == "QR") {
    if(control.solver$rowaug == TRUE){
      if(issparse(X)){
        coefs <- solve_sparse_QR(A = row_augX(Xw, lambdas),
                                 b = augY(yw, lambdas, ncol(Xw)))$coefficients
      } else {
        coefs <- solve_dense_QR(A = row_augX(Xw, lambdas),
                                b = augY(yw, lambdas, ncol(Xw)), pivot = control.solver$pivot)$coefficients
      }

    } else {
      D <- diag(1/pen_weights)
      if (issparse(Xw)){
        Dsp <- Matrix(D, sparse=T)
        coefs <- solve_sparse_QR(A = Dsp %*% crossprod(Xw) + diaglambdas(lambdas, issparse(Xw), ncol(Xw)),
                                 b = as.vector(Dsp %*% crossprod(Xw,yw)))$coefficients
      } else {
        coefs <- solve_dense_QR(A = D %*% crossprod(Xw) + diaglambdas(lambdas, issparse(Xw), ncol(Xw)),
                                b = as.vector(D %*% crossprod(Xw,yw)),
                                pivot = control.solver$pivot)$coefficients
      }
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
    D <- diag(1/pen_weights)
    if(issparse(Xw)){
      Dsp <- Matrix(D, sparse=T)
      coefs <- solve_sparse_chol(A = Dsp %*% crossprod(Xw) + diaglambdas(lambdas, issparse(Xw), ncol(Xw)),
                                 b = as.vector(Dsp %*% crossprod(Xw,yw)), pivot = control.solver$pivot,
                                 ordering = control.solver$ordering)$coefficients
    } else {
      coefs <- solve_dense_chol(A = D %*% crossprod(Xw) + diaglambdas(lambdas, issparse(Xw), ncol(Xw)),
                                b = as.vector(D %*% crossprod(Xw,yw)),
                                pivot = control.solver$pivot,
                                orderind = control.solver$ordering)$coefficients
    }
  } else {
    stop("Invalid method.")
  }
  coefs <- pmax(lower, pmin(coefs, upper)) # clip coefficients within allowed range
  return(coefs)
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
  if (sum(lambdas) != 0) {
    return(rbind(Xw, diaglambdas(lambdas, issparse(Xw), ncol(Xw))))
  } else {
    return(Xw)
  }
}

augY <- function(yw, lambdas, p){
  return(if (sum(lambdas)>0) c(yw, rep(0, p)) else yw)
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
