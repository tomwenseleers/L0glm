#' Set parameters for L0glm2 algorithm
#'
#' @description
#'
#' The functions generate a list with all required parameters for the control
#' the adaptive ridge iteration (\code{control.l0.gen}), the iteratively
#' weighted least squares for GLM fitting (\code{control.iwls.gen2}), or the
#' least squares fitting (\code{control.fit.gen2}).
#'
#' @param maxit
#' The maximum number of iterations.
#' @param reltol
#' Coefficients are update until the relative difference between two iteration
#' is smaller than this threshold. The relative difference is defined as
#' \eqn{2 \frac{|coef_i^{(k-1)} - coef_i^{(k)}|}{coef_i^{(k-1)} + coef_i^{(k)}}}
#' @param delta
#' A small constant added to the estimated coefficients when computing the
#' adaptive ridge penalty weights.
#' @param gamma
#' The exponent used for computing the adaptive ridge penalty weights. It is
#' recommended to work with \code{gamma = 2} (Frommlet et Nuel 2016).
#' @param warn
#' Should warnings from the iteration be produced ?
#' @param thresh
#' A threshold under which the fitted coefficients are clipped to 0.
#' @param block.size
#' Number of covariates to fit simultaneously. If \code{block.size} is smaller
#' than \code{p}, the fit is performed in different blocks which is updated
#' cyclically.
#' @param tol
#' Stop iteration when the difference in coefficients between iterations is
#' smaller than this threshold. The coefficient difference is defined as:
#' \eqn{\frac{\|coefs^{(k-1)} - coefs^{(k)}}{\|coefs^{(k-1)}\|}}
#'
#' @return
#'
#' The functions return a list with the required parameters for computing an
#' \code{L0glm2} fit.
#'
#' @references
#'
#' Frommlet, F. and G., Nuel. (2016) An adaptive ridge procedure for L0 regularization. PloS one.
#'
#' @seealso
#'
#' \code{\link{L0glm2}}
#'
#' @examples
#'
#' control.l0.gen()
#' control.iwls.gen()
#' control.fit.gen()
#'
#' @export
control.l0.gen <- function(maxit = 100, reltol = 1E-5, delta = 1E-5, gamma = 1.5, warn = F){
  l <- list(maxit = maxit, reltol = reltol, delta = delta, gamma = gamma, warn = warn)
  return(l)
}


#' @describeIn control.l0.gen
#'
#' Generate parameters that control the IWLS iteration.
#'
#' @export
control.iwls.gen <- function(maxit = 1, reltol = 1E-5, thresh = 1E-7, warn = F){ # TODO review, reltol used to be 1E-4
  l <- list(maxit = maxit, reltol = reltol, thresh = thresh, warn = warn)
  return(l)
}


#' @describeIn control.l0.gen
#'
#' Generates parameters that control the least square fit.
#'
#' @export
control.fit.gen <- function(maxit = 10, block.size = NULL, tol = 1E-10){
  l <- list(block.size = block.size, maxit = maxit, tol = tol)
  return(l)
}


#' @describeIn control.l0.gen
#'
#' Generates parameters that control the solver (lmridge function).
#'
#' @export
control.solver.gen <- function(solver, method = NULL,
                               maxit = NULL,
                               blocksize = NULL,
                               reltol = 1E-8,
                               tol = NULL,
                               abstol = NULL,
                               verbose = FALSE,
                               rowaug = NULL){
  if(!is.null(method) & solver != "eigen" & verbose == TRUE){
    print("A method is only required when using an Eigen solver.")
  }
  
  if(solver == "osqp"){
    if(is.null(maxit)) maxit <- 10000
    if(is.null(abstol)) abstol <- 1e-8
    if(is.null(rowaug)) rowaug <- TRUE
    
    l <- list(verbose = verbose, eps_abs = abstol, eps_rel = reltol,
              max_iter = maxit, linsys_solver = 0L, polish = FALSE, 
              rowaug = rowaug)
    
  } else if (solver == "glmnet") {
    if(is.null(maxit)) maxit <- 1e5
    if(is.null(abstol)) abstol <- 1e-10
    if(is.null(rowaug)) rowaug <- TRUE
    if(rowaug == FALSE) print("Only row augmentation is allowed with this method.")
    
    l <- list(thresh = abstol, maxit = maxit)
    
  } else if (solver == "nnls"){
    l <- list()
    
  } else if (solver == "eigen"){
    if (is.null(method)){
      if (verbose) print("No method for Eigen solver was provided. Method is set to 1.")
    }
    if(is.null(tol)) tol <- 1E-10
    
    if (method == 1){
      if(is.null(rowaug)) rowaug <- TRUE
      if(rowaug == FALSE & verbose == TRUE) print("Row augmentation is required with your chosen method. Using row augmentation.") 
      # TODO : review this? If cant use normal equations, remove rowaug as a return
      l <- list(method = "QR", pivot = 0L, rowaug = TRUE)
      
    } else if (method == 2){
      if(is.null(rowaug)) rowaug <- TRUE
      l <- list(method = "QR", pivot = 1, rowaug = rowaug)
      
    } else if (method == 3){
      if(is.null(rowaug)) rowaug <- FALSE
      if(rowaug == TRUE & verbose == TRUE) print("Row augmentation is not allowed with your chosen method. Using normal equations.")
      
      l <- list(method = "LSCG", tol = tol, maxit = 0L, precond = 0L)
      
    } else if (method == 4){
      if(is.null(rowaug)) rowaug <- FALSE
      
      if(rowaug == TRUE & verbose == TRUE) print("Row augmentation is not allowed with your chosen method. Using normal equations.")
      
      l <- list(method = "LSCG", tol = tol, maxit = 0L, precond = 1L)
      
    } else if (method == 5){
      l <- list(method = "CG", tol = tol, maxit = 0L, precond = 0L)
      
    } else if (method == 6){
      l <- list(method = "CG", tol = tol, maxit = 0L, precond = 1L)
      
    } else if (method == 7){
      if(is.null(rowaug)) rowaug <- FALSE
      if(rowaug == TRUE & verbose == TRUE) print("Row augmentation is not allowed with your chosen method. Using normal equations.")
      l <- list(method = "chol", pivot = 0L, ordering = 0L)
      
    } else if (method == 8){
      if(is.null(rowaug)) rowaug <- FALSE
      if(rowaug == TRUE & verbose == TRUE) print("Row augmentation is not allowed with your chosen method. Using normal equations.")
      l <- list(method = "chol", pivot = 1L, ordering = 0L)
    } else {
      stop("Chosen method for Eigen solvers is unavailable. Please choose a method (1-8).")
    }
  } else {
    stop("Chosen solver is not available.")
  }
  return(l)
}




#' @describeIn control.l0.gen
#'
#' Generates parameters that control the block fitting algorithm.
#'
#' @export
control.block.gen <- function(maxit = 70, block.size = NULL, tol = 1E-10){
  l <- list(block.size = block.size, maxit = maxit, tol = tol)
  return(l)
}


