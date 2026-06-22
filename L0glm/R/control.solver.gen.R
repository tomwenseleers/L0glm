#' @name control.solver.gen
#' @title Set Options of Constrained Ridge Regularised Weighted Least Squares Regression Solvers
#' 
#' @importFrom RhpcBLASctl get_num_cores get_num_procs
#' @description
#'
#' This functions generates a list with all the required parameters for the lmridge 
#' solver (\code{control.solver.gen}) & allows fine control over its options.
#'
#' @param solver
#' XXX
#' @param method
#' XXX
#' @param maxit
#' The maximum number of iterations.
#' @param standardize
#' XXX
#' @param reltol
#' XXX
#' @param tol
#' XXX
#' @param abstol
#' XXX
#' @param block_size
#' Optional block size used by solvers that support blocked updates.
#' @param precond
#' Optional integer preconditioner code for Eigen iterative and Jacobi SVD
#' methods. Defaults are method-specific.
#' @param pivot
#' Optional integer pivoting code for Eigen QR and Cholesky methods. Defaults
#' are method-specific.
#' @param ordering
#' Optional integer sparse ordering code for Eigen sparse QR and Cholesky
#' methods. Defaults are method-specific.
#' @param threshold
#' Optional threshold for Eigen SVD rank decisions. If \code{NULL}, Eigen's
#' default threshold is used.
#' @param verbose
#' XXX
#' @param rowaug
#' XXX
#'
#' @return
#'
#' The function return a list with the required options to solve a linear
#' weighted least squares \code{lmridge} ridge regression using the given solver
#' and method.
#'
#' @seealso
#'
#' \code{\link{L0glm.fit}}
#'
#' @example examples/control.solver.gen.R
#' @export
control.solver.gen <- function(solver, 
                               method = NULL,
                               maxit = NULL,
                               standardize = NULL,
                               reltol = NULL,
                               tol = NULL,
                               abstol = NULL,
                               block_size = NULL,
                               precond = NULL,
                               pivot = NULL,
                               ordering = NULL,
                               threshold = NULL,
                               verbose = FALSE,
                               rowaug = NULL){
  
  if(!(solver %in% c("bigGlm", "bigglm", "glmnet")) & !is.null(standardize) & verbose) {
    print("Selected solver does not support standardization. To enable standardization, please set the solver to bigGlm.")
  }
  
  if(solver == "osqp"){
    if(is.null(maxit)) {
      maxit <- 1e5
    } else if(maxit <= 0){
      maxit <- 1e6
      if(verbose) print("Control.solver$maxit needs to be a positive integer fot the chosen solver. Value was set to solver default (10000).")
    }
    
    if(is.null(abstol)) abstol <- 1e-7
    if(is.null(reltol)) reltol <- 1e-7
    if(is.null(rowaug)) rowaug <- TRUE
    
    l <- list(verbose = verbose, eps_abs = abstol, eps_rel = reltol,
              max_iter = maxit, 
              polish = FALSE, # polish ADMM solution?
              rowaug = rowaug)
    
  } else if (solver %in% c("bigGlm", "bigglm", "glmnet")) {
    if(is.null(maxit)) {
      maxit <- 1e6
    } else if(maxit <= 0){
      maxit <- 1e7
      if(verbose) print("Control.solver$maxit needs to be a positive integer fot the chosen solver. Value was set to solver default (1e5).")
    }
    
    if(is.null(abstol)) abstol <- 1e-10 # largest value that gave half decent results
    if(is.null(rowaug)) rowaug <- TRUE
    if(rowaug == FALSE && verbose) print("Row augmentation is required with chosen solver.")
    
    if(is.null(standardize)) standardize <- TRUE
    
    
    l <- list(thresh = abstol, maxit = maxit, standardize = standardize)
    
  # } else if (solver == "parglm"){ # removed this solver
  #   if(is.null(tol)) tol <- 1e-8
  #   if(is.null(block_size)) block_size <- 1000
  #   l <- list(epsilon = tol, 
  #             nthreads = RhpcBLASctl::get_num_cores(), # didn't see much improvement if I tried to use more
  #             block_size = block_size,
  #             method = "FAST"  # "FAST"=Armadillo Cholesky LLT or "LINPACK" or "LAPACK"
  #             )
    
  } else if (solver == "nnls"){
    l <- list()
    
  } else if (solver == "eigen"){
    if (is.null(method)){
      if (verbose) print("No method for Eigen solver was provided. Method is set to 7.")
      method <- 7L
    }
    tol_was_null <- is.null(tol)
    if(is.null(tol)) tol <- 1E-7
    if(!is.null(maxit) && maxit < 0) {
      maxit <- NULL
      if(verbose) print("Control.solver$maxit needs to be a non-negative integer for the chosen solver. Value was set to solver default.")
    }
    
    if (method == 1){
      if(is.null(rowaug)) rowaug <- FALSE # TODO: check why QR seems to be only solver that now works well with rowaug <- FALSE
      if(is.null(pivot)) pivot <- 0L
      if(is.null(ordering)) ordering <- 1L
      l <- list(method = "QR", pivot = as.integer(pivot), ordering = as.integer(ordering), rowaug = rowaug)
      
    } else if (method == 2){
      if(is.null(rowaug)) rowaug <- FALSE # TODO: check why QR seems to be only solver that now works well with rowaug <- FALSE
      if(is.null(pivot)) pivot <- 1L
      if(is.null(ordering)) ordering <- 1L
      l <- list(method = "QR", pivot = as.integer(pivot), ordering = as.integer(ordering), rowaug = rowaug)
      
    } else if (method == 3){
      if(is.null(rowaug)) rowaug <- TRUE
      if(is.null(maxit)) maxit <- 5000L
      if(is.null(precond)) precond <- 0L
      l <- list(method = "LSCG", tol = tol, maxit = as.integer(maxit), precond = as.integer(precond), rowaug = rowaug)
      
    } else if (method == 4){
      if(is.null(rowaug)) rowaug <- TRUE
      if(is.null(maxit)) maxit <- 5000L
      if(is.null(precond)) precond <- 1L
      l <- list(method = "LSCG", tol = tol, maxit = as.integer(maxit), precond = as.integer(precond), rowaug = rowaug)
      
    } else if (method == 5){
      if(is.null(rowaug)) rowaug <- TRUE # TODO check why rowaug <- FALSE gives bad results
      if(is.null(maxit)) maxit <- 5000L
      if(is.null(precond)) precond <- 0L
      l <- list(method = "CG", tol = tol, maxit = as.integer(maxit), precond = as.integer(precond), rowaug = rowaug)
      
    } else if (method == 6){
      if(is.null(rowaug)) rowaug <- TRUE # TODO check why rowaug <- FALSE gives bad results
      if(is.null(maxit)) maxit <- 5000L
      if(is.null(precond)) precond <- 1L
      l <- list(method = "CG", tol = tol, maxit = as.integer(maxit), precond = as.integer(precond), rowaug = rowaug)
      
    } else if (method == 7){
      if(is.null(rowaug)) rowaug <- TRUE # TODO check why rowaug <- FALSE gives bad results
      if(is.null(pivot)) pivot <- 0L
      if(is.null(ordering)) ordering <- 0L
      l <- list(method = "chol", pivot = as.integer(pivot), ordering = as.integer(ordering), rowaug = rowaug)
      
    } else if (method == 8){
      if(is.null(rowaug)) rowaug <- TRUE # TODO check why rowaug <- FALSE gives bad results
      if(is.null(pivot)) pivot <- 1L
      if(is.null(ordering)) ordering <- 0L
      l <- list(method = "chol", pivot = as.integer(pivot), ordering = as.integer(ordering), rowaug = rowaug)
    } else if (method == 9){
      if(is.null(rowaug)) rowaug <- TRUE # TODO check why rowaug <- FALSE gives bad results
      if(is.null(maxit)) maxit <- 5000L
      if(is.null(precond)) precond <- 2L
      l <- list(method = "CG", tol = tol, maxit = as.integer(maxit), precond = as.integer(precond), rowaug = rowaug)
    } else if (method == 10){
      if(is.null(rowaug)) rowaug <- TRUE
      if(is.null(precond)) precond <- 0L
      if(is.null(threshold)) threshold <- if (tol_was_null) 0 else tol
      l <- list(method = "SVD", svd = 0L, precond = as.integer(precond), threshold = threshold, rowaug = rowaug)
    } else if (method == 11){
      if(is.null(rowaug)) rowaug <- TRUE
      if(is.null(precond)) precond <- 0L
      if(is.null(threshold)) threshold <- if (tol_was_null) 0 else tol
      l <- list(method = "SVD", svd = 1L, precond = as.integer(precond), threshold = threshold, rowaug = rowaug)
    } else {
      stop("Chosen method for Eigen solvers is unavailable. Please choose a method (1-11).")
    }
  } else {
    stop("Chosen solver is not available.")
  }
  return(l)
}
