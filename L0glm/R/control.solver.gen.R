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
                               verbose = FALSE,
                               rowaug = NULL){
  
  if(solver != "glmnet" & !is.null(standardize) & verbose) {
    print("Selected solver does not support standardization. To enable standardization, please set the solver to glmnet.")
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
    
  } else if (solver == "glmnet") {
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
    if(is.null(tol)) tol <- 1E-10
    
    if (method == 1){
      if(is.null(rowaug)) rowaug <- FALSE # TODO: check why QR seems to be only solver that now works well with rowaug <- FALSE
      l <- list(method = "QR", pivot = 0L, rowaug = rowaug)
      
    } else if (method == 2){
      if(is.null(rowaug)) rowaug <- FALSE # TODO: check why QR seems to be only solver that now works well with rowaug <- FALSE
      l <- list(method = "QR", pivot = 1L, rowaug = rowaug)
      
    } else if (method == 3){
      if(is.null(rowaug)) rowaug <- TRUE
      l <- list(method = "LSCG", tol = 1E-7, maxit = 0L, precond = 0L, rowaug = rowaug)
      
    } else if (method == 4){
      if(is.null(rowaug)) rowaug <- TRUE
      l <- list(method = "LSCG", tol = 1E-7, maxit = 0L, precond = 1L, rowaug = rowaug)
      
    } else if (method == 5){
      if(is.null(rowaug)) rowaug <- TRUE # TODO check why rowaug <- FALSE gives bad results
      l <- list(method = "CG", tol = 1E-7, maxit = 0L, precond = 0L, rowaug = rowaug)
      
    } else if (method == 6){
      if(is.null(rowaug)) rowaug <- TRUE # TODO check why rowaug <- FALSE gives bad results
      l <- list(method = "CG", tol = 1E-7, maxit = 0L, precond = 1L, rowaug = rowaug)
      
    } else if (method == 7){
      if(is.null(rowaug)) rowaug <- TRUE # TODO check why rowaug <- FALSE gives bad results
      l <- list(method = "chol", pivot = 0L, ordering = 0L, rowaug = rowaug)
      
    } else if (method == 8){
      if(is.null(rowaug)) rowaug <- TRUE # TODO check why rowaug <- FALSE gives bad results
      l <- list(method = "chol", pivot = 1L, ordering = 0L, rowaug = rowaug)
    } else {
      stop("Chosen method for Eigen solvers is unavailable. Please choose a method (1-8).")
    }
  } else {
    stop("Chosen solver is not available.")
  }
  return(l)
}
