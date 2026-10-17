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
#' Solver-specific method. For \code{solver="eigen"}, the exposed methods are
#' \code{1L}=Cholesky LLT, \code{2L}=Cholesky LDLT, \code{3L}=CG with
#' incomplete Cholesky preconditioning for sparse systems or diagonal
#' preconditioning for dense systems, \code{4L}=LSCG with least-squares
#' diagonal preconditioner, \code{5L}=column-pivoted QR, and
#' \code{6L}=BDCSVD for dense systems or Jacobi SVD for sparse systems.
#' @param maxit
#' The maximum number of iterations.
#' @param standardize
#' XXX
#' @param reltol
#' XXX
#' @param tol
#' Relative convergence tolerance for Eigen iterative solvers. For
#' LSCG, Eigen tests the normal-equation residual in the original coefficient
#' coordinates; no additional external column equilibration is performed.
#' @param abstol
#' XXX
#' @param block_size
#' Optional block size used by solvers that support blocked updates.
#' @param precond
#' Optional integer preconditioner code for Eigen Jacobi SVD. For exposed
#' Eigen iterative methods this is fixed by the chosen method: method
#' \code{3L} uses incomplete Cholesky CG for sparse systems and diagonal CG
#' for dense systems; method \code{4L} uses LSCG with
#' least-squares diagonal preconditioning.
#' @param pivot
#' Optional integer pivoting code for Eigen QR and Cholesky methods. Defaults
#' are method-specific.
#' @param ordering
#' Optional integer sparse ordering code for Eigen sparse QR and Cholesky
#' methods. For sparse QR, \code{0L}=AMD, \code{1L}=COLAMD and
#' \code{2L}=NaturalOrdering, with \code{1L} the default. For sparse Cholesky,
#' \code{0L}=AMD and \code{1L}=NaturalOrdering, with \code{0L} the default.
#' @param NaturalOrdering
#' Optional logical flag for Eigen sparse Cholesky ordering. If \code{TRUE},
#' sparse Cholesky LLT/LDLT use Eigen's \code{NaturalOrdering}; if
#' \code{FALSE}, they use \code{AMDOrdering}. The default is \code{FALSE}.
#' The legacy integer \code{ordering} argument is still accepted and takes
#' precedence when supplied.
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
                               NaturalOrdering = NULL,
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
    
  } else if (solver == "eigen"){
    if (is.null(method)){
      if (verbose) print("No method for Eigen solver was provided. Method is set to 1.")
      method <- 1L
    }
    if(is.null(NaturalOrdering)) NaturalOrdering <- FALSE
    NaturalOrdering <- isTRUE(NaturalOrdering)
    tol_was_null <- is.null(tol)
    if(is.null(tol)) tol <- 1E-7
    if(!is.null(maxit) && maxit < 0) {
      maxit <- NULL
      if(verbose) print("Control.solver$maxit needs to be a non-negative integer for the chosen solver. Value was set to solver default.")
    }
    
    if (method == 1){
      if(is.null(rowaug)) rowaug <- TRUE # TODO check why rowaug <- FALSE gives bad results
      if(is.null(pivot)) pivot <- 0L
      if(is.null(ordering)) ordering <- if (NaturalOrdering) 1L else 0L
      NaturalOrdering <- as.integer(ordering[1L]) == 1L
      l <- list(method = "chol", pivot = as.integer(pivot),
                ordering = as.integer(ordering),
                NaturalOrdering = NaturalOrdering,
                rowaug = rowaug)
      
    } else if (method == 2){
      if(is.null(rowaug)) rowaug <- TRUE # TODO check why rowaug <- FALSE gives bad results
      if(is.null(pivot)) pivot <- 1L
      if(is.null(ordering)) ordering <- if (NaturalOrdering) 1L else 0L
      NaturalOrdering <- as.integer(ordering[1L]) == 1L
      l <- list(method = "chol", pivot = as.integer(pivot),
                ordering = as.integer(ordering),
                NaturalOrdering = NaturalOrdering,
                rowaug = rowaug)
    } else if (method == 3){
      if(is.null(rowaug)) rowaug <- TRUE # TODO check why rowaug <- FALSE gives bad results
      if(is.null(maxit)) maxit <- 5000L
      if(!is.null(precond) && as.integer(precond[1L]) != 2L && verbose) {
        print("Eigen CG method 3 uses incomplete Cholesky preconditioning; precond was set to 2.")
      }
      precond <- 2L
      l <- list(method = "CG", tol = tol, maxit = as.integer(maxit), precond = as.integer(precond), rowaug = rowaug)
    } else if (method == 4){
      if(is.null(rowaug)) rowaug <- TRUE
      if(is.null(maxit)) maxit <- 5000L
      if(!is.null(precond) && as.integer(precond[1L]) != 1L && verbose) {
        print("Eigen LSCG method 4 uses least-squares diagonal preconditioning; precond was set to 1.")
      }
      precond <- 1L
      l <- list(method = "LSCG", tol = tol, maxit = as.integer(maxit), precond = as.integer(precond), rowaug = rowaug)
    } else if (method == 5){
      if(is.null(rowaug)) rowaug <- FALSE
      if(!is.null(pivot) && as.integer(pivot[1L]) != 1L && verbose) {
        print("Eigen QR method 5 uses column-pivoted QR; pivot was set to 1.")
      }
      pivot <- 1L
      if(is.null(ordering)) ordering <- 1L
      l <- list(method = "QR", pivot = as.integer(pivot), ordering = as.integer(ordering), rowaug = rowaug)
    } else if (method == 6){
      if(is.null(rowaug)) rowaug <- TRUE
      if(is.null(precond)) precond <- 0L
      if(is.null(threshold)) threshold <- if (tol_was_null) 0 else tol
      l <- list(method = "SVD", svd = 0L, sparse.svd = 1L, precond = as.integer(precond), threshold = threshold, rowaug = rowaug)
    } else {
      stop("Chosen method for Eigen solvers is unavailable. Please choose a method (1-6).")
    }
  } else {
    stop("Chosen solver is not available.")
  }
  return(l)
}
