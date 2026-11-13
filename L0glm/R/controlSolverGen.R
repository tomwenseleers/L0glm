#' @name controlSolverGen
#' @title Set Options of Constrained Ridge Regularised Weighted Least Squares Regression Solvers
#' 
#' @importFrom RhpcBLASctl get_num_cores get_num_procs
#' @description
#'
#' This functions generates a list with all the required parameters for the lmRidge
#' solver (\code{controlSolverGen}) & allows fine control over its options.
#'
#' @param solver
#' XXX
#' @param method
#' Solver-specific method. For \code{solver="eigen"}, the exposed methods are
#' \code{1L}=Cholesky LLT, \code{2L}=Cholesky LDLT, \code{3L}=CG with
#' incomplete Cholesky preconditioning by default, \code{4L}=LSCG with least-squares
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
#' CG the default is \code{1e-10}; LSCG defaults to \code{1e-12}. LSCG
#' checks the normal-equation residual in the original coefficient coordinates
#' by recomputing it from the returned coefficients. A solve that does not
#' reach the requested accuracy uses the Cholesky LLT fallback.
#' @param abstol
#' XXX
#' @param block_size
#' Optional block size used by solvers that support blocked updates.
#' @param precond
#' For Eigen CG (\code{method=3L}), \code{0L} selects no preconditioning,
#' \code{1L} diagonal preconditioning and \code{2L} incomplete Cholesky
#' (the default for both dense and sparse systems). Dense CG retains its dense
#' matrix operator; only the incomplete Cholesky preconditioner uses a sparse
#' copy containing every nonzero entry. For fully populated systems this
#' factorization can cost more time and memory than diagonal preconditioning.
#' LSCG (\code{method=4L}) uses least-squares diagonal preconditioning.
#' This argument also specifies the preconditioner code for Eigen Jacobi SVD.
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
#' The integer \code{ordering} argument takes
#' precedence when supplied.
#' @param threshold
#' Optional threshold for Eigen SVD rank decisions. If \code{NULL}, Eigen's
#' default threshold is used.
#' @param verbose
#' XXX
#' @param rowaug
#' Whether to append diagonal ridge-penalty rows before solving the working
#' least-squares problem. For Eigen LLT, LDLT and CG, \code{NULL} chooses direct
#' Gram construction for dense designs and sparse augmentation for sparse
#' designs, including dense-stored matrices with at least one million entries
#' and fewer than 10 percent nonzeros. Explicit \code{TRUE} or \code{FALSE}
#' overrides this choice. Existing rectangular QR, SVD and LSCG defaults are
#' unchanged. Eligible wide Eigen systems use the Woodbury route instead.
#'
#' @return
#'
#' The function return a list with the required options to solve a linear
#' weighted least squares \code{lmRidge} ridge regression using the given solver
#' and method.
#'
#' @seealso
#'
#' \code{\link{l0glmFit}}
#'
#' @example examples/controlSolverGen.R
#' @export
controlSolverGen <- function(solver,
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
      # NULL selects Gram construction from the working design's sparsity.
      if(is.null(pivot)) pivot <- 0L
      if(is.null(ordering)) ordering <- if (NaturalOrdering) 1L else 0L
      NaturalOrdering <- as.integer(ordering[1L]) == 1L
      l <- list(method = "chol", pivot = as.integer(pivot),
                ordering = as.integer(ordering),
                NaturalOrdering = NaturalOrdering,
                rowaug = rowaug)
      
    } else if (method == 2){
      # Resolve the storage-dependent default in lmridge.eigen.
      if(is.null(pivot)) pivot <- 1L
      if(is.null(ordering)) ordering <- if (NaturalOrdering) 1L else 0L
      NaturalOrdering <- as.integer(ordering[1L]) == 1L
      l <- list(method = "chol", pivot = as.integer(pivot),
                ordering = as.integer(ordering),
                NaturalOrdering = NaturalOrdering,
                rowaug = rowaug)
    } else if (method == 3){
      if (tol_was_null) tol <- 1e-10
      # Resolve the storage-dependent default in lmridge.eigen.
      if(is.null(maxit)) maxit <- 5000L
      if(is.null(precond)) precond <- 2L
      if(!is.numeric(precond) || length(precond) != 1L ||
         !is.finite(precond) || !precond %in% 0:2) {
        stop("Eigen CG precond must be 0 (none), 1 (diagonal), or 2 (incomplete Cholesky).")
      }
      l <- list(method = "CG", tol = tol, maxit = as.integer(maxit), precond = as.integer(precond), rowaug = rowaug)
    } else if (method == 4){
      if(tol_was_null) tol <- 1e-12
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
