#' @name controlSolverGen
#' @title Set Options of Constrained Ridge Regularised Weighted Least Squares Regression Solvers
#' 
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
#' incomplete Cholesky preconditioning by default, \code{4L}=LSCG with
#' incomplete Cholesky preconditioning by default, \code{5L}=column-pivoted QR, and
#' \code{6L}=BDCSVD for dense systems or Jacobi SVD for sparse systems.
#' @param maxit
#' The maximum number of iterations.
#' @param standardize
#' Whether \code{bigGlm} standardizes columns of the augmented working design.
#' Returned coefficients retain the original scale.
#' @param reltol
#' Relative convergence tolerance for OSQP. The default is \code{1e-9}.
#' @param tol
#' Relative convergence tolerance for Eigen iterative solvers. For
#' CG the default is \code{1e-10}; LSCG defaults to \code{1e-12}. LSCG
#' checks the normal-equation residual for every right-hand side by
#' recomputing it from the returned coefficients. A solve that does not
#' reach the requested accuracy uses the Cholesky LLT fallback.
#' @param abstol
#' Absolute tolerance for OSQP (default \code{1e-9}); coordinate-update stopping threshold for
#' \code{bigGlm}. Tall and square \code{bigGlm} working systems also receive
#' a Cholesky accuracy polish with a fresh projected-gradient check; see
#' \code{\link{lmRidge}}.
#' @param block_size
#' Optional block size used by solvers that support blocked updates.
#' @param precond
#' For Eigen CG and LSCG (\code{method=3L} or \code{method=4L}),
#' \code{0L} selects no preconditioning,
#' \code{1L} diagonal preconditioning and \code{2L} incomplete Cholesky
#' (the default for both dense and sparse systems). Dense CG retains its dense
#' matrix operator; only the incomplete Cholesky preconditioner uses a sparse
#' copy containing every nonzero entry. For fully populated systems this
#' factorization can cost more time and memory than diagonal preconditioning.
#' For LSCG, incomplete Cholesky factors the normal matrix of the working
#' least-squares design; \code{1L} uses its least-squares diagonal instead.
#' This argument also specifies the preconditioner code for Eigen Jacobi SVD.
#' @param pivot
#' Optional integer pivoting code for Eigen QR and Cholesky methods. Defaults
#' are method-specific.
#' @param ordering
#' Optional integer sparse ordering code for Eigen sparse QR and Cholesky
#' methods. For sparse QR, \code{0L}=AMD, \code{1L}=COLAMD and
#' \code{2L}=NaturalOrdering, with \code{1L} the default. Rectangular sparse QR
#' uses COLAMD when AMD is requested. For sparse Cholesky,
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
#' @param max.memory
#' Optional Gram/factor workspace budget in GiB, shared across concurrent native
#' channels. \code{NULL} uses 25 percent of currently available physical memory;
#' \code{Inf} disables the guard. Estimates include working copies and factors;
#' sparse fill depends on the design. Oversized incomplete Cholesky uses diagonal
#' preconditioning. CG/LSCG can use matrix-free ridge products; direct solves and
#' LLT fallbacks require their estimated workspace to fit.
#' A failed coarse sparse check is refined using the unique structural Gram
#' entries and symbolic Cholesky fill, with allowances for work arrays and
#' copies. The probe is itself bounded by the budget. Repeated overlaps in a
#' banded convolution design therefore do not count as separate Gram entries.
#' The Poisson backend also excludes zero-scaled unpenalized columns from its
#' Woodbury Gram estimate and includes their separate Schur-complement work.
#' Allocation checks also cover dispersion pilots and serial special-family
#' solver adapters.
#' @param nonnegative.method
#' Nonnegative working-solve strategy: \code{"clamp"} solves without coefficient
#' bounds and then clips to them; \code{"constrained"} enforces bounds within
#' the weighted ridge solve. \code{NULL} preserves solver-specific defaults:
#' \code{"clamp"} for Eigen and \code{"constrained"} for OSQP and \code{bigGlm}.
#' For Eigen, constrained mode uses a checked active-set Cholesky polish.
#' LLT, LDLT and CG reuse eligible predictor Gram matrices. LSCG also corrects
#' rectangular and wide proposals, forming a predictor Gram matrix on demand
#' when bounds are violated. This requires room under \code{max.memory}; an
#' insufficient budget is an error. Other rectangular and Woodbury routes
#' retain clamping. An unsuccessful numerical polish retains the clamped
#' proposal. Unrestricted intercepts remain free.
#' Unconstrained OSQP steps use infeasibility tolerances of \code{1e-12}, unless
#' explicitly supplied in the control list.
#'
#' @details
#' For Eigen with \code{nonnegative.method="constrained"}, the polish solves the current weighted
#' ridge quadratic on successive feasible coefficient sets. Bounds, prior means
#' and column normalization are retained, including an unrestricted binomial
#' intercept. Batched bound updates use Cholesky on the free set; cycling or an
#' unsuccessful batch uses a sequential active-set solve. A fresh equilibrated
#' box KKT residual must be at most \code{1e-10}, and the quadratic must not
#' increase relative to the clipped proposal. An unsuccessful correction keeps
#' that proposal. Dense and sparse Gram matrices retain their storage type.
#' All native channel work remains compatible with OpenMP.
#'
#' LLT, LDLT and CG reuse a Gram matrix from a tall or square active working
#' design. LSCG retains its least-squares operator for the initial solve and,
#' only when needed, constructs the equivalent weighted ridge Gram system for
#' the constrained correction. Its allocation is checked against the workspace
#' budget. A feasible unconstrained proposal needs no correction. Rectangular
#' QR/SVD and the other observation-space Woodbury paths keep their usual solves.
#' This optional correction can recover
#' peaks lost through clipping in correlated nonnegative designs, but can also
#' change the selected support and increase fitting time. Eigen defaults to \code{"clamp"}.
#'
#' OSQP working systems are diagonally equilibrated, with bounds and starting
#' coefficients transformed consistently. Returned coefficients retain their
#' original units. Tall and square working designs also receive the same checked
#' constrained Cholesky polish as \code{bigGlm}; an unsuccessful polish retains
#' the OSQP coefficients. See \code{\link{lmRidge}}.
#'
#' LSCG retains the working least-squares operator \eqn{A}, including any ridge
#' rows, and uses incomplete Cholesky of \eqn{A^T A} as its default
#' preconditioner. Every nonzero normal-matrix entry is retained when constructing
#' the sparse preconditioner. For already formed Gram or Woodbury systems
#' \eqn{G\beta=h}, this setting uses CG directly on \eqn{G}, avoiding a second
#' squaring of the condition number.
#'
#' Before each incomplete-Cholesky LSCG solve, the previous-solution guess is
#' compared with zero. Zero is used when
#' \eqn{\|A^T(b-A\beta_{\mathrm{previous}})\|_2 > \|A^T b\|_2}.
#' This guards against poorly scaled inner guesses when adaptive ridge weights
#' change. The public starting coefficients still determine model initialization
#' and feature eligibility; penalties and prior means retain their specified values.
#'
#' Each LSCG right-hand side is checked using the freshly computed relative
#' normal residual \eqn{\|A^T(b-A\beta)\|_2/\|A^T b\|_2}. If recursive
#' convergence fails this check, incomplete-Cholesky LSCG can restart once from
#' its returned coefficients, reusing the same factorization and the remaining
#' iteration budget. A zero normal right-hand side passes only with a zero
#' normal residual. Inexact or nonfinite results use the Cholesky LLT fallback.
#' Forming the normal matrix for incomplete Cholesky can require quadratic
#' storage for fully populated designs; diagonal preconditioning avoids that
#' allocation in rectangular LSCG solves.
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
                               rowaug = NULL,
                               nonnegative.method = NULL, max.memory = NULL){
  L0glm_validate_memory(max.memory)
  if (is.null(nonnegative.method))
    nonnegative.method <- if (solver == "eigen") "clamp" else "constrained"
  if (!is.character(nonnegative.method) || length(nonnegative.method) != 1L ||
      is.na(nonnegative.method) || !nonnegative.method %in% c("clamp", "constrained"))
    stop('nonnegative.method must be "clamp" or "constrained".', call. = FALSE)

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
    
    if(is.null(abstol)) abstol <- 1e-9
    if(is.null(reltol)) reltol <- 1e-9
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
    } else if (method %in% c(3, 4)){
      type <- if (method == 3) "CG" else "LSCG"
      if (tol_was_null) tol <- if (method == 3) 1e-10 else 1e-12
      if (method == 4 && is.null(rowaug)) rowaug <- TRUE
      if(is.null(maxit)) maxit <- 5000L
      if(is.null(precond)) precond <- 2L
      if(!is.numeric(precond) || length(precond) != 1L ||
         !is.finite(precond) || !precond %in% 0:2) {
        stop(sprintf("Eigen %s precond must be 0 (none), 1 (diagonal), or 2 (incomplete Cholesky).", type))
      }
      l <- list(method = type, tol = tol, maxit = as.integer(maxit), precond = as.integer(precond), rowaug = rowaug)
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
  l$nonnegative.method <- nonnegative.method
  l$max.memory <- max.memory
  return(l)
}
