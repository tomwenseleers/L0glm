# Resolve a fit-wide Gram/factor budget once, before native channel threads start.
L0glm_validate_memory <- function(x) {
  if (!is.null(x) && (!is.numeric(x)||length(x)!=1L||is.na(x)||x<0))
    stop("max.memory must be NULL or one nonnegative number in GiB.",call.=FALSE)
  invisible(x)
}
L0glm_solver_memory <- function(control) {
  L0glm_validate_memory(control$max.memory)
  if (is.null(control$.workspace.total))
    control$.workspace.total <- solver_workspace_budget_cpp(
      if(is.null(control$max.memory))-1 else control$max.memory*1024^3)
  control
}
L0glm_workspace_fits <- function(X,control,gram=TRUE,rows=FALSE,nthreads=1L) {
  control <- L0glm_solver_memory(control)
  X <- if(issparse(X))as(as(X,'generalMatrix'),'CsparseMatrix')else as.matrix(X)
  budget <- control$.workspace.total/max(1L,nthreads)
  L0glm_workspace_estimate(X,budget,gram,rows) <= budget
}
L0glm_workspace_estimate <- function(X,budget,gram=TRUE,rows=FALSE) {
  estimate <- solver_workspace_estimate_cpp(X,gram,rows)
  if(gram && issparse(X) && estimate>budget && budget>0)
    estimate <- solver_sparse_workspace_cpp(X,rows,budget)
  estimate
}
L0glm_require_workspace <- function(X,control,gram=TRUE,rows=FALSE,nthreads=1L,operation="Factorization",
                                    unpenalized=integer()) {
  control <- L0glm_solver_memory(control)
  budget <- control$.workspace.total/max(1L,nthreads)
  X <- if(issparse(X))as(as(X,'generalMatrix'),'CsparseMatrix')else as.matrix(X)
  if(gram && rows && length(unpenalized)) {
    # Woodbury zeroes these columns before forming its observation Gram;
    # their separate Schur complement needs only narrow dense work arrays.
    projected <- X[,-unpenalized,drop=FALSE]
    nuisance <- 48.*nrow(X)*(length(unpenalized)+1L)
    estimate <- L0glm_workspace_estimate(projected,max(0,budget-nuisance),TRUE,TRUE)+nuisance
    if(ncol(projected)<=nrow(X))
      estimate <- max(estimate,L0glm_workspace_estimate(X,budget,TRUE,FALSE))
  } else estimate <- L0glm_workspace_estimate(X,budget,gram,rows)
  if(estimate>budget)
    stop(operation," exceeds the Gram/factor workspace budget; use screening, Eigen CG/LSCG ",
      "with diagonal preconditioning, or increase controlSolverGen(max.memory=...).",call.=FALSE)
  invisible(NULL)
}
L0glm_require_bytes <- function(estimate,control,operation="Factorization") {
  control <- L0glm_solver_memory(control)
  if(estimate>control$.workspace.total)
    stop(operation," exceeds the Gram/factor workspace budget; use screening, a matrix-free solver, ",
      "or increase controlSolverGen(max.memory=...).",call.=FALSE)
  invisible(NULL)
}
L0glm_matrixfree_ridge <- function(X,y,penalty,center,start,control,dual=FALSE) {
  X <- if(issparse(X))as(as(X,'generalMatrix'),'CsparseMatrix')else as.matrix(X)
  result <- ridge_matrixfree_cpp(X,as.numeric(y),penalty,center,start,
    if(is.null(control$tol)||control$tol<=0)1e-10 else control$tol,
    if(is.null(control$maxit))0L else control$maxit,dual)
  if(!isTRUE(result$success))stop("Matrix-free ridge did not reach its residual tolerance.",call.=FALSE)
  as.numeric(result$coefficients)
}
