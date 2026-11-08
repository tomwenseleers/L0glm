# A lazy augmented design stores the input once and builds requested products.
L0glm_design_view <- function(X,terms=list(),intercept=FALSE,names=NULL,check=FALSE) {
  if(is.null(names)) {
    original <- colnames(X)
    if(is.null(original))original <- paste0('x',seq_len(ncol(X)))
    names <- make.unique(c(original,vapply(terms,function(i)paste(original[i],collapse=':'),'')))
  }
  structure(list(data=X,terms=terms,intercept=intercept,columns=names,check=check),
    class='L0glm_design')
}

#' @export
#' @method dim L0glm_design
#' @noRd
dim.L0glm_design <- function(x)c(nrow(x$data),length(x$columns)+as.integer(x$intercept))

#' @export
#' @method dimnames L0glm_design
#' @noRd
dimnames.L0glm_design <- function(x)list(rownames(x$data),
  if(x$intercept)c('Intercept',x$columns)else x$columns)

#' @export
#' @method print L0glm_design
#' @noRd
print.L0glm_design <- function(x,...) {
  cat('Lazy training design:',nrow(x),'rows,',ncol(x$data),'original predictors and',
    length(x$terms),'interaction terms.\n')
  invisible(x)
}

#' @export
#' @method [ L0glm_design
#' @noRd
`[.L0glm_design` <- function(x,i,j,...,drop=TRUE) {
  rows <- if(missing(i))seq_len(nrow(x))else {
    index <- seq_len(nrow(x));names(index) <- rownames(x$data);index[i]
  }
  columns <- if(missing(j))seq_len(ncol(x))else {
    index <- seq_len(ncol(x));names(index) <- colnames(x);index[j]
  }
  if(anyNA(columns))stop('Unknown augmented-design columns.',call.=FALSE)
  p <- ncol(x$data);added <- as.integer(x$intercept)
  original <- columns-added
  product <- original[original>p]-p
  needed <- sort(unique(c(original[original>0 & original<=p],unlist(x$terms[product]))))
  A <- L0glm_data_matrix(x$data[rows,needed,drop=FALSE])
  sparse <- inherits(A,'sparseMatrix')
  if(!sparse)A <- as.matrix(A)
  values <- lapply(original,function(k) {
    if(k==0L) {
      value <- matrix(1,length(rows),1L)
      if(sparse)Matrix::Matrix(value,sparse=TRUE)else value
    } else if(k<=p)A[,match(k,needed),drop=FALSE] else {
      factors <- match(x$terms[[k-p]],needed)
      Reduce(`*`,lapply(factors,function(q)A[,q,drop=FALSE]))
    }
  })
  result <- if(length(values))do.call(cbind,values)else L0glm_data_matrix(x$data[rows,integer(),drop=FALSE])
  dimnames(result) <- list(rownames(x$data)[rows],colnames(x)[columns])
  if(x$check && any(!is.finite(if(inherits(result,'sparseMatrix'))result@x else result)))
    stop('Selected predictors and their interaction products must be finite.',call.=FALSE)
  if(drop)base::drop(result)else result
}

L0glm_interaction_newdata <- function(X,recipe) {
  p <- recipe$p.original
  if(is.data.frame(X))X <- as.matrix(X)
  if(is.null(dim(X)) && is.numeric(X))X <- if(p==1L)matrix(X,ncol=1L)else
    matrix(X,nrow=1L,dimnames=list(NULL,names(X)))
  if(!(is.matrix(X)||inherits(X,'Matrix')||inherits(X,'L0glm_data')) || ncol(X)!=p)
    stop('Supply all original main-effect predictors for interaction prediction.',call.=FALSE)
  expected <- recipe$input.names;provided <- colnames(X)
  if(!is.null(expected) && !is.null(provided)) {
    if(!anyNA(expected) && all(nzchar(expected)) && !anyDuplicated(expected)) {
      if(anyNA(provided)||anyDuplicated(provided)||!setequal(expected,provided))
        stop('Prediction column names must match the original main effects.',call.=FALSE)
      X <- X[,match(expected,provided),drop=FALSE]
    } else if(!identical(expected,provided))stop('Ambiguous predictor names require training order.',call.=FALSE)
  }
  if(is.matrix(X) && !(is.numeric(X)||is.logical(X)))stop('Predictors must be numeric.',call.=FALSE)
  L0glm_design_view(X,recipe$terms,names=recipe$columns,check=TRUE)
}

L0glm_criterion_size <- function(p,control) {
  value <- control$.ic.p
  if(is.null(value))return(p)
  if(length(value)!=1L || !is.numeric(value) || !is.finite(value) ||
     value<p || value!=floor(value))stop('Invalid original criterion dimension.',call.=FALSE)
  value
}
