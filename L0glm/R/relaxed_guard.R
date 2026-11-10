# Structural information checks for special-family unpenalized refits.
# Built-in ordinary-family/native passes share the C++ counterpart.
L0glm_relaxed_gram <- function(G,factor=NULL) {
 if(!ncol(G))return(TRUE)
 norm <- sqrt(pmax(0,diag(G)))
 if(any(!is.finite(norm)|norm==0))return(FALSE)
 if(is.null(factor)) {
  G <- as.matrix(G)/(norm%o%norm)
  factor <- tryCatch(chol(G),error=function(e)NULL)
  norm[] <- 1
 }
 if(is.null(factor))return(FALSE)
 # U / column_norm is the Cholesky factor of the normalized Gram matrix.
 pivots <- diag(factor)/norm
 all(is.finite(pivots))&&min(pivots)>sqrt(ncol(G)*.Machine$double.eps)*max(pivots)
}
L0glm_relaxed_standard <- function(X,weights,allowed,estimated=FALSE,nuisance_df=0L) {
 m <- ncol(allowed);estimated <- rep_len(estimated,m)
 last_rows <- last_columns <- NULL;last_value <- FALSE
 vapply(seq_len(m),function(k) {
  w <- if(is.null(weights))rep(1,nrow(X))else if(is.null(dim(weights)))rep_len(as.numeric(weights),nrow(X))else
    L0glm_observation_column(weights,k)
  rows <- which(w>0&is.finite(w));columns <- which(allowed[,k])
  if(length(columns)>length(rows)||estimated[k]&&length(columns)+nuisance_df>=length(rows))return(FALSE)
  if(identical(rows,last_rows)&&identical(columns,last_columns))return(last_value)
  value <- L0glm_relaxed_gram(Matrix::crossprod(X[rows,columns,drop=FALSE]))
  last_rows <<- rows;last_columns <<- columns;last_value <<- value;value
 },TRUE)
}
L0glm_relaxed_check <- function(data,B,allowed) {
 if(data$kind=='glm')return(rep(TRUE,data$m)) # checked by the native driver
 if(data$kind=='cox') {
  columns <- which(allowed[,1L]);stats <- L0glm_ext_cox_stats(data,B)
  rows <- which(stats$q>0)
  if(length(columns)>length(rows)-1L)return(FALSE)
  if(!length(columns))return(TRUE)
  A <- data$Z[rows,columns,drop=FALSE];w <- stats$q[rows]
  sum_x <- as.numeric(Matrix::crossprod(A,w))
  # All event risk sets are nested. Information is nonsingular exactly when
  # selected predictors are independent modulo a constant on their union.
  G <- Matrix::crossprod(A,Matrix::Diagonal(x=w)%*%A)-tcrossprod(sum_x)/sum(w)
  return(L0glm_relaxed_gram(G))
 }
 if(data$kind=='multinomial') {
  rows <- which(data$weights[,1L]>0 & rowSums(data$y)>0)
  columns <- which(rowSums(allowed)>=2L);C <- data$m
  if(!length(columns))return(TRUE)
  A <- data$Z[rows,columns,drop=FALSE];counts <- rowSums(allowed[columns,,drop=FALSE])-1L
  if(sum(counts)>length(rows)*(C-1L))return(FALSE)
  G <- as.matrix(Matrix::crossprod(A))
  if(all(counts==C-1L))return(L0glm_relaxed_gram(G))
  # Entrywise CLR support has a different zero-sum plane in each row.
  basis <- lapply(columns,function(j) {
   keep <- which(allowed[j,]);Q <- matrix(0,C,length(keep)-1L)
   Q[keep,] <- L0glm_ext_contrasts(length(keep));Q
  })
  ends <- cumsum(counts);starts <- c(1L,head(ends,-1L)+1L)
  H <- matrix(0,sum(counts),sum(counts))
  for(j in seq_along(columns))for(k in seq_len(j)) {
   a <- starts[j]:ends[j];b <- starts[k]:ends[k]
   H[a,b] <- G[j,k]*crossprod(basis[[j]],basis[[k]])
   H[b,a] <- t(H[a,b,drop=FALSE])
  }
  return(L0glm_relaxed_gram(H))
 }
 L0glm_relaxed_standard(data$Z,data$weights,allowed,data$response$estimated)
}
