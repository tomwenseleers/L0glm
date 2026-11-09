L0glm_benchmark_biglasso_matrix <- function(X) {
  if(!inherits(X,'L0glm_data') || X$format!='bigmemory' || X$missing!='error' ||
     anyDuplicated(X$rows) || !identical(X$columns,seq_len(X$dimensions[2])))return(NULL)
  M <- X$object
  if(bigmemory::is.nil(M@address))M <- bigmemory::attach.big.matrix(X$descriptor)
  if(bigmemory::typeof(M)=='double')M else NULL
}

L0glm_benchmark_biglasso <- function(data,grid,options,grouped) {
  X <- data$X
  M <- L0glm_benchmark_biglasso_matrix(X)
  if(!is.null(M))rows <- X$rows else {
    A <- as.matrix(X)
    M <- bigmemory::as.big.matrix(A,type='double');rows <- seq_len(nrow(A))
  }
  y <- data$y
  gaussian <- data$family$family=='gaussian'
  kind <- data$kind
  if(gaussian)y <- y-data$offset
  keep <- data$weights[,1L]>0
  response <- matrix(NA_real_,nrow(M),ncol(y));response[rows,] <- y
  family <- if(kind=='cox')'cox'else if(gaussian && grouped)'mgaussian'else
    if(gaussian)'gaussian'else 'binomial'
  args <- list(X=M,y=if(ncol(y)==1L)as.numeric(response)else response,
    row.idx=rows[keep],penalty='lasso',family=family,nlambda=data$nlambda,
    lambda.min=if(nrow(X)<ncol(X)).01 else 1e-4,ncores=1L,eps=1e-7,max.iter=10000L)
  if(!is.null(options$native.screen)) {
    args$screen <- options$native.screen;options$native.screen <- NULL
  }
  args <- utils::modifyList(args,options)
  if(!is.null(grid))args$lambda <- grid
  fit <- do.call(biglasso::biglasso,args)
  if(any(!is.finite(fit$center))||any(!is.finite(fit$scale)))
    stop('biglasso requires finite predictors without missing values.',call.=FALSE)
  coefficients <- lapply(seq_along(fit$lambda),function(j) {
    B <- if(is.list(fit$beta))do.call(cbind,lapply(fit$beta,function(z)z[,j,drop=FALSE]))else
      fit$beta[,j,drop=FALSE]
    B
  })
  L0glm_benchmark_align(list(grid=fit$lambda,coefficients=coefficients),grid)
}
