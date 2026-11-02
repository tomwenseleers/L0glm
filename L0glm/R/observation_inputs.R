# A shared observation column stays n x 1, including in native multichannel fits.
L0glm_observation_input <- function(value,n,m,default,name) {
  if(is.null(value))value<-default
  if(length(value)==1L)value<-rep(as.numeric(value),n)
  if(is.null(dim(value))) {
    if(length(value)!=n)stop(name," must be scalar, length n, n x 1 or n x m.",call.=FALSE)
    value<-matrix(as.numeric(value),n,1L)
  } else {
    if(nrow(value)!=n||!ncol(value)%in%c(1L,m))
      stop(name," must be scalar, length n, n x 1 or n x m.",call.=FALSE)
    if(!inherits(value,"Matrix")) {value<-as.matrix(value);storage.mode(value)<-"double"}
  }
  entries<-if(inherits(value,"sparseMatrix")) {
    if("x"%in%methods::slotNames(value))value@x else 1
  }else value
  if(any(!is.finite(entries))||(name=="weights"&&any(entries<0)))
    stop(name," must be finite",if(name=="weights")" and nonnegative"else"",".",call.=FALSE)
  value
}
L0glm_observation_column <- function(value,column,rows=seq_len(nrow(value)))
  as.numeric(value[rows,if(ncol(value)==1L)1L else column])
L0glm_observation_subset <- function(value,columns)
  value[,if(ncol(value)==1L)1L else columns,drop=FALSE]
L0glm_observation_broadcast <- function(value)
  if(!is.null(dim(value))&&ncol(value)==1L)as.numeric(value)else
    if(inherits(value,"Matrix"))as.matrix(value)else value
