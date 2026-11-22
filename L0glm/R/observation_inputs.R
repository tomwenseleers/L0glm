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

# Keep the complete linear predictor (including offsets) for the first IRLS step.
L0glm_working_start <- function(mustart, etastart, family, n, m) {
  if (is.null(mustart) && is.null(etastart)) return(NULL)
  on_link <- !is.null(etastart)
  name <- if (on_link) "etastart" else "mustart"
  value <- if (on_link) etastart else mustart
  if (!is.numeric(value) && !inherits(value, "dMatrix"))
    stop(name, " must be numeric.", call. = FALSE)
  value <- as.matrix(L0glm_observation_input(value, n, m, NULL, name))
  kind <- L0glm_ext_kind(family)
  if (identical(kind, "multinomial")) {
    if (ncol(value) != m)
      stop(name, " must have a column for each multinomial category.", call. = FALSE)
    if (!on_link) {
      if (any(value <= 0) || any(abs(rowSums(value) - 1) > 1e-8))
        stop("Multinomial mustart rows must be positive probabilities summing to one.", call. = FALSE)
      value <- log(value)
    }
    return(value - rowMeans(value))
  }
  if (identical(kind, "cox")) {
    if (!on_link && any(value <= 0)) stop("Cox mustart must be a positive relative risk.", call. = FALSE)
    return(if (on_link) value else log(value))
  }
  valid <- function(x, check) all(is.finite(x)) &&
    (!is.function(check) || isTRUE(all(check(x))))
  if (!on_link && !valid(value, family$validmu))
    stop("mustart is outside the family mean domain.", call. = FALSE)
  eta <- if (on_link) value else suppressWarnings(family$linkfun(value))
  mu <- if (on_link) suppressWarnings(family$linkinv(eta)) else value
  if (!valid(eta, family$valideta) || !valid(mu, family$validmu))
    stop(name, " gives an invalid starting mean or linear predictor.", call. = FALSE)
  dim(eta) <- dim(value)
  eta
}

L0glm_start_column <- function(eta, channel)
  if (!is.null(eta)) L0glm_observation_column(eta, channel)

L0glm_start_rows <- function(value, rows) {
  if (is.null(value) || length(value) == 1L) return(value)
  if (is.null(dim(value))) value[rows] else value[rows, , drop = FALSE]
}
