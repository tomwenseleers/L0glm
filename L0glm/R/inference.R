L0glm_inference_control <- function(control) {
  defaults <- list(maxit=100L,tol=1e-10,level=.95,max.parameters=1000L,
    rank.tol=1e-10,nboot=199L,seed=NULL,max.cells=5e7,refit.args=list(),refit=NULL,
    nsim=2000L,burnin=500L,thin=1L)
  if(!is.list(control) || length(control) && (is.null(names(control)) ||
     any(!names(control)%in%names(defaults)) || anyDuplicated(names(control))))
    stop('Use named control.inference entries.',call.=FALSE)
  defaults[names(control)] <- control
  for(name in c('maxit','tol','max.parameters','rank.tol','nboot','max.cells'))
    if(length(defaults[[name]])!=1L || !is.numeric(defaults[[name]]) ||
       !is.finite(defaults[[name]]) || defaults[[name]]<=0)
      stop('Inference limits and tolerances must be positive finite values.',call.=FALSE)
  if(length(defaults$level)!=1L || !is.finite(defaults$level) || defaults$level<=0 || defaults$level>=1 ||
     defaults$nboot<2L || defaults$nboot!=floor(defaults$nboot) || !is.list(defaults$refit.args) ||
     !is.null(defaults$refit) && !is.function(defaults$refit))
    stop('Use 0 < level < 1, at least two bootstraps and a refit.args list.',call.=FALSE)
  if(!is.null(defaults$seed) && (length(defaults$seed)!=1L || !is.numeric(defaults$seed) ||
    !is.finite(defaults$seed) || abs(defaults$seed)>.Machine$integer.max))
    stop('seed must be NULL or a finite integer-range number.',call.=FALSE)
  for(name in c('nsim','burnin','thin')) {
    value <- defaults[[name]]
    minimum <- switch(name,nsim=2L,burnin=0L,thin=1L)
    if(length(value)!=1L || !is.numeric(value) || !is.finite(value) || value<minimum ||
       value!=floor(value) || value>.Machine$integer.max)
      stop('Use nsim >= 2, burnin >= 0 and thin >= 1 as integer controls.',call.=FALSE)
  }
  defaults
}

L0glm_inference_inverse <- function(H,tol) {
  if(!nrow(H))return(H)
  H <- (H+t(H))/2
  if(any(!is.finite(H)))return(matrix(NA_real_,nrow(H),ncol(H)))
  scale <- sqrt(pmax(diag(H),0))
  if(any(scale==0))return(H*NA_real_)
  G <- H/tcrossprod(scale)
  R <- tryCatch(chol(G),error=function(e)NULL)
  if(is.null(R) || rcond(G)<tol)return(H*NA_real_)
  chol2inv(R)/tcrossprod(scale)
}

L0glm_inference_data <- function(object,control) {
  B <- coef(object);m <- ncol(B)
  rows <- sort(unique(c(which(Matrix::rowSums(B!=0)>0),object$always.include)))
  A <- object$X[,rows,drop=FALSE];B <- as.matrix(B[rows,,drop=FALSE])
  family <- object$family;kind <- L0glm_ext_kind(family)
  included <- match(object$always.include,rows)
  allowed <- B!=0;allowed[included,] <- TRUE
  if(identical(kind,'multinomial') && object$grouped)allowed[] <- TRUE
  dimension <- sum(allowed)-if(identical(kind,'multinomial'))nrow(B)else 0L
  if(dimension>control$max.parameters)
    stop('Selected-support inference exceeds control.inference$max.parameters.',call.=FALSE)
  nonnegative <- identical(object$engine,'l0nnPois') || isTRUE(object$specification$nonnegative)
  if(identical(kind,'multinomial'))nonnegative <- FALSE
  weights <- if(is.null(kind))object$prior.weights else object$weights
  data <- L0glm_ext_prepare(A,object$y,family,weights,object$offset,FALSE,
    included,nonnegative,'none',B,NULL,if(is.null(object$control.family))list()else object$control.family,
    if(is.null(kind))'glm'else kind)
  data$allowed <- allowed
  phi <- rep_len(object$dispersion,m)
  estimated <- rep_len(family$family%in%c('gaussian','Gamma','inverse.gaussian',
    'quasi','quasibinomial','quasipoisson'),m)
  if(!is.null(object$dispersion.estimated))estimated <- rep_len(object$dispersion.estimated,m)
  state <- if(is.null(kind)) {
    initial <- L0glm_group_glm_initialize(data,ifelse(estimated,-phi,phi),0)
    data <- initial$data;initial$state
  } else if(kind=='nb')list(theta=object$theta,estimated=object$theta.estimated)else list()
  list(data=data,rows=rows,B=B,allowed=allowed,state=state,kind=kind)
}

L0glm_inference_multinomial <- function(data,B,allowed,control) {
  D <- matrix(Inf,nrow(B),ncol(B));D[allowed] <- 0
  map <- L0glm_ext_multinomial_basis(D,B,B*0)$map
  probabilities <- L0glm_ext_softmax(as.matrix(data$Z%*%B)+
    L0glm_observation_broadcast(data$offset))$p
  p <- nrow(B);C <- ncol(B);H <- matrix(0,p*C,p*C)
  weight <- data$weights[,1L]*rowSums(data$y)
  for(a in seq_len(C))for(b in seq_len(C)) {
    w <- weight*probabilities[,a]*(as.integer(a==b)-probabilities[,b])
    H[(a-1L)*p+seq_len(p),(b-1L)*p+seq_len(p)] <-
      as.matrix(Matrix::crossprod(data$Z,Matrix::Diagonal(x=w)%*%data$Z))
  }
  reduced <- as.matrix(Matrix::crossprod(map,H%*%map))
  V <- L0glm_inference_inverse(reduced,control$rank.tol)
  as.matrix(map%*%V%*%Matrix::t(map))
}

L0glm_conditional_inference <- function(object,control) {
  prepared <- L0glm_inference_data(object,control)
  data <- prepared$data;B <- prepared$B;allowed <- prepared$allowed;kind <- prepared$kind
  args <- list(solver='eigen',method=1L,control.solver=controlSolverGen('eigen',1L),validation=NULL)
  fit <- if(length(B))L0glm_ext_pass(data,B,allowed,prepared$state,0,FALSE,0,
    as.integer(control$maxit),control$tol,args,TRUE)else
      list(B=B,state=prepared$state,converged=TRUE,iterations=0L)
  if(length(B) && any(data$lower>-1e100))fit <- L0glm_ext_boundary(data,fit,NULL)
  B <- fit$B;m <- ncol(B);p <- nrow(B)
  se <- matrix(NA_real_,p,m);covariance <- vector('list',if(identical(kind,'multinomial'))1L else m)
  df <- rep(Inf,m)
  if(identical(kind,'multinomial')) {
    covariance[[1L]] <- L0glm_inference_multinomial(data,B,allowed,control)
    se[] <- sqrt(pmax(0,diag(covariance[[1L]])))
  } else if(identical(kind,'cox')) {
    stats <- L0glm_ext_cox_stats(data,B)
    H <- L0glm_ext_cox_information(data$Z,stats)
    covariance[[1L]] <- L0glm_inference_inverse(H,control$rank.tol)
    se[] <- sqrt(pmax(0,diag(covariance[[1L]])))
  } else for(j in seq_len(m)) {
    keep <- which(allowed[,j]);A <- data$Z[,keep,drop=FALSE]
    eta <- as.numeric(data$Z%*%B[,j])+L0glm_observation_column(data$offset,j)
    family <- if(identical(kind,'nb'))stats::make.link(data$family$link)else data$family
    mu <- family$linkinv(eta);w <- L0glm_observation_column(data$weights,j)
    if(is.function(data$family$score)) {
      loss <- function(beta) {
        mu <- data$family$linkinv(as.numeric(A%*%beta)+L0glm_observation_column(data$offset,j))
        -sum(w[w>0]*data$family$log_density(data$y[,j],mu)[w>0])
      }
      H <- tryCatch(stats::optimHess(B[keep,j],loss),
        error=function(e)matrix(NA_real_,length(keep),length(keep)))
    } else {
      variance <- if(identical(kind,'nb'))mu+mu^2/fit$state$theta[j]else family$variance(mu)
      phi <- if(identical(kind,'nb'))1 else fit$state$phi[j]
      weight <- w*family$mu.eta(eta)^2/variance/phi;weight[w==0] <- 0
      H <- as.matrix(Matrix::crossprod(A,Matrix::Diagonal(x=weight)%*%A))
    }
    covariance[[j]] <- L0glm_inference_inverse(H,control$rank.tol)
    se[keep,j] <- sqrt(pmax(0,diag(covariance[[j]])))
    if(is.null(kind) && fit$state$estimated[j])df[j] <- sum(w>0)-length(keep)
  }
  boundary <- allowed & matrix(data$lower>-1e100,p,m) & B==0
  se[boundary] <- NA_real_
  index <- which(allowed,arr.ind=TRUE)
  row_names <- rownames(coef(object));if(is.null(row_names))row_names <- as.character(seq_len(nrow(coef(object))))
  unnamed <- is.na(row_names)|!nzchar(row_names)
  row_names[unnamed] <- paste0('x',which(unnamed)-as.integer(object$added_intercept))
  channel_names <- colnames(coef(object));if(is.null(channel_names))channel_names <- as.character(seq_len(m))
  if(identical(kind,'multinomial')) {
    names <- as.vector(outer(row_names[prepared$rows],channel_names,paste,sep=':'))
    dimnames(covariance[[1L]]) <- list(names,names)
  } else for(j in seq_len(m)) {
    names <- row_names[prepared$rows[allowed[,j]]]
    dimnames(covariance[[j]]) <- list(names,names)
  }
  degrees <- df[index[,2L]]
  critical <- rep(NA_real_,length(degrees));valid <- degrees>0
  critical[valid] <- stats::qt((1+control$level)/2,degrees[valid])
  statistic <- B[index]/se[index]
  table <- data.frame(term=row_names[prepared$rows[index[,1L]]],channel=channel_names[index[,2L]],
    selected.estimate=as.numeric(coef(object)[cbind(prepared$rows[index[,1L]],index[,2L])]),
    estimate=B[index],std.error=se[index],statistic,
    p.value=2*stats::pt(-abs(statistic),degrees),df=degrees,
    conf.low=B[index]-critical*se[index],conf.high=B[index]+critical*se[index],
    boundary=boundary[index],row=prepared$rows[index[,1L]],column=index[,2L])
  entries <- which(B!=0,arr.ind=TRUE)
  beta <- Matrix::sparseMatrix(i=prepared$rows[entries[,1L]],j=entries[,2L],x=B[entries],
    dims=dim(coef(object)),dimnames=dimnames(coef(object)))
  list(method='conditional',coefficients=table,refit=beta,covariance=covariance,
    covariance.rows=prepared$rows,converged=fit$converged,iterations=fit$iterations,
    level=control$level,dispersion=fit$state$phi,theta=fit$state$theta,conditional.on.support=TRUE)
}
