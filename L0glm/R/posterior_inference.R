# Conditional inference only: fitting and model selection never call this file.
L0glm_posterior_prior <- function(object,prepared) {
  B <- prepared$B;data <- prepared$data;rows <- prepared$rows;m <- ncol(B)
  scale <- object$preprocessing$scale[rows]
  if(!length(scale))scale <- rep(1,nrow(B))
  working <- B*scale
  lambda <- object$init.lambda
  if(is.null(lambda)||length(lambda)!=1L||!is.finite(lambda)||lambda<0)
    stop('Posterior inference requires a finite scalar fitted penalty.',call.=FALSE)
  flat <- isTRUE(object$support.refit) || lambda==0
  source <- if(flat)'flat prior on the refitted support'else'one inverse-square update at the fitted coefficients'
  if(flat)lambda <- 0
  if(is.null(object[['group']])) {
    grouped <- !is.null(prepared$kind) && isTRUE(object$grouped)
    D <- L0glm_ext_penalty(working,prepared$allowed,data$included,lambda,grouped)
  } else {
    labels <- object[['group']][rows];levels <- unique(labels)
    costs <- numeric(length(levels));location <- match(names(object$group.multiplier),levels)
    costs[location[!is.na(location)]] <- object$group.multiplier[!is.na(location)]
    D <- L0glm_group_penalty(working,prepared$allowed,data$included,lambda,object$grouped,
      list(id=match(labels,levels),multiplier=costs))
  }
  D <- D*scale^2
  prior <- object$specification$prior.mean
  center <- matrix(0,nrow(B),m)
  if(length(prior)==1L)center[] <- prior else if(length(prior)) {
    prior <- as.matrix(prior)
    added <- as.integer(object$added_intercept)
    if(nrow(prior)==nrow(coef(object))-added && added)prior <- rbind(0,prior)
    if(nrow(prior)<nrow(coef(object)) && !is.null(object$interactions))
      prior <- rbind(prior,matrix(0,nrow(coef(object))-nrow(prior),ncol(prior)))
    if(nrow(prior)!=nrow(coef(object)) || !ncol(prior)%in%c(1L,m))
      stop('The stored prior mean does not match posterior coefficient coordinates.',call.=FALSE)
    center[] <- as.numeric(prior[rows,,drop=FALSE])
  }
  if(identical(prepared$kind,'multinomial'))center <- center-rowMeans(center)
  list(precision=D,center=center,source=source)
}

L0glm_posterior_setup <- function(object,control) {
  prepared <- L0glm_inference_data(object,control)
  data <- prepared$data;allowed <- prepared$allowed
  index <- which(allowed);k <- length(index)
  map <- Matrix::Diagonal(k)
  if(identical(prepared$kind,'multinomial')) {
    i <- j <- x <- numeric();column <- 0L
    for(row in seq_len(nrow(allowed))) {
      categories <- which(allowed[row,]);s <- length(categories)
      if(s<2L)next
      contrasts <- L0glm_ext_contrasts(s)
      position <- match(row+(categories-1L)*nrow(allowed),index)
      i <- c(i,rep(position,s-1L));j <- c(j,rep(column+seq_len(s-1L),each=s))
      x <- c(x,as.numeric(contrasts));column <- column+s-1L
    }
    map <- Matrix::sparseMatrix(i=i,j=j,x=x,dims=c(k,column))
  }
  if(ncol(map)>control$max.parameters || max(as.double(k)^2,as.double(k)*control$nsim)>control$max.cells)
    stop('Posterior inference exceeds the parameter or draw-memory limit.',call.=FALSE)
  prior <- L0glm_posterior_prior(object,prepared)
  if(any(!is.finite(prior$precision[index])))stop('A retained coefficient has an undefined adaptive prior.',call.=FALSE)
  prepared$state$phi <- as.numeric(object$dispersion)
  prepared$state$estimated <- rep(FALSE,data$m)
  data$prior <- prior$center
  list(prepared=prepared,data=data,index=index,map=map,prior=prior,
    constrained=if(identical(prepared$kind,'multinomial'))rep(FALSE,ncol(map))else
      rep(data$lower,ncol(allowed))[index]>-1e100)
}

# Negative likelihood score and Fisher information on admitted entries. Cox
# and multinomial use their exact partial/joint likelihood information.
L0glm_posterior_information <- function(context,B) {
  d <- context$data;p <- nrow(B);m <- ncol(B);k <- length(context$index)
  score <- numeric(length(B));H <- matrix(0,k,k)
  positions <- lapply(seq_len(m),function(j)which((context$index-1L)%/%p+1L==j))
  rows <- lapply(seq_len(m),function(j)(context$index[positions[[j]]]-1L)%%p+1L)
  if(d$kind=='cox') {
    stats <- L0glm_ext_cox_stats(d,B)
    score[] <- -as.numeric(Matrix::crossprod(d$Z,stats$score))
    H[,] <- L0glm_ext_cox_information(d$Z,stats)
  } else if(d$kind=='multinomial') {
    probability <- L0glm_ext_softmax(as.matrix(d$Z%*%B)+L0glm_observation_broadcast(d$offset))$p
    total <- rowSums(d$y);w <- d$weights[,1L]
    score[] <- as.numeric(Matrix::crossprod(d$Z,(total*probability-d$y)*w))
    for(a in seq_len(m))for(b in seq_len(m))if(length(rows[[a]])&&length(rows[[b]])) {
      weight <- w*total*probability[,a]*(as.integer(a==b)-probability[,b])
      H[positions[[a]],positions[[b]]] <- as.matrix(Matrix::crossprod(d$Z[,rows[[a]],drop=FALSE],
        Matrix::Diagonal(x=weight)%*%d$Z[,rows[[b]],drop=FALSE]))
    }
  } else for(j in seq_len(m)) {
    eta <- as.numeric(d$Z%*%B[,j])+L0glm_observation_column(d$offset,j)
    family <- d$family;mu <- family$linkinv(eta)
    w <- L0glm_observation_column(d$weights,j);derivative <- family$mu.eta(eta)
    phi <- if(d$kind=='nb')1 else context$prepared$state$phi[j]
    positive <- family$family%in%c('poisson','quasipoisson','Gamma','inverse.gaussian') ||
      d$kind=='nb' || startsWith(family$family,'Negative Binomial(')
    safe <- if(positive)pmax(mu,1e-12)else if(family$family%in%c('binomial','quasibinomial'))
      pmin(1-1e-12,pmax(mu,1e-12))else mu
    if(is.function(family$score)) {
      terms <- family$score(d$y[,j],safe)
      g <- w*terms[,1L]*derivative/phi
      weight <- w*terms[,2L]*derivative^2/phi
    } else {
      variance <- if(d$kind=='nb')safe+safe^2/context$prepared$state$theta[j]else family$variance(safe)
      residual <- (if(positive)safe else mu)-d$y[,j]
      g <- w*residual*derivative/(phi*variance)
      weight <- w*derivative^2/(phi*variance)
    }
    g[w==0] <- 0;weight[w==0] <- 0
    score[(j-1L)*p+seq_len(p)] <- as.numeric(Matrix::crossprod(d$Z,g))
    A <- d$Z[,rows[[j]],drop=FALSE]
    H[positions[[j]],positions[[j]]] <- as.matrix(Matrix::crossprod(A,Matrix::Diagonal(x=weight)%*%A))
  }
  D <- context$prior$precision[context$index];center <- context$prior$center[context$index]
  gradient <- score[context$index]+D*(B[context$index]-center)
  H <- H+diag(D,k)
  list(gradient=as.numeric(Matrix::crossprod(context$map,gradient)),
    precision=as.matrix(Matrix::crossprod(context$map,H%*%context$map)))
}

L0glm_posterior_mode <- function(context,control) {
  prepared <- context$prepared;map <- context$map
  expand <- function(u) {
    B <- prepared$B*0;B[context$index] <- as.numeric(map%*%u);B
  }
  u <- as.numeric(Matrix::crossprod(map,prepared$B[context$index]))
  objective <- function(u,eta=NULL) {
    B <- expand(u)
    value <- suppressWarnings(L0glm_ext_loss(context$data,B,prepared$state))/2+
      sum(context$prior$precision[context$index]*(B[context$index]-context$prior$center[context$index])^2)/2
    if(is.finite(value))value else Inf
  }
  if(!length(u))return(list(mode=u,mean=u,precision=matrix(0,0,0),covariance=matrix(0,0,0),
    iterations=0L,converged=TRUE,score=0))
  converged <- FALSE
  for(iteration in seq_len(control$maxit)) {
    info <- L0glm_posterior_information(context,expand(u));Q <- info$precision;g <- info$gradient
    scale <- 1/sqrt(diag(Q));projected <- g
    projected[context$constrained & u<=0 & g>0] <- 0
    if(max(abs(projected*scale))<max(1e-8,sqrt(control$tol)*.01)) {converged<-TRUE;break}
    proposal <- as.numeric(nuisance_quadratic_cpp(Q*tcrossprod(scale),
      (as.numeric(Q%*%u)-g)*scale,context$constrained))*scale
    step <- L0glm_step(proposal,u,objective,baseline=objective(u),maxit=30L,rtol=1e-12)
    if(!step$accepted)break
    u <- as.numeric(step$beta)
  }
  info <- L0glm_posterior_information(context,expand(u))
  V <- L0glm_inference_inverse(info$precision,control$rank.tol)
  if(any(!is.finite(V)))stop('The conditional posterior information is singular or ill-conditioned.',call.=FALSE)
  projected <- info$gradient;projected[context$constrained & u<=0 & projected>0] <- 0
  score <- max(abs(projected)/sqrt(diag(info$precision)))
  converged <- converged || score<max(1e-8,sqrt(control$tol)*.01)
  if(!converged)stop('The frozen-prior posterior mode did not converge; increase the inference iteration limit.',call.=FALSE)
  # At a constrained mode the gradient is not zero: retain its linear term.
  list(mode=u,mean=u-as.numeric(V%*%info$gradient),precision=info$precision,covariance=V,
    iterations=iteration,converged=converged,score=score)
}

L0glm_posterior_sample <- function(context,mode,nsim,burnin,thin) {
  k <- length(mode$mean)
  if(!k)return(structure(matrix(0,0L,nsim),sampler='empty support'))
  scale <- 1/sqrt(diag(mode$precision));Q <- mode$precision*tcrossprod(scale)
  factor <- chol(Q)
  normal <- function(n)backsolve(factor,matrix(rnorm(k*n),k,n))*scale+mode$mean
  sampler <- 'independent Gaussian'
  if(!any(context$constrained))draws <- normal(nsim)else {
    batch <- min(512L,nsim);probe <- normal(batch)
    valid <- colSums(probe[context$constrained,,drop=FALSE]<0)==0L
    if(mean(valid)>=.1) {
      # Independent rejection avoids slow Gibbs mixing when bounds remove
      # little Gaussian probability, including highly correlated PSF columns.
      draws <- matrix(0,k,nsim);saved <- 0L
      repeat {
        accepted <- probe[,valid,drop=FALSE];take <- min(nsim-saved,ncol(accepted))
        if(take)draws[,saved+seq_len(take)] <- accepted[,seq_len(take),drop=FALSE]
        saved <- saved+take
        if(saved==nsim)break
        probe <- normal(batch);valid <- colSums(probe[context$constrained,,drop=FALSE]<0)==0L
      }
      sampler <- 'independent truncated-Gaussian rejection'
    } else {
      draws <- posterior_truncated_cpp(mode$mean/scale,Q,context$constrained,
        mode$mode/scale,as.integer(nsim),as.integer(burnin),as.integer(thin))*scale
      sampler <- 'truncated-Gaussian Gibbs'
    }
  }
  structure(as.matrix(context$map%*%draws),sampler=sampler)
}
