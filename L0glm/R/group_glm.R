# Ordinary-family adapter for the shared grouped-penalty IRLS orchestration.
L0glm_group_glm_initialize <- function(data,dispersion,lambda) {
  phi <- L0glm_dispersion(dispersion,data$family,data$m)
  estimate <- if(is.null(phi))rep(TRUE,data$m)else phi<0
  mu <- matrix(0,data$n,data$m)
  B <- matrix(0,data$p,data$m)
  prior <- L0glm_ext_original_coefficients(data,data$prior)
  for(j in seq_len(data$m)) {
    initialize <- new.env(parent=environment())
    initialize$y <- data$y[,j];initialize$weights <- L0glm_observation_column(data$weights,j)
    initialize$nobs <- data$n;initialize$nvars <- data$p
    initialize$mustart <- initialize$etastart <- NULL
    eval(data$family$initialize,initialize)
    mu[,j] <- initialize$mustart
    if(!data$explicit) {
      beta <- L0glm_initial_coefficients(data$X,data$y[,j],mu[,j],
        L0glm_observation_column(data$weights,j),L0glm_observation_column(data$offset,j),
        prior[,j],data$family,lambda,data$included,
        data$lower/data$preprocessing$scale,data$preprocessing)
      B[,j] <- L0glm_ext_working_coefficients(data,matrix(beta,data$p,1L))[,1L]
    }
  }
  if(is.null(phi))phi <- vapply(seq_len(data$m),function(j) {
    allowed <- which(data$allowed[,j])
    if(!length(allowed))return(L0glm_pearson_scale(data$y[,j],
      data$family$linkinv(L0glm_observation_column(data$offset,j)),
      L0glm_observation_column(data$weights,j),data$family,0))
    L0glm_initial_dispersion(data$X[,allowed,drop=FALSE],data$y[,j],mu[,j],
      L0glm_observation_column(data$weights,j),L0glm_observation_column(data$offset,j),
      data$family,match(data$included,allowed),L0glm_subset_preprocessing(data$preprocessing,allowed))
  },1)
  phi <- abs(phi);phi[!is.finite(phi)|phi<=0] <- 1
  data$initial.mu <- mu
  cacheable <- !is.null(L0glm_native_family(data$family))
  if((data$n+data$p)*data$m<=1e6 && cacheable)
    data$prediction.cache <- new.env(parent=emptyenv())
  if(data$family$family=='gaussian' && data$family$link=='identity' &&
     is.matrix(data$Z) && cacheable)
    data$gram <- new.env(parent=emptyenv())
  list(data=data,B=pmax(B,data$lower),state=list(phi=phi,estimated=estimate,stable=integer(data$m)))
}

L0glm_group_glm_working <- function(data,B,D,state,args) {
  prediction <- L0glm_ext_nb_mean(data,B)
  answer <- B; working <- L0glm_working_model(data$family)
  for(j in seq_len(data$m)) {
    mu <- prediction$mu[,j];eta <- prediction$eta[,j]
    if(any(!is.finite(mu)) || !isTRUE(data$family$validmu(mu)) ||
       is.function(data$family$valideta) && !isTRUE(data$family$valideta(eta))) {
      mu <- data$initial.mu[,j];eta <- data$family$linkfun(mu)
    }
    point <- working(data$y[,j],eta,mu,L0glm_observation_column(data$weights,j),
      L0glm_observation_column(data$offset,j))
    good <- point$good;options <- args
    active <- which(is.finite(D[,j]))
    if(!is.null(data$gram) && args$solver=='eigen' && args$method==1L &&
       !isTRUE(args$control.solver$rowaug) && length(active)<=length(good) &&
       length(active)<=2000L && length(active)^2*data$m<=8e6 && length(active)>0L) {
      key <- as.character(j);cached <- data$gram[[key]]
      if(is.null(cached) || !all(active%in%cached$active)) {
        Xw <- Xweights(data$Z[good,active,drop=FALSE],L0glm_observation_column(data$weights,j,good))
        yw <- Yweights(data$y[good,j]-L0glm_observation_column(data$offset,j,good),
          L0glm_observation_column(data$weights,j,good))
        cached <- list(xx=crossprod(Xw),xy=crossprod(Xw,yw),active=active)
        data$gram[[key]] <- cached
      }
      k <- match(active,cached$active)
      options$control.solver$.crossprod <- list(xx=cached$xx[k,k,drop=FALSE]/state$phi[j],
        xy=cached$xy[k,,drop=FALSE]/state$phi[j])
    }
    answer[,j] <- L0glm_ext_ridge(data$Z[good,,drop=FALSE],point$response[good],
      point$weights[good]/state$phi[j],D[,j],B[,j],data$lower,data$prior[,j],options)
  }
  answer
}

L0glm_ext_loss <- function(data,B,state) {
  if(data$kind!='glm')return(-2*L0glm_ext_loglik(data,B,state))
  eta <- L0glm_ext_eta(data,B)
  sum(vapply(seq_len(data$m),function(j)L0glm_irls_deviance(data$y[,j],eta[,j],
    L0glm_observation_column(data$weights,j),data$family)/state$phi[j],1))
}

L0glm_group_scale_update <- function(data,B,previous,allowed,state,iteration,relaxed) {
  if(!any(state$estimated))return(state)
  mu <- L0glm_ext_nb_mean(data,B)$mu
  for(j in which(state$estimated)) {
    changed <- any((B[,j]!=0)!=(previous[,j]!=0))
    if(changed)state$stable[j] <- 0L
    if(state$stable[j]>=3L && !relaxed)next
    weights <- L0glm_observation_column(data$weights,j)
    k <- sum(allowed[,j]);N <- sum(weights>0)
    if((iteration>=3L && k<.7*N)||relaxed) {
      target <- L0glm_pearson_scale(data$y[,j],mu[,j],weights,data$family,k)
      if(is.finite(target) && target>0) {
        state$stable[j] <- if(!changed && abs(target/state$phi[j]-1)<.01)state$stable[j]+1L else 0L
        state$phi[j] <- max(state$phi[j]/2,min(2*state$phi[j],target))
      }
    }
  }
  state
}

L0glm_group_glm_result <- function(data,result,passes,lambda,grouped,relaxed,args) {
  beta <- L0glm_ext_original_coefficients(data,result$B)
  dimnames(beta) <- list(colnames(data$X),colnames(data$y))
  convergence <- do.call(rbind,lapply(names(passes),function(stage)
    data.frame(stage,iterations=passes[[stage]]$iterations,converged=passes[[stage]]$converged)))
  prediction <- L0glm_ext_nb_mean(data,result$B)
  W <- z <- matrix(0,data$n,data$m);working <- L0glm_working_model(data$family)
  for(j in seq_len(data$m)) {
    mu <- prediction$mu[,j]
    if(!isTRUE(data$family$validmu(mu)))mu <- data$initial.mu[,j]
    point <- working(data$y[,j],prediction$eta[,j],mu,
      L0glm_observation_column(data$weights,j),L0glm_observation_column(data$offset,j))
    W[,j] <- point$weights;z[,j] <- point$response
  }
  W[is.infinite(W)] <- 1e11
  answer <- structure(list(coefficients=Matrix::Matrix(beta,sparse=TRUE),X=data$X,y=data$y,
    prior.weights=data$weights,offset=data$offset,family=data$family,
    weights=W,z=z,
    dispersion=result$state$phi,dispersion.estimated=result$state$estimated,
    init.lambda=lambda,always.include=data$included,preprocessing=data$preprocessing,
    p=data$p_selectable,added_intercept=data$added_intercept,grouped=grouped,relaxed=relaxed,
    history=lapply(passes,`[[`,'history'),convergence=convergence),class='L0glm_fit')
  if(!is.null(args$validation)) {
    answer$validation <- lapply(passes,`[[`,'validation')
    answer$validation$deviance <- L0glm_validation_loss(args$validation,beta)
  }
  answer
}
