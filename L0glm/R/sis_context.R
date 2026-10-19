# Null fits and sufficient statistics for native efficient-score screening.
L0glm_sis_inverse <- function(G) {
  if (!nrow(G)) return(G)
  scale <- sqrt(pmax(diag(G),0));active <- which(scale>0);answer <- G*0
  if (!length(active)) return(answer)
  norm <- tcrossprod(scale[active])
  eig <- eigen(G[active,active,drop=FALSE]/norm,symmetric=TRUE)
  keep <- eig$values>max(abs(eig$values))*1e-12
  if (any(keep)) answer[active,active] <-
    tcrossprod(sweep(eig$vectors[,keep,drop=FALSE],2L,sqrt(eig$values[keep]),"/"))/norm
  answer
}

L0glm_sis_context <- function(y,family,n,weights,offset,nuisance,nonnegative=FALSE) {
  kind <- L0glm_ext_kind(family)
  response <- if (is.null(kind)) NULL else switch(kind,
    nb=L0glm_ext_nb_response(y,family,n),
    multinomial=L0glm_ext_multinomial_response(y,family,n),cox=L0glm_ext_cox_response(y,n))
  Y <- if (is.null(kind)) as.matrix(y) else response$y
  if (!identical(kind,"cox") && (nrow(Y)!=n || !ncol(Y) ||
      (!is.numeric(Y) && !is.logical(Y)) || any(!is.finite(Y))))
    stop("y must have finite numeric responses with one row per observation.",call.=FALSE)
  m <- if (identical(kind,"cox")) 1L else ncol(Y)
  joint <- identical(kind,"multinomial") || identical(kind,"cox")
  w <- L0glm_ext_matrix(weights,n,if (joint) 1L else m,1,"weights")
  off <- L0glm_ext_matrix(offset,n,m,0,"offset")
  if (any(!is.finite(w)|w<0) || any(colSums(w)==0) || any(!is.finite(off)))
    stop("Each response needs positive finite weights; offsets must be finite.",call.=FALSE)
  if (identical(kind,"multinomial") && any(colSums(Y*w[,1L])==0))
    stop("Each multinomial category needs a positive-weight observation.",call.=FALSE)
  N <- if (is.null(nuisance)) matrix(if (identical(kind,"cox")) numeric() else 1,n,
    if (identical(kind,"cox")) 0L else 1L) else as.matrix(nuisance)
  if (nrow(N)!=n || !is.numeric(N) || any(!is.finite(N)))
    stop("nuisance must be a finite numeric matrix with one row per observation.",call.=FALSE)
  if (!is.double(N)) storage.mode(N) <- "double"
  if (identical(kind,"cox") && ncol(N)) N <- N[,!vapply(seq_len(ncol(N)),
    function(j) all(N[,j]==N[1L,j]),logical(1)),drop=FALSE]
  q <- ncol(N)
  context <- list(score=matrix(0,n,m),weight=matrix(0,n,m),nuisance=N,
    inverse=numeric(),correction=matrix(0,q,m),probability=matrix(0,0,0),
    means=matrix(0,0,0),basis=matrix(0,0,0),rows=seq_len(n),groups=integer(),
    e=numeric(),prefix=numeric(),cox=identical(kind,"cox"),
    multinomial=identical(kind,"multinomial"),one_sided=nonnegative,theta=numeric())
  if (identical(kind,"cox")) {
    if (!any(response$status*w[,1L]>0)) stop("Cox screening needs an observed event.",call.=FALSE)
    eta <- off[,1L]
    if (q) {
      fit <- L0glm.fit(N,y,family=family,weights=w,offset=off,lambda=0,
        intercept=FALSE,normalize="auto",maxit=100,tol=1e-10)
      eta <- as.numeric(N%*%coef(fit))+eta
    }
    risk <- L0glm_ext_cox_eta(eta,response,w[,1L])
    Z <- N[risk$order,,drop=FALSE];HZ <- matrix(0,n,q)
    for (j in seq_len(q)) {
      a <- ifelse(risk$S0>0,risk$d*cumsum(risk$e*Z[,j])[risk$ends]/risk$S0^2,0)
      HZ[,j] <- risk$q*Z[,j]-risk$e*rev(cumsum(rev(a)))[risk$group]
    }
    G <- L0glm_sis_inverse(crossprod(Z,HZ))
    context$score <- matrix(risk$score,n,1);context$weight <- matrix(risk$q,n,1)
    context$nuisance <- HZ;context$inverse <- as.numeric(G)
    context$correction <- matrix(G%*%crossprod(Z,risk$score),q,1)
    context$rows <- as.integer(risk$order);context$groups <- as.integer(risk$group)
    context$e <- risk$e
    context$prefix <- c(0,cumsum(ifelse(risk$S0>0,risk$d/risk$S0^2,0)))
    return(context)
  }
  if (identical(kind,"multinomial")) {
    eta <- off
    if (q) {
      fit <- L0glm.fit(N,y,family=family,weights=w,offset=off,lambda=0,
        intercept=FALSE,normalize="auto",maxit=100,tol=1e-10)
      eta <- as.matrix(N%*%coef(fit))+off
    }
    P <- L0glm_ext_softmax(eta)$p;Q <- L0glm_ext_contrasts(m);r <- m-1L
    M <- P%*%Q;mass <- w[,1L]*rowSums(Y)
    score <- ((Y-rowSums(Y)*P)*w[,1L])%*%Q
    G <- matrix(0,q*r,q*r)
    for (a in seq_len(r)) for (b in seq_len(a)) {
      value <- mass*(as.numeric(P%*%(Q[,a]*Q[,b]))-M[,a]*M[,b])
      rows <- seq_len(q)+(a-1L)*q;cols <- seq_len(q)+(b-1L)*q
      G[rows,cols] <- crossprod(N,value*N);G[cols,rows] <- t(G[rows,cols,drop=FALSE])
    }
    G <- L0glm_sis_inverse(G)
    context$score <- score;context$weight <- matrix(mass,n,1)
    context$inverse <- as.numeric(G)
    context$correction <- matrix(G%*%as.numeric(crossprod(N,score)),q,r)
    context$probability <- P;context$means <- M;context$basis <- Q
    context$one_sided <- FALSE
    return(context)
  }
  # Independent channels need q-by-q blocks, not a full channel Hessian.
  inverse <- array(0,c(q,q,m))
  if (identical(kind,"nb") && q) {
    fit <- L0glm.fit(N,Y,family=family,weights=w,offset=off,lambda=0,
      intercept=FALSE,normalize="auto",maxit=150,tol=1e-9)
    context$theta <- fit$theta
  }
  for (j in seq_len(m)) {
    f <- family;good <- w[,j]>0
    if (identical(kind,"nb")) {
      if (q) {
        eta <- as.numeric(N%*%as.matrix(coef(fit))[,j])+off[,j]
        f <- fit$families[[j]]
      } else {
        eta <- off[,j];mu <- family$linkinv(eta)
        theta <- if (family$n.theta>0) exp(stats::optimize(function(a)
          -sum(w[good,j]*stats::dnbinom(Y[good,j],mu=mu[good],size=exp(a),log=TRUE)),
          log(c(1e-6,1e7)))$minimum) else family$getTheta(TRUE)
        f <- unserialize(serialize(family,NULL));f$putTheta(log(theta));context$theta[j] <- theta
      }
      mu <- f$linkinv(eta);phi <- 1
    } else {
      fit0 <- stats::glm.fit(N,Y[,j],weights=w[,j],offset=off[,j],family=f,
        intercept=any(vapply(seq_len(q),function(k) all(N[,k]==1),logical(1))),
        control=stats::glm.control(maxit=100,epsilon=1e-10))
      eta <- fit0$linear.predictors;mu <- fit0$fitted.values
      phi <- if (f$family%in%c("poisson","binomial") || grepl("negative binomial",f$family,ignore.case=TRUE)) 1 else
        sum(w[good,j]*(Y[good,j]-mu[good])^2/f$variance(mu[good]))/max(1,sum(good)-fit0$rank)
    }
    derivative <- f$mu.eta(eta)[good];variance <- f$variance(mu[good])
    s <- w[good,j]*(Y[good,j]-mu[good])*derivative/variance/max(phi,.Machine$double.eps)
    h <- w[good,j]*derivative^2/variance/max(phi,.Machine$double.eps)
    boundary <- variance==0 & Y[good,j]==mu[good]
    s[boundary] <- h[boundary] <- 0
    if (any(!is.finite(s)) || any(!is.finite(h)|h<0))
      stop("Undefined null score/information; check family, offsets and nuisance fit.",call.=FALSE)
    context$score[good,j] <- s;context$weight[good,j] <- h
    G <- L0glm_sis_inverse(crossprod(N,context$weight[,j]*N))
    inverse[,,j] <- G;context$correction[,j] <- G%*%crossprod(N,context$score[,j])
  }
  context$inverse <- as.numeric(inverse)
  context
}
