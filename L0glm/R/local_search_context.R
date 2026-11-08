# Native SIS scoring at a fitted model, rather than at an intercept-only null.
# Only the small active/nuisance information blocks are formed in R.
L0glm_local_context <- function(model, d, adjust.max=24L) {
  pattern<-L0glm_local_pattern(model);B<-as.matrix(pattern$matrix);active<-sort(unique(pattern$rows))
  rows <- sort(unique(c(active,d$mandatory)))
  if (length(rows)>adjust.max || length(rows)*d$m>64L)
    rows <- d$mandatory
  if (length(rows)>adjust.max) rows <- integer()
  N <- as.matrix(d$D[,rows,drop=FALSE]);q <- ncol(N);n <- d$n;m <- d$m
  eta <- L0glm_local_eta(B,d)
  ctx <- list(score=matrix(0,n,m),weight=matrix(0,n,m),nuisance=N,
    inverse=numeric(),correction=matrix(0,q,m),probability=matrix(0,0,0),
    means=matrix(0,0,0),basis=matrix(0,0,0),rows=seq_len(n),groups=integer(),
    e=numeric(),prefix=numeric(),cox=identical(d$kind,"cox"),
    multinomial=identical(d$kind,"multinomial"),one_sided=d$nonnegative,
    theta=model$theta)
  inverse <- L0glm_sis_inverse
  if (ctx$cox) {
    r <- L0glm_ext_cox_eta(as.numeric(eta),d$response,d$W[,1L],d$cox.cache)
    Z <- N[r$order,,drop=FALSE];HZ <- matrix(0,n,q)
    for(j in seq_len(q)) {
      a <- ifelse(r$S0>0,r$d*cumsum(r$e*Z[,j])[r$ends]/r$S0^2,0)
      HZ[,j] <- r$q*Z[,j]-r$e*rev(cumsum(rev(a)))[r$group]
    }
    G <- inverse(crossprod(Z,HZ))
    ctx$score <- matrix(r$score,n,1L);ctx$weight <- matrix(r$q,n,1L)
    ctx$nuisance <- HZ;ctx$inverse <- as.numeric(G)
    ctx$correction <- matrix(G%*%crossprod(Z,r$score),q,1L)
    ctx$rows <- as.integer(r$order);ctx$groups <- as.integer(r$group)
    ctx$e <- r$e;ctx$prefix <- c(0,cumsum(ifelse(r$S0>0,r$d/r$S0^2,0)))
    return(list(context=ctx,adjusted=rows,inverse=list(G)))
  }
  if (ctx$multinomial) {
    P <- L0glm_ext_softmax(eta)$p
    Q <- L0glm_ext_contrasts(m);r <- m-1L
    M <- P%*%Q;mass <- d$W[,1L]*rowSums(d$Y)
    score <- ((d$Y-rowSums(d$Y)*P)*d$W[,1L])%*%Q
    G <- matrix(0,q*r,q*r)
    for(a in seq_len(r))for(b in seq_len(a)) {
      h <- mass*(as.numeric(P%*%(Q[,a]*Q[,b]))-M[,a]*M[,b])
      aa <- seq_len(q)+(a-1L)*q;bb <- seq_len(q)+(b-1L)*q
      G[aa,bb] <- crossprod(N,h*N);G[bb,aa] <- t(G[aa,bb,drop=FALSE])
    }
    G <- inverse(G)
    ctx$score <- score;ctx$weight <- matrix(mass,n,1L);ctx$inverse <- as.numeric(G)
    ctx$correction <- matrix(G%*%as.numeric(crossprod(N,score)),q,r)
    ctx$probability <- P;ctx$means <- M;ctx$basis <- Q;ctx$one_sided <- FALSE
    return(list(context=ctx,adjusted=rows,inverse=list(G)))
  }
  inv <- array(0,c(q,q,m));phi <- L0glm_local_phi(model,d)
  for(j in seq_len(m)) {
    f <- if(identical(d$kind,"nb")) model$families[[j]] else d$family
    mu <- as.numeric(f$linkinv(eta[,j]))
    positive <- f$family%in%c("poisson","quasipoisson","Gamma","inverse.gaussian") ||
      identical(d$kind,"nb") || !is.null(f$scalar)
    if(positive) mu <- pmax(mu,1e-8)
    derivative <- f$mu.eta(eta[,j]);variance <- f$variance(mu)
    s <- d$W[,j]*(d$Y[,j]-mu)*derivative/variance/phi[j]
    h <- d$W[,j]*derivative^2/variance/phi[j]
    if(is.function(f$score)) {
      z <- f$score(d$Y[,j],mu)
      s <- -d$W[,j]*z[,1L]*derivative;h <- d$W[,j]*z[,2L]*derivative^2
    }
    s[d$W[,j]==0] <- h[d$W[,j]==0] <- 0
    s[!is.finite(s)] <- 0;h[!is.finite(h)|h<0] <- 0
    ctx$score[,j] <- s;ctx$weight[,j] <- h
    G <- inverse(crossprod(N,h*N));inv[,,j] <- G
    ctx$correction[,j] <- G%*%crossprod(N,s)
  }
  ctx$inverse <- as.numeric(inv)
  list(context=ctx,adjusted=rows,inverse=lapply(seq_len(m),function(j) matrix(inv[,,j],q,q)))
}

L0glm_local_scan <- function(model,d,control,exclude=integer()) {
  clock <- proc.time()[[3L]]
  context <- L0glm_local_context(model,d,control$adjust.max)
  # Scan the original matrix without copying it when nearly all columns qualify.
  pool <- d$pool-d$added;pool <- pool[pool>0L]
  if(!length(pool)) return(list(indices=integer(),stat=numeric(),context=context,seconds=0))
  if(length(pool)>.8*ncol(d$X)) {
    X <- d$X;map <- seq_len(ncol(X))+d$added
    count <- min(ncol(X),control$top+length(exclude)+ncol(X)-length(pool))
  } else {X <- d$X[,pool,drop=FALSE];map <- pool+d$added
    count <- min(length(map),control$top+length(exclude))}
  X <- L0glm_sis_design(X)
  result <- sis_score_cpp(X,context$context,as.integer(max(1L,count)),FALSE,
    as.integer(control$nthreads))
  chosen <- map[result$indices]
  keep <- which(d$in.pool[chosen] & !chosen%in%exclude)
  keep <- head(keep,control$top)
  list(indices=chosen[keep],stat=result$top_statistics[keep],context=context,
    seconds=proc.time()[[3L]]-clock)
}

L0glm_local_drop_rank <- function(model,d,context) {
  pattern<-L0glm_local_pattern(model);B<-as.matrix(pattern$matrix)
  select <- setdiff(sort(unique(pattern$rows)),d$mandatory)
  if(!length(select)) return(integer())
  value <- numeric(length(select))
  if(!identical(d$kind,"multinomial")) for(j in seq_len(d$m)) {
    X <- d$D[,select,drop=FALSE]
    if(context$context$cox) X <- X[context$context$rows,,drop=FALSE]
    h <- as.numeric(Matrix::colSums(X*(context$context$weight[,j]*X)))
    i <- match(select,context$adjusted);valid <- which(!is.na(i))
    if(length(valid)) {
      v <- diag(context$inverse[[j]])[i[valid]]
      h[valid] <- ifelse(v>0,1/v,h[valid])
    }
    value <- value+B[select,j]^2*pmax(h,0)
  } else value <- rowSums(B[select,,drop=FALSE]^2)
  select[order(value,select)]
}
