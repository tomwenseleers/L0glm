#' Nonnegative compulsory terms and smooth backgrounds
#' @name compulsory_terms
#' @md
#' @description
#' Supply unpenalized columns through `always.include`, for example a low-df
#' B-spline basis for a smooth chromatographic background. With `nonnegative=TRUE`,
#' `L0glm.fit()` and `L0nnpois()` refine these coefficients in a small constrained
#' block, separately from sparse predictor selection.
#'
#' @details
#' Write the linear predictor as
#' \deqn{\eta_j=X\beta_j+B\alpha_j+o_j,\qquad \alpha_j\geq0,}
#' where \eqn{B} contains the compulsory columns. For an identity link,
#' \eqn{B\alpha_j} is an additive background in response units. For other links
#' it is additive on the link scale. Binomial and quasibinomial intercepts retain
#' their unrestricted bounds; Cox models have no identifiable intercept.
#'
#' Conditional on \eqn{\beta}, the compulsory block minimizes the likelihood
#' loss subject to its bounds. Gaussian identity models use an exact weighted
#' nonnegative quadratic solve. Poisson identity models with nonnegative designs
#' and fixed signals use the observed likelihood curvature. Other scalar families
#' and links use safeguarded Fisher scoring:
#' \deqn{\alpha^{+}=\arg\min_{a\in\mathcal C}
#' (z-X\beta-o-Ba)^T W(z-X\beta-o-Ba).}
#' An active-set solver handles both unrestricted and nonnegative coefficients.
#' Backtracking keeps the family's mean domain and checks the actual likelihood
#' loss. Cox uses its partial-likelihood score and information. Negative-binomial
#' fits also update estimated shape parameters. Multinomial models retain their
#' unconstrained CLR formulation.
#'
#' The fitter initializes jointly, then alternates the compulsory solve with
#' sparse fitting conditional on its contribution. At most two outer cycles and
#' a final compulsory-block polish are used. Each sparse refinement uses two
#' grouped passes unless `group.maxit` was supplied explicitly. Updates are
#' accepted using the original-data likelihood plus the requested complexity
#' penalty; validation stopping instead uses the held-out deviance. Quasi families
#' use scaled deviance plus complexity. The original predictor count determines
#' the preset penalty throughout. Compulsory parameters retain their usual
#' contribution to information-criterion degrees of freedom and estimated scale.
#'
#' For positive-mean identity-link distributions (Poisson, Gamma, inverse Gaussian
#' and negative binomial) with nonnegative designs, automatic starts use a positive
#' signal backprojection above a low channel-wise background. The pilot accounts
#' for predictor scales, observation weights and known offsets, avoiding excessive
#' initial shrinkage of peaks relative to unpenalized background terms.
#' Gaussian models retain their ridge pilot, rescaled to the response units.
#' Supplied starts keep their meaning;
#' zero starts exclude selectable coefficients while compulsory terms remain
#' eligible. `maxit=1` retains a single ridge-pilot iteration.
#'
#' Both dense and sparse designs retain their storage type. Native Poisson and
#' Gaussian compulsory solves use OpenMP across channels; native Poisson sparse
#' fitting also retains its OpenMP backend. With `precluster=TRUE`, screening uses
#' summed Poisson channels and every refinement is scored on the original channels.
#'
#' The fitted object contains `nuisance`, with compulsory indices, candidate rows,
#' initialization, accepted/rejected scores, conditional KKT residuals and timings.
#'
#' @seealso [L0glm.fit()], [L0nnpois()], [preclustering]
NULL

# Conditional handling of nonnegative unpenalized terms.
# A private control entry prevents recursive orchestration of primitive passes.
L0glm_nuisance_coefficients <- function(value, Z, added, mandatory, lower, normalize, m, start=FALSE) {
  if(is.null(value)) return(NULL)
  p <- ncol(Z)
  if(start && is.null(dim(value)) && length(value)==1L && isTRUE(value==1))
    value <- L0glm_unit_start(value,Z,L0glm_preprocessing(Z,normalize,lower,mandatory),m)
  else if(length(value)==1L) value <- matrix(as.numeric(value),p,m) else {
    value <- as.matrix(value)
    if(added && nrow(value)==p-1L) value <- rbind(0,value)
    if(ncol(value)==1L) value <- matrix(as.numeric(value),nrow(value),m)
  }
  if(nrow(value)!=p || ncol(value)!=m || any(!is.finite(value)))
    stop("Starting coefficients and prior.mean must match the original coefficient dimensions.",call.=FALSE)
  if(start) value <- pmax(value,lower)
  value
}

L0glm_nuisance_conditional <- function(B,Y,signal,W,alpha,lower,families,threads=1L,cox=NULL) {
  m <- ncol(alpha);q <- nrow(alpha);bounded <- lower==0
  first <- families[[1L]]
  common <- all(vapply(families,function(f) identical(f$family,first$family)&&identical(f$link,first$link),logical(1)))
  if(is.null(cox) && common && first$family=="gaussian" && first$link=="identity")
    return(nuisance_gaussian_cpp(B,Y,signal,W,alpha,bounded,threads))
  if(is.null(cox) && common && first$family%in%c("poisson","quasipoisson") &&
     first$link=="identity" && all(bounded) && all(B>=0) && all(signal>=0))
    return(nuisance_poisson_cpp(B,Y,signal,W,alpha,60L,1e-8,threads))
  iterations <- integer(m);kkt <- rep(NA_real_,m)
  for(j in seq_len(m)) {
    a <- alpha[,j];eta <- as.numeric(signal[,j]+B%*%a);f <- families[[j]]
    loss <- if(is.null(cox)) function(e) L0glm_irls_deviance(Y[,j],e,W[,j],f) else
      function(e) -2*L0glm_ext_cox_eta(e,cox$response,W[,1L],cox$cache)$loglik
    old <- loss(eta)
    if(!is.finite(old)) next
    for(iter in seq_len(60L)) {
      if(is.null(cox)) {
        mu <- suppressWarnings(f$linkinv(eta))
        work <- L0glm_working_model(f)(Y[,j],eta,mu,W[,j],signal[,j])
        good <- work$good
        A <- B[good,,drop=FALSE];ww <- work$weights[good]
        H <- crossprod(A,A*ww)
        g <- as.numeric(crossprod(A,ww*(eta[good]-signal[good,j]-work$response[good])))
      } else {
        stats <- L0glm_ext_cox_eta(eta,cox$response,W[,1L],cox$cache)
        A <- B[stats$order,,drop=FALSE]
        H <- L0glm_ext_cox_information(A,stats)
        g <- -as.numeric(crossprod(A,stats$score))
      }
      if(any(!is.finite(H)) || any(!is.finite(g))) break
      proposal <- as.numeric(nuisance_quadratic_cpp(H,as.numeric(H%*%a-g),bounded))
      pred <- as.numeric(signal[,j]+B%*%proposal)
      step <- L0glm_step(proposal,a,function(value,e) loss(e),pred,eta,old,maxit=35L,rtol=1e-12)
      if(!step$accepted) break
      change <- sqrt(sum((step$beta-a)^2))/max(1,sqrt(sum(a^2)))
      a <- step$beta;eta <- step$eta;old <- step$value;iterations[j] <- iter
      kkt[j] <- max(abs(ifelse(bounded & a==0,pmin(g,0),g)))/max(1,max(abs(H%*%a)))
      if(change<1e-10) break
    }
    if(is.null(cox)) {
      edge <- L0glm_boundary(B,Y[,j],a,eta,W[,j],f,lower,signal[,j])
      a <- edge$beta
    }
    if(is.null(cox)) {
      eta <- as.numeric(signal[,j]+B%*%a)
      mu <- suppressWarnings(f$linkinv(eta))
      work <- L0glm_working_model(f)(Y[,j],eta,mu,W[,j],signal[,j])
      A <- B[work$good,,drop=FALSE];ww <- work$weights[work$good]
      H <- crossprod(A,A*ww)
      g <- as.numeric(crossprod(A,ww*(eta[work$good]-signal[work$good,j]-work$response[work$good])))
    } else {
      stats <- L0glm_ext_cox_eta(as.numeric(signal[,j]+B%*%a),cox$response,W[,1L],cox$cache)
      A <- B[stats$order,,drop=FALSE];H <- L0glm_ext_cox_information(A,stats)
      g <- -as.numeric(crossprod(A,stats$score))
    }
    if(all(is.finite(g)) && all(is.finite(H)))
      kkt[j] <- max(abs(ifelse(bounded & a==0,pmin(g,0),g)))/max(1,max(abs(H%*%a)))
    alpha[,j] <- a
  }
  list(coefficients=alpha,background=B%*%alpha,iterations=iterations,kkt=kkt)
}

L0glm_nuisance_fit <- function(X,y,options,fit.fun,supplied=character()) {
  begun <- proc.time()[[3L]];f <- options$family;kind <- L0glm_ext_kind(f)
  if(is.data.frame(X)) X <- as.matrix(X)
  if(is.object(X) && !inherits(X,"Matrix")) class(X) <- NULL
  added <- isTRUE(options$intercept);Z <- if(added) cbind(Intercept=1,X) else X
  included <- L0glm_mandatory(Z,if(added) c(1L,options$always.include+1L) else options$always.include)
  p <- ncol(Z);m <- if(identical(kind,"cox")) 1L else NCOL(y);n <- nrow(Z)
  lower <- rep(0,p)
  if(f$family%in%c("binomial","quasibinomial")) lower[L0glm_intercepts(Z)] <- -.Machine$double.xmax
  options$control.solver$.nuisance <- FALSE
  primitive <- function(design,response,opt) do.call(fit.fun,c(list(X=design,y=response),
    opt[intersect(names(opt),names(formals(fit.fun)))]))
  if(!length(included) || !any(lower[included]==0))
    return(primitive(X,y,options))
  selectable <- setdiff(seq_len(p),included);B <- as.matrix(Z[,included,drop=FALSE]);q <- ncol(B)
  Y <- if(identical(kind,"cox")) y else as.matrix(y)
  storage.mode(Y) <- "double"
  O <- L0glm_ext_matrix(options$offset,n,m,0,"offset")
  S <- L0glm_nuisance_coefficients(options$start,Z,added,included,lower,options$normalize,m,TRUE)
  P <- L0glm_nuisance_coefficients(options$prior.mean,Z,added,included,lower,options$normalize,m)
  allowed <- if(is.null(S)) matrix(TRUE,p,m) else S!=0
  allowed[included,] <- TRUE
  original_w <- L0glm_ext_matrix(options$weights,n,m,1,"weights")
  automatic_signal <- FALSE
  screened <- !is.null(options$screen) &&
    L0glm_screen_method(options$screen,X,options$start,added,options$always.include)!="none"
  if(is.null(S) && L0glm_signal_family(f) && !any(lower[included]<0) && !screened && !identical(kind,"cox")) {
    S <- L0glm_signal_start(Z,Y,O,original_w,included,P)
    automatic_signal <- !is.null(S)
  }
  n_ic <- if(identical(kind,"cox")) sum(as.matrix(Y)[,2L]*original_w[,1L]) else
    if(identical(kind,"nb")) sum(rowSums(original_w)>0) else n
  lambda <- L0glm_resolve_lambda(options$lambda,max(1L,length(selectable)),max(2,n_ic))
  if(!"relaxed"%in%supplied && identical(kind,"cox")) options$relaxed <- FALSE
  if(!is.null(options$screen) && !isTRUE(options$precluster) &&
      L0glm_screen_method(options$screen,X,options$start,added,options$always.include)!="none") {
    if(!"lambda"%in%supplied) lambda <- L0glm_resolve_lambda("gic",max(1L,length(selectable)),max(2,n_ic))
    if(!"relaxed"%in%supplied) options$relaxed <- TRUE
    if(!"parallelize"%in%supplied) options$parallelize <- FALSE
  }
  opt <- options;opt$lambda <- lambda;opt$intercept <- FALSE;opt$start <- S;opt$prior.mean <- P;opt$always.include <- included
  if(isTRUE(options$precluster) && !"group.maxit"%in%supplied) opt$group.maxit <- NULL
  # A partition-of-unity basis admits a stable constant initializer.
  zero_basis <- all(B==0)
  collapsed <- !identical(kind,"cox") && f$link=="identity" && !isTRUE(f$validmu(0)) && !zero_basis && q>1L && is.null(S) && all(lower[included]==0) && all(B>=0) &&
    max(abs(rowSums(B)-rowSums(B)[1L]))<=1e-12*max(1,abs(rowSums(B)[1L]))
  if(collapsed) {
    seed_design <- cbind(Z[,selectable,drop=FALSE],nuisance_constant=rep(1,n))
    opt$always.include <- ncol(seed_design);opt$validation <- NULL
    opt$prior.mean <- if(is.null(P)) NULL else rbind(P[selectable,,drop=FALSE],0)
  } else if(zero_basis) {
    seed_design <- Z[,selectable,drop=FALSE];opt$always.include <- NULL;opt$validation <- NULL
    opt$start <- if(is.null(S)) NULL else S[selectable,,drop=FALSE]
    opt$prior.mean <- if(is.null(P)) NULL else P[selectable,,drop=FALSE]
  } else seed_design <- Z
  init_clock <- proc.time()[[3L]]
  fit <- primitive(seed_design,Y,opt)
  init_seconds <- proc.time()[[3L]]-init_clock
  background_seconds <- peak_seconds <- 0
  W <- if(is.null(fit$prior.weights)) as.matrix(fit$weights) else as.matrix(fit$prior.weights)
  storage.mode(W) <- "double"
  beta <- matrix(0,p,m)
  seedmap <- if(collapsed) c(selectable,NA_integer_) else if(zero_basis) selectable else seq_len(p)
  put <- function(C) {
    out <- matrix(0,p,ncol(C))
    good <- which(!is.na(seedmap));out[seedmap[good],] <- as.matrix(C[good,,drop=FALSE])
    if(collapsed) out[included,] <- matrix(rep(as.numeric(C[nrow(C),])/rowSums(B)[1L],each=q),q,ncol(C))
    out
  }
  beta <- put(coef(fit));alpha <- beta[included,,drop=FALSE]
  candidate <- if(isTRUE(fit$preclustering$used)) fit$preclustering$candidates else
    if(!is.null(fit$screening$candidates)) fit$screening$candidates else which(Matrix::rowSums(coef(fit)!=0)>0)
  if(collapsed || zero_basis) candidate <- seedmap[candidate]
  candidate <- intersect(selectable,candidate[!is.na(candidate)])
  candidate <- candidate[rowSums(allowed[candidate,,drop=FALSE])>0]
  seed <- beta
  init <- fit$preclustering$initialization
  if(!is.null(init)) {
    C <- matrix(0,ncol(seed_design),m)
    C[init$rows,] <- if(is.null(dim(init$values))) matrix(init$values,length(init$rows),m) else init$values
    seed <- put(C)
  } else if(isTRUE(options$grouped)) seed[selectable,] <- rowMeans(beta[selectable,,drop=FALSE])
  seed[!allowed] <- 0
  # Save all diagnostics before the small peak-only fits replace their template.
  preclustering <- fit$preclustering;screening <- fit$screening
  if(collapsed && isTRUE(preclustering$used)) {
    remap <- function(rows) sort(unique(c(seedmap[rows[rows<length(seedmap)]],
      if(length(seedmap)%in%rows) included)))
    preclustering$candidates <- remap(preclustering$candidates)
    if(!is.null(preclustering$screen.coefficients))
      preclustering$screen.coefficients <- Matrix::Matrix(put(preclustering$screen.coefficients),sparse=TRUE)
    preclustering$initialization <- list(method="nuisance_constant_seed",rows=which(rowSums(seed!=0)>0),
      values=seed[rowSums(seed!=0)>0,,drop=FALSE])
  }
  phi_spec <- if(is.null(kind)) L0glm_dispersion(options$dispersion,f,m) else rep(1,m)
  estimated_phi <- if(is.null(phi_spec)) rep(TRUE,m) else phi_spec<0
  families <- function(model) if(identical(kind,"nb")) model$families else rep(list(f),m)
  cox <- if(identical(kind,"cox")) list(response=fit$response,
    cache=L0glm_ext_cox_cache(fit$response,W[,1L])) else NULL
  eta <- function(coefs) as.matrix(Z%*%coefs)+O
  poisson_constant <- if(f$family=="poisson" && all(Y==round(Y)))
    -2*sum(W*stats::dpois(Y,Y,log=TRUE)) else NULL
  score <- function(coefs,model) {
    pred <- eta(coefs)
    df <- sum(coefs[selectable,,drop=FALSE]!=0)+q*m
    if(identical(kind,"cox")) return(-2*L0glm_ext_cox_eta(as.numeric(pred),cox$response,W[,1L],cox$cache)$loglik+lambda*df)
    mu <- matrix(f$linkinv(pred),n,m)
    if(identical(kind,"nb")) {
      ll <- sum(vapply(seq_len(m),function(j) {good<-W[,j]>0
        sum(W[good,j]*stats::dnbinom(Y[good,j],mu=mu[good,j],size=model$theta[j],log=TRUE))},numeric(1)))
      return(-2*ll+lambda*(df+sum(model$theta.estimated)))
    }
    if(f$family=="poisson" && f$link=="identity" && !is.null(poisson_constant))
      return(sum(nuisance_poisson_deviance_cpp(Y,mu,W,1L))+poisson_constant+lambda*df)
    phi <- model$dispersion;phi[estimated_phi] <- -abs(phi[estimated_phi])
    ll <- L0glm_loglik_mean(Y,mu,W,f,phi)
    if(is.finite(ll)) return(-2*ll+lambda*(df+sum(estimated_phi)))
    # Quasi families use their penalized deviance rather than likelihood IC.
    sum(vapply(seq_len(m),function(j) L0glm_irls_deviance(Y[,j],pred[,j],W[,j],f)/
      max(.Machine$double.eps,model$dispersion[j]),numeric(1)))+lambda*df
  }
  validation <- options$validation
  merit <- function(coefs,model) if(is.null(validation)) score(coefs,model) else
    L0glm_validation_loss(validation,coefs,state=list(theta=model$theta))
  value <- merit(beta,fit)
  trace <- list(data.frame(cycle=0L,step="initialization",score=value,accepted=TRUE))
  threads <- if(isTRUE(options$parallelize)) max(1L,as.integer(options$ncores)) else 1L
  if(!is.finite(threads)) threads <- 1L
  details <- NULL;cycles <- 0L
  accept <- function(coefs,model,cycle,step) {
    next_value <- merit(coefs,model)
    keep <- is.finite(next_value) && (!is.finite(value)||next_value<=value+1e-10*max(1,abs(value)))
    trace[[length(trace)+1L]] <<- data.frame(cycle,step,score=next_value,accepted=keep)
    if(keep) value <<- next_value
    keep
  }
  polish <- function(cycle) {
    signal <- as.matrix(Z[,selectable,drop=FALSE]%*%beta[selectable,,drop=FALSE])+O
    clock <- proc.time()[[3L]]
    model <- fit
    update <- L0glm_nuisance_conditional(B,Y,signal,W,alpha,lower[included],families(model),threads,cox)
    trial <- beta;trial[included,] <- update$coefficients
    if(identical(kind,"nb") && any(model$theta.estimated)) for(shape_step in seq_len(6L)) {
      state <- list(theta=model$theta,estimated=model$theta.estimated)
      data <- list(m=m,Z=matrix(0,n,0),offset=eta(trial),family=f,y=Y,weights=W,
        control=model$control.family)
      for(inner in seq_len(10L)) {
        next_state <- L0glm_ext_nb_update(data,matrix(0,0,m),state)
        change <- max(abs(log(next_state$theta/state$theta)))
        state <- next_state
        if(change<model$control.family$theta.tol) break
      }
      change <- max(abs(log(state$theta/model$theta)));model$theta <- state$theta
      model$families <- lapply(model$theta,function(theta) {
        family <- unserialize(serialize(f,NULL));family$putTheta(log(theta));family
      })
      if(m==1L) model$family <- model$families[[1L]]
      if(change<model$control.family$theta.tol) break
      update <- L0glm_nuisance_conditional(B,Y,signal,W,update$coefficients,lower[included],families(model),threads,cox)
      trial[included,] <- update$coefficients
    }
    background_seconds <<- background_seconds+proc.time()[[3L]]-clock
    if(accept(trial,model,cycle,"conditional nuisance")) {
      beta <<- trial;alpha <<- update$coefficients;details <<- update;fit <<- model
    }
  }
  for(cycle in 1:2) {
    cycles <- cycle;polish(cycle)
    if(!length(candidate)) break
    pk <- options;pk$lambda <- lambda;pk$precluster <- FALSE;pk$intercept <- FALSE;pk$always.include <- NULL;pk$screen <- "none"
    pk$weights <- W;pk$offset <- O+B%*%alpha
    pk$prior.mean <- if(is.null(P)) NULL else P[candidate,,drop=FALSE]
    pk$start <- if(cycle==1L) seed[candidate,,drop=FALSE] else if(isTRUE(options$grouped))
      matrix(rowMeans(beta[candidate,,drop=FALSE]),length(candidate),m) else beta[candidate,,drop=FALSE]
    pk$start[!allowed[candidate,,drop=FALSE]] <- 0
    pk$group.maxit <- if("group.maxit"%in%supplied) options$group.maxit else 2L
    pk$control.solver$.nuisance.df <- q
    if(any(estimated_phi)) pk$dispersion <- ifelse(estimated_phi,-abs(fit$dispersion),fit$dispersion)
    if(identical(kind,"nb")) {pk$dispersion <- 1;pk$control.solver$.theta.start <- fit$theta}
    if(!is.null(validation)) {
      v <- L0glm_subset_validation(validation,candidate)
      v$offset <- validation$offset+as.matrix(validation$X[,included,drop=FALSE]%*%alpha)
      v$baseline <- v$refit_baseline <- beta[candidate,,drop=FALSE]
      v$baseline_dispersion <- fit$dispersion;v$baseline_state <- list(theta=fit$theta,estimated=fit$theta.estimated)
      pk$validation <- v
    }
    clock <- proc.time()[[3L]]
    proposal <- primitive(Z[,candidate,drop=FALSE],Y,pk)
    peak_seconds <- peak_seconds+proc.time()[[3L]]-clock
    trial <- beta;trial[selectable,] <- 0;trial[candidate,] <- as.matrix(coef(proposal))
    if(!accept(trial,proposal,cycle,"peak refinement")) break
    beta <- trial;fit <- proposal;details <- NULL
  }
  if(length(candidate)) polish(cycles)
  # Reconstruct the public original-dimension model and its working quantities.
  fit$X <- Z;fit$coefficients <- Matrix::Matrix(beta,sparse=TRUE)
  dimnames(fit$coefficients) <- list(colnames(Z),if(identical(kind,"cox")) "linear_predictor" else colnames(Y))
  fit$always.include <- included;fit$p <- length(selectable);fit$added_intercept <- added
  fit$offset <- Matrix::Matrix(O,sparse=TRUE);fit$init.lambda <- lambda
  fit$grouped <- options$grouped;fit$relaxed <- options$relaxed
  fit$preprocessing <- L0glm_preprocessing(Z,options$normalize,lower,included)
  fit$preclustering <- preclustering;fit$screening <- screening
  if(is.null(kind)) {
    fit$prior.weights <- Matrix::Matrix(W,sparse=TRUE)
    fit$y <- Y;pred <- eta(beta);mu <- matrix(f$linkinv(pred),n,m)
    wrk <- zr <- matrix(0,n,m)
    for(j in seq_len(m)) {
      if(estimated_phi[j]) {
        target <- L0glm_pearson_scale(Y[,j],mu[,j],W[,j],f,sum(beta[selectable,j]!=0)+q)
        if(is.finite(target) && target>0) fit$dispersion[j] <- target
      }
      work <- L0glm_working_model(f)(Y[,j],pred[,j],mu[,j],W[,j],O[,j])
      wrk[,j] <- work$weights;zr[,j] <- work$response
    }
    fit$weights <- if(m==1L) as.numeric(wrk) else Matrix::Matrix(wrk,sparse=TRUE)
    fit$z <- Matrix::Matrix(zr,sparse=TRUE)
  } else {
    fit$weights <- W
    if(identical(kind,"cox")) {fit$cox.cache <- cox$cache;fit$baseline <- L0glm_ext_cox_baseline(fit)}
    fit$metrics <- L0glm.metrics(fit,mu_floor=0);fit$objective <- fit$metrics$objective
  }
  penalties <- matrix(Inf,p,m);scale <- matrix(rep(fit$dispersion,each=p),p,m)
  penalties[beta!=0] <- lambda*scale[beta!=0]/4/pmax(beta[beta!=0]^2,.Machine$double.xmin)
  penalties[included,] <- 0;fit$lambdas <- penalties
  fit$selected <- which(rowSums(beta!=0)>0)
  if(!is.null(validation)) fit$validation$deviance <- merit(beta,fit)
  if(isTRUE(fit$preclustering$used)) {
    fit$preclustering$seconds["refit"] <- fit$preclustering$seconds["refit"]+peak_seconds
    fit$preclustering$seconds["nuisance"] <- background_seconds
    fit$preclustering$seconds["total"] <- proc.time()[[3L]]-begun
  }
  fit$nuisance <- list(terms=included,cycles=cycles,initializer=if(automatic_signal) "signal_backprojection" else if(collapsed) "constant_in_spline_space" else "joint",
    candidates=candidate,history=do.call(rbind,trace),kkt=if(is.null(details)) NULL else details$kkt,
    criterion=if(!is.null(validation)) "validation deviance" else "likelihood + lambda * df",
    seconds=c(initialization=init_seconds,conditional=background_seconds,refit=peak_seconds,
      total=proc.time()[[3L]]-begun))
  fit
}
