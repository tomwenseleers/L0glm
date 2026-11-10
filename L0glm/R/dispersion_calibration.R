#' Cross-channel calibration of unknown dispersion
#'
#' @name dispersion_calibration
#' @md
#' @description
#' Grouped multichannel models can use independent response channels to calibrate
#' unknown noise levels before a final adaptive-ridge fit. This supplements the
#' inexpensive ridge-pilot and bounded Pearson updates in [l0glmFit()].
#'
#' @details
#' `control.family=list(dispersion.calibration="auto")` is the default. It checks
#' models with more selectable columns than positive-weight observations and at
#' most ten times as many columns as observations. Calibration requires at least
#' two grouped response channels with unknown noise parameters, automatic starts
#' and zero ridge centers. A fourfold increase relative to the permissive initial
#' scale, or an unsupported requested relaxed refit, triggers a residual-score
#' check for remaining shared signal. `FALSE` disables calibration. `TRUE`
#' bypasses the width and trigger checks, while retaining the reliability checks.
#' Supplied positive or negative dispersions, starts and ridge centers are
#' respected. Screening, preclustering, validation and conditional nuisance fits
#' retain their existing estimators. Fixed-scale families require no calibration.
#'
#' Channels are reproducibly split into two balanced sets, preserving the caller's
#' RNG. A permissive, norm-based grouped pilot in one set ranks predictors for
#' unpenalized refitting in the other set. Each calibration model retains at most
#' 70% of its positive-weight observations in coefficient dimensions, including
#' mandatory terms. It must retain at least eight and 20% residual degrees of
#' freedom. The original family, link, weights, offsets and coefficient constraints
#' are used. The family variance function gives
#' \deqn{\widehat\phi_c=\frac{\sum_{i:w_{ic}>0}w_{ic}
#' (y_{ic}-\widehat\mu_{ic})^2/V(\widehat\mu_{ic})}{N_c-r_c}.}
#' Here \eqn{r_c} is the calibration model's coefficient dimension. For Gaussian
#' identity models with independent channel errors, selection on other channels
#' avoids selection on the target channel's noise. If the selected span contains
#' its mean, this is the usual unbiased residual-variance estimate. Omitted signal
#' biases the estimate upward. For other GLMs the same argument applies to the
#' local working model; the Pearson calculation is an approximation.
#'
#' Estimated negative-binomial shape uses a bounded NB2 moment equation,
#' \deqn{\sum_i\frac{w_{ic}(y_{ic}-\widehat\mu_{ic})^2}
#' {\widehat\mu_{ic}+\widehat\mu_{ic}^2/\theta_c}=\nu_c(\theta_c).}
#' Unit case weights give \eqn{\nu_c=N_c-r_c}. Unequal NB case weights use
#' \eqn{\nu_c=\sum_i w_{ic}(1-h_{ic})}, where \eqn{h_{ic}} is the diagonal
#' of the weighted working-model projection. This accounts for the power-
#' likelihood interpretation of NB weights. Shape bounds follow `theta.min`
#' and `theta.max` in `control.family`.
#'
#' Rank-deficient or inaccurate calibration fits, invalid scales and substantial
#' residual dependence across channel sets retain the original fit. A repair
#' requires at least a twofold reduction in estimated variance (or increase in NB
#' shape) and an acceptable final refit. A full-response pilot supplies its start;
#' calibrated noise parameters remain fixed during its IRLS, relaxed and local-search passes,
#' while retaining their data-estimated parameter-count flags. Calibration refits
#' use default iteration budgets of 60 and 100, respectively, respecting a larger
#' supplied `maxit`. The repair is kept
#' only if its penalized likelihood, or quasi-deviance kernel, improves when both
#' models are compared at the common calibrated noise level.
#'
#' Single-channel fits retain the guarded ridge initialization, bounded updates
#' and relaxed-refit information checks. They do not sacrifice observations for
#' a noise-calibration holdout. This fallback also applies whenever independent
#' calibration is unavailable. [l0glmTune()] performs calibration within each
#' full-data candidate and each CV training fit using only that fit's responses.
#'
#' The reported likelihood and IC methods retain their conventions in
#' [model_criteria], including profiled Gaussian variance. The common-scale repair
#' score in `dispersionCalibration$objective` is a separate fitting diagnostic.
#' The channel-split adaptation implements the independent-screen/refit principle
#' of Fan, Guo and Hao's
#' [refitted cross-validation](https://arxiv.org/abs/1004.5178).
#'
#' @return Applicable fits include `dispersionCalibration`, with `applied`, a
#'   fallback `reason`, or calibrated values, channel folds, calibration supports,
#'   residual degrees of freedom, timings and common-noise objective comparisons.
#'   `initialDispersion` and `initialTheta` record initialization values for the
#'   ordinary and NB drivers, respectively.
#' @seealso [l0glmFit()], [l0glmTune()], [model_criteria]
NULL

# Independent channel screening supplies unshrunk calibration models.
L0glm_calibration_mode <- function(mode) {
 if(!(isTRUE(mode)||isFALSE(mode)||identical(mode,"auto")))
   stop('control.family$dispersion.calibration must be TRUE, FALSE or "auto".',call.=FALSE)
 mode
}
L0glm_calibrated_fit <- function(fit,args) {
 mode <- L0glm_calibration_mode(if(is.null(args$control.family$dispersion.calibration))
   "auto" else args$control.family$dispersion.calibration)
 if(isFALSE(mode)||isFALSE(args$control.solver$.dispersion.calibration)||
    !is.null(fit$dispersionCalibration))return(fit)
 nb <- identical(fit$kind,"nb")
 estimated <- if(nb)fit$theta.estimated else fit$dispersion.estimated
 if(!length(estimated)||!all(estimated))return(fit)
 fallback <- function(reason,details=NULL) {
  fit$dispersionCalibration <- list(method="cross-channel",applied=FALSE,reason=reason,details=details)
  fit
 }
 m <- NCOL(fit$y);n <- nrow(fit$X);p <- ncol(fit$X)
 if(m<2L)return(fallback("single_channel"))
 if(!isTRUE(fit$grouped))return(fallback("independent_supports"))
 if(!is.null(args$dispersion)||!is.null(args$start)||!is.null(args$prior.mean))
   return(fallback("supplied_initialization"))
 if(isTRUE(args$maxit==1L)||isTRUE(fit$init.lambda==0))return(fit)
 if(inherits(fit$X,"L0glm_data")||!is.null(fit$screening)||!is.null(fit$preclustering)||
    !is.null(fit$nuisance)||!is.null(args$validation)||isTRUE(args$interactions>0L)||
    !is.null(args$control.solver$.theta.start)||isTRUE(args$control.solver$.nuisance.df>0))
   return(fallback("conditional_fit"))
 if(!is.matrix(fit$X)&&!inherits(fit$X,"Matrix"))return(fallback("design"))
 W <- L0glm_observation_input(if(nb)fit$weights else fit$prior.weights,n,m,1,"weights")
 N <- vapply(seq_len(m),function(j)sum(L0glm_observation_column(W,j)>0),1L)
 q <- p-length(fit$always.include)
 if(identical(mode,"auto") && (q<=min(N)||q>10*n))return(fit)
 old <- if(nb)fit$theta else fit$dispersion
 initial <- if(nb)fit$initialTheta else fit$initialDispersion
 if(length(initial)!=m||any(!is.finite(initial)|initial<=0))
   return(fallback("initial_scale"))
 inflation <- if(nb)median(initial/old)else median(old/initial)
 if(identical(mode,"auto") && !any(fit$relaxationSkipped) &&
    (!is.finite(inflation)||inflation<=4))return(fit)
 started <- proc.time()[[3L]]
 answer <- tryCatch(withCallingHandlers(
   L0glm_channel_calibration(fit,args,W,N,initial,nb,mode),
   warning=function(w)invokeRestart("muffleWarning")),error=function(e)e)
 if(inherits(answer,"error"))return(fallback("calibration_failed",conditionMessage(answer)))
 if(!is.null(answer$reason)) {
  result <- fallback(answer$reason,answer$details)
  result$dispersionCalibration$seconds <- proc.time()[[3L]]-started
  return(result)
 }
 result <- answer$fit
 result$dispersionCalibration <- c(list(method="cross-channel",applied=TRUE,
   parameter=if(nb)"theta"else"dispersion",baseline=old,initial=initial,
   calibrated=answer$scale,inflation=inflation,seconds=proc.time()[[3L]]-started),
   answer[c("fold","support","residualDf","objective","residualCorrelation")])
 result
}

L0glm_calibration_fold <- function(m) {
 had <- exists(".Random.seed",.GlobalEnv,inherits=FALSE)
 old <- if(had)get(".Random.seed",.GlobalEnv)else NULL
 on.exit(if(had)assign(".Random.seed",old,.GlobalEnv)else rm(".Random.seed",envir=.GlobalEnv))
 set.seed(1701L);sample(rep(1:2,length.out=m))
}

L0glm_calibration_converged <- function(fit) {
 z <- fit$convergence
 if(is.null(z)||!nrow(z))return(FALSE)
 if(!is.null(z$channel))z <- z[!duplicated(z$channel,fromLast=TRUE),,drop=FALSE]
 else z <- tail(z,1L)
 if(all(z$converged%in%TRUE))return(TRUE)
 if(identical(fit$kind,"nb"))return(FALSE)
 # A flat objective can reach floating-point precision before the coefficient
 # tolerance. Check the likelihood score without changing main-fit convergence.
 if(any(as.matrix(coef(fit))<=0) && isTRUE(fit$specification$nonnegative))return(FALSE)
 mu <- as.matrix(fitted(fit));family <- fit$family;X <- as.matrix(fit$X)
 W <- L0glm_observation_input(fit$prior.weights,nrow(X),ncol(mu),1,"weights")
 failed <- if(is.null(z$channel))seq_len(ncol(mu))else z$channel[!z$converged]
 for(j in failed) {
  w <- L0glm_observation_column(W,j);good <- w>0
  derivative <- family$mu.eta(family$linkfun(mu[,j]))
  variance <- family$variance(mu[,j]);A <- X[good,,drop=FALSE]
  G <- crossprod(A*sqrt(w[good]*derivative[good]^2/variance[good]))
  U <- tryCatch(chol(G),error=function(e)NULL)
  if(is.null(U))return(FALSE)
  score <- crossprod(A,w[good]*(fit$y[good,j]-mu[good,j])*derivative[good]/variance[good])
  decrement <- sum(backsolve(U,score,transpose=TRUE)^2)
  dev <- sum(family$dev.resids(fit$y[,j],mu[,j],w))
  if(!is.finite(decrement)||decrement>1e-6*max(1,dev))return(FALSE)
 }
 TRUE
}

L0glm_channel_calibration <- function(base,args,W,N,initial,nb,mode) {
 X <- base$X;Y <- as.matrix(base$y);m <- ncol(Y);n <- nrow(Y);p <- ncol(X)
 included <- base$always.include;selectable <- setdiff(seq_len(p),included)
 O <- L0glm_observation_input(base$offset,n,m,0,"offset")
 control <- args$control.solver
 if(is.null(control))control <- controlSolverGen(solver="eigen",method=1L)
 control$.dispersion.calibration <- FALSE
 control$.criterion.p <- base$p
 family <- base$family
 pilot_family <- if(nb)MASS::negative.binomial(min(1e7,10*median(initial)),family$link)else family
 threads <- if(isFALSE(args$parallelize))1L else
   if(is.null(args$ncores))L0glm_default_ncores(TRUE)else args$ncores
 common <- list(X=X,intercept=FALSE,always.include=included,normalize=base$preprocessing$mode,
   solver=if(is.null(args$solver))"eigen"else args$solver,
   method=if(is.null(args$method))1L else args$method,
   parallelize=threads>1L,ncores=threads,nthreads=args$nthreads,grouped=TRUE,
   maxit=if(is.null(args$maxit))25L else args$maxit,
   tol=if(is.null(args$tol))1e-10 else args$tol,
   thresh=if(is.null(args$thresh))1e-5 else args$thresh,
   nonnegative=isTRUE(args$nonnegative),relaxed=FALSE,local.search=FALSE,screen="none",
   control.solver=control,control.family=if(is.null(args$control.family))list()else args$control.family)
 # A joint row-norm pilot avoids cancellation of signed channel coefficients.
 pilot <- function(cols) {
  a <- common;a$y <- Y[,cols,drop=FALSE];a$family <- pilot_family
  a$weights <- L0glm_observation_subset(W,cols);a$offset <- L0glm_observation_subset(O,cols)
  a$dispersion <- if(nb)1 else initial[cols]*.1
  a$lambda <- base$init.lambda
  a$group <- seq_len(p);a$group.multiplier <- as.numeric(length(cols))
  attr(a$group.multiplier,"dimension_based") <- TRUE
  do.call(l0glmFit,a)
 }
 if(identical(mode,"auto")) {
  gate_offset <- if(any(base$relaxationSkipped))O else family$linkfun(as.matrix(fitted(base)))
  score <- sisGLM(X,Y,family=pilot_family,weights=W,offset=gate_offset,nsis=1L,
    return_stats=FALSE,nuisance=X[,included,drop=FALSE],
    nonnegative=isTRUE(args$nonnegative),nthreads=1L)$top_statistics[1L]
  threshold <- stats::qchisq(1-.05/max(1,length(selectable)),df=m)
  if(!is.finite(score)||score<=threshold)return(list(reason="no_shared_signal"))
 }
 fold <- L0glm_calibration_fold(m);scale <- numeric(m);df <- integer(m)
 residual <- matrix(NA_real_,n,m);sets <- vector("list",2L)
 for(k in 1:2) {
  train <- which(fold!=k);target <- which(fold==k);screen_fit <- pilot(train)
  B <- as.matrix(coef(screen_fit))*screen_fit$preprocessing$scale
  energy <- rowSums(B^2)
  rows <- selectable[energy[selectable]>0]
  cap <- max(0L,min(length(rows),floor(.7*min(N[target]))-length(included)))
  S <- sort(unique(c(included,head(rows[order(energy[rows],decreasing=TRUE)],cap))))
  if(any(N[target]-length(S)<max(8L,ceiling(.2*min(N[target])))))
    return(list(reason="insufficient_information"))
  weights <- L0glm_observation_subset(W,target)
  if(any(!L0glm_relaxed_standard(X[,S,drop=FALSE],weights,
      matrix(TRUE,length(S),length(target)),TRUE)))
    return(list(reason="rank_deficient"))
  a <- common;a$X <- X[,S,drop=FALSE];a$y <- Y[,target,drop=FALSE]
  a$weights <- weights;a$offset <- L0glm_observation_subset(O,target)
  a$always.include <- match(included,S);a$grouped <- FALSE
  a$family <- if(nb)MASS::negative.binomial(median(initial),family$link)else family
  a$lambda <- 0;a$dispersion <- 1;a$maxit <- max(60L,common$maxit);a$tol <- 1e-8
  cal <- do.call(l0glmFit,a)
  if(!L0glm_calibration_converged(cal))return(list(reason="calibration_convergence"))
  mu <- as.matrix(fitted(cal));df[target] <- N[target]-length(S)
  for(t in seq_along(target)) {
   j <- target[t];w <- L0glm_observation_column(W,j);good <- w>0
   scale[j] <- if(nb)L0glm_calibration_theta(Y[,j],mu[,t],w,X[,S,drop=FALSE],
      family,df[j],c(base$control.family$theta.min,base$control.family$theta.max))else
      L0glm_pearson_scale(Y[,j],mu[,t],w,family,length(S))
   variance <- if(nb)mu[,t]+mu[,t]^2/scale[j]else family$variance(mu[,t])
   residual[good,j] <- (Y[good,j]-mu[good,t])/sqrt(variance[good]) * if(nb)1 else sqrt(w[good])
  }
  sets[[k]] <- S
 }
 if(any(!is.finite(scale)|scale<=0))return(list(reason="invalid_scale"))
 check <- unique(round(seq(1,m,length.out=min(m,64L))))
 correlation <- suppressWarnings(stats::cor(residual[,check,drop=FALSE],use="pairwise.complete.obs"))
 opposite <- outer(fold[check],fold[check],`!=`)
 maximum <- if(any(opposite))max(abs(correlation[opposite]),na.rm=TRUE)else 0
 limit <- tanh(stats::qnorm(1-.01/(2*max(1,sum(opposite))))/sqrt(max(1,min(df)-3)))
 if(!is.finite(maximum)||maximum>limit)
   return(list(reason="dependent_channels",details=c(correlation=maximum,limit=limit)))
 old <- if(nb)base$theta else base$dispersion
 improvement <- if(nb)median(scale/old)else median(old/scale)
 if(improvement<2)return(list(reason="consistent_scale"))
 mature <- pilot(seq_len(m));start <- as.matrix(coef(mature));original <- as.matrix(coef(base))
 missing <- rowSums(start!=0)==0 & rowSums(original!=0)>0 & !any(base$relaxationSkipped)
 start[missing,] <- original[missing,,drop=FALSE]
 seeds <- L0glm_group_start(start,isTRUE(args$nonnegative))$start
 zero <- start==0 & matrix(rowSums(start!=0)>0,p,m)
 start[zero] <- .1*seeds[zero]
 a <- common;a$family <- family;a$y <- Y;a$weights <- W;a$offset <- O
 a$start <- start;a$lambda <- base$init.lambda
 a$group.maxit <- if(is.null(args$group.maxit))6L else args$group.maxit
 a$relaxed <- if(is.null(args$relaxed))TRUE else args$relaxed
 a$maxit <- max(100L,common$maxit)
 a$control.solver$.dispersion.fixed <- TRUE
 if(nb) {a$control.solver$.theta.start <- scale;a$dispersion <- NULL}
 else a$dispersion <- -scale
 if(!is.null(base[["group"]])) {
  a$group <- base[["group"]];a$group.multiplier <- base$group.multiplier
 }
 result <- do.call(l0glmFit,a)
 result$added_intercept <- base$added_intercept
 if(!L0glm_calibration_converged(result))return(list(reason="refit_convergence"))
 before <- L0glm_calibration_score(base,scale,nb)
 after <- L0glm_calibration_score(result,scale,nb)
 if(!is.finite(after)||after>before+1e-8*max(1,abs(before)))
   return(list(reason="objective",details=c(before=before,after=after)))
 list(fit=result,scale=scale,fold=fold,support=sets,residualDf=df,
   residualCorrelation=maximum,objective=c(before=before,after=after))
}

L0glm_calibration_theta <- function(y,mu,w,X,family,df,bounds=c(1e-6,1e7)) {
 good <- w>0;y <- y[good];mu <- mu[good];w <- w[good]
 if(any(!is.finite(mu)|mu<=0))return(NA_real_)
 derivative <- family$mu.eta(family$linkfun(mu))
 A <- as.matrix(X[good,,drop=FALSE]);weight_df <- length(unique(w))>1L
 equation <- function(a) {
  variance <- mu+mu^2/exp(a)
  nu <- df*mean(w)
  if(weight_df && ncol(A)) {
   Z <- A*sqrt(w*derivative^2/variance);U <- tryCatch(chol(crossprod(Z)),error=function(e)NULL)
   if(is.null(U))return(NA_real_)
   leverage <- rowSums(t(backsolve(U,t(Z),transpose=TRUE))^2)
   nu <- sum(w*(1-leverage))
  }
  sum(w*(y-mu)^2/variance)/nu-1
 }
 lower <- log(bounds[1L]);upper <- log(bounds[2L])
 lo <- equation(lower);hi <- equation(upper)
 if(!is.finite(lo)||!is.finite(hi))return(NA_real_)
 if(hi<=0)return(exp(upper))
 if(lo>=0)return(exp(lower))
 exp(stats::uniroot(equation,c(lower,upper),tol=1e-6)$root)
}

L0glm_calibration_score <- function(fit,scale,nb) {
 model <- fit
 if(nb) {model$theta <- scale;model$theta.estimated[] <- FALSE}
 else {model$dispersion <- scale;model$dispersion.estimated[] <- FALSE}
 metrics <- l0glmMetrics(model,mu_floor=1e-8)
 if(!is.null(fit[["group"]]) && !identical(fit$group.cost.type,"dimensions") &&
    is.finite(metrics$group_objective))
   return(metrics$group_objective)
 likelihood <- metrics$loglik
 if(!is.finite(likelihood)) {
  mu <- as.matrix(fitted(model));W <- model$prior.weights
  likelihood <- -.5*sum(vapply(seq_len(ncol(mu)),function(j)
    sum(model$family$dev.resids(model$y[,j],mu[,j],L0glm_observation_column(W,j)))/scale[j],1))
 }
 -2*likelihood+fit$init.lambda*metrics$df
}
