#' Select groups of predictor variables
#'
#' @name predictor_groups
#' @md
#' @description
#' Select non-overlapping groups of predictors using adaptive ridge and a group
#' L0 penalty. Groups can represent factor contrasts, spline coefficients,
#' related measurements or other terms that should enter a model together.
#'
#' @param group Optional vector with one label per column of the supplied `X`,
#'   excluding an automatically added intercept. Factors, integers and character
#'   labels are accepted; ordering and consecutive numbering are unnecessary.
#'   Label `0` or `"0"` identifies unpenalized terms. Intercepts and
#'   `always.include` terms remain unpenalized. `NULL` uses individual predictors.
#' @param group.multiplier Optional positive penalty multipliers, one per
#'   penalized predictor group in first-appearance order, or a named vector
#'   matching the nonzero group labels. A scalar applies to every group.
#'   By default, named IC presets charge each group's coefficient dimension;
#'   a numeric `lambda` charges once per group. Returned `group.multiplier`
#'   values reproduce the preset when reusing its numeric lambda. Tuning uses
#'   dimension-based multipliers throughout its path and training folds.
#'
#' @details
#' For working coefficients \eqn{b_{jc}=s_j\beta_{jc}} and predictor group
#' \eqn{G_g}, independent-channel selection (`grouped=FALSE`) uses
#' \deqn{r_{gc}^2=\sum_{j\in G_g}b_{jc}^2.}
#' Joint selection (`grouped=TRUE`) uses
#' \deqn{r_g^2=\sum_c\sum_{j\in G_g}b_{jc}^2.}
#' The norm prevents opposing signs from cancelling. Each penalized coefficient
#' in a block receives ridge diagonal
#' \deqn{D_{jc}=\frac{\lambda a_g}{4r_g^2},}
#' with the appropriate channel-specific norm for independent selection. Scalar
#' GLM working weights are divided by the channel's dispersion, so this expression
#' is on the likelihood scale. Equivalently, an unscaled weighted least-squares
#' system multiplies the diagonal by dispersion. Adaptive weights use coefficients
#' relative to zero; `prior.mean` remains the fixed ridge target.
#'
#' The approximate selection target is
#' \deqn{-2\ell(B,\psi)+\lambda\sum_g a_g\,1(r_g>0),}
#' with separate blocks per response when `grouped=FALSE`. A numeric lambda uses
#' \eqn{a_g=1} by default. IC presets use the number of coefficients in the group,
#' multiplied by the number of responses for joint selection. Joint multinomial
#' groups use \eqn{(C-1)|G_g|} dimensions for \eqn{C} CLR categories. Multinomial
#' predictor groups require `grouped=TRUE`. Noise parameters and compulsory terms
#' contribute to reported information criteria; their complexity charges are
#' constant across group-support comparisons.
#' Reported GIC, AIC and BIC count free nonzero coefficient entries, with the
#' usual mandatory-term, noise-parameter and CLR adjustments. Group size is
#' accounted for through these entries, without an additional size multiplier.
#' When a retained group contains exact zeros, its reported parameter count can
#' be smaller than the nominal group cost used in the adaptive-ridge update.
#'
#' With zero ridge targets, the inverse-norm update is the quadratic
#' majorization of a log-norm penalty. For an orthogonal block with curvature
#' \eqn{h}, its nonzero fixed point obeys
#' \eqn{r^2-\widehat r r+\lambda a_g/(4h)=0}; its existence boundary
#' \eqn{h\widehat r^2=\lambda a_g} matches the block-L0 selection threshold.
#' The relaxed refit removes shrinkage on selected blocks. Coefficients can be
#' zero inside a retained group, for example at a nonnegative boundary.
#'
#' Hard thresholding removes whole blocks using `thresh` and normalized norms.
#' Explicit starting-zero exclusions remain binding for individual entries;
#' intercepts and compulsory terms remain eligible. Structural zero columns do
#' not enter solves. The returned sparse coefficient matrix retains the original
#' predictor dimensions. `group_support_size`, `group_df` and `group_objective`
#' in [l0glmMetrics()] describe the group fit, alongside coefficient-based scores.
#'
#' Group models use the shared structured IRLS driver and the existing ridge
#' solvers for ordinary GLMs, negative binomial, joint multinomial and Cox models.
#' Predictor-group screening and outcome preclustering are unavailable; use
#' `screen="none"` to fit a large group design directly. Local refinement tests
#' additions, deletions and exchanges of complete predictor groups, accepting
#' only improvements in the group objective (or supplied validation deviance).
#' Named presets and tuning refine the reported parameter-count IC. Numeric
#' penalties refine the weighted group-count objective; `control.local$target`
#' can explicitly choose `"IC"` or `"L0"`.
#' Automatic refinement uses the selected coefficient-dimension limit described
#' in [local_search]. Cholesky is used for regular systems; a rank-deficient
#' unpenalized group refit can fall back to pivoted QR.
#'
#' @return Fits retain `group`, `group.multiplier`, eligible-entry masks and
#'   selected-block masks, together with the usual fit fields. Use the stored
#'   multipliers with a numeric lambda to reproduce an IC-preset or tuned fit.
#' @seealso [l0glmFit()], [l0nnPois()], [l0glmTune()], [l0glmMetrics()]
#' @examples
#' set.seed(41)
#' X <- matrix(rnorm(200 * 12), 200, 12)
#' group <- rep(seq_len(6), each=2)
#' group[11:12] <- 0
#' beta <- c(1.5, -1.5, 0, 0, 1, -1, rep(0, 6))
#' y <- 1 + as.numeric(X %*% beta) + rnorm(200)
#' fit <- l0glmFit(X, y, group=group, lambda="gic", intercept=TRUE,
#'   relaxed=TRUE, local.search=FALSE)
#' l0glmMetrics(fit)[c("gic", "group_support_size", "group_df", "group_objective")]
#' same <- l0glmFit(X, y, group=group, lambda=fit$init.lambda,
#'   group.multiplier=fit$group.multiplier, intercept=TRUE,
#'   relaxed=TRUE, local.search=FALSE)
#' max(abs(coef(fit) - coef(same)))
NULL

# All norm calculations use working units.
L0glm_group_labels <- function(group, p) {
  if (!is.atomic(group) || length(group) != p || anyNA(group))
    stop("group needs one nonmissing label per input predictor.", call.=FALSE)
  label <- as.character(group)
  if (any(!nzchar(label))) stop("group labels must be nonempty.", call.=FALSE)
  label
}

L0glm_group_blocks <- function(B, group, grouped) {
  if (grouped) matrix(rep(group, ncol(B)), nrow(B), ncol(B)) else
    matrix(rep(group, ncol(B)) + rep((seq_len(ncol(B))-1L)*max(c(0,group)), each=nrow(B)),
      nrow(B), ncol(B))
}

L0glm_group_penalty <- function(B, allowed, included, lambda, grouped, groups) {
  block <- L0glm_group_blocks(B, groups$id, grouped)
  eligible <- allowed; eligible[included,] <- FALSE
  index <- which(eligible)
  D <- matrix(Inf,nrow(B),ncol(B))
  if (length(index)) {
    norm2 <- rowsum(as.vector(B^2)[index], as.vector(block)[index], reorder=FALSE)
    key <- match(as.vector(block)[index], as.numeric(rownames(norm2)))
    cost <- rep(groups$multiplier[groups$id],ncol(B))[index]
    D[index] <- if(lambda==0)0 else lambda*cost/(4*norm2[key,1L])
  }
  D[included,] <- 0
  D
}

L0glm_group_threshold <- function(B, allowed, included, thresh, grouped, groups) {
  block <- L0glm_group_blocks(B,groups$id,grouped)
  eligible <- allowed; eligible[included,] <- FALSE
  index <- which(eligible)
  if(length(index)) {
    norm <- sqrt(rowsum(as.vector(B^2)[index],as.vector(block)[index],reorder=FALSE)[,1L])
    drop <- as.numeric(names(norm))[norm < max(norm)*thresh | norm < 1e-20]
    remove <- block %in% drop & eligible
    B[remove] <- 0; allowed[remove] <- FALSE
  }
  allowed[included,] <- TRUE
  list(B=B,allowed=allowed)
}

L0glm_group_active <- function(B, allowed, included, grouped, groups) {
  block <- L0glm_group_blocks(B,groups$id,grouped)
  active <- matrix(block %in% block[B!=0],nrow(B),ncol(B)) & allowed
  active[included,] <- TRUE
  active
}

L0glm_group_cost <- function(label,m,grouped,kind,dimension,multiplier=NULL,estimable=rep(TRUE,length(label))) {
  labels <- unique(label[label!='0'])
  cost <- if(dimension)pmax(1,as.numeric(table(factor(label[estimable],levels=labels)))) *
    if(grouped) m-as.integer(identical(kind,'multinomial')) else 1 else rep(1,length(labels))
  if(!is.null(multiplier)) {
    if(!is.numeric(multiplier) || any(!is.finite(multiplier)|multiplier<=0) ||
       !length(multiplier)%in%c(1L,length(labels)))
      stop('group.multiplier needs positive finite values, one per penalized group.',call.=FALSE)
    if(!is.null(names(multiplier))) {
      if(anyDuplicated(names(multiplier)) || !setequal(names(multiplier),labels))
        stop('Named group.multiplier values must match penalized group labels.',call.=FALSE)
      multiplier <- multiplier[labels]
    }
    cost <- rep_len(multiplier,length(labels))
  }
  stats::setNames(as.numeric(cost),labels)
}

L0glm_group_tune_cost <- function(X,y,family,args) {
  if(is.null(args[["group"]]))return(args)
  if(is.null(args$control.local$target))args$control.local$target <- 'IC'
  label <- L0glm_group_labels(args[["group"]],ncol(X))
  included <- union(L0glm_mandatory(X,args$always.include),which(label=='0'))
  label[included] <- '0';args$always.include <- included
  if(!is.null(args$group.multiplier))return(args)
  kind <- L0glm_ext_kind(family)
  m <- if(identical(kind,'multinomial'))family$nlp+1L else NCOL(y)
  grouped <- if(is.null(args$grouped))TRUE else args$grouped
  args$group.multiplier <- L0glm_group_cost(label,m,grouped,kind,TRUE)
  attr(args$group.multiplier,'dimension_based') <- TRUE
  args
}

L0glm_group_metrics <- function(object,scores) {
  B <- coef(object);mandatory <- object$always.include
  dimension <- sum(B!=0)+sum(B[mandatory,,drop=FALSE]==0)
  if(identical(object$kind,'multinomial')) {
    selectable <- setdiff(seq_len(nrow(B)),mandatory)
    rows <- sum(Matrix::rowSums(B[selectable,,drop=FALSE]!=0)>0)
    dimension <- if(object$grouped)(rows+length(mandatory))*(ncol(B)-1L) else
      sum(B[selectable,,drop=FALSE]!=0)-rows+length(mandatory)*(ncol(B)-1L)
  }
  id <- match(object[["group"]],unique(object[["group"]]))
  block <- L0glm_group_blocks(B,id,object$grouped)
  penalized <- as.matrix(B!=0)
  if(length(mandatory))penalized[mandatory,] <- FALSE
  units <- unique(block[penalized])
  labels <- unique(object[["group"]])[(units-1L)%%max(id)+1L]
  charge <- sum(object$group.multiplier[labels])
  scores$group_support_size <- length(units)
  scores$group_df <- dimension-length(mandatory)*(ncol(B)-as.integer(identical(object$kind,'multinomial')))
  scores$group_objective <- -2*scores$loglik+object$init.lambda*charge
  scores
}

L0glm_group_fit <- function(args) {
  mode <- L0glm_local_mode(args$local.search)
  if(is.null(args$control.local$target))args$control.local$target <-
    if(is.character(args$lambda) || isTRUE(attr(args$group.multiplier,'dimension_based')))'IC'else'L0'
  if(!isFALSE(mode))return(L0glm_local_fit(l0glmFit,args,mode,args$control.local))
  family <- args$family
  if(!inherits(family,'family'))family <- family()
  label <- L0glm_group_labels(args[["group"]],ncol(args$X))
  kind <- L0glm_ext_kind(family)
  if(is.null(kind))kind <- 'glm'
  validate_thresh(args$thresh)
  if(length(args$maxit)!=1L || !is.numeric(args$maxit) || !is.finite(args$maxit) ||
     args$maxit<1 || args$maxit!=floor(args$maxit) || length(args$tol)!=1L ||
     !is.numeric(args$tol) || !is.finite(args$tol) || args$tol<=0)
    stop('maxit must be a positive integer and tol a positive finite number.',call.=FALSE)
  for(name in c('grouped','relaxed','intercept','nonnegative'))
    if(!is.logical(args[[name]]) || length(args[[name]])!=1L || is.na(args[[name]]))
      stop(name,' must be TRUE or FALSE.',call.=FALSE)
  if(kind=='multinomial' && !args$grouped)
    stop('Multinomial predictor groups require grouped=TRUE to respect joint CLR dimensions.',call.=FALSE)
  if(kind=='cox')args$intercept <- FALSE
  if(args$precluster)stop('Predictor groups with outcome preclustering are not implemented.',call.=FALSE)
  if(!tolower(args$screen)%in%c('auto','none'))stop('Group screening is not implemented yet.',call.=FALSE)
  include <- union(L0glm_mandatory(args$X,args$always.include),which(label=='0'))
  if(L0glm_screen_method(args$screen,args$X,args$start,args$intercept,include)!='none')
    stop('Group screening is not implemented yet; use screen="none" to fit all groups.',call.=FALSE)
  if(!is.null(args$validation) && !inherits(args$validation,'L0glm_validation'))
    return(L0glm_validated_fit(args))
  if(args$intercept)include <- include+1L
  data <- L0glm_ext_prepare(args$X,args$y,family,args$weights,args$offset,args$intercept,
    include,args$nonnegative,args$normalize,args$start,args$prior.mean,args$control.family,kind=kind)
  data$p_selectable <- L0glm_criterion_size(data$p_selectable,args$control.solver)
  structural <- which(L0glm_moments(data$Z)$ss==0)
  data$allowed[structural,] <- FALSE;data$B[structural,] <- 0
  data$included <- setdiff(data$included,structural)
  if(args$intercept)label <- c('0',label)
  label[data$included] <- '0'
  multiplier <- L0glm_group_cost(label,data$m,args$grouped,kind,is.character(args$lambda),args$group.multiplier,
    !seq_along(label)%in%structural)
  id <- match(label,unique(label));cost <- numeric(length(unique(label)))
  cost[match(names(multiplier),unique(label))] <- multiplier
  data$groups <- list(id=id,multiplier=cost)
  if(kind=='cox' && data$control$cox=='pcg' &&
     (args$solver!='eigen' || any(data$lower > -1e100)))
    stop('Cox cox="pcg" requires solver="eigen" and unrestricted coefficients.',call.=FALSE)
  lambda <- L0glm_resolve_lambda(args$lambda,max(1,data$p_selectable),max(2,data$n_ic))
  options <- list(solver=args$solver,method=args$method,control.solver=args$control.solver,
    validation=args$validation,predictor.groups=TRUE)
  state <- if(kind=='nb')list(theta=data$response$theta,estimated=data$response$estimated)else list()
  B <- data$B
  if(kind=='glm') {
    initial <- L0glm_group_glm_initialize(data,args$dispersion,lambda)
    data <- initial$data; state <- initial$state
    if(!data$explicit)B <- initial$B
    if(!is.null(options$validation)) {
      options$validation$baseline_state$phi <- options$validation$baseline_dispersion
      options$validation$baseline_state$estimated <- state$estimated
      options$validation$baseline_state$stable <- state$stable
    }
  } else if(!data$explicit || kind=='multinomial' && all(B==0)) {
    D <- matrix(1e-3,data$p,data$m);D[!data$allowed] <- Inf;D[data$included,] <- 0
    B <- switch(kind,nb=L0glm_ext_nb_working(data,B,D,state,options),
      cox=L0glm_ext_cox_initial(data),multinomial=L0glm_ext_multinomial_working(data,B,D,state,options))
  }
  B[!data$allowed] <- 0
  if(kind%in%c('glm','nb') && !is.finite(suppressWarnings(L0glm_ext_loss(data,B,state))))
    B <- L0glm_group_feasible(data,B,state)
  selection <- L0glm_ext_pass(data,B,data$allowed,state,lambda,args$grouped,args$thresh,
    as.integer(args$maxit),args$tol,options)
  result <- selection;passes <- list(selection=selection)
  if(args$relaxed) {
    active <- L0glm_group_active(selection$B,data$allowed,data$included,args$grouped,data$groups)
    result <- L0glm_ext_pass(data,selection$B,active,selection$state,0,args$grouped,0,
      as.integer(args$maxit),args$tol,options,TRUE)
    passes$relaxed <- result
  }
  if(args$relaxed && args$nonnegative) {
    result <- L0glm_ext_boundary(data,result,options$validation)
    passes$relaxed <- result
  }
  fit <- L0glm_ext_result(data,result,passes,lambda,args$grouped,args$relaxed,options)
  fit$coefficients <- methods::as(methods::as(fit$coefficients,'generalMatrix'),'CsparseMatrix')
  fit[["group"]] <- label
  fit$group.multiplier <- multiplier
  fit$group.cost.type <- if(is.character(args$lambda) ||
    isTRUE(attr(args$group.multiplier,'dimension_based')))'dimensions'else'specified'
  attr(fit$group.multiplier,'dimension_based') <- fit$group.cost.type=='dimensions'
  fit$group.allowed <- data$allowed
  fit$group.active <- L0glm_group_active(result$B,data$allowed,data$included,args$grouped,data$groups)
  fit$metrics <- l0glmMetrics(fit,mu_floor=0)
  fit$objective <- fit$metrics$group_objective
  fit
}

L0glm_group_feasible <- function(data,B,state) {
  anchor <- matrix(0,data$p,data$m)
  if(length(data$included))for(j in seq_len(data$m)) {
    family <- data$family
    if(data$kind=='nb') {
      family <- MASS::negative.binomial(state$theta[j],data$family$link)
    }
    null <- suppressWarnings(stats::glm.fit(as.matrix(data$Z[,data$included,drop=FALSE]),data$y[,j],
      weights=L0glm_observation_column(data$weights,j),offset=L0glm_observation_column(data$offset,j),
      family=family,intercept=FALSE,control=stats::glm.control(maxit=100L)))$coefficients
    null[!is.finite(null)] <- 0
    anchor[data$included,j] <- pmax(data$lower[data$included],null)
  }
  for(halving in seq_len(60L)) {
    trial <- anchor+2^(-halving)*(B-anchor)
    if(is.finite(suppressWarnings(L0glm_ext_loss(data,trial,state))))return(trial)
  }
  stop('Cannot initialize finite group-model means; supply a feasible start or an intercept/offset.',call.=FALSE)
}

L0glm_group_ridge <- function(X,y,w,D,start,lower,prior,args) {
  options <- list(X=X,y=y,weights=w,lambdas=D,start=start,lower=lower,prior.mean=prior,
    solver=args$solver,method=args$method,control.solver=args$control.solver)
  failed <- FALSE
  value <- withCallingHandlers(do.call(lmRidge,options),warning=function(warning) {
    if(grepl('Eigen Cholesky failed',conditionMessage(warning),fixed=TRUE)) {
      failed <<- TRUE;invokeRestart('muffleWarning')
    }
  })
  if(failed) {
    options$method <- 5L;options$control.solver <- controlSolverGen('eigen',5L)
    value <- do.call(lmRidge,options)
  }
  value
}
