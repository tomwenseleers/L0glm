#' Final local support refinement
#'
#' @name local_search
#' @md
#' @description
#' Optional finite support refinement after [l0glmFit()] or [l0nnPois()].
#' [l0glmTune()] applies the same procedure independently in each full-data
#' candidate and each cross-validation training fit.
#'
#' @details
#' Set `local.search=TRUE` to run additions, deletions, single exchanges and
#' small batch exchanges after the initial adaptive-ridge fit. SIS/ISIS,
#' compulsory-term fitting and any original-channel preclustering refit finish
#' first. Their internal fits do not run local search. `FALSE` is the default.
#' `local.search="auto"` runs refinement when the initial support has at most
#' 50 free coefficient coordinates, including mandatory terms. A grouped
#' multichannel row contributes `m` coordinates; a multinomial row contributes
#' `m-1`. Ungrouped fits count active entries, subtracting the CLR row constraint
#' for multinomial. Noise parameters are excluded from this computational size.
#' The decision is taken after screening and refitting, separately for each
#' training model. Increase `control.local$size.max` to change the limit.
#'
#' The default target is
#' \deqn{Q(B,\psi)=-2\ell(B,\psi)+\lambda d(B,\psi),}
#' where `lambda` is the fitted model's initial penalty and `d` uses the
#' parameter-count convention in [model_criteria]. Every accepted step strictly
#' improves its target on all training observations, with the original predictor
#' count. A supplied `validation` set instead selects moves by its held-out
#' deviance, consistently with validation-based stopping. CV validation folds
#' never enter the training searches.
#'
#' Current-model and removal-conditioned native efficient scores shortlist
#' candidates. Small supports use full restricted refits. Larger supports use
#' small conditional blocks, holding other fitted terms in a link-scale offset,
#' followed by a full-support refit. Independent response channels only need new
#' coefficient solves when their support changes. Sparse observations outside
#' a conditional block can be omitted from its coefficient solve for fixed
#' noise parameters; full-data likelihood scoring and Cox risk sets are retained.
#'
#' Explicit starting-zero exclusions and mandatory terms retain their meaning.
#' Multinomial refits respect CLR zero sums and entry masks. Ordinary and NB
#' refits use the existing family engines; nonnegative identity-Poisson fits
#' retain the native [l0nnPois()] backend. Free Gaussian variance is profiled for
#' comparison, NB theta is re-estimated when requested, and quasi scale is held
#' common across candidates. The original-dimension sparse coefficients are
#' returned.
#'
#' `control.local` is a named list. Main controls are:
#' \describe{
#' \item{`size.max`}{Maximum initial coefficient dimension for `"auto"`,
#' default 50. This limits automatic refinement independently of IC counts.}
#' \item{`maxit`, `max.evals`}{At most 12 search rounds and 160 candidate solves
#' by default. The final full-support polish may add one solve.}
#' \item{`top`, `drop`, `max.swap`}{Two shortlisted additions, removals and
#' simultaneous exchange terms by default.}
#' \item{`scope`}{`"auto"` reuses the retained SIS/ISIS pool when available;
#' otherwise all eligible predictors. `"screened"` explicitly requests that
#' pool and `"all"` permits a full eligible-predictor scan.}
#' \item{`block.max`}{32 coefficient rows per conditional trial for supports
#' exceeding 128 rows. `Inf` requests full-support trials.}
#' \item{`ncores`, `nthreads`}{Outcome-refitting and native score-screening thread
#' counts. Defaults inherit supplied fitting/screening counts, with one screening
#' thread otherwise.}
#' \item{`refit.maxit`, `refit.tol`, `tol`}{75 refit iterations, `1e-8` refit
#' tolerance and `1e-8` relative target-improvement tolerance.}
#' \item{`mu.floor`}{Scoring-only positive-mean floor, default `1e-4`.
#' Tuning inherits its `mu_floor` unless this control is supplied explicitly.
#' Custom scalar-family likelihoods retain their exact boundary conventions.}
#' \item{`target`, `criterion`}{`target="IC"` uses parameter counts; `"L0"`
#' charges selected rows for grouped multichannel fits and nonzero entries
#' otherwise. `criterion=NULL` retains the fit's numeric penalty; an IC name
#' such as `"bic"` requests that preset for refinement.}
#' \item{`conditional`, `adjust.max`, `entry.moves`}{Removal-conditioned scans
#' are enabled, small active information blocks are limited to 24 rows, and
#' multichannel entry moves are enabled. `conditional=FALSE, adjust.max=0`
#' requests a cheaper diagonal-score search.}
#' \item{`trace`}{Print accepted search steps, default `FALSE`.}
#' }
#'
#' The returned `local.search` component contains initial and final scores,
#' accepted history, candidate counts, timings, stopping reason and compact
#' `initial.fit` coefficients, adaptive penalties and noise parameters. These
#' preserve the adaptive-ridge state separately from the refined estimates.
#' A bounded search reports its tested neighbourhood and budget; it does not
#' certify a global optimum. Lower training IC can trade additional false
#' positives for better prediction, so refinement remains optional.
#'
#' @seealso [l0glmFit()], [l0nnPois()], [l0glmTune()], [model_criteria]
#' @examples
#' set.seed(42)
#' X <- matrix(rnorm(120 * 12), 120, 12)
#' y <- 2 * X[, 1] - X[, 2] + rnorm(120)
#' fit <- l0glmFit(X, y, lambda="gic", local.search=TRUE,
#'                 parallelize=FALSE)
#' fit$local.search$history
#' automatic <- l0glmFit(X, y, local.search="auto", parallelize=FALSE)
#' automatic$local.search[c("applied", "initial.size", "stopping")]
NULL

L0glm_local_mode <- function(value) {
  if (is.logical(value) && length(value)==1L && !is.na(value)) return(value)
  if (is.character(value) && length(value)==1L && !is.na(value) &&
      tolower(value)=="auto") return("auto")
  stop('local.search must be TRUE, FALSE or "auto".', call.=FALSE)
}

L0glm_local_size <- function(fit) {
  beta <- coef(fit)
  m <- ncol(beta)-as.integer(identical(fit$kind,"multinomial"))
  selected <- setdiff(which(Matrix::rowSums(beta!=0)>0),fit$always.include)
  slopes <- if(isTRUE(fit$grouped)) length(selected)*m else
    sum(beta[selected,,drop=FALSE]!=0)-
      if(identical(fit$kind,"multinomial"))length(selected)else 0L
  slopes+length(fit$always.include)*m
}

# The outer call owns refinement. Capture only supplied arguments so the
# existing missing-argument defaults (especially automatic ISIS) remain intact.
L0glm_local_fit <- function(fit.fun, arguments, mode, control, nonnegative=NULL) {
  control <- L0glm_local_control(control, arguments)
  arguments$local.search <- FALSE
  started <- proc.time()[[3L]]
  fit <- do.call(fit.fun, arguments)
  elapsed <- proc.time()[[3L]]-started
  if (isFALSE(mode)) return(fit)
  size <- L0glm_local_size(fit)
  if (identical(mode,"auto") && size>control$size.max) {
    fit$local.search <- list(requested=mode,applied=FALSE,initial.size=size,
      size.max=control$size.max,stopping="size_limit",initial.seconds=elapsed,
      seconds=0,total.seconds=elapsed)
    return(fit)
  }
  original <- fit
  if (!is.null(nonnegative)) arguments$nonnegative <- nonnegative
  if (identical(fit$kind,"multinomial")) arguments$nonnegative <- FALSE
  arguments <- L0glm_local_arguments(fit, arguments)
  fit <- L0glm_local_refine(fit, arguments, control)
  fit$local.search$requested <- mode
  fit$local.search$applied <- TRUE
  fit$local.search$initial.size <- size
  fit$local.search$size.max <- control$size.max
  fit$local.search$initial.seconds <- elapsed
  fit$local.search$initial.fit <- original[c("coefficients","lambdas","dispersion",
    "dispersion.estimated","theta","theta.estimated","init.lambda","preprocessing")]
  fit <- L0glm_local_preprocessing(fit, arguments)
  fit$local.search$total.seconds <- proc.time()[[3L]]-started
  fit
}

L0glm_local_arguments <- function(fit, args) {
  # An aggregate-only result has its own response/offset/noise units. Preserve
  # structural start exclusions by projecting their eligibility, not their signs.
  cl <- fit$preclustering
  if (isTRUE(cl$used) && !isTRUE(cl$refitted)) {
    P <- Matrix::sparseMatrix(i=seq_along(cl$clusterids),j=cl$clusterids,x=1,
      dims=c(cl$original_channels,cl$clusters))
    if (identical(fit$family$family,"gaussian"))
      P <- P %*% Matrix::Diagonal(x=1/tabulate(cl$clusterids,cl$clusters))
    if (!is.null(args$start) && length(args$start)>1L) {
      S <- as.matrix(args$start)
      if (ncol(S)==1L) S <- matrix(S,nrow(S),cl$original_channels)
      args$start <- 1*(as.matrix((S!=0)%*%P)>0)
    }
    args$validation <- L0glm_precluster_validation(args$validation,group=cl$clusterids,
      groups=P,family=fit$family,dispersion=fit$dispersion)
  }
  args
}

L0glm_local_preprocessing <- function(fit, args) {
  rows <- sort(unique(c(fit$always.include,L0glm_local_pattern(fit)$rows)))
  if (!length(rows)) return(fit)
  lower <- rep(if(isTRUE(args$nonnegative))0 else -.Machine$double.xmax,length(rows))
  if (fit$family$family%in%c("binomial","quasibinomial"))
    lower[L0glm_intercepts(fit$X[,rows,drop=FALSE])] <- -.Machine$double.xmax
  spec <- L0glm_preprocessing(fit$X[,rows,drop=FALSE],
    if(is.null(args$normalize))"auto" else args$normalize, lower,
    match(fit$always.include,rows))
  for (name in c("center","scale")) fit$preprocessing[[name]][rows] <- spec[[name]]
  fit
}

L0glm_local_control <- function(control, args=list()) {
  defaults <- list(size.max=50L,top=2L,drop=2L,maxit=12L,max.evals=160L,max.swap=2L,
    nthreads=if(is.null(args$screen.nthreads))1L else args$screen.nthreads,
    ncores=if(is.null(args$ncores))L0glm_default_ncores(!isFALSE(args$parallelize))else args$ncores,
    adjust.max=24L,scope="auto",target="IC",criterion=NULL,
    refit.maxit=75L,refit.tol=1e-8,tol=1e-8,mu.floor=1e-4,
    conditional=TRUE,entry.moves=TRUE,channel.refit=TRUE,refit.every.round=FALSE,
    native.score=TRUE,incremental=TRUE,block.max=32L,final.refit=TRUE,prune.rows=TRUE,trace=FALSE)
  if (!is.list(control) || (length(control) && (is.null(names(control)) ||
      any(!nzchar(names(control))) || anyDuplicated(names(control)))) ||
      any(!names(control)%in%names(defaults))) stop("Use named control.local entries; see ?local_search.",call.=FALSE)
  defaults[names(control)] <- control
  for (name in c("top","drop","maxit","max.evals","max.swap","nthreads","ncores","refit.maxit","adjust.max","size.max")) {
    value <- defaults[[name]]
    if (!is.numeric(value) || length(value)!=1L || !is.finite(value) ||
        value<as.integer(!name%in%c("adjust.max","size.max")) || value!=floor(value))
      stop("Invalid control.local integer: ",name,call.=FALSE)
  }
  for (name in c("tol","mu.floor","refit.tol")) {
    value <- defaults[[name]]
    if (!is.numeric(value) || length(value)!=1L || !is.finite(value) ||
        value<0 || name=="refit.tol" && value==0)
      stop("Invalid control.local tolerance: ",name,call.=FALSE)
  }
  for (name in c("conditional","entry.moves","channel.refit","refit.every.round",
      "native.score","incremental","final.refit","prune.rows","trace"))
    if (!is.logical(defaults[[name]]) || length(defaults[[name]])!=1L || is.na(defaults[[name]]))
      stop("control.local$",name," must be TRUE or FALSE.",call.=FALSE)
  b <- defaults$block.max
  if (!is.numeric(b) || length(b)!=1L || is.na(b) || b<1 || is.finite(b) && b!=floor(b))
    stop("control.local$block.max must be a positive integer or Inf.",call.=FALSE)
  defaults$scope <- match.arg(defaults$scope,c("auto","screened","all"))
  defaults$target <- match.arg(defaults$target,c("IC","L0"))
  if (!is.null(defaults$criterion)) L0glm_resolve_lambda(defaults$criterion,100L,100L)
  defaults
}
