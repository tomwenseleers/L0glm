#' Efficient-score screening for GLMs and special families
#' @md
#'
#' @description
#' Screen predictors using one null/nuisance fit and native column reductions,
#' without fitting a separate model for every predictor. Dense matrices stay
#' dense and sparse matrices stay sparse.
#' @param X Numeric predictor matrix, dense or sparse, with observations in rows.
#' @param y Response as in [L0glm.fit()]: a numeric vector or independent-channel
#'   matrix for ordinary/NB models; factors, class labels or count matrices for
#'   multinomial; or `survival::Surv(time, status)` for right-censored Cox models.
#' @param family Any standard GLM family/link, a fixed-theta negative-binomial
#'   family, or `mgcv::nb()`, `mgcv::multinom(K)` or `mgcv::cox.ph()`.
#' @param weights Nonnegative prior weights: a scalar, observation vector, or
#'   observation-by-channel matrix for ordinary/NB models. Multinomial and Cox
#'   use shared row weights.
#' @param offset Known linear-predictor offsets: a scalar, observation vector,
#'   or observation-by-channel/category matrix. Cox uses a vector.
#' @param nsis Number of predictors retained, capped at `ncol(X)`. Defaults to
#'   `min(p, max(1, floor(n/log(n))))`. No size guarantees true-feature retention.
#' @param return_stats Return all statistics and their full ranking, default
#'   `TRUE`. Use `FALSE` for large designs to retain only the top predictors.
#' @param nuisance Optional numeric nuisance design. `NULL` uses an intercept,
#'   except for Cox, which has none. An explicit zero-column matrix requests no
#'   nuisance terms. Include an intercept explicitly when supplying other
#'   nuisance predictors. Constant Cox nuisance columns are omitted.
#' @param nonnegative Use one-sided scores for nonnegative candidate effects.
#'   Nuisance terms remain unrestricted. Ignored for multinomial CLR effects.
#' @param iterative Use offset-based iterative SIS (ISIS), default `FALSE`.
#' @param niter Maximum ISIS rounds, default 10. Stops when the candidate
#'   set repeats, since refitting that same subproblem would give the same fit.
#' @param nthreads Positive integer OpenMP thread count for predictor screening,
#'   `NULL` (default) uses the physical core count, falling back to one if
#'   detection is unavailable. Capped at logical processors and predictor count.
#'   Builds without OpenMP use one thread. Does not change BLAS threads or
#'   outcome-fitting parallelism; use `screen.nthreads` in [L0glm.fit()].
#' @param ... Named [L0glm.fit()] arguments for ISIS, such as `lambda`, `start`,
#'   `intercept`, `always.include`, `grouped`, `normalize` and solver controls.
#'   ISIS defaults to `lambda="gic"`, `start=NULL`, `relaxed=TRUE`
#'   and `parallelize=FALSE`. Coefficient-shaped starts/prior means use the
#'   original predictor dimensions, optionally including an added intercept;
#'   they are subset internally. Scalar `start=1` retains its normalized-coordinate
#'   meaning. Initial zero starts remain excluded, except mandatory terms.
#' @details
#' For an ordinary GLM let \eqn{s_i=v_i(y_i-\mu_i)\dot\mu_i/
#' (\phi V(\mu_i))} and \eqn{h_i=v_i\dot\mu_i^2/(\phi V(\mu_i))}, evaluated
#' at the nuisance fit. For candidate \eqn{x} and nuisance design \eqn{N},
#' \deqn{U=x^Ts-x^THN(N^THN)^+N^Ts,\qquad
#' I=x^THx-x^THN(N^THN)^+N^THx,}
#' with \eqn{H=\operatorname{diag}(h)}. Rank \eqn{U^2/I}, or
#' \eqn{\max(U,0)^2/I} for a one-sided candidate. Unidentifiable columns get
#' zero; ties use original column order. Noncanonical derivatives retain their
#' sign. Free dispersion uses the null-fit Pearson statistic/residual degrees
#' of freedom. NB theta is fitted or fixed according to the supplied family.
#'
#' Independent channels contribute separate statistics, preventing cancellation
#' of opposite effects. Multinomial uses the exact joint score/information in
#' orthonormal contrasts of the CLR plane; grouped ranking is invariant to
#' category relabeling. Cox uses exact weighted Breslow risk sets, ties and
#' nuisance adjustment. Quasi families do not require a likelihood.
#'
#' Native kernels avoid `X^2` and scaled design copies. Independent channels
#' use diagonal information; multinomial accumulates curvature from
#' probabilities without an observation-by-category-squared array. Without
#' full statistics, a bounded heap selects top predictors. Screening can miss
#' weak/suppressor effects and does not certify a global L0 optimum.
#' Screening threads process separate predictors with unchanged within-column
#' summation order. Thread-local top-candidate heaps are merged deterministically,
#' including ties in original column order. Errors and interrupts are handled on
#' the main R thread; no design-sized copies are made for threading.
#'
#' With `iterative=TRUE`, add the fitted linear predictor to the original
#' offset for screening, merge new candidates with current nonzero/mandatory
#' rows, and refit on the original responses and original offsets. Each round
#' starts afresh; this does not reintroduce coefficients within an individual
#' BAR fit. Preset lambda and GIC use the original selectable predictor count.
#' The lowest-GIC round is returned, using the metrics scoring mean floor of
#' `1e-8`. For quasi families without a likelihood/GIC, return the last round.
#' `nuisance` must be `NULL` for ISIS: put mandatory covariates in `X` and pass
#' `always.include` instead. Intercept defaults match [L0glm.fit()].
#' `converged` describes repetition of the outer candidate set; inner-fit
#' convergence diagnostics remain in `fit$convergence`.
#' [L0glm.fit()] also exposes this workflow through `screen="SIS"` or
#' `screen="ISIS"`, with automatic ISIS for ultra-wide eligible designs.
#' Its small default ISIS batch and other screened-fit defaults are documented
#' in that function; pass `nsis` explicitly to match a standalone call exactly.
#' @return A list with decreasing-score `indices`, `top_statistics`, `nsis`,
#'   `family`, `link`, fitted `theta` for NB, and `nthreads` (screening threads
#'   actually used). With `return_stats=TRUE`, `ranking` lists all predictor
#'   indices and `stats` contains `feature`,
#'   `stat`, `sign`, `rank` and `keep`. Signs are `NA` for multichannel/joint
#'   scores. Full statistics/ranking are otherwise `NULL`.
#'   ISIS additionally returns a full-dimension `L0glm_fit` in `fit`, its sparse
#'   `coefficients`, the best round's `candidates`, resolved `lambda`, `history`
#'   (including GIC, deviance and support size), `best_iteration`, `iterations`,
#'   `converged` and elapsed `seconds`. Its `indices` are retained nonzero or
#'   mandatory predictors in original column order; `nsis` counts these and
#'   `screen_size` records the screening batch size. Optional statistics/ranking
#'   describe the screen preceding the chosen fit, with `keep` updated to its
#'   retained rows. `top_statistics` is `NULL` for ISIS. The returned `fit`
#'   accepts new data with all original columns in [stats::predict()].
#' @seealso [L0glm.fit()]
#' @example examples/sis_glm.R
#' @export
sis_glm <- function(X, y, family=gaussian(), weights=NULL, offset=NULL,
                    nsis=NULL, return_stats=TRUE, nuisance=NULL, nonnegative=FALSE,
                    iterative=FALSE, niter=10L, nthreads=NULL, ...) {
  X <- L0glm_sis_design(X)
  if (!inherits(family,"family")) family <- family()
  if (!is.logical(return_stats) || length(return_stats)!=1L || is.na(return_stats) ||
      !is.logical(nonnegative) || length(nonnegative)!=1L || is.na(nonnegative) ||
      !is.logical(iterative) || length(iterative)!=1L || is.na(iterative))
    stop("return_stats, nonnegative and iterative must be TRUE or FALSE.",call.=FALSE)
  n <- nrow(X);p <- ncol(X)
  if (is.null(nsis)) nsis <- max(1L,floor(n/log(max(n,3))))
  if (!is.numeric(nsis) || length(nsis)!=1L || !is.finite(nsis) ||
      nsis<1 || nsis!=floor(nsis)) stop("nsis must be a positive integer.",call.=FALSE)
  nsis <- as.integer(min(nsis,p))
  nthreads <- L0glm_sis_threads(nthreads)
  if (iterative) {
    if (!is.null(nuisance)) stop("For ISIS, put nuisance covariates in X and use always.include.",call.=FALSE)
    return(L0glm_isis(X,y,family,weights,offset,nsis,return_stats,nonnegative,niter,list(...),nthreads))
  }
  if (length(list(...))) stop("Additional fitting arguments require iterative=TRUE.",call.=FALSE)
  context <- L0glm_sis_context(y,family,n,weights,offset,nuisance,nonnegative)
  result <- sis_score_cpp(X,context,nsis,return_stats,nthreads)
  ranking <- stats <- NULL
  if (return_stats) {
    ranking <- order(result$statistics,decreasing=TRUE,method="radix")
    stats <- data.frame(feature=seq_len(p),stat=result$statistics,sign=result$sign,
      rank=match(seq_len(p),ranking),keep=seq_len(p)%in%result$indices)
  }
  list(indices=result$indices,top_statistics=result$top_statistics,
    ranking=ranking,stats=stats,nsis=length(result$indices),family=family$family,
    link=family$link,theta=context$theta,nthreads=result$nthreads)
}

L0glm_sis_threads <- function(nthreads) {
  if (is.null(nthreads)) {
    nthreads <- L0glm_detect_physical_cores()
    if (is.na(nthreads) || nthreads < 1L) nthreads <- 1L
  }
  if (!is.numeric(nthreads) || length(nthreads)!=1L || !is.finite(nthreads) ||
      nthreads<1 || nthreads!=floor(nthreads) || nthreads>.Machine$integer.max)
    stop("Screening nthreads must be a positive integer.",call.=FALSE)
  as.integer(nthreads)
}

L0glm_sis_design <- function(X) {
  if (inherits(X,"sparseMatrix"))
    X <- methods::as(methods::as(methods::as(X,"dMatrix"),"generalMatrix"),"CsparseMatrix")
  else {
    X <- as.matrix(X)
    if (!is.numeric(X) && !is.logical(X)) stop("X must be numeric.",call.=FALSE)
    if (!is.double(X)) storage.mode(X) <- "double"
  }
  if (length(dim(X))!=2L || nrow(X)<2L || ncol(X)<1L)
    stop("X needs at least two observations and one predictor.",call.=FALSE)
  X
}
