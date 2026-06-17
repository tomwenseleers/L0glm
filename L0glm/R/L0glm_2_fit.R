#' @name L0glm.fit
#' @aliases L0glm.fit
#' @title Fit generalized linear model with \eqn{L_0} or grouped
#'   \eqn{L_{0,\infty}} regularization subject to optional nonnegativity or box
#'   constraints
#'
#' @description Fits a generalized linear model with \eqn{L_0} norm or grouped
#' \eqn{L_{0,\infty}} mixed norm regularization to approximate Best Subset
#' Selection (for univariate data) or Best Subset of Group Selection (for
#' multivariate data) subject to optional nonnegativity or box constraints on
#' model coefficients. The objective functions that can be approximated are :
#' \itemize{
#'
#' \item for univariate data, using an \eqn{L_0} norm penalty to achieve Best Subset Selection:
#'
#' \deqn{-2.LL + lambda.||beta||_0}
#'
#' where \eqn{LL} is the log-likelihood of the GLM model and
#' \eqn{|beta||_0} is the number of non-zero elements in
#' coefficient vector beta, i.e. the \eqn{L_0} norm of beta.
#'
#' \item for multivariate data, using an \eqn{L_{0,0}} norm penalty, i.e. using independent \eqn{L_0}  penalized fits on each outcome:
#'
#' \deqn{-2.LL + lambda.||beta||_{0,0}}
#'
#' \item for multivariate data, using an \eqn{L_{0,\infty}} norm penalty, i.e. using grouped variable selection in multitask learning settings to achieve Best Subset of Groups Selection:
#'
#' \deqn{-2.LL + lambda.||beta||_{0,\infty}}
#'
#' }
#'
#' in all 3 cases subject to optional box constraints on the model coefficient
#' vector or matrix beta
#'
#' \deqn{lower <= beta <= upper}
#'
#' These objectives are approximated via an iterative adaptive ridge regression
#' procedure, see Details. The \eqn{L_0} and \eqn{L_{0,0}} norm penalties can be
#' used to approximate Best Subset Selection, while for multivariate
#' (multitask/multichannel) data, the \eqn{L_{0,\infty}} mixed norm penalty can be
#' used to approximate Best Subset of Group Selection (via option
#' \code{grouped=TRUE}). Fits on multichannel data can be parallelized. Option
#' \code{relaxed} allows one to finish with a (nearly) unpenalized fit on the
#' \code{L0glm} selected features.
#'
#' The level of regularisation can be manually chosen (if a numeric value is
#' passed to \code{lambda}) or preset to optimise given information criteria (if
#' set to \code{"aic"}, \code{"bic"}, \code{"gic"}, \code{"ebic"} or
#' \code{"mbic"}). For example, a lambda of 2 is expected to optimise AIC, a
#' lambda of log(\emph{n}) BIC, etc.
#'
#' The \code{L0glm.fit} function has various advanced features, including
#' support for both dense and sparse design matrices \emph{\bold{X}}, high
#' dimensional data (where \emph{p} >> \emph{n}), single or multichannel
#' outcomes \emph{\bold{y}} (with either independent or grouped variable
#' selection), support for any of the standard R GLM families & link functions,
#' prior weights (e.g. total nr. of trials for binomial proportions), model
#' offsets, starting coefficient values or dispersion (if prior estimate is
#' available, otherwise optimised), nonnegativity or box constraints,
#' unpenalized terms (e.g. intercept) and shrinkage towards given prior means
#' (via nonzero centered ridge regression updates).
#'
#' A choice of several fast C++ coded least squares solvers can be used, which
#' means this function scales well to solve fairly large tall or wide problems
#' with millions of observations or variables. As solvers, we currently
#' implemented the \code{Eigen} sparse or dense Cholesky, QR, PCG & LSCG
#' solvers, plus two solvers that natively support box constraints:
#' \code{glmnet} (\code{bigGLM}), which uses coordinate descent, and
#' \code{osqp}, a fast quadratic programming solver.
#'
#' @param X Design matrix with predictor variables of dimension \code{n} x
#'   \code{p}, with \code{n} corresponding to the number of observations and
#'   \code{p} being the number of columns in \code{y} (variables plus
#'   intercept). It should include a column of ones for an intercept term is
#'   desired.
#' @param y Response variable(s). Can be either a vector of length \code{n} or a
#'   matrix of dimension \code{n} x \code{m}, where \code{n} is the number of
#'   observations and \code{m} the number of response variables
#'   (tasks/channels).
#' @param family a description of the error distribution (gaussian, poisson,
#'   binomial, Gamma, inverse.gaussian) and link function (identity, log, sqrt,
#'   logit, cloglog, probit, cauchit, inverse, 1/mu^2) to be to be used for the
#'   GLM model. See \code{?stats::family} for details of family functions.
#' @param weights Optional vector of prior observation weights ("prior weights")
#'   of length \code{n}, or a matrix of weights with dimension \code{n} x
#'   \code{m}, with \code{n} corresponding to the number of observations and
#'   \code{m} the number of response variables (columns in \code{y}). For
#'   binomial proportion data this would be the total number of trials per
#'   observed proportion (passed in \code{y}).
#' @param dispersion Optional dispersion estimate; set to \code{NULL} by
#'   default, in which case it is initialised at 1 and then estimated during the
#'   IRLS loop; fixed at 1 for family binomial & Poisson.
#' @param start Optional vector of coefficient starting values of length
#'   \code{p}, or a matrix of starting values with dimension \code{p} x
#'   \code{m}, with \code{p} corresponding to the number of predictors and
#'   \code{m} to the number of columns in \code{y}. These will determine the
#'   level of penalisation of the different coefficients during the first IRLS
#'   ridge regression update. The ridge penalty used equals
#'   \eqn{(\lambda/4) \phi / |\beta|^\gamma}. Note that coefficients
#'   with absolute value smaller than \code{thresh} will be clipped to zero and
#'   are removed from the model. Set to \code{NULL} by default, in which case it
#'   is auto-initialised.
#' @param offset Optional vector of offset values of length \code{p}, or a
#'   matrix of offset values with dimensions \code{c(p, m)}. The offset is a
#'   known quantity that is added to the linear predictor of the model.
#' @param lambda Either a scalar with the \eqn{L_0} or \eqn{L_{0,0}} penalty
#'   (when \code{grouped=FALSE}) or \eqn{L_{0,\infty}} penalty (when
#'   \code{grouped=TRUE}), or a character with the name of an information
#'   criterion (\code{"aic"}, \code{"bic"}, \code{"gic"}, \code{"ebic"} or
#'   \code{"mbic"}), in which case \code{lambda} will be preset with a value
#'   that would be expected to optimise that particular information criterion.
#'   Set \code{lambda = 0} to disable penalization.
#' @param intercept Logical: Whether to include an intercept term in the model.
#'   If `TRUE`, a column of ones will be added to the design matrix `X`, and
#'   this term will **not** be penalized. If `FALSE`, no intercept column will
#'   be added, and the user must ensure that `X` already contains an intercept
#'   if needed. Defaults to `FALSE`.
#' @param always.include A vector of coefficient indices that should always be
#'   estimated & included in the model. This is done by updating \code{prior.mean}
#'   in each IRLS iteration for these coefficients, instead of using sceptical
#'   zero-centered prior updates for the other coefficients. Coefficients with 
#'   a name that includes \code{"Intercept"} in their name are always included. 
#'   Specify 0 if all terms should be allowed to be
#'   shrunk to zero. If `intercept = TRUE`, the first column (intercept) will
#'   always be estimated and included in the model.
#' @param nonnegative A logical indicating whether nonnegativity constraints
#'   should be applied. If \code{TRUE}, \code{lower} will be set to \code{rep(0,
#'   ncol(X))}.
#' @param normalize Logical indicating whether the columns of the design matrix
#'   \strong{X} should be scaled to have unit \eqn{L_2} norm. It is advised to
#'   use normalization unless the covariates were measured on a similar scale.
#'   Defaults to \code{TRUE}
#' @param solver Solver used for the ridge regression step in the IRLS algorithm
#'   (see details). It should be one of \code{"osqp"} (osqp quadratic
#'   programming solver), \code{"glmnet"} (\code{"bigGlm"}, using coordinate
#'   descent) or \code{"eigen"} (one of the \code{RcppEigen Eigen} solvers).
#'   \code{"osqp"} and \code{"glmnet"} support box or nonnegativity constraints
#'   and scale well to larger problems, with especially \code{"glmnet"} being
#'   suitable for very large problems. The \code{"eigen"} solvers achieve box or
#'   nonnegativity constraints via clipping, which in combination with the ridge
#'   penalties in practice work almost as well, and some of those solvers are
#'   often among the fastest (e.g. \code{method} \code{3L}, \code{5L} or
#'   \code{7L}, i.e. LSCG, CG or Cholesky LLT).
#' @param method When \code{solver="eigen"}, method used for the eigen solver
#'   (see details). It should be an integer from 1 to 8 :
#'
#'   \code{1L} = QR solver, Householder QR solver if \code{X} is dense
#'
#'   \code{2L} = QR solver, QR solver with column pivoting if \code{X} is dense
#'
#'   \code{3L} = Least Squares Conjugate Gradient (LSCG) with identity
#'   preconditioner (i.e. not preconditioned)
#'
#'   \code{4L} = Least Squares Conjugate Gradient (LSCG) with least squares
#'   diagonal preconditioner
#'
#'   \code{5L} = Conjugate Gradient (CG) with identity preconditioner (i.e. not
#'   preconditioned)
#'
#'   \code{6L} = Preconditioned Conjugate Gradient (PCG) with diagonal
#'   preconditioner
#'
#'   \code{7L} = Cholesky LLT solver
#'
#'   \code{8L} = Cholesky LDLT solver
#'
#'   Defaults to \code{7L} (Cholesky LLT), which is often the fastest.
#'
#' @param lower Optional lower limit on model coefficients to specify box
#'   constraints. It should either be a scalar, a vector of length \code{p} or a
#'   matrix with dimensions \code{c(p, m)}.
#' @param upper Optional upper limit on model coefficients to specify box
#'   constraints. It should either be a scalar, a vector of length \code{p} or a
#'   matrix with dimensions \code{c(p, m)}.
#' @param parallelize A logical value indicating whether in the multivariate
#'   case, the columns of the response matrix \code{y} should be fitted in
#'   parallel. If \code{TRUE}, each column of \code{y} is fitted concurrently,
#'   potentially speeding up computation. If \code{FALSE}, the columns are
#'   fitted sequentially.
#' @param ncores Optional integer specifying the number of cores to use when
#'   \code{parallelize=TRUE}. If \code{parallelize=TRUE}, defaults to the
#'   detected number of physical cores when that can be determined, falling back
#'   to half the number of logical cores. The requested value is capped at the
#'   detected number of physical cores when that can be determined.
#' @param nthreads Optional integer specifying the number of BLAS/OpenMP threads
#'   to use per fit process. The default \code{NULL} leaves the main R session's
#'   thread settings unchanged and uses one BLAS/OpenMP thread per parallel
#'   worker. Setting \code{nthreads} can be useful to benchmark serial fits or
#'   tune per-worker threading on systems with multithreaded BLAS/OpenMP.
#' @param maxit Integer specifying the maximum number of IRLS iterations to
#'   perform.
#' @param grouped A logical value indicating whether a grouped fit should be
#'   performed in the multivariate case or not. If \code{TRUE}, the algorithm
#'   will return a group penalized fit, otherwise independent L0glm fits will be
#'   performed on each response. These procedures approximate the \eqn{L_{0,\infty}}
#'   and \eqn{L_{0,0}} mixed norm penalties, and can be used to do Best Subset
#'   of Group Selection or plain Best Subset Selection.
#' @param group.maxit Integer  specifying the maximum number of iterations to
#'   perform for the group penalized fitting procedure when \code{grouped =
#'   TRUE}. In each iteration, the columns of the response matrix \code{y} will
#'   be refitted using the average of all columns of \code{y} as the starting
#'   value, which then imposes a group penalty.
#' @param relaxed Logical indicating if the algorithm should conduct a nearly
#'   unpenalized fit on the L0glm selected variables (using
#'   \code{lambda=1E-10}). This can eliminate any remaining bias on the
#'   coefficients and increase sensitivity, especially for the multivariate case
#'   with grouped variable selection. Defaults to \code{TRUE} when
#'   \code{grouped=TRUE} & multichannel data are passed in \code{y}.
#' @param warm.start A numeric value in the range of [0, 1] that represents the
#'   level of confidence in the initial start values when provided. This
#'   parameter is used in a weighted average to calculate the 'start' values.
#'   Default is set at 0.9. The formula used is: \deqn{start = \text{warm.start}
#'   \cdot start + (1 - \text{warm.start}) \cdot
#'   \text{rep}(\text{sd}(start)^{\frac{1}{\gamma}}, \text{ncol}(X))} If
#'   `warm.start` is closer to 1, more weight is given to the original 'start'
#'   value. Using this weighted average avoids the problem that coefficients
#'   initialized at zero would otherwise always be removed from the model.
#' @param prior.mean An optional numeric vector of length \code{p} specifying
#'   the center towards which coefficients should be shrunk in each IRLS ridge
#'   regression iteration. This is done via nonzero centered ridge regression.
#'   This can be useful if there is good prior information available that some
#'   coefficients should have a particular nonnegative value. The default
#'   \code{prior.mean = NULL} will shrink to zero using a default zero centered
#'   ridge regression.
#' @param gamma Scalar with the exponent to use when calculating the ridge
#'   penalty \code{(lambda/4)*dispersion/abs(beta)^gamma)} in each iteration.
#'   Defaults to \code{gamma = 2}, as in Frommlet & Nuel 2016.
#' @param thresh Coefficients below this scalar threshold are clipped to 0 in
#'   each of the IRLS ridge regression iterations and removed from the model.
#'   Defaults to \code{1E-5}.
#' @param tol Convergence tolerance. Iterations are stopped when the relative
#'   change in the \eqn{L_2} norm of the coefficients drops below this
#'   threshold.
#' @param control.solver List of parameters controlling the optimization
#'   solvers. Auto-initialised with \code{control.solver.gen} when \code{NULL}.
#'   See \code{?control.solver.gen} for more details on the available parameters
#'   for the different solvers.
#' @param verbose Logical indicating if the algorithm progression should be
#'   printed to the console.
#'
#' @details L0glm approximates the \eqn{L_0} or \eqn{L_{0,\infty}} mixed norm
#' penalized GLM objectives via an iterative adaptive ridge regression
#' procedure.
#'
#' \strong{Algorithm for univariate, single-channel data}
#'
#' For univariate data, the approach is inspired by earlier work on L0 penalized
#' regression by de Rooi & Eilers (2011), Frommlet & Nuel (2016), Liu et al.
#' (2017), Dai et al. (2018 & 2020) & Li et al. (2021) and is based in a simple
#' modification of the IRLS (Iteratively Reweighed Least Squares) algorithm used
#' to fit GLMs (cf. McCullach & Nelder 1989, Hastie & Pregibon 1992, Hardin &
#' Hilbe 2007). In this algorithm, each iteration of the IRLS algorithm updates
#' the coefficient estimates using a quadratic approximation of the log
#' likelihood of the GLM objective. This is done using a weighted least squares
#' step. In our adapted algorithm, this step is replaced by a (zero or nonzero
#' centered) adaptive ridge regularised regression with optional box or
#' nonnegativity constraints on the coefficients. In each iteration, this uses
#' the following ridge penalties for the penalized terms (excluding the
#' intercept, if included):
#' \deqn{\Lambda = (1/4)\lambda\sigma^2 w}
#' with adaptive penalty weights
#' \deqn{w = 1 / |\beta|^\gamma}
#' where \eqn{\beta} are the coefficient estimates
#' in each iteration and coefficients smaller than \code{thresh} are set to zero
#' and taken out of the model at the start of each iteration, and as in Frommlet
#' & Nuel (2016), \eqn{\gamma} is set to a default value of 2. Box or nonnegativity
#' constraints are natively supported if the \code{glmnet} (\code{bigGlm}) or
#' the \code{osqp} solver is used, or are obtained via clipping if the
#' \code{eigen} solvers are used. The factor of 1/4 was shown to be required to
#' approximate the L0 penalty objective and has been proven for the Gaussian
#' case and under orthogonality (Frommlet & Nuel 2016), but in practice also
#' work well to the non-orthogonal case and other noise families. In contrast to
#' Frommlet & Nuel (2016) and the algorithms given in Liu et al. (2017) & Dai et
#' al. (2018, 2020) (implemented in the \pkg{l0ara} package of Guo et al. 2022),
#' however, we do not require the dispersion (variance in case of Gaussian
#' linear models) to be known and standardized to 1, but instead allow it to be
#' estimated and updated during the IRLS algorithm (when
#' \code{dispersion=NULL}). We do this by first setting it to 1 and then after a
#' few iterations, after the estimated dispersion has converged to a stable
#' value (when the relative change between iterations is < 0.1), updating the
#' dispersion multiplier that controls the ridge penalties. If the dispersion is
#' known from prior data or has been estimated before, it can be also be passed
#' as a fixed value. For binomial & Poisson models the dispersion is fixed at 1.
#' In contrast to Frommlet & Nuel (2016) we use Fisher scoring instead of
#' Newton-Raphson iterations (i.e. using the expected Hessian instead of the
#' observed one). If option \code{relaxed} is set to \code{TRUE} the algorithm
#' will refit the \code{L0glm} selected variables using a nearly unpenalized
#' model (with \code{lambda=1E-10}). This can reduce any remaining bias and
#' increase sensitivity. If an appropriate level of regularisation lambda is chosen,
#' the L0 norm penalized solution approximates Best Subset selection. To achieve
#' variable selection consistency, the \code{"bic"} \code{lambda} preset is
#' recommended for low dimensional settings, whilst \code{"gic"}, \code{"ebic"}
#' or \code{"mbic"} may work better for high dimensional settings. Setting it to
#' \code{"aic"} is recommended if predictive performance is of key interest, as
#' the AIC criterion asymptotically approximates the leave-one-out
#' cross-validation error.
#'
#' \strong{Algorithm for multivariate/multichannel data}
#'
#' For the multivariate ("multitask learning") case there are two options. First,
#' if \code{grouped=FALSE}, an independent \code{L0glm} model is fit on each of
#' the responses, which comes down to using an \eqn{L_{0,0}} mixed norm penalty.
#' Second, if \code{grouped=TRUE}, grouped variable selection will be performed
#' via an \eqn{L_{0,\infty}} mixed norm penalty. This is achieved by taking the
#' \eqn{L_{0,0}} mixed norm penalized solution, and then, for all variables that
#' have nonzero coefficients in at least on of the channels, repeating the
#' \(L_0) penalized fits but using the average of the coefficients across all
#' channels (response variables) as starting values so as to achieve a group
#' penalty and perform grouped feature selection. This is then further iterated
#' until convergence. Finally, as for single channel data, if
#' \code{relaxed=TRUE}, the algorithm can finish with a final step whereby the
#' selected variables are refit using a nearly unpenalized model (with
#' \code{lambda=1E-10}). Compared to other group penalized regression methods,
#' such as the group LASSO, which uses an \(L_0/L_2) mixed norm penalty, this
#' approach produces a sparser solution that better approximates the Best Subset
#' of Groups solution if lambda is adequately chosen.
#'
#' @return Returns a list with the following elements:
#' \item{\code{coefficients}}{a sparse dgCmatrix with the fitted coefficients
#' (single column for the univariate case). If `intercept = TRUE`, the first 
#' row is for the intercept.}
#' \item{\code{weights}}{a vector or matrix of observation weights.}
#' \item{\code{init.lambda}}{intial lambda.} \item{\code{lambdas}}{a vector or
#' matrix with the adaptive penalties lambdas.} \item{\code{family}}{The error
#' structure used to fit the data.}
#' \item{\code{z}}{a vector or matrix containing the working response in the
#' final iteration of \code{L0glm.fit}}
#' \item{\code{X}}{design matrix X, including the intercept column if added}
#' \item{\code{y}}{response vector or matrix y.}
#'
#' @author Tom Wenseleers, Laura Lanzas & Christophe Vanderaa
#'
#' @references
#' \itemize{
#' \item J. de Rooi & P. Eilers, \emph{Deconvolution of pulse trains with the L0 penalty}. Analytica chimica acta 705, 218-226 (2011).
#' \item F. Frommlet & G. Nuel, \emph{An adaptive ridge procedure for L0 regularization}. PloS one 11, e0148620 (2016).
#' \item Z. Liu, F. Sun & D.P. McGovern, \emph{Sparse generalized linear model with L0 approximation for feature selection and prediction with big omics data}. BioData mining 10, 1-12 (2017).
#' \item L. Dai, K. Chen, Z. Sun, Z. Liu & G. Li, \emph{Broken adaptive ridge regression and its asymptotic properties}. Journal of multivariate analysis 168, 334-351 (2018).
#' \item L. Dai, K. Chen & G. Li, \emph{The broken adaptive ridge procedure and its applications}. Statistica Sinica 30, 1069-1094 (2020).
#' \item N. Li, X. Peng, E. Kawaguchi, M.A. Suchard & G. Li, \emph{A scalable surrogate L0 sparse regression method for generalized linear models with applications to large scale data}. Journal of Statistical Planning and Inference 213, 262-281 (2021).
#' \item W. Guo, S. Ma & Z. Liu, \emph{l0ara: Sparse Generalized Linear Model with L0 Approximation for Feature Selection}. R package version 0.1.6 (2022). Available at \url{https://CRAN.R-project.org/package=l0ara}
#' \item J.W. Hardin & J.M. Hilbe, \emph{Generalized linear models and extensions}. (Stata press, 2007).
#' \item T.J. Hastie & D. Pregibon, \emph{Generalized linear models}. Chapter 6 of \emph{Statistical Models in S} eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole (1992).
#' \item P. McCullagh & J.A. Nelder, \emph{Generalized Linear Models}. London: Chapman and Hall (1989).
#' }
#'
#' @importFrom parallel detectCores makeCluster clusterCall clusterApply
#'   stopCluster
#'
#' @export
#' @author
#' Tom Wenseleers \email{tom.wenseleers@kuleuven.be} \\
#' Laboratory of Socioecology and Social Evolution, KU Leuven \\
#' \url{https://bio.kuleuven.be/eeb/tw/research} \\
#' ORCID: \href{https://orcid.org/0000-0002-1434-861X}{0000-0002-1434-861X}
L0glm.fit <- function(X, # n x p design matrix (without intercept)
                      y, # n x m outcome matrix
                      weights       = NULL,
                      family        = gaussian(identity),
                      dispersion    = NULL,
                      start         = NULL,
                      offset        = NULL,
                      lambda        = 1,
                      intercept     = (family$family=="binomial"), # always fit an intercept with binomial family
                      always.include        = NULL,
                      nonnegative   = FALSE,
                      normalize     = FALSE,
                      solver        = "eigen",
                      method        = 7L,
                      lower         = NULL,
                      upper         = NULL,
                      parallelize   = TRUE,
                      ncores        = L0glm_default_ncores(parallelize),
                      nthreads      = NULL,
                      maxit         = 25,
                      grouped       = TRUE,
                      group.maxit   = 6,
                      relaxed       = grouped & ncol(y)>1,
                      warm.start    = 0.9,
                      prior.mean    = NULL,
                      gamma         = 2,
                      thresh        = 1e-5,
                      tol           = 1e-10,
                      control.solver = do.call("control.solver.gen", list(solver = solver, method = method)),
                      verbose       = FALSE) {
  
  # Helpers ----------------------------------------------------------------
  issparse  <- function(M) inherits(M, "sparseMatrix")
  as_list <- function(d)
    if (is.null(d)) rep(list(NULL), m) else as.list(d)
  thread_state <- L0glm_set_thread_count(nthreads)
  if (!is.null(thread_state)) on.exit(L0glm_restore_thread_count(thread_state), add = TRUE)
  worker_nthreads <- if (is.null(nthreads)) 1L else L0glm_normalize_thread_count(nthreads)
  
  # Check dimensions ------------ -------------------------------------------
  y = as.matrix(y) # ensure y is a matrix (n x m) & coerce y to dense
  n   = nrow(X)
  p   = ncol(X)
  m   = ncol(y)
  
  # Coerce X to dgCMatrix if sparse -----------------------------------------
  if (issparse(X)) {
    X <- as(Matrix(X, sparse = TRUE), "dgCMatrix")
  } else {
    X <- Matrix(X, sparse = FALSE)
  }
  
  # Intercept handling & set lower and upper ---------------------------------
  lower <- if (nonnegative) { rep(0, ncol(X)) } else { rep(-.Machine$double.xmax, ncol(X)) }
  upper <- rep(.Machine$double.xmax, ncol(X))
    
  if (intercept) {
    X <- cbind(Intercept = 1, X)
    intercept.col <- 1  # remember its position
    always.include <- c(1, always.include+1)   # always include intercept (indexing for always.include is for non-intercept coefs)
    lower <- c(-.Machine$double.xmax, lower) # fit intercept unconstrained (important e.g. for binomial)
    upper <- c(.Machine$double.xmax, upper)
    p <- ncol(X)               # update p after intercept addition
    # TODO also check dimensions of prior.mean, weights, start & offset
    # later in L0glm_1_fit_internal.R when colnames(X)[[1]]=="Intercept"
  } else { intercept.col <- NULL }
  
  # ----- preset lambda ------------------------------------------------------
  if (lambda %in% c("aic","bic","gic","ebic","mbic")) {
    lambda = preset.lambda(IC = lambda, p = ncol(X), n = nrow(X))
  }  
  
  # - Some argument checking & conversion to matrices -----------------------
  # 1  Observation weights --------------------------------------------------
  if (is.null(weights)) {
    weights <- if (m == 1) rep(1, n) else matrix(1, n, m)
  } else {
    weights <- as.matrix(weights)
    if (m > 1 && ncol(weights) == 1L) {
      weights <- matrix(rep(weights, m), n, m)
    }
  }
  
  # 2  Offset ---------------------------------------------------------------
  if (is.null(offset)) {
    offset <- Matrix::Matrix(0, n, m, sparse = TRUE)
  } else {
    if (!inherits(offset, "Matrix"))
      offset <- Matrix::Matrix(offset, n, 1, sparse = TRUE)
    
    if (m > 1 && ncol(offset) == 1L)
      offset <- Matrix::Matrix(rep(offset, m), n, m, sparse = TRUE)
  }
  
  # 3  Starting coefficient values --------------------------------------------
  if (!is.null(start)) {
    start <- as.matrix(start)
    if (m > 1 && ncol(start) == 1L)
      start <- matrix(rep(start, m), nrow = nrow(start), ncol = m)
    if (!identical(dim(start), c(p, m)))
      stop("start must have dimension p x m (or length p for univariate)")
    # if (normalize) { # normalisation is done in L0glm.fit.internal
    #   Xscaled <- normalizeX(X, no.norm = intercept.col)
    #   scaling <- attr(Xscaled, "scaling")
    #   # start   <- start * scaling
    # }
    start <- warm.start * start + (1 - warm.start) * rep(sd(start)^(1/gamma), p)
  } else {
    start <- NULL # leave initialisation to internal function
  }
  
  # 4  Prior.mean ------------------------------------------------------------
  if (is.null(prior.mean)) {
    prior.mean <- Matrix(0, p, m, sparse = TRUE)
  } else {
    prior.mean <- as.matrix(prior.mean)
    if (m == 1) {
      if (length(prior.mean) != p)
        stop("prior.mean must have length p for univariate y")
      prior.mean <- matrix(prior.mean, p, 1)
    } else {
      if (ncol(prior.mean) == 1L) {
        prior.mean <- matrix(rep(prior.mean, m), p, m)
      } else if (!identical(dim(prior.mean), c(p, m))) {
        stop("prior.mean must be p x m when y has m columns")
      }
    }
    prior.mean <- Matrix(prior.mean, sparse = TRUE)
  }
  
  # 5 dispersion (vector per per-channel or NULL -------------
  if (!is.null(dispersion)) {
    dispersion <- as.numeric(dispersion)
    if (length(dispersion) == 1) dispersion <- rep(dispersion, m)
    if (length(dispersion) != m)
      stop("dispersion must have length 1 or length m")
  }
  disp_list <- as_list(dispersion)   # <- list(NULL) or list(dispersion vector)
  
  # 6 misc --------------------------------------------------
  if (nonnegative) lower[lower < 0] <- 0
  always.include.start <- always.include
  init.lambda <- lambda
  
  
  # =======================================================================
  #  MULTIVARIATE RESPONSE  (m > 1) ***************************************
  # =======================================================================
  if (m > 1) {
    
    if (!grouped) group.maxit <- 1L
    ncores <- suppressWarnings(as.integer(ncores[1L]))
    if (is.na(ncores) || ncores < 1L) ncores <- 1L
    physical_cores <- L0glm_detect_physical_cores()
    logical_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
    if (!is.na(physical_cores) && physical_cores >= 1L) {
      max_ncores <- as.integer(physical_cores)
      max_ncores_label <- "physical core count"
    } else if (!is.na(logical_cores) && logical_cores >= 1L) {
      max_ncores <- max(1L, as.integer(logical_cores) - 2L)
      max_ncores_label <- "detectCores(logical = TRUE) - 2"
    } else {
      max_ncores <- NA_integer_
      max_ncores_label <- NULL
    }
    if (!is.na(max_ncores)) {
      if (ncores > max_ncores) {
        warning("ncores capped at ", max_ncores_label, " = ",
                max_ncores, call. = FALSE)
        ncores <- max_ncores
      }
    }
    
    fit_multichannel <- function(do_parallel) {
      if (do_parallel) {
        nworkers <- min(ncores, m)
        cluster_type <- if (Sys.info()["sysname"] == "Windows") "PSOCK" else "FORK"
        cl <- if (cluster_type == "PSOCK") {
          parallel::makeCluster(nworkers, type = cluster_type, useXDR = FALSE)
        } else {
          parallel::makeCluster(nworkers, type = cluster_type)
        }
        on.exit(parallel::stopCluster(cl), add = TRUE)
        L0glm_dll <- getLoadedDLLs()[["L0glm"]]
        L0glm_dll_path <- if (!is.null(L0glm_dll)) L0glm_dll[["path"]] else NULL
        parallel::clusterCall(cl, L0glm_worker_init, .libPaths(), L0glm_dll_path,
                              worker_nthreads)
      }
      
      make_chunks <- function(vec, k) {
        k <- min(as.integer(k), length(vec))
        cs  <- rep(length(vec) %/% k, k); extra <- length(vec) %% k
        if (extra) cs[1:extra] <- cs[1:extra] + 1
        split(vec, rep(seq_len(k), cs))
      }
      make_jobs <- function(chunks, vars, start.mat) {
        lapply(chunks, function(cols) {
          list(y          = y[, cols, drop = FALSE],
               weights    = weights[, cols, drop = FALSE],
               prior.mean = prior.mean[vars, cols, drop = FALSE],
               dispersion = disp_list[cols],
               start      = if (is.null(start.mat)) NULL else start.mat[vars, cols, drop = FALSE],
               offset     = offset[, cols, drop = FALSE])
        })
      }
      collect_fit <- function(fit, chunked) {
        if (chunked) {
          list(coefficients = do.call(cbind, lapply(fit, `[[`, "coefficients")),
               weights      = do.call(cbind, lapply(fit, `[[`, "weights")),
               dispersion   = unlist(lapply(fit, `[[`, "dispersion"), use.names = FALSE),
               lambdas      = do.call(cbind, lapply(fit, `[[`, "lambdas")),
               z            = do.call(cbind, lapply(fit, `[[`, "z")))
        } else {
          list(coefficients = do.call(cbind, lapply(fit, `[[`, "coefficients")),
               weights      = do.call(cbind, lapply(fit, `[[`, "weights")),
               dispersion   = as.numeric(Matrix(do.call(cbind, lapply(fit, `[[`, "dispersion")),
                                                sparse = FALSE)[1, ]),
               lambdas      = do.call(cbind, lapply(fit, `[[`, "lambdas")),
               z            = do.call(cbind, lapply(fit, `[[`, "z")))
        }
      }
      
      for (iter in seq_len(group.maxit)) {
        if (iter > 1) {
          avgc   <- Matrix::rowMeans(coefficients)
          vars   <- which(avgc != 0)
          always.include <- match(always.include.start[always.include.start %in% vars], vars)
          start  <- matrix(rep(avgc, m), length(avgc), m)
        } else {
          vars <- seq_len(p)
        }
        
        if (!length(vars)) {
          if (verbose) warning("All coefficients zero - use a smaller lambda")
          break
        }
        
        coefficients <- lambdas.mat <- Matrix(0, p, m, sparse = TRUE)
        weights.mat  <- Matrix(0, n, m, sparse = TRUE)
        
        if (do_parallel) {
          chunks <- make_chunks(seq_len(m), length(cl))
          jobs <- make_jobs(chunks, vars, start)
          fit <- parallel::clusterApply(cl, jobs, L0glmfun_chunk,
                                        X             = X[, vars, drop = FALSE],
                                        family        = family,
                                        lambda        = lambda,
                                        gamma         = gamma,
                                        maxit         = maxit,
                                        tol           = tol,
                                        always.include        = always.include,
                                        thresh        = thresh,
                                        solver        = solver,
                                        warm.start    = warm.start,
                                        lower         = lower[vars],
                                        upper         = upper[vars],
                                        method        = method,
                                        nonnegative   = nonnegative,
                                        normalize     = normalize,
                                        control.solver = control.solver,
                                        verbose       = verbose)
        } else {
          fit <- lapply(seq_len(m), function(col) {
            L0glm.fit.internal(X          = X[, vars, drop = FALSE],
                               y          = y[, col, drop = FALSE],
                               weights     = weights[, col],
                               family     = family,
                               lambda     = lambda,
                               gamma      = gamma,
                               normalize   = normalize,
                               offset      = offset[, col],
                               warm.start  = warm.start,
                               prior.mean  = prior.mean[vars, col],
                               dispersion  = disp_list[[col]],          # slice
                               start       = start[vars, col],
                               lower       = lower[vars],
                               upper       = upper[vars],
                               maxit       = maxit,
                               thresh      = thresh,
                               solver      = solver,
                               method      = method,
                               tol         = tol,
                               always.include      = always.include,
                               control.solver = control.solver,
                               verbose     = verbose)
          })
        }
        
        fit <- collect_fit(fit, chunked = do_parallel)
        if (length(vars) == p && identical(vars, seq_len(p))) {
          coefficients <- Matrix(fit$coefficients, sparse = TRUE)
        } else {
          coefficients[vars, ] <- Matrix(fit$coefficients, sparse = TRUE)
        }
        weights.mat          <- Matrix(fit$weights,      sparse = TRUE)
        dispersion           <- fit$dispersion
        disp_list            <- as_list(dispersion)
        if (length(vars) == p && identical(vars, seq_len(p))) {
          lambdas.mat <- Matrix(fit$lambdas, sparse = TRUE)
        } else {
          lambdas.mat[vars, ]  <- Matrix(fit$lambdas, sparse = TRUE)
        }
        z.mat                <- Matrix(fit$z,            sparse = TRUE)

      }
      
      # ----- relaxed - refit selected variables nearly unpenalised -----------
      if (relaxed) {
        disp_list <- as_list(dispersion)   # refresh helper
        avgc   <- Matrix::rowMeans(coefficients)
        vars   <- which(avgc != 0)
        if (length(vars)) {
          always.include <- match(always.include.start[always.include.start %in% vars], vars)
          start  <- matrix(rep(avgc, m), length(avgc), m)
          coefficients <- lambdas.mat <- Matrix(0, p, m, sparse = TRUE)
          weights.mat  <- Matrix(0, n, m, sparse = TRUE)
          
          if (do_parallel) {
            chunks <- make_chunks(seq_len(m), length(cl))
            jobs <- make_jobs(chunks, vars, start)
            fit <- parallel::clusterApply(cl, jobs, L0glmfun_chunk,
                                          X             = X[, vars, drop = FALSE],
                                          family        = family,
                                          lambda        = 1e-30, # using a very small lambda here
                                          gamma         = gamma,
                                          maxit         = maxit,
                                          tol           = tol,
                                          always.include        = always.include,
                                          thresh        = thresh,
                                          solver        = solver,
                                          warm.start    = warm.start,
                                          lower         = lower[vars],
                                          upper         = upper[vars],
                                          method        = method,
                                          nonnegative   = nonnegative,
                                          normalize     = normalize,
                                          control.solver = control.solver,
                                          verbose       = verbose)
          } else {
            fit <- lapply(seq_len(m), function(col) {
              L0glm.fit.internal(X          = X[, vars, drop = FALSE],
                                 y          = y[, col, drop = FALSE],
                                 weights     = weights[, col], 
                                 family     = family,
                                 lambda     = 1e-30,
                                 gamma      = gamma,
                                 normalize   = normalize,
                                 offset      = offset[, col], 
                                 warm.start  = warm.start,
                                 prior.mean  = prior.mean[vars, col], 
                                 dispersion  = disp_list[[col]],
                                 start       = start[vars, col],
                                 lower       = lower[vars],
                                 upper       = upper[vars],
                                 maxit       = maxit,
                                 thresh      = thresh,
                                 solver      = solver,
                                 method      = method,
                                 tol         = tol,
                                 always.include      = always.include,
                                 control.solver = control.solver,
                                 verbose     = verbose)
            })
          }
          
          fit <- collect_fit(fit, chunked = do_parallel)
          if (length(vars) == p && identical(vars, seq_len(p))) {
            coefficients <- Matrix(fit$coefficients, sparse = TRUE)
          } else {
            coefficients[vars, ] <- Matrix(fit$coefficients, sparse = TRUE)
          }
          weights.mat          <- Matrix(fit$weights,      sparse = TRUE)
          dispersion           <- fit$dispersion
          disp_list            <- as_list(dispersion)      # keep the helper in sync
          if (length(vars) == p && identical(vars, seq_len(p))) {
            lambdas.mat <- Matrix(fit$lambdas, sparse = TRUE)
          } else {
            lambdas.mat[vars, ] <- Matrix(fit$lambdas, sparse = TRUE)
          }
          z.mat                <- Matrix(fit$z,            sparse = TRUE)
          
        } else if (verbose) warning("Relaxed fit skipped - all coeffs zero")
      }
      
      # return everything needed by the caller
      list(coefficients = coefficients,
           weights      = weights.mat,
           dispersion   = dispersion,
           lambdas      = lambdas.mat,
           z            = z.mat)
    } # ------- end fit_multichannel() --------------------------------------
    
    res_multi <- fit_multichannel(do_parallel = parallelize && ncores > 1)
    
    coefficients <- res_multi$coefficients
    weights.ret  <- res_multi$weights
    dispersion   <- res_multi$dispersion
    lambdas      <- res_multi$lambdas
    z            <- res_multi$z
    
  } else {
    # =======================================================================
    #  UNIVARIATE RESPONSE  (m == 1) ****************************************
    # =======================================================================
    fit <- L0glm.fit.internal(X          = X,
                              y          = as.vector(y),
                              weights    = weights,
                              family     = family,
                              normalize   = normalize,
                              offset      = offset,
                              dispersion  = dispersion,
                              start       = start,
                              lower       = lower,
                              upper       = upper,
                              lambda      = lambda,
                              warm.start  = warm.start,
                              prior.mean  = prior.mean,
                              gamma       = gamma,
                              thresh      = thresh,
                              solver      = solver,
                              method      = method,
                              maxit       = maxit,
                              tol         = tol,
                              always.include      = always.include,
                              control.solver = control.solver,
                              verbose     = verbose)
    
    coefficients <- fit$coefficients
    dispersion   <- fit$dispersion
    lambdas      <- fit$lambdas
    weights.ret  <- fit$weights
    z            <- fit$z
    
    # --- relaxed - refit selected variables nearly unpenalised -------------
    if (relaxed) { 
      vars <- which(coefficients != 0)
      if (length(vars)) {
        always.include <- match(always.include.start[always.include.start %in% vars], vars)
        fit     <- L0glm.fit.internal(X          = X[, vars, drop = FALSE],
                                      y          = y,
                                      weights    = weights,
                                      family     = family,
                                      normalize   = normalize,
                                      offset      = offset,
                                      dispersion  = dispersion,
                                      start       = coefficients[vars],
                                      lower       = lower[vars],
                                      upper       = upper[vars],
                                      lambda      = 1e-30,
                                      warm.start  = warm.start,
                                      prior.mean  = prior.mean[vars],
                                      gamma       = gamma,
                                      thresh      = thresh,
                                      solver      = solver,
                                      method      = method,
                                      maxit       = maxit,
                                      tol         = tol,
                                      always.include      = always.include,
                                      control.solver = control.solver,
                                      verbose     = verbose)
        coefficients[]     <- 0
        lambdas[]          <- 0
        coefficients[vars] <- fit$coefficients
        lambdas[vars]      <- fit$lambdas
        dispersion         <- fit$dispersion
        weights.ret        <- fit$weights
        z                  <- fit$z
      }
    }
  }
  
  
  # ------------------------------------------------------------------------
  #  Final clean-up  &  return
  # ------------------------------------------------------------------------
  weights.ret[is.infinite(weights.ret)] <- 1e11
  lambdas    [is.infinite(lambdas)]     <- 1e11
  
  ## ensure coefficients are returned as sparse coeff-matrix also for univariate case
  if (m == 1L) {
    coefficients <- Matrix::Matrix(
      coefficients,
      ncol      = 1L,
      sparse    = TRUE,
      dimnames  = list(colnames(X), NULL)   # keep predictor names
    )
  }
  
  list(coefficients = coefficients,
                 weights      = weights.ret,
                 dispersion   = dispersion,
                 init.lambda  = init.lambda,
                 lambdas      = lambdas,
                 family       = family,
                 z            = z,
                 X            = X,
                 y            = y)
}


# ---------------------------------------------------------------------------
#  PARALLEL HELPER (called on the cluster)
# ---------------------------------------------------------------------------
L0glm_detect_physical_cores <- function() {
  physical_cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if ((is.na(physical_cores) || physical_cores < 1L) &&
      requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    physical_cores <- suppressWarnings(RhpcBLASctl::get_num_cores())
  }
  if (is.na(physical_cores) || physical_cores < 1L) {
    NA_integer_
  } else {
    as.integer(physical_cores)
  }
}

L0glm_default_ncores <- function(parallelize = TRUE) {
  if (!isTRUE(parallelize)) return(1L)
  physical_cores <- L0glm_detect_physical_cores()
  if (!is.na(physical_cores)) return(max(1L, physical_cores))
  logical_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
  if (!is.na(logical_cores) && logical_cores >= 1L) {
    return(max(1L, as.integer(round(logical_cores / 2))))
  }
  1L
}

L0glm_thread_env_vars <- function() {
  c("OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS")
}

L0glm_normalize_thread_count <- function(nthreads) {
  nthreads <- suppressWarnings(as.integer(nthreads[1L]))
  if (is.na(nthreads) || nthreads < 1L) {
    stop("nthreads must be NULL or a positive integer", call. = FALSE)
  }
  nthreads
}

L0glm_set_thread_count <- function(nthreads) {
  if (is.null(nthreads)) return(NULL)
  nthreads <- L0glm_normalize_thread_count(nthreads)
  env_vars <- L0glm_thread_env_vars()
  old_env <- Sys.getenv(env_vars, unset = NA_character_)
  blas_threads <- omp_threads <- NA_integer_
  has_blasctl <- requireNamespace("RhpcBLASctl", quietly = TRUE)
  if (has_blasctl) {
    blas_threads <- tryCatch(RhpcBLASctl::blas_get_num_procs(),
                             error = function(e) NA_integer_)
    omp_threads <- tryCatch(RhpcBLASctl::omp_get_max_threads(),
                            error = function(e) NA_integer_)
  }
  Sys.setenv(
    OMP_NUM_THREADS = as.character(nthreads),
    OPENBLAS_NUM_THREADS = as.character(nthreads),
    MKL_NUM_THREADS = as.character(nthreads),
    VECLIB_MAXIMUM_THREADS = as.character(nthreads)
  )
  if (has_blasctl) {
    try(RhpcBLASctl::blas_set_num_threads(nthreads), silent = TRUE)
    try(RhpcBLASctl::omp_set_num_threads(nthreads), silent = TRUE)
  }
  list(env = old_env,
       blas_threads = blas_threads,
       omp_threads = omp_threads,
       has_blasctl = has_blasctl)
}

L0glm_restore_thread_count <- function(state) {
  if (is.null(state)) return(invisible(NULL))
  for (nm in names(state$env)) {
    if (is.na(state$env[[nm]])) {
      Sys.unsetenv(nm)
    } else {
      do.call(Sys.setenv, structure(list(state$env[[nm]]), names = nm))
    }
  }
  if (isTRUE(state$has_blasctl) &&
      requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    if (!is.na(state$blas_threads)) {
      try(RhpcBLASctl::blas_set_num_threads(state$blas_threads), silent = TRUE)
    }
    if (!is.na(state$omp_threads)) {
      try(RhpcBLASctl::omp_set_num_threads(state$omp_threads), silent = TRUE)
    }
  }
  invisible(NULL)
}

L0glm_worker_init <- function(lib_paths, dll_path, nthreads = 1L) {
  .libPaths(lib_paths)
  nthreads <- suppressWarnings(as.integer(nthreads[1L]))
  if (is.na(nthreads) || nthreads < 1L) nthreads <- 1L
  Sys.setenv(
    OMP_NUM_THREADS = as.character(nthreads),
    OPENBLAS_NUM_THREADS = as.character(nthreads),
    MKL_NUM_THREADS = as.character(nthreads),
    VECLIB_MAXIMUM_THREADS = as.character(nthreads)
  )
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    try(RhpcBLASctl::blas_set_num_threads(nthreads), silent = TRUE)
    try(RhpcBLASctl::omp_set_num_threads(nthreads), silent = TRUE)
  }
  if (!is.null(dll_path) && nzchar(dll_path) &&
      file.exists(dll_path) &&
      !("L0glm" %in% names(getLoadedDLLs()))) {
    dyn.load(dll_path)
  }
  if (!("L0glm" %in% loadedNamespaces())) {
    requireNamespace("L0glm", quietly = TRUE)
  }
  NULL
}

L0glmfun_chunk <- function(job, ...) {
  fit <- L0glmfun(cols       = seq_len(ncol(job$y)),
                  y          = job$y,
                  weights    = job$weights,
                  prior.mean = job$prior.mean,
                  dispersion = job$dispersion,
                  start      = job$start,
                  offset     = job$offset,
                  ...)

  list(coefficients = do.call(cbind, lapply(fit, `[[`, "coefficients")),
       weights      = do.call(cbind, lapply(fit, `[[`, "weights")),
       dispersion   = vapply(fit, `[[`, numeric(1), "dispersion"),
       lambdas      = do.call(cbind, lapply(fit, `[[`, "lambdas")),
       z            = do.call(cbind, lapply(fit, `[[`, "z")))
}

L0glmfun <- function(cols, lambda, gamma, thresh, solver, method,
                     nonnegative, normalize,
                     X, y, family,
                     warm.start = 0.9,
                     prior.mean = NULL, weights = NULL, dispersion = NULL,
                     start = NULL, lower, upper,
                     verbose = FALSE, maxit, tol, offset, always.include = 0,
                     control.solver = do.call("control.solver.gen",
                                              list(solver = solver, method = method))) {
  
  lapply(cols, function(col) {
    disp_col <- if (is.list(dispersion)) dispersion[[col]] else dispersion
    L0glm.fit.internal(
      X            = X,
      y            = as.vector(y[, col]),
      weights      = as.vector(weights[, col]),
      family       = family,
      lambda       = lambda,
      gamma        = gamma,
      dispersion   = disp_col,          # <- per channel list element (may be NULL)
      start        = as.vector(start[, col]),
      maxit        = maxit,
      tol          = tol,
      always.include       = always.include,
      lower        = lower,
      offset       = as.vector(offset[, col]),
      upper        = upper,
      prior.mean   = prior.mean[, col],
      thresh       = thresh,
      solver       = solver,
      normalize    = normalize,
      method       = method,
      verbose      = verbose,
      control.solver = control.solver)
  })
}
  
