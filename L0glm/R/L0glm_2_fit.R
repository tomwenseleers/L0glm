#' @name L0glm.fit
#' @aliases L0glm.fit
#' @title Fit generalized linear model with \eqn{L_0} or grouped row-\eqn{L_0}
#'   regularization subject to optional nonnegativity constraints
#'
#' @description Fits a generalized linear model with an \eqn{L_0} penalty or,
#' for multichannel responses, a grouped row-\eqn{L_0} penalty, equivalent
#' to an \eqn{L_{2,0}} penalty, to approximate Best Subset Selection (for
#' univariate data) or Best Subset of Groups Selection (for multivariate data)
#' subject to optional nonnegativity constraints on model coefficients.
#' The objective functions that can be approximated are :
#' \itemize{
#'
#' \item for univariate data, using an \eqn{L_0} norm penalty to achieve Best Subset Selection:
#'
#' \deqn{-2\,LL(\beta) + \lambda \|\beta\|_0}
#'
#' where \eqn{LL} is the log-likelihood of the GLM model and
#' \eqn{\|\beta\|_0} is the number of non-zero elements in
#' coefficient vector beta, i.e. the \eqn{L_0} norm of beta.
#'
#' \item for multivariate data, using independent elementwise \eqn{L_0}
#' penalized fits on each outcome:
#'
#' \deqn{-2\,LL(B) + \lambda \sum_{j=1}^{p}\sum_{k=1}^{m} I(B_{jk} \neq 0)}
#'
#' \item for multivariate data, using grouped variable selection in multitask
#' learning settings to achieve Best Subset of Groups Selection:
#'
#' \deqn{-2\,LL(B) + \lambda \sum_{j=1}^{p} I(\|B_{j\cdot}\|_2 > 0)}
#'
#' }
#'
#' Here \eqn{B} is the \eqn{p \times m} multichannel coefficient matrix and the
#' grouped penalty counts active rows, i.e. predictors with a non-zero
#' coefficient in at least one channel. Any positive within-row norm gives the
#' same row support, so this is equivalently group-\eqn{L_0}, row-\eqn{L_0}, or
#' \eqn{L_{2,0}} regularization.
#'
#' In all three cases, \code{nonnegative = TRUE} restricts the coefficients to
#' be nonnegative, including any intercept, except that binomial and
#' quasibinomial intercepts remain unrestricted:
#'
#' \deqn{\beta_j \geq 0 \quad\hbox{or}\quad B_{jk} \geq 0.}
#' With \code{nonnegative = FALSE}, coefficients may have either sign.
#'
#' These objectives are approximated via an iterative adaptive ridge regression
#' procedure, see Details. The univariate and independent multichannel
#' elementwise \eqn{L_0} penalties can be used to approximate Best Subset
#' Selection, while for multivariate (multitask/multichannel) data, the
#' group-\eqn{L_0}/row-\eqn{L_0}/\eqn{L_{2,0}} penalty can be used to approximate
#' Best Subset of Groups Selection (via option \code{grouped=TRUE}). Fits on
#' multichannel data can be parallelized. Option \code{relaxed} allows one to
#' finish with a (nearly) unpenalized fit on the \code{L0glm} selected features.
#'
#' The level of regularisation can be manually chosen (if a numeric value is
#' passed to \code{lambda}) or preset to approximate an information criterion.
#' Criterion names are case-insensitive. See \link{information_criteria} for
#' the penalty formulas, guidance on choosing a criterion and original references.
#'
#' The \code{L0glm.fit} function has various advanced features, including
#' support for both dense and sparse design matrices \emph{\bold{X}}, high
#' dimensional data (where \emph{p} >> \emph{n}), single or multichannel
#' outcomes \emph{\bold{y}} (with either independent or grouped variable
#' selection), support for any of the standard R GLM families & link functions,
#' prior weights (e.g. total nr. of trials for binomial proportions), model
#' offsets, starting coefficient values or dispersion (if prior estimate is
#' available, otherwise optimised), optional nonnegativity constraints,
#' unpenalized terms (e.g. intercept) and shrinkage towards given prior means
#' (via targeted nonzero-centered ridge regression updates).
#'
#' A choice of several fast C++ coded least squares solvers can be used, which
#' means this function scales well to solve fairly large tall or wide problems
#' with millions of observations or variables. As solvers, we currently
#' implemented the \code{Eigen} sparse or dense Cholesky, QR, PCG & LSCG plus 
#' BDCSVD and Jacobi SVD solvers, as well as two solvers that natively support 
#' nonnegativity constraints: \code{bigGlm} from \pkg{glmnet}, and \code{osqp}, a fast
#' quadratic programming solver.
#'
#' @param X Design matrix with predictor variables of dimension \code{n} x
#'   \code{p}, with \code{n} corresponding to the number of observations and
#'   \code{p} being the number of predictor columns. Use \code{intercept = TRUE}
#'   to add an intercept, or include a column of ones in \code{X} and set
#'   \code{intercept = FALSE}.
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
#' @param dispersion \code{NULL} (default) automatically initializes and
#'   estimates dispersion. A positive number fixes it at that value. A negative
#'   number initializes it at its absolute value, then estimates it; the supplied
#'   magnitude is used directly, without the automatic permissive scaling.
#'   For multichannel data, supply one number or a vector with one number per
#'   channel; fixed and estimated channels may be mixed. Zero and non-finite
#'   values are invalid. Poisson, binomial and negative binomial with known
#'   \eqn{\theta} have dispersion fixed at 1; other fixed-scale families respect
#'   \code{family$dispersion}. See Details for initialization and updating.
#' @param start Optional scalar, vector of coefficient starting values of length
#'   \code{p}, or matrix of starting values with dimension \code{p} x
#'   \code{m}, with \code{p} the number of design columns, including any added
#'   intercept, and \code{m} the number of response columns. A scalar is repeated
#'   for every coefficient in every response channel.
#'   Numeric scalar \code{start = 1} initializes all coefficients to one in
#'   normalized working coordinates, including the intercept. With a free
#'   intercept this refers to the centered design. All other scalars and all
#'   explicit vectors/matrices use original input coefficient units. To request
#'   original-coordinate ones, pass a vector or a \code{p x m} matrix; use a
#'   matrix to distinguish this from the scalar shorthand in a one-term model.
#'   Starting values are used directly, subject to the coefficient
#'   constraints. Exact zeros exclude those coefficients throughout selection,
#'   grouped passes and relaxed refits, separately in each channel. Rows zero
#'   across every channel are screened before fitting. Intercepts and
#'   \code{always.include} terms remain eligible even with zero starts.
#'   Only \code{NULL} requests automatic initialization; an all-zero start
#'   fits just the mandatory terms, if any. For binomial
#'   and quasibinomial models with an intercept, automatic slope starts account
#'   for its cross-products with the predictors.
#'   These values determine the
#'   level of penalisation of the different coefficients during the first IRLS
#'   ridge regression update. The ridge penalty used equals
#'   \eqn{(\lambda/4) \phi / |\beta|^\gamma}. After each
#'   selection update, coefficients below the relative cutoff or the absolute
#'   numerical floor described in \code{thresh} are excluded from later solves.
#'   Set to \code{NULL} by default, in which case it
#'   is auto-initialised.
#' @param offset Optional vector of offset values of length \code{n}, or a
#'   matrix of offset values with dimensions \code{c(n, m)}. The offset is a
#'   known quantity that is added to the linear predictor of the model.
#' @param lambda Either a scalar with the \eqn{L_0} penalty for univariate or
#'   independent multichannel fits (when \code{grouped=FALSE}) or the grouped
#'   row-\eqn{L_0}/\eqn{L_{2,0}} penalty (when \code{grouped=TRUE}), or a
#'   case-insensitive criterion name: \code{"aic"}, \code{"bic"},
#'   \code{"gic"}, \code{"ebic"}, or \code{"mbic"}.
#'   For example, \code{"GIC"} and \code{"gic"} are equivalent.
#'   The preset is a per-coefficient penalty on the \eqn{-2\log L} scale:
#'   AIC uses 2, BIC uses \eqn{\log n}, and GIC uses \eqn{\log p\,\log\log n}.
#'   Here \eqn{p} counts the original selectable predictors, excluding intercepts
#'   and \code{always.include} terms, before screening by explicit starts.
#'   A pool with no selectable predictors uses \eqn{p=1} in logarithms.
#'   See \link{information_criteria} for
#'   all formulas, intended uses and references. Set \code{lambda = 0} to
#'   disable penalization. Presets approximate the chosen criterion; they do
#'   not guarantee its global minimizer.
#' @param intercept Logical: Whether to include an intercept term in the model.
#'   If `TRUE`, a column of ones will be added to the design matrix `X`, and
#'   this term will **not** be penalized. If `FALSE`, no intercept column will
#'   be added, and the user must ensure that `X` already contains an intercept
#'   if needed. Defaults to `TRUE` for binomial models and `FALSE` otherwise
#'   in `L0glm.fit()`, and to `FALSE` in `L0nnpois()`. With
#'   `nonnegative = TRUE`, the intercept is also constrained to be nonnegative,
#'   except for binomial and quasibinomial models, where it may have either sign.
#' @param always.include A vector of coefficient indices that should always be
#'   estimated & included in the model. These terms are unpenalized and remain
#'   eligible when their current estimate or explicit starting value is zero.
#'   Coefficients with
#'   a name that includes \code{"Intercept"} in their name are always included. 
#'   Use \code{NULL} or 0 to request no additional mandatory terms;
#'   recognized intercepts remain mandatory. If `intercept = TRUE`, the first column (intercept) will
#'   always be estimated and included in the model.
#' @param nonnegative A single logical value. If \code{TRUE}, all coefficients,
#'   including any intercept and \code{always.include} terms, must be
#'   nonnegative. For binomial and quasibinomial models only, intercept columns
#'   (columns of ones, whether supplied in \code{X} or added by
#'   \code{intercept = TRUE}) remain unrestricted. If \code{FALSE}, all
#'   coefficients may have either sign. Defaults to \code{FALSE} in
#'   \code{L0glm.fit()}. The specialized \code{L0nnpois()} function always
#'   enforces nonnegativity and has no \code{nonnegative} argument.
#' @param normalize Preprocessing convention: \code{"auto"} (default) uses
#'   \code{"rms"} for signed fits or a free intercept, and \code{"none"} for
#'   fully nonnegative fits. Other choices are \code{FALSE} or \code{"none"},
#'   \code{TRUE} or \code{"l2"} for unit L2 norm, \code{"rms"} for unit
#'   root-mean-square (\eqn{\sum_i X_{ij}^2/n=1}), and \code{"sd"} for unit
#'   sample standard deviation. With a free intercept, column means are
#'   accounted for and the norms use centered columns. The original \code{X}
#'   and \code{y} are unchanged. Dense working columns are centered once;
#'   sparse working columns are shifted only when this adds no stored entries,
#'   with the other means handled by the intercept. A constrained intercept
#'   prevents this centering, including in \code{L0nnpois()}. Mandatory terms
#'   and zero-scale columns use scale one. Penalties, starts, ridge centers and
#'   bounds are transformed consistently. Returned coefficients and explicit
#'   vector/matrix starts retain their input units; numeric scalar \code{start=1}
#'   is the normalized-coordinate shorthand described above. Non-\code{"none"} modes also use
#'   these scales for automatic initialization, thresholding and convergence.
#'   With \code{gamma = 2}, automatic signed fits are therefore invariant to
#'   changes of predictor units, up to numerical precision. Sparse designs
#'   remain sparse. The training centers and scales are returned in
#'   \code{preprocessing}; predictions use the original columns directly.
#' @param solver Solver used for the ridge regression step in the IRLS algorithm
#'   (see details). It should be one of \code{"osqp"} (osqp quadratic
#'   programming solver), \code{"bigGlm"} (using \code{\link[glmnet]{bigGlm}}
#'   coordinate descent; legacy alias \code{"glmnet"}) or \code{"eigen"} (one
#'   of the \code{RcppEigen Eigen} solvers).
#'   \code{"osqp"} and \code{"bigGlm"} support nonnegativity constraints
#'   and scale well to larger problems, with especially \code{"bigGlm"} being
#'   suitable for very large problems. The \code{"eigen"} solvers apply
#'   nonnegativity constraints via clipping, which in combination with the ridge
#'   penalties in practice work almost as well, and some of those solvers are
#'   often among the fastest (especially \code{method} \code{1L}, Cholesky
#'   LLT, for PSF-like sparse deconvolution problems).
#' @param method When \code{solver="eigen"}, method used for the eigen solver
#'   (see details). It should be an integer from 1 to 6 :
#'
#'   \code{1L} = Cholesky LLT solver, with sparse AMD ordering by default
#'
#'   \code{2L} = Cholesky LDLT solver, with sparse AMD ordering by default
#'
#'   \code{3L} = Preconditioned Conjugate Gradient (PCG) with incomplete
#'   Cholesky preconditioner
#'
#'   \code{4L} = Least Squares Conjugate Gradient (LSCG) with least squares
#'   diagonal preconditioner
#'
#'   \code{5L} = Column-pivoted QR solver. Sparse QR uses numerical column
#'   pivoting and COLAMD fill-reducing ordering by default.
#'
#'   \code{6L} = BDCSVD solver for dense systems, or Jacobi SVD for sparse
#'   systems.
#'
#'   Defaults to \code{1L} (Cholesky LLT), which is often the fastest.
#'
#' @param openmp Logical; defaults to \code{TRUE}. For nonnegative
#'   identity-link Poisson fits with \code{solver="eigen"} and
#'   \code{method=1L}, use the compiled Rcpp OpenMP \code{\link{L0nnpois}} backend.
#'   Other families, links and solvers use the general R fitting path.
#'   Set \code{FALSE} to use that path explicitly. Solver controls requesting
#'   another factorization or \code{rowaug=FALSE} also retain the R path.
#'   This option selects the backend; \code{parallelize} controls whether
#'   its outcome channels run concurrently.
#' @param parallelize A logical value indicating whether in the multivariate
#'   case, the columns of the response matrix \code{y} should be fitted in
#'   parallel. If \code{TRUE}, each column of \code{y} is fitted concurrently,
#'   potentially speeding up computation. If \code{FALSE}, the columns are
#'   fitted sequentially.
#' @param ncores Optional integer specifying the number of cores to use when
#'   \code{parallelize=TRUE}. If \code{parallelize=TRUE}, defaults to half the
#'   detected number of physical cores when that can be determined, falling back
#'   to half the number of logical cores. The requested value is capped at the
#'   detected number of physical cores when that can be determined.
#'   With the \code{L0nnpois} backend, \code{ncores} instead controls OpenMP
#'   threads sharing the design matrix and is capped at the logical core count.
#' @param nthreads Optional integer specifying the number of BLAS threads
#'   to use per fit process. The default \code{NULL} leaves the main R session's
#'   thread settings unchanged and uses one BLAS thread per parallel
#'   worker. Setting \code{nthreads} can be useful to benchmark serial fits or
#'   tune per-worker threading on systems with a multithreaded BLAS.
#'   The \code{L0nnpois} backend uses \code{ncores}; it ignores \code{nthreads}.
#' @param maxit Integer specifying the maximum number of IRLS iterations to
#'   perform per channel and fitting pass (default 25). Fits stop earlier
#'   when the convergence criterion is met; inspect \code{convergence} to
#'   distinguish this from reaching the iteration cap.
#' @param grouped A logical value indicating whether a grouped fit should be
#'   performed in the multivariate case or not. If \code{TRUE}, the algorithm
#'   will return a group penalized fit, otherwise independent L0glm fits will be
#'   performed on each response. These procedures approximate grouped
#'   row-\eqn{L_0}/\eqn{L_{2,0}} and independent elementwise \eqn{L_0} penalties,
#'   and can be used to do Best Subset of Groups Selection or plain Best Subset
#'   Selection.
#' @param group.maxit Integer  specifying the maximum number of iterations to
#'   perform for the group penalized fitting procedure when \code{grouped =
#'   TRUE}. In each iteration, the columns of the response matrix \code{y} will
#'   be refitted on the active row support using shared row-level starting
#'   values derived from the current coefficients, which then imposes a group
#'   penalty.
#' @param relaxed Logical indicating if the algorithm should conduct a nearly
#'   unpenalized fit on the L0glm selected variables (using
#'   \code{lambda=1E-30}). Hard thresholding is disabled in this final refit,
#'   while adaptive ridge weights remain unsmoothed. This can
#'   reduce coefficient bias and increase sensitivity, especially for the multivariate case
#'   with grouped variable selection. Defaults to \code{TRUE} when
#'   \code{grouped=TRUE} & multichannel data are passed in \code{y}.
#' @param prior.mean An optional numeric vector of length \code{p} specifying
#'   the center towards which coefficients should be shrunk in each IRLS ridge
#'   regression iteration. This is done via nonzero centered ridge regression.
#'   This can be useful if there is good prior information available that some
#'   coefficients should have a particular nonnegative value. The default
#'   \code{prior.mean = NULL} will shrink to zero using a default zero centered
#'   ridge regression.
#' @param gamma Scalar with the exponent to use when calculating the ridge
#'   penalty \code{(lambda/4)*dispersion/abs(beta)^gamma}.
#'   Defaults to \code{gamma = 2}, as in Frommlet & Nuel 2016.
#' @param thresh Finite nonnegative scalar, default \code{1e-5}, controlling
#'   hard thresholding during selection. After each update, coefficients with
#'   \code{abs(beta)*scale < max(abs(beta)*scale)*thresh} or
#'   \code{abs(beta)*scale < 1e-20} are set to zero. With
#'   \code{normalize = "none"}, \code{scale = 1}. The maximum is computed
#'   separately for each response channel, excluding unpenalized mandatory
#'   terms; these terms are never hard thresholded.
#'   Exact zeros remain excluded within that IRLS run. Setting \code{thresh=0}
#'   disables only the relative cutoff; the absolute numerical floor remains.
#'   The final relaxed refit applies neither cutoff; nonnegativity clipping
#'   still applies when requested. Adaptive ridge weights are unsmoothed,
#'   \code{1/abs(beta)^gamma}, independently of \code{thresh}.
#' @param tol Convergence tolerance (default \code{1e-10}). Iterations stop when the relative
#'   change in the \eqn{L_2} norm of the coefficients drops below this
#'   threshold. With preprocessing, the test uses scaled selectable slopes
#'   and, when mandatory terms exist, changes in predictions, so that the
#'   predictor or response origin cannot hide slope changes.
#'   When dispersion is estimated, selection also waits until its
#'   update has stabilized as described in Details, subject to \code{maxit}.
#' @param control.solver List of parameters controlling the optimization
#'   solvers. Auto-initialised with \code{control.solver.gen} when \code{NULL}.
#'   See \code{?control.solver.gen} for more details on the available parameters
#'   for the different solvers. For sparse Eigen Cholesky LLT/LDLT,
#'   \code{control.solver$NaturalOrdering} selects Eigen natural ordering
#'   (\code{TRUE}) or AMD ordering (\code{FALSE}, default). Sparse QR
#'   defaults to Eigen \code{COLAMDOrdering} (\code{ordering=1L});
#'   \code{ordering=0L} uses AMD and \code{ordering=2L} uses NaturalOrdering.
#' @param verbose Logical indicating if the algorithm progression should be
#'   printed to the console.
#'
#' @details L0glm approximates the \eqn{L_0}, independent multichannel
#' elementwise \eqn{L_0}, or grouped row-\eqn{L_0}/\eqn{L_{2,0}} penalized GLM
#' objectives via an iterative adaptive ridge regression procedure.
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
#' centered) adaptive ridge regularised regression with optional
#' nonnegativity constraints on the coefficients. In each iteration, this uses
#' the following ridge penalties for the penalized terms (excluding the
#' intercept, if included):
#' \deqn{\Lambda = (1/4)\lambda\sigma^2 w}
#' with adaptive penalty weights
#' \deqn{w = 1 / |\beta|^\gamma}
#' where \eqn{\beta} are the coefficient estimates
#' in each iteration. During selection, the relative cutoff and absolute
#' numerical floor specified by \code{thresh} remove coefficients from later solves.
#' The final relaxed refit does not apply this hard cutoff. As in Frommlet
#' & Nuel (2016), \eqn{\gamma} is set to a default value of 2. Nonnegativity
#' constraints are natively supported if the \code{bigGlm} or
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
#' \code{dispersion = NULL} or a negative initial value is supplied).
#'
#' Automatic initialization uses one weak ridge fit to the initial working
#' response, with columns scaled to unit weighted norm and ridge constant
#' \code{1e-6}. Four reproducible random-sign probes share this factorization
#' to estimate the residual noise degrees of freedom,
#' \eqn{\operatorname{tr}\{(I-H)^T(I-H)\}}, where \eqn{H} is the pilot hat matrix.
#' The caller's random-number state is preserved. The pilot residual scale is
#' divided by \code{max(1, preset.lambda("gic", p, n))} for a permissive start
#' independent of the fitted lambda; its support is never
#' imposed on the main fit. A supplied negative value bypasses this pilot.
#' With preprocessing enabled, mandatory terms are profiled out and column
#' information is conditional on them, preserving changes of predictor origin
#' when an unrestricted intercept is present.
#' If this weak pilot estimates a variance below \code{1e-6} times the
#' nuisance-adjusted working-response variance, one retry uses pilot ridge
#' \code{0.1}. This guards against a nearly interpolating pilot freezing a
#' dense fit at an extremely small dispersion. It is an initialization
#' safeguard, not a statistical lower bound on the true dispersion.
#'
#' Starting at iteration 3, when the number \eqn{k} of active coefficients is
#' below \eqn{0.7N}, dispersion is updated toward
#' \deqn{\widetilde\phi =
#'   \frac{\sum_{i:w_i>0}w_i(y_i-\mu_i)^2/V(\mu_i)}{N-k}.}
#' Here \eqn{N} counts positive-weight observations, \eqn{k} includes the
#' intercept, and \eqn{V} is \code{family$variance}, evaluated at the newly
#' fitted mean. Each update is bounded to \eqn{[\phi/2,2\phi]}. After three
#' consecutive updates with unchanged support and a target within 1\% of the
#' current dispersion, it is held fixed until the support changes again.
#' A later support loss reactivates these bounded updates in the same IRLS loop.
#' Each subsequent grouped selection pass may update it again on the new support.
#' A relaxed refit holds dispersion fixed during its coefficient solves and
#' reports the final Pearson residual scale for estimated channels; user-supplied
#' positive dispersions remain fixed throughout.
#'
#' These are inexpensive moment updates, not exact dispersion maximum likelihood
#' for every family. The active coefficient count approximates model degrees of
#' freedom and does not correct for selection. The pilot, permissive start and
#' freezing safeguards are heuristics; they do not guarantee a global optimum or
#' an unbiased noise estimate, especially in nearly interpolating models. If the
#' active model remains too large for the update, the starting scale is retained.
#' The pilot normally requires one factorization (two with the rare retry);
#' subsequent scale updates require only
#' \eqn{O(N)} arithmetic, without nested fitting or further factorizations.
#' Supplying a dispersion value avoids the pilot factorization, which can be
#' appreciable for very large dense designs even when the main solver is iterative.
#'
#' In contrast to Frommlet & Nuel (2016) we use Fisher scoring instead of
#' Newton-Raphson iterations (i.e. using the expected Hessian instead of the
#' observed one). If option \code{relaxed} is set to \code{TRUE} the algorithm
#' will refit the \code{L0glm} selected variables using a nearly unpenalized
#' model (with \code{lambda=1E-30}) without hard thresholding. This can reduce bias and
#' increase sensitivity. If an appropriate level of regularisation lambda is chosen,
#' the L0 norm penalized solution approximates Best Subset selection. To achieve
#' variable selection consistency, the \code{"bic"} \code{lambda} preset is
#' recommended for low dimensional settings, whilst \code{"gic"},
#' \code{"ebic"} or \code{"mbic"} may work better for high
#' dimensional settings, subject to the conditions discussed in
#' \link{information_criteria}. Setting it to
#' \code{"aic"} is recommended if predictive performance is of key interest, as
#' the AIC criterion asymptotically approximates the leave-one-out
#' cross-validation error.
#'
#' \strong{Algorithm for multivariate/multichannel data}
#'
#' For the multivariate ("multitask learning") case there are two options. First,
#' if \code{grouped=FALSE}, an independent \code{L0glm} model is fit on each of
#' the responses, which corresponds to an elementwise \eqn{L_0} penalty on the
#' coefficient matrix.
#' Second, if \code{grouped=TRUE}, grouped variable selection will be performed
#' via a row-\eqn{L_0} penalty, conventionally denoted \eqn{L_{2,0}} in the
#' multitask/MMV literature. This is achieved by taking the independent
#' elementwise \eqn{L_0} solution, and then, for all variables that have nonzero
#' coefficients in at least one of the channels, repeating the \eqn{L_0}
#' penalized fits but using a shared row-level start derived from the
#' coefficients across all channels (response variables) so as to achieve a
#' group penalty and perform grouped feature selection. Active rows are detected
#' by a row-wise nonzero count, equivalently by any positive within-row norm.
#' With signed coefficients, the shared start uses mean absolute magnitude
#' with the sign of the row mean. If that mean is exactly zero, the sign of the
#' first largest-magnitude entry is used. Thus opposite signs across channels
#' cannot reduce the starting magnitude by cancellation. Nonnegative fits retain
#' arithmetic averaging.
#' This is then further iterated until convergence. Finally, as for single
#' channel data, if \code{relaxed=TRUE}, the algorithm can finish with a final
#' step whereby the selected variables are refit using a nearly unpenalized model
#' (with \code{lambda=1E-10}). Compared to other group penalized regression
#' methods, such as the group LASSO or multitask lasso, which use an
#' \eqn{L_{2,1}} penalty, this approach produces a sparser solution that better
#' approximates the Best Subset of Groups solution if lambda is adequately
#' chosen.
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
#' \item{\code{p}}{original number of selectable predictors.}
#' \item{\code{always.include}}{indices of unpenalized mandatory coefficients,
#' including recognized intercepts, in the returned coefficient matrix.}
#' \item{\code{preprocessing}}{training column centers and scales, resolved
#' normalization mode, centering eligibility, and original selectable count.
#' Coefficients are already on the original input scale.}
#' \item{\code{convergence}}{one row per channel and selection or relaxed pass,
#' containing its stage, pass, channel, iteration count, convergence flag,
#' number of dispersion updates and number of reactivations after a freeze.
#' Convergence refers to the iterative stopping test, not a certificate of
#' constrained or global optimality.}
#'
#' @examples
#' x <- seq(-1, 1, length.out = 30)
#' X <- matrix(x, ncol = 1)
#'
#' # Nonnegativity includes the Gaussian intercept.
#' constrained <- L0glm.fit(X, -2 + 3 * x, family = gaussian(),
#'   intercept = TRUE, nonnegative = TRUE, lambda = 1e-8,
#'   dispersion = 1, parallelize = FALSE)
#' as.numeric(constrained$coefficients)
#'
#' # Estimate Gaussian variance from an automatic or user-supplied start.
#' set.seed(1)
#' y_gaussian <- 2 * x + rnorm(length(x), sd = 0.5)
#' automatic <- L0glm.fit(X, y_gaussian, family = gaussian(),
#'   lambda = "gic", dispersion = NULL, maxit = 100, parallelize = FALSE)
#' initialized <- L0glm.fit(X, y_gaussian, family = gaussian(),
#'   lambda = "gic", dispersion = -0.25, maxit = 100, parallelize = FALSE)
#' # dispersion = 0.25 instead keeps the variance fixed at 0.25.
#'
#' # A binomial intercept may remain negative with a nonnegative slope.
#' y <- round(1000 * plogis(-2 + x)) / 1000
#' logistic <- L0glm.fit(X, y, family = binomial(),
#'   weights = rep(1000, length(y)), nonnegative = TRUE,
#'   lambda = 1e-8, parallelize = FALSE)
#' as.numeric(logistic$coefficients)
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
#' @param control.family Family-specific IRLS controls. A named list with \code{theta.min=1e-6}, \code{theta.max=1e7}, \code{theta.tol=1e-6}. See the dedicated family help pages.
#' @section Additional families:
#' Use \code{mgcv::nb()} for negative binomial fitting with fixed or estimated
#' shape. See \link{L0glm_negative_binomial} for mathematics and an example.
#' These family-specific loops are serial. Ordinary GLMs retain their existing paths.
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
                      normalize     = "auto",
                      solver        = "eigen",
                      method        = 1L,
                      parallelize   = TRUE,
                      ncores        = L0glm_default_ncores(parallelize),
                      nthreads      = NULL,
                      maxit         = 25,
                      grouped       = TRUE,
                      group.maxit   = 6,
                      relaxed       = grouped & NCOL(y)>1,
                      prior.mean    = NULL,
                      gamma         = 2,
                      thresh         = 1e-5,
                      tol           = 1e-10,
                      control.solver = do.call("control.solver.gen", list(solver = solver, method = method)),
                      verbose       = FALSE,
                      openmp        = TRUE,
                      control.family = list()) {
  validate_thresh(thresh)
  L0glm_normalization(normalize)
  if (!solver %in% c("eigen", "bigGlm", "bigglm", "glmnet", "osqp"))
    stop("Invalid solver choice; use eigen, bigGlm, or osqp.", call. = FALSE)
  
  if (!is.logical(nonnegative) || length(nonnegative) != 1L || is.na(nonnegative)) {
    stop("nonnegative must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(openmp) || length(openmp) != 1L || is.na(openmp)) {
    stop("openmp must be TRUE or FALSE.", call. = FALSE)
  }
  if (!inherits(family, "family")) family <- family()

  kind <- L0glm_ext_kind(family)
  if (!is.null(kind)) return(L0glm_ext_fit(
    X, y, weights, family, dispersion, start, offset, lambda,
    if (missing(intercept)) FALSE else intercept,
    always.include, nonnegative, normalize, solver, method, nthreads, maxit,
    grouped, if (missing(relaxed)) grouped && NCOL(y) > 1L else relaxed,
    prior.mean, gamma, thresh, tol, control.solver, control.family))

  # Dispatch before R-specific coercion, initialization or worker creation.
  # Forward the public arguments unchanged so defaults and warm starts are
  # applied exactly once by L0nnpois.
  if (openmp && nonnegative && identical(family$family, "poisson") &&
      identical(family$link, "identity") && identical(solver, "eigen") &&
      isTRUE(method == 1L) &&
      (is.null(control.solver$method) || identical(control.solver$method, "chol")) &&
      (is.null(control.solver$pivot) || isTRUE(control.solver$pivot == 0L)) &&
      (is.null(control.solver$rowaug) || isTRUE(control.solver$rowaug))) {
    return(L0nnpois(
      X = X, y = y, weights = weights, family = family,
      dispersion = dispersion, start = start, offset = offset,
      lambda = lambda, intercept = intercept, always.include = always.include,
      normalize = normalize, solver = solver, method = method,
      parallelize = parallelize, ncores = ncores, nthreads = nthreads,
      maxit = maxit, grouped = grouped, group.maxit = group.maxit,
      relaxed = relaxed, prior.mean = prior.mean,
      gamma = gamma, thresh = thresh, tol = tol,
      control.solver = control.solver, verbose = verbose
    ))
  }

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
  
  # Add the intercept before deriving the solver's coefficient bounds --------
  if (intercept) {
    X <- cbind(Intercept = 1, X)
    intercept.col <- 1  # remember its position
    always.include <- c(1, always.include+1)   # always include intercept (indexing for always.include is for non-intercept coefs)
    p <- ncol(X)               # update p after intercept addition
    # TODO also check dimensions of prior.mean, weights, start & offset
    # later in L0glm_1_fit_internal.R when colnames(X)[[1]]=="Intercept"
  } else { intercept.col <- NULL }

  lower <- rep(if (nonnegative) 0 else -.Machine$double.xmax, p)
  upper <- rep(.Machine$double.xmax, p)
  always.include <- L0glm_mandatory(X, always.include)
  p_selectable <- p - length(always.include)
  if (nonnegative && family$family %in% c("binomial", "quasibinomial")) {
    # Recognize intercepts supplied explicitly in X as well as the added column.
    lower[L0glm_intercepts(X)] <- -.Machine$double.xmax
  }
  preprocessing <- L0glm_preprocessing(X, normalize, lower, always.include)
  
  # ----- preset lambda ------------------------------------------------------
  lambda <- L0glm_resolve_lambda(lambda, p = max(1L, p_selectable), n = nrow(X))
  
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
      offset <- if (is.null(dim(offset))) Matrix::Matrix(offset, n, 1, sparse = TRUE) else
        Matrix::Matrix(offset, sparse = TRUE)
    
    if (m > 1 && ncol(offset) == 1L)
      offset <- Matrix::Matrix(rep(offset, m), n, m, sparse = TRUE)
  }
  
  # 3  Starting coefficient values --------------------------------------------
  if (!is.null(start) && length(start) != 1L) {
    start <- as.matrix(start)
    if (m > 1 && ncol(start) == 1L)
      start <- matrix(rep(start, m), nrow = nrow(start), ncol = m)
    if (!identical(dim(start), c(p, m)))
      stop("start must be a scalar, have dimension p x m, or have length p for univariate y")
    # if (normalize) { # normalisation is done in L0glm.fit.internal
    #   Xscaled <- normalizeX(X, no.norm = intercept.col)
    #   scaling <- attr(Xscaled, "scaling")
    #   # start   <- start * scaling
    # }
  }
  if (!is.null(start) && any(!is.finite(start)))
    stop("start must contain finite coefficients.", call. = FALSE)
  start <- L0glm_unit_start(start, X, preprocessing, m)
  allowed <- if (is.null(start)) NULL else matrix(as.numeric(start) != 0, p, m)
  if (!is.null(allowed)) allowed[always.include, ] <- TRUE
  
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
  
  # 5 dispersion: positive values are fixed, negative values are free starts.
  dispersion <- L0glm_dispersion(dispersion, family, m)
  estimate_dispersion <- if (is.null(dispersion)) rep(TRUE, m) else dispersion < 0
  disp_list <- as_list(dispersion)
  next_dispersion <- function(d) ifelse(estimate_dispersion, -d, d)
  
  # 6 misc --------------------------------------------------
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
      convergence <- list()
      coefficients <- lambdas.mat <- Matrix(0, p, m, sparse = TRUE)
      weights.mat <- z.mat <- Matrix(0, n, m, sparse = TRUE)
      if (do_parallel) {
        nworkers <- min(ncores, m)
        cluster_type <- if (Sys.info()["sysname"] == "Windows") "PSOCK" else "FORK"
        cl <- if (cluster_type == "PSOCK") {
          parallel::makeCluster(nworkers, type = cluster_type, useXDR = FALSE,
                                master = "127.0.0.1", rscript_args = "--vanilla")
        } else {
          parallel::makeCluster(nworkers, type = cluster_type)
        }
        on.exit(parallel::stopCluster(cl), add = TRUE)
        L0glm_dll <- getLoadedDLLs()[["L0glm"]]
        L0glm_dll_path <- if (!is.null(L0glm_dll)) L0glm_dll[["path"]] else NULL
        L0glm_pkg_path <- getNamespaceInfo("L0glm", "path")
        L0glm_is_dev <- requireNamespace("pkgload", quietly = TRUE) &&
          pkgload::is_dev_package("L0glm")
        worker_init <- L0glm_worker_init
        environment(worker_init) <- baseenv()
        parallel::clusterCall(cl, worker_init, .libPaths(), L0glm_dll_path,
                              worker_nthreads, L0glm_pkg_path, L0glm_is_dev)
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
               start      = if (is.null(start.mat) || length(start.mat) == 1L)
                 start.mat else start.mat[vars, cols, drop = FALSE],
               offset     = offset[, cols, drop = FALSE])
        })
      }
      collect_fit <- function(fit, chunked) {
        if (chunked) {
          list(coefficients = do.call(cbind, lapply(fit, `[[`, "coefficients")),
               weights      = do.call(cbind, lapply(fit, `[[`, "weights")),
               dispersion   = unlist(lapply(fit, `[[`, "dispersion"), use.names = FALSE),
               lambdas      = do.call(cbind, lapply(fit, `[[`, "lambdas")),
               convergence  = do.call(rbind, lapply(fit, `[[`, "convergence")),
               z            = do.call(cbind, lapply(fit, `[[`, "z")))
        } else {
          list(coefficients = do.call(cbind, lapply(fit, `[[`, "coefficients")),
               weights      = do.call(cbind, lapply(fit, `[[`, "weights")),
               dispersion   = as.numeric(Matrix(do.call(cbind, lapply(fit, `[[`, "dispersion")),
                                                sparse = FALSE)[1, ]),
               lambdas      = do.call(cbind, lapply(fit, `[[`, "lambdas")),
               convergence  = do.call(rbind, lapply(fit, `[[`, "convergence")),
               z            = do.call(cbind, lapply(fit, `[[`, "z")))
        }
      }
      group_support_start <- function(coefs) {
        support <- L0glm_group_start(coefs, nonnegative)
        support$vars <- sort(unique(c(support$vars, always.include.start)))
        if (!is.null(allowed)) support$start[!allowed] <- 0
        support
      }
      
      for (iter in seq_len(group.maxit)) {
        if (iter > 1) {
          support <- group_support_start(coefficients)
          vars   <- support$vars
          always.include <- match(always.include.start[always.include.start %in% vars], vars)
          start  <- support$start
        } else {
          vars <- if (is.null(allowed)) seq_len(p) else which(rowSums(allowed) > 0)
          always.include <- match(always.include.start, vars)
        }
        
        if (!length(vars) && iter > 1L) {
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
                                        thresh         = thresh,
                                        solver        = solver,
                                        lower         = lower[vars],
                                        upper         = upper[vars],
                                        method        = method,
                                        normalize     = normalize,
                                        preprocessing = L0glm_subset_preprocessing(preprocessing, vars),
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
                               preprocessing = L0glm_subset_preprocessing(preprocessing, vars),
                               offset      = offset[, col],
                               prior.mean  = prior.mean[vars, col],
                               dispersion  = disp_list[[col]],          # slice
                               start       = if (is.null(start) || length(start) == 1L)
                                 start else start[vars, col, drop = FALSE],
                               lower       = lower[vars],
                               upper       = upper[vars],
                               maxit       = maxit,
                               thresh       = thresh,
                               solver      = solver,
                               method      = method,
                               tol         = tol,
                               always.include      = always.include,
                               control.solver = control.solver,
                               verbose     = verbose)
          })
        }
        
        fit <- collect_fit(fit, chunked = do_parallel)
        convergence[[length(convergence) + 1L]] <-
          cbind(stage = "selection", pass = iter, channel = seq_len(m), fit$convergence)
        if (length(vars) == p && identical(vars, seq_len(p))) {
          coefficients <- Matrix(fit$coefficients, sparse = TRUE)
        } else {
          coefficients[vars, ] <- Matrix(fit$coefficients, sparse = TRUE)
        }
        weights.mat          <- Matrix(fit$weights,      sparse = TRUE)
        dispersion           <- fit$dispersion
        disp_list            <- as_list(next_dispersion(dispersion))
        if (length(vars) == p && identical(vars, seq_len(p))) {
          lambdas.mat <- Matrix(fit$lambdas, sparse = TRUE)
        } else {
          lambdas.mat[vars, ]  <- Matrix(fit$lambdas, sparse = TRUE)
        }
        z.mat                <- Matrix(fit$z,            sparse = TRUE)

      }
      
      # ----- relaxed - refit selected variables nearly unpenalised -----------
      if (relaxed) {
        disp_list <- as_list(next_dispersion(dispersion))
        support <- group_support_start(coefficients)
        vars   <- support$vars
        if (length(vars)) {
          always.include <- match(always.include.start[always.include.start %in% vars], vars)
          start  <- support$start
          coefficients <- lambdas.mat <- Matrix(0, p, m, sparse = TRUE)
          weights.mat  <- Matrix(0, n, m, sparse = TRUE)
          
          if (do_parallel) {
            chunks <- make_chunks(seq_len(m), length(cl))
            jobs <- make_jobs(chunks, vars, start)
            fit <- parallel::clusterApply(cl, jobs, L0glmfun_chunk,
                                          X             = X[, vars, drop = FALSE],
                                          family        = family,
                                          lambda        = 1e-30, # using a very small lambda here
                                          hard_threshold = FALSE,
                                          gamma         = gamma,
                                          maxit         = maxit,
                                          tol           = tol,
                                          always.include        = always.include,
                                          thresh         = thresh,
                                          solver        = solver,
                                          lower         = lower[vars],
                                          upper         = upper[vars],
                                          method        = method,
                                          normalize     = normalize,
                                          preprocessing = L0glm_subset_preprocessing(preprocessing, vars),
                                          control.solver = control.solver,
                                          verbose       = verbose)
          } else {
            fit <- lapply(seq_len(m), function(col) {
              L0glm.fit.internal(X          = X[, vars, drop = FALSE],
                                 y          = y[, col, drop = FALSE],
                                 weights     = weights[, col], 
                                 family     = family,
                                 lambda     = 1e-30,
                                 hard_threshold = FALSE,
                                 gamma      = gamma,
                                 normalize   = normalize,
                                 preprocessing = L0glm_subset_preprocessing(preprocessing, vars),
                                 offset      = offset[, col], 
                                 prior.mean  = prior.mean[vars, col], 
                                 dispersion  = disp_list[[col]],
                                 start       = start[vars, col, drop = FALSE],
                                 lower       = lower[vars],
                                 upper       = upper[vars],
                                 maxit       = maxit,
                                 thresh       = thresh,
                                 solver      = solver,
                                 method      = method,
                                 tol         = tol,
                                 always.include      = always.include,
                                 control.solver = control.solver,
                                 verbose     = verbose)
            })
          }
          
          fit <- collect_fit(fit, chunked = do_parallel)
          convergence[[length(convergence) + 1L]] <-
            cbind(stage = "relaxed", pass = 1L, channel = seq_len(m), fit$convergence)
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
           convergence  = do.call(rbind, convergence),
           z            = z.mat)
    } # ------- end fit_multichannel() --------------------------------------
    
    res_multi <- fit_multichannel(do_parallel = parallelize && ncores > 1)
    
    coefficients <- res_multi$coefficients
    weights.ret  <- res_multi$weights
    dispersion   <- res_multi$dispersion
    lambdas      <- res_multi$lambdas
    z            <- res_multi$z
    convergence  <- res_multi$convergence
    
  } else {
    # =======================================================================
    #  UNIVARIATE RESPONSE  (m == 1) ****************************************
    # =======================================================================
    fit <- L0glm.fit.internal(X          = X,
                              y          = as.vector(y),
                              weights    = weights,
                              family     = family,
                              normalize   = normalize,
                              preprocessing = preprocessing,
                              offset      = offset,
                              dispersion  = dispersion,
                              start       = start,
                              lower       = lower,
                              upper       = upper,
                              lambda      = lambda,
                              prior.mean  = prior.mean,
                              gamma       = gamma,
                              thresh       = thresh,
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
    convergence <- cbind(stage = "selection", pass = 1L, channel = 1L, fit$convergence)
    
    # --- relaxed - refit selected variables nearly unpenalised -------------
    if (relaxed) { 
      vars <- sort(unique(c(which(coefficients != 0), always.include.start)))
      if (length(vars)) {
        always.include <- match(always.include.start[always.include.start %in% vars], vars)
        fit     <- L0glm.fit.internal(X          = X[, vars, drop = FALSE],
                                      y          = y,
                                      weights    = weights,
                                      family     = family,
                                      normalize   = normalize,
                                      preprocessing = L0glm_subset_preprocessing(preprocessing, vars),
                                      offset      = offset,
                                      dispersion  = next_dispersion(dispersion),
                                      start       = matrix(coefficients[vars], ncol = 1L),
                                      lower       = lower[vars],
                                      upper       = upper[vars],
                                      lambda      = 1e-30,
                                      hard_threshold = FALSE,
                                      prior.mean  = prior.mean[vars],
                                      gamma       = gamma,
                                      thresh       = thresh,
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
        convergence <- rbind(convergence,
          cbind(stage = "relaxed", pass = 1L, channel = 1L, fit$convergence))
      }
    }
  }
  
  
  # ------------------------------------------------------------------------
  #  Final clean-up  &  return
  # ------------------------------------------------------------------------
  if (any(is.infinite(weights.ret))) {
    weights.ret[is.infinite(weights.ret)] <- 1e11
  }
  if (!is.null(allowed)) {
    if (m == 1L) lambdas[!as.vector(allowed)] <- Inf else lambdas[!allowed] <- Inf
  }
  if (any(is.infinite(lambdas))) {
    lambdas[is.infinite(lambdas)] <- 1e11
  }
  
  ## ensure coefficients are returned as sparse coeff-matrix also for univariate case
  coefficients <- Matrix::Matrix(matrix(as.numeric(coefficients), p, m), sparse = TRUE,
                                  dimnames = list(colnames(X), if (m > 1L) colnames(y) else NULL))
  coefficients <- as(as(coefficients, "generalMatrix"), "CsparseMatrix")
  
  list(coefficients = coefficients,
                 weights      = weights.ret,
                 dispersion   = dispersion,
                 init.lambda  = init.lambda,
                 lambdas      = lambdas,
                 family       = family,
                 p            = p_selectable,
                 always.include = always.include.start,
                 preprocessing = preprocessing,
                 convergence = convergence,
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
  if (!is.na(physical_cores)) return(max(1L, round(physical_cores/2)))
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

L0glm_quiet_try <- function(expr) {
  try(
    suppressMessages(suppressWarnings(utils::capture.output(force(expr)))),
    silent = TRUE
  )
  invisible(NULL)
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
    L0glm_quiet_try(RhpcBLASctl::blas_set_num_threads(nthreads))
    L0glm_quiet_try(RhpcBLASctl::omp_set_num_threads(nthreads))
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
      L0glm_quiet_try(RhpcBLASctl::blas_set_num_threads(state$blas_threads))
    }
    if (!is.na(state$omp_threads)) {
      L0glm_quiet_try(RhpcBLASctl::omp_set_num_threads(state$omp_threads))
    }
  }
  invisible(NULL)
}

L0glm_worker_init <- function(lib_paths, dll_path, nthreads = 1L,
                              package_path = NULL, is_dev_package = FALSE) {
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
    try(
      suppressMessages(suppressWarnings(
        utils::capture.output(RhpcBLASctl::blas_set_num_threads(nthreads))
      )),
      silent = TRUE
    )
    try(
      suppressMessages(suppressWarnings(
        utils::capture.output(RhpcBLASctl::omp_set_num_threads(nthreads))
      )),
      silent = TRUE
    )
  }
  if (!is.null(dll_path) && nzchar(dll_path) &&
      file.exists(dll_path) &&
      !("L0glm" %in% names(getLoadedDLLs()))) {
    dyn.load(dll_path)
  }
  if (isTRUE(is_dev_package) && !is.null(package_path) &&
      nzchar(package_path) && file.exists(file.path(package_path, "DESCRIPTION")) &&
      requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(package_path, quiet = TRUE)
  } else if (!("L0glm" %in% loadedNamespaces())) {
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
       convergence  = do.call(rbind, lapply(fit, `[[`, "convergence")),
       z            = do.call(cbind, lapply(fit, `[[`, "z")))
}

L0glmfun <- function(cols, lambda, gamma, thresh, solver, method,
                     normalize,
                     X, y, family,
                     prior.mean = NULL, weights = NULL, dispersion = NULL,
                     start = NULL, lower, upper,
                     verbose = FALSE, maxit, tol, offset, always.include = 0,
                     control.solver = do.call("control.solver.gen",
                                              list(solver = solver, method = method)),
                     hard_threshold = TRUE, preprocessing = NULL) {
  
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
      start        = if (is.null(start) || length(start) == 1L)
        start else start[, col, drop = FALSE],
      maxit        = maxit,
      tol          = tol,
      always.include       = always.include,
      lower        = lower,
      offset       = as.vector(offset[, col]),
      upper        = upper,
      prior.mean   = prior.mean[, col],
      thresh        = thresh,
      hard_threshold = hard_threshold,
      solver       = solver,
      normalize    = normalize,
      preprocessing = preprocessing,
      method       = method,
      verbose      = verbose,
      control.solver = control.solver)
  })
}
  
