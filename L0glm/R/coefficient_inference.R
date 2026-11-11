#' Coefficient inference after sparse selection
#'
#' @name coefficient_inference
#' @description
#' Optional posterior, selected-support or case-bootstrap uncertainty summaries
#' for ordinary GLMs, negative binomial, multinomial and Cox models.
#' Request these calculations through \code{summary(fit, inference=...)}.
#' \code{simulate(fit)} and \code{simulate(posterior)} return posterior
#' coefficient draws in the original coefficient dimensions.
#' @param object A fitted \code{L0glm_fit} model, or an inference object of class
#'   \code{L0glm_posterior} returned by \code{summary(fit,inference="posterior")$inference}.
#' @param control.inference Named inference controls described below.
#'
#' @details
#' With \code{inference="posterior"}, the selected support and fitted noise
#' parameters are held fixed. When an accepted full support refit supplies the
#' returned estimates, or \code{lambda=0}, the calculation uses a flat prior on
#' the admitted coefficients. Otherwise it freezes one adaptive-ridge prior
#' update at the fitted coefficients. In original units its precision is
#' \eqn{P_{jj}=s_j^2\lambda/(4u_j^2)}, where \eqn{u_j=s_j\widehat\beta_j}
#' is the normalized coefficient. Shared row or predictor-group penalties use
#' their corresponding squared block norms and group multipliers. Mandatory
#' terms have zero prior precision; supplied \code{prior.mean} values give the
#' prior centres. Starting exclusions and the selected support remain fixed.
#'
#' The conditional target is
#' \deqn{p(\beta_J\mid y,J,P,\widehat\psi)\ \propto
#'   \exp\{\ell(\beta_J,\widehat\psi)
#'   -\tfrac12(\beta_J-c)^T P(\beta_J-c)\}\,1_{\mathcal C}(\beta_J).}
#' A separate fixed-prior Fisher-scoring solve finds its constrained mode.
#' The original fitted coefficients and predictions are preserved. The local
#' Gaussian uses precision \eqn{Q=\mathcal I+P}, where \eqn{\mathcal I} is
#' the IRLS information; multinomial and Cox use joint likelihood information.
#' At a constrained mode \eqn{b_*}, the negative log-target gradient \eqn{g_*}
#' can be nonzero. The approximating normal therefore has centre
#' \eqn{b_*-Q^{-1}g_*} and covariance \eqn{Q^{-1}}, truncated to the coefficient
#' constraints. Multinomial draws lie in the admitted zero-sum CLR subspace.
#'
#' Gaussian identity models with fixed dispersion give the exact frozen-prior
#' Gaussian target, subject to any truncation. Other families use a local
#' likelihood approximation. Quasi families use the corresponding working
#' posterior. Estimated dispersion and NB shape are plug-in values. The method
#' reports posterior means, standard deviations and equal-tail credible intervals,
#' together with the selected estimates and conditional modes. These intervals
#' condition on data-dependent selection and prior estimation; their frequentist
#' coverage after selection can be substantially below the nominal level.
#' Case bootstrapping repeats selection when that variability is of interest.
#'
#' Sampling uses Cholesky Gaussian draws, independent rejection when the
#' nonnegative region has sufficient probability, and native truncated-normal
#' Gibbs updates otherwise. The full precision is retained in the Gibbs
#' conditionals. Nonnegative draws are sampled from the constrained distribution.
#' Strong posterior correlations can require longer Gibbs runs; use
#' \code{burnin}, \code{thin} and repeated seeds to assess Monte Carlo stability.
#'
#' With \code{inference="conditional"}, the selected coefficients and compulsory
#' terms are refitted without a sparsity penalty, keeping the model family,
#' observation weights, offsets and coefficient constraints. The original fitted
#' object and its coefficients are preserved. Reported estimates and standard
#' errors refer to this selected-support refit. Its original penalized estimates
#' are also retained in the coefficient table.
#'
#' For ordinary GLMs, the model-based covariance is
#' \deqn{\widehat{\mathrm{Var}}(\widehat\beta_J)
#'   =\widehat\phi(X_J^T W X_J)^{-1},\qquad
#'   W_{ii}=w_i\{\mu'(\eta_i)\}^2/V(\mu_i).}
#' Supplied fixed dispersion is respected. Estimated dispersion uses the
#' selected-model Pearson scale and residual degrees of freedom. Unknown-scale
#' ordinary-family intervals use a Student t reference; known-scale, negative
#' binomial, multinomial and Cox intervals use a normal reference. NB shape is
#' reestimated when it was estimated in the fit, using plug-in expected
#' coefficient information. Fixed-noise scalar families use observed likelihood
#' information. Cox uses weighted Breslow partial-likelihood information.
#'
#' Multinomial inference works on the selected zero-sum coefficient subspace.
#' If \eqn{M} maps independent contrasts into the displayed CLR coefficients,
#' \deqn{\widehat{\mathrm{Var}}(\mathrm{vec}\widehat B)
#'   =M(M^T\mathcal I M)^{-1}M^T.}
#' Grouped selection refits the complete CLR rows. Ungrouped selection preserves
#' individual selected CLR entries and the row-sum constraint. Structurally
#' excluded coefficients remain zero. Cholesky inversion uses diagonal
#' equilibration; singular information and coefficients on nonnegative boundaries
#' receive unavailable ordinary Wald standard errors.
#'
#' These are fixed-model intervals after selection: their nominal interpretation
#' treats the selected model as given. Selection itself contributes uncertainty.
#' Strong-signal settings with reliable support recovery approach the usual
#' fixed-model setting; a finite bootstrap cannot certify perfect support recovery.
#'
#' With \code{inference="bootstrap"}, observation rows are sampled with
#' replacement, together with their weights and offsets. Every replicate reruns
#' the saved fitting specification, including its screening, interactions, adaptive-ridge
#' fitting and local refinement. Independent validation data remain fixed.
#' Explicit starting exclusions and externally supplied screening or clustering
#' settings retain their specified meanings. A fit extracted from tuning keeps
#' its chosen lambda; a custom \code{refit} callback can repeat tuning as well.
#'
#' Bootstrap standard errors and percentile intervals include zero estimates in
#' replicates where a coefficient is absent. Selection frequencies are reported
#' for all coefficients selected in the original fit or any replicate. Predictors
#' are retained in every search according to the original fitting specification;
#' there is no pruning based on nonselection in a limited number of replicates.
#' The bootstrap describes variability and selection stability. Discontinuities
#' at selection thresholds can affect nominal interval coverage; see Leeb and
#' Poetscher (2008).
#'
#' @section Controls:
#' \code{control.inference} is a named list with:
#' \describe{
#' \item{\code{level=.95}}{Interval probability.}
#' \item{\code{maxit=100L}, \code{tol=1e-10}}{Selected-support refit controls.}
#' \item{\code{max.parameters=1000L}, \code{rank.tol=1e-10}}{Information-matrix
#' dimension guard and numerical rank tolerance.}
#' \item{\code{nboot=199L}, \code{seed=NULL}}{Bootstrap replicate count and
#' optional seed. A supplied seed preserves the caller's random-number state.}
#' \item{\code{nsim=2000L}, \code{burnin=500L}, \code{thin=1L}}{Posterior-summary
#' draw count and Gibbs burn-in/thinning. The \code{simulate()} methods use their
#' \code{nsim} argument. Independent sampling needs no burn-in or thinning.}
#' \item{\code{max.cells=5e7}}{Memory guard for row resampling and dense bootstrap
#' summary blocks. Coefficient draws are stored sparsely.}
#' \item{\code{refit.args=list()}}{Overrides to the saved fitting options for
#' bootstrap replicates.}
#' \item{\code{refit=NULL}}{Optional function receiving a fitting-argument list
#' and returning a fitted L0glm model, or a tuning result with \code{fit_best}.}
#' }
#'
#' @section Returned quantities:
#' \code{summary(fit)$coefficients} retains the original sparse point estimates.
#' With inference requested, \code{summary(fit)$inference$coefficients} is the
#' coefficient table, including estimates, standard errors, interval limits and
#' either credible limits, Wald tests or bootstrap selection frequencies. Original row and channel
#' indices identify each coefficient. Conditional inference also returns the
#' selected-support refit, covariance matrices and convergence diagnostics.
#' Posterior output contains \code{prior.source}, \code{prior.precision},
#' \code{prior.mean}, frozen dispersion/shape, \code{sampler}, covariance and
#' sparse coefficient \code{draws}. \code{gaussian.covariance} is the covariance
#' before truncation; constrained \code{covariance} is estimated from the draws.
#' It stores a compact sampling state for \code{simulate(posterior)}, excluding
#' the training design and response. Draw rows flatten the original coefficient
#' matrix in column-major order; \code{coefficient.dim},
#' \code{coefficient.dimnames} and \code{coefficient.index} attributes identify
#' the original coordinates. Excluded coefficients have zero draws.
#'
#' Bootstrap output retains sparse \code{draws}, a sparse original-dimension
#' \code{selection.frequency} matrix (extended by the union of generated
#' interaction terms), any generated interaction terms and
#' per-replicate \code{diagnostics} with failures, convergence, warnings and
#' timings. Draw rows flatten coefficient matrices in column-major order.
#' \code{replicates} identifies successful columns; failed fits are reported and
#' excluded explicitly. At least two successful fits are required.
#' Covariance matrices follow the selected terms in each channel's table;
#' multinomial covariance uses column-major CLR order on \code{covariance.rows}.
#' Fitted models retain evaluated options in \code{specification}, excluding
#' copies of the design and response; \code{engine} identifies the fitter.
#'
#' @seealso \code{\link{model_summary}}, \code{\link{l0glmFit}}
#' @references
#' Leeb, H. and Poetscher, B. M. (2008). Sparse estimators and the oracle property,
#' or the return of Hodges' estimator. Journal of Econometrics, 142, 201-211.
#' \doi{10.1016/j.jeconom.2007.05.017}.
#' Robert, C. P. (1995). Simulation of truncated normal variables.
#' Statistics and Computing, 5, 121-125. \doi{10.1007/BF00143942}.
#' @examples
#' set.seed(51)
#' X <- matrix(rnorm(200*10),200,10)
#' y <- 1+2*X[,1]-2*X[,2]+rnorm(200)
#' fit <- l0glmFit(X,y,intercept=TRUE,lambda="gic")
#' details <- summary(fit,inference="conditional")
#' details$inference$coefficients
#' posterior <- summary(fit,inference="posterior",
#'   control.inference=list(nsim=1000L,seed=1))$inference
#' posterior$coefficients
#' draws <- simulate(posterior,nsim=1000L,seed=2)
#' matrix(draws[,1L],nrow(coef(fit)),ncol(coef(fit)))
#' \dontrun{
#' uncertainty <- summary(fit,inference="bootstrap",
#'   control.inference=list(nboot=199L,seed=1))
#' uncertainty$inference$coefficients
#' }
NULL
