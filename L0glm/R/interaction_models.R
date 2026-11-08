#' Staged selection of interaction terms
#'
#' Add products of selected predictors to an L0glm model. Each stage uses the
#' same family, constraints, weights, offsets and sparse-selection options as
#' the main-effect fit. Products are constructed only for the selected parents,
#' keeping the working problem small even when the original design is wide.
#'
#' @name interaction_models
#' @param interactions Nonnegative integer controlling staged interaction
#'   generation. The default \code{0L} fits main effects.
#'   \code{1L} also considers pairwise products of selected main effects;
#'   \code{2L} then considers three-way products, and so on. Each product contains distinct
#'   original predictors. Intercepts, \code{always.include} and predictor
#'   group-0 terms remain as compulsory main effects. See Details.
#' @param control.interactions Named advanced controls: \code{max.terms=10000L}
#'   bounds the total number of generated products; \code{max.cells=5e7}
#'   bounds rows times columns in a working design; \code{tol=1e-8} is the
#'   relative tolerance for stopping when validation deviance worsens; and
#'   \code{criterion=NULL} chooses the information criterion used to select
#'   the final stage. A named \code{lambda} supplies that criterion; otherwise
#'   GIC is used. With \code{l0glmTune}, the tuning criterion is used, or GIC
#'   inside each CV training fit. Limits are checked before forming products.
#'
#' @details
#' The main-effect stage first runs the requested fit, including SIS/ISIS when
#' enabled. At stage \eqn{d-1}, products of \eqn{d} currently selected original
#' main effects are generated. Compulsory terms are retained and previously
#' generated products remain candidates. Selection is repeated on this reduced
#' design. Parents determine which products are considered; an interaction can
#' remain selected after a parent main effect is discarded. Variable powers
#' require explicitly supplied predictor columns.
#'
#' For a named criterion preset, the stage penalty uses the original selectable
#' predictor count plus the products generated so far. For example,
#' \deqn{\lambda_s=\log(p_0+q_s)\log\log(n)}
#' for GIC, where \eqn{p_0} excludes compulsory terms and \eqn{q_s} counts
#' generated products. Numeric \code{lambda} stays fixed across stages.
#' Final stage scores use the common candidate count \eqn{p_0+q_{\rm final}}
#' and the usual free-coefficient degrees of freedom. With predictor groups,
#' products are grouped by the tuple of their parent group labels; named
#' presets use coefficient-dimension group costs, as in
#' \link{predictor_groups}. With independent validation data, its deviance
#' selects the stage and generation stops when this deviance worsens.
#' For quasi families with unavailable absolute likelihood criteria, stage
#' selection uses scaled quasi-deviance plus the criterion's parameter penalty,
#' holding the main-effect dispersion fixed. Their public likelihood criteria
#' retain the availability reported by the supplied family.
#'
#' Scalar starts retain their usual interpretation at every stage. Explicit
#' original-coordinate vectors or matrices supply main-effect starts; new
#' products are seeded at one in normalized coefficient coordinates, or a
#' centered contrast for multinomial models. A product is excluded in a channel
#' if any parent was explicitly excluded by a starting zero. Multinomial starts
#' use row exclusions. New products have prior mean zero.
#'
#' Cross-validation regenerates main-effect selection and interactions within
#' each training fold. Prediction takes the original main-effect design and
#' constructs the required products. Numeric and \code{NULL} starts are supported
#' for interaction tuning. Returned path coefficients share the union of the
#' full-data fits' generated products, with zeros for absent terms; path IC
#' scores use this common candidate count. Fold-specific products are retained
#' in \code{cv$interaction_terms}.
#'
#' The stored training design is a lightweight view of the original matrix and
#' product recipe. Subsetting \code{fit$X[, j, drop=FALSE]} materializes only
#' requested columns. Dense inputs and sparse inputs retain their storage type.
#' Preclustering uses original-channel refinement (\code{cluster.refit=TRUE});
#' an aggregate-only analysis can pass aggregated responses directly.
#'
#' @section Returned recipe:
#' \code{fit$interactions} contains \code{order}, \code{selected_order},
#' \code{terms} (lists of original column indices), predictor names, the fitting
#' \code{specification}, and a \code{history} table with each stage's lambda,
#' degrees of freedom and selection score. Sparse coefficients contain every
#' original main effect and generated product, including excluded terms as zero.
#' To reproduce the fitted pipeline, pass its specification together with the
#' original \code{X} and \code{y} to \code{do.call(l0glmFit, ...)}. A tuning
#' path also includes passive zero rows from other lambda fits for display.
#'
#' @seealso \code{\link{l0glmFit}}, \code{\link{l0glmTune}},
#'   \link{predictor_groups}, \link{model_predictions}
#' @examples
#' set.seed(31)
#' X <- matrix(rnorm(250 * 8), 250, 8)
#' y <- 1 + 2 * X[, 1] - 2 * X[, 2] +
#'   2 * X[, 1] * X[, 2] + rnorm(250)
#' fit <- l0glmFit(X, y, family=gaussian(), intercept=TRUE,
#'   lambda="gic", interactions=1L, relaxed=TRUE)
#' fit$interactions$history
#' predict(fit, X[1:5, ])
#' \dontrun{
#' cv <- l0glmTune(X, y, family=gaussian(), intercept=TRUE,
#'   interactions=1L, tune="cv", nfolds=5L, relaxed=TRUE)
#' predict(cv$fit_1se, X[1:5, ])
#' }
NULL
