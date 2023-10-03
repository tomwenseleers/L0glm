#' @name L0glm.fit
#' @aliases L0glm.fit
#' @title Fit generalized linear model with \eqn{L_0} or grouped \eqn{L_{0,∞}} regularization subject to optional 
#'   nonnegativity or box constraints
#'
#' @description
#' Fits a generalized linear model with \eqn{L_0} norm or grouped \eqn{L_{0,∞}} mixed norm
#' regularization to approximate Best Subset Selection (for univariate data) or
#' Best Subset of Group Selection (for multivariate data) subject to optional
#' nonnegativity or box constraints on model coefficients. The objective
#' functions that can be approximated are :
#' \itemize{
#' 
#' \item for univariate data, using an \eqn{L_0} norm penalty to achieve Best Subset Selection:
#' 
#' \deqn{-2.LL + λ.||β||_0}
#' 
#' where \eqn{LL} is the log-likelihood of the GLM model and
#' \eqn{|β||_0} is the number of non-zero elements in
#' coefficient vector β, i.e. the \eqn{L_0} norm of β.
#' 
#' \item for multivariate data, using an \eqn{L_{0,0}} norm penalty, i.e. using independent \eqn{L_0}  penalized fits on each outcome:
#' 
#' \deqn{-2.LL + λ.||β||_{0,0}}
#' 
#' \item for multivariate data, using an \eqn{L_{0,∞}} norm penalty, i.e. using grouped variable selection in multitask learning settings to achieve Best Subset of Groups Selection:
#' 
#' \deqn{-2.LL + λ.||β||_{0,∞}}
#' 
#' }
#' 
#' in all 3 cases subject to optional box constraints on the model coefficient vector or matrix β 
#' 
#' \deqn{lower ≤ β ≤ upper}
#' 
#' These objectives are approximated via an iterative adaptive ridge regression
#' procedure, see Details. The \eqn{L_0} and \eqn{L_{0,0}} norm penalties can
#' be used to approximate Best Subset Selection, while for multivariate
#' (multitask/multichannel) data, the \eqn{L_{0,∞}} mixed norm penalty can be
#' used to approximate Best Subset of Group Selection (via option
#' \code{grouped=TRUE}). Fits on multichannel data can be parallelized. Option
#' \code{relaxed} allows one to finish with a (nearly) unpenalized fit on the
#' \code{L0glm} selected features.
#' 
#' The level of regularisation can be manually chosen (if a numeric value is
#' passed to \code{lambda}) or preset to optimise given information criteria (if
#' set to \code{"aic"}, \code{"bic"}, \code{"gic"}, \code{"ebic"} or
#' \code{"mbic"}). For example, a λ of 2 is expected to optimise AIC, a
#' λ of log(\emph{n}) BIC, etc. 
#'
#' The \code{L0glm.fit} function has various advanced features, including
#' support for both dense and sparse design matrices \emph{\bold{X}}, high
#' dimensional data (where \emph{p} >> \emph{n}), single or multichannel outcomes \emph{\bold{y}} 
#' (with either independent or grouped variable selection),
#' support for any of the standard R GLM families & link functions, prior
#' weights (e.g. total nr. of trials for binomial proportions), model offsets,
#' starting coefficient values or dispersion (if prior estimate is available,
#' otherwise optimised), nonnegativity or box constraints, unpenalized terms
#' (e.g. intercept) and shrinkage towards given prior means (via nonzero
#' centered ridge regression updates).
#'
#' A choice of several fast C++ coded least squares solvers can be used, which
#' means this function scales well to solve fairly large tall or wide problems with
#' millions of observations or variables. 
#' As solvers, we currently implemented
#' the \code{Eigen} sparse or dense Cholesky, QR, PCG & LSCG solvers, plus two
#' solvers that natively support box constraints: \code{glmnet} (\code{bigGLM}), which
#' uses coordinate descent, and \code{osqp}, a fast quadratic programming
#' solver.
#' 
#' @param X
#' Design matrix with predictor variables of dimension \code{n} x \code{p}, with \code{n}
#' corresponding to the number of observations and \code{p} being the number of
#' columns in \code{y} (variables plus intercept). It should include a column of 
#' ones for an intercept term is desired. 
#' @param y
#' Response variable(s). Can be either a vector of length \code{n} or a matrix
#' of dimension \code{n} x \code{m}, where \code{n} is the number of observations and 
#' \code{m} the number of response variables (tasks/channels). 
#' @param family
#' a description of the error distribution (gaussian, poisson, binomial,
#' Gamma, inverse.gaussian) and link function (identity, log, sqrt, logit, cloglog, probit, cauchit, inverse, 1/mu^2) to be
#' to be used for the GLM model. See \code{?stats::family} for details of family functions.
#' @param weights 
#' Optional vector of prior observation weights ("prior weights") of length \code{n}, 
#' or a matrix of weights with dimension \code{n} x \code{m}, with \code{n} corresponding to 
#' the number of observations and \code{m} the number of response variables (columns in \code{y}).
#' For binomial proportion data this would be the total number of trials per
#' observed proportion (passed in \code{y}).
#' @param dispersion
#' Optional dispersion estimate; set to \code{NULL} by default, in which case it is initialised
#' at 1 and then estimated during the IRLS loop; fixed at 1 for family binomial & Poisson.
#' @param start
#' Optional vector of coefficient starting values of length \code{p}, or a
#' matrix of starting values with dimension \code{p} x \code{m}, with \code{p}
#' corresponding to the number of predictors and \code{m} to the number of
#' columns in \code{y}. These will determine the level of penalisation of the
#' different coefficients during the first IRLS ridge regression update (the
#' ridge penalty used equals \code{(lambda/4)*dispersion/abs(start)^gamma)}).
#' Note that coefficients with absolute value smaller than \code{thresh} will be
#' clipped to zero and are removed from the model. Set to \code{NULL} by default, in
#' which case it is auto-initialised.
#' @param offset 
#' Optional vector of offset values of length \code{p}, or a matrix of offset values
#' with dimensions \code{c(p, m)}. The offset is a known quantity that is 
#' added to the linear predictor of the model.
#' @param lambda
#' Either a scalar with the \eqn{L_0} or \eqn{L_{0,0}} penalty (when \code{grouped=FALSE}) 
#' or \eqn{L_{0,∞}} penalty (when \code{grouped=TRUE}), or a character with the name of an
#' information criterion (\code{"aic"}, \code{"bic"}, \code{"gic"}, \code{"ebic"} or \code{"mbic"}), 
#' in which case \code{lambda} will be preset with a value that would be expected to
#' optimise that particular information criterion. Set \code{lambda = 0} 
#' to disable penalization.
#' @param no.pen
#' A vector of coefficient indices that should 
#' not be penalized. Coefficients with a name that includes \code{"Intercept"}
#' in their name are not penalized by default. Specify 0 if all terms
#' should be allowed to be shrunk to zero.
#' @param nonnegative
#' A logical indicating whether nonnegativity constraints should be applied.
#' If \code{TRUE}, \code{lower} will be set to \code{rep(0, ncol(X))}.
#' @param normalize
#' Logical indicating whether the columns of the design matrix \strong{X} 
#' should be scaled to have unit \eqn{L_2} norm.
#' It is advised to use normalization unless the covariates were measured on a
#' similar scale.
#' @param solver
#' Solver used for the ridge regression step in the IRLS algorithm (see
#' details). It should be one of \code{"osqp"} (osqp quadratic programming
#' solver), \code{"glmnet"} (\code{"bigGlm"}, using coordinate descent) or
#' \code{"eigen"} (one of the \code{RcppEigen Eigen} solvers). 
#' \code{"osqp"} and \code{"glmnet"} support box or nonnegativity constraints
#' and scale well to larger problems, with especially \code{"glmnet"} being suitable
#' for very large problems. The \code{"eigen"} solvers achieve box or nonnegativity
#' constraints via clipping, which in combination with the ridge penalties in
#' practice work almost as well, and some of those solvers are often among
#' the fastest (e.g. \code{method} \code{3L}, \code{5L} or \code{7L}, i.e.
#' LSCG, CG or Cholesky LLT).
#' @param method
#' When \code{solver="eigen"}, method used for the eigen solver (see details). 
#' It should be an integer from 1 to 8 :
#' 
#' \code{1L} = QR solver, Householder QR solver if \code{X} is dense
#' 
#' \code{2L} = QR solver, QR solver with column pivoting if \code{X} is dense
#' 
#' \code{3L} = Least Squares Conjugate Gradient (LSCG) with identity preconditioner (i.e. not preconditioned)
#' 
#' \code{4L} = Least Squares Conjugate Gradient (LSCG) with least squares diagonal preconditioner
#' 
#' \code{5L} = Conjugate Gradient (CG) with identity preconditioner (i.e. not preconditioned)    
#' 
#' \code{6L} = Preconditioned Conjugate Gradient (PCG) with diagonal preconditioner
#' 
#' \code{7L} = Cholesky LLT solver
#' 
#' \code{8L} = Cholesky LDLT solver
#' 
#' Defaults to \code{7L} (Cholesky LLT), which is often the fastest.
#' 
#' @param lower
#' Optional lower limit on model coefficients to specify box constraints. 
#' It should either be a scalar, a vector of length 
#' \code{p} or a matrix with dimensions \code{c(p, m)}.
#' @param upper
#' Optional upper limit on model coefficients to specify box constraints.  
#' It should either be a scalar, a vector of length 
#' \code{p} or a matrix with dimensions \code{c(p, m)}.
#' @param parallelize
#' A logical value indicating whether in the multivariate case,
#' the columns of the response matrix \code{y} 
#' should be fitted in parallel. If \code{TRUE}, each 
#' column of \code{y} is fitted concurrently, potentially speeding up computation.
#' If \code{FALSE}, the columns are fitted sequentially.
#' @param ncores
#' Optional integer specifying the number of cores to use when \code{parallelize=TRUE}.
#' If \code{parallelize=TRUE}, defaults to the number of available logical cores
#' divided by 2, which usually corresponds to the number of physical cores available.
#' @param maxit
#' Integer specifying the maximum number of IRLS iterations to perform.
#' @param grouped
#' A logical value indicating whether a grouped fit should be performed in the
#' multivariate case or not. 
#' If \code{TRUE}, the algorithm will return a group penalized fit, otherwise
#' independent L0glm fits will be performed on each response. These procedures
#' approximate the \eqn{L_{0,∞}} and \eqn{L_{0,0}} mixed norm penalties, and can be
#' used to do Best Subset of Group Selection or plain Best Subset Selection. 
#' @param group.maxit
#' Integer  specifying the maximum number of iterations to perform for the group 
#' penalized fitting procedure when \code{grouped = TRUE}. In each iteration,
#' the columns of the response matrix \code{y} will be refitted using the average 
#' of all columns of \code{y} as the starting value, which then imposes a group penalty. 
#' @param relaxed
#' Logical indicating if the algorithm should conduct a nearly unpenalized fit 
#' on the L0glm selected variables (using \code{lambda=1E-10}).
#' This can eliminate any remaining bias on the coefficients and increase
#' sensitivity, especially for the multivariate case with grouped variable
#' selection. Defaults to \code{TRUE} when \code{grouped=TRUE}.
#' @param warm.start 
#' A numeric value in the range of [0, 1] that represents the level of confidence 
#' in the initial start values when provided. This parameter is used in a weighted average to 
#' calculate the 'start' values. Default is set at 0.1. 
#' The formula used is: 
#' \deqn{start = \text{warm.start} \cdot start + (1 - \text{warm.start}) \cdot \text{rep}(\text{sd}(start)^{\frac{1}{\gamma}}, \text{ncol}(X))}
#' If `warm.start` is closer to 1, more weight is given to the original 'start' value.
#' Using this weighted average avoids the problem that coefficients 
#' initialized at zero would otherwise always be removed from the model.
#' @param prior.mean
#' An optional numeric vector of length \code{p} specifying the center towards
#' which coefficients should be shrunk in each IRLS ridge regression iteration.
#' This is done via nonzero centered ridge regression. This can be useful if
#' there is good prior information available that some coefficients should have
#' a particular nonnegative value. The default \code{prior.mean = NULL} will
#' shrink to zero using a default zero centered ridge regression.
#' @param gamma
#' Scalar with the exponent to use when calculating the ridge penalty 
#' \code{(lambda/4)*dispersion/abs(beta)^gamma)} in each iteration.
#' Defaults to \code{gamma = 2}, as in Frommlet & Nuel 2016.
#' @param thresh
#' Coefficients below this scalar threshold are clipped to 0 in each of the
#' IRLS ridge regression iterations and removed from the model. Defaults to
#' \code{1E-5}.
#' @param tol
#' Convergence tolerance. Iterations are stopped when the relative change in the
#' \eqn{L_2} norm of the coefficients drops below this threshold. 
#' @param control.solver
#' List of parameters controlling the optimization solvers. Auto-initialised with  
#' \code{control.solver.gen} when \code{NULL}. See \code{?control.solver.gen} 
#' for more details on the available parameters for the different solvers.
#' @param verbose
#' Logical indicating if the algorithm progression should be printed to the console.
#'
#' @details
#' L0glm approximates the \eqn{L_0} or \eqn{L_{0,∞}} mixed norm 
#' penalized GLM objectives via an iterative adaptive ridge regression procedure.
#' 
#' \strong{Algorithm for univariate, single-channel data}
#' 
#' For univariate data, the approach is inspired by earlier work on L0 penalized
#' regression by de Rooi & Eilers (2011), Frommlet & Nuel (2016), Liu et al.
#' (2017), Dai et al. (2018 & 2020) & Li et al. (2021) and is based in a simple
#' modification of the IRLS (Iteratively Reweighed Least Squares)
#' algorithm used to fit GLMs (cf. McCullach &
#' Nelder 1989, Hastie & Pregibon 1992, Hardin & Hilbe 2007). In this algorithm,
#' each iteration of the IRLS algorithm updates the coefficient estimates using
#' a quadratic approximation of the log likelihood of the GLM objective. This is done
#' using a weighted least squares step. In our adapted algorithm, this step is
#' replaced by a (zero or nonzero centered) adaptive ridge regularised regression
#' with optional box or nonnegativity constraints on the coefficients. 
#' In each iteration, this uses the following ridge penalties for the penalized terms (excluding the
#' intercept, if included):
#' \[
#' Λ = (1/4).λ.σ^2.w 
#' \]
#' with adaptive penalty weights 
#' \[
#' w = 1/{{| β |^γ}}
#' \]
#' where \( β \) are the coefficient estimates in each iteration and
#' coefficients smaller than \code{thresh} are set to zero and taken out of the
#' model at the start of each iteration, and as in Frommlet & Nuel (2016), \(
#' γ \) is set to a default value of 2. Box or nonnegativity constraints
#' are natively supported if the \code{glmnet} (\code{bigGlm}) or the
#' \code{osqp} solver is used, or are obtained via clipping if the \code{eigen}
#' solvers are used. The factor of 1/4 was shown to be required to
#' approximate the L0 penalty objective and has been proven for the
#' Gaussian case and under orthogonality (Frommlet & Nuel 2016), but in practice
#' also work well to the non-orthogonal case and other noise families. In
#' contrast to Frommlet & Nuel (2016) and the algorithms given in Liu et al.
#' (2017) & Dai et al. (2018, 2020) (implemented in the \pkg{l0ara} package of
#' Guo et al. 2022), however, we do not require the dispersion (variance in case
#' of Gaussian linear models) to be known and standardized to 1, but instead
#' allow it to be estimated and updated during the IRLS algorithm (when
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
#' increase sensitivity. If an appropriate level of regularisation λ
#' is chosen, the L0 norm penalized solution approximates Best Subset
#' selection. To achieve variable selection consistency, the \code{"bic"} \code{lambda} 
#' preset is recommended for low dimensional settings, whilst \code{"gic"}, \code{"ebic"} or 
#' \code{"mbic"} may work better for high dimensional settings. Setting it to \code{"aic"} is 
#' recommended if predictive performance is of key interest, as
#' the AIC criterion asymptotically approximates the leave-one-out cross-validation
#' error.
#' 
#' \strong{Algorithm for multivariate/multichannel data}
#' 
#' For the multivariate case there are two options. First, if \code{grouped=FALSE},
#' an independent \code{L0glm} model is fit on each of the responses, which comes 
#' down to using an \eqn{L_{0,0}} mixed norm penalty. Second, if \code{grouped=TRUE},
#' grouped variable selection will be performed via an \eqn{L_{0,∞}} mixed norm penalty.
#' This is achieved by taking the \eqn{L_{0,0}} mixed norm penalized solution, and then,
#' for all variables that have nonzero coefficients in at least on of the channels, 
#' repeating the \(L_0) penalized fits but using the average of the coefficients across
#' all channels (response variables) as starting values so as to achieve a group 
#' penalty and perform grouped feature selection. This is then further iterated
#' until convergence. Finally, as for single channel data, if \code{relaxed=TRUE},
#' the algorithm can finish with a final step whereby the selected variables are 
#' refit using a nearly unpenalized model(with \code{lambda=1E-10}). 
#' Compared to other group penalized regression methods, such
#' as the group LASSO, which uses an \(L_0/L_2) mixed norm penalty, this approach
#' produces a sparser solution that better approximates the Best Subset of Groups solution
#' if λ is adequately chosen.
#' 
#' @return
#' Returns a list with the following elements:
#' \item{\code{coefficients}}{a vector or matrix of fitted coefficients.}
#' \item{\code{weights}}{a vector or matrix of observation weights.}
#' \item{\code{init.lambda}}{intial lambda.}
#' \item{\code{lambdas}}{a vector or matrix with the adaptive penalties lambdas.}
#' \item{\code{family}}{The error structure used to fit the data.}
#' \item{\code{z}}{a vector or matrix containing the working response in the 
#' final iteration of \code{L0glm.fit}}
#' \item{\code{X}}{design matrix X}
#' \item{\code{y}}{response vector or matrix y.}
#' 
#' @author
#' Tom Wenseleers, Laura Lanzas & Christophe Vanderaa
#' 
#' @references
#' \itemize{
#' \item J. de Rooi & P. Eilers, \emph{Deconvolution of pulse trains with the L0 penalty}. Analytica chimica acta 705, 218-226 (2011).
#' \item F. Frommlet & G. Nuel, \emph{An adaptive ridge procedure for L0 regularization}. PloS one 11, e0148620 (2016).
#' \item Z. Liu, F. Sun & D.P. McGovern, \emph{Sparse generalized linear model with L0 approximation for feature selection and prediction with big omics data}. BioData mining 10, 1-12 (2017).
#' \item L. Dai, K. Chen, Z. Sun, Z. Liu & G. Li, \emph{Broken adaptive ridge regression and its asymptotic properties}. Journal of multivariate analysis 168, 334-351 (2018).
#' \item L. Dai, K. Chen & G. Li, \emph{The broken adaptive ridge procedure and its applications}. Statistica Sinica 30, 1069-1094 (2020).
#' \item N. Li, X. Peng, E. Kawaguchi, M.A. Suchard & G. Li, \emph{A scalable surrogate L0 sparse regression method for generalized linear models with applications to large scale data}. Journal of Statistical Planning and Inference 213, 262-281 (2021).
#' \item W. Guo, S. Ma & Z. Liu, \emph{l0ara: Sparse Generalized Linear Model with L0 Approximation for Feature Selection}. R package version 0.1.6 (2022). Available at \url{https://cran.r-project.org/web/packages/l0ara/index.html}
#' \item J.W. Hardin & J.M. Hilbe, \emph{Generalized linear models and extensions}. (Stata press, 2007).
#' \item T.J. Hastie & D. Pregibon, \emph{Generalized linear models}. Chapter 6 of \emph{Statistical Models in S} eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole (1992).
#' \item P. McCullagh & J.A. Nelder, \emph{Generalized Linear Models}. London: Chapman and Hall (1989).
#' }
#'
#' @importFrom parallel detectCores makeCluster clusterExport clusterEvalQ clusterApply stopCluster
#'
#' @export
L0glm.fit = function(X, 
                     y, 
                     family = gaussian(identity),
                     weights = NULL,
                     dispersion = NULL,
                     start = NULL,
                     offset = NULL,
                     lambda = 1,
                     no.pen = if(!is.null(colnames(X))) which(grepl("Intercept", colnames(X))), # if there is an Intercept do not penalize it
                     nonnegative = FALSE,
                     normalize = FALSE,
                     solver = "eigen", 
                     method = 7L,
                     lower = if (nonnegative) { rep(0, ncol(X)) } else { rep(-.Machine$double.xmax, ncol(X)) },
                     upper = rep(.Machine$double.xmax, ncol(X)),
                     parallelize = TRUE,
                     ncores = if (parallelize) { round(parallel::detectCores()/2) } else { 1 },
                     maxit = 25,
                     grouped = TRUE, 
                     group.maxit = 6, 
                     relaxed = grouped,
                     warm.start = 0.1,
                     prior.mean = NULL,
                     gamma = 2,
                     thresh = 1E-5, 
                     tol = 1E-10,
                     control.solver = do.call("control.solver.gen", list(solver = solver, method = method)),
                     verbose = FALSE) {
  if (lambda %in% c("aic","bic","gic","ebic","mbic")) { lambda = preset.lambda(IC=lambda, p=ncol(X), n=nrow(X)) }
  y = as.matrix(y)
  if (issparse(X)) { X = as(Matrix(X, sparse = TRUE), "generalMatrix") } else { X = Matrix(X, sparse = FALSE) }
  no.pen.start <- no.pen
  init.lambda <- lambda
  # Make weights, offset and start into a matrix if they are provided as a vector and y is multichannel.
  if(!is.null(weights)){
    if(ncol(y) > 1 & ncol(as.matrix(weights)) == 1){
      weights = matrix(rep(weights, times = ncol(y)), nrow = length(weights), ncol = ncol(y))
    }
  } else { # If no weights are passed
    if(ncol(y) > 1){ # Multichannel
      weights = matrix(1, nrow = nrow(y), ncol = ncol(y))
    } else { # Single channel
      weights = rep(1, nrow(X))
    }
  }
  
  if(!is.null(offset)){
    if(ncol(y) > 1 & ncol(as.matrix(offset)) == 1){
      offset <- matrix(rep(offset, times = ncol(y)), nrow = length(offset), ncol = ncol(y))
    }
  } else { # If no offset is passed
    if(ncol(y) > 1){ # Multichannel
      offset <- Matrix(0, nrow = nrow(X), ncol=ncol(y), sparse=T)
    } else { # Single channel
      offset = rep(0, nrow(X))
    }
  }
  if(!is.null(start)){
    if(ncol(y) > 1 & ncol(as.matrix(start)) == 1){
      start <- matrix(rep(start, times = ncol(y)), nrow = length(start), ncol = ncol(y))
    } else if(nrow(as.matrix(start)) != ncol(X) | ncol(as.matrix(start)) != ncol(y)){ # nrow(as.matrix(start)) != )
      stop("Start needs to be of dimensions ncol(X), ncol(y).")
    } 
    if(normalize == TRUE){
      Xscaled = normalizeX(X)
      scaling = attr(Xscaled, "scaling")
      start = start*scaling
    }
    start = warm.start*start + (1-warm.start)*rep(sd(start)^(1/gamma), ncol(X))
  } else {
    start <- matrix(0, nrow = ncol(X), ncol = ncol(y)) # passing all zero start values is interpreted as NULL & will result in proper initialisation later on # matrix(sqrt(3), nrow = ncol(X), ncol = ncol(y)) # corresponding in Bayesian sense to assuming initial Gaussian prior with SD of 3, which is the most commonly chosen prior (lambda=1/start^2 and lambda=1/SD prior^2)
  }
  
  if(nonnegative) lower[lower < 0] <- 0
  if(ncol(y) > 1){
    if(!grouped) group.maxit <- 1 # Set group max iteration to 1 if group selection is negative -> channels will be fit independently

    if(parallelize & ncores > 1){
      # Determine the type of cluster based on the operating system
      cluster_type <- if (Sys.info()["sysname"] == "Windows") "PSOCK" else "FORK"
      
      # Create a cluster using the number of cores
      cl <- parallel::makeCluster(ncores, type = cluster_type) # on Mac & Ubuntu type="FORK" is faster
      
      # Export the functions to the workers
      parallel::clusterExport(cl, list("L0glmfun", "L0glm.fit.internal"))
      
      # Load the required libraries and functions on each worker
      # (this only needs to be done once)
      cleval <- parallel::clusterEvalQ(cl, {
        # require(L0glm)
        devtools::load_all(".") # TODO remove this? this should not be necessary once all callable functions are correctly exported by the L0glm package
      })
      
      # Create the chunks for parallelization
      cols <- 1:ncol(y)
      chunk.size <- length(cols) %/% ncores
      extra <- length(cols) %% ncores
      chunk_sizes <- rep(chunk.size, ncores)
      if(extra > 0) {
        chunk_sizes[1:extra] <- chunk_sizes[1:extra] + 1
      }
      
      split_points <- c(0, cumsum(chunk_sizes))
      split_vector <- rep(1:ncores, chunk_sizes)
      chunks <- split(cols, split_vector)
      
      
      for(i in 1L:group.maxit){
        # print(i)
        if(i > 1){
          avgcoefs = Matrix::rowMeans(coefficients)
          variable.indices = which(avgcoefs > 0)
          no.pen.new.indices <- no.pen.start[no.pen.start %in% variable.indices] # Recalculate no.pen indices
          no.pen <- match(no.pen.new.indices, variable.indices)
          start <- matrix(rep(avgcoefs, times = ncol(y)), nrow = length(avgcoefs), ncol = ncol(y))
        } else {
          variable.indices <- 1:ncol(X)
        }
        coefficients <- lambdas.mat <- Matrix(0, nrow = ncol(X), ncol = ncol(y), sparse = TRUE)
        weights.mat <- Matrix(0, nrow = nrow(X), ncol = ncol(y), sparse = TRUE)
        
        if(length(variable.indices) == 0) {
          break
          if(verbose) warning("Zero coefficients produced. Try using a lower lambda.")
        } else if (length(variable.indices) == 1) { # If only 1 variable selected => don't parallelize
          fit <- do.call(cbind, lapply(1:ncol(y), function(col) {
            L0glm.fit.internal(X = X[,variable.indices,drop=F],
                               y = y[,col,drop=F],
                               family=family,
                               lambda = lambda,
                               gamma = gamma,
                               nonnegative = nonnegative,
                               normalize = normalize,
                               offset = as.vector(offset[col, drop = F]),
                               weights = as.vector(weights[,col,drop=F]),
                               warm.start = warm.start,
                               prior.mean = prior.mean[variable.indices],
                               dispersion = dispersion,
                               start = as.vector(start[variable.indices, col, drop = F]),
                               lower = lower[variable.indices],
                               upper = upper[variable.indices],
                               maxit = maxit,
                               thresh = thresh,
                               solver = solver,
                               method = method,
                               tol = tol,
                               no.pen = no.pen,
                               control.solver = control.solver,
                               verbose = verbose)
          }))
          coefficients[variable.indices,] <-  Matrix(do.call(cbind, fit[1, ]), sparse = TRUE)
          weights.mat <-  Matrix(do.call(cbind, fit[2, ]), sparse = TRUE)
          # pen.weights.mat[variable.indices,] <-  Matrix(do.call(cbind, fit[3, ]), sparse = TRUE)
          lambdas.mat[variable.indices,] <-  Matrix(do.call(cbind, fit[3, ]), sparse = TRUE)
          z.mat <-  Matrix(do.call(cbind, fit[4, ]), sparse = TRUE)
          
        } else { # If length(variable.indices) > 1 (more than 1 variable selected)
          fit <- parallel::clusterApply(cl, chunks, L0glmfun,
                              X = X[,variable.indices, drop = F], 
                              y = y, 
                              family = family, 
                              lambda = lambda, 
                              gamma = gamma,
                              maxit = maxit,
                              tol = tol,
                              no.pen = no.pen,
                              thresh = thresh, 
                              solver = solver, 
                              warm.start = warm.start,
                              prior.mean = prior.mean[variable.indices],
                              weights = weights, 
                              dispersion = dispersion,
                              start = start[variable.indices, , drop = F], 
                              lower = lower[variable.indices], 
                              upper = upper[variable.indices], 
                              offset = offset,
                              method = method, 
                              nonnegative = nonnegative, 
                              normalize = normalize,
                              control.solver = control.solver, 
                              verbose = verbose)
          
          # Extract results
          coefficients[variable.indices,] <- do.call(cbind, lapply(fit, function(x) sapply(x, function(y) y$coefficients)))
          weights.mat <- do.call(cbind, lapply(fit, function(x) sapply(x, function(y) y$weights)))
          # pen.weights.mat[variable.indices,] <- do.call(cbind, lapply(fit, function(x) sapply(x, function(y) y$pen.weights)))
          dispersion.mat <- do.call(cbind, lapply(fit, function(x) sapply(x, function(y) y$dispersion)))
          lambdas.mat[variable.indices,] <- do.call(cbind, lapply(fit, function(x) sapply(x, function(y) y$lambdas)))
          z.mat <- do.call(cbind, lapply(fit, function(x) sapply(x, function(y) y$z)))
        }
      }
      
      if(relaxed){
        # Apply relaxed group L0 feature selection
        avgcoefs = rowMeans(coefficients)
        variable.indices = which(avgcoefs != 0)
        no.pen.new.indices <- no.pen.start[no.pen.start %in% variable.indices] # Recalculate no.pen indices
        no.pen <- match(no.pen.new.indices, variable.indices)
        
        if(length(variable.indices) > 1){
          start <- matrix(rep(avgcoefs, times = ncol(y)), nrow = length(avgcoefs), ncol = ncol(y))
          coefficients <- lambdas.mat <- Matrix(0, nrow = ncol(X), ncol = ncol(y), sparse = TRUE)
          weights.mat <- Matrix(0, nrow = nrow(X), ncol = ncol(y), sparse = TRUE)
          
          fit <- parallel::clusterApply(cl, chunks, L0glmfun,
                              X = X[,variable.indices, drop = F], 
                              y = y, 
                              family = family, 
                              lambda = 1e-30, 
                              gamma = gamma,
                              maxit = maxit,
                              tol = tol,
                              no.pen = no.pen,
                              thresh = thresh, 
                              solver = solver, 
                              warm.start = warm.start,
                              prior.mean = prior.mean[variable.indices],
                              weights = weights, 
                              dispersion = dispersion,
                              start = start[variable.indices, , drop = F], 
                              lower = lower[variable.indices], 
                              upper = upper[variable.indices], 
                              offset = offset,
                              method = method, 
                              nonnegative = nonnegative, 
                              normalize = normalize,
                              control.solver = control.solver, 
                              verbose = verbose)
          # Extract results
          coefficients[variable.indices,] <- do.call(cbind, lapply(fit, function(x) sapply(x, function(y) y$coefficients)))
          weights.mat <- do.call(cbind, lapply(fit, function(x) sapply(x, function(y) y$weights)))
          # pen.weights.mat[variable.indices,] <- do.call(cbind, lapply(fit, function(x) sapply(x, function(y) y$pen.weights)))
          dispersion.mat <- do.call(cbind, lapply(fit, function(x) sapply(x, function(y) y$dispersion)))
          lambdas.mat[variable.indices,] <- do.call(cbind, lapply(fit, function(x) sapply(x, function(y) y$lambdas)))
          z.mat <- do.call(cbind, lapply(fit, function(x) sapply(x, function(y) y$z)))
          
        } else if (length(variable.indices) == 1){
          fit <- do.call(cbind, lapply(1:ncol(y), function(col) {
            L0glm.fit.internal(X = X[,variable.indices,drop=F],
                               y = y[,col,drop=F],
                               family=family,
                               lambda = lambda,
                               gamma = gamma,
                               nonnegative = nonnegative,
                               normalize = normalize,
                               offset = as.vector(offset[col, drop = F]),
                               weights = as.vector(weights[,col,drop=F]),
                               warm.start = warm.start,
                               prior.mean = prior.mean[variable.indices],
                               dispersion = dispersion,
                               start = as.vector(start[variable.indices, col, drop = F]),
                               lower = lower[variable.indices],
                               upper = upper[variable.indices],
                               maxit = maxit,
                               thresh = thresh,
                               solver = solver,  
                               method = method,  
                               tol = tol,
                               no.pen = no.pen,
                               control.solver = control.solver,
                               verbose = verbose)
          }))
          coefficients[variable.indices,] <-  Matrix(do.call(cbind, fit[1, ]), sparse = TRUE)
          weights.mat <-  Matrix(do.call(cbind, fit[2, ]), sparse = TRUE)
          # pen.weights.mat[variable.indices,] <-  Matrix(do.call(cbind, fit[3, ]), sparse = TRUE)
          # residuals.mat <- Matrix(do.call(cbind, fit[3, ]), sparse = TRUE)
          dispersion.mat[variable.indices,] <- Matrix(do.call(cbind, fit[3, ]), sparse = FALSE)
          lambdas.mat[variable.indices,] <-  Matrix(do.call(cbind, fit[4, ]), sparse = TRUE)
          z.mat <-  Matrix(do.call(cbind, fit[5, ]), sparse = TRUE)
        } else { # If variable indices == 0
          if(verbose) warning("Relaxed fit wasn't possible due to all coefficients being 0.")
        }
      }
      parallel::stopCluster(cl)
    } else { # If parallelize = FALSE
      for(i in 1L:group.maxit){
        # print(i)
        if(i > 1){
          avgcoefs <- rowMeans(coefficients)
          variable.indices = which(avgcoefs > 0) # Only refit nonzero indices
          no.pen.new.indices <- no.pen.start[no.pen.start %in% variable.indices] # Recalculate no.pen indices
          no.pen <- match(no.pen.new.indices, variable.indices)
          start <- matrix(rep(avgcoefs, times = ncol(y)), nrow = length(avgcoefs), ncol = ncol(y))
        } else {
          variable.indices <- 1:ncol(X)
        }
        
        coefficients <- lambdas.mat <- Matrix(0, nrow = ncol(X), ncol = ncol(y), sparse = TRUE)
        weights.mat <- Matrix(0, nrow = nrow(X), ncol = ncol(y), sparse = TRUE)
        
        if(length(variable.indices) == 0) break
        
        
        fit <- do.call(cbind, lapply(1:ncol(y), function(col) {
          L0glm.fit.internal(X = X[,variable.indices,drop=F],
                             y = y[,col,drop=F],
                             family=family,
                             lambda = lambda,
                             gamma = gamma,
                             nonnegative = nonnegative,
                             normalize = normalize,
                             offset = as.vector(offset[col, drop = F]),
                             weights = as.vector(weights[,col,drop=F]),
                             warm.start = warm.start,
                             prior.mean = prior.mean[variable.indices],
                             dispersion = dispersion,
                             start = as.vector(start[variable.indices, col, drop = F]),
                             lower = lower[variable.indices],
                             upper = upper[variable.indices],
                             maxit = maxit,
                             thresh = thresh,
                             solver = solver,  
                             method = method,  
                             tol = tol,
                             no.pen = no.pen,
                             control.solver = control.solver,
                             verbose = verbose)
        }))
        coefficients[variable.indices,] <-  Matrix(do.call(cbind, fit[1, ]), sparse = TRUE)
        weights.mat <-  Matrix(do.call(cbind, fit[2, ]), sparse = TRUE)
        # pen.weights.mat[variable.indices,] <-  Matrix(do.call(cbind, fit[3, ]), sparse = TRUE)
        dispersion.mat <-  Matrix(do.call(cbind, fit[3, ]), sparse = TRUE)
        lambdas.mat[variable.indices,] <-  Matrix(do.call(cbind, fit[4, ]), sparse = TRUE)
        z.mat <-  Matrix(do.call(cbind, fit[5, ]), sparse = TRUE)
      }
      
      if(relaxed){
        # Apply relaxed group L0 feature selection
        # using 1/ avg coefficients over channels from previous fit as starting values for another L0 GLM fit but now using lambda=0, 
        # i.e. using a nonnegative Poisson GLM, using variables selected earlier via iterated L0 GLMs
        avgcoefs = rowMeans(coefficients)
        variable.indices = which(avgcoefs!=0)
        no.pen.new.indices <- no.pen.start[no.pen.start %in% variable.indices] # Recalculate no.pen indices
        no.pen <- match(no.pen.new.indices, variable.indices)
        
        if(length(variable.indices) != 0){
          start <- matrix(rep(avgcoefs, times = ncol(y)), nrow = length(avgcoefs), ncol = ncol(y))
          coefficients <- lambdas.mat <- Matrix(0, nrow = ncol(X), ncol = ncol(y), sparse = TRUE)
          weights.mat <- Matrix(0, nrow = nrow(X), ncol = ncol(y), sparse = TRUE)
          
          fit <- do.call(cbind, lapply(1:ncol(y), function(col) {
            L0glm.fit.internal(X = X[,variable.indices,drop=F],
                               y = as.vector(y[,col,drop=F]),
                               family = family,
                               lambda = 1e-30,
                               gamma = gamma,
                               nonnegative = nonnegative,
                               normalize = normalize,
                               weights = as.vector(weights[,col,drop=F]),
                               dispersion = dispersion,
                               start = as.vector(start[variable.indices, col, drop = F]),
                               lower = lower[variable.indices],
                               upper = upper[variable.indices],
                               maxit = maxit,
                               warm.start = warm.start,
                               prior.mean = prior.mean[variable.indices],
                               thresh = thresh,
                               solver = solver,
                               method = method,
                               tol = tol,
                               no.pen = no.pen,
                               control.solver = control.solver,
                               verbose = verbose)}))
          
          coefficients[variable.indices,] <-  Matrix(do.call(cbind, fit[1, ]), sparse = TRUE)
          weights.mat <-  Matrix(do.call(cbind, fit[2, ]), sparse = TRUE)
          # pen.weights.mat[variable.indices,] <-  Matrix(do.call(cbind, fit[3, ]), sparse = TRUE)
          dispersion.mat <- Matrix(do.call(cbind, fit[3, ]), sparse = TRUE)
          lambdas.mat[variable.indices,] <-  Matrix(do.call(cbind, fit[4, ]), sparse = TRUE)
          z.mat <-  Matrix(do.call(cbind, fit[5, ]), sparse = TRUE)
        }
      }
    }
    weights.ret <- weights.mat
    # pen.weights <- pen.weights.mat
    lambdas <- lambdas.mat
    z <- z.mat
  } else { # Single variable
    fit <- L0glm.fit.internal(X = X, 
                              y = as.vector(y),
                              weights = weights,
                              family = family,
                              nonnegative = nonnegative,
                              normalize = normalize,
                              offset = offset,
                              dispersion = dispersion,
                              start = start,
                              lower = lower,
                              upper = upper,
                              lambda = lambda, 
                              warm.start = warm.start,
                              prior.mean = prior.mean,
                              gamma = gamma,
                              thresh = thresh,
                              solver = solver, 
                              method = method, 
                              maxit = maxit,
                              tol = tol,
                              no.pen = no.pen,
                              control.solver = control.solver,
                              verbose = verbose)
    
    coefficients <- fit$coefficients
    dispersion <- fit$dispersion
    # pen.weights <- fit$pen.weights
    lambdas <- fit$lambdas
    weights.ret <- fit$weights
    z <- fit$z
    
    if(relaxed){
      # Apply relaxed group L0 feature selection
      # using 1/ avg coefficients over channels from previous fit as starting values for another L0 GLM fit but now using lambda=0, 
      # i.e. using a nonnegative Poisson GLM, using variables selected earlier via iterated L0 GLMs
      variable.indices = which(coefficients!=0)
      no.pen.new.indices <- no.pen.start[no.pen.start %in% variable.indices] # Recalculate no.pen indices
      no.pen <- match(no.pen.new.indices, variable.indices)
      start <- coefficients
      
      fit <- L0glm.fit.internal(X = X[,variable.indices, drop = F],
                                y = y,
                                family = family,
                                weights = weights,
                                nonnegative = nonnegative,
                                normalize = normalize,
                                offset = offset,
                                dispersion = dispersion,
                                start = start[variable.indices],
                                lower = lower[variable.indices],
                                upper = upper[variable.indices],
                                lambda = 1e-30,
                                warm.start = warm.start,
                                prior.mean = prior.mean[variable.indices],
                                gamma = gamma,
                                thresh = thresh,
                                solver = solver,
                                method = method,
                                maxit = maxit,
                                tol = tol,
                                no.pen = no.pen,
                                control.solver = control.solver,
                                verbose = verbose)
      
      coefficients <- lambdas <- rep(0, ncol(X))
      weights.ret <- rep(0, nrow(X))
      
      coefficients[variable.indices] <- fit$coefficients
      # pen.weights[variable.indices] <- fit$pen.weights
      lambdas[variable.indices] <- fit$lambdas
      dispersion <- fit$dispersion
      z <- fit$z
      weights.ret <- fit$weights
    }
    
    
    # TODO: Possible option to add:
    # Thresholding the solution to achieve optimal IC
    # might work well especially with nonnegativity constraints,
    # cf https://proceedings.neurips.cc/paper/2011/hash/d6723e7cd6735df68d1ce4c704c29a04-Abstract.html
    # and would be faster than cross validation & maybe as effective
    # in presence of nonnegativity constraints
    
    # # Option 1: with optimize function. 0.49s for SC p=n=100
    # thresh.post <- optimize(optimize.threshold, interval = c(0, max(coefficients)), coefficients = coefficients, Xmat = X, y = y, IC = IC, family = family, maximum = FALSE)
    # coefficients[abs(coefficients) < thresh.post$minimum] = 0
    
    # # Option 2: with optim package. https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/optim
    # # method "SANN" => VERY slow (136.78s for n=p=100).
    # # method "Nelder-Mead" => unreliable for 1 dimension.
    # thresh.post <- optim(par = thresh,
    #                      fn = optimize.threshold,
    #                      coefficients = coefficients,
    #                      Xmat = X,
    #                      y = y,
    #                      IC = IC,
    #                      family = family,
    #                      method = "Brent"  # "SANN" or "Brent"
    # )
    # # Gives a negative threshold
    # coefficients[abs(coefficients) < thresh.post$par] = 0
    
    # # Option 3: package DEoptim. https://cran.r-project.org/web/packages/DEoptim/DEoptim.pdf. 27.83s
    # thresh.post <- DEoptim(fn = optimize.threshold, lower = thresh, upper = max(coefficients),
    #                        coefficients = coefficients,
    #                        Xmat = X,
    #                        y = y,
    #                        IC = IC,
    #                        family = family)
    # coefficients[abs(coefficients) < thresh.post$optim$bestmem] = 0
    
    # All three threshold optimizations give a threshold that's too high and only leaves 1 nonzero coefficient except for optim method "SANN", which gives a negative threshold
    
  }
  
  weights.ret[is.infinite(weights)] <- 1e11
  # pen.weights[is.infinite(pen.weights)] <- 1e11
  lambdas[is.infinite(lambdas)] <- 1e11
  
  return(list(coefficients = coefficients, weights = weights.ret, 
              # pen.weights = pen.weights, 
              dispersion = dispersion,
              init.lambda = lambda, 
              lambdas = lambdas, family = family, z = z, X = X, y = y))
}


# Internal helper function
#
# Helper function used in the \code{L0glm.fit} function for parallel computation. It takes a subset of column indices cols along with 
# other parameters such as regularization parameters, solver options, and data matrices X and y. It performs L0 regularization using 
# the \code{L0glm.fit.internal} function on each column independently and returns the results as a matrix. The function supports weighting 
# of observations and nonnegative constraints if specified. 
#
L0glmfun <- function(cols, lambda, gamma, thresh, solver, method, nonnegative, normalize, X, y, family, warm.start = 0.9, 
                prior.mean = NULL, weights = NULL, dispersion = 1, start = NULL, lower, upper, verbose = FALSE, maxit, tol, offset, no.pen = 0, 
                control.solver = do.call("control.solver.gen", list(solver = solver, method = method))) { 
  fit <- lapply(cols, function(col) {
    y = as.vector(y[,col,drop=F])
    
    # Call L0glm.fit.internal
    L0glm.fit.internal(X = X,
                       y = y,
                       family = family,
                       lambda = lambda,
                       gamma = gamma,
                       weights = as.vector(weights[,col,drop=F]),
                       dispersion = dispersion,
                       start = as.vector(start[, col, drop = F]),
                       maxit = maxit,
                       tol = tol,
                       no.pen = no.pen,
                       lower = lower,
                       offset = as.vector(offset[,col, drop = F]),
                       upper = upper,
                       prior.mean = prior.mean,
                       thresh = thresh,
                       solver = solver,
                       nonnegative = nonnegative,
                       normalize = normalize,
                       method = method,
                       verbose = verbose)
  })
  
  return(fit)
}