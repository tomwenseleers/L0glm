#' L0-penalized GLM fit using an iterative adaptive ridge regression approach
#' 
#' This function fits an L0-penalized GLM models using an iterative adaptive ridge
#' regression approach under the specified noise model & link function, with or
#' without nonnegativity or box constraints (upper & lower bounds) on the fitted
#' coefficients.
#' 
#' @param X  the covariate (design) matrix (can be dense or sparse, i.e. of class dgCMatrix)
#' @param y  vector of values of the dependent variable
#' @param weights  observation weights
#' @param family  desired noise model and link function; all default R families
#'   are supported, e.g. gaussian(identity), poisson(identity), Gamma(identity),
#'   inv.gaussian(identity) or binomial(logit)
#' @param lam regularization parameter lambda, which will be multiplied with the
#' adaptive penalty weights during each iteration of the L0 adaptive ridge algorithm
#' to give a per-coefficient adaptive ridge penalty
#' @param nonnegative a vector of TRUE/FALSE representing whether
#' to apply nonnegative constraints on each of the model's coefficients
#' @param algo the algorithm used to apply the constraints. Can be 
#' osqp (quadratic programming), clipping or scaling
#' @param miniter the minimum number of iterations of the iterative ridge after
#'   which the observation weights of the IRLS GLM algorithm will be updated
#' @param maxit the maximum number of iterations of the iterative adaptive ridge algorithm
#' @param tol  convergence tolerance
#' @param thresh coefficients with absolute value below this threshold will be
#'   set to zero
#' 
#' @return list with elements beta (L0 adaptive ridge coefficients) and weights (the optimized weights of the IRLS GLM fitting algorithm)
#' 
#' @examples 
# TODO: make proper reproducible example, using simulate_spike_train
#' fit <- l0adri(X = sim$X, y = sim$y, family = poisson(identity),
#'                    algo = "osqp", lam = 0,
#'                    nonnegative = rep(TRUE,nrow(X)))
#' l0araR_wgaus_fit <- l0adri(X = sim$X, y = sim$y, weights = 1/(sim$y+1),
#'                            family = gaussian(identity), lam = 1, 
#'                            nonnegative = TRUE, maxit = 25, tol = 1E-8)
#' 
#' 
#' @export
l0adri = function(X, y, weights = rep( 1, nrow( X)), family = c("gaussian(identity)", "poisson(identity)", "Gamma(identity)", "inv.gaussian(identity)", "binomial(logit)"), 
                  lam = 0, lower=rep(-.Machine$double.xmax,ncol(X)),upper=rep(.Machine$double.xmax,ncol(X)), 
                  algo=c('osqp','linearsolver'), 
                  miniter=5, 
                  maxit=30, tol=1E-8, thresh=1E-3) {
  
  # based on basic GLM algo with weighted LS step replaced by (constrained) weighted adaptive ridge LS regression
  
  if (!is(family, "family")) family = family()
  variance = family$variance
  linkinv = family$linkinv
  mu.eta = family$mu.eta
  etastart = NULL
  
  nobs = nrow(X)    # nobs & nvars are needed by the initialize expression below
  nvars = ncol(X)
  eval(family$initialize) # initializes n and mustart
  eta = family$linkfun(mustart) # we initialize eta based on this
  dev.resids = family$dev.resids #not needed
  dev = sum(dev.resids(y, linkinv(eta), weights)) #not needed
  devold = 0 #not needed
  
  iter = 1
  # ybar = mean(y) # this was initialization used by l0ara, but I use initialization of GLMs as in R's glm code
  # beta = rep(0, nvars)
  # beta[1] = family$linkfun(ybar) # assumes 1st column of X is intercept term
  
  beta = rep(1,nvars)
  while (iter < maxit) {
    # print(iter)
    mu = linkinv(eta)
    varg = variance(mu)
    gprime = mu.eta(eta)
    z = eta + (y - mu) / gprime # potentially -offset
    
    if (iter==1|iter>miniter) W = weights * as.vector(gprime^2 / varg) # TODO: add miniter as argument
    nonzero = (beta!=0)
    beta_old = beta
    
   
    beta = lm_adridge(X, y=z, weights=W, lambda=lam, beta=beta_old, thresh=thresh, 
                      lower=lower, upper=upper, algo=algo) # TODO allow lower and upper instead of nonnegative to allow box constraints, cf https://stats.stackexchange.com/questions/136563/linear-regression-with-individual-constraints-in-r
    
    eta = as.matrix(X %*% beta) # potentially +offset
    
    dev = sum(dev.resids(y, mu, weights)) #notneeded
    iter = iter + 1
    if(iter >= maxit) { warning("Algorithm did not converge. Increase maxit.") }
    if(sum(nonzero)==0){ break}
    if(NaN %in% beta){  break}
    if((sum(beta))==0){break}
    
    if(norm(beta-beta_old, "2") < tol){break}
    
    devold = dev #delete
  }
  
   mu = linkinv(eta)
  return(list(beta=beta,weights=W, family=family, fitted.values=as.vector(mu), lambda = lam)) # weights are needed to calculate weighted mean square error
}
