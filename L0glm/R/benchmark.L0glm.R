# Functions for benchmarking L0glm algo against other algorithms

#' Fit Sparse Nonnegative Identity Link Generalized Linear Models using Various Methods
#'
#' Fits a sparse generalized linear model to simulated data using different
#' methods specified by the user subject to nonnegativity constraints. 
#' This function utilizes different algorithms
#' and packages to fit sparse nonnegative identity link GLMs, and various aspects of the fitting process can
#' be controlled using the \code{options} parameter.
#'
#' @param sim Simulated data generated using the \code{\link{simul}}
#'   function.
#' @param family An optional \code{family} object specifying the GLM family. Currently (nonnegative) identity link \code{gaussian()} and \code{poisson()} are supported.
#' @param method A numeric value indicating the method to use for fitting the
#'   GLM. The possible values are:
#' \describe{
#'   \item{1}{\strong{Nonnegative L0glm}: Uses the \code{L0glm.fit} function to produce a nonnegative L0-penalized GLM fit.}
#'   \item{2}{\strong{Nonnegative LASSO using glmnet}: Employs the \code{glmnet} package for LASSO regression with nonnegativity constraints.}
#'   \item{3}{\strong{Adaptive nonnegative LASSO using glmnet}: Performs adaptive LASSO with weights from an initial nonnegative LASSO fit.}
#'   \item{4}{\strong{Nonnegative L0Learn}: Uses the \code{L0Learn} package for L0-penalized regression with nonnegativity constraints.}
#'   \item{5}{\strong{Nonnegative L0Learn after glmnet LASSO}: Applies \code{L0Learn} after selecting variables with a LASSO fit from \code{glmnet}.}
#'   \item{6}{\strong{Nonnegative abess}: Utilizes the \code{abess} package for best subset selection with nonnegativity constraints.}
#'   \item{7}{\strong{Nonnegative abess after nonnegative LASSO}: Performs best subset selection with \code{abess} after a nonnegative LASSO.}
#'   \item{8}{\strong{Nonnegative weighted MCP in ordinis}: Fit an MCP penalized GLM model with nonnegativity constraints using the \code{ordinis} package.}
#'   \item{9}{\strong{ordinis SCAD}: Fit a SCAD penalized GLM model with nonnegativity constraints using the \code{ordinis} package.}
#' }
#' @param options A list of additional options specific to the method being
#'   used. These options are passed to the underlying functions from the
#'   respective packages. Refer to the documentation of those packages for
#'   details.
#' @param plot Logical; if \code{TRUE}, plots related to the fitting process are
#'   displayed.
#' @param IC Character string specifying the information criterion to use for
#'   model selection. Possible values are \code{"gic"}, \code{"aic"}, \code{"bic"}, \code{"mbic"} or \code{"ebic"}.
#' @param internal Logical; for internal calls (e.g., when \code{method = 3})
#'   where a starting fit is needed. If \code{TRUE}, the function returns the
#'   full fitted model object instead of just the coefficients.
#'
#' @return If \code{internal = FALSE}, returns a vector or matrix of
#' coefficients from the fitted model. If \code{internal = TRUE}, returns the
#' full fitted model object from the underlying fitting function.
#'
#' @details This function supports various methods for fitting GLMs, including
#' L0 penalization (using \code{L0glm}, \code{abess} or \code{L0Learn}), LASSO or adaptive LASSO (using \code{glmnet}) and MCP or SCAD (using \code{ordinis}). 
#' The \code{options} parameter allows passing additional arguments to these underlying functions.
#'
#' **Method Details:**
#'
#' - **Method 1 (Nonnegative L0glm):** Fits a GLM using L0 penalization via the \code{L0glm.fit} function.
#' - **Method 2 (Nonnegative LASSO using glmnet):** Fits a LASSO regression model using the \code{glmnet} package.
#' - **Method 3 (Adaptive nonnegative LASSO using glmnet):** Performs an adaptive LASSO by weighting coefficients from an initial LASSO fit.
#' - **Method 4 (Nonnegative L0Learn):** Fits a nonnegative L0-penalized regression using \code{L0Learn}.
#' - **Method 5 (Nonnegative L0Learn after glmnet LASSO):** Applies \code{L0Learn} on variables selected by a previous LASSO fit.
#' - **Method 6 (Nonnegative abess):** Uses \code{abess} for best subset selection with nonnegativity constraints.
#' - **Method 7 (Nonnegative abess after nonnegative LASSO):** Performs best subset selection after a nonnegative LASSO.
#' - **Method 8 (Nonnegative MCP):** Fits a model using the MCP penalty with nonnegativity constraints via \code{ordinis}.
#' - **Method 9 (Nonnegative SCAD):** Fits a model using the SCAD penalty with nonnegativity constraints via via \code{ordinis}.
#'
#' @seealso \code{\link{simul}}, \code{\link[glmnet]{glmnet}},
#' \pkg{L0Learn}, \code{\link[abess]{abess}}, \pkg{ordinis}
#'
#' @example examples/fit.GLM.R
#' @export
fit.GLM <- function(sim, family = NULL, method = 1, options = NULL, plot = TRUE, 
                    IC = "gic", internal = FALSE){
  X <- sim$X
  y <- sim$y
  p <- ncol(X)
  n <- nrow(X)
  if(family$family == "gaussian"){
    W <- rep(1,ncol(X))
  } else if(family$family == "poisson"){
    W <- as.vector(1/(y+0.1))
  }
  
  if(method == 1) { # L0glm
    if(is.null(family)){
      fam.fun <- get(sim$family)
      family <- fam.fun(link = "identity")
    }
    
    if(is.null(options)){
      result <- L0glm.fit(X, y, 
                          family = family, 
                          lambda = 1, 
                          prior.mean = NULL,
                          nonnegative = TRUE,
                          solver = "eigen", 
                          method = 1L)
    } else {
      args <- c(list(X=X, y=y, family=family), options)
      result <- do.call(L0glm.fit, args)
    }
    if(plot){
      plot(x = result$coefficients, y=sim$beta_true, pch=16, col='blue',
           main=paste0("L0glm identity link Poisson with lambda preset to optimise ", toupper(IC)),
           xlab="estimated coefficients", ylab="true coefficients")
    }
    
    coefs <- result$coefficients
    
    
  } else if(method == 2) { # glmnet
    if(is.null(options)){
      result <- glmnet(x = X, 
                       y = y, 
                       weights = W,
                       family = "gaussian", # poisson(identity) should be supported but gives an error
                       nlambda = 20L,
                       intercept = F, 
                       standardize = F, 
                       alpha = 1, 
                       lower.limits = 0)
    } else {
      args <- c(list(x = X, y = y, family = "gaussian"), options)
      result <- do.call(glmnet, args)
    }
    LLs <- logLik(result, weights = W, X = X, y = y)
    edfs = attr(LLs, "df")
    IC.val = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,IC]
    
    bestlambda_lasso_glmnet = result$lambda[which.min(IC.val)]
    coefs = coef(result, s = bestlambda_lasso_glmnet)
    
    if(plot){
      plot(result, xvar='lambda',
           main = "Regularisation path glmnet nonnegative weighted LASSO")
      abline(h=as.vector(sim$beta_true[sim$beta_true!=0]), col=alpha(1:6,0.2))
      
      plot(x = coefs[-1], y=sim$beta_true, pch=16, col='blue',
           main="glmnet nonnegative weighted LASSO",
           xlab="estimated coefficients", ylab="true coefficients")
    }
    
  } else if(method == 3) { # adaptive LASSO using glmnet 
    # since glmnet is so much faster than ordinis & supports sparse
    # covariate matrices, let's use that to fit an adaptive LASSO model, cf
    # https://rpubs.com/kaz_yos/alasso
    # note that for the first step we use LASSO rather than ridge coefficients, i.e. we use
    # glmnet.fit and then in a second step use a glmnet LASSO fit with
    # adaptive penalty weights penalty.factor = 1 / (abs(glmnet.fit))+delta)
    # with delta=1E-5
    # NOTE: If this would be iterated further this would be a way to approximate L0 norm
    # penalized regression (cf. Nuel 2017 talk)
    
    # Start with regular LASSO fit
    glmnet.fit <- fit.GLM(sim = sim, method = 2, family = family, plot = F, internal = TRUE)
    
    LLs <- logLik(glmnet.fit, weights=W, X=X, y=y)
    edfs = attr(LLs, "df")
    IC.val = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,IC]
    bestlambda_lasso_glmnet = glmnet.fit$lambda[which.min(IC.val)]
    coefs.glmnet.fit = coef(glmnet.fit, s = bestlambda_lasso_glmnet)
    
    # in second step we do another LASSO fit with adaptive penalty weights
    delta = 1E-5
    if(is.null(options)){
      result <- glmnet(x = X, 
                       y = y, 
                       weights = W,
                       family = "gaussian", # poisson(identity) should be supported but gave an error
                       nlambda = 20L,
                       penalty.factor = 1 / (abs(coefs.glmnet.fit[-1])+delta), # adaptive penalty weights
                       intercept = F, 
                       standardize = F, 
                       alpha = 1, # lasso (lower value, ie elastic net, made solution worse)
                       lower.limits = 0) # for nonnegativity constraints
    } else {
      args <- c(list(x = X, y = y, weights = W, family = "gaussian", penalty.factor = 1/(abs(coefs.glmnet.fit[-1])+delta)), options)
      do.call(glmnet, args)
    }
    LLs <- logLik(result, weights=W, X=X, y=y)
    edfs = attr(LLs, "df")
    IC.val = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,IC]
    bestlambda_alasso_glmnet = result$lambda[which.min(IC.val)]
    coefs = coef(result, s = bestlambda_lasso_glmnet)
    
    if(plot){
      plot(result, xvar='lambda',
           main = "Regularisation path glmnet nonnegative weighted adaptive LASSO")
      abline(h=as.vector(sim$beta_true[sim$beta_true!=0]), col=alpha(1:6,0.2))
      
      plot(x=coefs[-1], y=sim$beta_true, pch=16, col='blue',
           main="True vs. estimated coefficients \n glmnet nonnegative weighted adaptive LASSO",
           xlab="estimated coefficients", ylab="true coefficients")
    }
    
  } else if(method == 4) { # nonnegative L0 learn
    if(sim$family == "poisson"){
      X = Diagonal(x=sqrt(W)) %*% X
      y = as.vector(y*sqrt(W))
    } else {
      y = as.vector(y)
    }
    if(is.null(options)){
      result <- L0Learn.fit(x = X,  
                            y = y,
                            loss = "SquaredError",
                            penalty = "L0", algorithm="CD",
                            maxSuppSize = round(p/2),
                            intercept = FALSE,
                            nLambda = 100, # NOTE: to make timings fair we should probably use 100 or 20 lambda values everywhere?
                            autoLambda=TRUE,
                            lows=0) # nonnegativity constraints
    } else {
      args = c(list(x = X, y = y), options)
      result <- do.call(L0Learn.fit, args)
    }
    
    LLs = logLik(result, weights = W, X, y)
    edfs = attr(LLs, "df")
    LLs = as.vector(LLs)
    IC.val = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,IC]
    coefs = coef(result, lambda = result$lambda[[1]][which.min(IC.val)], gamma = 0)
    
    if(plot){
      plot(result)
      
      plot(x=coefs, y=sim$beta_true, pch=16, col='blue',
           main="True vs. estimated coefficients\n nnL0Learn weighted least squares best subset",
           xlab="estimated coefficients", ylab="true coefficients")
    }
    
  } else if(method == 5) { # nonnegative L0Learn after glmnet LASSO
    # Start with regular LASSO fit
    glmnet.fit <- fit.GLM(sim = sim, method = 2, family = family, plot = F, internal = TRUE)
    LLs <- logLik(glmnet.fit, weights = W, X = X, y = y)
    edfs = attr(LLs, "df")
    IC.val = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,IC]
    
    bestlambda = glmnet.fit$lambda[which.min(IC.val)]
    glmnet.coefs = coef(glmnet.fit, s = bestlambda)
    
    selcoefs = which(glmnet.coefs[-1]!=0)
    p_selcoefs = length(selcoefs)
    
    if(sim$family == "poisson"){
      X = as.matrix(Diagonal(x=sqrt(W)) %*% X)
      y = as.vector(y*sqrt(W))
    } else {
      X = as.matrix(X)
      y = as.vector(y)
    }
    
    if(is.null(options)){
      fit <- L0Learn.fit(x = X[,selcoefs,drop=F],  # note: no support for sparse matrices
                         y = y,
                         loss = "SquaredError",
                         penalty = "L0", algorithm = "CD",
                         maxSuppSize = (p_selcoefs-1), 
                         intercept = FALSE, 
                         nLambda = 100,
                         autoLambda=TRUE)
      
    } else {
      args <- c(list(x = X[,selcoefs,drop=F], y = y, maxSuppSize = (p_selcoefs-1)), options)
      fit <- do.call(L0Learn.fit, args)
    }
    if(plot) plot(fit)
    
    LLs = logLik(fit, weights=W, X[,selcoefs,drop=F], y)
    edfs = attr(LLs, "df")
    LLs = as.vector(LLs)
    IC.val = get.IC(loglik = LLs, n = n, p = p_selcoefs, edf = edfs)[,IC]
    coefs.sel = coef(fit,
                     lambda = fit$lambda[[1]][which.min(IC.val)],
                     gamma = 0)
    coefs = Matrix(0, p)
    coefs[selcoefs] = coefs.sel
    if(plot){
      plot(fit)
      plot(x=coefs, y=sim$beta_true, pch=16, col='blue',
           main="True vs. estimated coefficients\nweighted nnL0Learn after nonnegative LASSO",
           xlab="estimated coefficients", ylab="true coefficients")
    } 
    
  } else if(method == 6) { # abess
    # fast for small problems, but becomes slow for n=p=10000
    # if(IC != "ebic") warning("Ebic is recommended with method == 6 for better performance")
    if(is.null(options)){
      fit <- abess(x = X, 
                   y = y,
                   weight = W, 
                   family = "gaussian", # identity link poisson not supported
                   tune.path = "sequence",
                   support.size = c(1:round(p*0.5)),
                   fit.intercept = FALSE,
                   beta.low = 0, # nonnegativity constraints
                   lambda = 0,
                   warm.start = TRUE,
                   tune.type = "bic", # trick to support any IC that L0glm supports
                   ic.scale = preset.lambda(IC = IC, # options aic, bic, gic, ebic, mbic 
                                            p = p, 
                                            n = n)/log(n)) # trick to support any IC that L0glm supports
    } else {
      args <- c(list(x = X, y = y, weight = W, 
                     family = "gaussian", # identity link poisson not supported
                     tune.type = "bic", # trick to support any IC that L0glm supports
                     ic.scale = preset.lambda(IC = IC, # options aic, bic, gic, ebic, mbic 
                                              p = p, 
                                              n = n)/log(n) ), 
                options)
      fit <- do.call(abess, args)
    }
    
    # coefs = as.matrix(extract(fit)$beta) # coefficient matrix for best subset
    coefs = coef(fit, 
                 support.size = fit[["best.size"]], 
                 sparse = FALSE)[-1,] # coefficient matrix for best subset
    
    if(plot){
      plot(fit, type="tune")
      
      plot(x=coefs, y=sim$beta_true, pch=16, col='blue',
           main="True vs. estimated coefficients abess best\nsubset nonnegative weighted least squares",
           xlab="estimated coefficients", ylab="true coefficients")
    }
    
  } else if(method == 7) { # nonnegative abess after nonnegative LASSO
    # very good solution, and fast for small problems, but
    # becomes very slow for large problems, e.g. n=p=10000
    glmnet.fit <- fit.GLM(sim = sim, method = 2, family = family, plot = F, internal = TRUE)
    LLs <- logLik(glmnet.fit, weights = W, X = X, y = y)
    edfs = attr(LLs, "df")
    IC.val = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,IC]
    
    bestlambda = glmnet.fit$lambda[which.min(IC.val)]
    glmnet.coefs = coef(glmnet.fit, s = bestlambda)
    
    selcoefs = which(glmnet.coefs[-1]!=0)
    p_selcoefs = length(selcoefs)
    
    if(is.null(options)){
      fit <- abess(x = X[,selcoefs,drop=F], 
                   y = y, 
                   weight = W,
                   family = "gaussian", # identity link poisson not supported
                   tune.path = "sequence", 
                   support.size = c(1:(p_selcoefs-1)),
                   fit.intercept = FALSE,
                   beta.low = 0, # nonnegativity constraints
                   lambda = 0, 
                   warm.start = TRUE, # FALSE worked better here, but much slower
                   tune.type = "bic", # trick to support any IC that L0glm supports
                   ic.scale = preset.lambda(IC = IC, # options aic, bic, gic, ebic, mbic 
                                            p = p, 
                                            n = n)/log(n)) # trick to support any IC that L0glm supports
    } else {
      args <- c(list(x = X[,selcoefs,drop=F], y = y, 
                     weight = W,
                     family = "gaussian",
                     tune.type = "bic", # trick to support any IC that L0glm supports
                     ic.scale = preset.lambda(IC = IC, # options aic, bic, gic, ebic, mbic 
                                              p = p, 
                                              n = n)/log(n),
                     support.size = c(1:(p_selcoefs-1))), options)
      fit <- do.call(abess, args)
    }
    
    plot(fit, type="tune")  # optimal support size 
    # coefs.selcoefs = as.matrix(extract(fit)$beta) # coefficient matrix for best subset
    coefs.selcoefs = coef(fit, 
                 support.size = fit[["best.size"]], 
                 sparse = FALSE)[-1,] # coefficient matrix for best subset
    coefs = Matrix(rep(0, p), sparse=T)
    coefs[selcoefs] = coefs.selcoefs
    
    if(plot){
      plot(fit, type="tune")  # optimal support size 
      
      plot(x=coefs, y=sim$beta_true, pch=16, col='blue',
           main="True vs estimated coefficients abess best subset\nweighted least squares after nonnegative LASSO",
           xlab="estimated coefficients", ylab="true coefficients")
    }
  } else if(method == 8) { # nonnegative weighted MCP penalized regression in ordinis 
    require_ordinis()
    if(is.null(options)){
      fit <- ordinis_fit(X, y, 
                     weights = W,
                     nlambda = 100L, 
                     penalty = "mcp", # using MCP penalty
                     intercept = FALSE,
                     lower.limits = rep(0, p), # force all coefficients to be positive
                     standardize = FALSE,
                     alpha = 1) 
    } else {
      args <- c(list(X = X, y = y, weights = W), options)
      fit <- do.call(ordinis_fit, args)
    }
    
    LLs = logLik2(fit, weights= as.vector(W), X, as.vector(y))
    edfs = attr(LLs, "df")
    LLs = as.vector(LLs)
    IC.val = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,IC]
    coefs = fit$beta[,which.min(IC.val[-1])+1]
    
    if(plot){
      plot(fit)
      abline(h=as.vector(sim$beta_true[sim$beta_true!=0]), col=alpha(1:6,0.2))
      
      plot(x=coefs[-1], y=sim$beta_true, pch=16, col='blue',
           main="True vs. estimated coefficients\nordinis nonnegative MCP",
           xlab="estimated coefficients", ylab="true coefficients")
    }
  } else if(method == 9) { # ordinis scad
    require_ordinis()
    if(is.null(options)){
      fit <- ordinis_fit(X, y, 
                     weights = W,
                     nlambda = 100L, 
                     penalty = "scad", # using scad penalty
                     intercept = FALSE,
                     lower.limits = rep(0, p), # force all coefficients to be positive
                     standardize = FALSE,
                     alpha = 1) 
    } else {
      args <- c(list(X = X, y = y, weights = W), options)
      fit <- do.call(ordinis_fit, args)
    }
    
    LLs = logLik2(fit, weights=W, X, y)
    edfs = attr(LLs, "df")
    LLs = as.vector(LLs)
    IC.val = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,IC]
    coefs = fit$beta[,which.min(IC.val[-1])+1]
    
    if(plot){
      plot(fit)
      abline(h=as.vector(sim$beta_true[sim$beta_true!=0]), col=alpha(1:6,0.2))
      
      plot(x=coefs[-1], y=sim$beta_true, pch=16, col='blue',
           main="True vs. estimated coefficients\nordinis nonnegative SCAD",
           xlab="estimated coefficients", ylab="true coefficients")
    }
  }
  if(internal) return(result)
  return(coefs)
}





#' Fit Sparse Multichannel Nonnegative Gaussian or Identity Link Poisson Generalized Linear Models using Various Methods
#'
#' Fits a multichannel nonnegative identity link generalized linear model (GLM) to simulated data using different methods specified by the user. 
#' This function supports fitting multivariate response models, utilizing various algorithms and packages. The fitting process can be customized using the \code{options} parameter.
#'
#' @param sim Simulated data generated using the \code{\link{simul}} function. It should be a list containing at least \code{X} (design matrix) and \code{y} (response matrix).
#' @param family An optional \code{family} object specifying the GLM family. Currently supports (nonnegative) identity link \code{gaussian()} and \code{poisson()}. Default is \code{NULL}, which tries to infer from \code{sim}.
#' @param method A numeric value indicating the method to use for fitting the multichannel GLM. The possible values are:
#' \describe{
#'   \item{1}{\strong{abess Multitask Learning}: Uses the \code{\link[abess]{abess}} package for best subset selection with multivariate Gaussian responses.}
#'   \item{2}{\strong{glmnet Multivariate Group LASSO}: Employs the \code{\link[glmnet]{glmnet}} package for multivariate LASSO regression with nonnegativity constraints.}
#' }
#' @param options A list of additional options specific to the method being used. These options are passed to the underlying functions from the respective packages. Refer to the documentation of those packages for details.
#' @param plot Logical; if \code{TRUE}, plots related to the fitting process are displayed. Default is \code{TRUE}.
#' @param IC Character string specifying the information criterion to use for model selection. Possible values are \code{"gic"} (default), \code{"aic"}, \code{"bic"}, \code{"ebic"} and \code{"mbic"}.
#' @param internal Logical; for internal calls where a starting fit is needed. If \code{TRUE}, the function returns the full fitted model object instead of just the coefficients. Default is \code{FALSE}.
#'
#' @return A matrix of coefficients from the fitted model. Each column corresponds to a response variable (channel).
#'
#' @details
#' This function supports fitting multichannel GLMs using different methods. The methods utilize different packages such as \code{\link[abess]{abess}} and \code{\link[glmnet]{glmnet}}. The \code{options} parameter allows passing additional arguments to these underlying functions.
#'
#' **Method Details:**
#'
#' - **Method 1 (abess Best Subset of Group Selection):**
#'   - Uses the \code{\link[abess]{abess}} package to perform best subset of group selection for multivariate Gaussian responses (\code{family = "mgaussian"}) subject to nonnegativity constraints.
#'   - Nonnegativity constraints are applied via the \code{beta.low} parameter.
#'
#' - **Method 2 (glmnet Multivariate Group LASSO):**
#'   - Utilizes the \code{\link[glmnet]{glmnet}} package to fit a multivariate group LASSO regression model with nonnegativity constraints.
#'   - Performs cross-validation to select the best lambda value.
#'
#' @seealso
#' \code{\link{simulate}},
#' \code{\link[abess]{abess}},
#' \code{\link[glmnet]{glmnet}}
#'
#' @example examples/fit.GLM.mc.R
#' @export

# TODO: method 1 should be a group L0glm fit, method 2 abess best subset of groups
# method 3 glmnet multivariate group LASSO
# TODO: add option relaxed=T/F which if T should (for all methods) after initial multitask feature selection
# refit the model unregularised for the selected variables only

fit.GLM.mc <- function(sim, family = NULL, method = 1, options = NULL, plot = TRUE, IC = "gic", internal = FALSE){
  X <- sim$X
  y <- sim$y
  p <- ncol(X)
  n <- nrow(X)
  m <- ncol(sim$y)
  
  if (family$family == "poisson") { 
    W <- 1/(y+0.1) } else {  # approx 1/variance weights to approximate identity link poisson
      W <- NULL } 
  if(method == 1){ # best subset multitask learning using family="mgaussian" using abess
    # TODO check solution quality & timings of multiple response Poisson case using
    # https://skscope.readthedocs.io/en/latest/gallery/GeneralizedLinearModels/poisson-identity-link.html
    # (can be called from R via reticulate package)
    # for tuning sparsity based on IC see
    # https://skscope.readthedocs.io/en/latest/gallery/SparsityLevelSelection/supported-information-criteria-and-cross-validation.html
    if (is.null(options)) {
      fit <- abess(x = X, 
                   y = y, 
                   # weight = W,  # right now not allowed, but I raised a ticket about this on github: https://github.com/abess-team/abess/issues/510 
                   family = "mgaussian",
                   tune.path = "sequence", 
                   support.size = c(1:round(p*0.5)),
                   fit.intercept = FALSE,
                   beta.low = 0, # nonnegativity constraints
                   lambda=0, 
                   warm.start = TRUE,
                   tune.type = "aic" #, # trick to support any IC that L0glm supports
                   # ic.scale = preset.lambda(IC = IC, # options aic, bic, gic, ebic, mbic 
                   #                         p = p, 
                   #                         n = n)/log(n)
                   ) 
    } else {
      fit <- abess(x = X, 
                   y = y, 
                   # weight = W, # right now not allowed, but I raised a ticket about this on github: https://github.com/abess-team/abess/issues/510 
                   family = "mgaussian",
                   tune.path = options$tune.path, 
                   support.size = c(1:round(p*0.5)),
                   fit.intercept = FALSE,
                   beta.low = 0, # nonnegativity constraints
                   lambda=options$lambda, 
                   warm.start = options$warm.start,
                   tune.type = "bic", # trick to support any IC that L0glm supports
                   ic.scale = preset.lambda(IC = IC, # options aic, bic, gic, ebic, mbic 
                                             p = p, 
                                             n = n)/log(n)
                   ) 
    }
    
    if (plot) {
      plot(fit, type="tune")  # optimal support size 
      # coefs = as.matrix(abess::extract(fit)$beta) # coefficient matrix for best subset
      coefs = coef(fit, 
                   support.size = fit[["best.size"]], 
                   sparse = FALSE)[[1]][-1,] # coefficient matrix for best subset
      # coefs[abs(coefs)<1E-5] = 0 # apply same thresholding as in L0glm
      image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=coefs^0.1, col = topo.colors(255),
            useRaster=TRUE,
            xlab="Time", ylab="Channel", main="abess multitask weighted mgaussian\n(red = true peaks, cyan = estimated peaks)")
      abline(v=(1:nrow(sim$X))[as.vector(qlcMatrix::rowMax(sim$beta_true)!=0)], col="red")
      abline(v=(1:nrow(sim$X))[as.vector(qlcMatrix::rowMax(abs(coefs))!=0)], col="cyan")
      
      R2 = round(100*summary(lm(as.vector(rowMeans(sim$beta_true) ~ 0+as.vector(rowMeans(coefs)))))$r.squared,2) # Calculate R2
      title = paste("True vs estimated coefficients\nabess multitask mgaussian\nR2=", R2)
      plot(as.vector(rowMeans(coefs)), as.vector(rowMeans(sim$beta_true)), pch=16, 
           xlab="estimated coefficients",
           ylab="true coefficients", col="steelblue",
           main=title)
    }
    
  } else if(method ==2){ # glmnet nonnegative mgaussian multivariate group LASSO fit. Currently only support nonnegativity unweighted identity link Gaussian. Normally it should also allow using family=poisson(identity) but in practice that one returns an error.
    if(is.null(options)){
      cvfit <- cv.glmnet(x = X, 
                         y = y, 
                         # weights = W, # passing matrix of weights not currently supported
                         alpha = 1,  # lasso
                         family = "mgaussian", 
                         nlambda = 100,
                         nfolds = 5, 
                         standardize = FALSE,
                         standardize.response = FALSE,
                         intercept = FALSE, 
                         relax = FALSE,
                         lower.limits = rep(0, ncol(X)+1)) # for nonnegativity constraints
      
      # Fit the nonnegative group LASSO glmnet model over a sequence of lambda values
      fit <- glmnet(x = X, 
                    y = y,
                    # weights = W, # passing matrix of weights not currently supported
                    alpha = 1, # lasso
                    family = "mgaussian", 
                    nlambda = 100,
                    standardize = FALSE,
                    standardize.response = FALSE,
                    intercept = FALSE, 
                    relax = FALSE, 
                    lower.limits = rep(0, ncol(X)+1))
    } else {
      cvfit <- cv.glmnet(x = X, 
                         y = y, 
                         # weights = W, # passing matrix of weights not currently supported
                         alpha = 1,  # lasso
                         family = "mgaussian", 
                         nlambda = options$nlambda,
                         nfolds = 5, 
                         standardize = options$standardize,
                         standardize.response = options$standardize.response,
                         intercept = options$intercept, 
                         relax = options$relax,
                         lower.limits = rep(0, ncol(X)+1)) # for nonnegativity constraints
      
      # Fit the nonnegative group LASSO glmnet model over a sequence of lambda values
      fit <- glmnet(x = X, 
                    y = y,
                    alpha = 1, # lasso
                    # weights = W, # passing matrix of weights not currently supported
                    family = "mgaussian", 
                    nlambda = options$nlambda,
                    standardize = options$standardize,
                    standardize.response = options$standardize.response,
                    intercept = options$intercept, 
                    relax = options$relax, 
                    lower.limits = rep(0, ncol(X)+1))
    }
    best_lambda <- cvfit$lambda.1se
    # get the coefficients for each task for the best lambda value
    coefs <- coef(fit, s = best_lambda)
    beta <- do.call(cbind, lapply(seq_len(m), function (channel) as.matrix(coefs[[channel]][-1,,drop=F])))
    
    colnames(beta) = 1:ncol(beta)
    # beta[abs(beta)<1E-5] = 0 # slight amount of thresholding of small coefficients using same thresh as used in L0glm
    
    coefs <- beta
    if(plot){
      # plot(cvfit) # TODO remove?
      # plot(fit, xvar="lambda")
      
      # image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=coefs^0.1, col = topo.colors(255),
      #      useRaster=TRUE,
      #      xlab="Time", ylab="Channel", main="nonnegative group Lasso glmnet\n(red = true peaks, cyan = estimated peaks)")
      # abline(v=(1:nrow(sim$X))[as.vector(qlcMatrix::rowMax(sim$beta_true)!=0)], col="red")
      # abline(v=(1:nrow(sim$X))[as.vector(qlcMatrix::rowMax(abs(coefs))!=0)], col="cyan")
      
      R2 = round(100*summary(lm(as.vector(rowMeans(sim$beta_true) ~ 0+as.vector(rowMeans(coefs)))))$r.squared,2) # Calculate R2
      title = paste("True vs estimated coefficients\nnonnegative group Lasso glmnet\nR2=", R2)
      plot(as.vector(rowMeans(coefs)), as.vector(rowMeans(sim$beta_true)), pch=16, 
           xlab="estimated coefficients",
           ylab="true coefficients", col="steelblue",
           main=title)
      
    }
    
    
  }
  
  
  return(coefs)
}



# UTILITY FUNCTIONS ####

# logLik method for gaussian glmnet
# cf logLik.ordinis in https://github.com/jaredhuling/ordinis/blob/master/R/ordinis_methods.R

#' @title Compute Log-Likelihood for glmnet Objects
#' @description Calculates the log-likelihood for models fitted using `glmnet`.
#' @param object A `glmnet` model object.
#' @param ... Additional arguments. Named arguments \code{weights}, \code{X},
#'   and \code{y} are required.
#' @return The log-likelihood of the model.
#' @export
logLik.glmnet = function(object, ...) {
  args <- get_logLik_args(...)
  weights <- args$weights
  X <- args$X
  y <- args$y
  n  = as.numeric(object$nobs)
  df = Matrix::colSums(object$beta!=0)
  # df = nr of nonzero coefficients (incl intercept) (NOTE: more correct for the LASSO would be to calculate this as the trace of the hat matrix)

  # Compute y_pred for all columns of object$beta at once
  y_pred = X %*% object$beta

  # Compute squared residuals, multiply by weights, and sum over rows to get resid.ss
  resid.ss = Matrix::colSums(weights * (matrix(rep(y, ncol(y_pred)), ncol = ncol(y_pred) ) - y_pred) ^ 2)

  rdf = n - df
  logL = -0.5 * n * (log(2 * pi) - log(rdf) + log(resid.ss)) - 0.5 * rdf
  attr(logL, "df") = df  # nr of nonzero coefficients
  attr(logL, "nobs") = n
  attr(logL, "p") = object$dim[[2]]
  return(logL)
}



# log-likelihood method for gaussian L0Learn fit,
# cf logLik.ordinis in https://github.com/jaredhuling/ordinis/blob/master/R/ordinis_methods.R
#' @title Compute Log-Likelihood for L0Learn Objects
#' @description Calculates the log-likelihood for models fitted using `L0Learn`.
#' @param object An `L0Learn` model object.
#' @param ... Additional arguments. Named arguments \code{weights}, \code{X},
#'   and \code{y} are required.
#' @return The log-likelihood of the model.
#' @export
logLik.L0Learn = function(object, ...) {
  args <- get_logLik_args(...)
  weights <- args$weights
  X <- args$X
  y <- args$y
  n  = object$n
  p = object$p
  df = object$suppSize[[1]]
  # df = nr of nonzero coefficients (incl intercept)
  
  # Compute y_pred for all columns of object$beta at once
  y_pred = X %*% object$beta[[1]]
  
  # Compute squared residuals, multiply by weights, and sum over rows to get resid.ss
  resid.ss = colSums(weights * (matrix(rep(y, ncol(y_pred)), ncol = ncol(y_pred) ) -
                                  y_pred) ^ 2)
  
  rdf = n - df
  logL = -0.5 * n * (log(2 * pi) - log(rdf) + log(resid.ss)) - 0.5 * rdf
  attr(logL, "df") = df  # nr of nonzero coefficients
  attr(logL, "nobs") = n
  attr(logL, "p") = object$dim[[2]]
  return(logL)
}

get_logLik_args <- function(...) {
  args <- list(...)
  weights <- args$weights
  X <- args$X
  y <- args$y
  matched <- match(c("weights", "X", "y"), names(args), nomatch = 0L)
  positional <- args[setdiff(seq_along(args), matched)]
  next_pos <- 1L
  if (is.null(weights) && length(positional) >= next_pos) {
    weights <- positional[[next_pos]]
    next_pos <- next_pos + 1L
  }
  if (is.null(X) && length(positional) >= next_pos) {
    X <- positional[[next_pos]]
    next_pos <- next_pos + 1L
  }
  if (is.null(y) && length(positional) >= next_pos) {
    y <- positional[[next_pos]]
  }
  if (is.null(weights) || is.null(X) || is.null(y)) {
    stop("logLik methods require weights, X, and y.", call. = FALSE)
  }
  list(weights = weights, X = X, y = y)
}

require_ordinis <- function() {
  pkg <- paste0("ordi", "nis")
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package 'ordinis' is required for this method. Install it to use MCP or SCAD benchmark fits.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

ordinis_fit <- function(...) {
  require_ordinis()
  getExportedValue(paste0("ordi", "nis"), paste0("ordi", "nis"))(...)
}



# log-likelihood method for gaussian ordinis fit, 
# cf logLik.ordinis in https://github.com/jaredhuling/ordinis/blob/master/R/ordinis_methods.R
# but correctly taking into account weights
# NOTE: ordinis is currently not taking into account the weights in
# the calculation of the log likelihood as loss=residual SS do not take into 
# account weights here https://github.com/jaredhuling/ordinis/blob/master/R/ordinis_call.R
# hence the redefinition of the log likelihood here
logLik2 = function(object, weights, X, y) {  
  n  = as.numeric(object$nobs)
  p = object$nvars
  df = Matrix::colSums(object$beta!=0)  
  # df = nr of nonzero coefficients (incl intercept)
  
  # Compute y_pred for all columns of object$beta at once
  y_pred =  cbind(1, X) %*% object$beta
  
  # Compute squared residuals, multiply by weights, and sum over rows to get resid.ss
  resid.ss = colSums(weights * (matrix(rep(y, ncol(y_pred)), ncol = ncol(y_pred) ) - 
                                  y_pred) ^ 2)
  
  rdf = n - df 
  logL = -0.5 * n * (log(2 * pi) - log(rdf) + log(resid.ss)) - 0.5 * rdf
  attr(logL, "df") = df  # nr of nonzero coefficients
  attr(logL, "nobs") = n
  attr(logL, "p") = p
  return(logL)
}
