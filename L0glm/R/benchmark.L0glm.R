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
