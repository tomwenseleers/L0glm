#' Inference on L0glm Models
#'
#' @import boot
#'
#' @description
#'
#' Performs inference on the estimated coefficients of an \code{L0glm} object.
#'
#' @param fit
#' an \code{L0glm} object (see \code{\link{L0glm}}).
#' @param level
#' the confidence level of the confidence interval. \code{level} must be in
#' \eqn{[0,1]} (e.g. \code{level = 0.95} means a 95\% confidence interval).
#' @param boot.repl
#' the number of bootstrap replicates. Only used if non parametric bootstrapping
#' is performed (see Details).
#' @param control.l0
#' a list of parameters controling the L0 penalty loop (see \code{\link{control.l0.gen}}
#' and \code{\link{L0glm}}). This should be the same parameters as used
#' to compute \code{fit}.
#' @param control.iwls
#' a list of parameters controling the L0 penalty loop (see \code{\link{control.iwls.gen}}
#' and \code{\link{L0glm}}).This should be the same parameters as used
#' to compute \code{fit}.
#' @param control.fit
#' a list of parameters controling the L0 penalty loop (see \code{\link{control.fit.gen}}
#' and \code{\link{L0glm}}).This should be the same parameters as used
#' to compute \code{fit}.
#' @param verbose
#' print algorithm progression to console ?
#'
#' @details
#'
#' The inference is performed using two distinct methods according to the
#' constraint and the type of penalty used for fitting the \code{L0glm} object.
#' If no constraint and no penalty or ridge penalty were applied, the inference
#' is performed based on the Fischer information matrix. Otherwise, inference
#' is performed using non parametric percentile bootstrapping.
#'
#' In the case of unconstrained fit, the test hypotheses are:
#' \itemize{
#' \item H0: coef = 0
#' \item H1: coef != 0
#' }
#' In the case of nonnegativity constraint, the test hypotheses are:
#' \itemize{
#' \item H0: coef = 0
#' \item H1: coef > 0
#' }
#'
#' \strong{Fischer information matrix}
#'
#' The returned confidence intervals are Wald confidence intervals.
#'
#' \strong{Non parametric percentile boostrapping}
#'
#'
#' @return
#' \code{L0glm.inference} returns a list containing the following elements:
#' \item{\code{CI.lower}}{The lower bound of the confidence interval for every
#' coefficient.}
#' \item{\code{CI.upper}}{The lower upper of the confidence interval for every
#' coefficient.}
#' \item{\code{constraint}}{The constraint used to compute \code{fit}. "nonneg"
#' for nonnegativity constraint; "none" if there was no constraint.}
#' \item{\code{estimates}}{The estimated coefficients (see coefficients in \code{\link{L0glm}}).}
#' \item{\code{family}}{The error structure used to fit the data.}
#' \item{\code{method}}{The inference method used (see details).}
#' \item{\code{p.value}}{The infered type I error rate for every estimate.}
#' \item{\code{penalty}}{The type of penalty used to compute \code{fit}.}
#' \strong{If inference is performed using the Fischer information matrix}
#' \item{\code{SE}}{The standard error on the estimates}
#' \item{\code{vcov}}{The variance-covariance matrix of the estimates.}
#' \item{\code{z}}{The z-score associated to every }
#' \strong{If inference is performed using bootstrapping}
#' \item{\code{boot.result}}{The output of the bootstrapping. See \code{\link{boot}}
#' for more details.}
#'
#' @seealso
#' \code{\link{L0glm}}, \code{\link{control.l0.gen}}, \code{\link{control.iwls.gen}},
#' \code{\link{control.fit.gen}}, \code{\link{boot}}
#'
#' @example examples/L0glm_inference_examples.R
#'
#' @export
L0glm.inference <- function(fit, level = 0.95, boot.repl = 200,
                            control.l0, control.iwls, control.fit,
                            verbose = TRUE){

  # Extract data structure from the formula (code taken from stats::glm)
  formula <- eval(fit$call$formula)
  data <- eval(fit$call$data)
  mf <- fit$call
  if (is.null(data)) data <- environment(formula)
  m <- match(c("formula", "data", "weights"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  # Get the response variable
  y <- model.response(mf, "any")
  if (length(dim(y)) == 1L) {
    nm <- rownames(y)
    dim(y) <- NULL
    names(y) <- nm
  }
  # Get the design matrix
  X <- if (!is.empty.model(mt)){
    model.matrix(mt, mf, contrasts)
  } else {
    matrix(NA, length(y), 0L)
  }
  intercept <- as.logical(attr(mt, "intercept"))

  control.l0 <- do.call("control.l0.gen", control.l0)
  control.iwls <- do.call("control.iwls.gen", control.iwls)
  control.fit <- do.call("control.fit.gen", control.fit)

  if(!inherits(fit, "L0glm")) stop("Function supports only 'L0glm' class objects. Make sure 'fit' is an object returned from 'L0glm(...)'.")
  if(is.matrix(fit$coefficients)) stop("Lambda tuning should be performed first. See lambda tuning section in '?L0glm'.")

  out <- list()
  fam <- out$family <- fit$family
  out$level <- level
  alpha <- 1 -level
  out$constraint <- fit$constraint
  out$penalty <- if(fit$iter.l0 > 1){
    "L0"
  } else if (!is.null(fit$prior.coefficients)){
    "adaptive ridge"
  } else if (sum(fit$lambda != 0)){
    "ridge"
  } else {
    "none"
  }
  out$method <- ifelse(out$constraint == "none" && out$penalty %in% c("ridge", "none"), "FIM", "bootstrap")
  if(verbose) cat(paste0("==== Inference ====\n\n",
                         "Nonnegativity constraint: ", ifelse(out$constraint == "none", "no", "yes"), "\n",
                         "Penalty type: ", out$penalty, "\n",
                         "Method: ", ifelse(out$method == "FIM",
                                            "Inference based on the Fischer information matrix",
                                            "Inference based on non-parametric percentile bootstrapping"), "\n"))

  if(out$method == "FIM"){ # For ordinary regression or (non adaptive) ridge regression, compute pvalues numerically using the Fischer information matrix
    out$estimates <- fit$coefficients
    lambda <- fit$lambda * fit$lambda.w

    # Estimate dispersion (1 if noise is poisson or binomial)
    if(fam$family %in% c("poisson","binomial")){
      dispersion <- 1
    } else {
      dispersion <- sum(fit$weights * fit$residuals^2)/fit$df$rdf
    }
    Xw <- if(sum(lambda) != 0){
      sqrt(fit$weights) * rbind(X, diag(sqrt(lambda),ncol(X))) # the weighted row augmented matrix subseted for positive covariates
    } else {
      sqrt(fit$weights) * X
    }
    FIM <- (1/dispersion) * crossprod(Xw)
    # Variance-covariance matrix of coefs
    out$vcov <- tryCatch(solve(FIM),
                         error = function(e) return(e$message))
    if(is.character(out$vcov) && grepl(out$vcov, pattern = "condition number")){
      stop("Cannot compute Fischer information matrix because system is computationally singular. Maybe refit data with a small ridge penalty, or increase lambda.")
    }
    # The approximate standard error
    out$SE <- sqrt(diag(out$vcov))
    # Confidence intervals on back tranformed scale to impose non-negativity
    out$CI.lower <- out$estimates - qnorm(alpha/2, lower.tail = FALSE) * out$SE
    out$CI.upper <- out$estimates + qnorm(alpha/2, lower.tail = FALSE) * out$SE
    # Compute the p-values
    out$z <- abs(out$estimates/out$SE)
    out$p.value <- pnorm(out$z, lower.tail = FALSE)*2
  } else { # Otherwise Coefficient inference using bootstraping
    # Normalize data
    normalize <- fit$call$normalize
    if(is.null(normalize)) normalize <- TRUE
    if(normalize){
      X.n <- apply(X, 2, norm, type = "2")
      X <- sweep(X, 2, X.n, "/")
    }

    lambda <- max(unique(fit$lambda))
    no.pen <- which(fit$lambda == 0)
    if(verbose){
      cat(paste0("\nFitting ", boot.repl, " bootstrap replicates\n",
                 "Estimated time: ", round(fit$timing * boot.repl), " secs...\n"))
    }
    out$boot.result <- boot(data = y, X = X, wts = fit$prior.weights, family = fam,
                            lambda = fit$lambda, start = fit$prior.coefficients,
                            nonnegative = ifelse(out$constraint == "none", FALSE, TRUE),
                            control.l0 = control.l0, control.iwls = control.iwls,
                            control.fit = control.fit, post.filter.fn = fit$post.filter.fn,
                            # boot arguments
                            statistic = L0glm.bfun, R = boot.repl, stype = "i")
    out$estimates <- out$boot.result$t0
    # Reassign norms
    if(normalize){
      out$estimates <- out$estimates * X.n
      out$boot.result$t <- sweep(out$boot.result$t, 2, X.n, "*")
    }
    # Confidence intervals
    # Check https://stats.stackexchange.com/questions/20701/computing-p-value-using-bootstrap-with-r
    # or https://tolstoy.newcastle.edu.au/R/e6/help/09/04/11096.html
    # https://stats.stackexchange.com/questions/231074/confidence-intervals-on-predictions-for-a-non-linear-mixed-model-nlme
    # Perform inference using percentile bootstrap
    out$CI.lower <- apply(out$boot.result$t, 2, quantile, probs = alpha/2)
    out$CI.upper <- apply(out$boot.result$t, 2, quantile, probs = 1-alpha/2)
    if(out$constraint == "none"){
      out$p.value <- sapply(1:length(out$estimates), function(i){
        # First exclude the 0
        is.zero <- out$boot.result$t[,i] != 0
        theta <- out$boot.result$t[!is.zero,i]
        # Center the coefficient distribution. We assume that the coefficient
        # distribution does not depend on the true coefficient value !
        theta.c <- theta - mean(theta)
        # Reinsert the zero
        p.val <- (sum(theta.c >= abs(out$estimates[i]))+sum(is.zero)+1)/(boot.repl+2)  # abs for 2-sided; constants are added as in Laplace's rule of succession
        return(p.val)
      })
    } else { # if nonnegative
      out$p.value <- (colSums(out$boot.result$t == 0)+1)/(boot.repl+2)  # constants are added as in Laplace's rule of succession
    }
  }
  out <- out[sort(names(out))]
  return(out)
}

# NOTES:
# SEs & post-selection inference on nonzero coefficients:
# see https://stats.stackexchange.com/questions/373253/calculating-the-p-values-in-a-constrained-non-negative-least-squares/405606#405606
# and https://stats.stackexchange.com/questions/14471/how-do-you-calculate-standard-errors-for-a-transformation-of-the-mle/14472
# and  stat.psu.edu/~sesa/stat504/Lecture/lec2part2.pdf
# Method is based on calculating Fisher information matrix (negative of 2nd
# derivate of log-likelihood) for covariates with nonzero coefs on log scale
# to make log-likelihood surface more symmetric & enforce nonneg constraint
#
# Observed Fisher information matrix for nonzero coefs is the negative of the
# 2nd derivative (Hessian) of the log likelihood at parameter estimates. We
# assume the non zero coefficients follow a normal distribution for which
# the associated log likelihood is:
# l(mu, sig^2; y) = -n/2*log(2*pi) - n/2*log(sig^2) - 1/(2*sig^2)*sum((y - mu)^2)
# where sig^2 is the dispersion factor, mu = X_pos %*% beta_pos, _pos
# indicates the subset for which beta is positive.
#
# First derivative wrt to beta_pos = delta l/delta beta_pos = 1/sig^2 * (t(y)%*%X_pos - t(X_pos)%*%X_pos%*%beta_pos)
# Second derivative wrt to beta_pos = delta^2 l/delta beta_pos^2 = -1/(2*sig^2)* t(X_pos) %*% X_pos
# Fischer information matrix FIM = - delta^2 l/delta beta_pos^2 = 1/(2*sig^2)* t(X_pos) %*% X_pos
# So,
# FIM <- (1/dispersion) * crossprod(Xw)
#
# Let's now calculate the information matrix on a log transformed Y scale to
# take into account the nonnegativity constraints on the parameters
# see stat.psu.edu/~sesa/stat504/Lecture/lec2part2.pdf, slide 20 eqn 8 & Table 1
# and https://stats.stackexchange.com/questions/14471/how-do-you-calculate-standard-errors-for-a-transformation-of-the-mle/14472
# and https://stats.stackexchange.com/questions/373253/calculating-the-p-values-in-a-constrained-non-negative-least-squares/405606#405606
# FIM(Phi) = FIM/(Phi'(beta_pos)^2) where in our case Phi(x)=log(x) => Phi'(x) = 1/x
# So,
# FIM_phi <- FIM/(1/beta_pos)^2
#
# Variance-covariance matrix of coefs on log scale, i.e. of log(beta) is 1/FIM_phi
# vcov_log <- solve(FIM_phi)
#
# The approximate standard error of the positive coefficients (on log scale) is
# SE_logbeta_X = sqrt(diag(VCOV))[X]
#
# Confidence intervals on back tranformed scale
# CI = exp( log(beta_pos) +- q_level * SE_logbeta
#
# Compute the p-values
# One-sided p-values for log(coefs) being greater than 0, ie coefs being > 1 (since log(1) = 0)
# z = log(beta_pos)/SE_logbetapos

