#' @name L0glm.tune
#' @aliases L0glm.tune
#' @title Tunes the regularisation parameter lambda of an L0glm fit.
#'
#' @importFrom qlcMatrix colMax colMin
#' 
#' @description
#' Tunes the lambda parameter for the \code{\link{L0glm.fit}} 
#' function. It uses a sequence of penalty parameter values (\code{lambda}) to fit 
#' L0-penalized GLMs, and selects the best  lambda value according to a specified 
#' information criterion. It can handle multiple output variables (channels), 
#' and can produce a plots of the regularization path, true vs. estimated coefficients (if true coefficients are passed, e.g. from a simulation),
#' estimated support size, deviance, log likelihood, mean absolute deviation, 
#' accuracy, precision, recall and the IC of interest, all in function of the level of
#' regularization (L0 penalty lambda). Some plots are only produced if the user 
#' passes the true coefficients \code{coefs.true} (e.g. from a simulation).
#' 
#' @param X 
#' A matrix of predictor variables.
#' @param y 
#' A vector or matrix of response variables.
#' @param family 
#' A description of the error distribution and link function to be
#' used in the model. This should be a call to a family function. (See
#' \code{stats::family} for details of family functions.)
#' @param weights 
#' A matrix of observation weights.
#' @param offset 
#' A matrix to be added to the linear predictor.
#' @param dispersion 
#' A dispersion parameter for the family function. If not NULL, it is used to 
#' set the lambda sequence.
#' @param IC A character string specifying the information criterion to use for
#' model selection. It should be one of \code{"aic"}, 
#' \code{"bic"}, \code{"gic"}, \code{"ebic"}, \code{"mbic"}.
#' @param nlambda Number of lambda values to use for tuning (default 10). 
#' If \code{nlambda} is even, half of the values lie above 
#' the lambda value that would be expected to optimise
#' the chosen IC (\code{lambda_IC}) and half below;  
#' if odd, the upper branch receives the extra point so that
#' \code{lambda_IC} is included exactly once.
#' @param prop_range Numeric in \((0,\,1]\).  
#' Proportion of the log10 distance between \code{lambda_IC} and
#' \code{lambda_max} (that would penalise all coefficients to zero
#' or the chosen prior.mean) that is *retained on each side* of 
#' \code{lambda_IC}. \code{prop_range = 1} (default) keeps the 
#' entire range, whilst smaller values "zoom‑in", yielding grids 
#' that are denser around \code{lambda_IC}. For example, 
#' \code{prop_range = 0.1} preserves only 10% of the original 
#' log10 span above and below the preset lambda value that 
#' a priori would be expected to optimise the chosen IC.
#' @param warm.start A numeric value indicating the ratio of warm start 
#' coefficients to use. Default 0.01.
#' @param coefs.true A vector or matrix of true coefficient values, for comparison 
#' with the estimated coefficients.
#' @param printbest A logical value indicating whether to print the best lambda 
#' value according to the selected IC (default is TRUE).
#' @param plot A logical value indicating whether to produce the plots mentioned
#' in the description, all in function of regularization (lambda).
#' @param ... Other arguments to pass to the `L0glm.fit` function except for \code{start}.
#' 
#' @return
#' \code{L0glm.tune} returns a list containing the following elements:
#' \item{\code{IC}}{chosen \code{IC}.}
#' \item{\code{metrics}}{a data frame containing metrics evaluating the models fit
#' with different regularization parameters.}
#' \item{\code{lambda_best}}{best lambda according to the \code{IC} 
#' criterion specified.}
#' \item{\code{lambda_IC}}{lambda that was a priori expected to approx. optimise
#'the chosen \code{IC}.}
#' \item{\code{lambdas}}{a vector containing all the lambdas that were used.}
#' 
#' @export
L0glm.tune <- function(X, y,
                       family,
                       weights      = matrix(1, nrow = nrow(as.matrix(y)),
                                             ncol = ncol(as.matrix(y))),
                       offset       = Matrix::Matrix(0, nrow = nrow(X),
                                                     ncol = ncol(as.matrix(y)),
                                                     sparse = TRUE),
                       dispersion   = NULL,
                       IC           = c("gic", "aic", "bic", "ebic", "mbic"),
                       nlambda      = 10,
                       prop_range   = 1,
                       warm.start   = 0.01,
                       coefs.true   = NULL,
                       printbest    = TRUE,
                       plot         = TRUE,
                       ...)
{
  ## ------------------------------------------------------------------ set-up
  IC <- match.arg(tolower(IC),
                  choices = c("gic", "aic", "bic", "ebic", "mbic"))
  
  y <- as.matrix(y)
  n <- nrow(X); p <- ncol(X); m <- ncol(y)
  
  ## ---- warm-up fit (used to build lambda-grid and optional warm start) ----------
  lambda_preset  <- preset.lambda(IC = IC, p = p, n = n)
  
  warmup <- L0glm.fit(X = X, y = y,
                      lambda     = lambda_preset,
                      family     = family,
                      weights    = weights,
                      offset     = offset,
                      dispersion = dispersion,
                      ...)
  
  beta0   <- warmup$coefficients          # p × m (sparse)
  Z    <- warmup$z
  W    <- warmup$weights
  
  if (all(beta0 == 0))
    stop("warm-up produced only zero coefficients — cannot build lambda sequence.")
  
  ## ---- lambda grid (ceiling(nlambda/2) points between lambda_IC and lambda_max on log10 scale, remaining below lambda_IC with same log10 step size)
  lambda_max <- maxlambda(X, weights = W, z = Z,
                     coefs = beta0, family = family) # max lambda where all coefs would be penalised to zero
  lambda_IC <- preset.lambda(IC=IC, p, n) # lambda that would be expected to approx optimise the chosen IC
  n_up   <- ceiling(nlambda / 2)        # includes lambda_IC
  n_down <- nlambda - n_up              # below lambda_IC
  # ---- log10 limits, shrunk by prop_range 
  g_IC  <- log10(lambda_IC)
  g_max <- log10(lambda_max)
  g_up_max <- g_IC + prop_range * (g_max - g_IC)   # new “effective” lambda_max
  log_step <- (g_up_max - g_IC) / (n_up - 1)
  log_up   <- g_IC + log_step * seq(0, n_up - 1)
  log_down <- g_IC - log_step * seq_len(n_down)
  lambda_grid <- sort(10^c(log_up, log_down))
  
  ## ---- optional warm-start coefficients for every lambda ------------------------
  if (warm.start > 0) {
    beta_const <- matrix(start_coefs(X, family = family), nrow = p, ncol = m)
    beta_start <- beta0 * warm.start + beta_const * (1 - warm.start)
    beta_start[beta0 == 0] <- 0                      # keep strict zeros
  } else {
    beta_start <- NULL
  }
  
  ## containers ---------------------------------------------------------------
  coefmat  <- Matrix::Matrix(0, nrow = length(lambda_grid), ncol = p,
                             sparse = TRUE)         # row-means of bet̂a(lambda)
  metric_rows <- vector("list", length(lambda_grid))     # each element numeric
  
  ## ------------------------------------------------------------------ loop --
  for (i in seq_along(lambda_grid)) {
    lambda <- lambda_grid[i]
    
    fit <- L0glm.fit(X, y,
                     family     = family,
                     weights    = weights,
                     lambda     = lambda,
                     dispersion = dispersion,
                     start      = beta_start,
                     offset     = offset,
                     ...)
    
    beta_hat <- fit$coefficients                           # p × m (sparse)
    
    ## store coefficient path (row mean for multichannel)
    coefmat[i, ] <- if (m > 1) Matrix(Matrix::rowMeans(beta_hat), sparse=TRUE) else beta_hat
    
    ## quality metrics --------------------------------------------------------
    met <- get.metrics(beta = beta_hat,
                       X = X, y = y, family = family,
                       beta_true = coefs.true)

    ## numeric row *not* data.frame
    metric_rows[[i]] <- c(lambda = unname(lambda), unlist(met))
  }
  
  ## ----------------------------- build tidy data-frame ----------------------
  metrics <- as.data.frame(do.call(rbind, metric_rows),
                           check.names = FALSE,
                           stringsAsFactors = FALSE)

  ## ------------------------------------------------------------- best lambda -----
  best.index  <- which.min(metrics[[IC]])
  lambda_best <- metrics$lambda[best.index]
  
  if (printbest)
    message(sprintf("Best lambda by %s = %.3g", toupper(IC), lambda_best))
  
  ## ---------------------------------------------------------------- plotting
  if (plot) {
    old_par <- par(no.readonly = TRUE); on.exit(par(old_par))
    layout(matrix(1:2, nrow = 1))          # side-by-side
    ## 1. coefficient paths
    matplot(x  = log10(metrics$lambda),
            y  = as.matrix(coefmat),
            type = "l", lty = 1, col = "grey40",
            xlab = expression(log[10](lambda)),
            ylab = "Coefficient estimate",
            main = "Regularisation path")
    abline(v = log10(lambda_max), col = "black", lty = 3)
    abline(v = log10(lambda_best), col = "red", lty = 2)
    abline(v = log10(lambda_IC), col = "red", lty = 3)
    if (!is.null(coefs.true)) {
      if (m > 1) coefs.true <- Matrix::rowMeans(coefs.true)
      abline(h = as.vector(coefs.true[coefs.true!=0]),
             col = scales::alpha("blue", 0.1)) 
    }
    
    ## 2. chosen information criterion
    plot(log10(metrics$lambda), metrics[[IC]],
         type = "b", pch = 19,
         xlab = expression(log[10](lambda)), ylab = toupper(IC),
         main = paste("Tuning –", toupper(IC)))
    abline(v = log10(lambda_best), col = "red", lty = 2)
    abline(v = log10(lambda_IC), col = "red", lty = 3)
  }
  
  ## ----------------------------------------------------------------- return
  list(IC           = IC, # chosen IC
       metrics      = metrics,
       lambda_best  = lambda_best, # 
       lambda_IC    = lambda_IC, # lambda preset expected to optimise chosen IC
       lambdas      = lambda_grid)
}

  # TODO; maybe take the lambda sequence from 100x lower to 100x higher than the preset
  # instead of up to maxlambda which is unrealistically high probably
  # TODO: add more returns
  # - the optimal values of lambda & the index in lambdas of that optimal value based on given metrics - so this could also 
  # be a dataframe with columns lambda, lambda_index, metric (e.g. different information criteria plus accuracy & F1 score for 
  # support recovery & metrics for bias in coefficients  if beta_true was given
  # - a sparse matrix with the coefficient estimates for those lambda values (rowMeans of the coefficients for multichannel) (with the right optimal threshold applied for each)
  # - a list with the originally estimated sparse coefficient matrices for the different lambda values (but with the right optimal threshold applied for each)




tune.plots <- function(metrics, coefficients, lambdas, IC = "aic", coefs.true = NULL){
  IC = tolower(IC)
  if(!(IC %in% c("aic", "bic", "ebic", "gic", "mbic"))){
    stop("Invalid information criteria.")
  }
  bestlambda = lambdas[which.min(metrics[,IC])]
  bestlambda = lambdas[which.min(abs(lambdas-bestlambda))] 
  beta_relaxedgroupL0 = as.matrix(coefficients[which(metrics$lambda==bestlambda),]) 
  IC = toupper(IC) # better for plot titles
  # PLOT TRUE VS ESTIMATED COEFFICIENTS 
  if(!is.null(coefs.true)){
    R2 = round(cor(as.vector(coefs.true), as.vector(beta_relaxedgroupL0)), digits = 4)^2 # Calculate R2
    title = paste("True vs estimated coefficients\nnonnegative L0 group GLM\nR2=", R2)
    plot(as.vector(beta_relaxedgroupL0), as.vector(coefs.true), pch=16, 
         xlab="estimated coefficients nonnegative L0 group GLM",
         ylab="true coefficients", col="steelblue",
         main=title)
  }
  
  # PLOT ESTIMATED SUPPORT SIZE IN FUNCTION OF REGULARISATION
  if(!is.null(coefs.true)){
    title = paste("Estimated support size in function of regularisation\n(red = estimate at best", IC, ", cyan = true support size)")
  } else {
    title = paste("Estimated support size in function of regularisation\n(red = estimate at best", IC, ")")
  }
  
  
  plot(x = log10(lambdas), y = metrics$support_size_estimated,
       ylab="estimated support size",
       main = title,
       pch=16, type="o")
  abline(v = log10(bestlambda), col="red")
  if(!is.null(coefs.true)) abline(h = metrics$support_size_true, col="cyan")
  
  # PLOT DEVIANCE
  title = paste("Theoretical deviance relative to true underlying signal without noise\n(red=estimate at best", IC,", cyan=best deviance)")
  plot(x=log10(lambdas), y=metrics$deviance, pch=16, type="o", 
       ylim=c(0,max(unlist(metrics$deviance))), 
       ylab="deviance",
       main = title) # theoretical deviance relative to true underlying signal without noise
  lambda_bestdeviance = lambdas[which.min(metrics$deviance)]
  abline(v=log10(bestlambda), col="red")
  abline(v=log10(lambda_bestdeviance), col="cyan")
  par(oma=c(0,2,0,0))
  
  # PLOT LOG LIKELIHOOD IN FUNCTION OF REGULARISATION
  title = paste("Log likelihood in function of regularisation\n(red = estimate at best", IC,", cyan = estimate at best log likelihood)")
  plot(x=log10(lambdas), y=metrics$loglik,
       main=title,
       ylab="log likelihood",
       pch=16, type="o") 
  lambda_bestloglik = lambdas[which.max(metrics$loglik)]
  abline(v=log10(bestlambda), col="red")
  abline(v=log10(lambda_bestloglik), col="cyan")
  
  # PLOT MAD
  # TODO plot MD too?
  title = paste("Mean absolute deviation of coefficients vs true coefs\n(red = estimate at best", IC,")")
  plot(x=log10(lambdas), y=metrics$MAD_coefs, 
       main=title,
       ylab="Mean absolute deviation coefficients",
       pch=16, type="o")
  abline(v=log10(bestlambda), col="red")
  
  # PLOT ACCURACY
  title = paste("Accuracy to recover true support\n(red = estimate at best", IC,", cyan = estimate at best accuracy)")
  plot(x=log10(lambdas), y=metrics$accuracy,
       main=title,
       ylab="Accuracy",
       pch=16, type="o") 
  lambda_best_accuracy = lambdas[which.max(metrics$accuracy)] 
  abline(v=log10(bestlambda), col="red")
  abline(v=log10(lambda_best_accuracy), col="cyan")
  
  # PLOT PRECISION
  title = paste("Precision in function of regularisation\n(red = estimate at best", IC, ", blue = estimate at best precision)")
  plot(x=log10(lambdas), y=metrics$precision,
       main=title,
       ylab="Precision",
       pch=16, type="o") 
  lambda_best_precision = lambdas[which.max(metrics$precision)] 
  abline(v=log10(bestlambda), col="red")
  abline(v=log10(lambda_best_precision), col="cyan")
  
  # PLOT RECALL
  # TODO should it be precision and recall separately or F1 score?
  title = paste("Recall in function of regularisation\n(red = estimate at best", IC, ", blue = estimate at best recall)")
  plot(x=log10(lambdas), y=metrics$recall,
       main=title,
       ylab="Recall",
       pch=16, type="o") 
  lambda_best_recall = lambdas[which.max(metrics$recall)] 
  abline(v=log10(bestlambda), col="red")
  abline(v=log10(lambda_best_recall), col="cyan")
  
  # PLOT SPECIFICITY
  title = paste("Specificity in function of regularisation\n(red = estimate at best", IC, ", blue = estimate at best specificity)")
  plot(x=log10(lambdas), y=metrics$specificity,
       main=title,
       ylab="Specificity",
       pch=16, type="o") 
  lambda_best_specficity = lambdas[which.max(metrics$specificity)] 
  abline(v=log10(bestlambda), col="red")
  abline(v=log10(lambda_best_specficity), col="cyan")
  
  
  # R2
  title = paste("R2 value in function of regularisation\n(red = estimate at best", IC,")")
  plot(x=log10(lambdas), y = metrics$R2_coefs,
       main = title, 
       ylab = "R2", 
       pch = 16, type="o")
  abline(v=log10(bestlambda), col="red")
  
  # TP, FP, TN, FN
  color = c("steelblue", "purple")
  title = paste("Effect of regularisation on True and False Positives\n(red = estimate at best", IC, ")")
  plot(x=log10(lambdas), y = metrics$TP,
       main = title,
       ylab = "Value",
       col = color[1],
       pch = 16,
       type = "o")
  lines(log10(unlist(metrics$lambda)), metrics$FP, type="o", pch=16, col=color[2])
  # lines(log10(metrics$lambda), metrics$TN, type="o", pch=16, col=color[1])
  # lines(log10(metrics$lambda), metrics$FN, type="o", pch=16, col=color[4])
  
  legend("topright", legend=c("TP", "FP"), col=color, pch=16, lty=1)
  abline(v=log10(bestlambda), col="red")
  
  # PLOT IC OF INTEREST
  # TODO should it also plot other ICs?
  title = paste("Effect of regularisation on", IC, "\n(red = estimate at best", IC, ")")
  plot(x = log10(lambdas), y = metrics[[paste0(tolower(IC))]],
       main = title,
       ylab = IC, 
       pch = 16,
       type = "o")
  abline(v=log10(bestlambda), col="red")
}

optimize.threshold <- function(thresh, IC = "bic", coefficients, Xmat, y, family){
  coefficients[abs(coefficients) < thresh[1]] = 0
  IC.val <- get.metrics(beta = coefficients, 
                        X = Xmat, y = y, 
                        family = family)[[paste0(IC)]]
  return(IC.val)
}


