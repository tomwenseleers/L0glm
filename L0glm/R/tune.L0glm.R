#' @importFrom qlcMatrix colMax
#' 
#' @describeIn L0glm
#' Selects optimal value of regularization parameter lambda
#' @export
tune.L0glm <- function(X, y, family, 
                       weights = matrix(1, nrow = nrow(y), ncol = ncol(y)), 
                       offset = Matrix(0, nrow = nrow(X), ncol=ncol(y), sparse=T), 
                       dispersion = NULL, 
                       IC = "aic", 
                       warm.start.ratio = 0.9, 
                       coefs.true = NULL,
                       ...){

  # Checking if the IC is one of the valid criteria
  IC = tolower(IC)
  if(!(IC %in% c("aic", "bic", "ebic", "gic", "mbic"))){
    stop("Invalid information criteria.")
  }

  if(!is.null(dispersion)){
    lambda.preset <- preset.lambda(IC = IC, dispersion = dispersion, family = family)
  } else {
    lambda.preset = preset.lambda(IC = IC, y = y[,which.max(colMeans(abs(y)))],
                                      X = X, weights=rep(1, nrow(y)),
                                  family = family)
  }
  
  
  # TODO: Calculate Xw and yw propperly
  Xw <- Xweights(X, weights[,1]) # TODO: talk with Tom about this
  Yw <- sapply(seq_len(ncol(Y)), function(col) {
    yw <- Yweights(Y[, col, drop = FALSE], weights[, col, drop = FALSE])
    return(yw)
  })
  
  # for the lambda sequence below we use lambda_preset_aic 
  # (round(log10(lambda_preset_aic))-3) as the lower bound and the theoretical
  # lambda that would penalize all variables to zero with an L0 penalty as maximum
  # see https://stats.stackexchange.com/questions/416144/minimum-and-maximum-regularization-in-l0-pseudonorm-penalized-regression
  maxlambda = max((1/colSums(X^2)) * (crossprod(X, Y)^2))
  lambdas = 10^seq((round(log10(lambda.preset))-3), round(log10(maxlambda))+1, by=1)

  # Get warm start
  warm_start <- L0glm.fit_new(X = X, y = y,
                              lambda = min(lambdas), 
                              family = family,
                              offset = offset, 
                              weights = weights,
                              ...) # using provided / default params
  
  
  no.warm.start = matrix(1, nrow=nrow(warm_start), ncol=ncol(warm_start))
  coefficients <- lapply(lambdas, function (lambda) {
    new_coefs <- L0glm.fit_new(X = X, y = y, 
                               family = family,
                               start = warm.start.ratio*warm_start+(1-warm.start.ratio)*no.warm.start, # warm.start.ratio= 0 => no warm starts,warm.start.ratio= 1 => full warm starts
                               lambda = lambda, # needs to be larger with larger sd_noise for family=gaussian, e.g. 0.5 with sd_noise=1 & 5 with sd_noise=5 
                               weights = weights,
                               offset = offset,
                               ...)
    warm_start <<- new_coefs
  }) # 82s, TODO check how warm starts would affect solution quality & stability
  
   metrics <- data.frame(do.call(rbind, lapply(coefficients, function (coefs) get.metrics(beta = coefs, X = X, y = y, 
                                                                                         weights = weights, 
                                                                                         offset = offset, 
                                                                                         family = family, 
                                                                                         beta_true = coefs.true))))

  # metrics <- structure(metrics, class = "tune.L0glm") # Deformats the data frame. TODO: should I make a constructor for the tune.L0glm class?
  return(metrics)
}

abessfit <- function(X, y, ...){
  fit <- abess(x = X, 
        y = y, 
        ...)
  
  
  plot(fit)
  return(fit)
}

glmnetfit <- function(X, y, ...){
  
}
