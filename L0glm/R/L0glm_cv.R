# L0glm.trainval ####
# Function that tunes a sequence of
L0glm.trainval <- function(X, y,
                           weights=rep(1,length(y)),
                           family = poisson(identity),
                           start=NULL,
                           lambdas = 0,
                           no.pen = 0,
                           nonnegative = FALSE,
                           control.l0 = list(maxit = 100, rel.tol = 1E-4, delta = 1E-5, gamma = 2, warn = FALSE),
                           control.iwls = list(maxit = 100, rel.tol = 1E-4, thresh = 1E-5, warn = FALSE),
                           control.fit = list(maxit = 10, block.size = NULL, tol = 1E-7),
                           post.filter.fn = function(u) return(u),
                           tune.crit = "bic",
                           seed = NULL,
                           verbose = FALSE){
  n <- nrow(X)
  p <- ncol(X)
  # Subset the train and validation set
  if(!is.null(seed)) set.seed(seed)
  if(n < 10) stop("At least 10 observations are needed to split data in training and
                  validation set.")
  train <- sample(1:n, size = round(n*0.9))
  fits <- lapply(lambdas, function(lambda){
    fit <- L0glm.qual(lambda = lambda, crit = "all", no.pen = no.pen,
                      X = X[train,], y = y[train], weights = weights[train],
                      x.val = X[-train,,drop=F], y.val = y[-train],
                      w.val = weights[-train], family = family,start = start,
                      nonnegative = nonnegative, control.l0 = control.l0,
                      control.iwls = control.iwls, control.fit = control.fit,
                      post.filter.fn = post.filter.fn)
    if(verbose) print.progress(which(lambda == lambdas), length(lambdas))
    return(fit)
  })
  res <- list()
  res$lambdas <- lambdas
  res$IC <- do.call(rbind, lapply(fits, "[[", "IC"))
  res$coefficients.lam <- t(sapply(fits, "[[", "coefficients.lam"))
  rownames(res$IC) <- rownames(res$coefficients.lam) <- res$lambdas
  if(tune.crit %in% c("loglik", "R2", "adjR2")){
    best <- which.max(res$IC[,tune.crit])
  } else {
    best <- which.min(res$IC[,tune.crit])
  }
  res$best.lam <- lambdas[best]

  # TODO delete
  # plot(x = res$lambdas, y = res$IC[,tune.crit], log = "xy", type = "l")
  # abline(v = res$best.lam)
  # stop("ierjgruioeghe")

  return(res)
}


L0glm.cv <- function(X, y,
                     weights=rep(1,length(y)),
                     family = poisson(identity),
                     start=NULL,
                     lambdas = 0,
                     no.pen = 0,
                     nonnegative = FALSE,
                     control.l0 = list(maxit = 100, rel.tol = 1E-4, delta = 1E-5, gamma = 2, warn = FALSE),
                     control.iwls = list(maxit = 100, rel.tol = 1E-4, thresh = 1E-5, warn = FALSE),
                     control.fit = list(maxit = 10, block.size = NULL, tol = 1E-7),
                     post.filter.fn = function(u) return(u),
                     tune.crit = "bic",
                     k = 5,
                     seed = NULL,
                     verbose = FALSE){
  n <- nrow(X)
  p <- ncol(X)
  # Subset the train and validation set
  if(!is.null(seed)) set.seed(seed)
  foldid <- sample(rep(1:k, length = n))
  fits <- t(sapply(lambdas, function(lambda){
    fits.fold <- lapply(1:k, function(fold){
      val.ind <- foldid == fold
      fit <- L0glm.qual(lambda = lambda, crit = "all", no.pen = no.pen,
                        X = X[!val.ind,,drop=F], y = y[!val.ind], weights = weights[!val.ind],
                        x.val = X[val.ind,,drop=F], y.val = y[val.ind], w.val = weights[val.ind],
                        family = family, start = start, nonnegative = nonnegative,
                        control.l0 = control.l0, control.iwls = control.iwls,
                        control.fit = control.fit, post.filter.fn = post.filter.fn)
      return(fit$IC)
    })
    out <- rowMeans(do.call(cbind, fits.fold)) # compute the average criterion over the folds
    # TODO delete
    # cat(paste0("Lambda = ", lambda, " | ", tune.crit, " = ", round(out[tune.crit], 2), "\n"))
    if(verbose) print.progress(which(lambdas == lambda), length(lambdas))
    return(out)
  }))
  res <- list()
  res$lambdas <- lambdas
  res$IC <- fits
  res$coefficients.lam <- NA
  rownames(res$IC) <- res$lambdas
  if(tune.crit %in% c("loglik", "R2", "adjR2")){
    best <- which.max(res$IC[,tune.crit])
  } else {
    best <- which.min(res$IC[,tune.crit])
  }
  res$best.lam <- lambdas[best]

  # plot(x = res$lambdas, y = res$IC[,tune.crit], log = "xy", type = "l")
  # abline(v = res$best.lam)

  return(res)
}

L0glm.IC <- function(X, y,
                     weights=rep(1,length(y)),
                     family = poisson(identity),
                     start=NULL,
                     lambdas = 0,
                     no.pen = 0,
                     nonnegative = FALSE,
                     control.l0 = list(maxit = 100, rel.tol = 1E-4, delta = 1E-5, gamma = 2, warn = FALSE),
                     control.fit = list(maxit = 10, block.size = NULL, tol = 1E-7),
                     control.iwls = list(maxit = 100, rel.tol = 1E-4, thresh = 1E-5, warn = FALSE),
                     post.filter.fn = function(u) return(u),
                     tune.crit = "bic",
                     verbose = FALSE){

  fits <- lapply(lambdas, function(lambda){
    fit <- L0glm.qual(lambda = lambda, crit = "all", no.pen = no.pen,
                      X = X, y = y, weights = weights, x.val = X, y.val = y,
                      w.val = weights, family = family, start = start,
                      nonnegative = nonnegative, control.l0 = control.l0,
                      control.iwls = control.iwls, control.fit = control.fit,
                      post.filter.fn = post.filter.fn)
    # TODO delete
    # cat(paste0("Lambda = ", lambda, " | ", tune.crit, " = ", round(fit$IC[tune.crit], 2), "\n"))
    if(verbose) print.progress(which(lambda == lambdas), length(lambdas))
    return(fit)
  })
  res <- list()
  res$lambdas <- lambdas
  res$IC <- do.call(rbind, lapply(fits, "[[", "IC"))
  res$coefficients.lam <- t(sapply(fits, "[[", "coefficients.lam"))
  rownames(res$IC) <- rownames(res$coefficients.lam) <- res$lambdas
  if(tune.crit %in% c("loglik", "R2", "adjR2")){
    best <- which.max(res$IC[,tune.crit])
  } else {
    best <- which.min(res$IC[,tune.crit])
  }
  res$best.lam <- lambdas[best]

  # TODO delete
  # plot(x = res$lambdas, y = res$IC[,tune.crit], log = "xy", type = "l")
  # abline(v = res$best.lam)
  # stop("ierjgruioeghe")

  return(res)
}


# L0glm.qual ####
# Functions to optimize optimal lambda
L0glm.qual <- function(lambda, crit, no.pen,  X, y, weights, x.val, y.val,
                       w.val, nonnegative, family, start, control.l0,
                       control.iwls, control.fit, post.filter.fn){
  lambda <- rep(lambda, ncol(X))
  lambda[no.pen] <- 0
  fit <- L0glm.fit(X = X, y = y, weights = weights, family = family,
                   lambda = lambda, start = start, nonnegative = nonnegative,
                   control.l0 = control.l0, control.iwls = control.iwls,
                   control.fit = control.fit, post.filter.fn = post.filter.fn)
  val <- vector()
  # Compute effective degrees of freedom
  fit$df <- compute.df(fit = fit, X = X)
  nonzero <- fit$coefficients != 0
  eta.hat <- as.vector(x.val[,nonzero,drop=F] %*% fit$coefficients[nonzero])
  y.hat <- family$linkinv(eta.hat)

  # Compute the selection criteria
  if(crit %in% c("all", "loocv")){
    z.res <- (y.val - y.hat)/family$mu.eta(eta.hat) # the residuals on the adjusted scale = z.res = z - z.fit = eta + (y - mu)/g'(eta) - eta = (y - mu)/g'(eta)
    w <- w.val * as.vector(family$mu.eta(eta.hat)^2 / family$variance(y.hat)) # IWLS weights from last iteration
    H.diag <- compute.diagH(fit = fit, X = x.val, w = w)
    val["loocv"] <- 1/length(y.val) * sum(w * (z.res/(1 - H.diag))^2) # see p.7 in https://scholarworks.gsu.edu/cgi/viewcontent.cgi?referer=https://scholar.google.com/&httpsredir=1&article=1100&context=math_theses
  }
  if(crit %in% c("all", "rss", "bic", "aic", "loglik", "aicc", "ebic", "hq",
                 "ric", "mric", "cic",  "bicg", "bicq")){ # TODO check all IC are present
    val <- c(val, compute.ic(y = y.val, y.hat = y.hat, w = w.val, fit = fit))
    if(crit != "all") val <- val[crit]
  }
  if(length(val) == 0) stop("Invalid IC supplied.")

  # TODO delete
  # plot(x = 1:length(y), y = y, main = paste0("Lambda = ", lambda))
  # lines(x = seq(1, length(y), length.out = length(fit$coefficients)), y = fit$coefficients, col = 2, type = "h")
  # lines(x = 1:length(y), y = fit$fitted.values)

  # Return the value to minimize
  return(list(IC = val, coefficients.lam = fit$coefficients, lambda = lambda))
}


