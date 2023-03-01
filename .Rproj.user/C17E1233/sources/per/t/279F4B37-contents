###########################
#### GENERIC FUNCTIONS ####
###########################

# TODO include generic functions

# print.L0glm <- function(x){
#  # TODO
# }
#
# summary.L0glm <- function(x){
#   # TODO
# }
#
# predict.L0glm <- function(x, newdata){
#   if (missing(newdata)) return(x$fitted.values)
#   # TODO if (x$intercept) newdata <- cbind(1, newdata) # TODO make sure this still works
#   return(newdata %*% x$coefficients)
# }
#
# residuals.L0glm <- function(x){
#   # TODO
#   print(stats:::residuals.lm)
#   return(x$residuals)
# }
#
# coefficients.L0glm <- function(x, lambda){
#   if(missing(lambda)) return(x$coefficients)
#   else stop("to implement") # TODO
#   return(x$coefficients)
# }


###############################
#### HIDDEN CORE FUNCTIONS ####
###############################


# compute.df ####
# Calculate appropriate degrees of freedom with ridge penalty.
#   - Effective degrees of freedom (edf) taken up by the fit = trace(H). This
#     can be seen as the Cp correction for the average squared residual
#   - Degrees of freedom for error, aka residual degrees of freedom
#     (rdf) = n - 2 * trace(H) + trace(H %*% t(H))
# The hat matrix is not computed explicitely but we take a shortcut using QR
# decomposition.
#
# REF
#   Cule et al., Significance testing in ridge regression for genetic data, 2011
#
# INPUT:
#   X:    the model matrix
#   fit:  the fitted model
#   return.diagH:   if TRUE, the diagonal of the hat matrix is returned
# OUTPUT:
#   edf: effective degrees of freedom
#   rdf: residual degrees of freedom
compute.df <- function(X, fit){
  n <- sum(fit$weights != 0)
  lambda <- fit$lambda.w * fit$lambda
  if(sum(lambda) != 0){
    X <- rbind(X, diag(sqrt(lambda),ncol(X)))
    rw <- c(sqrt(fit$weights), rep(1, nrow(X) - length(fit$weights)))
  }

  # Calculate trace of the hat matrix using QR decomposition:
  # see https://stackoverflow.com/questions/20562177/get-hat-matrix-from-qr-decomposition-for-weighted-least-square-regression
  QR <- qr.default(X, LAPACK = TRUE)
  if(sum(lambda) == 0){
    edf <- QR$rank # p
    rdf <- n - edf # n - p
    return(list(edf = edf, rdf = rdf))
  }
  if (QR$rank != 0) {
    Q <- qr.qy(QR, diag(1, nrow = nrow(QR$qr), ncol = QR$rank))
    Q1 <- (1 / rw) * Q
    Q2 <- rw * Q
    d <- rowSums(Q1 * Q2)[1:n] # diag of H
    edf <- sum(d) # effective df = trace(H)
    if(!fit$family$family %in% c("poisson","binomial")) edf <- edf + 1 # estimated dispersion eats up one more df
    rdf <- n - 2 * sum((Q1[1:n,] * Q2[1:n,])) + sum(crossprod(Q1[1:n,]) * crossprod(Q2[1:n,])) # residual df = n-tr(2H-HH')
    # Note: in case sum(lambda)!=0, edf != n - rdf
  } else {
    edf <- 0
    rdf <- n
  }
  return(list(edf = edf, rdf = rdf))
}


# compute.diagH ####
# Calculate the diagonal of the hat matrix. The hat matrix is not computed
# explicitely but we take a shortcut using QR decomposition.
#
# INPUT:
#   X:    the model matrix
#   fit:  the fitted model
# OUTPUT:
#   d:    the diagonal of the hat matrix returned as a vector
compute.diagH <- function(X, w, fit){
  n <- sum(w != 0)
  lambda <- fit$lambda.w * fit$lambda
  if(sum(lambda) != 0){
    X <- rbind(X, diag(sqrt(lambda),ncol(X)))
    rw <- c(sqrt(w), rep(1, nrow(X) - length(w)))
  }
#  QR decomposition
  QR <- qr.default(X, LAPACK = TRUE)
  if (QR$rank != 0) {
    Q <- qr.qy(QR, diag(1, nrow = nrow(QR$qr), ncol = QR$rank))
    Q1 <- (1 / rw) * Q
    Q2 <- rw * Q
    d <- rowSums(Q1 * Q2)[1:n] # diag of H
  } else {
    d <- rep(n, 0)
  }
  return(d)
}

# compute.ic ####
# Function that computes a variety of information criteria
# INPUT:
#   y:      the response signal
#   y.hat:  the fitted response
#   w:      the observation weights associated to y
#   fit:    the L0glm object that predicted y.hat
# OUTPUT
#   ICs:  loglik, rss, aic, aicc, bic, ebic, hq, ric, mric, cic, bicg, bicq
compute.ic <- function(y, y.hat, w, fit){
  fam <- fit$family
  if(sum(fit$coefficients != 0) == 0){
    k <- 0
    rss <- mse <- min2LL <- aic <- bic <- aicc <- ebic <- hq <- ric <- mric <- cic <- bicg <- bicq <- Inf
    loglik <- -Inf
  } else {
    if(fam$family == "quasipoisson") fam <- poisson(identity) # Fix to let loglik be calculated for quasi families as for regular families
    n <- sum(w != 0)
    p <- length(fit$coefficients) # pn = model size

    # Get the number of parameters
    # k <- fit$df$edf
    k <- sum(fit$coefficients != 0)
    # print(k) # TODO delete

    # Residual sums of square
    rss <- compute.ss(mu = y.hat, y = y, w = w, fam = fam)

    # Mean squared error
    mse <- sum((y - y.hat)^2)/length(y)

    # TODO implement R2, adj R2, D2 here

    # Get -2 * log-likelihood
    min2LL <- fam$aic(y = y, mu = y.hat, wt = w, n = sum(w != 0),
                      dev = sum(fam$dev.resids(y = y, mu = y.hat, wt = w)))

    # Compute the information criteria
    # See Zhang & Chen 2010
    aic <- min2LL + 2*k # Akaike information criterion
    aicc <- ifelse(n > k, min2LL + 2*k*n/(n-k), NA) # corrected AIC
    bic <- min2LL + log(n)*k # Baysesian information criterion
    xi <- 1-log(n)/(2*log(p))
    ebic <- bic + 2*xi*k*log(p) # extenved BIC
    c <- 2 # TODO default, add as a parameter
    hq <- min2LL + c*log(log(n))*k # Hannan and Quinnn information criterion
    ric <- min2LL + 2 * log(p) # risk inflation criterion
    mric <- min2LL + 2 * sum(log(p/(1:k))) # modified risk inflation criterion
    cic <- min2LL + 4 * sum(log(p/(1:k))) # covariance inflation criterion
    # See Xu thesis 2010, Model Selection with Information Criteria
    g <- 1 # TODO Defautl add as argument
    bicg =  min2LL + log(n)*k + 2*g*lchoose(p,round(k))
    q <- 0.25 # TODO Defautl add as argument
    bicq <-  min2LL + log(n)*k - 2*k*log(q/(1-q))
  }
  return(c(k = k, loglik = -1/2 * min2LL, rss = rss, mse = mse, aic = aic, aicc = aicc, bic = bic, ebic = ebic,
           hq = hq, ric = ric, mric = mric, cic = cic, bicg = bicg, bicq = bicq))
}

# compute.ss ####
# Compute sums of square on the adjusted IWLS scale
# INPUT
#   mu:   the fitted response
#   y:    the response signal
#   w:    the observation weights associated to y
#   fam:  the error family used to fit mu
# OUTPUT
#   ss:   the sums of squares computed on the adjsuted IWLS scale
compute.ss <- function(mu, y, w, fam){
  eta <- fam$linkfun(mu)
  w <- w * as.vector(fam$mu.eta(eta)^2 / fam$variance(mu)) # IWLS weights from last iteration
  z.res <- (y - mu)/fam$mu.eta(eta) # the residuals on the adjusted scale = z.res = z - z.fit = eta + (y - mu)/g'(eta) - eta = (y - mu)/g'(eta)
  ss <- sum(w * z.res^2) # sums of squares
  return(ss)
}

# AR.weights ####
# adaptive ridge penalty weights to approximate L0 norm penalty:
# INPUT
#   beta:   a vector of coefficents
#   delta:  a small constant
#   gamma:  the power value
# OUTPUT
#   AR.weights: a vector with adaptive weights associated to beta
# REF
#   Zou, 2006; Candes et al., 2008 & Nuel talk
#   de Rooi & Eilers 2011, Osborne et al. 2000 & Frommlet & Nuel 2016 (eqn 5)
AR.weights <- function(beta, delta = 1E-5, gamma = 2){
  return(1/(abs(beta)^gamma + delta^gamma))
}

# AL.weights ####
# adaptive lasso penalty weights to approximate L0 norm penalty:
# INPUT
#   beta:   a vector of coefficents
#   delta:  a small constant
# OUTPUT
#   AL.weights: a vector with adaptive weights associated to beta
# REF
#   Zou, 2006; Candes et al., 2008 & Nuel talk
AL.weights <- function(beta, delta){
  return(1/(abs(beta) + delta))
}


# block.fit ####
# Function to fit the (NN)LS problem. Large data is fitted using a block
# coordinate descent algorithm, where sparse matrix (from Matrix) is supported.
# INPUT
#   y:  the response vector
#   X:  the design matrix (can be sparse, see Matrix)
#   coefs:  vector of initial coefficients
#   nonnegative:  a logical indicating whether nonnegative constraints should
#     be applied.
#   block.id: a vector indicating to which block a covariate belongs to.
#   control: a list with control parameters of the algorithm
#     block.size: amount of covariate that are fit simulatenously in each block
#                 of the coordinate descent procedure
#     maxit:  maximum number of iteration of the coordinate descent loop
#     tol:    convergence threshold. The loop stops when
#             ||coefs_{k-1} - coefs_k||_2 / ||coefs_{k-1}||_2 < tol
# OUTPUT
#  coefs:   the fitted coefficients
block.fit <- function(y, X, coefs = rep(0, ncol(X)), nonnegative = FALSE,
                      block.id = NULL, control = list(block.size = NULL,
                                                      maxit = 10, tol = 1E-7)){
  n <- length(y)
  p <- ncol(X)
  if(is.null(control$block.size) || control$block.size > p) control$block.size <- p
  if(is.null(block.id)) block.id <- ceiling((1:p)/control$block.size)

  interc.ind <- which(grepl(colnames(X), pattern = "intercept", ignore.case = TRUE))
  if(max(block.id) == 1){
    control$maxit <- 1
  } else if(length(interc.ind) != 0){ # Rearrange blocks to account for intercepts
    if(nonnegative){ # fit intercepts at the end
      block.id[interc.ind] <- max(block.id) + 1
    } else { # start with fitting intercept
      block.id <- block.id + 1
      block.id[interc.ind] <- 1
    }
  }
  blocks <- sort(unique(block.id))

  iter <- 1
  while(iter <= control$maxit){
    coefs0 <- coefs
    for(i in blocks){
      # Compute residuals
      if(length(blocks) > 1){
        A <- X[, block.id == i,drop=F]
        is.zero <- Matrix::rowSums(A) == 0
        A <- A[!is.zero,,drop=F]
        b <- y[!is.zero]
        b <- b - X[!is.zero, block.id != i, drop = F] %*% coefs[block.id != i]
      } else {
        A <- X
        b <- y
      }

      # Compute coefficients
      if(nonnegative){
        coefs[block.id == i] <- nnls(A = A, b = b)$x
      } else {
        coefs[block.id == i] <- lm.fit(x = as.matrix(A), y = as.vector(b), tol = 1E-07)$coefficients
        coefs[is.na(coefs) | is.infinite(coefs)] <- 0
      }

      # TODO delete
      # plot(coefs, type = "h", main = i)
      # Sys.sleep(0.1)
    }
    conv <- sqrt(sum((coefs0 - coefs)^2))/sqrt(sum((coefs0)^2)) < control$tol || all(coefs == 0)  # ||beta_{k-1} - beta_k||_2 / ||beta_{k-1}||_2
    if(conv) break
    iter <- iter + 1
  }
  if(control$maxit > 1 && !conv) warning("Algorithm did not converge, maybe increase control.fit$maxit or control.fit$block.size")

  return(coefs)
}


# L0glm.bfun ####
# Wrapper function around L0glm.fit to be used with boot for bootstrapping
# see L0glm.fit documentatio for more details
L0glm.bfun <- function(data, indices, X, family, lambda, start, wts, nonnegative,
                       normalize, control.l0, control.iwls, control.fit, post.filter.fn){
  beta <- L0glm.fit(y = data[indices], X = X[indices,], weights = wts[indices],
                    family = family, lambda = lambda, start = start, nonnegative = nonnegative,
                    normalize = normalize, control.l0 = control.l0, control.iwls = control.iwls,
                    control.fit = control.fit, post.filter.fn = post.filter.fn)$coefficients
  return(beta)
}

# preset.lambda ####
# Function to compute the best lambda given the IC. For some families, the
# dispersion is estimated using glmnet (+lambda tuning) which is faster than
# tuning lambda with L0glm.
preset.lambda <- function(IC, y, X, weights, family, start){
  N <- n <- nrow(X)
  # AR is an approximation of L0 penalty so the given lambda and lambda from
  # the L0 penalty objective do not exactly match (see Frommlet et al, 2016).
  # We will use the factor 4 the authors derived for orthogonal case, assuming
  # this approximately holds for any data
  lambda.hat <- 4 # The estimated true lambda
  # Lambda is dependent on the dispersion of the model
  # The issue is that the dispersion is not known beforehand, so we use glmnet
  # which is very fast to optimize but has no L0 penalty and estimate the
  # dispersion from the fit
  if(family$family %in% c("poisson", "binomial")){ # constant dispersion
    dispersion <- 1
  } else {
    # Perform an initial fit of the data
    if (!is.null(start)){
      eta <- X %*% start
      mu <- family$linkinv(eta)
    } else {
      # Initial estimate
      nobs <- n # needed for evaluating the initialization
      etastart <- mustart <- NULL # needed for evaluating the initialization
      eval(family$initialize) # generates mustart
      n <- N
      mu <- mustart
      eta <- family$linkfun(mustart)
    }
    # Transform data on the link scale
    if(family$link == "identity"){ # simplified formula and faster than formula in else statement
      W <- as.vector(weights/family$variance(mu))
      z <- y
    } else {
      mu.prime <- family$mu.eta(eta)
      z <- as.vector(eta + (y - mu) / mu.prime) # TODO Remove offset here if needed.
      W <- as.vector(weights * mu.prime^2 / family$variance(mu))
    }
    glmnet.cv.out <- cv.glmnet(x = X, y = z, weights = W, lambda = 10^seq(-10, 10, length.out = 100),
                               family = "gaussian")
    res <- y - predict(glmnet.cv.out, newx = X, s = glmnet.cv.out$lambda.min)
    dispersion <- sum(weights * res^2)/n # To avoid overhead of computing rdf or having negative variance, we use just use n instead of n-p
  }
  lambda.hat <- lambda.hat * dispersion / 2 # estimated lambda factor. This is additionally scaled to the desired IC
  lambda <- switch(IC,
                   aic = 2 * lambda.hat,
                   bic = log(n) * lambda.hat,
                   hq = 2*log(log(n)) * lambda.hat, # Hannan and Quinnn information criterion,
                   bicq = log(n) - 2*log(0.25/(1-0.25)) * lambda.hat, # TODO q = 0.25 should be a constant to optimize...
                   stop("invalid information criterion for initializing lambda. It must be one of 'aic', 'bic', 'hq', 'bicq'."))
  return(lambda)
}


##########################
#### UTILITY FUNCTION ####
##########################


#' Plot benchmark examples
#'
#' @importFrom graphics plot legend lines segments par
#'
#' @description
#'
#' Function to plot the comparison between the fitted model and the known true
#' model. It used in the \code{\link{L0glm}} examples and should not be expected
#' to useful in another use.
#'
#' @param x
#' a vector of length \code{n} with time labels.
#' @param y
#' a vector of length \code{n} with signal.
#' @param fit
#' an \code{L0glm} object
#' @param inference
#' the output of \code{L0glm.inference}
#' @param a.true
#' the true coefficients
#' @param ...
#' further arguments passed to \code{plot}
#'
#' @return
#'
#' The function returns \code{NULL}, but plots \code{y} with the expected
#' coefficients and the fitted values from \code{fit} with the estimated
#' coefficients in a way that gives an intuitive view on the quality of the fit.
#'
#' @export
plot_L0glm_benchmark <- function(x, y, fit, inference = NULL, a.true, ...){
  a.fit <- fit$coefficients
  if(!is.null(names(a.fit))){
    int.ind <- grepl("Intercept", names(a.fit))
    inference$CI.lower <- inference$CI.lower[!int.ind]
    inference$CI.upper <- inference$CI.upper[!int.ind]
    inference$p.value <- inference$p.value[!int.ind]
    a.fit <- a.fit[!int.ind]
  }
  # transf <- function(x) sqrt(x) # transform data for better rendering
  transf <- function(x) return(x)
  plot(x, transf(y), type="l", ylab="Signal", xlab="Time", col = "grey40",
       ylim=c(-transf(max(y)), transf(max(y))), ...)
  lines(x,-transf(y), col = "grey40")
  lines(transf(a.true), type = "h", col = "red", lwd = 2)
  lines(x, -transf(a.fit), type="h", col="blue", lwd = 2)
  lines(x, -transf(fit$fitted.values), col="orange2", lwd = 1)
  if(!is.null(inference)){
    # Change color of significant covariates
    sign <- inference$p.value < 0.05
    lines(x[sign], -transf(a.fit[sign]), type="h", col="green4", lwd=2)
    # Draw CIs
    segments(x0 = x, x1 = x, y0 = -transf(inference$CI.lower), y1 = -transf(inference$CI.upper), col = "grey40")
    segments(x0 = x-0.35, x1 = x+0.35, y0 = -transf(inference$CI.lower), y1 = -transf(inference$CI.lower), col = "grey40")
    segments(x0 = x-0.35, x1 = x+0.35, y0 = -transf(inference$CI.upper), y1 = -transf(inference$CI.upper), col = "grey40")
    # Legend
    legend("topleft", lty = c(1,1,0,0,0), col = c("grey40", "orange2", "red", "blue", "green2"),
           pch = c("", "", "|", "|", "|"),
           legend = c("Input signal", "Fitted signal", "Ground truth",
                      "Non significant estimates (p-value >= 0.05)", "Significant estimates (p-value < 0.05)"))
  } else {
    legend("topleft", lty = c(1,1,0,0), col = c("grey40", "orange2", "red", "blue"),
           pch = c("", "", "|", "|"),
           legend = c("Input signal", "Fitted signal", "Ground truth", "L0glm estimates"))
  }
}


# print.progress ####
# Prints progress of a loop to the console
print.progress <- function(current, total, before = "Progress: ", after = "", ...){
  if(length(current)!=1 && length(total) != 1) stop("Arguments \'current\' and \'total\' must have length 1.")
  perc <- format(current/total*100, digits = 0,...)
  cat(paste0("\r", before, perc, " %", after, "   "))
  if(current == total) cat ("\n")
  flush.console()
}


#' Simulate example data
#'
#' @description
#'
#' The function generates synthetic data consisting of a time series made of
#' Gaussian shaped spike trains, as could be observed in Gas chromatography/Mass
#' spectrometry data. The objective is to recover the spike train used to
#' generate the data by fitting a banded matrix of shifted Gaussian shapes to
#' the signal. A subset of this matrix is sampled and used to generate the
#' response signal. The signal contains Poisson noise, and hence is nonnegative
#' with integer counts.
#'
#' The function is used as a test case in the \code{\link{L0glm}} examples.
#'
#' @param n
#' the problem size.
#' @param npeaks
#' the amount of non zero coefficients. Smaller values lead to increased
#' sparsity in the coefficients. It should be smaller than \code{n}.
#' @param peakhrange
#' the range of the non zero coefficients. The coefficients are linearly
#' distributed on logarithmic scale within this range
#' @param seed
#' the seed used to sample the non zero coefficients and to generate the
#' response signal.
#' @param Plot
#' should the generated data be plotted?
#'
#' @return
#'
#' The function return a list with the following elements:
#' \item{x}{: a vector of length \code{n} containing time labels.}
#' \item{y}{: a vector of lenght \code{n} containing the response signal.}
#' \item{X}{: a \code{n} by \code{n} covariate matrix. The matrix is a banded
#' matrix with shifted Gaussian shapes from which a subset has associated non
#' zero (true) coefficients.}
#' \item{a}{: a vector of lenght \code{p = n} containing the true coefficients.
#' This is the spike train used to generate the data.}
#'
#' @export
simulate_spike_train <- function(n = 200, npeaks = 20, peakhrange = c(10,1E3),
                                 seed = 123, Plot = TRUE){
  set.seed(seed)
  x = 1:n
  # unkown peak locations
  if(npeaks > n) npeaks <- n
  u = sample(x, npeaks, replace=FALSE)
  # unknown peak heights
  h = 10^runif(npeaks, min=log10(min(peakhrange)), max=log10(max(peakhrange)))
  # locations of spikes of simulated spike train, which are assumed to be unknown
  # here, and which needs to be estimated from the measured total signal
  a = rep(0, n)
  a[u] = h
  # peak shape function
  gauspeak = function(x, u, w, h=1) h*exp(((x-u)^2)/(-2*(w^2)))
  # banded matrix with peak shape measured beforehand
  X = do.call(cbind, lapply(1:n, function (u) gauspeak(x, u=u, w=5, h=1) ))
  # colnames(X) <- paste("Var", 1:ncol(X))
  # noiseless simulated signal = linear convolution of spike train with peak shape function
  y_nonoise = as.vector(X %*% a)
  # simulated signal with random poisson noise
  y <- rpois(n, y_nonoise)

  # Plot the data
  if(Plot){
    par(mfrow=c(1,1))
    plot(y, type = "l", ylab = "Signal", xlab = "x",
         main = "Simulated spike train (red) to be estimated given known blur kernel & with Poisson noise")
    lines(a, type = "h", col = "red")
  }

  return(list(x = x, y = y, X = X, a = a))
}






