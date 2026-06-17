# quick try of HALS type GNMF with unimodality constraints
# see also Li et al. 2012 with HALS for Bregman divergence NMF using Taylor expansion
# Cichocki et al. 2008 HALS type algos under beta divergence
# didn't work though

# what would probably work though is:
# initialisation of algo as in L0BestRankGNMF with grouped_Hupdate=T

# --- rescale row sums of H to unit L2 norm & jointly rescale W -------------
s <- sqrt(rowSums(H^2))
H <- sweep(H, 1, s, FUN="/")
W <- sweep(W, 2, s, FUN="*")

# cluster factors with high overlap & join (sum) their elution profiles
# and upon joining only keep compound spectrum in H with highest
# peak height (colMax(W)) (or average them)
# (maybe based on overlap in outer products between outer(W[,k], H[,k]))
# where k=compound, maybe using Xmeans clustering, but not sure of 
# exact distance metric to use - given separability, it may also be possible
# to do two separate clusterings in function of overlap in elution profiles
# and in function of spectral similarity & join compounds only if in both
# clusterings they appear in the same group)
# Then make new W based on aggregating (summing) elution profiles in W within clusters

# update k

# loop until convergence
# force log-concavity on columns of W using log-concave spline fit, fit using cobs
# update objective

# update H via multivariate group L0 penalized GLM / L0glm.fit using starting values provided
# by updated H with updated k & log-concavity constraints

# update W via L0 penalized GLM / L0glm.fit using starting values
# provided by previous iteration W






# a component-wise rank-1 deflation approach maybe also possible:

# helper functions
dev_fun <- function(mu) {
  sum(family$dev.resids(as.vector(D), as.vector(mu), rep_len(1, length(D))))
}

var_fun <- family$variance           # Var(Y|µ)

# calculate
mu <- W %*% H + eps               # current mean
B  <- 1 / pmax(var_fun(mu), eps)  # n × p 1/variance working weights

# sort factors in W and H from most to least important based on colMax(W)

# set remaining residual data matrix Dres = D, i.e. all signal remains to be fitted

# set remaining_factors to be fitted equal to 1:k (with most important one now first)
# remaining_factors = 1:k

# set residual matrix with signal that remains to be fitted Dres = D

# loop over remaining elution profiles to locally update elution profiles:
# while (length(remaining_factors)>0) do
# k = 1 # largest remaining factor
# join together identical compounds :
# iteratively decide which other compounds should belong together with the
# first largest peak & whose elution profiles should be added together
# (e.g. for those with broader peaks than the assumed PSF) - do this gradually starting from
# the closest peaks (i.e. with overlapping elution profile in the columns of
# W) when (1) elution profiles overlap, (2) they have very high spectral similarity 
# (based on cosine similarity between spectra in rows of H) and (3) adding together the
# elution profiles would retain approximate log-concavity
# each time compound elution profiles are joined/summed update W and remaining_factors
# accordingly & adjust the matrix W accordingly (remove the rows that proved to belong
# to the same factor than the 1st largest one)

# calculate for the local window=which(W[,1,drop=F]>0):
#  notk = 2:length(remaining_factors)
#  Dk = D[window,,drop=F]-W[window,notk,drop=F] %*% H[notk,,drop=F] 
#  (local window in data matrix deflated by all other factors than the one to be considered)
# update 

# outer loop
#   inner loop looping over k compounds in order of importance
#       - calculate residual-(k-1) data matrix D(k-1) with contributions of all compounds except kth one removed
#       - update kth column of H with univariate quantile regressions that are constrained to
#       be log-concave via a cobs log-concave spline fit on obtained coefficients
#       (this can be fit on local window of D corresponding to nonzero entries in column of H
#       but window can then be dynamically expanded to region where predictions of cobs
#       fit remain above small threshold)

# cluster columns of H with high correlation & make new H based on aggregating
# (averaging) clusters, maybe using Xmeans clustering
# also calculate 

# update W by regressing D on H with grouped penalty 

# repeat until convergence


# set Hnew = H
# loop over k compounds
#   update columns of Hnew with 
#    calculating residual-(k-1) data matrix D(k-1) with contributions of all compounds except k removed
#    for positive part of data matrix D(k-1) 


irls_hals_gnmf <- function(D, rank = NULL,
                           family     = poisson("identity"),
                           max_outer  = 50,
                           max_inner  = 3,
                           tol        = 1e-4,
                           Wstart     = NULL,
                           unimodal   = FALSE,
                           eps        = 1e-10) {
  if (!is.matrix(D)) D <- as.matrix(D)
  n <- nrow(D); p <- ncol(D)
  if (any(D < 0)) stop("D must be non‑negative")
  
  if (unimodal && !requireNamespace("cobs", quietly = TRUE))
    stop("Package 'cobs' is required for unimodality projection")
  if (!requireNamespace("NMF",  quietly = TRUE))
    stop("Package 'NMF' is required for NNDSVD initialisation")
  
  # ---- initial W, H ---------------------------------------------------------
  if (is.null(Wstart)) {
    init <- suppressWarnings(NMF:::.nndsvd.internal(A = D, k = rank))
    W <- pmax(init$W, eps)
    H <- pmax(init$H, eps)
    rank <- ncol(W)
  } else {
    rank <- ncol(W)
    W <- pmax(Wstart, eps)
    max_indices <- apply(W, 2, which.max) # initialise H with data matrix at peak maxima
    H <- Matrix(D[max_indices,], sparse=TRUE) 
    # H <- matrix(runif(rank * p), rank, p)  # random non‑negative
  }
  
  # ---- helpers --------------------------------------------------------------
  dev_fun <- function(mu) {
    sum(family$dev.resids(as.vector(D), as.vector(mu), rep_len(1, length(D))))
  }
  
  var_fun <- family$variance           # Var(Y|µ)
  linkinv <- family$linkinv            # should be identity ⇒ µ = η
  
  dev_old <- dev_fun(W %*% H + eps)
  
  # ---- IRLS outer loop ------------------------------------------------------
  for (outer in 1:max_outer) {
    mu      <- W %*% H + eps          # current mean
    B       <- 1 / pmax(var_fun(mu), eps)   # n x p working weights with 1/Var

    # ---- HALS inner sweeps -------------------------------------------------
    for (inner in 1:max_inner) {
      for (k in 1:rank) {
        # residual after taking out component k-------------------------------
        Rk <- D - mu + tcrossprod(W[, k], H[k, ])
        
        # -- update H row ----------------------------------------------------
        numer_h <- crossprod(W[, k],  B * Rk)            # 1 × p
        denom_h <- crossprod(W[, k]^2, B)                # 1 × p
        H[k, ]   <- pmax(0, as.numeric(numer_h) / (as.numeric(denom_h) + eps))
        
        # -- update W column --------------------------------------------------
        Rk <- D - mu + tcrossprod(W[, k], H[k, ])        # refresh residual
        numer_w <- (B * Rk) %*% H[k, ]                   # n × 1
        denom_w <- (B %*% (H[k, ]^2))                    # n × 1
        W[, k]  <- pmax(0, as.numeric(numer_w) / (as.numeric(denom_w) + eps))
        
        # -- optional unimodality projection ---------------------------------
        if (unimodal && any(W[, k] > 0)) {
          wlog <- log(pmax(W[, k], eps))
          fit  <- cobs::cobs(x = seq_len(n), y = wlog,
                             constraint = "concave", lambda = 0,
                             print.warn = FALSE, print.mesg = FALSE)
          W[, k]  <- pmax(0, exp(fit$fitted))
          Ww[, k] <- sqrtB * W[, k]
        }
      }
    }
    
    # ---- deviance & convergence --------------------------------------------
    mu_new <- W %*% H + eps
    dev_new <- dev_fun(mu_new)
    if (abs(dev_new - dev_old) / (abs(dev_old) + eps) < tol) break
    dev_old <- dev_new
  }
  
  list(W = W, H = H, dev = dev_old, iter = outer)
}




library(Matrix)
# devtools::load_all(".")
# library(L0glm)

n <- 100 # timepoints
k_true <- 10 # true rank (nr of compounds)
peak_width <- 2.5 # width / sigma of Gaussian peak shapes 
psf_fun <- function(x, u, w=peak_width) { exp(((x-u)^2)/(-2*(w^2))) } # gaussian point spread / peak shape function
psf <- psf_fun(x=0:20, u=10) # vector with PSF
plot(psf, type="o")
sim.pois.mc <- simulate_spike_train(n = n, k = k_true, family = "poisson", 
                                    seed = 123,
                                    multichannel = TRUE,
                                    mean_beta = 10000,
                                    sparse = TRUE,
                                    sd_logbeta = 1,
                                    psf = psf_fun)
D = sim.pois.mc$y # n x p with n timepoints and p mass channels
H_true = sim.pois.mc$beta_true[rowSums(sim.pois.mc$beta_true)!=0,,drop=F]
image(D)
plot(rowSums(D), type="l")
image(H_true)

# initialise elution profiles (W matrix)
# we initialise W with an n x n banded matrix with 
# time-shifted copies of the point spread / peak shape function (PSF) 
# (this can be estimated by running standards or tuned based on the IC of different L0glm fits 
# with banded matrices with different peak widths)
Wstart = sim.pois.mc$X  

# this matrix can subsequently be pre-subsetted by regressing the row sums of D
# (=total ion chromatogram) on our banded covariate matrix Wstart using an
# L0-penalised nonnegative identity link Poisson fit, fit using L0glm.fit  and with
# lambda=0, after which we can retain only the columns of Wstart with nonzero
# coefficients and use this as initialisation for W
keep <- (L0glm.fit(X = Wstart, 
                   y = rowSums(D), # n x p data matrix, now regressing on row sums of D (=total ion chromatogram)
                   intercept = FALSE,
                   family = poisson(link=identity),
                   lambda = 0,  
                   nonnegative = TRUE,
                   solver = "osqp" # using unregularised nonnegative identity link Poisson fit using osqp solver and with lambda=0
)$coefficients != 0)
sum(keep) # 21
Wstart_rowsums = as.matrix(Wstart[, keep, drop=F])
Wstart = Wstart_rowsums

IC = "bic" # information criterion to be optimised, we choose BIC since after
# initial variable preselection we are working in a low-dimensional regime


uninmf <- irls_hals_gnmf(D, family = poisson("identity"), 
                             max_outer  = 50,
                             max_inner  = 3,
                             tol = 1e-4,
                             Wstart = Wstart)

dim(W) # 100 x 21  
matplot(W, type="l") # all profiles are identical though
dim(H) # 21 x 132
image(H) # all rows are identical though






uninmf <- hals_unimodal_gnmf(D, family = poisson("identity"), 
                               max_iter = 50, 
                               tol = 1e-4,
                               IC = "bic",
                               Wstart = Wstart)

dim(W) # 100 x 21  
matplot(W, type="l") # all profiles are identical though
dim(H) # 21 x 132
image(H) # all rows are identical though







#### pickLogConPeak ####
# Function that approximate the given signal (y) by a log-concave smooth signal. The
# function uses the cobs algorithm which is based on quatile regression
# INPUT
#   y:    the signal (vector) to approximate
#   tau:  the quantile of the signal to fit
#   weights:    vector of weights for the quantile regression
#   local.fit:  a logical indicating if the fitting range be bound by zero values
#   mode: a numeric giving the index of the expected mode. Used only when local.fit==TRUE.
#         The fit is performed on the part of the signal that contains the mode.
#   basfun: a function that fits a baseline that is subtracted before fitting.
#   Plot: a logical indicating if the fitting should be plotted
# OUTPUT
#   y.new:  the fitted signal
# TODO make this a generic function to fit different kind of peak shape constraints for instance: gaussian, exponetially modified gaussian, or any user defined paramteric shape
pickLogConPeak <- function(y, tau = 0.1, weights = NULL, local.fit = TRUE,
                           mode = NULL, basfun = function(x)return(0),
                           Plot = FALSE, ...){
  y <- pmax(y-basfun(y),0)
  y.new <- rep(0, length(y))
  if(sum(y) == 0) return(y.new)
  
  # Local Fit
  if(local.fit){
    starts <- pmax(which(diff(sign(c(0,y))) == 1),1) # before: 0, after: pos value
    stops <- pmin(which(diff(sign(c(y,0))) == -1),length(y)) # before: pos value, after: 0
    if(is.null(mode)){ # if no mode is supplied, take the chunk that has highest norm
      norm.ysub <- sapply(1:length(starts), function(l){
        return(sqrt(sum(y[starts[l]:stops[l]]^2)))
      })
      start <- starts[which.max(norm.ysub)]
      stop <- stops[which.max(norm.ysub)]
      
    } else {
      if(y[mode] == 0) return(y.new)
      start <- tail(starts[starts <= mode],1)
      stop <- head(stops[stops >= mode],1)
      if(length(start) * length(stop) == 0) return(y.new)
    }
  } else {
    start <- 1
    stop <- length(y)
  }
  y.sub <- y[start:stop]
  x <- start:stop
  
  # Weights
  if(is.null(weights)){
    offs <- max(y.sub)*1E-4
    weights <- sqrt(y.sub+offs)
  } else {
    weights <- weights[start:stop]
  }
  weights <- weights/sum(weights)*length(weights)
  
  # Knots
  nknots <- round(length(y.sub)/3)
  if(nknots <= 3) return(y.new)
  
  # Log-concave fit using cobs
  require(cobs)
  offs <- max(y.sub)*1E-6
  fit <- suppressWarnings(cobs(x=x,y=log10(y.sub+offs),
                               constraint= "concave", lambda=0,
                               knots = seq(min(x),max(x),length.out = nknots),
                               nknots=nknots, knots.add = FALSE, repeat.delete.add = FALSE,
                               w=weights,
                               tau=tau, print.warn = F, print.mesg = F,
                               ...))
  y.fit <- fit$fitted
  y.fit <- pmax(10^y.fit - offs,0)
  y.fit[is.infinite(y.fit)] <- 0 # Happens when cobs does not converge and that the fit is exponentiated
  y.new[start:stop] <- y.fit
  
  if(Plot){
    plot(y+1, log = "y", type = "l")
    lines(y.new+1, col = 4)
    lines(x = start:stop, y = weights*max(y.sub)/max(weights)+1, col = "grey")
    abline(v = c(start, stop), lty = 2)
    abline(v = mode, col = 2)
  }
  return(y.new)
}