#' L0-Penalised Generalised Best-Rank NMF
#'
#' Factorises a non-negative data matrix \eqn{D} into \eqn{W\,H} under any
#' exponential-family noise model while **simultaneously selecting the
#' best rank** via hard (L0 or grouped L0) best-subset penalties. This NMF
#' can for example be used to separate (deconvolute) overlapping peaks in
#' GC/MS (gas chromatography - mass spectrometry) datasets, allowing one
#' to simultaneously estimate compound elution profiles (W), electron impact
#' mass spectra (H) and the true number of compounds that are present in the
#' chromatogram (the best rank).
#'
#' The routine alternates two calls to `L0glm.fit()` and works via cyclical 
#' block coordinate descent:
#' \describe{
#'   \item{H-update}{multivariate GLM with L0 or optional **grouped** L0 penalty
#'                   (selects or discards an entire latent component
#'                   across all outcome channels).}
#'   \item{W-update}{plain L0 GLM row-wise, warm-started from the previous
#'                   iterate.}
#' }
#' After each full sweep factors that are shrunk exactly to 0 are pruned,
#' which means one can start with \code{k} larger than the true rank.
#' The optimal best rank is then typically very well recovered.
#'
#' The *penalised objective* minimised is
#' \deqn{
#'   \sum_{ij}\frac12\,D\_{F}\!\bigl(D\_{ij}\,\|\,\mu\_{ij}\bigr) \;+\;
#'   \lambda_{W}\!\sum_{k}\mathbf 1[\|W_{\cdot k}\|_2>0] \;+\;
#'   \lambda_{H}\!\sum_{k}\mathbf 1[\|H_{k\cdot}\|_2>0],
#' }
#' where \eqn{D\_{F}} is the Bregman (half-deviance) loss of the chosen
#' family and \eqn{\mu = W\,H}.
#'
#' @section Implementation notes:
#' * Dispersion \eqn{\phi} is re-estimated after every sweep for
#'   Gaussian, Gamma and inverse-Gaussian models (or fixed to 1 for Poisson).
#' * Rows of \eqn{H} are rescaled to unit \eqn{L^2} norm at the end;
#'   corresponding columns of \eqn{W} are scaled back to keep
#'   \eqn{W\,H} unchanged.
#'
#' @param D      Numeric (or \link[Matrix]{Matrix}) \eqn{n\times p} input
#'               matrix with non-negative entries.
#' @param k      Initial rank of the non-negative double SVD (NNDSVD) that is
#'               used for initialisation if \code{Wstart=NULL}. Allowed to 
#'               over-estimate the true rank.
#' @param Wstart Optional \eqn{n\times k_0} non-negative matrix with
#'               starting values for \eqn{W}.  If supplied, \code{k} is
#'               ignored and set to \code{ncol(Wstart)}.
#' @param grouped_Hupdate Logical; if \code{TRUE} a **grouped** best-subset
#'               penalty is used in the \eqn{H} update (one latent
#'               component on/off for all outcome channels, i.e. best subset of
#'               group selection).
#' @param IC     Information criterion used to set the penalty size
#'               (\code{"bic"}, \code{"gic"}, \code{"ebic"}, …); passed to
#'               \code{preset.lambda()}. Defaults to \code{"gic"}, which works
#'               well for both low- and high-dimensional settings.
#' @param max_iter Maximum number of outer iterations.
#' @param tol    Relative change in penalised objective to determine convergence.
#' @param method  Eigen solver to use, as in `L0glm.fit()` (e.g. 4L=least squares
#'               conjugate gradient, 7L=Cholesky LDL).
#' @param parallelize  Use parallelisation? May only be faster for larger problems.
#' @param ncores  Number of cores to use when \code{"parallelize=TRUE"}.
#' @param family A \link[stats]{family} object defining the exponential
#'               family and link; identity-link Poisson provides
#'               Lee-Seung KL-NMF but with additional L0 or grouped L0
#'               regularisation.
#'
#' @return A list with components
#' \describe{
#'   \item{\code{W}}{Estimated \eqn{n\times \hat k} W matrix.}
#'   \item{\code{H}}{Estimated \eqn{\hat k\times p} H matrix.}
#'   \item{\code{dev}}{Final penalised objective value, i.e. Bregman (half-deviance) 
#'   loss plus imposed L0 or grouped L0 penalties).}
#'   \item{\code{k}}{Discovered rank \eqn{\hat k}.}
#'   \item{\code{iter}}{Number of outer iterations run.}
#'   \item{\code{dev_trace}}{Vector of penalised objective values over iterations.}
#'   \item{\code{lambdaW}}{Lambda penalty used for W.}
#'   \item{\code{lambdaH}}{Lambda penalty used for H.}
#'   
#' }
#'
#' @examples
#' ## ---------- GC/MS deconvolution toy example with Poisson noise ------------
#' if (requireNamespace("Matrix", quietly = TRUE)) {
#'   library(Matrix)
#' }
#' n  <- 200                 # time points
#' k0 <- 50                  # true rank/ number of compounds
#' set.seed(123)
#'
#' peak_width <- 2.5         # unknown, but can be estimated from standards or tuned from chromatograms
#' psf_fun <- function(x,u,w = peak_width) exp(-0.5*((x-u)/w)^2)
#'
#' sim <- simulate_spike_train(n       = n,
#'                             k       = k0,
#'                             family  = "poisson",
#'                             multichannel = TRUE,
#'                             mean_beta    = 1e4,
#'                             sparse       = TRUE,
#'                             psf          = psf_fun)
#'
#' D <- sim$y                # observed GC/MS chromatogram
#'
#' ## build a design matrix consisting of time-shifted copies of the peak shape function (PSF)
#' Wband <- sim$X            # n × n banded matrix with time-shifted copies of PSF
#' ## pre-select putative peaks using unregularised nonnegative identity link Poisson fit on total ion chromatogram
#' sel   <- L0glm.fit(Wband, rowSums(D),
#'                    lambda = 0, nonnegative = TRUE,
#'                    family = poisson(link=identity),
#'                    solver = "osqp")$coefficients != 0
#' Winit <- as.matrix(Wband[, sel, drop = FALSE])
#'
#' ## sparse Poisson Best-Rank NMF with grouped-L0 penalty to infer mass spectra, elution profiles and number of compounds present
#' nmf <- L0BestRankGNMF(D,
#'                       Wstart          = Matrix(Winit, sparse = TRUE),
#'                       grouped_Hupdate = TRUE,
#'                       IC              = "gic",
#'                       max_iter        = 5,
#'                       family          = poisson("identity"))
#'
#' cat(" True rank (nr. of compounds present):", k0, "\n", "Rank recovered by Poisson L0BestRankGNMF:", nmf$k, "\n")
#' matplot(nmf$W, type="l", xlab="Time", ylab="Signal", 
#' main="Elution profiles inferred by Poisson L0BestRankGNMF")
#'
#' ## ---------- Comparison with MCR-ALS (alternating non-negative least squares) -
#' if (requireNamespace("ALS", quietly = TRUE)) {
#'   library(ALS)
#'   als_fit <- als(CList   = list(Winit),
#'                  S       = matrix(1, ncol(D), ncol(Winit)),
#'                  PsiList = list(as.matrix(D)),
#'                  uniC    = TRUE,  # impose unimodality constraints
#'                  maxiter = 8,
#'                  normS   = 1)
#'
#'   cat(" Rank recovered by MCR-ALS:", nrow(t(als_fit$S)), "\n")
#' }
#'
#' @import Matrix
#' @export
L0BestRankGNMF <- function(D, k, Wstart = NULL,
                     grouped_Hupdate = FALSE, # use plain L0 or grouped L0 penalty in H update?
                     IC = "gic", # information criterion to be optimised
                     max_iter = 5, tol = 1e-4, 
                     method = 4L, parallelize = FALSE,
                     ncores = if (parallelize) { round(parallel::detectCores()/2) } else { 1 },
                     family = gaussian(link="identity")) {
  
  D <- Matrix(D, sparse = TRUE)
  Dt <- t(D)
  n <- nrow(D); p <- ncol(D)
  if (!is.null(Wstart)) { k <- ncol(Wstart) }
  
  ## --- translate target IC to numeric lambda values --------------------------
  lambdaH <- preset.lambda(IC = IC, p = k, n = n)
  lambdaW <- preset.lambda(IC = IC, p = k, n = p)
  
  ## --- initialise W and H ----------------------------------------------------
  if (!is.null(Wstart)) { W <- Wstart # n x k_init matrix
         suppressWarnings(H <- L0glm.fit(X = W, 
                                         y = D, # n x p data matrix, now regressing on columns of D
                                         intercept = FALSE,
                                         family = family,
                                         lambda = lambdaH,
                                         nonnegative = TRUE,
                                         parallelize = parallelize, 
                                         ncores = ncores,
                                         maxit = 10, 
                                         group.maxit = 3, 
                                         grouped = grouped_Hupdate, # apply L0 or grouped L0 penalty
                                         relaxed = !grouped_Hupdate, # refit selected groups of coefficients unpenalized at the end if grouped=T
                                         solver = "eigen", # method 4L: least squares conjugate gradient
                                         method = method)$coefficients) # k_init x p matrix
  } else {
     # NNDSVD initialisation
     suppressWarnings(nndsvd <- NMF:::.nndsvd.internal(A=D, k=k))
     W <- Matrix(nndsvd$W, sparse=T) # n x k
     H <- Matrix(nndsvd$H, sparse=T) # k x p
  } 
  
  # prune empty factors
  keep = which(rowSums(H)!=0 & colSums(W)!=0)
  H_old <- H <- H[keep, , drop=F] 
  W_old <- W <- W[, keep, drop=F] 
  k <- nrow(H)
  
  ## --- calculate objective ---------------------------------------------------
  obj <- penalised_nll(W, H, D, family, lambdaW, lambdaH)
  trace <- obj
  
  for (iter in seq_len(max_iter)) {
    ## --- update W -------------------------------------------------------
    suppressWarnings(W <- t(L0glm.fit(
                              X = t(H), # p x k matrix
                              y = Dt, # p x n data matrix, now regressing on rows of D
                              intercept = FALSE,
                              start = t(W), # use warm starts to speed up algo
                              family = family,
                              lambda = lambdaW,
                              nonnegative = TRUE,
                              parallelize = parallelize,
                              maxit = 10, 
                              grouped = FALSE,
                              solver = "eigen", # method 4L: least squares conjugate gradient
                              method = method)$coefficients)) # n x k matrix
    if (sum(range(W) == 0)==2) { 
      message("All factors were penalized to zero in the NMF")
      return(list(W = NULL, H = NULL, dev = obj, k = 0, iter = iter, 
                  dev_trace = trace, lambdaW = lambdaW, lambdaH = lambdaH)) }
    
    ## --- update H -------------------------------------------------------
    suppressWarnings(H <- L0glm.fit(
                   X = W, # n x k matrix
                   y = D, # n x p data matrix, now regressing on columns of D
                   intercept = FALSE,
                   start = H, # use warm starts to speed up algo
                   family = family,
                   lambda = lambdaH,
                   nonnegative = TRUE,
                   parallelize = parallelize,
                   ncores = ncores,
                   maxit = 10, 
                   group.maxit = 3, 
                   grouped = grouped_Hupdate, # apply L0 or grouped L0 penalty
                   relaxed = !grouped_Hupdate, # refit selected groups of coefficients unpenalized at the end if grouped=T
                   solver = "eigen", # method 4L: least squares conjugate gradient
                   method = method)$coefficients) # k x p matrix
    
    ## --- prune dead factors --------------------------------------------
    keep <- which(rowSums(H)!=0 & colSums(W)!=0) # only keep factors with nonzero row sums in H & nonzero column sums in W
    H <- H[keep, , drop=FALSE] 
    W <- W[, keep, drop=FALSE] 
    k <- nrow(H)

    ## --- convergence check---------------------------------------------
    obj_new <- penalised_nll(W, H, D, family, lambdaW, lambdaH)
    if (obj_new > obj) { H <- H_old
                         W <- W_old
                         k <- nrow(H)
                         obj_new <- obj } # fit deteriorated so return previous iteration estimates
    trace <- c(trace, obj_new)
    if (abs(obj_new - obj)/(abs(obj) + 1e-12) < tol) break
    obj <- obj_new
  }
  
  ## --- rescale row sums of H to unit L2 norm --------------------------
  s <- sqrt(rowSums(H^2))
  H <- sweep(H, 1, s, FUN="/")
  W <- sweep(W, 2, s, FUN="*")
  
  list(W = W, H = H, dev = obj, k = k, iter = iter, dev_trace = trace, lambdaW = lambdaW, lambdaH = lambdaH)
}

# a few internal helper functions (not exported):

# un-penalised negative log-likelihood objective = summed Bregman divergence
# (beta divergence is a special case of Bregman divergence, and for beta divergence
# squared Euclidean, Kullback-Leibler, Itakura-Saito and beta=-1 divergences are
# special cases that are often used in the NMF literature to deal with
# Gaussian, Poison, Gamma or inverse Gaussian noise in NMF problems, 
# but Bregman divergence is more general & can readily be calculated from
# the inbuilt deviance residuals). 
neg_loglik <- function(D, mu, fam, disp = 1) {
  dev <- fam$dev.resids(as.vector(D), 
                        as.vector(mu), 
                        wt = rep_len(1, length(D)))
  good <- !is.infinite(dev)
  return(sum(dev[good]) / (2 * disp))
}

# L0 penalized negative log-likelihood objective
penalised_nll <- function(W, H, D, fam, lambdaW, lambdaH) {
  mu  <- W %*% H
  disp <- estimate_phi(D, mu, fam, # estimate dispersion if need be
                       df = sum(colSums(W)!=0) + sum(rowSums(H)!=0))
  loss <- neg_loglik(D, mu, fam, disp)
  pen  <- lambdaW * sum(colSums(W) > 0) + lambdaH * sum(rowSums(H) > 0)
  return(loss + pen)
}

# helper function to estimate dispersion coefficient phi
estimate_phi <- function(D, mu, fam, df) {
  if (fam$family %in% c("binomial","poisson")) { return(1) } else {
    dev <- fam$dev.resids(as.vector(D), as.vector(mu), wt = rep_len(1, length(D)))
    phi <- sum(dev) / (2 * (length(dev) - df))
    return(max(phi, 1e-12)) }             # avoid division by 0 later
}