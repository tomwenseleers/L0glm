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
#' The objective minimised after each full W-H sweep is
#' \deqn{
#'   \underbrace{\frac{1}{\phi}\sum_{ij} D_{F}\!\bigl(D_{ij}\,\|\,\mu_{ij}\bigr)}
#'   _{-2\,\mathrm{logLik}}
#'   \;+\;
#'   \lambda_W \sum_{k}\mathbf 1[\|W_{\cdot k}\|_2>0]
#'   \;+\;
#'   \begin{cases}
#'     \lambda_H\sum_{k}\mathbf 1[\|H_{k\cdot}\|_2>0],
#'         &\text{if } grouped\_Hupdate=TRUE\\[6pt]
#'     \lambda_H\!\sum_{k,j}\mathbf 1[H_{kj}\neq 0],
#'         &\text{otherwise.}
#'   \end{cases}
#' }
#'
#' The grouped penalty is a group-\eqn{L_0}/row-\eqn{L_0} penalty, equivalently
#' an \eqn{L_{2,0}} penalty, counting whole latent components as on/off. The
#' ungrouped penalty is the plain elementwise \eqn{L_0} penalty used by
#' \code{L0glm.fit}.
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
#' @inheritParams L0glm.tune
#' @param max_iter Maximum number of outer iterations.
#' @param tol    Relative change in penalised objective to determine convergence.
#' @param method  Eigen solver to use, as in `L0glm.fit()` (e.g. 1L=Cholesky
#'               LLT, 4L=least squares conjugate gradient).
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
#'   \item{\code{dev}}{Final penalised objective value, i.e. scaled Bregman (deviance) 
#'   loss (=-2*log-likelihood), with added elementwise \eqn{L_0} or grouped
#'   row-\eqn{L_0}/\eqn{L_{2,0}} penalties).}
#'   \item{\code{k}}{Discovered rank \eqn{\hat k}.}
#'   \item{\code{iter}}{Number of outer iterations run.}
#'   \item{\code{dev_trace}}{Vector of penalised objective values over iterations.}
#'   \item{\code{lambdaW}}{Lambda penalty used for W.}
#'   \item{\code{lambdaH}}{Lambda penalty used for H.}
#'   
#' }
#'
#' @example examples/L0BestRankGNMF.R
#' @import Matrix
#' @export
L0BestRankGNMF <- function(D, k, Wstart = NULL,
                           grouped_Hupdate = FALSE, # use plain L0 or grouped L0 penalty in H update?
                           IC = "gic", # information criterion to be optimised
                           max_iter = 5, tol = 1e-4, 
                           method = 1L, parallelize = FALSE,
                           ncores = if (parallelize) { round(parallel::detectCores()/2) } else { 1 },
                           family = gaussian(link="identity")) {
  
  # D <- as.matrix(D) # coerce to dense
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
    nndsvd <- nndsvd_init(D, k)
    W <- Matrix(nndsvd$W, sparse=T) # n x k
    H <- Matrix(nndsvd$H, sparse=T) # k x p
  } 
  
  # prune empty factors
  keep = which(active_rows(H) & active_cols(W))
  H_old <- H <- H[keep, , drop=F] 
  W_old <- W <- W[, keep, drop=F] 
  k <- nrow(H)
  
  ## --- calculate objective ---------------------------------------------------
  obj <- penalised_nll(W, H, D, family, lambdaW, lambdaH, grouped_Hupdate)
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
    keep <- which(active_rows(H) & active_cols(W)) # keep factors with nonzero rows in H and nonzero columns in W
    H <- H[keep, , drop=FALSE] 
    W <- W[, keep, drop=FALSE] 
    k <- nrow(H)
    
    ## --- convergence check---------------------------------------------
    obj_new <- penalised_nll(W, H, D, family, lambdaW, lambdaH, grouped_Hupdate)
    if (obj_new > obj) { H <- H_old
    W <- W_old
    k <- nrow(H)
    obj_new <- obj } # fit deteriorated so return previous iteration estimates
    trace <- c(trace, obj_new)
    if (abs(obj_new - obj)/(abs(obj) + 1e-12) < tol) break
    obj <- obj_new
  }
  
  ## --- rescale row sums of H to unit L2 norm & jointly rescale W -------------
  s <- sqrt(rowSums(H^2))
  H <- sweep(H, 1, s, FUN="/")
  W <- sweep(W, 2, s, FUN="*")
  
  list(W = W, H = H, dev = obj, k = k, iter = iter, dev_trace = trace, lambdaW = lambdaW, lambdaH = lambdaH)
}

# a few internal helper functions (not exported):

# un-penalised -2*negative log-likelihood objective = scaled Bregman / deviance loss
# (i.e. -2*log-likelihood)
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
  return(sum(dev[good]) / disp) # no /2 to get L0 penalties on same scale as in L0glm.fit, i.e. use -2*LL+lambda*L0 penalty
}

# L0 penalized negative log-likelihood objective
penalised_nll <- function(W, H, D, fam, lambdaW, lambdaH, grouped_Hupdate) {
  mu  <- W %*% H
  disp <- estimate_phi(D, mu, fam,
                       df = sum(active_cols(W)) +
                         if (grouped_Hupdate)
                           sum(active_rows(H)) else sum(H!=0))
  loss <- neg_loglik(D, mu, fam, disp)
  
  penW <- lambdaW * sum(active_cols(W))
  penH <- if (grouped_Hupdate)
    lambdaH * sum(active_rows(H))      # row-L0 / L2,0
  else
    lambdaH * sum(H != 0)              # elementwise L0
  return(loss + penW + penH)
}

active_rows <- function(M) {
  as.vector(Matrix::rowSums(abs(M) > 0) > 0)
}

active_cols <- function(M) {
  as.vector(Matrix::colSums(abs(M) > 0) > 0)
}

# helper function to estimate dispersion coefficient phi
estimate_phi <- function(D, mu, fam, df) {
  if (fam$family %in% c("binomial","poisson")) { return(1) } else {
    dev <- fam$dev.resids(as.vector(D), as.vector(mu), wt = rep_len(1, length(D)))
    phi <- sum(dev) / (2 * (length(dev) - df))
    return(max(phi, 1e-12)) }             # avoid division by 0 later
}

nndsvd_init <- function(A, k) {
  A <- as.matrix(A)
  if (any(A < 0, na.rm = TRUE)) {
    stop("NNDSVD initialisation requires non-negative data.", call. = FALSE)
  }
  nr <- nrow(A)
  nc <- ncol(A)
  rank <- min(k, nr, nc)
  decomp <- svd(A, nu = rank, nv = rank)
  W <- matrix(0, nr, k)
  H <- matrix(0, k, nc)
  if (rank == 0L || length(decomp$d) == 0L) {
    return(list(W = W, H = H))
  }

  W[, 1L] <- sqrt(decomp$d[1L]) * abs(decomp$u[, 1L])
  H[1L, ] <- sqrt(decomp$d[1L]) * abs(decomp$v[, 1L])
  if (rank == 1L) {
    return(list(W = W, H = H))
  }

  for (j in 2L:rank) {
    u <- decomp$u[, j]
    v <- decomp$v[, j]
    u_pos <- pmax(u, 0)
    u_neg <- pmax(-u, 0)
    v_pos <- pmax(v, 0)
    v_neg <- pmax(-v, 0)

    u_pos_norm <- sqrt(sum(u_pos^2))
    u_neg_norm <- sqrt(sum(u_neg^2))
    v_pos_norm <- sqrt(sum(v_pos^2))
    v_neg_norm <- sqrt(sum(v_neg^2))
    pos_norm <- u_pos_norm * v_pos_norm
    neg_norm <- u_neg_norm * v_neg_norm

    if (pos_norm >= neg_norm) {
      basis_u <- if (u_pos_norm > 0) u_pos / u_pos_norm else u_pos
      basis_v <- if (v_pos_norm > 0) v_pos / v_pos_norm else v_pos
      sigma <- pos_norm
    } else {
      basis_u <- if (u_neg_norm > 0) u_neg / u_neg_norm else u_neg
      basis_v <- if (v_neg_norm > 0) v_neg / v_neg_norm else v_neg
      sigma <- neg_norm
    }
    scale <- sqrt(decomp$d[j] * sigma)
    W[, j] <- scale * basis_u
    H[j, ] <- scale * basis_v
  }
  list(W = W, H = H)
}
