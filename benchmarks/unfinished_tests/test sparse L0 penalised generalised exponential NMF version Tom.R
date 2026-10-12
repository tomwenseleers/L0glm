# L0 or grouped L0 penalised generalised exponential-amily NMF with automatic best rank discovery
# with applications for the deconvolution of overlapping peaks in GC/MS chromatography

# a few helper functions (not exported):

# un-penalised negative log-likelihood objective = summed Bregman divergence
# (beta divergence is a special case of Bregman divergence, and for beta divergence
# squared Euclidean, Kullback-Leibler, Itakura-Saito and beta=-1 divergences are
# special cases that are often used in the NMF literature to deal with
# Gaussian, Poison, Gamma or inverse Gaussian noise in NMF problems, 
# but Bregman divergence is more general & can readily be calculated from
# the inbuilt deviance residuals)
neg_loglik <- function(D, mu, fam, disp = 1) {
  dev <- fam$dev.resids(as.vector(D), 
                        as.vector(mu), 
                        wt = rep_len(1, length(D)))
  good <- !is.infinite(dev)
  return(sum(dev[good]) / (2 * disp))
}

# L0 penalized negative log-likelihood objective
penalised_nll <- function(W, H, D, fam, lambdaW, lambdaH) {
  disp <- estimate_phi(D, mu = W %*% H, family, # estimate dispersion if need be
               df = sum(colSums(W)!=0) + sum(rowSums(H)!=0))
  mu  <- W %*% H
  loss <- neg_loglik(D, mu, fam, disp)
  pen  <- lambdaW * sum(colSums(W) > 0) + lambdaH * sum(rowSums(H) > 0)
  return(loss + pen)
}

# helper function to estimate dispersion coefficient phi
estimate_phi <- function(D, mu, fam, df) {
  if (family$family %in% c("binomial","poisson")) { return(1) } else {
  dev <- fam$dev.resids(as.vector(D), as.vector(mu), wt = rep_len(1, length(D)))
  phi <- sum(dev) / (2 * (length(dev) - df))
  return(max(phi, 1e-12)) }             # avoid division by 0 later
}

# L0 or grouped L0 penalised generalised Best Rank NMF
# i.e. sparse NMF for any exponential family distribution (gaussian, poisson, gamma, inverse gamma)
# with automatic best rank recovery

# decomposes data matrix D as W %*% H subject to exponential family noise 
# with inbuilt best rank discovery 
# (we assume an identity link & impose nonnegativity constraints)
# function to be exported:
L0BestRankGNMF <- function(D, 
                     k,
                     Wstart = NULL,
                     grouped_Hupdate = FALSE, # use plain L0 or grouped L0 penalty in H update?
                     IC = "gic", # information criterion to be optimised
                     max_iter = 5, tol = 1e-4,
                     family = gaussian(link="identity")) {
  
  D <- Matrix(D, sparse = TRUE) # coerce to sparse format
  Dt <- t(D)  # keep for later
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
                                         parallelize = FALSE, # TODO : fix: with parallelize=T returns error "object 'intercept' not found"
                                         maxit = 10,
                                         group.maxit = 3,
                                         grouped = grouped_Hupdate, # apply L0 or grouped L0 penalty
                                         relaxed = !grouped_Hupdate, # refit selected groups of coefficients unpenalized at the end if grouped=T
                                         solver = "eigen", # use Cholesky solver
                                         method=1L)$coefficients) # k_init x p matrix
  } else {
     # NNDSVD initialisation
     require(NMF)
     suppressWarnings(nndsvd <- NMF:::.nndsvd.internal(A=D, k=k))
     W <- Matrix(nndsvd$W, sparse=T) # n x k
     H <- Matrix(nndsvd$H, sparse=T) # k x p
  } 
  
  ## --- drop factors that were fully penalized to zero ------------------------
  keep = which(rowSums(H)!=0 & colSums(W)!=0) # only keep factors with nonzero row sums in H & nonzero column sums in W
  H_old <- H <- H[keep, , drop=F] 
  W_old <- W <- W[, keep, drop=F] 
  k <- nrow(H)
  
  ## --- calculate objective ---------------------------------------------------
  obj <- penalised_nll(W, H, D, family, lambdaW, lambdaH)
  obj_vect <- obj
  
  for (iter in 1:max_iter) {
    ## --- update W -------------------------------------------------------
    suppressWarnings(W <- t(L0glm.fit(
                     X = t(H), # p x k matrix
                     y = Dt, # p x n data matrix, now regressing on rows of D
                     intercept = FALSE,
                     family = family,
                     lambda = lambdaW,
                     nonnegative = TRUE,
                     parallelize = FALSE, # TODO : fix: with parallelize=T returns error "object 'intercept' not found"
                     maxit = 10,
                     grouped = FALSE,
                     solver = "eigen", # use Cholesky solver
                     method = 1L)$coefficients)) # n x k matrix
    if (sum(range(W) == 0)==2) { 
      message("All factors were penalized to zero in the NMF")
      return(list(W = NULL, H = NULL, dev = obj, k = 0, iter = iter, obj = obj_vect)) }
    
    ## --- update H -------------------------------------------------------
    suppressWarnings(H <- L0glm.fit(
                   X = W, # n x k matrix
                   y = D, # n x p data matrix, now regressing on columns of D
                   intercept = FALSE,
                   family = family,
                   lambda = lambdaH,
                   nonnegative = TRUE,
                   parallelize = FALSE, # TODO : fix: with parallelize=T returns error "object 'intercept' not found"
                   maxit = 10,
                   group.maxit = 3,
                   grouped = grouped_Hupdate, # apply L0 or grouped L0 penalty
                   relaxed = !grouped_Hupdate, # refit selected groups of coefficients unpenalized at the end if grouped=T
                   solver = "eigen", # use Cholesky solver
                   method = 1L)$coefficients) # k x p matrix
    
    ## --- drop factors that were fully penalized to zero ------------------------
    keep <- which(rowSums(H)!=0 & colSums(W)!=0) # only keep factors with nonzero row sums in H & nonzero column sums in W
    H <- H[keep, , drop=F] 
    W <- W[, keep, drop=F] 
    k <- nrow(H)

    ## --- convergence check---------------------------------------------
    obj_new <- penalised_nll(W, H, D, family, lambdaW, lambdaH)
    if (obj_new > obj) { H <- H_old
                         W <- W_old
                         k <- nrow(H)
                         obj_new <- obj } # fit deteriorated so use previous iteration estimates
    obj_vect <- c(obj_vect, obj_new)
    if (abs(obj_new - obj)/(abs(obj) + 1e-12) < tol) break
    obj <- obj_new
  }
  
  ## --- rescale row sums of H to unit L2 norm --------------------------
  s <- sqrt(rowSums(H^2))
  H <- sweep(H, 1, s, FUN="/") # rescale rows of H
  W <- sweep(W, 2, s, FUN="*") # correspondingly rescale columns of W
  
  list(W = W, H = H, dev = obj, k = k, iter = iter, obj = obj_vect)
}



# EXAMPLE: deconvolute overlapping peaks in multichannel GC/MS data (blurred
# mutichannel spike train with Poisson noise) with n = 100 & true rank k = 10,
# i.e. with 10 compounds being present using L0BestRankGNMF
n <- 100 # timepoints
k_true <- 20 # true rank (nr of compounds)
psf_fun <- function(x, u, w=2.5) { exp(((x-u)^2)/(-2*(w^2))) } # gaussian point spread / peak shape function
psf <- psf_fun(x=1:n/5, u=n/2) # vector with PSF
sim.pois.mc <- simulate_spike_train(n = n, k = k_true, family = "poisson", 
                                    seed = 123,
                                    multichannel = TRUE,
                                    mean_beta = 10000,
                                    sparse = TRUE,
                                    sd_logbeta = 1,
                                    psf = psf_fun)
D = sim.pois.mc$y # n x p with n timepoints and p mass channels
image(D)
sim.pois.mc
IC = "gic" # information criterium to use 

# initialise elution profiles (W matrix)
# we initialise W with an n x n banded matrix with 
# time-shifted copies of the point spread / peak shape function (PSF) 
# (this can be estimated by running standards or tuned based on the IC of different L0glm fits 
# with banded matrices with different peak widths)
Wstart = sim.pois.mc$X  

# this matrix can subsequently be pre-subsetted by 
# regressing the row sums of D (=total ion chromatogram) on
# our banded covariate matrix Wstart using an unpenalised nnls fit, after which 
# we can retain only the columns of Wstart with nonzero coefficients 
# and use this as initialisation for W
peak_coefs <- L0glm.fit(X = Wstart, 
                        y = rowSums(D), # n x p data matrix, now regressing on row sums of D (=total ion chromatogram)
                        intercept = FALSE,
                        family = family,
                        lambda = 0, # use unregularised nnls
                        nonnegative = TRUE,
                        solver = "osqp" # use unregularised nnls via osqp solver
)$coefficients
keep = peak_coefs!=0
sum(keep) # 32
Wstart_rowsums = as.matrix(Wstart[, keep, drop=F])

system.time(nmf <- L0BestRankGNMF(D, 
                Wstart = Wstart_rowsums,  # initially selected elution profiles
                grouped_Hupdate = TRUE, # here in the H=t(S) with S=mass spectra we use a grouped L0 penalty 
                IC = "gic", # information criterion to be optimised
                max_iter = 5, tol = 1e-4,
                family = poisson(link="identity"))) # nonnegative identity link Poisson here
# 8s non-parallelized
dim(nmf$W)
dim(nmf$H)
nmf$iter # 1
nmf$k # 20 out of 20 components recovered
nmf$obj # objective over subsequent iterations

FP = (nmf$k - k_true) 
FP # 0 false positives
nmf$k / k_true # 1, i.e. exact best rank recovery

matplot(nmf$W, type="l") # obtained elution profiles
image(sim.pois.mc$beta_true[rowSums(sim.pois.mc$beta_true)!=0,, drop=F]) # true spectra
image(nmf$H) # recovered spectra
image(nmf$H-sim.pois.mc$beta_true[rowSums(sim.pois.mc$beta_true)!=0,, drop=F]) # residual against ground truth for spectra

image(as.matrix(nmf$W %*% nmf$H))    # sparse NMF reconstruction
image(as.matrix(sim.pois.mc$y))      # observed data matrix (with Poisson noise)
image(as.matrix(sim.pois.mc$y_true)) # ground truth (without Poisson noise)
image(as.matrix(nmf$W %*% nmf$H-sim.pois.mc$y_true)) # residual against ground truth
range(as.matrix(nmf$W %*% nmf$H-sim.pois.mc$y_true)) # -87.43010  91.38636
mean(abs(as.matrix(nmf$W %*% nmf$H-sim.pois.mc$y_true))) # MAE against ground truth = 1.98
image(as.matrix(nmf$W %*% nmf$H-sim.pois.mc$y)) # residual against observed signal
range(as.matrix(nmf$W %*% nmf$H-sim.pois.mc$y)) # -142.0393  181.0895

# TODO: calculate maximum pairwise cosine similarity between each recovered spectrum (rows in nmf$H)
# and the best matching true spectrum ()

# ALS::matchFactor(S[,1],test0$S[,1]) (dot product match factor)

# comparison with solution provided by ALS package (which uses plain alternating nonnegative least squares, i.e. ignoring Poisson noise and using a fixed rank decomposition, but including optional unimodality constraints)
library(ALS)
system.time(nmf_als <- als(CList = list(Wstart_rowsums), # starting values for elution profiles
               S = matrix(1, nrow=ncol(D), ncol=ncol(Wstart_rowsums)), # dummy matrix for S
               PsiList = list(as.matrix(D)), # data matrix 
               optS1st = TRUE,
               uniC = TRUE, 
               normS = 2)) # 1.6s
nmf_als$W = nmf_als$CList[[1]]
nmf_als$H = t(nmf_als$S)
nmf_als$k = nrow(nmf_als$H)
FP = (nmf_als$k - k_true) 
FP # 12 false positives
nmf_als$k / k_true # 1.6, i.e. estimated rank 60% too high

matplot(nmf_als$W, type="l") # obtained elution profiles
image(sim.pois.mc$beta_true[rowSums(sim.pois.mc$beta_true)!=0,, drop=F]) # true spectra
image(Matrix(nmf_als$H, sparse=T)) # recovered spectra

image(as.matrix(nmf_als$W %*% nmf_als$H))    # sparse NMF reconstruction
image(as.matrix(sim.pois.mc$y))      # observed data matrix (with Poisson noise)
image(as.matrix(sim.pois.mc$y_true)) # ground truth (without Poisson noise)
image(as.matrix(nmf_als$W %*% nmf_als$H-sim.pois.mc$y_true)) # residual against ground truth
range(as.matrix(nmf_als$W %*% nmf_als$H-sim.pois.mc$y_true)) # -114.8259  149.3363
mean(abs(as.matrix(nmf_als$W %*% nmf_als$H-sim.pois.mc$y_true))) # MAE against ground truth = 3.38
image(as.matrix(nmf_als$W %*% nmf_als$H-sim.pois.mc$y)) # residual against observed signal
range(as.matrix(nmf_als$W %*% nmf_als$H-sim.pois.mc$y)) # -82.60136  73.13906
