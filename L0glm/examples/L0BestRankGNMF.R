\donttest{
## ---------- GC/MS deconvolution toy example with Poisson noise ------------
if (requireNamespace("Matrix", quietly = TRUE)) {
  library(Matrix)
}
n  <- 200                 # time points
k0 <- 50                  # true rank/ number of compounds
set.seed(123)

peak_width <- 2.5         # unknown, but can be estimated from standards or tuned from chromatograms
psf_fun <- function(x,u,w = peak_width) exp(-0.5*((x-u)/w)^2)

sim <- simulate_spike_train(n       = n,
                            k       = k0,
                            family  = "poisson",
                            multichannel = TRUE,
                            mean_beta    = 1e4,
                            sparse       = TRUE,
                            psf          = psf_fun)

D <- sim$y                # observed GC/MS chromatogram

## build a design matrix with time-shifted copies of the PSF
Wband <- sim$X            # n x n banded matrix with time-shifted copies of PSF
## pre-select putative peaks using an unregularised nonnegative
## identity-link Poisson fit on the total ion chromatogram
sel   <- L0glm.fit(Wband, rowSums(D),
                   lambda = 0, nonnegative = TRUE,
                   family = poisson(link=identity),
                   solver = "osqp")$coefficients != 0
Winit <- as.matrix(Wband[, which(sel), drop = FALSE])

## sparse Poisson Best-Rank NMF with grouped-L0 penalty
nmf <- L0BestRankGNMF(D,
                      Wstart          = Matrix(Winit, sparse = TRUE),
                      grouped_Hupdate = TRUE,
                      IC              = "bic",
                      max_iter        = 5,
                      family          = poisson("identity"))

cat("True rank:", k0, "\n",
    "Rank recovered by Poisson L0BestRankGNMF:", nmf$k, "\n")
matplot(nmf$W, type="l", xlab="Time", ylab="Signal", 
main="Elution profiles inferred by Poisson L0BestRankGNMF")

## ---------- Comparison with MCR-ALS (alternating non-negative least squares) -
if (requireNamespace("ALS", quietly = TRUE)) {
  library(ALS)
  als_fit <- als(CList   = list(Winit),
                 S       = matrix(1, ncol(D), ncol(Winit)),
                 PsiList = list(as.matrix(D)),
                 uniC    = TRUE,  # impose unimodality constraints
                 maxiter = 8,
                 normS   = 1)

  cat(" Rank recovered by MCR-ALS:", nrow(t(als_fit$S)), "\n")
}
}
