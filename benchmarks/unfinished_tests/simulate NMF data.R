# simulate some data to benchmark our sparse NMF algo
# right now only Gaussian with identity link supported, could do other identity link noise families too
simulate_nmf_data <- function(n = 500, p = 300, k = 10,
                              rho_W = 0.15,  # sparsity = proportion of non-zeros
                              rho_H = 0.10,
                              shape = 2,     # Gamma(shape, rate = shape) ⇒ mean 1
                              target_SNR = 5,
                              seed = 123) {
  
  set.seed(seed)
  
  ## 1.  sparse positive factors ------------------------------------------
  draw_sparse_gamma <- function(rows, cols, rho, shape) {
    mat <- matrix(0, rows, cols)
    keep <- matrix(runif(rows * cols) < rho, rows, cols)
    nnz  <- rgamma(sum(keep), shape = shape, rate = shape)  # mean = 1
    mat[keep] <- nnz
    mat
  }
  
  W <- draw_sparse_gamma(n, k, rho_W, shape)
  H <- draw_sparse_gamma(k, p, rho_H, shape)
  
  ## 2.  scale columns / rows to unit 2-norm ------------------------------
  if (any(colSums(W^2) == 0) || any(rowSums(H^2) == 0))
    stop("rank or sparsity level too high: some factors vanished completely")
  W <- sweep(W, 2, sqrt(colSums(W^2)),  "/")
  H <- sweep(H, 1, sqrt(rowSums(H^2)),  "/")
  
  ## 3.  noiseless mean matrix & noise variance ---------------------------
  mu  <- W %*% H
  sigma <- norm(mu, "F") / (target_SNR * sqrt(n * p))
  
  ## 4.  Gaussian noise ---------------------------------------------------
  E <- matrix(rnorm(n * p, sd = sigma), n, p)
  
  D <- mu + E
  
  D[D<0] <- 0
  
  list(D = D, W = W, H = H, sigma = sigma,
       SNR_empirical = norm(mu, "F")/norm(E, "F"))
}


# sim = simulate_nmf_data(n = 500, p = 300, k = 10, target_SNR = 50)
# image(sim$D)

# EXAMPLE: factorise gaussian NMF dataset with simulated example where correct rank k = 10
# using L0BestRankGNMF - this didn't work yet for some reason...
sim <- simulate_nmf_data(n = 200, p = 100, k = 10)
D <- sim$D
image(D)
IC = "gic" # information criterium to use 
nmf <- L0BestRankGNMF(D, 
                      k = 15,  # choose rank > true rank & let L0 penalty eliminate some of the factors to zero
                      grouped_Hupdate = FALSE, 
                      IC = IC, # information criterion to be optimised
                      max_iter = 5, tol = 1e-4,
                      family = gaussian(link="identity"))
# returns error Error in if (sum(scaling) != length(scaling)) { : 
# missing value where TRUE/FALSE needed - check source of error?
dim(nmf$W)
dim(nmf$H)
