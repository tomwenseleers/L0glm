# Install instructions OpenBLAS for Windows
# see https://github.com/david-cortes/R-openblas-in-windows?tab=readme-ov-file

# timings with OpenBLAS installed and using 16 threads 
# (after setting system environment variable OPENBLAS_NUM_THREADS to 16)

library(microbenchmark)
set.seed(1)
m <- 1e5
n <- 1e3
k <- 3e2
X <- matrix(rnorm(m*k), nrow=m)
Y <- matrix(rnorm(n*k), ncol=n)
microbenchmark({
  Z <- X %*% Y
}, times=10L) # mean 384 ms (appears multithreaded)

