devtools::load_all(".")

# timings for some torch & tensorflow solvers available via GPUmatrix

# see https://cran.r-project.org/web/packages/GPUmatrix/vignettes/vignette.html
library(remotes)

# install torch compiled binaries

options(timeout = 6000)
# install_torch() # see https://torch.mlverse.org/docs/articles/installation.html
kind <- "cu117" 
version <- "0.10.0"
options(repos = c(
  torch = sprintf("https://storage.googleapis.com/torch-lantern-builds/packages/%s/%s/", kind, version),
  CRAN = "https://cloud.r-project.org" # or any other from which you want to install the other R dependencies.
))
# install.packages("torch")
library(torch)

# remotes::install_github("ceslobfer/GPUmatrix", ref="main")
library(GPUmatrix)
# install_tensorflow(version="tf-nightly") # to install GPU version
# install_tensorflow(version="nightly-cpu") # to install CPU version
library(tensorflow)

set.seed(1)
sim = simulate_spike_train(n=10000, p=10000, 
                           k=1000, sparse=TRUE) # blurred spike train with 10000 timepoints & 200 blurred peaks
library(Matrix)
Xs = sim$X # sparse covariate matrix Xs of class dgCMatrix (time-shifted copies of the point spread function)
y = sim$y
weights <- 1/(y+0.1) # observation weights (approx 1/Poisson noise variance)
yw = y*sqrt(weights) # weighted y
ywaug = c(yw, rep(0, ncol(Xs))) # augmented weighted y
lam = 1 # 1/max(sim$beta_true) # 1
lambdas = rep(lam, ncol(Xs)) # ridge penalties on individual coefficients
Xsw = Diagonal(x=sqrt(weights)) %*% Xs # weighted X
Xswaug = rbind(Xsw, Diagonal(x=sqrt(lambdas))) # row augmented covariate matrix

# diagonal/Jacobi preconditioner initialisation, tailored to a weighted least squares ridge regression
diag_crossprod_Xswaug = colSums(Xswaug^2) # = diag(crossprod(Xswaug)))
P = Diagonal(x=1/diag_crossprod_Xswaug) # tailored diagonal preconditioner for weighted least square ridge regression
beta_init_jacobi = as.vector(P %*% crossprod(Xswaug, ywaug)) 
hist(as.vector(beta_init_jacobi)) 

library(microbenchmark)
microbenchmark(beta_init_jacobi <- as.vector(crossprod(Xswaug, ywaug) / colSums(Xswaug^2)), unit="s") # 0.001s, colSums(Xswaug^2) = diag(crossprod(Xswaug)))



# Cholesky decomposition based solver using base R functions
# (fast with IntelMKL, very slow without - better to use Eigen's solvers)
chol.solve = function (X, y, tol=.Machine$double.eps^(2/3)) {
  ch <- base::chol(crossprod(X)) # solve using Cholesky decomposition
  backsolve(ch, forwardsolve(ch, crossprod(X, y), upper = TRUE, trans = TRUE))
}

# least squares conjugate gradient solver, using pureR function from https://github.com/ceslobfer/GPUmatrix/issues/2
cg_solve2 = function(A, b, 
                      x=b[1:ncol(A)], # coefficient starting values
                      tol=1e-8) {
  # We will solve AtA x = Atb without explicitly computing AtA
  At = t(A)
  Atb = At %*% b
  r = Atb - At %*% (A %*% x)
  p = r
  rsold = sum(r^2)
  
  for (i in 1:length(Atb)) {
    Ap = At %*% (A %*% p)
    alpha = drop(rsold / crossprod(p,Ap))
    x = x + alpha * p
    r = r - alpha * Ap
    rsnew <- sum(r^2)
    if (sqrt(rsnew) < tol) break
    p = r + (rsnew / rsold) * p
    rsold = rsnew
  }
  return(x)
}


# timings Eigen sparse Cholesky solve, using CPU - VERY FAST
library(microbenchmark)
microbenchmark(s1 <- solve_sparse_chol(A = as(crossprod(Xswaug),"generalMatrix"),
                                        b = as.matrix(crossprod(Xswaug, ywaug)), 
                                        pivot = 1L,
                                        ordering = 0L)$coefficients, unit="s") # 0.04s using Eigen sparse Cholesky solver
cor(as.vector(s1),as.vector(sim$beta_true)) # 0.1547611
plot(as.vector(s1),as.vector(sim$beta_true))
hist(s1)


# timings Eigen dense Cholesky solve, using CPU - quite slow
system.time(s1_dense <- solve_dense_chol(A = as(crossprod(Xswaug), "matrix"),
                                          b = as.matrix(crossprod(Xswaug, ywaug)), 
                                          pivot = 1L)$coefficients) # 62s using Eigen dense Cholesky solver
max(abs(s1_dense-s1)) # 4.377387e-12



# timings on Intel Core i9-11950H / Nvidia RTX A2000

# # timings sparse Cholesky decomposition based solve of weighted least squares problem using base R
# system.time(s1 <- chol.solve(Xswaug, ywaug)) # 1.4s when using Intel MKL (MRO 4.0.2), 135s without Intel MKL
# cor(as.vector(s1),as.vector(sim$beta_true)) # 0.155
# plot(as.vector(s1),as.vector(sim$beta_true))
# 
# # timings dense Cholesky decomposition based solve of weighted least squares problem using base R
# system.time(s1_dense <- chol.solve(as(Xswaug, "matrix"), ywaug)) # 8.6s when using Intel MKL (MRO 4.0.2), 1017s without Intel MKL
# # cor(as.vector(s1),as.vector(sim$beta_true)) # 0.155
# # plot(as.vector(s1),as.vector(sim$beta_true))



# timings sparse least squares conjugate gradient pure R implementation, using CPU
system.time(s2 <- cg_solve2(Xswaug, ywaug)) # 0.7s
max(abs(s2-s1)) # 1.346478e-10
# cor(as.vector(s2),as.vector(sim$beta_true)) # 0.155
# plot(as.vector(s2),as.vector(sim$beta_true))

# torch & tensorflow solver timings
Xswaug_gpu <- gpu.matrix(Xswaug, dtype = "float32", device="cuda", sparse = T)
ywaug_gpu <- gpu.matrix(ywaug, dtype = "float32", device="cuda", sparse = F)  
Xswaug_cpu <- gpu.matrix(Xswaug, dtype = "float32", device="cpu", sparse = T)
ywaug_cpu <- gpu.matrix(ywaug, dtype = "float32", device="cpu", sparse = F)  
Xswaug_dense_gpu <- gpu.matrix(Xswaug, dtype = "float32", device="cuda", sparse = F)
Xswaug_dense_cpu <- gpu.matrix(Xswaug, dtype = "float32", device="cpu", sparse = F)
Xswaug_tensorflow_gpu <- gpu.matrix(Xswaug, type = "tensorflow", dtype = "float32", sparse = T)
ywaug_tensorflow_gpu <- gpu.matrix(ywaug, type = "tensorflow", dtype = "float32", sparse = F) 

# timings torch sparse Cholesky solve, using GPU
system.time(s3 <- chol_solve(crossprod(Xswaug_gpu, ywaug_gpu), chol(crossprod(Xswaug_gpu)))) # 0.71s, faster than Base R Cholesky solve with IntelMKL, but slower than Eigen's sparse Cholesky solver
max(abs(s3-s1)) # 8

# timings torch sparse Cholesky solve, using CPU
system.time(s4 <- chol_solve(crossprod(Xswaug_cpu, ywaug_cpu), chol(crossprod(Xswaug_cpu)))) # 3.06s, slightly faster than Base R Cholesky solve with IntelMKL, but slower than Eigen's sparse Cholesky solver
max(abs(s4-s1)) # 8

# timings torch dense Cholesky solve, using GPU - FASTEST DENSE SOLVER
system.time(s5 <- chol_solve(crossprod(Xswaug_dense_gpu, ywaug_gpu), 
                             chol(crossprod(Xswaug_dense_gpu)))) # 1.13s, faster than torch Cholesky solve using CPU, faster than Base R Cholesky solve with IntelMKL & faster than Eigen's dense Cholesky solver
max(abs(s5-s1)) # 8

# timings torch dense Cholesky solve, using CPU
system.time(s6 <- chol_solve(crossprod(Xswaug_dense_cpu, ywaug_cpu), 
                             chol(crossprod(Xswaug_dense_cpu)))) # 6.35s
max(abs(s6-s1)) # 8

# timings tensorflow sparse least squares conjugate gradient pure R implementation, using GPU
system.time(s7 <- cg_solve2(Xswaug_tensorflow_gpu, ywaug_tensorflow_gpu)) # 19s
max(abs(s7-s1)) # 1.073286e-05

# timings torch sparse least squares conjugate gradient pure R implementation, using GPU
system.time(s8 <- cg_solve2(Xswaug_gpu, ywaug_gpu)) # 1.68s
max(abs(s8-s1)) # 1.317928e-06

# timings torch sparse least squares conjugate gradient pure R implementation, using CPU
system.time(s9 <- cg_solve2(Xswaug_cpu, ywaug_cpu)) # 4.92s
max(abs(s9-s1)) # 6.344871e-05

# timings torch dense least squares conjugate gradient pure R implementation, using GPU
system.time(s10 <- cg_solve2(Xswaug_dense_gpu, ywaug_gpu)) # 116.53s
max(abs(s10-s1)) # 1.489769e-06

# timings torch dense least squares conjugate gradient pure R implementation, using CPU
system.time(s11 <- cg_solve2(Xswaug_dense_cpu, ywaug_cpu)) # 176.11s
max(abs(s11-s1)) # 2.542638e-05

# timings Eigen sparse QR solver, using CPU - slow
library(microbenchmark)
microbenchmark({s12 <- solve_sparse_QR(A = Xswaug,
                                       b = ywaug,
                                       ordering = 1L)$coefficients}, unit="s") # 510s using Eigen QR solver
max(abs(s12-s1)) # 1.17506e-11
# cor(as.vector(s12),as.vector(sim$beta_true)) # 0.47
# plot(as.vector(s12),as.vector(sim$beta_true))


# timings Eigen sparse conjugate gradient, using normal equations, using CPU
microbenchmark({A=as(crossprod(Xswaug),"generalMatrix");
                b=as.matrix(crossprod(Xswaug, ywaug));
                s15 <- solve_sparse_CG(A = A,
                       b = b,
                       x0 = beta_init_jacobi, # b[1:ncol(A)], # or even better: use beta_init below in section diagonal/jacobi preconditioner
                       tol = 1E-5,
                       maxit = 0L,
                       precond = 1L)$coefficients}, unit="s") # 0.03s using Eigen sparse conjugate gradient & using normal equations
max(abs(s15-s1)) # 0.001054265


# timings Eigen dense conjugate gradient, using normal equations, using CPU - FASTEST DENSE SOLVER ON CPU
microbenchmark({A=as(crossprod(Xswaug),"matrix");
                b=as.matrix(crossprod(Xswaug, ywaug));
s16 <- solve_dense_CG(A = A,
                       b = b,
                       x0 = beta_init_jacobi,
                       tol = 1E-5,
                       maxit = 0L,
                       precond = 1L)$coefficients}, unit="s") # 4.3s using Eigen dense conjugate gradient & using normal equations
max(abs(s16-s1)) # 0.001054265


# timings Eigen sparse least squares conjugate gradient, without forming normal equations & calculating XtX, using CPU - FASTEST
library(microbenchmark)
microbenchmark({ x0 <- as.vector(crossprod(Xswaug, ywaug) / colSums(Xswaug^2)); # colSums(Xswaug^2) = diag(crossprod(Xswaug)))
                 s17 <- solve_sparse_LSCG(A = Xswaug,
                                         b = ywaug,
                                         x0 = x0, # ywaug[1:ncol(Xswaug)], # or even better: use beta_init below in section diagonal/jacobi preconditioner
                                         tol = 1E-7,
                                         maxit = 0L,
                                         precond = 1L)$coefficients }, unit="s") # 0.04s using Eigen least squares conjugate gradient
max(abs(s17-s1)) # 4.748111e-06



# timings Eigen dense least squares conjugate gradient, without forming normal equations & calculating XtX, using CPU
system.time({x0 <- as.vector(crossprod(Xswaug, ywaug) / colSums(Xswaug^2));
             s18 <- solve_dense_LSCG(A = as(Xswaug, "matrix"),
                                     b = as.matrix(ywaug),
                                     x0 = as.matrix(x0), # as.matrix(ywaug[1:ncol(Xswaug)]),
                                     tol = 1E-7,
                                     maxit = 0L,
                                     precond = 1L)$coefficients}) # 18.9s using Eigen least squares conjugate gradient
max(abs(s18-s1)) # 6.044192e-06









