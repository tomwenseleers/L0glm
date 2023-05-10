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
# install_tensorflow(version="nightly-cpu") # to install CPU version, installs OK
# install_tensorflow(version = "nightly-gpu") # to install GPU version - this one fails with an error
library(tensorflow)

sim = simulate_spike_train(n=10000, p=10000, k=200, sparse=TRUE) # blurred spike train with 10000 timepoints & 200 blurred peaks
library(Matrix)
Xs = sim$X # sparse covariate matrix Xs of class dgCMatrix (time-shifted copies of the point spread function)
y = sim$y
weights <- 1/(y+0.1) # observation weights (approx 1/Poisson noise variance)
yw = y*sqrt(weights) # weighted y
ywaug = c(yw, rep(0, ncol(Xs))) # augmented weighted y
lam = 1 # 1/max(sim$beta_true) # 1
lambdas = rep(lam, ncol(Xs)) # ridge penalties on individual coefficients
Xsw = Xs*Diagonal(x=sqrt(weights)) # weighted X
Xswaug = rbind(Xsw, Diagonal(x=sqrt(lambdas))) # row augmented covariate matrix

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


# timings on Intel Core i9-11950H / Nvidia RTX A2000

# timings sparse Cholesky decomposition based solve of weighted least squares problem using base R
system.time(s1 <- chol.solve(Xswaug, ywaug)) # 1.4s when using Intel MKL (MRO 4.0.2), 135s without Intel MKL
max(s1) # 0.9999187
cor(as.vector(s1),as.vector(sim$beta_true)) # 0.47
plot(as.vector(s1),as.vector(sim$beta_true))

# timings dense Cholesky decomposition based solve of weighted least squares problem using base R
system.time(s1_dense <- chol.solve(as(Xswaug, "matrix"), ywaug)) # 8.6s when using Intel MKL (MRO 4.0.2), 1017s without Intel MKL
max(abs(s1_dense-s1)) # 0
# cor(as.vector(s1),as.vector(sim$beta_true)) # 0.47
# plot(as.vector(s1),as.vector(sim$beta_true))

# timings sparse least squares conjugate gradient pure R implementation, using CPU
system.time(s2 <- cg_solve2(Xswaug, ywaug)) # 0.03s
max(abs(s2-s1)) # 1.658031e-10
# cor(as.vector(s2),as.vector(sim$beta_true)) # 0.47
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
system.time(s3 <- chol_solve(crossprod(Xswaug_gpu, ywaug_gpu), chol(crossprod(Xswaug_gpu)))) # 0.58s, faster than Base R Cholesky solve with IntelMKL, but slower than Eigen's sparse Cholesky solver
max(abs(s3-s1)) # 3.131536e-07

# timings torch sparse Cholesky solve, using CPU
system.time(s4 <- chol_solve(crossprod(Xswaug_cpu, ywaug_cpu), chol(crossprod(Xswaug_cpu)))) # 1.09s, slightly faster than Base R Cholesky solve with IntelMKL, but slower than Eigen's sparse Cholesky solver
max(abs(s4-s1)) # 2.365705e-07

# timings torch dense Cholesky solve, using GPU
system.time(s5 <- chol_solve(crossprod(Xswaug_dense_gpu, ywaug_gpu), chol(crossprod(Xswaug_dense_gpu)))) # 1.17s, faster than torch Cholesky solve using CPU, faster than Base R Cholesky solve with IntelMKL & faster than Eigen's dense Cholesky solver
max(abs(s5-s1)) # 3.131536e-07

# timings torch dense Cholesky solve, using CPU
system.time(s6 <- chol_solve(crossprod(Xswaug_dense_cpu, ywaug_cpu), chol(crossprod(Xswaug_dense_cpu)))) # 6.08s
max(abs(s6-s1)) # 2.365705e-07

# timings tensorflow sparse least squares conjugate gradient pure R implementation, using GPU
system.time(s7 <- cg_solve2(Xswaug_tensorflow_gpu, ywaug_tensorflow_gpu)) # 0.61s
max(abs(s7-s1)) # 1.073286e-05

# timings torch sparse least squares conjugate gradient pure R implementation, using GPU
system.time(s8 <- cg_solve2(Xswaug_gpu, ywaug_gpu)) # 0.06s
max(abs(s_torch_gpu-s1)) # 1.011147e-07

# timings torch sparse least squares conjugate gradient pure R implementation, using CPU
system.time(s9 <- cg_solve2(Xswaug_cpu, ywaug_cpu)) # 0.03s
max(abs(s9-s1)) # 1.011147e-07

# timings torch dense least squares conjugate gradient pure R implementation, using GPU
system.time(s10 <- cg_solve2(Xswaug_dense_gpu, ywaug_gpu)) # 4.07s
max(abs(s10-s1)) # 1.011147e-07

# timings torch dense least squares conjugate gradient pure R implementation, using CPU
system.time(s11 <- cg_solve2(Xswaug_dense_cpu, ywaug_cpu)) # 5.20s
max(abs(s11-s1)) # 1.011147e-07

# timings Eigen sparse QR solver, using CPU
library(microbenchmark)
microbenchmark({s12 <- solve_sparse_QR(A = Xswaug,
                                       b = ywaug,
                                       ordering = 1L)$coefficients}, unit="s") # 0.0263564s using Eigen QR solver
max(abs(s12-s1)) # 2.375877e-14
# cor(as.vector(s12),as.vector(sim$beta_true)) # 0.47
# plot(as.vector(s12),as.vector(sim$beta_true))

# timings Eigen sparse Cholesky solve, using CPU
microbenchmark(s13 <- solve_sparse_chol(A = as(crossprod(Xswaug),"dgCMatrix"),
                                        b = as.matrix(crossprod(Xswaug, ywaug)), 
                                        pivot = 1L,
                                        ordering = 0L)$coefficients, unit="s") # 0.00093359s using Eigen sparse Cholesky solver, 135/0.00093359=144603x faster than native R Cholesky when not using IntelMKL
max(abs(s13-s1)) # 3.330669e-16
cor(as.vector(s13),as.vector(sim$beta_true)) # 0.47
plot(as.vector(s13),as.vector(sim$beta_true))
hist(s13)
max(s13) # maximum value of obtained coefficients is 1/lam


# timings Eigen dense Cholesky solve, using CPU
system.time(s14_dense <- solve_dense_chol(A = as(crossprod(Xswaug), "matrix"),
                                          b = as.matrix(crossprod(Xswaug, ywaug)), 
                                          pivot = 1L)$coefficients) # 57s using Eigen dense Cholesky solver
max(abs(s14_dense-s1)) # 3.330669e-16

# timings Eigen sparse conjugate gradient, using normal equations, using CPU
microbenchmark({A=as(crossprod(Xswaug),"dgCMatrix");
                b=as.matrix(crossprod(Xswaug, ywaug));
                s15 <- solve_sparse_CG(A = A,
                       b = b,
                       x0 = b[1:ncol(A)], # or even better: use beta_init below in section diagonal/jacobi preconditioner
                       tol = 1E-5,
                       maxit = 0L,
                       precond = 1L)$coefficients}, unit="s") # 0.000690752s using Eigen sparse conjugate gradient & using normal equations
max(abs(s15-s1)) # 4.440892e-16


# timings Eigen sparse least squares conjugate gradient, using normal equations, using CPU
microbenchmark({A=as(crossprod(Xswaug),"dgCMatrix");
                b=as.matrix(crossprod(Xswaug, ywaug));
                s16 <- solve_sparse_LSCG(A = A,
                                         b = b,
                                         x0 = b[1:ncol(A)], # or even better: use beta_init below in section diagonal/jacobi preconditioner
                                         tol = 1E-5,
                                         maxit = 0L,
                                         precond = 1L)$coefficients}, unit="s") # 0.000688393s using Eigen sparse least squares conjugate gradient & using normal equations
max(abs(s16-s1)) # 4.440892e-16

# timings Eigen sparse least squares conjugate gradient, without forming normal equations & calculating XtX, using CPU - FASTEST
microbenchmark({s17 <- solve_sparse_LSCG(A = Xswaug,
                                         b = ywaug,
                                         x0 = ywaug[1:ncol(Xswaug)], # or even better: use beta_init below in section diagonal/jacobi preconditioner
                                         tol = 1E-5,
                                         maxit = 0L,
                                         precond = 1L)$coefficients}, unit="s") # 0.00016174s using Eigen least squares conjugate gradient, without forming normal equations, 135/0.00016174=834672.9x faster than native R Cholesky solve when not using IntelMKL
max(abs(s17-s1)) # 1.957323e-13

# timings Eigen dense least squares conjugate gradient, without forming normal equations & calculating XtX, using CPU
system.time({s18 <- solve_dense_LSCG(A = as(Xswaug, "matrix"),
                                     b = as.matrix(ywaug),
                                     x0 = as.matrix(ywaug[1:ncol(Xswaug)]),
                                     tol = 1E-5,
                                     maxit = 0L,
                                     precond = 1L)$coefficients}) # 1.42s using Eigen least squares conjugate gradient
max(abs(s18-s1)) # 1.957323e-13


# diagonal/Jacobi preconditioner initialisation, tailored to a weighted least squares ridge regression - BEST
diag_crossprod_Xswaug = colSums(Xswaug^2) # = diag(crossprod(Xswaug)))
P = Diagonal(x=1/diag_crossprod_Xswaug) # tailored diagonal preconditioner for weighted least square ridge regression
beta_hat = P %*% crossprod(Xswaug, ywaug) # these coefficients range between 0 and 1/lambda 
hist(as.vector(beta_hat)) 

library(microbenchmark)
microbenchmark(beta_hat <- crossprod(Xswaug, ywaug) / colSums(Xswaug^2), unit="s") # 0.0002s, colSums(Xswaug^2) = diag(crossprod(Xswaug)))

# obtained coefficients are nearly identical to the actual ridge regression solution :
range(beta_hat) # 0.0000000 0.9999187 - coefficients are bounded between 0 and 1/lambda
max(abs(beta_hat-s1)) # 3.330669e-16, where s1 is my actual ridge regression Cholesky solver solution
cor(as.vector(beta_hat), s1) # 1



deco_nnL0poisglm = function(X = sim$X, y = sim$y, 
                            lambda=1, 
                            delta=1E-10, 
                            gamma=2,
                            thresh=1E-10) {  
  print(sum(sim$beta_true!=0))
  
  weights = 1/(as.vector(y)+0.1) # observation weights (approx 1/Poisson noise variance)
  yw = y*sqrt(weights) # weighted y
  n = nrow(X)
  p = ncol(X)
  ywaug = c(yw, rep(0, p)) # augmented weighted y
  penweights = rep(1, p) # ridge regression penalty weights on individual coefficients (sum to p)
  lambdas = lambda*penweights # ridge penalties on individual coefficients
  # Xw = X * sqrt(weights) # weighted X 
  Xw = Diagonal(x=sqrt(weights)) %*% X
  Xwaug = rbind(Xw, Diagonal(x=sqrt(lambdas))) # row augmented covariate matrix
  # X.n = colSums(Xwaug^2) # squared L2 norm normalisation of columns of Xwaug
  # Xwaug_scaled = Xwaug / X.n
  # X.n = Diagonal(x=1/colSums(Xwaug^2)) # Jacobi/diagonal type preconditioning matrix
  # beta_hat = crossprod(Xwaug, ywaug) / colSums(Xwaug^2) # coefficient initialisation based on preconditioning matrix
  beta_hat = crossprod(Xwaug, ywaug) / colSums(Xwaug^2) # coefficient initialisation based on preconditioning matrix
  # beta_hat = lmridge(X=X, # /X.n,
  #         y=y,
  #         weights=weights,
  #         lambdas=lambdas,
  #         nonnegative=TRUE,
  #         solver="osqp") 
  hist(as.vector(beta_hat))
  # beta_hat = beta_hat_scaled# * X.n
  beta_hat[beta_hat<thresh] = 0
  # hist(as.vector(beta_hat)) # bounded between 0 and 1/lambda
  hist(lambdas)
  
  for (i in 1:100) {  
    nonzero = as.vector(beta_hat!=0) # nonzero coefficients
    print(sum(nonzero))
    penweights[nonzero] = 1/(1+abs(beta_hat[nonzero])^gamma+delta^gamma) # update L2 norm ridge regression penalty weights
    lambdas[nonzero] = lambda*penweights[nonzero] # updated ridge regression penalties lambdas
    hist(lambdas)
    # Xwaug = rbind(Xw[,nonzero,drop=F], Diagonal(x=sqrt(lambdas[nonzero])))  
    # y waug = c(yw, rep(0, sum(nonzero)))
    # P = Diagonal(x=1/colSums(Xwaug^2)) # Jacobi/diagonal type preconditioning matrix
    # beta_hat[nonzero] = lmridge(X=Xw[,nonzero,drop=F],
    #                             y=yw,
    #                             # weights=weights[nonzero],
    #                             lambdas=lambdas[nonzero],
    #                             nonnegative=TRUE,
    #                             solver="osqp") # * P # P %*% crossprod(Xwaug, ywaug) # updated coefficients
    Xwaug = rbind(Xw[,nonzero,drop=F], Diagonal(x=sqrt(lambdas[nonzero])))
    ywaug = c(yw, rep(0, sum(nonzero)))
    x0 = as.vector(crossprod(Xwaug, ywaug) / colSums(Xwaug^2))
    # beta_hat[nonzero] = lmridge(X=Xw[,nonzero,drop=F],
    #                             y=yw,
    #                             # weights=weights[nonzero],
    #                             lambdas=lambdas[nonzero],
    #                             nonnegative=TRUE,
    #                             solver="eigen",
    #                             method=4)
    X.n = sqrt(colSums(Xw[,nonzero,drop=F] ^ 2))
    beta_hat[nonzero] = solve_sparse_LSCG(A = Xw[,nonzero,drop=F] / X.n,
                      b = matrix(yw),
                      x0 = matrix(x0), # or even better: use beta_init below in section diagonal/jacobi preconditioner
                      tol = 1E-3,
                      maxit = 0L,
                      precond = 1L)$coefficients * X.n
    beta_hat[beta_hat<thresh] = 0
    hist(as.vector(beta_hat))
    plot(sim$beta_true, as.vector(beta_hat))  
    cor(as.vector(sim$beta_true), as.vector(beta_hat)) # 0.47  
  }

  sum(sim$beta_true!=0)
  sum(beta_hat!=0)
 

}