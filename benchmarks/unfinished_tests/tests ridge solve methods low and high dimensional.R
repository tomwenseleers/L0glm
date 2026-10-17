# fix to set working package repository when using Microsoft R Open with Intel MKL
options (repos = c(CRAN="https://cloud.r-project.org/"))
# install.packages("Rcpp")
library(Matrix)
library(nnls)
library(microbenchmark)
library(sanic)
library(osqp)
library(devtools)
library(glmnet)
devtools::load_all()
library(RhpcBLASctl)

blas_threads = get_num_cores() # 8
omp_threads = get_num_procs() # 16

blas_get_num_procs() # 8
blas_set_num_threads(blas_threads)
omp_get_num_procs() # 16
omp_get_max_threads()
omp_set_num_threads(omp_threads)


# couple of solves (also: QR decomposition based solve = solve.qr(a=X, b=y))
# SVD based solve (good when p > n & lambda=0) (SLOW)
library(pracma)
pinv.solve = function (X, y, tol=.Machine$double.eps^(2/3)) pracma::pinv(X, tol=tol) %*% y # Moore-Penrose / SVD-based solve for p>n case

# Cholesky decomposition based solve (FAST)
chol.solve = function (X, y, tol=.Machine$double.eps^(2/3)) {
  ch <- base::chol(crossprod(X)) # solve using Cholesky decomposition
  backsolve(ch, forwardsolve(ch, crossprod(X, y), upper = TRUE, trans = TRUE))
}

# Simulate data
set.seed(1)
n <- 500 # 10000 # 10000 # 500
p <- 50000 # 10000 # 20000 # 50000
s <- 0.1 # sparsity as proportion of p that have nonzero coefficients // p*s == #of peaks
k <- round(p*s) # nr of nonzero covariates
sim <- simulate_spike_train(n=n, p=p, k=k, mean_beta = 1000, sd_logbeta = 1, family="poisson")

system.time(Matrix::crossprod(sim$X)) # 0.02s for n=p=10000, 0.1s for n=1000 & p=10000

coefs_true <- as.vector(sim$beta_true)

Xsp <- sim$X # sparse covariate matrix X
y <- as.numeric(as.vector(sim$y))
weights <- 1/(y+0.1) # observation weights
Xwsp <-  Xweights(X=Xsp, weights=weights) # weighted X=Diagonal(x=sqrt(weights)) %*% Xsp
yw <- Yweights(y, weights) # weighted y = y*sqrt(weights)

lambdas <- rep(1, p) # = lambda*pen_weights
# NOTE: if you would like pen_weights to sum to length(lambdas) you would have to use
# lambdas <- lambdas / sum(lambdas) * length(lambdas)
lambda <- mean(lambdas+1E-100) # average ridge penalty
pen_weights <- (lambdas+1E-100) / lambda # penalization weights
D <- Diagonal(x=1/pen_weights) # cf notation Liu & Li 2016 # 1/pen_weights
prior.mean <- rep(0, p) # prior.mean = ridge center for nonzero-centered ridge regression

Xwaugsp <- row_augX(Xw=Xwsp, lambdas) # sparse augmented version
X <- as.matrix(Xsp) # dense version
Xw <- as.matrix(Xwsp) # dense version
Xwaug <- as.matrix(Xwaugsp) # dense version

ywaug <- augY(yw = yw, p = p, prior.mean = prior.mean, lambdas = lambdas)
lower <- rep(0, p)
upper <- rep(Inf, p)
nonnegative <- TRUE
start <- rep(1, p)
verbose <- TRUE


# # example of multichannel signal
# yvect <- c(Y)
# library(Matrix)
# # create the weight matrix
# W <- Diagonal(sqrt(1/(y+0.1)))
# # create the block diagonal matrix
# Xblock <- bdiag(lapply(1:nrow(Y), function(i) W %*% X))
#
# #
# library(spam)
# # create the weight matrix
# W <- Diagonal(sqrt(1/(yvect+0.1)))
# # create the block diagonal matrix
# Xblock <- bKMatrix(W, X)


# 1. TESTS WITH SOME DEFAULT SOLVE METHODS TO SOLVE ZERO OR NONZERO CENTERED RIDGE REGRESSION, EITHER WITH DIRECT SOLVE OR VIA ROW AUGMENTATION ####
# (using default BLAS that R is compiled against, can be Accelerate (Mac), OpenBLAS (Ubuntu/Win) or IntelMKL (Ubuntu/Win))

# 1.1 SOLVE VIA NORMAL EQUATIONS FOR n >= p, CF. Liu & Li 2016, here adapted to also cover nonzero centered ridge regression ####
# dense case
system.time(coefs_solve_largen_dense <- as.vector(Matrix::solve(
  as.matrix(D %*% Matrix::crossprod(Xwsp) + # FAST with Intel MKL
            Diagonal(x=rep(lambda, p))),
  as.matrix(D %*% Matrix::crossprod(Xwsp,yw) + lambda*Matrix::crossprod(D, prior.mean)),
            tol=2*.Machine$double.eps))) 
# 3.3s for n=p=10000 or 51s for n=10000 & p=20000 or 8s for n=1000 & p=10000 with MRO & Intel MKL, 1.42s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
# sparse case
system.time(coefs_solve_largen_sparse <- as.vector(Matrix::solve(
  D %*% Matrix::crossprod(Xwsp) + # FAST
  Diagonal(x=rep(lambda, p)),
  as.matrix(D %*% Matrix::crossprod(Xwsp,yw) + lambda*Matrix::crossprod(D, prior.mean))))) 
# 0.05s for n=p=10000 or XXs for n=10000 & p=20000 or XXs for n=1000 & p=10000 with Intel MKL, cf. Liu & Li 2016 eqn (10)

cor(coefs_solve_largen_sparse, coefs_true) # 0.1572276 for n=p=10000 / 0.1087041 for n=10000 & p=20000
cor(coefs_solve_largen_dense, coefs_solve_largen_sparse) # 1

# 1.2 SOLVE USING OPTIMISED FORMULA FOR p > n, CF. Liu & Li 2016 eqn (28), here adapted to also cover nonzero centered ridge regression ####
# dense case
system.time({
  # Precompute some values
  DXw <- Matrix::tcrossprod(D, Xwsp)
  y_star <- as.vector(yw - Xwsp %*% prior.mean)
  
  # Calculate β using Woodbury identity
  coefs_solve_largep_dense <- as.vector(prior.mean + DXw %*% Matrix::solve(
    as.matrix(Xwsp %*% DXw + Diagonal(x=rep(lambda, n))),
    y_star,
    tol=2*.Machine$double.eps ) )
}) # 3.1s for n=p=10000 or 3.3s for n=10000 and p=20000 or 0.06s for n=500 and p=50000
cor(coefs_solve_largen_dense, coefs_solve_largep_dense) # 1
cor(coefs_solve_largen_dense, coefs_true) # 0.1557878 / 0.1087041 for n=10000 and p=20000
cor(coefs_solve_largep_dense, coefs_true) # 0.1557878

# sparse case
system.time({
  # Precompute some values
  DXw <- Matrix::tcrossprod(D, Xwsp)
  y_star <- as.vector(yw - Xwsp %*% prior.mean)
  
  # Calculate β using Woodbury identity
  coefs_solve_largep_sparse <- as.vector(prior.mean + DXw %*% Matrix::solve(
    Xwsp %*% DXw + Diagonal(x=rep(lambda, n)),
    y_star,
    tol=2*.Machine$double.eps )) 
}) 
# 0.06s for n=p=10000 or 0.06s for n=500 and p=50000 or 
# XXs for n=10000 and p=20000 or 0.08s for n=500 and p=50000
cor(coefs_solve_largen_sparse, coefs_solve_largep_sparse) # 1
cor(coefs_solve_largen_sparse, coefs_true) # 0.1557878 / 0.1087041 for n=10000 and p=20000
cor(coefs_solve_largep_sparse, coefs_true) # 0.1557878


# 1.3 SOLVE VIA ROW AUGMENTATION ####

# sparse case 
# NOTE: if solved via row augmentation then always better to use sparse Xwaugsp
# for p > n case ridge & all solvers except LSCG best solved via Woodbury optimised solution method & normal equations
system.time(coefs_solve_rowaug_sparse_LSCG <- as.vector(solve_sparse_LSCG(A = Xwaugsp,
                                                                          b = as.matrix(ywaug),
                                                                          x0 = as.matrix(start))$coefficients)) 
# 3.99s for n=10000 & p=10000, 0.01s for n=500 & p=50000 - FAST FOR p > n

system.time(coefs_solve_rowaug_sparse_QR <- as.vector(solve_sparse_QR(A = Xwaugsp,
                                                                      b = as.matrix(ywaug))$coefficients)) 
# 2.62s for n=10000 & p=10000, SLOW for n=500 & p=50000

system.time(coefs_solve_rowaug_sparse_CHOL <- as.vector(solve_sparse_chol(A = as(Matrix::crossprod(Xwaugsp),"dgCMatrix"),
                                                                          b = as.matrix(Matrix::crossprod(Xwaugsp, ywaug)))$coefficients))
# 0.03s for n=10000 & p=10000, XXs for n=500 & p=50000

system.time(coefs_solve_rowaug_sparse_CG <- as.vector(solve_sparse_CG(A = as(Matrix::crossprod(Xwaugsp),"dgCMatrix"),
                                                                      b = as.matrix(Matrix::crossprod(Xwaugsp, ywaug)),
                                                                      x0 = as.matrix(start))$coefficients))
# 0.05s for n=10000 & p=10000, 9.7s for n=500 & p=50000


# TESTS OF IMPLEMENTED lmridge SOLVERS ####
# sparse case
system.time(coefs_eigen1 <- lmridge(X=Xsp, y=y, 
                                    lambdas = rep(1, p),
                                    weights = weights,
                                    lower = rep(0, p),
                                    upper = rep(.Machine$double.xmax,ncol(X)),
                                    nonnegative = TRUE,
                                    prior.mean = NULL,
                                    solver = "eigen",
                                    method = 5L,
                                    verbose = T,
                                    start = rep(1, ncol(X)),
                                    control.solver = do.call("control.solver.gen", 
                                                             list(solver = "eigen", method = 5L))))
# 2.63s for n=p=10000, 0.08s for n=500 and p=50000

system.time(coefs_eigen3 <- lmridge(X=Xsp, y=y, 
                    lambdas = rep(1, p),
                    weights = weights,
                    lower = rep(0, p),
                    upper = rep(.Machine$double.xmax,ncol(X)),
                    nonnegative = TRUE,
                    prior.mean = NULL,
                    solver = "eigen",
                    method = 4L,
                    start = rep(1, ncol(X)),
                    verbose = T,
                    control.solver = modifyList(do.call("control.solver.gen", 
                                                        list(solver = "eigen", 
                                                             method = 4L)),
                                                list(rowaug=TRUE)) ))  
# this one sets rowaug=TRUE
# 0.04s for n=p=10000, 0.06s for n=500 and p=50000
# TODO: check - according to row augmentation benchmark above this should only
# take 0.02s instead of 0.06s with this solver

system.time(coefs_eigen4 <- lmridge(X=Xsp, y=y, 
                                    lambdas = rep(1, p),
                                    weights = weights,
                                    lower = rep(0, p),
                                    upper = rep(.Machine$double.xmax,ncol(X)),
                                    nonnegative = TRUE,
                                    prior.mean = NULL,
                                    solver = "eigen",
                                    method = 4L,
                                    start = rep(1, ncol(X)),
                                    verbose = T,
                                    control.solver = modifyList(do.call("control.solver.gen", 
                                                                        list(solver = "eigen", 
                                                                             method = 4L)), 
                                                                list(rowaug=TRUE)) )) # this one sets rowaug=TRUE
# 0.04s for n=p=10000, 0.07s for n=500 and p=50000 - FASTEST FOR p > n for both DENSE & SPARSE
# TODO: check - according to row augmentation benchmark above this should only
# take 0.02s instead of 0.07s with this solver

system.time(coefs_eigen5 <- lmridge(X=Xsp, y=y, 
                                    lambdas = rep(1, p),
                                    weights = weights,
                                    lower = rep(0, p),
                                    upper = rep(.Machine$double.xmax,ncol(X)),
                                    nonnegative = TRUE,
                                    prior.mean = NULL,
                                    solver = "eigen",
                                    method = 3L,
                                    verbose = T,
                                    start = rep(1, ncol(X)),
                                    control.solver = modifyList(do.call("control.solver.gen", 
                                                                        list(solver = "eigen", 
                                                                             method = 3L)),
                                                                list(rowaug=FALSE)) )) 
# 0.03s for n=p=10000, 0.06s for n=500 and p=50000

system.time(coefs_eigen7 <- lmridge(X=Xsp, y=y, 
                                    lambdas = rep(1, p),
                                    weights = weights,
                                    lower = rep(0, p),
                                    upper = rep(.Machine$double.xmax,ncol(X)),
                                    nonnegative = TRUE,
                                    prior.mean = NULL,
                                    solver = "eigen",
                                    verbose =T,
                                    method = 1L,
                                    start = rep(1, ncol(X)),
                                    control.solver = modifyList(do.call("control.solver.gen", 
                                                                        list(solver = "eigen", 
                                                                             method = 1L)),
                                                                list(rowaug=TRUE)) )) 
# 0.05s for n=p=10000, 0.07s for n=500 and p=50000

system.time(coefs_osqp <- lmridge(X=Xsp, y=y, 
                                   # Xw=Xw, yw=yw,
                                   lambdas = rep(1, p),
                                   weights = weights,
                                   lower = rep(0, p),
                                   upper = rep(.Machine$double.xmax,ncol(X)),
                                   start = rep(1, ncol(X)),
                                   nonnegative = TRUE,
                                   prior.mean = NULL,
                                   solver = "osqp", #,
                                   control.solver = do.call("control.solver.gen", 
                                                            list(solver = "osqp")))) 
# 0.15s for n=p=10000, crashes for n=500 and p=50000

system.time(coefs_glmnet <- lmridge(X=Xsp, y=y, 
                                    # Xw=Xw, yw=yw,
                                    lambdas = rep(1, p),
                                    weights = weights,
                                    lower = rep(0, p),
                                    upper = rep(.Machine$double.xmax,ncol(X)),
                                    start = rep(1, ncol(X)),
                                    nonnegative = TRUE,
                                    prior.mean = NULL,
                                    solver = "glmnet", #,
                                    control.solver = do.call("control.solver.gen", 
                                                             list(solver = "glmnet")))) 
# 0.03s for n=p=10000, 0.1s for n=500 and p=50000


# dense case
system.time(coefs_eigen1 <- lmridge(X=X, y=y, 
                                    lambdas = rep(1, p),
                                    weights = weights,
                                    lower = rep(0, p),
                                    upper = rep(.Machine$double.xmax,ncol(X)),
                                    nonnegative = TRUE,
                                    prior.mean = NULL,
                                    solver = "eigen",
                                    verbose = T,
                                    method = 5L,
                                    start = rep(1, ncol(X)),
                                    control.solver = modifyList(do.call("control.solver.gen", 
                                                                        list(solver = "eigen", 
                                                                             method = 5L)),
                                                                list(rowaug=FALSE)) )) 
# 6.1s for n=p=10000, 0.32s for n=500 and p=50000

system.time(coefs_eigen3 <- lmridge(X=X, y=y, 
                                    lambdas = rep(1, p),
                                    weights = weights,
                                    lower = rep(0, p),
                                    upper = rep(.Machine$double.xmax,ncol(X)),
                                    nonnegative = TRUE,
                                    prior.mean = NULL,
                                    solver = "eigen",
                                    verbose = T,
                                    method = 4L,
                                    start = rep(1, ncol(X)),
                                    control.solver = modifyList(do.call("control.solver.gen", 
                                                                        list(solver = "eigen", 
                                                                             method = 4L)),
                                                                list(rowaug=FALSE)) )) # this one sets rowaug=TRUE
# 3.6s for n=p=10000, 0.2s for n=500 and p=50000

system.time(coefs_eigen4 <- lmridge(X=X, y=y, 
                                    lambdas = rep(1, p),
                                    weights = weights,
                                    lower = rep(0, p),
                                    upper = rep(.Machine$double.xmax,ncol(X)),
                                    nonnegative = TRUE,
                                    prior.mean = NULL,
                                    verbose = T,
                                    solver = "eigen",
                                    method = 4L,
                                    start = rep(1, ncol(X)),
                                    control.solver = modifyList(do.call("control.solver.gen", 
                                                                        list(solver = "eigen", 
                                                                             method = 4L)), 
                                                                list(rowaug=TRUE)) )) # this one sets rowaug=TRUE
# 3.6s for n=p=10000, 0.2s for n=500 and p=50000


system.time(coefs_eigen5 <- lmridge(X=X, y=y, 
                                    lambdas = rep(1, p),
                                    weights = weights,
                                    lower = rep(0, p),
                                    upper = rep(.Machine$double.xmax,ncol(X)),
                                    nonnegative = TRUE,
                                    prior.mean = NULL,
                                    solver = "eigen",
                                    verbose = T,
                                    method = 3L,
                                    start = rep(1, ncol(X)),
                                    control.solver = modifyList(do.call("control.solver.gen", 
                                                                        list(solver = "eigen", 
                                                                             method = 3L)),
                                                                list(rowaug=FALSE)) )) 
# 3.5s for n=p=10000, 0.27s for n=500 and p=50000

system.time(coefs_eigen7 <- lmridge(X=X, y=y, 
                                    lambdas = rep(1, p),
                                    weights = weights,
                                    lower = rep(0, p),
                                    upper = rep(.Machine$double.xmax,ncol(X)),
                                    nonnegative = TRUE,
                                    prior.mean = NULL,
                                    verbose = T,
                                    solver = "eigen",
                                    method = 1L,
                                    start = rep(1, ncol(X)),
                                    control.solver = modifyList(do.call("control.solver.gen", 
                                                                        list(solver = "eigen", 
                                                                             method = 1L)),
                                                                list(rowaug=FALSE)) )) 
# XXs for n=p=10000, 0.24s for n=500 and p=50000

system.time(coefs_osqp <- lmridge(X=X, y=y, 
                                  # Xw=Xw, yw=yw,
                                  lambdas = rep(1, p),
                                  weights = weights,
                                  lower = rep(0, p),
                                  upper = rep(.Machine$double.xmax,ncol(X)),
                                  start = rep(1, ncol(X)),
                                  nonnegative = TRUE,
                                  prior.mean = NULL,
                                  solver = "osqp", #,
                                  control.solver = do.call("control.solver.gen", 
                                                           list(solver = "osqp")) )) 
# 3.7s for n=p=10000, crashes for n=500 and p=50000

system.time(coefs_glmnet <- lmridge(X=X, y=y, 
                                    # Xw=Xw, yw=yw,
                                    lambdas = rep(1, p),
                                    weights = weights,
                                    lower = rep(0, p),
                                    upper = rep(.Machine$double.xmax,ncol(X)),
                                    start = rep(1, ncol(X)),
                                    nonnegative = TRUE,
                                    prior.mean = NULL,
                                    verbose = T,
                                    solver = "glmnet", #,
                                    control.solver = do.call("control.solver.gen", 
                                                             list(solver = "glmnet")))) 
# 3.6s for n=p=10000, 0.25s for n=500 and p=50000







# OLD CODE BELOW ####

system.time(coefs_qrsolve_dense <- as.vector(qr.solve(D %*% crossprod(Xw) +
                                                   diag(rep(lambda, p)),
                                                 D %*% crossprod(Xw,yw),
                                                 tol=2*.Machine$double.eps))) # 386.13s for n=p=10000, 27.47s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
system.time(coefs_qrsolve_sparse <- as.vector(qr.solve(D %*% crossprod(Xwsp) +
                                                    .sparseDiagonal(p, lambda),
                                                  D %*% crossprod(Xwsp,yw),
                                                  tol=2*.Machine$double.eps))) # 298s for n=p=3000, 28.35s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
cor(coefs_solve_sparse, coefs_true) # 0.1651061
cor(coefs_solve_dense, coefs_qrsolve_dense) # 1
cor(coefs_solve_dense, coefs_qrsolve_sparse) # 1

# eigen's QR methods: slower than default QR solve methods
system.time(coefs_eigenqrsolvepivot0_dense <- as.vector(solve_dense_QR(D %*% crossprod(Xw) +
                                                        diag(rep(lambda, p)),
                                                     D %*% crossprod(Xw,yw),
                                                      pivot = 0L)$coefficients)) # 40.62s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
system.time(coefs_eigenqrsolveord0_sparse <- as.vector(solve_sparse_QR(D %*% crossprod(Xwsp) +
                                                         .sparseDiagonal(p, lambda),
                                                       as.matrix(D %*% crossprod(Xwsp,yw)),
                                                       ordering = 0L)$coefficients)) # 35s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
cor(coefs_eigenqrsolveord0_sparse, coefs_qrsolve_sparse) # 1


# eigen's iterative preconditioned conjugate gradient method: FAST, but least squares conjugate gradient via row augmentation even faster
system.time(coefs_eigenCGsolve_dense <- as.vector(solve_dense_CG(D %*% crossprod(Xw) +
                                                                         diag(rep(lambda, p)),
                                                                       D %*% crossprod(Xw,yw),
                                                                 x0=rep(0,p),
                                                                 tol=0,
                                                                 maxit=0,
                                                                 precond = 0L)$coefficients)) # 29s for p=n=10000, 1.23s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
system.time(coefs_eigenCGsolve_sparse <- as.vector(solve_sparse_CG(D %*% crossprod(Xwsp) +
                                                                         .sparseDiagonal(p, lambda),
                                                                       as.matrix(D %*% crossprod(Xwsp,yw)),
                                                                       x0=rep(0,p),
                                                                       tol=0,
                                                                       maxit=0,
                                                                       precond = 0L)$coefficients)) # 0.10s for p=n=10000, 2.66s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
cor(coefs_solve_sparse, coefs_eigenCGsolve_sparse) # 1
cor(coefs_solve_dense, coefs_eigenCGsolve_dense) # 1


# default chol.solve: relatively fast, but not great
system.time(coefs_cholsolve_dense <- as.vector(chol.solve(D %*% crossprod(Xw) +
                                                        diag(rep(lambda, p)),
                                                      D %*% crossprod(Xw,yw),
                                                      tol=2*.Machine$double.eps))) # 23s for n=p=10000, 2.02s for n=p=5000, 2.16s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
system.time(coefs_cholsolve_sparse <- as.vector(chol.solve(D %*% crossprod(Xwsp) +
                                                         .sparseDiagonal(p, lambda),
                                                       D %*% crossprod(Xwsp,yw),
                                                       tol=2*.Machine$double.eps))) # 5s for n=p=10000, 48.75s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
cor(coefs_solve_sparse, coefs_true) # 0.1745552
cor(coefs_solve_dense, coefs_cholsolve_dense) # 1
cor(coefs_solve_dense, coefs_cholsolve_sparse) # 1

# eigen's solve_dense_chol & solve_sparse_chol: FAST for sparse only
system.time(coefs_cholsolve_dense <- as.vector(solve_dense_chol(D %*% crossprod(Xw) +
                                                            diag(rep(lambda, p)),
                                                          pivot=1L)$coefficients)) # 84s for n=p=10000, 7.96s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
system.time(coefs_cholsolve_sparse <- as.vector(solve_sparse_chol(D %*% crossprod(Xwsp) +
                                                             .sparseDiagonal(p, lambda),
                                                           as.matrix(D %*% crossprod(Xwsp,yw)),
                                                           pivot=0L,
                                                           ordering=0L)$coefficients)) # 0.08s for n=p=10000, 18.76s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
cor(coefs_cholsolve_sparse, coefs_true) # 0.1745552
cor(coefs_cholsolve_dense, coefs_cholsolve_sparse) # 1
cor(coefs_cholsolve_dense, coefs_cholsolve_sparse) # 1


# using row augmentation & default Cholesky decomposition -> FAST both for dense & sparse X BUT only square matrices (otherwise crashes)
system.time(coefs_cholsolve_rowaug_dense <- as.vector(chol.solve(Xwaug,
                                                                 ywaug,
                                                                 tol=2*.Machine$double.eps))) # 13s for n=p=10000, 0.90s for p=5000 & n=1000
system.time(coefs_cholsolve_rowaug_sparse <- as.vector(chol.solve(Xwaugsp,
                                                                  ywaug,
                                                           tol=2*.Machine$double.eps))) # 1.6s for n=p=10000, 2.36s for p=5000 & n=1000
cor(coefs_solve_dense, coefs_cholsolve_rowaug_dense) # 1
cor(coefs_solve_dense, coefs_cholsolve_rowaug_sparse) # 1


# using row augmentation & eigen's least squares conjugate gradient -> SUPER FAST BOTH FOR SPARSE & DENSE!! -> BEST FOR US?
system.time(coefs_eigenLSCG_rowaug_dense <- as.vector(solve_dense_LSCG(Xwaug,
                                                                         as.matrix(ywaug),
                                                                         x0=as.matrix(rep(0,ncol(Xwaug))),
                                                                         tol=0,
                                                                         maxit=0L,
                                                                         precond=1L)$coefficients)) # 53s for p=n=10000, 0.42s for p=5000 & n=1000
system.time(coefs_eigenLSCG_rowaug_sparse <- as.vector(solve_sparse_LSCG(Xwaugsp,
                                                                          as.matrix(ywaug),
                                                                          x0=as.matrix(rep(0,ncol(Xwaug))),
                                                                          tol=0,
                                                                          maxit=0L,
                                                                          precond=0L)$coefficients)) # 0.19s for p=n=10000, 0.03s for p=5000 & n=1000
cor(coefs_solve_dense, coefs_eigenLSCG_rowaug_dense) # 1
cor(coefs_solve_dense, coefs_eigenLSCG_rowaug_sparse) # 1



# SOLVE USING OPTIMISED FORMULA FOR p > n, CF. Liu & Li 2016 eqn (28) ####
# -> THIS IS MUCH FASTER WHEN p > n
system.time({DXT <- as(Matrix::tcrossprod(D, Xw), "dgCMatrix"); 
             coefs_solve_largep_dense <- as.vector(DXT %*% base::solve(
               as.matrix(Xw %*% DXT + .sparseDiagonal(nrow(X), lambda)),
               yw,
               tol=2*.Machine$double.eps ) )}) 
# 10.8s for n=p=10000 or 0.2s for n=1000 & p=10000 or 
# 0.78s for n=500 & p=50000
# 0.48s for n=1000 & p=50000 with Intel MKL, cf. Liu & Li 2016 eqn (28)
system.time({DXT <- Matrix::tcrossprod(D, Xwsp); 
             coefs_solve_largep_sparse <- as.vector(DXT %*% base::solve(
               Xwsp %*% DXT + .sparseDiagonal(nrow(X), lambda), 
               yw,
               tol=2*.Machine$double.eps ) )}) # 3.4s for n=p=10000 or 0.01s for n=1000 & p=10000 or 0.06s for n=1000 & p=50000 with Intel MKL, cf. Liu & Li 2016 eqn (28)
# 0.08s for n=500 & p=50000
# coefs_solve_largep_sparse <- pmax(lower, pmin(coefs_solve_largep_sparse, upper))
cor(coefs_solve_dense, coefs_solve_largep_dense) # 1
cor(coefs_solve_dense, coefs_solve_largep_sparse) # 1

# using Eigen Cholesky decomposition:
system.time({
  DXT <- Matrix::tcrossprod(D, Xw)
  
  coefs_solveChol_largep_dense <-
    as.vector(
      DXT %*% solve_dense_chol(
        A = as.matrix(Xw %*% DXT +
                        Diagonal(x=rep(lambda, nrow(X))) ),
        b = yw,
        pivot = 1L
      )$coefficients
    )
}) 
# 8.35s for n=500 & p=50000
# 0.50s for n=1000 & p=50000, cf. Liu & Li 2016 eqn (28)
system.time({
  DXT <-
    tcrossprod(D, Xwsp)
  coefs_solveChol_largep_sparse <-
    as.vector(
      DXT %*% solve_sparse_chol(
        A = Xwsp %*% DXT +
          Diagonal(x=rep(lambda, nrow(X))),
        b = yw,
        pivot = 1L,
        ordering = 0L)$coefficients)}) 
# 0.06s for n=500 & p=50000
# 0.06s for n=1000 & p=50000, cf. Liu & Li 2016 eqn (28)
cor(coefs_solve_largep_dense, coefs_solveChol_largep_dense) # 1
cor(coefs_solve_largep_sparse, coefs_solveChol_largep_sparse) # 1

# using Eigen conjugate gradient solver:
system.time({
  DXT <- Matrix::tcrossprod(D, Xw)
  
  coefs_solveCG_largep_dense <-
    as.vector(
      DXT %*% solve_dense_CG(
        A = as.matrix(Xw %*% DXT +
          Diagonal(x=rep(lambda, nrow(X))) ),
        b = yw,
        x0 =
          rep(0, nrow(X)),
        tol =
          0,
        maxit =
          0L,
        precond = 0L
      )$coefficients
    )
}) 
# 8.23s for n=500 & p=50000
# 0.45s for n=1000 & p=50000, cf. Liu & Li 2016 eqn (28)
system.time({
  DXT <-
    tcrossprod(D, Xwsp)
  coefs_solveCG_largep_sparse <-
    as.vector(
      DXT %*% solve_sparse_CG(
        A = Xwsp %*% DXT +
          Diagonal(x=rep(lambda, nrow(X))),
        yw,
        x0 =
          rep(0, nrow(X)),
        tol =
          0,
        maxit=0,
        precond = 0L)$coefficients)}) # 0.05s for n=1000 & p=50000, cf. Liu & Li 2016 eqn (28)
# 0.06s for n=500 & p=50000
cor(coefs_solve_largep_dense, coefs_solveCG_largep_dense) # 1
cor(coefs_solve_largep_sparse, coefs_solveCG_largep_sparse) # 1


# using Eigen QR solver:
system.time({
  DXT <- Matrix::tcrossprod(D, Xw)
  
  coefs_solveQR_largep_dense <-
    as.vector(
      DXT %*% solve_dense_QR(
        A = as.matrix(Xw %*% DXT +
                        Diagonal(x=rep(lambda, nrow(X))) ),
        b = as.matrix(yw),
        pivot = 1L
      )$coefficients
    )
}) 
# 8.22s for n=500 & p=50000
# 0.54s for n=1000 & p=50000, cf. Liu & Li 2016 eqn (28)
system.time({
  DXT <-
    tcrossprod(D, Xwsp)
  coefs_solveQR_largep_sparse <-
    as.vector(
      DXT %*% solve_sparse_QR(
        A = Xwsp %*% DXT +
          Diagonal(x=rep(lambda, nrow(X))),
        b = as.matrix(yw),
        ordering = 1L)$coefficients)}) # 0.07s for n=1000 & p=50000, cf. Liu & Li 2016 eqn (28)
# 0.06s for n=500 & p=50000
cor(coefs_solve_largep_dense, coefs_solveQR_largep_dense) # 1
cor(coefs_solve_largep_sparse, coefs_solveQR_largep_sparse) # 1


# using Eigen LS conjugate gradient solver:
system.time({
  DXT <- Matrix::tcrossprod(D, Xw)
  
  coefs_solveLSCG_largep_dense <-
    as.vector(
      DXT %*% solve_dense_LSCG(
        A = as.matrix(Xw %*% DXT +
                        Diagonal(x=rep(lambda, nrow(X))) ),
        b = as.matrix(yw),
        x0 =
          as.matrix(rep(0, nrow(X))),
        tol =
          0,
        maxit =
          0L,
        precond = 1L
      )$coefficients
    )
}) 
# 8.25s for n=500 & p=50000
# 0.50s for n=1000 & p=50000, cf. Liu & Li 2016 eqn (28)
system.time({
  DXT <-
    tcrossprod(D, Xwsp)
  coefs_solveLSCG_largep_sparse <-
    as.vector(
      DXT %*% solve_sparse_LSCG(
        A = Xwsp %*% DXT +
          Diagonal(x=rep(lambda, nrow(X))),
        b = as.matrix(yw),
        x0 =
          rep(0, nrow(X)),
        tol =
          0,
        maxit=0,
        precond = 1L)$coefficients)}) # 0.06s for n=1000 & p=50000, cf. Liu & Li 2016 eqn (28)
# 0.07s for n=500 & p=50000
cor(coefs_solve_largep_dense, coefs_solveLSCG_largep_dense) # 1
cor(coefs_solve_largep_sparse, coefs_solveLSCG_largep_sparse) # 1

# test that all solvers work
system.time(coefs_eigen1 <- lmridge.eigen(Xw=Xwsp, yw=yw, 
                              lambdas=lambdas, 
                              lower=lower,
                              upper=upper,
                              start=rep(1, ncol(Xwsp)),
                              control.solver=
                                do.call("control.solver.gen", 
                                        list(solver = "eigen", method = 5L))))
# 0.07s for n=500 and p=50000
system.time(coefs_eigen2 <- lmridge.eigen(Xw=Xwsp, yw=yw, 
                              lambdas=lambdas, 
                              lower=lower,
                              upper=upper,
                              start=rep(1, ncol(Xwsp)),
                              control.solver=
                                do.call("control.solver.gen", 
                                        list(solver = "eigen", method = 5L))))
# 0.06s for n=500 and p=50000
system.time(coefs_eigen3 <- lmridge.eigen(Xw=Xwsp, yw=yw, 
                              lambdas=lambdas, 
                              lower=lower,
                              upper=upper,
                              start=rep(1, ncol(Xwsp)),
                              control.solver=
                                do.call("control.solver.gen", 
                                        list(solver = "eigen", method = 4L))))
# 0.14s for n=500 and p=50000
system.time(coefs_eigen4 <- lmridge.eigen(Xw=Xwsp, yw=yw, 
                              lambdas=lambdas, 
                              lower=lower,
                              upper=upper,
                              start=rep(1, ncol(Xwsp)),
                              control.solver=
                                do.call("control.solver.gen", 
                                        list(solver = "eigen", method = 4L))))
# 0.02s for n=500 and p=50000 - BEST for sparse & dense case and p > n
system.time(coefs_eigen5 <- lmridge.eigen(Xw=Xwsp, yw=yw, 
                              lambdas=lambdas, 
                              lower=lower,
                              upper=upper,
                              start=rep(1, ncol(Xwsp)),
                              control.solver=
                                do.call("control.solver.gen", 
                                        list(solver = "eigen", method = 3L))))
# 7.79s for n=500 and p=50000
system.time(coefs_eigen6 <- lmridge.eigen(Xw=Xwsp, yw=yw, 
                              lambdas=lambdas, 
                              lower=lower,
                              upper=upper,
                              start=rep(1, ncol(Xwsp)),
                              control.solver=
                                do.call("control.solver.gen", 
                                        list(solver = "eigen", method = 3L))))
# 8.0s for n=500 and p=50000
system.time(coefs_eigen7 <- lmridge.eigen(Xw=Xwsp, yw=yw, 
                              lambdas=lambdas, 
                              lower=lower,
                              upper=upper,
                              start=rep(1, ncol(Xwsp)),
                              control.solver=
                                do.call("control.solver.gen", 
                                        list(solver = "eigen", method = 1L))))
# 120s for n=500 and p=50000
system.time(coefs_eigen8 <- lmridge.eigen(Xw=Xwsp, yw=yw, 
                              lambdas=lambdas, 
                              lower=lower,
                              upper=upper,
                              start=rep(1, ncol(Xwsp)),
                              control.solver=
                                do.call("control.solver.gen", 
                                        list(solver = "eigen", method = 2L))))
# 116s for n=500 and p=50000
system.time(coefs_eigend1 <- lmridge.eigen(Xw=Xw, yw=yw, 
                              lambdas=lambdas, 
                              lower=lower,
                              upper=upper,
                              start=rep(1, ncol(Xwsp)),
                              control.solver=
                                do.call("control.solver.gen", 
                                        list(solver = "eigen", method = 5L))))
# 7s for n=500 and p=50000
system.time(coefs_eigend2 <- lmridge.eigen(Xw=Xw, yw=yw, 
                              lambdas=lambdas, 
                              lower=lower,
                              upper=upper,
                              start=rep(1, ncol(Xwsp)),
                              control.solver=
                                do.call("control.solver.gen", 
                                        list(solver = "eigen", method = 5L))))
# 7s for n=500 and p=50000
system.time(coefs_eigend3 <- lmridge.eigen(Xw=Xw, yw=yw, 
                              lambdas=lambdas, 
                              lower=lower,
                              upper=upper,
                              start=rep(1, ncol(Xwsp)),
                              control.solver=
                                do.call("control.solver.gen", 
                                        list(solver = "eigen", method = 4L))))
# 0.06s for n=500 and p=50000 - very good for dense case and p > n
system.time(coefs_eigend4 <- lmridge.eigen(Xw=Xw, yw=yw, 
                              lambdas=lambdas, 
                              lower=lower,
                              upper=upper,
                              start=rep(1, ncol(Xwsp)),
                              control.solver=
                                do.call("control.solver.gen", 
                                        list(solver = "eigen", method = 4L))))
# 0.08s for n=500 and p=50000 - BEST for sparse & dense case and p > n
system.time(coefs_eigend5 <- lmridge.eigen(Xw=Xw, yw=yw, 
                              lambdas=lambdas, 
                              lower=lower,
                              upper=upper,
                              start=rep(1, ncol(Xwsp)),
                              control.solver=
                                do.call("control.solver.gen", 
                                        list(solver = "eigen", method = 3L))))
# 8s for n=500 and p=50000
system.time(coefs_eigend6 <- lmridge.eigen(Xw=Xw, yw=yw, 
                              lambdas=lambdas, 
                              lower=lower,
                              upper=upper,
                              start=rep(1, ncol(Xwsp)),
                              control.solver=
                                do.call("control.solver.gen", 
                                        list(solver = "eigen", method = 3L))))
# 8s for n=500 and p=50000
system.time(coefs_eigend7 <- lmridge.eigen(Xw=Xw, yw=yw, 
                              lambdas=lambdas, 
                              lower=lower,
                              upper=upper,
                              start=rep(1, ncol(Xwsp)),
                              control.solver=
                                do.call("control.solver.gen", 
                                        list(solver = "eigen", method = 1L))))
# 115s for n=500 and p=50000
system.time(coefs_eigend8 <- lmridge.eigen(Xw=Xw, yw=yw, 
                              lambdas=lambdas, 
                              lower=lower,
                              upper=upper,
                              start=rep(1, ncol(Xwsp)),
                              control.solver=
                                do.call("control.solver.gen", 
                                        list(solver = "eigen", method = 2L))))
# 116s for n=500 and p=50000



# TESTS WITH ALREADY IMPLEMENTED SOLVERS ####

# PS I think we best make 1 meta lmridge_solve_eigen R function that comprises all the eigen solvers with common options
# (some ignored for particular solvers) and that we would then have an argument solver="eigen" with
# extra option control.solvers with the desired method & a list with the desired options
# TO DO: Check which of the available solvers allow for the outcome y to be a matrix, ie allowing for multichannel data
# (would be good for multitask extension with grouped variable selection, though due to the weights being different
# we would not have the same covariate matrix for all channels, only the sparsity structure would be the same)
# Vectorizing

# CG
system.time(coefs_QR <- lmridge(X=Xsp, y, lambda, weights, lower, upper, nonnegative,
                                solver="eigen", method = 3L)) # CG, 3.0s for n=p=10000
lines(coefs_QR, col="blue", type="h")
cor(coefs_QR, coefs_solve_dense) # 0.9990464
plot(coefs_QR, coefs_solve_sparse)


# LSCG: VERY FAST!
system.time(coefs_LSCG <- lmridge_eigen(X=Xsp, y=y, weights=weights, lambdas=rep(1,p),
                          lower=lower, upper=upper, control.solver=control.solver.gen("eigen", 3))) # 0.18s for n=p=10000

system.time(coefs_LSCG <- lmridge(X=Xsp, y, lambda, weights, lower, upper, nonnegative,
                                  solver="eigen", method=4L)) # 0.17s TODO make sure that start values can be passed OK via arg start (x0 in Rcpp function)
lines(coefs_LSCG, col="blue", type="h")
cor(coefs_LSCG, coefs_solve_sparse) # 0.9990464
plot(coefs_LSCG, coefs_solve_sparse)


system.time(coefs_eigenLSCG_rowaug_sparse <- as.vector(solve_sparse_LSCG(Xwaugsp,
                                                                         as.matrix(ywaug),
                                                                         x0=as.matrix(rep(0,ncol(Xwaug))),
                                                                         tol=0,
                                                                         maxit=0L,
                                                                         precond=0L)$coefficients)) # 0.19s for p=n=10000, 0.03s for p=5000 & n=1000




# glmnet (unbounded): FAST, but LSCG still a bit faster if there is no box constraints??
# also allows option with desired maximum nr of variables

# model <- bigGlm(x = Xaug, y = yaug,
# alpha = 0, lambda = 0, intercept = FALSE, standardize = FALSE, standardize.response = FALSE,
# lower.limits = c(-Inf, lower), upper.limits = c(Inf, upper),
# family = "gaussian", path = FALSE)
# coefs <- coef(model)[-1]

system.time(coefs_glmnet_unbounded <- lmridge(X=Xsp, y, lambda, weights,
                                              lower=rep(-Inf,p), upper=rep(Inf,p), nonnegative=F,
                                              solver="glmnet")) # 0.56s with p=n=10000 & Intel MKL
lines(coefs_glmnet_unbounded, col="orange", type="h")
cor(coefs_glmnet_unbounded, coefs_solve_sparse) # 1


# solvers with box or nonnegativity constraints
# glmnet (bounded): FAST, but LSCG still a bit faster?? Fastest method with explicit box constraints though.
system.time(coefs_glmnet <- lmridge(X=Xsp, y, lambda = lambdas, weights,
                                    lower, upper, nonnegative,
                                    solver="glmnet")) # 0.03s with p=n=10000 & IntelMKL

cor(coefs_glmnet, coefs_true)
lines(coefs_glmnet, col="blue", type="h")

plot(coefs_glmnet, coefs_true)

system.time(coefs_glmnet_dense <- lmridge(X=X, y, lambda, weights,
                                          lower, upper, nonnegative,
                                          solver="glmnet")) # 8.3s with p=n=10000 with MRO/IntelMKL


# osqp: FASTER if  n > p & IntelMKL is installed, otherwise rather slow
control.solver = do.call("control.solver.gen", list(solver = "osqp", method = NULL, rowaug = FALSE))
system.time(coefs_osqp <- lmridge(X=Xsp, y, lambda=0, weights, lower, upper, nonnegative,
                                  solver="osqp")) # 0.21s with p=n=10000 with MRO/IntelMKL
cor(coefs_osqp, coefs_true)

lines(coefs_osqp, col="green", type="h")
cor(coefs_osqp, coefs_true) # 0.9999903
plot(coefs_osqp, coefs_glmnet) # correct...


# nnls: SLOW
system.time(coefs_nnls <- lmridge(X=Xsp, y, lambda, weights, lower, upper,
                                  nonnegative, solver="osqp")) # 1841s for n=p=10000
lines(coefs_nnls, col="yellow", type="h")
cor(coefs_nnls, coefs_true) # 0.98

# fastnnls:
system.time(coefs_fastnnls_dense <- fastnnls_dense(Xwaug, ywaug)) # 21s for p=n=10000


# tests with some other solvers
# unbounded ones
# sgd: currently no box constraints though, works if x is a big.matrix though (bigmemory package)
# also allows start values
library(sgd)
system.time(coefs_sgd <- sgd(x = Xwaug, y = ywaug, model="lm",
                             sgd.control=list(method="ai-sgd",
                                              lr="adagrad",
                                              reltol=1E-10,
                                              start=rep(0,p)))$coefficients) # 27s with p=n=10000
cor(coefs_sgd, coefs_solve_sparse) # 0.9953188


# some other options for GLMs with box constraints (returning an object of class glm, check code to see
# how that's done)
# holiglm, using conic solver, also allows option with desired maximum nr of variables
library(holiglm)

# glm.cons from zetadiv -> uses nnnpls from library(nnls) -> beter use nnnpls directly
library(zetadiv)
?glm.cons

# restriktor

# ConsReg




# Case I: fit nonnegative identity link Poisson GLM with no penalty
L0glm1 <- L0glm2(formula = "y ~ X + .", data = data.frame(X, y = y), solver = "osqp", family = poisson(identity),
                lambda = 0, tune.meth = "none", nonnegative = TRUE, normalize = FALSE,
                control.iwls = list(maxit = 100),
                control.l0 = list(maxit = 1, delta = 1E-2, gamma = 1.8))

plot_L0glm_benchmark(x = sim$x, y = y, fit = L0glm1, a.true = sim$a,
                     main="Ground truth vs L0glm estimates")


##################################
##########BLOCK LMRIDGE###########
##################################
block.size <- 500
block.ids <- ceiling((1:p)/block.size)
control.fit <- do.call("control.fit.gen2", list(block.size = block.size, maxit = 1000))
#OSQP
system.time(coefs_osqp <- lmridge(X = Xsp, y = y, lambda = lambda, weights = weights, lower = lower, upper = upper,
                                  nonnegative = nonnegative, solver = "osqp", method = 5L))

  #sparse
system.time(coefs_block_osqp <- block_lmridge(Xsp, y = y, lambda = lambda, weights = weights, lower = lower, upper = upper,
                      nonnegative = nonnegative, solver = "osqp", verbose = verbose,
                      control.fit = control.fit, block.ids = block.ids))

  #dense
system.time(coefs_block_osqp <- block_lmridge(X, y = y, lambda = lambda, weights = weights, lower = lower, upper = upper,
                                              nonnegative = nonnegative, solver = "osqp", verbose = verbose,
                                              control.fit = control.fit, block.ids = block.ids))
cor(coefs_osqp, coefs_block_osqp)

#EIGEN
system.time(coefs_eigen1 <- lmridge(X = Xsp, y = y, lambda = lambda, weights = weights, lower = lower, upper = upper,
                                  nonnegative = nonnegative, solver = "eigen", method = 1L))

system.time(coefs_block_eigen1 <- block_lmridge(Xsp, y = y, lambda = lambda, weights = weights, lower = lower, upper = upper,
                                              nonnegative = nonnegative, solver = "eigen", method = 1L, verbose = verbose,
                                              control.fit = control.fit, block.ids = block.ids))

cor(coefs_eigen1, coefs_block_eigen1)


