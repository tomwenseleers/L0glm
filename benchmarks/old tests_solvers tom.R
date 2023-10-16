# install.packages("Matrix", repos = "https://mran.microsoft.com")
library(nnls)
library(microbenchmark)
# devtools::install_github("nk027/sanic")
library(sanic)
library(osqp)
library(devtools)
library(glmnet)
devtools::load_all()
library(RhpcBLASctl) # to set nr of threads & cores to use in your BLAS
# fix to set working package repository when using Microsoft R Open with Intel MKL
options(repos = c(CRAN="https://cloud.r-project.org/"))

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

# Cholesky decomposition based solve, using normal equations (FAST)
chol.solve = function (X, y, tol=.Machine$double.eps^(2/3)) {
  ch <- base::chol(crossprod(X)) # solve using Cholesky decomposition
  backsolve(ch, forwardsolve(ch, crossprod(X, y), upper = TRUE, trans = TRUE))
}

library(Matrix)

# Simulate data
set.seed(1)
n <- 10000
p <- 10000
s <- 0.1 # sparsity as proportion of p that have nonzero coefficients // p*s == #of peaks
k <- round(p*s) # nr of nonzero covariates
sim <- simulate_spike_train(n=n, p=p, k=k, mean_beta = 1000, sd_logbeta = 1, family="poisson")

system.time(Matrix::crossprod(sim$X)) # 0.02s for n=p=10000

coefs_true <- as.vector(sim$beta_true)
Xsp <- sim$X
X.n <- apply(Xsp, 2, norm, type = "2") # slow, maybe use max(abs(columns of X)) ?
# X.n <- rep(1, ncol(Xsp))
Xsp <- sweep(Xsp, 2, X.n, "/") # slow
X <- as.matrix(Xsp)

lambda <- 1 # ridge penalty
y <- as.numeric(as.vector(sim$y))
weights <- 1/(y+0.1) # observation weights
lower <- rep(0, ncol(X))
upper <- rep(Inf, ncol(X))
nonnegative <- TRUE
thresh <- 1E-3 # coefficients below this value are set to zero in each iteration
beta <- start <- rep(1, ncol(X))

# # we use starting values starting from a nonnegative weighted LS fit (fit using osqp)
# system.time(start <- lmridge(X=Xsp, y, lambda = rep(1E-5, ncol(Xsp)), weights = weights,
#                              lower = lower, upper = upper, nonnegative = nonnegative,
#                              solver="osqp") / X.n) # 4.64s for n=p=10000
# start[abs(start)<thresh] <- 0
# nonzero <- (abs(start)>0)
# 
# cor(coefs_true, start) # 0.69
# plot(log10(coefs_true), log10(start))


delta <- 1E-4 # not needed if thresholding is used
gamma <- 2 # gamma=1.5 sometimes more stable
nonzero <- beta!=0
sum(nonzero)
pen_weights <- rep(0, ncol(X))
pen_weights[nonzero] <- (1/(start[nonzero]^gamma+delta)) # penalization weights of different variables (columns of X)
pen_weights[nonzero] <- pen_weights[nonzero] / (sum(pen_weights[nonzero])/sum(nonzero)) # normalize sum of pen weights to ncol(X[,nonzero,drop=F])
hist(log10(pen_weights[nonzero]))
sum(pen_weights[nonzero]) == sum(nonzero) # TRUE
lambdas <- rep(0, ncol(X))
lambdas[nonzero] <- lambda*pen_weights[nonzero]
D <- diag(1/pen_weights[nonzero]) # cf notation Liu & Li 2016
Dsp <- Matrix(D, sparse=T) # .sparseDiagonal(ncol(X), 1/pen_weights)  # sparse version of D, cf notation Liu & Li 2016
Xw <- as.matrix(X[,nonzero,drop=F]*sqrt(weights))
Xwaug <- rbind(Xw,diag(sqrt(lambdas[nonzero])))
Xwsp <- Matrix(Xw, sparse=TRUE)
Xwaugsp <- Matrix(Xwaug, sparse=TRUE)
yw <- y*sqrt(weights)
ywaug <- c(yw, rep(0, sum(nonzero)))
# start <- NULL
verbose <- FALSE







# 1. TESTS WITH SOME DEFAULT SOLVE METHODS TO SOLVE RIDGE REGRESSION, EITHER WITH DIRECT SOLVE OR VIA ROW AUGMENTATION ####
# 1.1 TEST WITH DIFFERENT SOLVE METHODS ####

# (using default BLAS that R is compiled against, 
# can be Accelerate (Mac), OpenBLAS (Ubuntu/Win) or IntelMKL (Ubuntu/Win))
# SOLVE USING BASIC FORMULA FOR n >= p, CF. Liu & Li 2016 ####

# with base R solve which uses LAPACK routines DGESV and ZGESV which use LU decomposition XXXX
system.time(coefs_solve_dense <- as.vector(base::solve(D %*% crossprod(Xw) + 
                                                    diag(rep(lambda, sum(nonzero))),
                                                    D %*% crossprod(Xw,yw),
                                                    tol=2*.Machine$double.eps))) 
# 13s for n=p=10000 with Intel MKL, cf. Liu & Li 2016 eqn (10)
# this is faster than eigen Cholesky for dense case below
system.time(coefs_solve_sparse <- as.vector(base::solve(Dsp %*% crossprod(Xwsp) + # FAST
                                                     .sparseDiagonal(sum(nonzero), lambda),
                                                     as.matrix(Dsp %*% crossprod(Xwsp,yw)),
                                                     tol=2*.Machine$double.eps))) 
# 3.59s for n=p=10000 with Intel MKL, cf. Liu & Li 2016 eqn (10)
cor(coefs_solve_dense, coefs_solve_sparse) # 1

# default chol.solve: FAST but only with Intel MKL & solve_sparse_CG faster
system.time(coefs_cholsolve_dense <- as.vector(chol.solve(D %*% crossprod(Xw) +
                                                            diag(rep(lambda, ncol(Xw))),
                                                          D %*% crossprod(Xw,yw),
                                                          tol=2*.Machine$double.eps))) 
# 15s for n=p=10000 with Intel MKL, cf. Liu & Li 2016 eqn (10)
system.time(coefs_cholsolve_sparse <- as.vector(chol.solve(Dsp %*% crossprod(Xwsp) +
                                                             .sparseDiagonal(sum(nonzero), lambda),
                                                           Dsp %*% crossprod(Xwsp,yw),
                                                           tol=2*.Machine$double.eps))) 
# 1.58s for n=p=10000 with Intel MKL, cf. Liu & Li 2016 eqn (10)
# Sometimes gives Error in chol.default(crossprod(X)) :
# the leading minor of order 3843 is not positive definite, so overfitted, 
# some eigenvectors near zero

# alt syntax (more stable) :
system.time(coefs_cholsolve_sparse <- as.vector(chol.solve(as(crossprod(Xwaugsp),"dgCMatrix"),
                                                           crossprod(Xwaugsp,ywaug),
                                                           tol=2*.Machine$double.eps))) 
# 1.5s for n=p=10000 
cor(coefs_cholsolve_sparse, coefs_solve_sparse) # 1

# eigen's solve_dense_chol & solve_sparse_chol: 
# fast only for sparse - check better solution for dense? (base LU solver faster then)
# very slow for dense case
system.time(coefs_cholsolve_dense <- as.vector(solve_dense_chol(A=D %*% crossprod(Xw) +
                                                                  diag(rep(lambda, ncol(Xw))),
                                                                b=as.matrix(D %*% crossprod(Xw,yw)),
                                                                pivot=0L)$coefficients)) 
# 71s for n=p=10000 with Intel MKL, cf. Liu & Li 2016 eqn (10)
cor(coefs_cholsolve_dense, coefs_solve_sparse) # 1
system.time(coefs_cholsolve_sparse <- as.vector(solve_sparse_chol(A=Dsp %*% crossprod(Xwsp) +
                                                                    .sparseDiagonal(sum(nonzero), lambda),
                                                                  b=as.matrix(Dsp %*% crossprod(Xwsp,yw)),
                                                                  pivot=0L,
                                                                  ordering=0L)$coefficients)) 
# 0.05s for n=p=10000 with Intel MKL, cf. Liu & Li 2016 eqn (10)
cor(coefs_cholsolve_sparse, coefs_solve_sparse) # 1


# more stable fit:
# FASTEST CHOLESKY SOLVER FOR SPARSE CASE : via row augmentation
system.time(coefs_cholsolve_sparse <- as.vector(solve_sparse_chol(as(crossprod(Xwaugsp),"dgCMatrix"),
                                                                  as.matrix(crossprod(Xwaugsp,ywaug)),
                                                                  pivot=0L,
                                                                  ordering=0L)$coefficients)) 
# 0.03s for n=p=10000, cf. Liu & Li 2016 eqn (10)
cor(coefs_cholsolve_sparse, coefs_solve_sparse) # 1

# using row augmentation & default Cholesky decomposition -> slower than lmridge eigen solver
system.time(coefs_cholsolve_rowaug_dense <- as.vector(chol.solve(Xwaug,
                                                                 ywaug,
                                                                 tol=2*.Machine$double.eps))) 
# 7.83s for n=p=10000 with Intel MKL, slower than lmridge eigen above
system.time(coefs_cholsolve_rowaug_sparse <- as.vector(chol.solve(Xwaugsp,
                                                                  ywaug,
                                                                  tol=2*.Machine$double.eps))) 
# 1.24s for n=p=10000 with Intel MKL, slower than lmridge eigen above
cor(coefs_cholsolve_rowaug_sparse, coefs_solve_sparse) # 1




# system.time(coefs_qrsolve_dense <- as.vector(qr.solve(D %*% crossprod(Xw) +
#                                                    diag(rep(lambda, ncol(X))),
#                                                  D %*% crossprod(Xw,yw),
#                                                  tol=2*.Machine$double.eps))) # 386.13s for n=p=10000, 27.47s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
# system.time(coefs_qrsolve_sparse <- as.vector(qr.solve(Dsp %*% crossprod(Xwsp) +
#                                                     .sparseDiagonal(sum(nonzero), lambda),
#                                                   Dsp %*% crossprod(Xwsp,yw),
#                                                   tol=2*.Machine$double.eps))) # 298s for n=p=3000, 28.35s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
# beta[nonzero] <- coefs_qrsolve_sparse
# cor(beta, coefs_true) # 0.8589209
# cor(coefs_qrsolve_sparse, coefs_solve_sparse) # 1
# cor(coefs_solve_dense, coefs_qrsolve_dense) # 1
# cor(coefs_solve_dense, coefs_qrsolve_sparse) # 1
# 
# 
# # eigen's QR methods: slower than default QR solve methods
# system.time(coefs_eigenqrsolvepivot0_dense <- as.vector(solve_dense_QR(D %*% crossprod(Xw) +
#                                                         diag(rep(lambda, ncol(X))),
#                                                      D %*% crossprod(Xw,yw),
#                                                       pivot = 0L)$coefficients)) # 40.62s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
# system.time(coefs_eigenqrsolveord0_sparse <- as.vector(solve_sparse_QR(Dsp %*% crossprod(Xwsp) +
#                                                          .sparseDiagonal(sum(nonzero), lambda),
#                                                        as.matrix(Dsp %*% crossprod(Xwsp,yw)),
#                                                        ordering = 0L)$coefficients)) # 35s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
# beta[nonzero] <- coefs_eigenqrsolveord0_sparse
# cor(beta, coefs_true) # 0.8589209
# cor(coefs_eigenqrsolveord0_sparse, coefs_solve_sparse) # 1
# cor(coefs_eigenqrsolveord0_sparse, coefs_qrsolve_sparse) # 1



# sparse Cholesky solver when called via lmridge function ####
# EXCELLENT TIMINGS BOTH FOR SPARSE & DENSE CASE
microbenchmark(coefs_cholsolve_sparse_lmridge <- lmridge(X=Xsp, y=y, lambdas=lambdas,
                                                         weights = weights,
                                                         lower = lower,
                                                         upper = upper,
                                                         nonnegative = nonnegative,
                                                         # scaling = rep(1, ncol(X)),
                                                         prior.mean = NULL,
                                                         solver = "eigen",
                                                         method = 7L,
                                                         start = rep(0, ncol(X)),
                                                         control.solver = 
                                                           modifyList(do.call("control.solver.gen", 
                                                                              list(solver = "eigen", 
                                                                                   method = 7L)),
                                                                      list(rowaug=TRUE)) )
)
# with n=p=10000 for sparse case: 0.058s with rowaug=Fv, 0.042s with rowaug=T
cor(coefs_cholsolve_sparse_lmridge, coefs_solve_sparse) # 0.9965996
# TODO
system.time(coefs_cholsolve_sparse <- as.vector(solve_sparse_chol(as(crossprod(Xwaugsp),"dgCMatrix"),
                                                                  as.matrix(crossprod(Xwaugsp,ywaug)),
                                                                  pivot=0L,
                                                                  ordering=0L)$coefficients))
# gets to 0.03s for n=p=10000, so should still be possible to get this slightly faster?
cor(coefs_cholsolve_sparse, coefs_solve_sparse) # 1

# dense Cholesky solver when called via lmridge function
system.time(coefs_cholsolve_dense_lmridge <- lmridge(X=X, y=y, lambdas=lambdas,
                                                     weights = weights,
                                                     lower = lower,
                                                     upper = upper,
                                                     nonnegative = nonnegative,
                                                     # scaling = rep(1, ncol(X)),
                                                     prior.mean = NULL,
                                                     solver = "eigen",
                                                     method = 7L,
                                                     start = rep(0, ncol(X)),
                                                     control.solver = 
                                                       modifyList(do.call("control.solver.gen", 
                                                                          list(solver = "eigen", 
                                                                               method = 7L)),
                                                                  list(rowaug=FALSE)) )
)
# with n=p=10000 for dense case: 4.08s with rowaug=TRUE, 4.03s with rowaug=FALSE
# TODO -> set rowaug=T only when n>=p & matrix is sparse?
cor(coefs_cholsolve_sparse_lmridge, coefs_cholsolve_sparse) # 0.9965996




# eigen's iterative preconditioned conjugate gradient method: 
# FAST, but least squares conjugate gradient via row augmentation even faster
system.time(coefs_eigenCGsolve_dense <- as.vector(solve_dense_CG(D %*% crossprod(Xw) +
                                                                         diag(rep(lambda, ncol(Xw))),
                                                                       D %*% crossprod(Xw,yw),
                                                                 x0=rep(0,ncol(Xw)),
                                                                 tol=0,
                                                                 maxit=0,
                                                                 precond = 0L)$coefficients)) 
# 15s for p=n=10000, 1.23s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
# slower than lmridge cholesky though
system.time(coefs_eigenCGsolve_sparse <- as.vector(solve_sparse_CG(Dsp %*% crossprod(Xwsp) +
                                                                         .sparseDiagonal(sum(nonzero), lambda),
                                                                       as.matrix(Dsp %*% crossprod(Xwsp,yw)),
                                                                       x0=rep(0,ncol(Xw)), # as.matrix(start[nonzero]), # rep(0,ncol(X)),
                                                                       tol=0,
                                                                       maxit=0,
                                                                       precond = 1L)$coefficients)) 
# 0.05s for p=n=10000, 2.66s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
# slower than lmridge cholesky
coefs_eigenCGsolve_sparse
# sometimes returns NULL
cor(coefs_eigenCGsolve_sparse, coefs_solve_sparse) # 1

# using crossproduct of row augmented cov matrix: SUPER FAST
# -> faster than using normal equations above for sparse case
system.time(coefs_eigenCGsolve_sparse <- as.vector(solve_sparse_CG(A=as(crossprod(Xwaugsp), "dgCMatrix"),
                                                                   b=as.matrix(crossprod(Xwaugsp,ywaug)),
                                                                   x0=rep(0,ncol(Xw)), # as.matrix(start[nonzero]), # rep(0,ncol(X)),
                                                                   tol=1E-10,
                                                                   maxit=0,
                                                                   precond = 1L)$coefficients)) 
# 0.03s for n=p=10000
cor(coefs_eigenCGsolve_sparse, coefs_solve_sparse) # 1

system.time(coefs_eigenCGsolve_dense <- as.vector(solve_dense_CG(A=as.matrix(crossprod(Xwaug)),
                                                                 b=as.matrix(crossprod(Xwaug,ywaug)),
                                                                   x0=rep(0,ncol(Xw)), # as.matrix(start[nonzero]), # rep(0,ncol(X)),
                                                                   tol=1E-10,
                                                                   maxit=0,
                                                                   precond = 0L)$coefficients)) 
# 9.8s for n=p=10000


# using row augmentation & eigen's least squares conjugate gradient 
# -> SUPER FAST BOTH FOR SPARSE & DENSE!! -> BEST FOR US?
system.time(coefs_eigenLSCG_rowaug_dense <- as.vector(solve_dense_LSCG(Xwaug,
                                                                         as.matrix(ywaug),
                                                                         x0=as.matrix(rep(0,ncol(Xwaug))),
                                                                         tol=0,
                                                                         maxit=0L,
                                                                         precond=0L)$coefficients)) 
# 21s for p=n=10000, CHOLESKY & CG FASTER
cor(coefs_eigenLSCG_rowaug_dense, coefs_solve_sparse) # 1

system.time(coefs_eigenLSCG_rowaug_dense2 <- as.vector(solve_sparse_LSCG(rbind(Xw, Diagonal(x=sqrt(lambdas))),
                                                                         as.matrix(ywaug),
                                                                         x0=as.matrix(rep(0,ncol(Xwaug))), # as.matrix(start[nonzero]), # as.matrix(rep(0,ncol(Xwaug))),
                                                                         tol=0,
                                                                         maxit=0L,
                                                                         precond=0L)$coefficients)) 
# 0.28s for p=n=10000, SUPER FAST!!!! - TODO code it like this in lmridge.eigen
cor(coefs_eigenLSCG_rowaug_dense2, coefs_solve_sparse) # 1

system.time(coefs_eigenLSCG_rowaug_sparse <- as.vector(solve_sparse_LSCG(Xwaugsp,
                                                                          as.matrix(ywaug),
                                                                          x0=as.matrix(rep(0,ncol(Xwaug))), # as.matrix(start[nonzero]), # as.matrix(rep(0,ncol(Xwaug))),
                                                                          tol=0,
                                                                          maxit=0L,
                                                                          precond=1L)$coefficients)) 
# 0.03s for p=n=10000 - SUPER FAST
# NOTE: sometimes returns NULL
cor(coefs_eigenLSCG_rowaug_sparse, coefs_solve_sparse) # 1



# 1.2. SOLVE USING OPTIMISED FORMULA FOR p > n, CF. Liu & Li 2016 eqn (28) ####
# -> THIS IS INDEED FASTER WHEN p > n, ALTHOUGH solve_dense_LSCG & solve_sparse_LSCG STILL FASTEST
system.time({DXT <- tcrossprod(D, Xw); coefs_solve_largep_dense <- as.vector(DXT %*% base::solve(Xw %*% DXT +
                                                                             diag(rep(lambda, nrow(X))),
                                                                             yw,
                                                                             tol=2*.Machine$double.eps))}) 
# 14s for p=n=10000 with IntelMKL, cf. Liu & Li 2016 eqn (28)
cor(coefs_solve_largep_dense, coefs_solve_sparse) # 1

system.time({DXT <- tcrossprod(Dsp, Xwsp); coefs_solve_largep_sparse <- as.vector(DXT %*% base::solve(Xwsp %*% DXT +
                                                                                    .sparseDiagonal(nrow(X), lambda),
                                                                                  yw,
                                                                                  tol=2*.Machine$double.eps))}) 
# 3.2s for p=n=10000 with IntelMKL, cf. Liu & Li 2016 eqn (28)
cor(coefs_solve_largep_sparse, coefs_solve_sparse) # 1

system.time({
  DXT <-
    tcrossprod(D, Xw)
  coefs_solveCG_largep_dense <-
    as.vector(
      DXT %*% solve_dense_CG(
        Xw %*% DXT +
          diag(rep(lambda, nrow(X))),
        yw,
        x0 =
          rep(0, nrow(X)),
        tol =
          0,
        maxit =
          0,
        precond = 0L)$coefficients)}) 
# 16.47s for p=n=10000, cf. Liu & Li 2016 eqn (28)
cor(coefs_solveCG_largep_dense, coefs_solve_sparse) # 1

system.time({
  DXT <- tcrossprod(Dsp, Xwsp)
  
  coefs_solveCG_largep_sparse <-
    as.vector(
      DXT %*% solve_sparse_CG(
        Xwsp %*% DXT +
          .sparseDiagonal(nrow(X), lambda),
        yw,
        x0 =
          rep(0, nrow(X)),
        tol =
          0,
        maxit =
          0,
        precond = 0L
      )$coefficients
    )
}) 
# 0.05s for p=n=10000 - FAST,- cf. Liu & Li 2016 eqn (28)
cor(coefs_solveCG_largep_sparse, coefs_solve_sparse) # 1




# 1.3. TESTS WITH ALREADY IMPLEMENTED SOLVERS ####

# PS I think we best make 1 meta lmridge_solve_eigen R function that comprises all the eigen solvers with common options
# (some ignored for particular solvers) and that we would then have an argument solver="eigen" with
# extra option control.solvers with the desired method & a list with the desired options

## 1.3.1 FOR SPARSE CASE ####
### 1.3.1.1 CG ####
system.time(coefs_CG <- lmridge(X=Xsp[,nonzero,drop=F], 
                                y=y, 
                                lambdas=lambdas[nonzero],
        weights = weights,
        lower = lower[nonzero],
        upper = upper[nonzero],
        nonnegative = nonnegative,
        # scaling = rep(1, ncol(X)),
        prior.mean = NULL,
        solver = "eigen",
        method = 5L,
        start = rep(0, ncol(X))[nonzero],
        control.solver = 
          modifyList(do.call("control.solver.gen", 
                             list(solver = "eigen", 
                                  method = 5L)),
                     list(rowaug=FALSE)) )) # best with rowaug=FALSE here for sparse case?
# 0.05s for n=p=10000
cor(coefs_CG, coefs_solve_sparse) # 0.9965962
# TODO make sure that start values can be passed OK via arg start (x0 in Rcpp function)
# TODO if solver fails (returns NULL as coefficients) make it use cholesky LLT instead

# TODO
system.time(coefs_cgsolve_sparse <- as.vector(solve_sparse_CG(A=as(crossprod(Xwaugsp),"dgCMatrix"),
                                                              b=as.matrix(crossprod(Xwaugsp,ywaug)),
                                                              x0=rep(0, ncol(X))[nonzero],
                                                              tol=0,
                                                              maxit=0L,
                                                              precond=0L)$coefficients))
# gets to 0.03s for n=p=10000, so should still be possible to get this slightly faster?



### 1.3.1.2 LSCG ####
system.time(coefs_LSCG <- lmridge(X=Xsp[,nonzero,drop=F], 
                                y=y, 
                                lambdas=lambdas[nonzero],
                                weights = weights,
                                lower = lower[nonzero],
                                upper = upper[nonzero],
                                nonnegative = nonnegative,
                                # scaling = rep(1, ncol(X)),
                                prior.mean = NULL,
                                solver = "eigen",
                                method = 3L,
                                start = rep(0, ncol(X))[nonzero],
                                control.solver = 
                                  modifyList(do.call("control.solver.gen", 
                                                     list(solver = "eigen", 
                                                          method = 3L)),
                                             list(rowaug=TRUE)) )) # fastest with rowaug=TRUE - TODO set rowaug always to TRUE with this solver
# 0.01s for n=p=10000
cor(coefs_LSCG, coefs_solve_sparse) # 0.9965962
# TODO make sure that start values can be passed OK via arg start (x0 in Rcpp function)
# TODO if solver fails (returns NULL as coefficients) make it use cholesky LLT instead



### 1.3.1.3 Cholesky ####
system.time(coefs_chol <- lmridge(X=Xsp[,nonzero,drop=F], 
                                  y=y, 
                                  lambdas=lambdas[nonzero],
                                  weights = weights,
                                  lower = lower[nonzero],
                                  upper = upper[nonzero],
                                  nonnegative = nonnegative,
                                  # scaling = rep(1, ncol(X)),
                                  prior.mean = NULL,
                                  solver = "eigen",
                                  method = 7L,
                                  start = rep(0, ncol(X))[nonzero],
                                  control.solver = 
                                    modifyList(do.call("control.solver.gen", 
                                                       list(solver = "eigen", 
                                                            method = 7L)),
                                               list(rowaug=TRUE)) )) # fastest with rowaug=TRUE here
# 0.05s for n=p=10000
cor(coefs_chol, coefs_solve_sparse) # 0.9965996
# TODO
system.time(coefs_cholsolve_sparse <- as.vector(solve_sparse_chol(as(crossprod(Xwaugsp),"dgCMatrix"),
                                                                  as.matrix(crossprod(Xwaugsp,ywaug)),
                                                                  pivot=0L,
                                                                  ordering=0L)$coefficients))
# gets to 0.03s for n=p=10000, so should still be possible to get this slightly faster?


### 1.3.1.4 QR ####
system.time(coefs_qr <- lmridge(X=Xsp[,nonzero,drop=F], 
                                  y=y, 
                                  lambdas=lambdas[nonzero],
                                  weights = weights,
                                  lower = lower[nonzero],
                                  upper = upper[nonzero],
                                  nonnegative = nonnegative,
                                  # scaling = rep(1, ncol(X)),
                                  prior.mean = NULL,
                                  solver = "eigen",
                                  method = 1L,
                                  start = rep(0, ncol(X))[nonzero],
                                  control.solver = 
                                    modifyList(do.call("control.solver.gen", 
                                                       list(solver = "eigen", 
                                                            method = 1L)),
                                               list(rowaug=FALSE)) )) # fastest with rowaug=FALSE for sparse case
# 0.1s for n=p=10000
cor(coefs_qr, coefs_solve_sparse) # 0.9965996

# TODO maybe still add an eigen LU & SVD based solver too?


### 1.3.1.5 osqp ####
system.time(coefs_osqp <- lmridge(X=Xsp[,nonzero,drop=F], 
                                y=y, 
                                lambdas=lambdas[nonzero],
                                weights = weights,
                                lower = lower[nonzero],
                                upper = upper[nonzero],
                                nonnegative = nonnegative,
                                # scaling = rep(1, ncol(X)),
                                prior.mean = NULL,
                                solver = "osqp",
                                start = rep(0, ncol(X))[nonzero],
                                control.solver = 
                                  modifyList(do.call("control.solver.gen", 
                                                     list(solver = "osqp")),
                                             list(rowaug=TRUE,
                                                  eps_abs=1E-8,
                                                  eps_rel=1E-8)) )) # fastest with rowaug=TRUE - always set to TRUE with this solver?
# 0.11s for n=p=10000
cor(coefs_osqp, coefs_solve_sparse) # 0.9844926


### 1.3.1.6 glmnet / bigglm ####
system.time(coefs_bigglm <- lmridge(X=Xsp[,nonzero,drop=F], 
                                y=y, 
                                lambdas=lambdas[nonzero],
                                weights = weights,
                                lower = lower[nonzero],
                                upper = upper[nonzero],
                                nonnegative = nonnegative,
                                # scaling = rep(1, ncol(X)),
                                prior.mean = NULL,
                                solver = "glmnet",
                                start = rep(0, ncol(X))[nonzero],
                                control.solver = 
                                  modifyList(do.call("control.solver.gen", 
                                                     list(solver = "glmnet")),
                                             list(thresh=1E-7)) )) # quite a bit faster with lower thresh
# 0.01s for n=p=10000
cor(coefs_bigglm, coefs_solve_sparse) # 0.9844928
cor(coefs_bigglm, coefs_osqp) # 1

# TO DO: solve by using inbuilt ridge in glmnet with alpha=0, lambda=lambda/nrow(X) & appropriate penalty weights
# instead of using row augmentation & lambda=0 to put in ridge penalty
# now it uses
system.time(model <- bigGlm(x = Xwaugsp, y = ywaug,
  alpha = 0, lambda = 0, intercept = FALSE, 
  standardize = FALSE, 
  standardize.response = FALSE,
  lower.limits = c(-Inf, lower), 
  upper.limits = c(Inf, upper),
  family = "gaussian", path = FALSE))
# 0.01s
coefs <- coef(model)[-1]
cor(coefs, coefs_solve_sparse) # 0.9844928


### 1.3.1.7 nnls ####
# SUPER SLOOOOW....
system.time(coefs_nnls <- lmridge(X=Xsp[,nonzero,drop=F], y,
                                  lambdas[nonzero],
                                  weights,
                                  lower=rep(0,sum(nonzero)),
                                  upper=rep(Inf,sum(nonzero)),
                                  nonnegative=T, solver="nnls")) 
# 1841s for n=p=10000 (31 mins), ie 184100x slower than LSCG
cor(coefs_nnls, coefs_osqp) # 0.9712636


## tests with some other solvers ####
# unbounded ones

# parglm ####
# https://htmlpreview.github.io/?https://github.com/boennecd/parglm/blob/master/vignettes/parglm.html
library(parglm)
system.time(parglm_dense <- parglm.fit(
  x = Xwaug, 
  y = ywaug, 
  intercept = F,
  family = gaussian(identity), 
  control = parglm.control(maxit = 1L, 
                           nthreads=1L, # using more threads hardly reduces timing
                           # block_size=1000, 
                           method="FAST")) )
# 9s for n=p=10000


# speedglm ####
# https://cran.r-project.org/web/packages/speedglm/speedglm.pdf
library(speedglm)
system.time(spglm_dense <- speedglm.wfit(
  X = Xwaug, 
  y = ywaug, 
  intercept = F,
  # row.chunk= XX,
  family = gaussian(identity), 
  chunksize = 1000, 
  sparse = FALSE,
  maxit = 1L))
# 19.91s for n=p=10000

system.time(spglm_sparse <- speedglm.wfit(
  X = Xwaugsp, 
  y = ywaug, 
  intercept = F,
  # row.chunk= XX,
  family = gaussian(identity), 
  chunksize = 1000, 
  sparse = TRUE,
  maxit = 1L))
# 12.28s for n=p=10000


# regular glm.fit ####
system.time(plainglm <- glm.fit(
  x = Xwaug, 
  y = ywaug, 
  family = gaussian(identity), 
  control = glm.control(maxit = 1)))
# 208s for n=p=10000


# fastglm (here using Eigen's Cholesky LLT solver) ####
# slow as only dense solvers have been implemented
# and appears to be slow even for dense case
library(fastglm)
system.time(fglm1 <- fastglmPure(
                            x = Xwaug, 
                            y = ywaug, 
                            family = gaussian(identity), 
                            method = 2L, # LLT Cholesky
                            maxit = 1L))
# 237s for n=p=10000
system.time(fglm2 <- fastglmPure(
  x = Xwaug, 
  y = ywaug, 
  family = gaussian(identity), 
  method = 3L, # LDLT Cholesky
  maxit = 1L))
# 347s for n=p=10000



# sgd: currently no box constraints; works if x is a big.matrix though (bigmemory package)
# also allows start values
library(sgd)
system.time(coefs_sgd <- sgd(x = Xwaug, y = ywaug, model="lm",
                             sgd.control=list(method="ai-sgd",
                                              lr="adagrad",
                                              reltol=1E-10,
                                              start=rep(0,ncol(X))))$coefficients) 
# 27s with p=n=10000 - SLOOOOW


# some other options for GLMs with box constraints (returning an object of class glm, check code to see
# how that's done)
# holiglm, using conic solver, also allows option with desired maximum nr of variables
# no support for sparse design matrix though & only formula interface
# library(holiglm)

# glm.cons from zetadiv -> uses nnnpls from library(nnls) -> beter use nnnpls directly
# library(zetadiv)
# ?glm.cons

# restriktor

# ConsReg




# TODO: fix formula interface
L0glm1 <- L0glm(formula = "y ~ X + .", 
                data = list(X = Xsp, # data should be passed as list & allow X to be sparse
                            y = y), 
                solver = "osqp", 
                family = poisson(identity),
                lambda = "gic", 
                nonnegative = TRUE, 
                normalize = FALSE)

# plot_L0glm_benchmark(x = sim$x, y = y, fit = L0glm1, a.true = sim$a,
#                      main="Ground truth vs L0glm estimates")
