# install.packages("Matrix", repos = "https://mran.microsoft.com")
library(nnls)
library(microbenchmark)
# devtools::install_github("nk027/sanic")
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

system.time(Matrix::crossprod(sim$X)) # 0.02s

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

# we use starting values starting from a nonnegative weighted LS fit (fit using osqp)
system.time(start <- lmridge(X=Xsp, y, lambda = rep(1E-5, ncol(Xsp)), weights,
                             lower, upper, nonnegative,
                             solver="osqp") / X.n) # 1.97s
start[abs(start)<thresh] <- 0
nonzero <- (abs(start)>0)
pen_weights <- rep(0, ncol(X))
lambdas <- rep(0, ncol(X))

cor(coefs_true, start) # 0.65
plot(log10(coefs_true), log10(start))


beta <- start
delta <- 1E-4 # not needed if thresholding is used
gamma <- 1.5 # gamma=2 sometimes returns in instability in the fits
pen_weights[nonzero] <- (1/(start[nonzero]^gamma+delta)) # penalization weights of different variables (columns of X)
pen_weights[nonzero] <- pen_weights[nonzero] / (sum(pen_weights[nonzero])/sum(nonzero)) # normalize sum of pen weights to ncol(X[,nonzero,drop=F])
hist(log10(pen_weights[nonzero]))
sum(pen_weights[nonzero]) == sum(nonzero) # TRUE
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
# (using default BLAS that R is compiled against, can be Accelerate (Mac), OpenBLAS (Ubuntu/Win) or IntelMKL (Ubuntu/Win))
# SOLVE USING BASIC FORMULA FOR n >= p, CF. Liu & Li 2016 ####
system.time(coefs_solve_dense <- as.vector(base::solve(D %*% crossprod(Xw) + # FAST
                                                    diag(rep(lambda, sum(nonzero))),
                                                    D %*% crossprod(Xw,yw),
                                                    tol=2*.Machine$double.eps))) # 13s for n=p=10000, 1.42s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
system.time(coefs_solve_sparse <- as.vector(base::solve(Dsp %*% crossprod(Xwsp) + # FAST
                                                     .sparseDiagonal(sum(nonzero), lambda),
                                                     as.matrix(Dsp %*% crossprod(Xwsp,yw)),
                                                     tol=2*.Machine$double.eps))) # 3.59s for n=p=10000, 3.41s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
beta[nonzero] <- coefs_solve_sparse
cor(beta, coefs_true) # 0.8589209
cor(coefs_solve_dense, coefs_solve_sparse) # 1



system.time(coefs_qrsolve_dense <- as.vector(qr.solve(D %*% crossprod(Xw) +
                                                   diag(rep(lambda, ncol(X))),
                                                 D %*% crossprod(Xw,yw),
                                                 tol=2*.Machine$double.eps))) # 386.13s for n=p=10000, 27.47s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
system.time(coefs_qrsolve_sparse <- as.vector(qr.solve(Dsp %*% crossprod(Xwsp) +
                                                    .sparseDiagonal(sum(nonzero), lambda),
                                                  Dsp %*% crossprod(Xwsp,yw),
                                                  tol=2*.Machine$double.eps))) # 298s for n=p=3000, 28.35s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
beta[nonzero] <- coefs_qrsolve_sparse
cor(beta, coefs_true) # 0.8589209
cor(coefs_qrsolve_sparse, coefs_solve_sparse) # 1
cor(coefs_solve_dense, coefs_qrsolve_dense) # 1
cor(coefs_solve_dense, coefs_qrsolve_sparse) # 1


# eigen's QR methods: slower than default QR solve methods
system.time(coefs_eigenqrsolvepivot0_dense <- as.vector(solve_dense_QR(D %*% crossprod(Xw) +
                                                        diag(rep(lambda, ncol(X))),
                                                     D %*% crossprod(Xw,yw),
                                                      pivot = 0L)$coefficients)) # 40.62s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
system.time(coefs_eigenqrsolveord0_sparse <- as.vector(solve_sparse_QR(Dsp %*% crossprod(Xwsp) +
                                                         .sparseDiagonal(sum(nonzero), lambda),
                                                       as.matrix(Dsp %*% crossprod(Xwsp,yw)),
                                                       ordering = 0L)$coefficients)) # 35s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
beta[nonzero] <- coefs_eigenqrsolveord0_sparse
cor(beta, coefs_true) # 0.8589209
cor(coefs_eigenqrsolveord0_sparse, coefs_solve_sparse) # 1
cor(coefs_eigenqrsolveord0_sparse, coefs_qrsolve_sparse) # 1


# eigen's iterative preconditioned conjugate gradient method: FAST, but least squares conjugate gradient via row augmentation even faster
system.time(coefs_eigenCGsolve_dense <- as.vector(solve_dense_CG(D %*% crossprod(Xw) +
                                                                         diag(rep(lambda, ncol(X))),
                                                                       D %*% crossprod(Xw,yw),
                                                                 x0=rep(0,ncol(X)),
                                                                 tol=0,
                                                                 maxit=0,
                                                                 precond = 0L)$coefficients)) # 29s for p=n=10000, 1.23s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
system.time(coefs_eigenCGsolve_sparse <- as.vector(solve_sparse_CG(Dsp %*% crossprod(Xwsp) +
                                                                         .sparseDiagonal(sum(nonzero), lambda),
                                                                       as.matrix(Dsp %*% crossprod(Xwsp,yw)),
                                                                       x0=as.matrix(start[nonzero]), # rep(0,ncol(X)),
                                                                       tol=0,
                                                                       maxit=0,
                                                                       precond = 1L)$coefficients)) # 0.10s for p=n=10000, 2.66s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
coefs_eigenCGsolve_sparse
# TODO returns NULL - check what's wrong

# this syntax works though :
system.time(coefs_eigenCGsolve_sparse <- as.vector(solve_sparse_CG(as(crossprod(Xwaugsp), "dgCMatrix"),
                                                                   as.matrix(crossprod(Xwaugsp,ywaug)),
                                                                   x0=as.matrix(start[nonzero]), # rep(0,ncol(X)),
                                                                   tol=0,
                                                                   maxit=0,
                                                                   precond = 1L)$coefficients))
coefs_eigenCGsolve_sparse

beta[nonzero] <- coefs_eigenCGsolve_sparse
cor(beta, coefs_true) # 0.8589209
cor(coefs_solve_sparse, coefs_eigenCGsolve_sparse) # 1
cor(coefs_solve_dense, coefs_eigenCGsolve_dense) # 1


# default chol.solve: relatively fast, but not great
system.time(coefs_cholsolve_dense <- as.vector(chol.solve(D %*% crossprod(Xw) +
                                                        diag(rep(lambda, ncol(X))),
                                                      D %*% crossprod(Xw,yw),
                                                      tol=2*.Machine$double.eps))) # 23s for n=p=10000, 2.02s for n=p=5000, 2.16s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
system.time(coefs_cholsolve_sparse <- as.vector(chol.solve(Dsp %*% crossprod(Xwsp) +
                                                         .sparseDiagonal(sum(nonzero), lambda),
                                                       Dsp %*% crossprod(Xwsp,yw),
                                                       tol=2*.Machine$double.eps))) # 5s for n=p=10000, 48.75s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
# Sometimes gives Error in chol.default(crossprod(X)) :
# the leading minor of order 3843 is not positive definite, so overfitted, some eigenvectors near zero

# alt syntax (more stable) :
system.time(coefs_cholsolve_sparse <- as.vector(chol.solve(as(crossprod(Xwaugsp),"dgCMatrix"),
                                                           crossprod(Xwaugsp,ywaug),
                                                           tol=2*.Machine$double.eps))) 
beta[nonzero] <- coefs_cholsolve_sparse
cor(beta, coefs_true) # 0.8948964
cor(coefs_solve_sparse, coefs_true) # 0.1745552
cor(coefs_solve_dense, coefs_cholsolve_dense) # 1
cor(coefs_solve_dense, coefs_cholsolve_sparse) # 1

# eigen's solve_dense_chol & solve_sparse_chol: FAST for sparse only
system.time(coefs_cholsolve_dense <- as.vector(solve_dense_chol(D %*% crossprod(Xw) +
                                                            diag(rep(lambda, ncol(X))),
                                                          pivot=1L)$coefficients)) # 84s for n=p=10000, 7.96s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
system.time(coefs_cholsolve_sparse <- as.vector(solve_sparse_chol(Dsp %*% crossprod(Xwsp) +
                                                             .sparseDiagonal(sum(nonzero), lambda),
                                                           as.matrix(Dsp %*% crossprod(Xwsp,yw)),
                                                           pivot=0L,
                                                           ordering=0L)$coefficients)) # 0.08s for n=p=10000, 18.76s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
# gives bad solution
# more stable fit:
system.time(coefs_cholsolve_sparse <- as.vector(solve_sparse_chol(as(crossprod(Xwaugsp),"dgCMatrix"),
                                                                  as.matrix(crossprod(Xwaugsp,ywaug)),
                                                                  pivot=0L,
                                                                  ordering=0L)$coefficients)) # 0.08s for n=p=10000, 18.76s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
beta[nonzero] <- coefs_cholsolve_sparse
cor(beta, coefs_true) # 0.8948963
cor(coefs_cholsolve_sparse, coefs_true) # 0.1745552
cor(coefs_cholsolve_dense, coefs_cholsolve_sparse) # 1
cor(coefs_cholsolve_dense, coefs_cholsolve_sparse) # 1

system.time(coefs_cholsolve_sparse_lmridge <- lmridge(X, y, lambda,
                                                weights = weights,
                                                lower = lower,
                                                upper = upper,
                                                nonnegative = nonnegative,
                                                solver = "eigen",
                                                method = 7))



# using row augmentation & default Cholesky decomposition -> FAST both for dense & sparse
system.time(coefs_cholsolve_rowaug_dense <- as.vector(chol.solve(Xwaug,
                                                                 ywaug,
                                                                 tol=2*.Machine$double.eps))) # 13s for n=p=10000, 0.90s for p=5000 & n=1000
system.time(coefs_cholsolve_rowaug_sparse <- as.vector(chol.solve(Xwaugsp,
                                                                  ywaug,
                                                           tol=2*.Machine$double.eps))) # 1.6s for n=p=10000, 2.36s for p=5000 & n=1000

beta[nonzero] <- coefs_cholsolve_rowaug_sparse
cor(beta, coefs_true) # 0.8948963
cor(coefs_cholsolve_rowaug_sparse, coefs_true) # 0.019
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
                                                                          x0=as.matrix(start[nonzero]), # as.matrix(rep(0,ncol(Xwaug))),
                                                                          tol=0,
                                                                          maxit=0L,
                                                                          precond=1L)$coefficients)) # 0.19s for p=n=10000, 0.03s for p=5000 & n=1000
# NOTE: sometimes returns NULL
beta[nonzero] <- coefs_eigenLSCG_rowaug_sparse
cor(beta, coefs_true) # 0.8589209
cor(coefs_eigenLSCG_rowaug_sparse, coefs_true) # 1
cor(coefs_solve_dense, coefs_eigenLSCG_rowaug_dense) # 1
cor(coefs_solve_dense, coefs_eigenLSCG_rowaug_sparse) # 1



# 1.2. SOLVE USING OPTIMISED FORMULA FOR p > n, CF. Liu & Li 2016 eqn (28) ####
# -> THIS IS INDEED FASTER WHEN p > n, ALTHOUGH solve_dense_LSCG & solve_sparse_LSCG STILL FASTEST
system.time({DXT <- tcrossprod(D, Xw); coefs_solve_largep_dense <- as.vector(DXT %*% base::solve(Xw %*% DXT +
                                                                             diag(rep(lambda, nrow(X))),
                                                                             yw,
                                                                             tol=2*.Machine$double.eps))}) # 24s for p=n=10000, 0.26s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (28)
system.time({DXT <- tcrossprod(Dsp, Xwsp); coefs_solve_largep_sparse <- as.vector(DXT %*% base::solve(Xwsp %*% DXT +
                                                                                    .sparseDiagonal(nrow(X), lambda),
                                                                                  yw,
                                                                                  tol=2*.Machine$double.eps))}) # 11s for p=n=10000, 0.44s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (28)
# coefs_solve_largep_sparse <- pmax(lower, pmin(coefs_solve_largep_sparse, upper))
beta[nonzero] <- coefs_solve_largep_sparse
cor(beta, coefs_true) # 0.8948963
cor(coefs_solve_dense, coefs_solve_largep_dense) # 1
cor(coefs_solve_dense, coefs_solve_largep_sparse) # 1

system.time({DXT <- tcrossprod(D, Xw); coefs_solveCG_largep_dense <- as.vector(DXT %*% solve_dense_CG(Xw %*% DXT +
                                                                                                   diag(rep(lambda, nrow(X))),
                                                                                                 yw,
                                                                                                 x0=rep(0,nrow(X)),
                                                                                                 tol=0,
                                                                                                 maxit=0,
                                                                                                 precond = 0L)$coefficients)}) # 33s for p=n=10000, 0.28s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (28)
system.time({DXT <- tcrossprod(Dsp, Xwsp); coefs_solveCG_largep_sparse <- as.vector(DXT %*% solve_sparse_CG(Xwsp %*% DXT +
                                                                                                        .sparseDiagonal(nrow(X), lambda),
                                                                                                      yw,
                                                                                                      x0=rep(0,nrow(X)),
                                                                                                      tol=0,
                                                                                                      maxit=0,
                                                                                                      precond = 0L)$coefficients)}) # 18s for p=n=10000, 0.45s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (28)
beta[nonzero] <- coefs_solveCG_largep_sparse
cor(beta, coefs_true) # 0.8948963
cor(coefs_solve_dense, coefs_solveCG_largep_dense) # 1
cor(coefs_solve_dense, coefs_solveCG_largep_sparse) # 1




# 1.3. TESTS WITH ALREADY IMPLEMENTED SOLVERS ####

# PS I think we best make 1 meta lmridge_solve_eigen R function that comprises all the eigen solvers with common options
# (some ignored for particular solvers) and that we would then have an argument solver="eigen" with
# extra option control.solvers with the desired method & a list with the desired options
# TO DO: Check which of the available solvers allow for the outcome y to be a matrix, ie allowing for multichannel data
# (would be good for multitask extension with grouped variable selection, though due to the weights being different
# we would not have the same covariate matrix for all channels, only the sparsity structure would be the same)
# Vectorizing

# CG
system.time(coefs_CG <- lmridge(X=Xsp[,nonzero,drop=F], y, lambdas[nonzero], weights, 
                                lower=rep(0,sum(nonzero)), 
                                upper=rep(Inf,sum(nonzero)), 
                                nonnegative=TRUE,
                                solver="eigen", method = 5)) # CG, 3.0s for n=p=10000
# RETURNS error
beta[nonzero] <- coefs_CG
cor(beta, coefs_true) # 0.8948963

lines(coefs_CG, col="blue", type="h")
cor(coefs_CG, coefs_true) # 0.9990464
plot(coefs_CG, coefs_solve_sparse)


# LSCG: VERY FAST!
system.time(coefs_LSCG <- lmridge_eigen(X=Xsp, y=y, weights=weights, lambdas=rep(1,ncol(X)),
                                        start = as.matrix(start),
                          lower=lower, upper=upper, control.solver=control.solver.gen("eigen", 3))) # 0.18s for n=p=10000

system.time(coefs_LSCG <- lmridge(X=Xsp, y, lambdas, weights, lower, upper, nonnegative,
                                  solver="eigen", method=3)) # 0.17s TODO make sure that start values can be passed OK via arg start (x0 in Rcpp function)
# RETURNS NULL
beta[nonzero] <- coefs_LSCG
cor(beta, coefs_true) #

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
# TO DO: solve by using inbuilt ridge in glmnet with alpha=0, lambda=lambda/nrow(X) & appropriate penalty weights
# instead of using row augmentation & lambda=0 to put in ridge penalty
# now it uses

# model <- bigGlm(x = Xaug, y = yaug,
# alpha = 0, lambda = 0, intercept = FALSE, standardize = FALSE, standardize.response = FALSE,
# lower.limits = c(-Inf, lower), upper.limits = c(Inf, upper),
# family = "gaussian", path = FALSE)
# coefs <- coef(model)[-1]

system.time(coefs_glmnet_unbounded <- lmridge(X=Xsp[,nonzero,drop=F], y, lambdas[nonzero], weights,
                                              lower=rep(-Inf,sum(nonzero)),
                                              upper=rep(Inf,sum(nonzero)),
                                              nonnegative=F,
                                              solver="glmnet")) # 0.56s with p=n=10000 & Intel MKL
beta[nonzero] <- coefs_glmnet_unbounded
cor(beta, coefs_true) # 0.8720732

lines(coefs_glmnet_unbounded, col="orange", type="h")
cor(coefs_glmnet_unbounded, coefs_solve_sparse) # 1


# solvers with box or nonnegativity constraints
# glmnet (bounded): FAST, but LSCG still a bit faster?? Fastest method with explicit box constraints though.
system.time(coefs_glmnet <- lmridge(X=Xsp[,nonzero,drop=F], y,
                                    lambdas[nonzero],
                                    weights,
                                    lower=rep(0,sum(nonzero)),
                                    upper=rep(Inf,sum(nonzero)),
                                    nonnegative=F,
                                    solver="glmnet")) # 0.03s with p=n=10000 & IntelMKL
beta[nonzero] <- coefs_glmnet
cor(beta, coefs_true) # 0.9553911
cor(coefs_glmnet, coefs_true)
lines(coefs_glmnet, col="blue", type="h")

plot(coefs_glmnet, coefs_true)

system.time(coefs_glmnet_dense <- lmridge(X=X, y, lambda, weights,
                                          lower, upper, nonnegative,
                                          solver="glmnet")) # 8.3s with p=n=10000 with MRO/IntelMKL


# osqp: FASTER if  n > p & IntelMKL is installed, otherwise rather slow
system.time(coefs_osqp <- lmridge(X=Xsp[,nonzero,drop=F], y,
                                    lambdas[nonzero],
                                    weights,
                                    lower=rep(0,sum(nonzero)),
                                    upper=rep(Inf,sum(nonzero)),
                                    nonnegative=T,
                                  solver="osqp")) # 0.21s with p=n=10000 with MRO/IntelMKL
beta[nonzero] <- coefs_osqp
cor(beta, coefs_true) # 0.9230191

cor(coefs_osqp, coefs_true)

lines(coefs_osqp, col="green", type="h")
cor(coefs_osqp, coefs_true) # 0.9999903
plot(coefs_osqp, coefs_glmnet) # correct...


# nnls: SLOW
system.time(coefs_nnls <- lmridge(X=Xsp[,nonzero,drop=F], y,
                                  lambdas[nonzero],
                                  weights,
                                  lower=rep(0,sum(nonzero)),
                                  upper=rep(Inf,sum(nonzero)),
                                  nonnegative=T, solver="nnls")) # 1841s for n=p=10000
beta[nonzero] <- coefs_nnls
cor(beta, coefs_true) # 0.9712636

lines(coefs_nnls, col="yellow", type="h")
cor(coefs_nnls, coefs_true) # 0.98


# fastnnls:
system.time(coefs_fastnnls_dense <- fastnnls_dense(Xwaug,
                                                   ywaug)) # 21s for p=n=10000
beta[nonzero] <- coefs_fastnnls_dense
cor(beta, coefs_true) # 0.005004153 # not OK, check where it is going wrong...



# tests with some other solvers
# unbounded ones

# fastglm (here using Eigen's Cholesky LLT solver) - slow as only dense solvers have been implemented
fglm <- fastglm(Xwaug, ywaug, family = gaussian(), method = 2)
coefs_LLT <- coef(fglm)
beta[nonzero] <- coefs_LLT
cor(beta, coefs_true) # 0.8948963



# sgd: currently no box constraints though, works if x is a big.matrix though (bigmemory package)
# also allows start values
library(sgd)
system.time(coefs_sgd <- sgd(x = Xwaug, y = ywaug, model="lm",
                             sgd.control=list(method="ai-sgd",
                                              lr="adagrad",
                                              reltol=1E-10,
                                              start=rep(0,ncol(X))))$coefficients) # 27s with p=n=10000
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
