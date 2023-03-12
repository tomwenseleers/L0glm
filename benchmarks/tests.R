library(Matrix)
library(nnls)
library(microbenchmark)
# devtools::install_github("nk027/sanic")
library(sanic)
library(osqp)
library(devtools)
library(glmnet)
devtools::load_all()

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
n <- 1000
p <- 1000
s <- 0.1 # sparsity as proportion of p that have nonzero coefficients // p*s == #of peaks
k <- round(p*s) # nr of nonzero covariates
sim <- simulate_spike_train(n=n, p=p, k=k, mean_beta = 1000, sd_logbeta = 1, family="poisson")

coefs_true <- sim$beta_true

X <- sim$X
lambda <- 1 # ridge penalty
pen_weights <- rep(1, ncol(X)) # penalization weights
lambdas <- lambda*pen_weights
D <- diag(1/pen_weights) # cf notation Liu & Li 2016
Dsp <- .sparseDiagonal(ncol(X), 1/pen_weights)  # sparse version of D, cf notation Liu & Li 2016
y <- sim$y
weights <- 1/(y+0.1) # observation weights
Xw <- X*sqrt(weights)
Xwaug <- rbind(Xw,diag(sqrt(lambdas)))
Xsp <- Matrix(X, sparse=TRUE)
Xwsp <- Matrix(Xw, sparse=TRUE)
Xwaugsp <- Matrix(Xwaug, sparse=TRUE)
yw <- y*sqrt(weights)
ywaug <- c(yw, rep(0, ncol(X)))
lower <- rep(0, ncol(X))
upper <- rep(Inf, ncol(X))
nonnegative <- TRUE
start <- NULL
verbose <- FALSE


family <- poisson(identity)
converged_set = rep(FALSE, ncol(X))
control.iwls = list(maxit = 1, reltol = 1E-4, thresh = 1E-5, warn = FALSE)
control.solver = NULL
solver = "osqp"

system.time(c <- glm.iwls2(X = X, y = y, weights = weights, family = poisson(identity),
          solver = "osqp",
          lambda = lambdas,
          nonnegative = nonnegative,
          lower = lower,
          upper = upper,
          verbose = verbose))

cor(b$coefficients, coefs_true)

plot_L0glm_benchmark(x = sim$x, y = y, fit = a, a.true = sim$beta_true,
                     main="Ground truth vs L0glm estimates")


system.time(b <- glm.iwls(X = X, y = y, weights = weights, family = poisson(identity),
               lambda = lambdas,
               nonnegative = nonnegative))


plot_L0glm_benchmark(x = sim$x, y = y, fit = b, a.true = sim$a,
                     main="Ground truth vs L0glm estimates")



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


# 1. TESTS WITH SOME DEFAULT SOLVE METHODS TO SOLVE RIDGE REGRESSION, EITHER WITH DIRECT SOLVE OR VIA ROW AUGMENTATION ####
# 1.1 TEST WITH DIFFERENT SOLVE METHODS ####
# (using default BLAS that R is compiled against, can be Accelerate (Mac), OpenBLAS (Ubuntu/Win) or IntelMKL (Ubuntu/Win))
# SOLVE USING BASIC FORMULA FOR n >= p, CF. Liu & Li 2016 ####
system.time(coefs_solve_dense <- as.vector(base::solve(D %*% crossprod(Xw) + # FAST
                                                    diag(rep(lambda, ncol(X))),
                                                    D %*% crossprod(Xw,yw),
                                                    tol=2*.Machine$double.eps))) # 13s for n=p=10000, 1.42s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
system.time(coefs_solve_sparse <- as.vector(base::solve(Dsp %*% crossprod(Xwsp) + # FAST
                                                     .sparseDiagonal(ncol(X), lambda),
                                                     as.matrix(Dsp %*% crossprod(Xwsp,yw)),
                                                     tol=2*.Machine$double.eps))) # 4.7s for n=p=10000, 3.41s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
cor(coefs_solve_sparse, coefs_true) # 0.1572276
cor(coefs_solve_dense, coefs_solve_sparse) # 1



system.time(coefs_qrsolve_dense <- as.vector(qr.solve(D %*% crossprod(Xw) +
                                                   diag(rep(lambda, ncol(X))),
                                                 D %*% crossprod(Xw,yw),
                                                 tol=2*.Machine$double.eps))) # 386.13s for n=p=10000, 27.47s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
system.time(coefs_qrsolve_sparse <- as.vector(qr.solve(Dsp %*% crossprod(Xwsp) +
                                                    .sparseDiagonal(ncol(X), lambda),
                                                  Dsp %*% crossprod(Xwsp,yw),
                                                  tol=2*.Machine$double.eps))) # 298s for n=p=3000, 28.35s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
cor(coefs_solve_sparse, coefs_true) # 0.1651061
cor(coefs_solve_dense, coefs_qrsolve_dense) # 1
cor(coefs_solve_dense, coefs_qrsolve_sparse) # 1

# eigen's QR methods: slower than default QR solve methods
system.time(coefs_eigenqrsolvepivot0_dense <- as.vector(solve_dense_QR(D %*% crossprod(Xw) +
                                                        diag(rep(lambda, ncol(X))),
                                                     D %*% crossprod(Xw,yw),
                                                      pivot = 0L)$coefficients)) # 40.62s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
system.time(coefs_eigenqrsolveord0_sparse <- as.vector(solve_sparse_QR(Dsp %*% crossprod(Xwsp) +
                                                         .sparseDiagonal(ncol(X), lambda),
                                                       as.matrix(Dsp %*% crossprod(Xwsp,yw)),
                                                       ordering = 0L)$coefficients)) # 35s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
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
                                                                         .sparseDiagonal(ncol(X), lambda),
                                                                       as.matrix(Dsp %*% crossprod(Xwsp,yw)),
                                                                       x0=rep(0,ncol(X)),
                                                                       tol=0,
                                                                       maxit=0,
                                                                       precond = 0L)$coefficients)) # 1.58s for p=n=10000, 2.66s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
cor(coefs_solve_sparse, coefs_eigenCGsolve_sparse) # 1
cor(coefs_solve_dense, coefs_eigenCGsolve_dense) # 1


# default chol.solve: relatively fast, but not great
system.time(coefs_cholsolve_dense <- as.vector(chol.solve(D %*% crossprod(Xw) +
                                                        diag(rep(lambda, ncol(X))),
                                                      D %*% crossprod(Xw,yw),
                                                      tol=2*.Machine$double.eps))) # 23s for n=p=10000, 2.02s for n=p=5000, 2.16s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
system.time(coefs_cholsolve_sparse <- as.vector(chol.solve(Dsp %*% crossprod(Xwsp) +
                                                         .sparseDiagonal(ncol(X), lambda),
                                                       Dsp %*% crossprod(Xwsp,yw),
                                                       tol=2*.Machine$double.eps))) # 9.69s for n=p=10000, 48.75s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
cor(coefs_solve_sparse, coefs_true) # 0.1745552
cor(coefs_solve_dense, coefs_cholsolve_dense) # 1
cor(coefs_solve_dense, coefs_cholsolve_sparse) # 1

# eigen's solve_dense_chol & solve_sparse_chol: FAST for sparse only
system.time(coefs_cholsolve_dense <- as.vector(solve_dense_chol(D %*% crossprod(Xw) +
                                                            diag(rep(lambda, ncol(X))),
                                                          D %*% crossprod(Xw,yw),
                                                          pivot=1L)$coefficients)) # 84s for n=p=10000, 7.96s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
system.time(coefs_cholsolve_sparse <- as.vector(solve_sparse_chol(Dsp %*% crossprod(Xwsp) +
                                                             .sparseDiagonal(ncol(X), lambda),
                                                           as.matrix(Dsp %*% crossprod(Xwsp,yw)),
                                                           pivot=0L,
                                                           ordering=0L)$coefficients)) # 1.7s for n=p=10000, 18.76s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (10)
cor(coefs_cholsolve_sparse, coefs_true) # 0.1745552
cor(coefs_cholsolve_dense, coefs_cholsolve_sparse) # 1
cor(coefs_cholsolve_dense, coefs_cholsolve_sparse) # 1


# using row augmentation & default Cholesky decomposition -> FAST both for dense & sparse X
system.time(coefs_cholsolve_rowaug_dense <- as.vector(chol.solve(Xwaug,
                                                                 ywaug,
                                                                 tol=2*.Machine$double.eps))) # 13s for n=p=10000, 0.90s for p=5000 & n=1000
system.time(coefs_cholsolve_rowaug_sparse <- as.vector(chol.solve(Xwaugsp,
                                                                  ywaug,
                                                           tol=2*.Machine$double.eps))) # 3s for n=p=10000, 2.36s for p=5000 & n=1000
cor(coefs_solve_dense, coefs_cholsolve_rowaug_dense) # 1
cor(coefs_solve_dense, coefs_cholsolve_rowaug_sparse) # 1


# using row augmentation & eigen's Cholesky decomposition -> crashes for the moment (sanic version also crashes), check why...
system.time(coefs_eigencholsolve_rowaug_dense <- as.vector(solve_dense_chol(Xwaug,
                                                                       as.matrix(ywaug),
                                                                       pivot=1L)$coefficients)) # CRASHES XXs for p=5000 & n=1000
system.time(coefs_eigencholsolve_rowaug_sparse <- as.vector(solve_sparse_chol(Xwaugsp,
                                                                  ywaug,
                                                                  pivot=0L,
                                                                  ordering=0L)$coefficients)) # CRASHES XXs for p=5000 & n=1000
cor(coefs_solve_dense, coefs_eigencholsolve_rowaug_dense) # 1
cor(coefs_solve_dense, coefs_eigencholsolve_rowaug_sparse) # 1



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
                                                                          precond=1L)$coefficients)) # 1.4s for p=n=10000, 0.03s for p=5000 & n=1000
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
                                                                                  tol=2*.Machine$double.eps))}) # 15s for p=n=10000, 0.44s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (28)
# coefs_solve_largep_sparse <- pmax(lower, pmin(coefs_solve_largep_sparse, upper))
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
                                                                                                      precond = 0L)$coefficients)}) # 21s for p=n=10000, 0.45s for p=5000 & n=1000, cf. Liu & Li 2016 eqn (28)
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

system.time(coefs_QR <- lmridge(X=Xsp, y, lambda, weights, lower, upper, nonnegative, solver="eigen", method = 6)) # 19.97s for p=5000 & n=1000
lines(coefs_QR, col="blue", type="h")
cor(coefs_QR, coefs_solve_dense) # 0.9990464
plot(coefs_QR, coefs_solve_sparse)

# LSCG: VERY FAST!
system.time(coefs_LSCG <- lmridge(X=Xsp, y, lambda, weights, lower, upper, nonnegative, solver="LSCG")) # 1.4s for p=n=10000, 0.06s for p=5000 & n=1000
lines(coefs_LSCG, col="blue", type="h")
cor(coefs_LSCG, coefs_solve_dense) # 0.9990464
plot(coefs_LSCG, coefs_solve_sparse)

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

system.time(coefs_glmnet_unbounded <- lmridge(X=Xsp, y, lambda, weights,
                                              lower=rep(-Inf,p), upper=rep(Inf,p), nonnegative=F,
                                              solver="glmnet")) # 3.6s with p=n=10000, 0.1s with p=5000 & n=1000 with MRO/IntelMKL
lines(coefs_glmnet_unbounded, col="orange", type="h")
cor(coefs_glmnet_unbounded, coefs_solve_dense) # 1


# solvers with box or nonnegativity constraints
# glmnet (bounded): FAST, but LSCG still a bit faster?? Fastest method with explicit box constraints though.
system.time(coefs_glmnet <- lmridge(X=Xsp, y, lambda = lambdas, weights,
                                    lower, upper, nonnegative,
                                    solver="glmnet")) # 0.11s with p=n=10000, 0.06s with p=5000 & n=1000 with MRO/IntelMKL


cor(coefs_glmnet, coefs_true)
lines(coefs_glmnet, col="blue", type="h")

plot(coefs_glmnet, coefs_glmnet2)

system.time(coefs_glmnet_dense <- lmridge(X=X, y, lambda, weights,
                                          lower, upper, nonnegative,
                                          solver="glmnet")) # 6.8s with p=n=10000, 1.26s with p=5000 & n=1000 with MRO/IntelMKL


# osqp: FASTER if  n > p & IntelMKL is installed, otherwise rather slow
system.time(coefs_osqp <- lmridge(X=Xsp, y, lambda=0.1, weights, lower, upper, nonnegative, solver="osqp")) # 3.6s with p=n=10000, 48.57s with p=5000 & n=1000, 3.7s for n=p=10000 with MRO/IntelMKL
# 6.78s when using Xaug, 7.09s when using Xw * diag(lambda), coefs are the same
cor(coefs_osqp, coefs_osqp0)

lines(coefs_osqp, col="green", type="h")
cor(coefs_osqp, coefs_true) # 0.9999903
plot(coefs_osqp, coefs_glmnet) # correct...
# TO DO check how to optimise for p > n

# nnls: SLOW
system.time(coefs_nnls <- lmridge(X=Xsp, y, lambda, weights, lower, upper, nonnegative, solver="nnls")) # 1841s for n=p=10000
lines(coefs_nnls, col="yellow", type="h")
cor(coefs_nnls, coefs_true) # 0.98

# fastnnls:
system.time(coefs_fastnnls_dense <- fastnnls(Xwaug, ywaug)) # 24s for p=n=10000


# tests with some other solvers
# unbounded ones
# sgd: currently no box constraints though, works if x is a big.matrix though (bigmemory package)
# also allows start values
library(sgd)
system.time(coefs_sgd <- sgd(x = Xwaug, y = ywaug, model="lm",
                             sgd.control=list(method="ai-sgd",
                                              lr="adagrad",
                                              reltol=1E-10,
                                              start=rep(0,ncol(X))))$coefficients) # 27s with p=n=10000, 3.3s with p=5000 & n=1000
cor(coefs_sgd, coefs_solve_dense) # 0.9953188


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
