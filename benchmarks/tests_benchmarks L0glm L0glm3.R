devtools::load_all()

# L0glm solver benchmarks ####
library(remotes)
# remotes::install_github("tomwenseleers/L0glm/L0glm")
library(L0glm)
library(microbenchmark)
library(Matrix)

# Simulate data
set.seed(1)
n <- 10000
p <- 10000
s <- 0.1 # sparsity as proportion of p that have nonzero coefficients // p*s == #of peaks
k <- round(p*s) # nr of nonzero covariates
sim <- simulate_spike_train(n=n, p=p, k=k, 
                            mean_beta = 1000, 
                            sd_logbeta = 1, 
                            family="poisson")
# sim <- simulate_spike_train()
coefs_true <- as.vector(sim$beta_true)
hist(coefs_true)
sum(coefs_true!=0) # 1000

Xsp <- sim$X
X <- as.matrix(Xsp)
y <- as.vector(sim$y)
lambda <- 1
pen_weights <- rep(1, ncol(X))
lambdas <- lambda*pen_weights
lower <- rep(0, ncol(X))
upper <- rep(Inf, ncol(X))
nonnegative <- TRUE
start <- NULL
verbose <- TRUE
weights <- 1/(y+0.1) # observation weights
Xspwaug = rbind(Xsp*sqrt(weights), Diagonal(x=sqrt(lambdas)))
ywaug = c(y*sqrt(weights), rep(0, ncol(Xsp)))

# colSums(Xsp) # 1st column for intercept is set at 1, rest sum to norm(Xsp[,100,drop=F], type="1")=sqrt(2*pi)/sqrt(1/w^2) with peak width w=2.5 = 6.266548 (max is 1)
# X.n <- sqrt(colSums(Xsp^2)) # L2 norm normalisation, scaling factors for columns of X = norm(Xsp[,2,drop=F], type="2") = sqrt(sum(Xsp[,2,drop=F]^2)) = 2.105026, PS much faster than apply(Xsp, 2, base::norm, type = "2")
X.n <- colSums(Xspwaug^2) # squared L2 norm normalisation, scaling factors for columns of X
# X.n <- colSums(abs(Xsp)) # L1 norm normalisation
# X.n <- colMax(abs(Xsp)) # normalisation by absolute value of column max
# X.n <- rep(1, ncol(Xsp))
# microbenchmark(Xspwaug <- sweep(Xspwaug, 2, X.n, "/"))
microbenchmark(Xspwaug <- Xspwaug / X.n) # 2.4 ms, faster than Diagonal(x=1/X.n) * Xspwaug) & much faster than Xspwaug <- sweep(Xspwaug, 2, X.n, "/")

thresh = 1E-10 # coefficients below this threshold will be set to 0

# couple of possible coefficient initialisations

# ridge regression
microbenchmark(beta_init_ridge <- solve_sparse_LSCG(A = Xspwaug,
                                                    b = ywaug,
                                                    x0 = ywaug[1:ncol(Xspwaug)], # or even better: use beta_init below in section diagonal/jacobi preconditioner
                                                    tol = 1E-5,
                                                    maxit = 0L,
                                                    precond = 1L)$coefficients, unit="s") 
# 0.000298331s, faster than Jacobi initialisation below even
cor(as.vector(beta_init_ridge)*X.n, coefs_true)
hist(as.vector(beta_init_ridge)*X.n) # coefficients are bounded between 0 and 1/lambda
range(as.vector(beta_init_ridge)*X.n) # 0.0000000 0.9999891
length(intersect(which(coefs_true>0), which(as.vector(beta_init_ridge*X.n)>thresh))) # 1000, all true coefficients that are >0 come out in inferred set
sum(as.vector(beta_init_ridge*X.n)>thresh) # 8683

# Jacobi/diagonal preconditioner initialisation - near identical to LSCG ridge regression coefs above & slightly slower
microbenchmark(beta_init_jacobi <- crossprod(Xspwaug, ywaug) * colSums(Xspwaug^2), unit="s") 
# 0.0014s, colSums(Xspwaug^2) = diag(crossprod(Xspwaug)))

max(abs(beta_init_jacobi-beta_init_ridge)) # 9.65894e-14, near identical to ridge regression coefficients
cor(as.vector(beta_init_jacobi)/X.n, coefs_true)
plot((as.vector(beta_init_jacobi)/X.n), coefs_true)
hist(as.vector(beta_init_jacobi)/X.n) # coefficients are bounded between 0 and 1/lambda
range(as.vector(beta_init_jacobi)/X.n, na.rm=T) # 0.0000000 0.9999891
length(intersect(which(coefs_true>0), which(as.vector(beta_init_jacobi/X.n)>thresh))) # 1000, all true coefficients that are >0 come out in inferred set
sum(as.vector(beta_init_jacobi/X.n)>thresh) # 8683


# osqp nonnegative weighted least square initialisation
system.time(beta_init_nnwls_osqp <- L0glm.fit(X=Xsp, y=y, weights = weights, family = gaussian(identity),
                                              start = as.vector(beta_init_ridge),
                                              solver = "osqp", lambda = 0, nonnegative = TRUE,
                                              control.iwls = do.call("control.iwls.gen", list(maxit = 0)),
                                              control.l0 = do.call("control.l0.gen", list(maxit = 0, delta = 1E-10, gamma = 2)),
                                              control.solver = do.call("control.solver.gen", 
                                                                       list(solver = "osqp", rowaug = FALSE)))$coefficients 
) # 2.16s
hist(log10(beta_init_nnwls_osqp/X.n))
length(intersect(which(coefs_true>0), which(as.vector(beta_init_nnwls_osqp/X.n)>thresh))) # 927, not all true coefficients that are >0 come out in inferred set
sum(as.vector(beta_init_nnwls_osqp)/X.n>thresh) # 2300
beta_init_nnwls_osqp[beta_init_nnwls_osqp/X.n < thresh] = 0
cor(beta_init_nnwls_osqp/X.n, coefs_true) # 0.9573048
plot(beta_init_nnwls_osqp/X.n, coefs_true)
mean(abs(beta_init_nnwls_osqp/X.n-coefs_true)) # 54.40869
sum(coefs_true!=0) # 1000
sum(beta_init_nnwls_osqp!=0) # 2300
length(intersect(which(coefs_true>0), which(as.vector(beta_init_nnwls_osqp/X.n)>0))) # 927


system.time(beta_init_nnwls_osqp2 <- lmridge.osqp(Xw = Xsp, yw = y, weights = weights, lambdas = rep(0, ncol(X)), 
                                                  start = as.vector(beta_init_ridge), # start=y or beta_init_ridge both OK
                                                  lower = lower, upper = upper, control.solver = do.call("control.solver.gen", 
                                                                                                         list(solver = "osqp", rowaug = FALSE)))) # 0.94s
# 0.78s, note: figure out why call above is so much slower
hist(log10(beta_init_nnwls_osqp2/X.n))
length(intersect(which(coefs_true>0), which(as.vector(beta_init_nnwls_osqp2/X.n)>thresh))) # 993, nearly all true coefficients that are >0 come out in inferred set
sum(as.vector(beta_init_nnwls_osqp2)/X.n>thresh) # 9361
beta_init_nnwls_osqp2[(beta_init_nnwls_osqp2/X.n) < thresh] = 0
cor(beta_init_nnwls_osqp2/X.n, coefs_true) # 0.7966744
plot(beta_init_nnwls_osqp2/X.n, coefs_true)
mean(abs(beta_init_nnwls_osqp2/X.n-coefs_true)) # 166.6863
sum(coefs_true!=0) # 1000
sum(beta_init_nnwls_osqp2!=0) # 9361
hist(log10(beta_init_nnwls_osqp2/X.n))
length(intersect(which(coefs_true>0), which(as.vector(beta_init_nnwls_osqp/X.n)>0))) # 850


# L0glm with OSQP solver
system.time(L0glm_osqp <- L0glm.fit(X=Xsp, y=y, 
                                    start=as.vector(beta_init_ridge), # beta_init_nnwls_osqp2, OR y,
                                    family = poisson(identity),
                                    solver = "osqp", lambda = 1, nonnegative = TRUE,
                                    control.iwls = do.call("control.iwls.gen", list(maxit = 1)),
                                    control.l0 = do.call("control.l0.gen", list(maxit = 1000, delta = 1E-10, gamma = 2)),
                                    control.solver = do.call("control.solver.gen", 
                                                             list(solver = "osqp", rowaug = FALSE))
)) # 14s, # 1.94s for n=p=10000
hist(log10(L0glm_osqp$coefficients/X.n))
length(intersect(which(coefs_true>0), which(as.vector(L0glm_osqp$coefficients/X.n)>1E-10))) # 850, not all true coefficients that are >0 come out in inferred set

cor(L0glm_osqp$coefficients/X.n, coefs_true) # 0.8279579, # 0.9563597
plot(L0glm_osqp$coefficients/X.n, coefs_true)
mean(abs(L0glm_osqp$coefficients/X.n-coefs_true)) # 113.027 # 43.75156
sum(coefs_true!=0) # 1000
sum(L0glm_osqp$coefficients!=0) # 1568
length(intersect(which(coefs_true>0), which(as.vector(L0glm_osqp$coefficients/X.n)>0))) # 850


# L0glm with Glmnet solver
system.time(L0glm_glmnet <- L0glm.fit(X=Xsp, y=y,
                                      start=as.vector(beta_init_ridge), # y, # as.vector(beta_init_ridge), # beta_init_nnwls_osqp2,
                                      family = poisson(identity),
                                      solver = "glmnet", lambda = 1, nonnegative = TRUE,
                                      control.iwls = do.call("control.iwls.gen", list(maxit = 1)),
                                      control.l0 = do.call("control.l0.gen", list(maxit = 1000, delta = 1E-10, gamma = 2)))
) # 7.3s, # 3.91s for n=p=10000
cor(L0glm_glmnet$coefficients/X.n, coefs_true) # 0.9356395
plot(L0glm_glmnet$coefficients/X.n, coefs_true)
mean(abs(L0glm_glmnet$coefficients/X.n-coefs_true)) # 44.82373
sum(coefs_true!=0) # 1000
sum(L0glm_glmnet$coefficients!=0) # 1078
length(intersect(which(coefs_true>0), which(as.vector(L0glm_glmnet$coefficients/X.n)>0))) # 794


# L0glm with Eigen Cholesky solver
system.time(L0glm_eigenCHOL <- L0glm.fit(X=Xsp, y=y, 
                                         start=as.vector(beta_init_ridge), # y, # as.vector(beta_init_ridge), # beta_init_nnwls_osqp2,
                                         family = poisson(identity),
                                         solver = "eigen", method = 7, 
                                         lambda = 1, nonnegative = TRUE,
                                         control.iwls = do.call("control.iwls.gen", list(maxit = 1)),
                                         control.l0 = do.call("control.l0.gen", list(maxit = 200, delta = 1E-10, gamma = 2)),
                                         control.solver = do.call("control.solver.gen", 
                                                                  list(solver = "eigen", method = 7,
                                                                       rowaug = FALSE))
)) # 0.67s for n=p=10000
cor(L0glm_eigenCHOL$coefficients/X.n, coefs_true) # 0.9484723
plot(L0glm_eigenCHOL$coefficients/X.n, coefs_true)
mean(abs(L0glm_eigenCHOL$coefficients/X.n-coefs_true)) # 39.49221
sum(coefs_true!=0) # 1000
sum(L0glm_eigenCHOL$coefficients!=0) # 1072
length(intersect(which(coefs_true>0), which(as.vector(L0glm_eigenCHOL$coefficients/X.n)>0))) # 803

# L0glm with Eigen CG solver
system.time(L0glm_eigenCG <- L0glm.fit(X=Xsp, y=y, 
                                       start=as.vector(beta_init_ridge), # y, # as.vector(beta_init_ridge), # beta_init_nnwls_osqp2,
                                       family = poisson(identity),
                                       solver = "eigen", method = 5, 
                                       lambda = 1, nonnegative = TRUE,
                                       control.iwls = do.call("control.iwls.gen", list(maxit = 1)),
                                       control.l0 = do.call("control.l0.gen", list(maxit = 200, delta = 1E-10, gamma = 2)),
                                       control.solver = do.call("control.solver.gen", 
                                                                list(solver = "eigen", method = 5, rowaug = FALSE))
)) # 2.29s for n=p=10000
cor(L0glm_eigenCG$coefficients/X.n, coefs_true) # 0.9484723
plot(L0glm_eigenCG$coefficients/X.n, coefs_true)
mean(abs(L0glm_eigenCG$coefficients/X.n-coefs_true)) # 39.49221
sum(coefs_true!=0) # 1000
sum(L0glm_eigenCG$coefficients!=0) # 1057
length(intersect(which(coefs_true>0), which(as.vector(L0glm_eigenCG$coefficients/X.n)>0))) # 803

# L0glm with QR solver
system.time(L0glm_eigenQR <- L0glm.fit(X=Xsp, y=y, 
                                       start=as.vector(beta_init_ridge), # y, # as.vector(beta_init_ridge), # beta_init_nnwls_osqp2,
                                       family = poisson(identity),
                                       solver = "eigen", method = 2, 
                                       lambda = 1, nonnegative = TRUE,
                                       control.iwls = do.call("control.iwls.gen", list(maxit = 1)),
                                       control.l0 = do.call("control.l0.gen", list(maxit = 100, delta = 1E-10, gamma = 2)),
                                       control.solver = do.call("control.solver.gen", 
                                                                list(solver = "eigen", method = 2, rowaug = TRUE))
)) # 589s for n=p=10000
cor(L0glm_eigenQR$coefficients/X.n, coefs_true) # 0.9441687
plot(L0glm_eigenQR$coefficients/X.n, coefs_true)
mean(abs(L0glm_eigenQR$coefficients/X.n-coefs_true)) # 39.87579
sum(coefs_true!=0) # 1000
sum(L0glm_eigenQR$coefficients!=0) # 1072
length(intersect(which(coefs_true>0), which(as.vector(L0glm_eigenQR$coefficients/X.n)>0))) # 803


# LSCG DOESN'T WORK YET
system.time(L0glm_eigenLSCG <- L0glm.fit(X=Xsp, y=y, 
                                         start=beta_init_nnwls_osqp2, # as.vector(beta_init_ridge), # y, # as.vector(beta_init_ridge), # beta_init_nnwls_osqp2,
                                         family = poisson(identity),
                                         solver = "eigen", method = 4, 
                                         lambda = 1, nonnegative = TRUE,
                                         control.iwls = do.call("control.iwls.gen", list(maxit = 1)),
                                         control.l0 = do.call("control.l0.gen", list(maxit = 100, delta = 1E-10, gamma = 2)),
                                         control.solver = do.call("control.solver.gen", 
                                                                  list(solver = "eigen", method = 4, rowaug = TRUE))
)) # XXXs for n=p=10000
cor(L0glm_eigenLSCG$coefficients/X.n, coefs_true) # 
plot(L0glm_eigenLSCG$coefficients/X.n, coefs_true)
mean(abs(L0glm_eigenLSCG$coefficients/X.n-coefs_true)) # 
sum(coefs_true!=0) # 
sum(L0glm_eigenLSCG$coefficients!=0) # 

# error in 
lmridge(X[good,beta != 0, drop = F], y = y[good], lambdas = lambda[beta != 0], weights = W[good],
        pen.weights = pen.weights[good], lower = lower[beta != 0], upper = upper[beta != 0], 
        nonnegative = nonnegative, solver = "eigen", method = 4, start = betaold[beta != 0], 
        control.solver = control.solver, verbose = verbose)