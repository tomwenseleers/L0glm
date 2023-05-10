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
weights <- 1/(as.vector(y)+0.1) # observation weights
Xspwaug = rbind(Diagonal(x=sqrt(weights)) %*% Xsp, Diagonal(x=sqrt(lambdas)))
ywaug = c(y*sqrt(weights), rep(0, ncol(Xsp)))

# colSums(Xsp) # 1st column for intercept is set at 1, rest sum to norm(Xsp[,100,drop=F], type="1")=sqrt(2*pi)/sqrt(1/w^2) with peak width w=2.5 = 6.266548 (max is 1)
X.n.L2 <- sqrt(colSums(Xsp^2)) # L2 norm normalisation, scaling factors for columns of X = norm(Xsp[,2,drop=F], type="2") = sqrt(sum(Xsp[,2,drop=F]^2)) = 2.105026, PS much faster than apply(Xsp, 2, base::norm, type = "2")
X.n.sqL2 <- colSums(Xspwaug^2) # squared L2 norm normalisation, scaling factors for columns of X
# X.n <- colSums(abs(Xsp)) # L1 norm normalisation
# X.n <- colMax(abs(Xsp)) # normalisation by absolute value of column max
# X.n <- rep(1, ncol(Xsp))
# microbenchmark(Xspwaug <- sweep(Xspwaug, 2, X.n, "/"))
# microbenchmark(Xspwaug <- Xspwaug / X.n) # 2.4 ms, faster than Diagonal(x=1/X.n) * Xspwaug) & much faster than Xspwaug <- sweep(Xspwaug, 2, X.n, "/")

thresh = 1E-10 # coefficients below this threshold will be set to 0

# couple of possible coefficient initialisations ####

# ridge regression initialisation ####
microbenchmark(beta_init_ridge <- solve_sparse_LSCG(A = Xspwaug / X.n.L2,
                                                    b = ywaug,
                                                    x0 = ywaug[1:ncol(Xspwaug)], # or even better: use beta_init below in section diagonal/jacobi preconditioner
                                                    tol = 1E-3,
                                                    maxit = 0L,
                                                    precond = 1L)$coefficients / X.n.L2, unit="s") 
# 0.02s
cor(as.vector(beta_init_ridge), coefs_true)
hist(as.vector(beta_init_ridge)) #
length(intersect(which(coefs_true>0), which(as.vector(beta_init_ridge)>thresh))) # 999, nearly all true coefficients that are >0 come out in inferred set
sum(as.vector(beta_init_ridge)>thresh) # 6662

# Jacobi/diagonal preconditioner initialisation ####
microbenchmark(beta_init_jacobi <- as.vector(crossprod(Xspwaug, ywaug) / colSums(Xspwaug^2)), unit="s") 
# 0.0012s, colSums(Xspwaug^2) = diag(crossprod(Xspwaug)))

cor(beta_init_jacobi, beta_init_ridge) # 0.96, highly correlated to true ridge solution
plot(beta_init_jacobi, beta_init_ridge)

cor(beta_init_jacobi, coefs_true)
plot(beta_init_jacobi, coefs_true)
length(intersect(which(coefs_true>0), which(beta_init_jacobi>thresh))) # 1000, all true coefficients that are >0 come out in inferred set
sum(beta_init_jacobi>thresh) # 9894


# osqp nonnegative weighted least square initialisation ####
system.time(beta_init_nnwls_osqp <- L0glm.fit(X=Xsp / X.n.L2, y=y, weights = weights, family = gaussian(identity),
                                              start = as.vector(beta_init_ridge),
                                              solver = "osqp", lambda = 0, nonnegative = TRUE,
                                              control.iwls = do.call("control.iwls.gen", list(maxit = 0)),
                                              control.l0 = do.call("control.l0.gen", list(maxit = 0, delta = 1E-10, gamma = 2)),
                                              control.solver = do.call("control.solver.gen", 
                                                                       list(solver = "osqp", rowaug = FALSE)))$coefficients / X.n.L2
) # 2.6s
hist(log10(beta_init_nnwls_osqp))
length(intersect(which(coefs_true>0), which(as.vector(beta_init_nnwls_osqp)>thresh))) # 923, not all true coefficients that are >0 come out in inferred set
sum(as.vector(beta_init_nnwls_osqp)>thresh) # 2047
beta_init_nnwls_osqp[beta_init_nnwls_osqp < thresh] = 0
cor(beta_init_nnwls_osqp, coefs_true) # 0.9629725
plot(beta_init_nnwls_osqp, coefs_true)
mean(abs(beta_init_nnwls_osqp-coefs_true)) # 51.99452
sum(coefs_true!=0) # 1000
sum(beta_init_nnwls_osqp!=0) # 2047
length(intersect(which(coefs_true>0), which(as.vector(beta_init_nnwls_osqp)>0))) # 923


system.time(beta_init_nnwls_osqp2 <- lmridge.osqp(Xw = Xsp / X.n.L2, yw = y, weights = weights, lambdas = rep(0, ncol(X)), 
                                                  start = as.vector(beta_init_ridge), # start=y or beta_init_ridge both OK
                                                  lower = lower, upper = upper, control.solver = do.call("control.solver.gen", 
                                                                                                         list(solver = "osqp", rowaug = FALSE))) / X.n.L2 ) # 0.94s
# 0.78s, note: figure out why call above is so much slower
# CHECK WHY SOLUTION BELOW IS SO BAD
hist(log10(beta_init_nnwls_osqp2))
length(intersect(which(coefs_true>0), which(as.vector(beta_init_nnwls_osqp2)>thresh))) # 990, nearly all true coefficients that are >0 come out in inferred set
sum(as.vector(beta_init_nnwls_osqp2)>thresh) # 9312
beta_init_nnwls_osqp2[(beta_init_nnwls_osqp2) < thresh] = 0
cor(beta_init_nnwls_osqp2, coefs_true) # 0.7931252
plot(beta_init_nnwls_osqp2, coefs_true)
mean(abs(beta_init_nnwls_osqp2-coefs_true)) # 167.8864
sum(coefs_true!=0) # 1000
sum(beta_init_nnwls_osqp2!=0) # 9312
hist(log10(beta_init_nnwls_osqp2))
length(intersect(which(coefs_true>0), which(as.vector(beta_init_nnwls_osqp)>0))) # 923


# L0glm with OSQP solver ####
system.time(L0glm_osqp <- L0glm.fit(X=Xsp / X.n.L2, y=y, 
                                    # start=as.vector(beta_init_ridge), # beta_init_nnwls_osqp2, OR y,
                                    family = poisson(identity),
                                    solver = "osqp", lambda = 1, nonnegative = TRUE,
                                    control.iwls = do.call("control.iwls.gen", list(maxit = 1)),
                                    control.l0 = do.call("control.l0.gen", 
                                                         list(maxit = 100, delta = 1E-10, gamma = 2)),
                                    control.solver = do.call("control.solver.gen", 
                                                             list(solver = "osqp", rowaug = FALSE)))$coefficients / X.n.L2 
            ) # 1.62s for n=p=10000
hist(log10(L0glm_osqp))
length(intersect(which(coefs_true>0), which(as.vector(L0glm_osqp)>thresh))) # 833, all true coefficients that are >0 come out in inferred set

cor(L0glm_osqp, coefs_true) # 0.9615766
plot(L0glm_osqp, coefs_true)
mean(abs(L0glm_osqp-coefs_true)) # 42.5446
sum(coefs_true!=0) # 1000
sum(L0glm_osqp!=0) # 1084
length(intersect(which(coefs_true>0), which(as.vector(L0glm_osqp)>0))) # 833


# L0glm with Glmnet solver ####
system.time(L0glm_glmnet <- L0glm.fit(X=Xsp, y=y,
                                      start=as.vector(beta_init_ridge), # y, # as.vector(beta_init_ridge), # beta_init_nnwls_osqp2,
                                      family = poisson(identity),
                                      solver = "glmnet", lambda = 1, nonnegative = TRUE,
                                      control.iwls = do.call("control.iwls.gen", list(maxit = 1)),
                                      control.l0 = do.call("control.l0.gen", list(maxit = 1000, delta = 1E-10, gamma = 2)))$coefficients
) # 8.4s for n=p=10000
cor(L0glm_glmnet, coefs_true) # 0.9547656
plot(L0glm_glmnet, coefs_true)
mean(abs(L0glm_glmnet-coefs_true)) # 43.6681
sum(coefs_true!=0) # 1000
sum(L0glm_glmnet!=0) # 1071
length(intersect(which(coefs_true>0), which(as.vector(L0glm_glmnet)>0))) # 834


# L0glm with Eigen Cholesky solver ####
system.time(L0glm_eigenCHOL <- L0glm.fit(X=Xsp / X.n.L2, y=y, 
                                         start=as.vector(beta_init_ridge), # y, # as.vector(beta_init_ridge), # beta_init_nnwls_osqp2,
                                         family = poisson(identity),
                                         solver = "eigen", method = 7, 
                                         lambda = 1, nonnegative = TRUE,
                                         control.iwls = do.call("control.iwls.gen", list(maxit = 1)),
                                         control.l0 = do.call("control.l0.gen", list(maxit = 200, delta = 1E-10, gamma = 2)),
                                         control.solver = do.call("control.solver.gen", 
                                                                  list(solver = "eigen", method = 7,
                                                                       rowaug = FALSE)))$coefficients / X.n.L2
) # 0.33s for n=p=10000
cor(L0glm_eigenCHOL, coefs_true) # 0.9313193
plot(L0glm_eigenCHOL, coefs_true)
mean(abs(L0glm_eigenCHOL-coefs_true)) # 73.16189
sum(coefs_true!=0) # 1000
sum(L0glm_eigenCHOL!=0) # 1000
length(intersect(which(coefs_true>0), which(as.vector(L0glm_eigenCHOL)>0))) # 640

# L0glm with Eigen CG solver ####
system.time(L0glm_eigenCG <- L0glm.fit(X=Xsp / X.n.sqL2, y=y, 
                                       start=as.vector(beta_init_ridge), # y, # as.vector(beta_init_ridge), # beta_init_nnwls_osqp2,
                                       family = poisson(identity),
                                       solver = "eigen", method = 5, 
                                       lambda = 1, nonnegative = TRUE,
                                       control.iwls = do.call("control.iwls.gen", list(maxit = 1)),
                                       control.l0 = do.call("control.l0.gen", list(maxit = 200, delta = 1E-10, gamma = 2)),
                                       control.solver = do.call("control.solver.gen", 
                                                                list(solver = "eigen", method = 5, rowaug = FALSE)))$coefficients / X.n.sqL2
) # 0.83s for n=p=10000
cor(L0glm_eigenCG, coefs_true) # 0.9410897
plot(L0glm_eigenCG, coefs_true)
mean(abs(L0glm_eigenCG-coefs_true)) # 78.45663
sum(coefs_true!=0) # 1000
sum(L0glm_eigenCG!=0) # 1026
length(intersect(which(coefs_true>0), which(as.vector(L0glm_eigenCG)>0))) # 644

# L0glm with QR solver ####
system.time(L0glm_eigenQR <- L0glm.fit(X=Xsp, y=y, 
                                       start=as.vector(beta_init_ridge), # y, # as.vector(beta_init_ridge), # beta_init_nnwls_osqp2,
                                       family = poisson(identity),
                                       solver = "eigen", method = 2, 
                                       lambda = 1, nonnegative = TRUE,
                                       control.iwls = do.call("control.iwls.gen", list(maxit = 1)),
                                       control.l0 = do.call("control.l0.gen", list(maxit = 100, delta = 1E-10, gamma = 2)),
                                       control.solver = do.call("control.solver.gen", 
                                                                list(solver = "eigen", method = 2, rowaug = TRUE)))$coefficients
) # 253s for n=p=10000
cor(L0glm_eigenQR, coefs_true) # 0.9444676
plot(L0glm_eigenQR, coefs_true)
mean(abs(L0glm_eigenQR-coefs_true)) # 67.84128
sum(coefs_true!=0) # 1000
sum(L0glm_eigenQR!=0) # 1033
length(intersect(which(coefs_true>0), which(as.vector(L0glm_eigenQR)>0))) # 679


# LSCG SOLVER ####
system.time(L0glm_eigenLSCG <- L0glm.fit(X=Xsp / X.n.sqL2, y=y, 
                                         start=as.vector(beta_init_ridge), # y, # as.vector(beta_init_ridge), # beta_init_nnwls_osqp2,
                                         family = poisson(identity),
                                         solver = "eigen", method = 4, 
                                         lambda = 1, nonnegative = TRUE,
                                         control.iwls = do.call("control.iwls.gen", list(maxit = 1)),
                                         control.l0 = do.call("control.l0.gen", 
                                                              list(maxit = 100, 
                                                                   delta = 1E-10, 
                                                                   gamma = 2)),
                                         control.solver = do.call("control.solver.gen", 
                                                                  list(solver = "eigen", method = 4, 
                                                                       rowaug = FALSE
                                                                       )))$coefficients / X.n.sqL2
) # 0.84s for n=p=10000
cor(L0glm_eigenLSCG, coefs_true) # 0.683617
plot(L0glm_eigenLSCG, coefs_true)
mean(abs(L0glm_eigenLSCG-coefs_true)) # 171
sum(coefs_true!=0) # 1000
sum(L0glm_eigenLSCG!=0) # 1260
