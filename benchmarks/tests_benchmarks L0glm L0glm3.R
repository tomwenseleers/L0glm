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

thresh = 1 # coefficients below this threshold will be set to 0

# couple of possible coefficient initialisations ####

# Jacobi/diagonal preconditioner initialisation ####
microbenchmark(beta_init_jacobi <- as.vector(crossprod(Xspwaug, ywaug) / colSums(Xspwaug^2)), unit="s") 
# 0.001s, colSums(Xspwaug^2) = diag(crossprod(Xspwaug)))

# ridge regression initialisation ####
microbenchmark(beta_init_ridge <- solve_sparse_LSCG(A = Xspwaug, # / X.n.L2,
                                                    b = ywaug,
                                                    x0 = beta_init_jacobi, # TODO check if this could be good default initialisation # ywaug[1:ncol(Xspwaug)], # or even better: use beta_init below in section diagonal/jacobi preconditioner
                                                    tol = 1E-3,
                                                    maxit = 0L,
                                                    precond = 1L)$coefficients, # / X.n.L2, 
               unit="s") 
# 0.012s
beta_init_ridge[beta_init_ridge<thresh] = 0
cor(as.vector(beta_init_ridge), coefs_true) # 0.201394
hist(as.vector(beta_init_ridge)) #
length(intersect(which(coefs_true>0), which(as.vector(beta_init_ridge)>thresh))) # 995, nearly all true coefficients that are >0 come out in inferred set
sum(as.vector(beta_init_ridge)>thresh) # 6386

cor(beta_init_jacobi, beta_init_ridge) # 0.74, highly correlated to true ridge solution, beta_init_jacobi all > 0 though
plot(beta_init_jacobi, beta_init_ridge)

cor(beta_init_jacobi, coefs_true) # 0.12
plot(beta_init_jacobi, coefs_true)
length(intersect(which(coefs_true>0), which(beta_init_jacobi>thresh))) # 1000, all true coefficients that are >0 come out in inferred set
sum(beta_init_jacobi>thresh) # 8058


# osqp nonnegative weighted least square initialisation ####
system.time(beta_init_nnwls_osqp <- L0glm.fit(X=Xsp / X.n.L2, y=y, weights = weights, family = gaussian(identity),
                                              start = as.vector(beta_init_ridge),
                                              solver = "osqp", lambda = 0, nonnegative = TRUE,
                                              control.iwls = do.call("control.iwls.gen", list(maxit = 0)),
                                              control.l0 = do.call("control.l0.gen", 
                                                                   list(maxit = 2, delta = 1E-10, gamma = 2)),
                                              control.solver = do.call("control.solver.gen", 
                                                                       list(solver = "osqp", rowaug = TRUE)))$coefficients 
                                     / X.n.L2
) # 0.78s
hist(log10(beta_init_nnwls_osqp))
length(intersect(which(coefs_true>0), which(as.vector(beta_init_nnwls_osqp)>thresh))) # 916, not all true coefficients that are >0 come out in inferred set
sum(as.vector(beta_init_nnwls_osqp)>thresh) # 1882
beta_init_nnwls_osqp[beta_init_nnwls_osqp < thresh] = 0
cor(beta_init_nnwls_osqp, coefs_true) # 0.9647845
plot(beta_init_nnwls_osqp, coefs_true)
mean(abs(beta_init_nnwls_osqp-coefs_true)) # 48.35177
sum(coefs_true!=0) # 1000
sum(beta_init_nnwls_osqp!=0) # 1882
length(intersect(which(coefs_true>0), which(as.vector(beta_init_nnwls_osqp)>0))) # 916


system.time(beta_init_nnwls_osqp2 <- lmridge.osqp(Xw = Xsp / X.n.L2, yw = y, weights = weights, lambdas = rep(0, ncol(X)), 
                                                  start = as.vector(beta_init_jacobi), # start=y or beta_init_ridge both OK
                                                  lower = lower, upper = upper, 
                                                  control.solver = do.call("control.solver.gen", 
                                                                   list(solver = "osqp", rowaug = FALSE))) / X.n.L2 ) # 0.94s
# 0.77s, note: figure out why call above is so much slower
# CHECK WHY SOLUTION BELOW IS SO BAD
hist(log10(beta_init_nnwls_osqp2))
length(intersect(which(coefs_true>0), which(as.vector(beta_init_nnwls_osqp2)>thresh))) # 942, nearly all true coefficients that are >0 come out in inferred set
sum(as.vector(beta_init_nnwls_osqp2)>thresh) # 2684
beta_init_nnwls_osqp2[(beta_init_nnwls_osqp2) < thresh] = 0
cor(beta_init_nnwls_osqp2, coefs_true) # 0.7963078
plot(beta_init_nnwls_osqp2, coefs_true)
mean(abs(beta_init_nnwls_osqp2-coefs_true)) # 166.8044
sum(coefs_true!=0) # 1000
sum(beta_init_nnwls_osqp2!=0) # 2684
hist(log10(beta_init_nnwls_osqp2))
length(intersect(which(coefs_true>0), which(as.vector(beta_init_nnwls_osqp)>0))) # 916


# L0glm with OSQP solver ####
system.time(L0glm_osqp <- L0glm.fit(X=Xsp / X.n.L2, y=y, 
                                    start=as.vector(beta_init_jacobi/ X.n.L2), # beta_init_jacobi, beta_init_nnwls_osqp2, OR y,
                                    family = poisson(identity),
                                    solver = "osqp", lambda = 1, nonnegative = TRUE,
                                    control.iwls = do.call("control.iwls.gen", list(maxit = 1)),
                                    control.l0 = do.call("control.l0.gen", 
                                                         list(maxit = 100, delta = 1E-10, gamma = 2)),
                                    control.solver = do.call("control.solver.gen", 
                                                             list(solver = "osqp", rowaug = FALSE)))$coefficients 
                            / X.n.L2 
            ) # 2.22s for n=p=10000
hist(log10(L0glm_osqp))
length(intersect(which(coefs_true>0), which(as.vector(L0glm_osqp)>thresh))) # 830, nearly all true coefficients that are >0 come out in inferred set

cor(L0glm_osqp, coefs_true) # 0.9569505
plot(L0glm_osqp, coefs_true)
mean(abs(L0glm_osqp-coefs_true)) # 45.4362
sum(coefs_true!=0) # 1000
sum(L0glm_osqp!=0) # 1080
length(intersect(which(coefs_true>0), which(as.vector(L0glm_osqp)>0))) # 830


# L0glm with Glmnet solver ####
system.time(L0glm_glmnet <- L0glm.fit(X=Xsp/X.n.L2, y=y,
                                      start=as.vector(beta_init_jacobi/X.n.L2), # y, # as.vector(beta_init_ridge), # beta_init_nnwls_osqp2,
                                      family = poisson(identity),
                                      solver = "glmnet", lambda = 1, nonnegative = TRUE,
                                      control.iwls = do.call("control.iwls.gen", list(maxit = 1)),
                                      control.l0 = do.call("control.l0.gen", list(maxit = 1000, delta = 1E-10, gamma = 2)))$coefficients 
                              / X.n.L2
) # 4.65s for n=p=10000
cor(L0glm_glmnet, coefs_true) # 0.960394
plot(L0glm_glmnet, coefs_true)
mean(abs(L0glm_glmnet-coefs_true)) # 40.80738
sum(coefs_true!=0) # 1000
sum(L0glm_glmnet!=0) # 1092
length(intersect(which(coefs_true>0), which(as.vector(L0glm_glmnet)>0))) # 842


# L0glm with Eigen Cholesky solver ####
system.time(L0glm_eigenCHOL <- L0glm.fit(X=Xsp, # / X.n.L2, 
                                         y=y, 
                                         start=as.vector(beta_init_ridge), # y, # as.vector(beta_init_ridge), # beta_init_nnwls_osqp2,
                                         family = poisson(identity),
                                         solver = "eigen", method = 7, 
                                         lambda = 1, nonnegative = TRUE,
                                         control.iwls = do.call("control.iwls.gen", list(maxit = 1)),
                                         control.l0 = do.call("control.l0.gen", list(maxit = 200, delta = 1E-10, gamma = 2)),
                                         control.solver = do.call("control.solver.gen", 
                                                                  list(solver = "eigen", method = 7,
                                                                       rowaug = FALSE)))$coefficients 
            # / X.n.L2
) # 0.55s for n=p=10000
cor(L0glm_eigenCHOL, coefs_true) # 0.9608598
plot(L0glm_eigenCHOL, coefs_true)
mean(abs(L0glm_eigenCHOL-coefs_true)) # 42.62425
sum(coefs_true!=0) # 1000
sum(L0glm_eigenCHOL!=0) # 1049
length(intersect(which(coefs_true>0), which(as.vector(L0glm_eigenCHOL)>0))) # 807


# L0glm with Eigen CG solver ####
system.time(L0glm_eigenCG <- L0glm.fit(X=Xsp, # / X.n.L2, 
                                       y=y, 
                                       start=as.vector(beta_init_ridge), # y, # as.vector(beta_init_ridge), # beta_init_nnwls_osqp2,
                                       family = poisson(identity),
                                       solver = "eigen", method = 5, 
                                       lambda = 1, nonnegative = TRUE,
                                       control.iwls = do.call("control.iwls.gen", list(maxit = 1)),
                                       control.l0 = do.call("control.l0.gen", list(maxit = 200, delta = 1E-10, gamma = 2)),
                                       control.solver = do.call("control.solver.gen", 
                                                                list(solver = "eigen", method = 5, rowaug = FALSE)))$coefficients 
            # / X.n.L2
) # 1.43s for n=p=10000
cor(L0glm_eigenCG, coefs_true) # 0.9608598
plot(L0glm_eigenCG, coefs_true)
mean(abs(L0glm_eigenCG-coefs_true)) # 42.62425
sum(coefs_true!=0) # 1000
sum(L0glm_eigenCG!=0) # 1049
length(intersect(which(coefs_true>0), which(as.vector(L0glm_eigenCG)>0))) # 807

# L0glm with QR solver ####
system.time(L0glm_eigenQR <- L0glm.fit(X=Xsp, # / X.n.L2, 
                                       y=y, 
                                       start=as.vector(beta_init_ridge), # y, # as.vector(beta_init_ridge), # beta_init_nnwls_osqp2,
                                       family = poisson(identity),
                                       solver = "eigen", method = 2, 
                                       lambda = 1, nonnegative = TRUE,
                                       control.iwls = do.call("control.iwls.gen", list(maxit = 1)),
                                       control.l0 = do.call("control.l0.gen", 
                                                            list(maxit = 100, delta = 1E-10, gamma = 2)),
                                       control.solver = do.call("control.solver.gen", 
                                                                list(solver = "eigen", method = 2, 
                                                                     rowaug = TRUE)))$coefficients 
                                # / X.n.L2
) # 151s for n=p=10000
cor(L0glm_eigenQR, coefs_true) # 0.9608598
plot(L0glm_eigenQR, coefs_true)
mean(abs(L0glm_eigenQR-coefs_true)) # 42.62425
sum(coefs_true!=0) # 1000
sum(L0glm_eigenQR!=0) # 1049
length(intersect(which(coefs_true>0), which(as.vector(L0glm_eigenQR)>0))) # 807


# L0glm with LSCG solver ####
system.time(L0glm_eigenLSCG <- L0glm.fit(X=Xsp, # / X.n.sqL2, 
                                         y=y, 
                                         start=as.vector(beta_init_ridge), # y, # as.vector(beta_init_ridge), # beta_init_nnwls_osqp2,
                                         family = poisson(identity),
                                         solver = "eigen", method = 4, 
                                         lambda = 1, nonnegative = TRUE,
                                         # block.size = 500,
                                         # control.block = list(block.size = 500, maxit = 1, tol = 1E-4),
                                         control.iwls = do.call("control.iwls.gen", list(maxit = 1)),
                                         control.l0 = do.call("control.l0.gen", 
                                                              list(maxit = 100, 
                                                                   delta = 1E-3, 
                                                                   gamma = 2)),
                                         control.solver = do.call("control.solver.gen", 
                                                                  list(solver = "eigen", method = 4, 
                                                                       rowaug = FALSE
                                                                       )))$coefficients 
                             # / X.n.sqL2
) # 3.2s for n=p=10000
cor(L0glm_eigenLSCG, coefs_true) # 0.9608598
plot(L0glm_eigenLSCG, coefs_true)
mean(abs(L0glm_eigenLSCG-coefs_true)) # 42.62425
sum(coefs_true!=0) # 1000
sum(L0glm_eigenLSCG!=0) # 1053
length(intersect(which(coefs_true>0), which(as.vector(L0glm_eigenLSCG)>0))) # 809


