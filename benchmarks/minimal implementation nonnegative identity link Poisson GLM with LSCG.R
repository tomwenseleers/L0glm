devtools::load_all(".")

# minimal implementation of nonnegative L0 penalized identity link Poisson GLM
nnL0poisglm = function(X, 
                       y, 
                       lambda = 1, 
                       delta = 1E-10, 
                       gamma = 2,
                       thresh = 1,
                       solver = "eigen", # or "osqp" or "glmnet"
                       method = 7L) { # 3(+4) = LSCG, 5+6 = CG, 7+8 = Cholesky  
  # print(sum(sim$beta_true!=0))
  
  weights = 1/(as.vector(y)+0.1) # observation weights (approx 1/Poisson noise variance)
  yw = y*sqrt(weights) # weighted y
  n = nrow(X)
  p = ncol(X)
  ywaug = c(yw, rep(0, p)) # augmented weighted y
  pen.weights = rep(1, p) # ridge regression penalty weights on individual coefficients (sum to p)
  lambdas = lambda*pen.weights # ridge penalties on individual coefficients
  Xw = Diagonal(x=sqrt(weights)) %*% X # = X * sqrt(weights)
  Xwaug = rbind(Xw, Diagonal(x=sqrt(lambdas))) # row augmented covariate matrix
  # require(microbenchmark)
  beta_hat <- lmridge(X = Xwaug, y = ywaug, # we use fast LSCG as initialisation
                                        pen.weights = pen.weights,
                                        lambdas = rep(0, p), # lambdas, # TODO really slow if lambdas are given & Xw and yw is passed
                                        start = as.vector(crossprod(Xwaug, ywaug) / colSums(Xwaug^2)),
                                        lower = rep(0, p), 
                                        upper = rep(Inf, p),
                                        solver = "eigen", # PS something wrong with solver="glmnet"
                                        method = 7L) # 0.03s, Eigen Cholesky was fastest, so we use that for initialisation
  # print(range(beta_hat))
  
  for (i in 1:20) {  
    nonzero = as.vector(beta_hat!=0) # nonzero coefficients
    # print(sum(nonzero))
    pen.weights[nonzero] = 1/(1+abs(beta_hat[nonzero])^gamma+delta^gamma) # update L2 norm ridge regression penalty weights
    # print(range(pen.weights[nonzero]))    
    lambdas[nonzero] = lambda*pen.weights[nonzero] # updated ridge regression penalties lambdas
    # print(range(lambdas[nonzero]))
    Xwaug = rbind(Xw[,nonzero,drop=F], Diagonal(x=sqrt(lambdas[nonzero])))
    ywaug = c(yw, rep(0, sum(nonzero)))
    
    # X.n = sqrt(colSums(Xwaug^2))
    # X.n = rep(1, sum(nonzero))
    
    beta_hat[nonzero] <- lmridge(X = Xwaug, y = ywaug, # we use fast LSCG as initialisation
                                 pen.weights = rep(1, ncol(Xwaug)),
                                 lambdas = rep(0, ncol(Xwaug)),
                                 start = as.vector(crossprod(Xwaug, ywaug) / colSums(Xwaug^2)),
                                 lower = rep(0, ncol(Xwaug)), 
                                 upper = rep(Inf, ncol(Xwaug)),
                                 solver = solver, # PS something wrong with solver="glmnet"
                                 method = method) # 0.03s, Eigen Cholesky was fastest, so we use that for initialisation
    # pred = as.vector(beta_hat %*% Xw)
    # weights = 1/pred # updated 1/variance weights for Poisson
    # Xw = Diagonal(x=sqrt(weights)) %*% X # updated Xw
    # yw = as.vector(y*sqrt(weights)) # updated yw
  }
  return(beta_hat)
}

# example
set.seed(1)
s = 0.1
p <- n <- 10000
k = round(p*s)
sim = simulate_spike_train(n=n, p=p, k=k, 
                           mean_beta = 1000, 
                           sd_logbeta = 1, 
                           family="poisson",
                           sparse=TRUE) # blurred spike train with 10000 timepoints & 200 blurred peaks
weights = 1/(sim$y+0.1) # approx 1/variance Poisson weights
Xw = Diagonal(x=sqrt(weights)) %*% sim$X # = X * sqrt(weights)
Xwaug = rbind(Xw, Diagonal(x=sqrt(rep(1,p)))) # row augmented covariate matrix
ywaug = as.vector(c(sim$y, rep(0, p))) # 


# couple of possible coefficient initialisations to use as starting values in L0glm.fit ####
# Jacobi/diagonal preconditioner initialisation ####
library(microbenchmark)
microbenchmark(beta_init_jacobi <- as.vector(crossprod(Xwaug, ywaug) / colSums(Xwaug^2)), unit="s") 
# 0.001s, colSums(Xwaug^2) = diag(crossprod(Xwaug)))

# ridge regression initialisation ####
microbenchmark(beta_init_ridge <- solve_sparse_chol(A = as(crossprod(Xwaug), "dgCMatrix"), # / X.n.L2,
                                                    b = as.vector(crossprod(Xwaug, ywaug)),
                                                    # x0 = beta_init_jacobi, # TODO check if this could be good default initialisation # ywaug[1:ncol(Xspwaug)], # or even better: use beta_init below in section diagonal/jacobi preconditioner
                                                    # tol = 1E-3,
                                                    pivot = 1L,
                                                    ordering = 0L
                                                    #maxit = 0L,
                                                    #precond = 1L
                                                    )$coefficients, # / X.n.L2, 
               unit="s") 
# 0.038s

library(microbenchmark)
system.time(coefs_lscg <- nnL0poisglm(X = sim$X,
            y = sim$y, 
            lambda = 1, 
            delta = 1E-10, 
            gamma = 2,
            thresh = 1,
            solver="eigen", 
            method=3L)) # 0.65s, LSCG
system.time(coefs_cg <- nnL0poisglm(X = sim$X,
                                 y = sim$y, 
                                 lambda = 1, 
                                 delta = 1E-10, 
                                 gamma = 2,
                                 thresh = 1,
                                 solver="eigen", 
                                 method=5L)) # 1.24s, CG
system.time(coefs_chol <- nnL0poisglm(X = sim$X,
                                 y = sim$y, 
                                 lambda = 1, 
                                 delta = 1E-10, 
                                 gamma = 2,
                                 thresh = 1,
                                 solver="eigen", 
                                 method=7L)) # 0.2s - Cholesky the fastest
# L0glm with Eigen Cholesky solver ####
system.time(coefs_chol2 <- L0glm.fit(X=sim$X, # / X.n.L2, 
                                     y=as.vector(sim$y),
                                     weights=1/(as.vector(sim$y)+0.1),
                                     # start=as.vector(beta_init_ridge) * X.n.L2, # y, # as.vector(beta_init_ridge), # beta_init_nnwls_osqp2,
                                     # family = poisson(identity),
                                     family = gaussian(identity),
                                     solver = "eigen", method = 7L, 
                                     lambda = 1, nonnegative = TRUE,
                                     control.iwls = do.call("control.iwls.gen", list(maxit = 1)),
                                     control.l0 = c(do.call("control.l0.gen", 
                                                          list(maxit = 20, delta = 1E-10, gamma = 2)),
                                                    thresh=1),
                                     control.solver = do.call("control.solver.gen", 
                                                                  list(solver = "eigen", method = 7L,
                                                                       rowaug = FALSE)))$coefficients 
            # / X.n.L2
) # 0.25s for n=p=10000

X.n.L2 <- sqrt(colSums(sim$X^2)) # L2 norm column normalisation
system.time(coefs_osqp <- nnL0poisglm(X = sim$X / X.n.L2,
                                 y = sim$y, 
                                 lambda = 1, 
                                 delta = 1E-10, 
                                 gamma = 2,
                                 thresh = 1,
                                 solver="osqp") 
            / X.n.L2
            ) # 2.4s, 9.4s in MRO

# L0glm with osqp solver ####
X.n.L2 <- sqrt(colSums(sim$X^2)) # L2 norm column normalisation
system.time(coefs_osqp2 <- L0glm.fit(X=sim$X / X.n.L2, 
                                     y=as.vector(sim$y),
                                     weights=1/(as.vector(sim$y)+0.1),
                                     # start=as.vector(beta_init_ridge)*X.n.L2, # y, # as.vector(beta_init_ridge), # beta_init_nnwls_osqp2,
                                     # family = poisson(identity),
                                     family = gaussian(identity),
                                     solver = "osqp", # method = 7L, 
                                     lambda = 1, nonnegative = TRUE,
                                     control.iwls = do.call("control.iwls.gen", list(maxit = 2)), # 2 gave better solution
                                     control.l0 = c(do.call("control.l0.gen", 
                                                            list(maxit = 20, delta = 1E-10, gamma = 2)),
                                                    thresh=1))$coefficients 
             / X.n.L2
) # 1.54s for n=p=10000

system.time(coefs_glmnet <- nnL0poisglm(X = sim$X,
                                      y = sim$y, 
                                      lambda = 1, 
                                      delta = 1E-10, 
                                      gamma = 2,
                                      thresh = 1,
                                      solver="glmnet")) # 1.71s

mean(abs(coefs_chol-coefs_lscg)) # 0.04224502 with tol=1E-5, 3.38143e-07 with tol=1E-10
mean(abs(coefs_chol-coefs_cg)) # 1.783757e-11
mean(abs(coefs_chol-coefs_osqp)) # 19.01243
mean(abs(coefs_chol-coefs_osqp2)) # 17.34572
mean(abs(coefs_chol-coefs_glmnet)) # 17.36323
mean(abs(coefs_chol-coefs_chol2)) # 0.227413

mean(abs(coefs_lscg-sim$beta_true)) # 37.55488
mean(abs(coefs_cg-sim$beta_true)) # 37.54006
mean(abs(coefs_chol-sim$beta_true)) # 37.54006
mean(abs(coefs_chol2-sim$beta_true)) # 37.59185
mean(abs(coefs_osqp-sim$beta_true)) # 41.85887
mean(abs(coefs_osqp2-sim$beta_true)) # 40.32613
mean(abs(coefs_glmnet-sim$beta_true)) # 44.86981

plot(sim$beta_true, coefs_chol, pch=16, col="steelblue")
plot(sim$beta_true, coefs_chol2, pch=16, col="steelblue")

thresh = 1

length(intersect(which(sim$beta_true>0), which(as.vector(coefs_lscg)>thresh))) # 859 (should be 1000)
length(intersect(which(sim$beta_true>0), which(as.vector(coefs_cg)>thresh))) # 859 (should be 1000)
length(intersect(which(sim$beta_true>0), which(as.vector(coefs_chol)>thresh))) # 859 (should be 1000)
length(intersect(which(sim$beta_true>0), which(as.vector(coefs_chol2)>thresh))) # 860 (should be 1000)
length(intersect(which(sim$beta_true>0), which(as.vector(coefs_osqp)>thresh))) # 851 (should be 1000)
length(intersect(which(sim$beta_true>0), which(as.vector(coefs_osqp2)>thresh))) # 855 (should be 1000)
length(intersect(which(sim$beta_true>0), which(as.vector(coefs_glmnet)>thresh))) # 851 (should be 1000)

sum(coefs_lscg>thresh) # 1079 (should be 1000)
sum(coefs_cg>thresh) # 1079 (should be 1000)
sum(coefs_chol>thresh) # 1079 (should be 1000)
sum(coefs_chol2>thresh) # 1082 (should be 1000)
sum(coefs_osqp>thresh) # 1090 (should be 1000)
sum(coefs_osqp2>thresh) # 1097 (should be 1000)
sum(coefs_glmnet>thresh) # 1098 (should be 1000)

