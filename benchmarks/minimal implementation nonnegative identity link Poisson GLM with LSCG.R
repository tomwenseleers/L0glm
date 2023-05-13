devtools::load_all(".")

# minimal implementation of nonnegative L0 penalized identity link Poisson GLM
nnL0poisglm = function(X, 
                       y, 
                       lambda = 1, 
                       delta = 1E-10, 
                       gamma = 2,
                       thresh = 1,
                       solver = "eigen", # or "osqp" or "glmnet"
                       method = 7L, # = Cholesky # 3(+4) = LSCG, 5(+6) = CG, 7(+8) = Cholesky  
                       maxit = 25,
                       reltol = 1E-5,
                       warmupiters = 4) { # nr of iterations before obs weights are updated in addition to penalty weights in each iteration, 3 or 4 gave best performance
  # print(sum(sim$beta_true!=0))
  
  weights = 1/(as.vector(y)+0.1) # observation weights (approx 1/Poisson noise variance)
  goodrows = (!is.na(weights)&!(weights==0)&!is.infinite(weights)) # good rows to include
  yw = y*sqrt(weights[goodrows]) # weighted y
  n = sum(goodrows)
  p = ncol(X)
  ywaug = c(yw, rep(0, p)) # augmented weighted y
  pen.weights = rep(1, p) # ridge regression penalty weights on individual coefficients (sum to p)
  lambdas = lambda*pen.weights # ridge penalties on individual coefficients
  Xw = Diagonal(x=sqrt(weights[goodrows])) %*% X[goodrows,,drop=F] # = X * sqrt(weights)
  Xwaug = rbind(Xw, Diagonal(x=sqrt(lambdas))) # row augmented covariate matrix
  # require(microbenchmark)
  beta_hat <- lmridge(X = Xwaug, y = ywaug, # we use Eigen Cholesky as initialisation
                                        pen.weights = pen.weights,
                                        lambdas = NULL,
                                        start = as.vector(crossprod(Xwaug, ywaug) / colSums(Xwaug^2)),
                                        lower = rep(0, p), 
                                        upper = rep(Inf, p),
                                        solver = "eigen",
                                        method = ifelse(p>=20000, 3L, 5L)) # 0.035s, we use Eigen Cholesky (for small p) or LSCG (for large p) for initialisation
  beta_hat[abs(beta_hat)<thresh] = 0
  nonzero = as.vector(beta_hat!=0) # nonzero coefficients
  # print(range(beta_hat))

  for (i in 1:maxit) {  
    # print(i)
    # print(sum(nonzero))
    pen.weights[nonzero] = 1/(1+abs(beta_hat[nonzero])^gamma+delta^gamma) # update L2 norm ridge regression penalty weights
    # print(range(pen.weights[nonzero]))    
    lambdas[nonzero] = lambda*pen.weights[nonzero] # updated ridge regression penalties lambdas
    # print(range(lambdas[nonzero]))
    
    Xwaug = rbind(Diagonal(x=sqrt(weights[goodrows])) %*% X[goodrows,nonzero,drop=F], # Xw
                  Diagonal(x=sqrt(lambdas[nonzero])))
    ywaug = c(y[goodrows] * sqrt(weights[goodrows]), rep(0, sum(nonzero)))
    
    beta_hat_old <- beta_hat
    beta_hat[nonzero] <- lmridge(X = Xwaug, y = ywaug,
                                 pen.weights = rep(1, ncol(Xwaug)),
                                 lambdas = rep(0, ncol(Xwaug)),
                                 start = beta_hat[nonzero], # as.vector(crossprod(Xwaug, ywaug) / colSums(Xwaug^2)), 
                                 lower = rep(0, ncol(Xwaug)), 
                                 upper = rep(Inf, ncol(Xwaug)),
                                 solver = solver,
                                 method = method)
    beta_hat[abs(beta_hat)<thresh] = 0
    nonzero = as.vector(beta_hat!=0) # nonzero coefficients
    if (i>=warmupiters) { # after nr of warmup iterations weights are also updated (could also be done if model predictions have nearly converged)
      pred = as.vector(X[,nonzero,drop=F] %*% beta_hat[nonzero])
      weights = 1/pred # updated 1/variance weights for Poisson
      goodrows = (!is.na(weights)&!(weights==0)&!is.infinite(weights)) # good rows to include
    }
        
    previous_norm <- sqrt(sum(beta_hat_old^2))
    current_norm <- sqrt(sum(beta_hat^2))
    relative_difference <- abs(current_norm - previous_norm) / previous_norm
    
    # print(relative_difference)
    if (relative_difference < reltol) break
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
X = sim$X
y = sim$y
weights = 1/(sim$y+0.1) # approx 1/variance Poisson weights
Xw = Diagonal(x=sqrt(weights)) %*% sim$X # = X * sqrt(weights)
Xwaug = rbind(Xw, Diagonal(x=sqrt(rep(1,p)))) # row augmented covariate matrix
ywaug = as.vector(c(sim$y, rep(0, p))) # 
lambda = 1 
delta = 1E-10
gamma = 2
thresh = 1
solver = "eigen" # or "osqp" or "glmnet"
method = 7L # = Cholesky
maxit = 25
reltol = 1E-5
warmupiters = 4 # 3 or 4 gave best performance



# couple of possible coefficient initialisations to use as starting values in L0glm.fit ####
# Jacobi/diagonal preconditioner initialisation ####
library(microbenchmark)
microbenchmark(beta_init_jacobi <- as.vector(crossprod(Xwaug, ywaug) / colSums(Xwaug^2)), unit="s") 
# 0.001s, colSums(Xwaug^2) = diag(crossprod(Xwaug)))

# ridge regression initialisation ####
microbenchmark(beta_init_ridge <- solve_sparse_chol(A = as(crossprod(Xwaug), "generalMatrix"),
                                                    b = as.vector(crossprod(Xwaug, ywaug)),
                                                    # x0 = beta_init_jacobi,
                                                    # tol = 1E-3,
                                                    pivot = 1L,
                                                    ordering = 0L
                                                    # maxit = 0L,
                                                    # precond = 1L
                                                    )$coefficients, 
               unit="s") 
# 0.037s
beta_init_ridge[abs(beta_init_ridge)<0] = 0

library(microbenchmark)
system.time(coefs_lscg <- nnL0poisglm(X = sim$X,
            y = sim$y, 
            lambda = 1, 
            delta = 1E-10, 
            gamma = 2,
            thresh = 1,
            solver="eigen", 
            method=3L,
            warmupiters = warmupiters)) # LSCG, 0.65s with tol=1E-5 or 1.30s with tol=1E-10
system.time(coefs_cg <- nnL0poisglm(X = sim$X,
                                 y = sim$y, 
                                 lambda = 1, 
                                 delta = 1E-10, 
                                 gamma = 2,
                                 thresh = 1,
                                 solver="eigen", 
                                 method=5L,
                                 warmupiters = warmupiters)) # 0.48s, CG
system.time(coefs_chol <- nnL0poisglm(X = sim$X,
                                 y = sim$y, 
                                 lambda = 1, 
                                 delta = 1E-10, 
                                 gamma = 2,
                                 thresh = 1,
                                 solver="eigen", 
                                 method=7L,
                                 warmupiters = warmupiters)) # 0.16s - Cholesky : FASTEST

 
# L0glm with Eigen Cholesky solver - NOTE weighted least squares gives better solution than Poisson identity link ####
system.time(coefs_chol2 <- L0glm.fit(X=sim$X, 
                                     y=as.vector(sim$y),
                                     weights=1/(as.vector(sim$y)+0.1),
                                     # start=as.vector(beta_init_ridge),
                                     # family = poisson(identity), # gives worse solution
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
) # 0.18s for n=p=10000

system.time(coefs_osqp <- nnL0poisglm(X = sim$X, 
                                      y = sim$y, 
                                      lambda = 1, 
                                      delta = 1E-10, 
                                      gamma = 2,
                                      thresh = 1,
                                      solver="osqp",
                                      warmupiters = warmupiters)
) # 2.68s

# L0glm with osqp solver ####
system.time(coefs_osqp2 <- L0glm.fit(X=sim$X, 
                                     y=as.vector(sim$y),
                                     weights=1/(as.vector(sim$y)+0.1),
                                     # start=as.vector(beta_init_ridge), 
                                     # family = poisson(identity), # gives worse solution
                                     family = gaussian(identity),
                                     solver = "osqp", 
                                     lambda = 1, nonnegative = TRUE,
                                     control.iwls = do.call("control.iwls.gen", list(maxit = 1)),
                                     control.l0 = c(do.call("control.l0.gen", 
                                                            list(maxit = 20, delta = 1E-10, gamma = 2)),
                                                    thresh=1))$coefficients 
) # 4.14s for n=p=10000

system.time(coefs_glmnet <- nnL0poisglm(X = sim$X,
                                      y = sim$y, 
                                      lambda = 1, 
                                      delta = 1E-10, 
                                      gamma = 2,
                                      thresh = 1,
                                      solver="glmnet",
                                      warmupiters = warmupiters)) # 1.67s

mean(abs(coefs_chol-coefs_lscg)) # 0.04224502 with tol=1E-5, 7.718954e-07 with tol=1E-10
mean(abs(coefs_chol-coefs_cg)) # 7.722587e-07
mean(abs(coefs_chol-coefs_osqp)) # 10.87642
mean(abs(coefs_chol-coefs_osqp2)) # 12.78663
mean(abs(coefs_chol-coefs_glmnet)) # 17.08276
mean(abs(coefs_chol-coefs_chol2)) # 5.520426

mean(abs(coefs_lscg-sim$beta_true)) # 36.83137
mean(abs(coefs_cg-sim$beta_true)) # 36.83137
mean(abs(coefs_chol-sim$beta_true)) # 36.83137
mean(abs(coefs_chol2-sim$beta_true)) # 37.59185
mean(abs(coefs_osqp-sim$beta_true)) # 39.19337
mean(abs(coefs_osqp2-sim$beta_true)) # 37.72411
mean(abs(coefs_glmnet-sim$beta_true)) # 42.75899

plot(sim$beta_true, coefs_chol, pch=16, col="steelblue")
plot(sim$beta_true, coefs_chol2, pch=16, col="steelblue")

thresh = 1

length(intersect(which(sim$beta_true>0), which(as.vector(coefs_lscg)>thresh))) # 859 (should be 1000)
length(intersect(which(sim$beta_true>0), which(as.vector(coefs_cg)>thresh))) # 859 (should be 1000)
length(intersect(which(sim$beta_true>0), which(as.vector(coefs_chol)>thresh))) # 859 (should be 1000)
length(intersect(which(sim$beta_true>0), which(as.vector(coefs_chol2)>thresh))) # 860 (should be 1000)
length(intersect(which(sim$beta_true>0), which(as.vector(coefs_osqp)>thresh))) # 850 (should be 1000)
length(intersect(which(sim$beta_true>0), which(as.vector(coefs_osqp2)>thresh))) # 856 (should be 1000)
length(intersect(which(sim$beta_true>0), which(as.vector(coefs_glmnet)>thresh))) # 852 (should be 1000)

sum(coefs_lscg>thresh) # 1072 (should be 1000)
sum(coefs_cg>thresh) # 1072 (should be 1000)
sum(coefs_chol>thresh) # 1072 (should be 1000)
sum(coefs_chol2>thresh) # 1082 (should be 1000)
sum(coefs_osqp>thresh) # 1073 (should be 1000)
sum(coefs_osqp2>thresh) # 1088 (should be 1000)
sum(coefs_glmnet>thresh) # 1084 (should be 1000)

