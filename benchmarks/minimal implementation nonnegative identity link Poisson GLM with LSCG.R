devtools::load_all(".")

# minimal implementation of nonnegative L0 penalized identity link Poisson GLM
nnL0poisglm = function(X, 
                       y, 
                       lambda = 1, 
                       delta = 1E-10, 
                       gamma = 2,
                       thresh = 1,
                       solver = "eigen", # or "osqp" or "glmnet"
                       method = ifelse(ncol(X)<20000, 7L, 3L), # = Cholesky # 3(+4) = LSCG, 5(+6) = CG, 7(+8) = Cholesky  
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
                                 lambdas = NULL, # rep(0, ncol(Xwaug)),
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
    # print(previous_norm)
    # print(current_norm)
    relative_difference <- abs(current_norm - previous_norm) / (previous_norm+0.1)
    
    # print(relative_difference)
    if (relative_difference < reltol) break
  }
  return(beta_hat)
}

# 1. EXAMPLE SINGLE CHANNEL ####
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
# orthogonal initialisation (assuming variables are orthogonal) ####
library(microbenchmark)
microbenchmark(beta_init_orthogonal <- as.vector(crossprod(Xwaug, ywaug) / colSums(Xwaug^2)), unit="s") 
# 0.001s, colSums(Xwaug^2) = diag(crossprod(Xwaug)))
cor(beta_init_orthogonal, as.vector(sim$beta_true)) # 0.336
plot(beta_init_orthogonal, sim$beta_true)
# PS original variables in columns of X are obviously very far from orthogonal, 
# but the augmentation with a diagonal matrix with sqrt(lambdas) along the diagonal
# makes the assumption hold better

# diagonal initialisation (assuming system is dominated by diagonal of X) ####
# for our problem this gives the best initial solution
# only works if X is square though
microbenchmark(beta_init_diagonal <- y*sqrt(weights) / diag(Xw))
# 104µs
cor(beta_init_diagonal, as.vector(sim$beta_true)) # 0.4709823
plot(beta_init_diagonal, sim$beta_true)

# ridge regression initialisation ####
microbenchmark(beta_init_ridge <- solve_sparse_chol(A = as(crossprod(Xwaug), "generalMatrix"),
                                                    b = as.vector(crossprod(Xwaug, ywaug)),
                                                    # x0 = beta_init_orthogonal,
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
                                      warmupiters = warmupiters)) # LSCG, 0.65s with tol=1E-5 or 1.30s with tol=1E-10 & 5.77 in other comp

system.time(coefs_cg <- nnL0poisglm(X = sim$X,
                                    y = sim$y, 
                                    lambda = 1, 
                                    delta = 1E-10, 
                                    gamma = 2,
                                    thresh = 1,
                                    solver="eigen", 
                                    method=5L,
                                    warmupiters = warmupiters)) # 0.48s, CG / 2.64 in other comp

system.time(coefs_chol <- nnL0poisglm(X = sim$X,
                                      y = sim$y, 
                                      lambda = 1, 
                                      delta = 1E-10, 
                                      gamma = 2,
                                      thresh = 1,
                                      solver="eigen", 
                                      method=7L,
                                      warmupiters = warmupiters)) # 0.16s - Cholesky : FASTEST // 0.42 in other comp


# L0glm with Eigen Cholesky solver - NOTE weighted least squares gives better solution than Poisson identity link ####
system.time(coefs_chol2 <- L0glm.fit(X=sim$X, 
                                     y=as.vector(sim$y),
                                     weights=1/(as.vector(sim$y)+0.1),
                                     # start=as.vector(beta_init_orthogonal),
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
) # 0.22s for n=p=10000 / 0.46 in other comp

system.time(coefs_osqp <- nnL0poisglm(X = sim$X, 
                                      y = sim$y, 
                                      lambda = 1, 
                                      delta = 1E-10, 
                                      gamma = 2,
                                      thresh = 1,
                                      solver="osqp",
                                      warmupiters = warmupiters)
) # 2.68s / 5.63 in other comp

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
) # 4.14s for n=p=10000 / 7.5 in other comp

system.time(coefs_glmnet <- nnL0poisglm(X = sim$X,
                                        y = sim$y, 
                                        lambda = 1, 
                                        delta = 1E-10, 
                                        gamma = 2,
                                        thresh = 1,
                                        solver="glmnet",
                                        warmupiters = warmupiters)) # 1.67s / 2.99 in other comp

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



# 2. EXAMPLE MULTICHANNEL ####

library(export)
library(qlcMatrix)

set.seed(1)
sim <- simulate_spike_train(n=100, p=100, k=10, multichannel=T, sparse=TRUE)
X <- sim$X
Y <- sim$y
true_coefs <- sim$beta_true
m <- ncol(Y)
n <- nrow(X)
p <- ncol(X)

# image(X, useRaster=TRUE)
# graph2png(file="simulated_chromatogram.png", width=7, height=5)



# ITERATED L0 penalized nonnegative Poisson GLMs FOR MULTICHANNEL DECONVOLUTION ####

# 2.1 USING BASIC APPLY, I.E. NOT PARALLELIZED ####

# L0glm fit applied separately to each of the columns of y :
# (see section 2.2 for a parallelized version of this)
# NOTE I think L0glm.fit should return this solution by default if you pass it
# a Y matrix, this solution could be called "independent L0" by analogy with https://hichamjanati.github.io/mutar/api.html

system.time(beta_L0_1 <- do.call(cbind, lapply(1:ncol(sim$y), function(col) {
  # print(col)
  nnL0poisglm(X = sim$X,
              y = as.vector(sim$y[,col,drop=F]),
              lambda = 1,
              delta = 1E-10,
              gamma = 2,
              thresh = 1,
              solver = "osqp",  # TODO check why glmnet & osqp crash or give error. Osqp crashes R, glmnet gives error because sum(Xwaug) == 0 before calling bigGLM in lmridge
              method = 7L,
              warmupiters = 4L)
  # L0glm.fit(X = sim$X,
  #           y = as.vector(sim$y[,col,drop=F]),
  #           )
}))) # 11.7s with eigen solver 7 (Cholesky) / 35s with other computer, TODO check why glmnet & osqp crash or give error
cor(as.vector(beta_L0_1), as.vector(sim$beta_true)) # 0.9979857
# NOTE I think L0glm.fit should return this solution by default if you pass it
# a Y matrix, this solution could be called "independent L0" by analogy with https://hichamjanati.github.io/mutar/api.html

image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=beta_L0_1^1E-10, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="L0 nonneg Poisson GLM coefficients")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim$X))[as.vector(rowMin(beta_L0_1)!=0)], col="cyan")
sum(rowMax(sim$beta_true)!=0) # 100 true peaks
sum(rowMin(beta_L0_1)!=0) # 299 peaks detected - still too high
length(intersect(which(rowMax(sim$beta_true)!=0), which(rowMin(beta_L0_1)!=0))) # 99 out of 100 peaks detected

# 2nd step using 1/ avg coefficients over channels from previous fit as starting values for another L0 GLM fit
avgcoefs = rowMeans(beta_L0_1) # PS as(rowMax(beta_L0_1),"vector") gives worse solution
sum(avgcoefs>0) # 299 selected variables in 1st pass
variable_indices = which(avgcoefs>0)
beta_L0_1_2 = Matrix(0, nrow=nrow(beta_L0_1), ncol=ncol(beta_L0_1), sparse=TRUE)
system.time(beta_L0_1_2[variable_indices,] <- do.call(cbind, lapply(1:ncol(sim$y), function(col) {
  # print(col)
  L0glm.fit(X=sim$X[,variable_indices,drop=F], 
            y=as.vector(sim$y[,col,drop=F]),
            weights=1/(as.vector(sim$y[,col,drop=F])+0.1),
            start=as.vector(avgcoefs[variable_indices]),
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
  
}))) # 3.2s with eigen solver 7 (Cholesky) / 12s with other computer
cor(as.vector(beta_L0_1_2), as.vector(sim$beta_true)) # 0.9993619

image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=as.matrix(beta_L0_1_2)^1E-10, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="L0 nonneg Poisson GLM coefficients")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim$X))[as.vector(rowMin(beta_L0_1_2)!=0)], col="cyan")
sum(rowMax(sim$beta_true)!=0) # 100 true peaks
sum(rowMin(beta_L0_1_2)!=0) # 111 peaks detected - still too high 
length(intersect(which(rowMax(sim$beta_true)!=0), which(rowMin(beta_L0_1)!=0))) # 99 out of 100 peaks detected


# 3d step using 1/ avg coefficients over channels from previous fit as starting values for another L0 GLM fit
avgcoefs = rowMeans(beta_L0_1_2)
sum(avgcoefs>0) # 111 selected variables in 2nd pass
variable_indices = which(avgcoefs>0)
beta_L0_1_3 = Matrix(0, nrow=nrow(beta_L0_1), ncol=ncol(beta_L0_1), sparse=TRUE)
system.time(beta_L0_1_3[variable_indices,] <- do.call(cbind, lapply(1:ncol(sim$y), function(col) {
  # print(col)
  L0glm.fit(X=sim$X[,variable_indices,drop=F], 
            y=as.vector(sim$y[,col,drop=F]),
            weights=1/(as.vector(sim$y[,col,drop=F])+0.1),
            start=as.vector(avgcoefs[variable_indices]),
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
  
}))) # 2.87s with eigen solver 7 (Cholesky) / 15.2S in other computer
cor(as.vector(beta_L0_1_3), as.vector(sim$beta_true)) # 0.9994057

image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=as.matrix(beta_L0_1_3)^1E-10, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="L0 nonneg Poisson GLM coefficients")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim$X))[as.vector(rowMin(beta_L0_1_3)!=0)], col="cyan")
sum(rowMax(sim$beta_true)!=0) # 100 true peaks
sum(rowMin(beta_L0_1_3)!=0) # 98 peaks detected
length(intersect(which(rowMax(sim$beta_true)!=0), which(rowMin(beta_L0_1_3)!=0))) # 98 out of 100 peaks detected
sum(sim$beta_true!=0) # 7645 true nonzero coefficients
sum(beta_L0_1_3!=0) # 3976 coefficients detected - spectra too sparse...
length(intersect(which(sim$beta_true!=0), which(beta_L0_1_3!=0))) # 3971 out of 7645 = 52% nonzero coefficients detected
length(intersect(which(sim$beta_true==0), which(beta_L0_1_3!=0))) # 5 out of 7645 = 0.07% false positive coefficients

# NOTE: I would call this solution an L0,Linfinity GLM, which we can maybe just call group L0 ?


# 4th step using 1/ avg coefficients over channels from previous fit as starting values for another L0 GLM fit but now using lambda=0, i.e. using a nonnegative Poisson GLM, using variables selected earlier via iterated L0 GLMs
avgcoefs = rowMeans(beta_L0_1_3)
sum(avgcoefs>0) # 98 selected variables in 3d pass
variable_indices = which(avgcoefs>0)
beta_L0_1_4 = Matrix(0, nrow=nrow(beta_L0_1), ncol=ncol(beta_L0_1), sparse=TRUE)
system.time(beta_L0_1_4[variable_indices,] <- do.call(cbind, lapply(1:ncol(sim$y), function(col) {
  # print(col)
  L0glm.fit(X=sim$X[,variable_indices,drop=F], 
            y=as.vector(sim$y[,col,drop=F]),
            # weights=1/(as.vector(sim$y[,col,drop=F])+0.1),
            start=as.vector(avgcoefs[variable_indices]),
            family = poisson(identity), # gives worse solution
            # family = gaussian(identity), 
            solver = "osqp", method = 7L, 
            lambda = 0, nonnegative = TRUE,
            control.iwls = do.call("control.iwls.gen", list(maxit = 20)),
            control.l0 = c(do.call("control.l0.gen", 
                                   list(maxit = 20, delta = 1E-10, gamma = 2)),
                           thresh=1E-5),
            control.solver = do.call("control.solver.gen", 
                                     list(solver = "osqp", method = 7L,
                                          rowaug = TRUE)))$coefficients 
  
}))) # 2.73s with osqp solver, eigen slightly faster / 8.84s with other computer
cor(as.vector(beta_L0_1_4), as.vector(sim$beta_true)) # 0.9998314

# NOTE: this solution I would call "relaxed group L0" by analogy of a relaxed lasso, https://cran.microsoft.com/snapshot/2023-03-06/web/packages/glmnet/vignettes/relax.pdf


image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=as.matrix(beta_L0_1_4)^1E-10, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="L0 nonneg Poisson GLM coefficients")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim$X))[as.vector(rowMin(beta_L0_1_4)!=0)], col="cyan")
sum(rowMax(sim$beta_true)!=0) # 100 true peaks
sum(rowMin(beta_L0_1_4)!=0) # 98 peaks detected
length(intersect(which(rowMax(sim$beta_true)!=0), which(rowMin(beta_L0_1_4)!=0))) # 98 out of 100 peaks detected
sum(sim$beta_true!=0) # 7645 true nonzero coefficients
sum(beta_L0_1_4!=0) # 6882 coefficients detected - very good solution!!!
length(intersect(which(sim$beta_true!=0), which(beta_L0_1_4!=0))) # 6882 out of 7645 = 90% nonzero coefficients detected
length(intersect(which(sim$beta_true==0), which(beta_L0_1_4!=0))) # 104 out of 7645 = 1.4% false positive coefficients

# TODO performance should be compared to methods in
# https://hichamjanati.github.io/mutar/index.html#
# https://hichamjanati.github.io/mutar/api.html
# https://scikit-learn.org/stable/auto_examples/linear_model/plot_multi_task_lasso_support.html






# 2.2 USING PARALLELISATION USING library(parallel) ####

# Load the parallel package
library(parallel)

# Define the number of cores to use
n_cores <- round(detectCores()/2) # 8 for me

# Determine the type of cluster based on the operating system
cluster_type <- if (Sys.info()["sysname"] == "Windows") "PSOCK" else "FORK"

# Create a cluster using the number of cores
cl <- makeCluster(n_cores, type=cluster_type) # on Mac & Ubuntu type="FORK" is faster

# Export the sim data and the function to the workers
clusterExport(cl, list("sim", "fun", "nnL0poisglm"))

# Load the required libraries and functions on each worker
# (this only needs to be done once)
cleval <- clusterEvalQ(cl, {
  require(L0glm)
  devtools::load_all(".") # this should not be necessary once all callable functions are
  # correctly exported by the L0glm package
})

system.time({
  # Define the function to be applied
  fun <- function(cols) {
    # require(L0glm)
    # devtools::load_all(".")
    
    results <- lapply(cols, function(col) {
      # Get the required data
      X = sim$X
      y = as.vector(sim$y[,col,drop=F])
      
      # Call nnL0poisglm
      nnL0poisglm(X = X,
                  y = y, 
                  lambda = 1, 
                  delta = 1E-10, 
                  gamma = 2,
                  thresh = 1,
                  solver = "eigen",
                  method = 7L,
                  warmupiters = 4L)
    })
    
    # Return the results as a matrix
    do.call(cbind, results)
  }
  
  # Create a list of column index chunks
  cols <- 1:ncol(sim$y)
  chunks <- split(cols, cut(cols, n_cores))
  
  # Use clusterApply to apply the function in parallel & combine results
  beta_L0_2 <- do.call(cbind, clusterApply(cl, chunks, fun))
  
})# 4.46s using 8 cores

beta <- L0glm.fit(X = sim$X, y = sim$y)
# Stop the cluster
stopCluster(cl)

cor(as.vector(beta_L0_2), as.vector(sim$beta_true)) # 0.9979857

image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=beta_L0_2^1E-10, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="L0 nonneg Poisson GLM coefficients")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim$X))[as.vector(rowMin(beta_L0_2)!=0)], col="cyan")

sum(rowMax(sim$beta_true)!=0) # 100 true peaks
sum(rowMin(beta_L0_2)!=0) # 299 peaks detected - still too high
length(intersect(which(rowMax(sim$beta_true)!=0), which(rowMin(beta_L0_2)!=0))) # 99 out of 100 peaks detected


system.time(par <- L0glm.fit(X = sim$X, y = sim$y, lambda = 1, nonnegative = T, parallellize = T, grouped = T, relaxed = F, solver = "eigen", family = poisson(identity))) # 34.88s for parallellize = F, 19.75s for parallellize = T
cor(as.vector(par$coefficients), as.vector(sim$beta_true)) # 0.994652 with variable indices
# TODO: check why if data is family = "gaussian" gives an error

length(intersect(which(sim$beta_true!=0), which(par$coefficients!=0))) # 4515 out of 7645
length(intersect(which(sim$beta_true==0), which(par$coefficients!=0))) # 4 out of 7645

length(intersect(which(rowMax(sim$beta_true)!=0), which(rowMin(par$coefficients)!=0)))  # 99 out of 100



