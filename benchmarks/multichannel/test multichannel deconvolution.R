# TESTS L0 GROUP PENALIZED DECONVOLUTION OF MULTICHANNEL DATA

library(export)
library(qlcMatrix)

sim <- simulate_spike_train(n=5, p=5, k=2, multichannel=TRUE, sparse=TRUE)
X <- sim$X
Y <- sim$y
true_coefs <- sim$beta_true
m <- ncol(Y)
n <- nrow(X)
p <- ncol(X)


# graph2png(file="simulated_chromatogram.png", width=7, height=5)

# test implementation of L0 group penalized nonnegative identity link Poisson GLM

# wnnls result fitted on each channel independently

# 1. using basic apply ####
system.time(beta_wnnls <- apply(sim$y, 2, function(ycol) {
             lmridge(X = sim$X,
                y = as.vector(ycol),
                lambda = rep(0, ncol(sim$X)),
                weights = 1/(as.vector(ycol)+0.1),
                lower = rep(0, ncol(sim$X)),
                upper = rep(Inf, ncol(sim$X)),
                nonnegative = TRUE,
                solver = "osqp", # note: nnls or osqp solution much better than glmnet
                method = NULL,
                start = rep(0, ncol(sim$X)),
                control.solver = NULL,
                verbose = TRUE) })) # 5s with nnls, 1.3s with osqp, 185s with glmnet, for n=p=500 & k=10
# TODO check why glmnet is so slow here and gives off warnings
image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=beta_wnnls^0.1, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="wnnls coefficients (nnls)")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")
# graph2png(file="simulated_chromatogram_coefficients_wnnls.png", width=7, height=5)


# 2. using parallelized apply using library(parallel) ####

library(parallel)
n_cores <- round(parallel::detectCores()/2)

system.time({
  # Set up parallel backend
  cl <- parallel::makeCluster(n_cores, type="PSOCK") # NOTE: on Ubuntu & Mac type="FORK" is faster
  clusterExport(cl, list("sim", "lmridge")) # Export the 'sim' object to the workers

  # Run the loop in parallel
beta_wnnls <- parallel::parApply(cl, sim$y, 2, FUN = function(ycol) {
  lmridge(X = sim$X,
          y = as.vector(ycol),
          lambda = rep(0, ncol(sim$X)),
          weights = 1/(as.vector(ycol)+0.1),
          lower = rep(0, ncol(sim$X)),
          upper = rep(Inf, ncol(sim$X)),
          nonnegative = TRUE,
          solver = "nnls",
          method = NULL,
          start = rep(0, ncol(sim$X)),
          control.solver = NULL,
          verbose = TRUE)
  })

  # Clean up
  stopCluster(cl)
}) # 3.4s with nnls, 2.9s with osqp, 52s with glmnet, for n=p=500 & k=10
# Gives error invalid 'times' argument


# 3. by working with a block diagonal matrix to make multichannel problem into single channel one: ####

# TODO finish this code - ask ChatGPT
# I think this can be done using the spam library & using a Kronecker product
library(spam)
library(Matrix)


yvect = as.vector(Y) # vectorized version of multichannel data

# create the weight matrix and the X weighted matrix
W <- diag(sqrt(1/(yvect+0.1)))
Xw <- lapply(1:m, function(i) X * W[i])


# create the block diagonal matrix
Xblock <- bdiag(Xw)


lambda <- 1
lower <- rep(0, ncol(Xblock))
upper <- rep(Inf, ncol(Xblock))
nonnegative <- TRUE
start <- NULL
verbose <- TRUE
X.n <- apply(Xblock, 2, norm, type = "2")

system.time(coefs_multi <- lmridge(X = Xblock, y = yvect, lambda = lambda, lower = lower, upper = upper,
                                  nonnegative = nonnegative, solver = "osqp", method = 7)) 
#0.09s glmnet, 0.61s osqp, 82s nnls, 0.83 eigen.1, 0.15 eigen.4, (n=p=150, m=131) 

system.time(coefs_multi <- L0glm.fit2(X=Xblock, y=yvect, family = poisson(identity),
                                      solver = "glmnet", method = 4, block = FALSE,
                                      lambda = 1, nonnegative = TRUE,
                                      control.iwls = do.call("control.iwls.gen2", list(maxit = 100)),
                                      control.l0 = do.call("control.l0.gen2", list(maxit = 100, delta = 1E-10, gamma = 2)),
                                      control.solver = do.call("control.solver.gen", list(solver = "glmnet", method = 4, rowaug = TRUE))))
#206s osqp, 2.5s glmnet,  (n=p=150, m=131, 100 maxit)

# Turn coefficient vector into matrix
original_dim <- c(p, m)

coefficients_matrix <- matrix(coefs_multi$coefficients, nrow = original_dim[1], ncol = original_dim[2], byrow = FALSE)
cor(coefs_multi$coefficients, as.numeric(true_coefs)) #bad correlation



image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=matrix(coefs_multi$coefficients)^0.1, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="Coefficients")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")












# nonnegative identity link Poisson GLM fitted on each channel independently
system.time(beta_nnpois <- apply(sim$y, 2, function(ycol) {
  L0glm.fit2(X = sim$X,
             y = as.vector(ycol),
             lambda = 0,
             # weights = 1/(as.vector(ycol)+0.1),
             family = poisson(identity),
             lower = rep(0, ncol(sim$X)),
             upper = rep(Inf, ncol(sim$X)),
             nonnegative = TRUE,
             solver = "nnls",
             method = NULL,
             start = rep(0, ncol(sim$X)), # prob faster to start from nnls coefficients?
             control.solver = NULL,
             control.iwls = do.call("control.iwls.gen2", list(maxit = 10)),
             control.l0 = do.call("control.l0.gen2", list(maxit = 1, delta = 1E-10, gamma = 2)))$coefficients }))
# 10.7s with nnls, 4.3s with osqp for n=p=500 & k=10
image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=beta_nnpois^0.1, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="nonnegative ident link Poisson GLM coefficients (nnls)")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")
# graph2png(file="simulated_chromatogram_coefficients_nnpois.png", width=7, height=5)


# nonnegative L0 penalized identity link weighted LS fitted on each channel independently
system.time(beta_nnL0pois <- apply(sim$y, 2, function(ycol) {
  L0glm.fit2(X = sim$X,
          y = as.vector(ycol),
          lambda = 1,
          # weights = 1/(as.vector(ycol)+0.1),
          family = poisson(identity),
          lower = rep(0, ncol(sim$X)),
          upper = rep(Inf, ncol(sim$X)),
          nonnegative = TRUE,
          solver = "osqp",
          method = NULL,
          # start = rep(0, ncol(sim$X)), # prob faster to start from nnls coefficients?
          # control.solver = NULL,
          control.iwls = do.call("control.iwls.gen2", list(maxit = 1)),
          control.l0 = do.call("control.l0.gen2", list(maxit = 10, delta = 1E-10, gamma = 2)))$coefficients }))
# 78s with osqp
image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=beta_nnL0pois^0.1, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="nonnegative L0 pen ident link Poisson GLM coefficients")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")
# graph2png(file="simulated_chromatogram_coefficients_nnL0pois.png", width=7, height=5)


# TODO: finish group penalized algo below
# (below just for nonnegative identity link Poisson)

lambda <- 1
beta <- beta_wnnls
thresh <- 1E-20
for (iter in 1:10) {
  beta[beta < thresh] <- 0
  beta_rowmax <- rowMax(abs(beta))
  nonzero <- as.vector(beta_rowmax>0) # variables that have nonzero coefficients in any of the channels
  lambdas <- lambda*1/beta_rowmax^2
  fitted <- sim$X[,nonzero,drop=F] %*% beta[nonzero,,drop=F]
  beta[nonzero,,drop=F] <- lmridge(X=X[,nonzero,drop=F],
                           y=y,
                           lambda=lambdas[nonzero],
                           weights = as.vector(1/(fitted+0.1)),
                           #lower = lower[nonzero],
                           #upper = upper[nonzero],
                           nonnegative = TRUE,
                           solver = "osqp",
                           method = NULL,
                           start = beta[nonzero],
                           control.solver = NULL,
                           verbose = TRUE)

  beta_nnL0pois <- do.call(rbind, lapply(which(nonzero),
                                         sim$y,
                                         2,
                                         function(ycol) {
    L0glm.fit2(X = sim$X,
               y = as.vector(ycol),
               lambda = 1,
               # weights = 1/(ycol+0.1),
               family = poisson(identity),
               lower = rep(0, ncol(sim$X)),
               upper = rep(Inf, ncol(sim$X)),
               nonnegative = TRUE,
               solver = "nnls",
               method = NULL,
               # start = rep(0, ncol(sim$X)), # prob faster to start from nnls coefficients?
               # control.solver = NULL,
               control.iwls = do.call("control.iwls.gen2", list(maxit = 1)),
               control.l0 = do.call("control.l0.gen2", list(maxit = 100, delta = 1E-10, gamma = 2)))$coefficients }))

  print(sum(beta>0))
}

plot(coefs_true, beta / X.n)


