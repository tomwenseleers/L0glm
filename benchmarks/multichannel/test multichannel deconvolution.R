# TESTS L0 GROUP PENALIZED DECONVOLUTION OF MULTICHANNEL DATA

library(export)
library(qlcMatrix)

sim = simulate_spike_train(n=500, p=500, k=10, multichannel=TRUE, sparse=TRUE)
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
                solver = "nnls", # note: nnls or osqp solution much better than glmnet
                method = NULL,
                start = rep(0, ncol(sim$X)),
                control.solver = NULL,
                verbose = TRUE) })) # 5s with nnls, 1.3s with osqp, 185s with glmnet, for n=p=500 & k=10
# TODO check why glmnet is so slow here
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



# 3. by working with a block diagonal matrix to make multichannel problem into single channel one: ####

# TODO finish this code - ask ChatGPT
# I think this can be done using the spam library & using a Kronecker product
library(spam)

yvect = as.vector(sim$y) # vectorized version of multichannel data

# create the weight matrix
W <- Diagonal(sqrt(1/(yvect+0.1)))

# create the block diagonal matrix
Xblock <- bKMatrix(W, X)


library(Matrix)
# create the weight matrix
W <- Diagonal(sqrt(1/(y+0.1)))

# create the block diagonal matrix
Xblock <- bdiag(lapply(1:nrow(Y), function(i) W %*% X))





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
             # start = rep(0, ncol(sim$X)), # prob faster to start from nnls coefficients?
             # control.solver = NULL,
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
                           lower = lower[nonzero],
                           upper = upper[nonzero],
                           nonnegative = TRUE,
                           solver = "nnls",
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
               control.l0 = do.call("control.l0.gen2", list(maxit = 100, delta = 1E-10, gamma = 2)))$coefficients })

  print(sum(beta>0))
}

plot(coefs_true, beta / X.n)


