# TESTS L0 GROUP PENALIZED DECONVOLUTION OF MULTICHANNEL DATA

library(export)
library(qlcMatrix)

set.seed(1)
sim <- simulate_spike_train(n=200, p=200, k=20, multichannel=TRUE, sparse=TRUE)
X <- sim$X
Y <- sim$y
true_coefs <- sim$beta_true
m <- ncol(Y)
n <- nrow(X)
p <- ncol(X)

image(X, useRaster=TRUE)
# graph2png(file="simulated_chromatogram.png", width=7, height=5)

# test implementation of L0 group penalized nonnegative identity link Poisson GLM

# wnnls result fitted on each channel independently

# 1. using basic apply ####
system.time(beta_wnnls <- apply(sim$y, 2, function(ycol, index) {
  lmridge(X = sim$X,
          y = as.vector(ycol),
          lambda = rep(1, ncol(sim$X)),
          #weights = 1/(as.vector(ycol)+0.1),
          lower = rep(0, ncol(sim$X)),
          upper = rep(Inf, ncol(sim$X)),
          nonnegative = TRUE,
          solver = "eigen", # note: nnls or osqp solution much better than glmnet
          method = 7,
          start = rep(0, ncol(sim$X)),
          control.solver = NULL,
          verbose = TRUE) })) 

cor(as.vector(beta_wnnls), as.vector(sim$beta_true))
# 5s with nnls, 1.3s with osqp, 185s with glmnet, for n=p=500 & k=10
# 95s osqp for n=p=10000 & k=100
# TODO check why glmnet is so slow here and gives off warnings
image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=beta_wnnls^0.1, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="wnnls coefficients (nnls)")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim$X))[as.vector(rowMin(beta_wnnls)!=0)], col="cyan")
# graph2png(file="simulated_chromatogram_coefficients_wnnls.png", width=7, height=5)


system.time(beta_wnnls <- apply(sim$y, 2, function(ycol) {
  L0glm.fit(X = sim$X,
             y = as.vector(ycol),
             lambda = 1,
             #weights = 1 / (as.vector(ycol) + 0.1),
             family = poisson(identity),
             nonnegative = TRUE,
             solver = "eigen",
             method = 7,
             #start = rep(1, ncol(sim$X)),
             control.iwls = do.call("control.iwls.gen", list(maxit = 100)),
             control.l0 = do.call("control.l0.gen", list(maxit = 100, delta = 1E-10, gamma = 2)))$coefficients })) 

#  for n=p=500 & k=10
cor(as.vector(beta_wnnls), as.vector(true_coefs))
image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=beta_wnnls^0.1, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="wnnls coefficients (nnls)")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim$X))[as.vector(rowMin(beta_wnnls)!=0)], col="cyan")


# TODO check why glmnet is so slow here and gives off warnings
image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=beta_wnnls^0.1, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="wnnls coefficients (nnls)")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim$X))[as.vector(rowMin(beta_wnnls)!=0)], col="pink")
abline(v=which(colSums(beta_wnnls != 0) >= 5), col="pink")
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

# create the weighted block diagonal matrix
W <- sqrt(1/(yvect+0.1))
Xblock <- kronecker(diag(m), Xsp)
Xwblock <- Xweights(Xblock, W)

image(Xwblock[1:1080, 1:1080], useRaster = TRUE)

# define other variables
lambda <- 1
lower <- rep(0, ncol(Xwblock))
upper <- rep(Inf, ncol(Xwblock))
nonnegative <- TRUE
start <- NULL
verbose <- TRUE
X.n <- apply(Xblock, 2, norm, type = "2")

system.time(coefs_multi <- lmridge(X = Xwblock, y = yvect, lambda = lambda, lower = lower, upper = upper,
                                   nonnegative = nonnegative, solver = "osqp", method = 7)) 
#0.09s glmnet, 0.61s osqp, 82s nnls, 0.83 eigen.1, 0.15 eigen.4, (n=p=150, m=131) 
cor(coefs_multi, as.numeric(true_coefs))


system.time(fit_multi <- L0glm.fit(X=Xwblock, y=yvect, family = poisson(identity),
                                    solver = "eigen", method = 7, block = FALSE,
                                    lambda = 1, nonnegative = TRUE,
                                    control.iwls = do.call("control.iwls.gen", list(maxit = 100)),
                                    control.l0 = do.call("control.l0.gen", list(maxit = 100, delta = 1E-10, gamma = 2)),
                                    control.solver = do.call("control.solver.gen", list(solver = "eigen", method = 7, rowaug = TRUE))))
#767.3s osqp,  14.3s glmnet, 767.4s eigen.1 (bad result), 17.8s eigen.56, TOO LONG nnls, (n=p=150, m=131, 100 maxit)

# Turn coefficient vector into matrix
original_dim <- c(p, m)

coefficients_matrix <- matrix(fit_multi$coefficients, nrow = original_dim[1], ncol = original_dim[2], byrow = FALSE)
cor(fit_multi$coefficients, as.numeric(true_coefs)) #0.98 glmnet, 0.98 ospq, -0,002 eigen.1??, 0.977 eigen.56, 



image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=matrix(coefs_multi$coefficients)^0.1, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="Coefficients")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")












# nonnegative identity link Poisson GLM fitted on each channel independently
system.time(beta_nnpois <- apply(sim$y, 2, function(ycol) {
  L0glm.fit(X = sim$X,
             y = as.vector(ycol),
             lambda = 0,
             # weights = 1/(as.vector(ycol)+0.1),
             family = poisson(identity),
             nonnegative = TRUE,
             solver = "nnls",
             method = NULL,
             start = rep(0, ncol(sim$X)), # prob faster to start from nnls coefficients?
             control.solver = NULL,
             control.iwls = do.call("control.iwls.gen", list(maxit = 10)),
             control.l0 = do.call("control.l0.gen", list(maxit = 1, delta = 1E-10, gamma = 2)))$coefficients }))
# 10.7s with nnls, 4.3s with osqp for n=p=500 & k=10
image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=beta_nnpois^0.1, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="nonnegative ident link Poisson GLM coefficients (nnls)")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")
# graph2png(file="simulated_chromatogram_coefficients_nnpois.png", width=7, height=5)


# nonnegative L0 penalized identity link weighted LS fitted on each channel independently
system.time(beta_nnL0pois <- apply(sim$y, 2, function(ycol) {
  L0glm.fit(X = sim$X,
             y = as.vector(ycol),
             lambda = 1,
             # weights = 1/(as.vector(ycol)+0.1),
             family = poisson(identity),
             nonnegative = TRUE,
             solver = "osqp",
             method = NULL,
             # start = rep(0, ncol(sim$X)), # prob faster to start from nnls coefficients?
             # control.solver = NULL,
             control.iwls = do.call("control.iwls.gen", list(maxit = 1)),
             control.l0 = do.call("control.l0.gen", list(maxit = 10, delta = 1E-10, gamma = 2)))$coefficients }))
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
                                           L0glm.fit(X = sim$X,
                                                      y = as.vector(ycol),
                                                      lambda = 1,
                                                      # weights = 1/(ycol+0.1),
                                                      family = poisson(identity),
                                                      nonnegative = TRUE,
                                                      solver = "nnls",
                                                      method = NULL,
                                                      # start = rep(0, ncol(sim$X)), # prob faster to start from nnls coefficients?
                                                      # control.solver = NULL,
                                                      control.iwls = do.call("control.iwls.gen", list(maxit = 1)),
                                                      control.l0 = do.call("control.l0.gen", list(maxit = 100, delta = 1E-10, gamma = 2)))$coefficients }))
  
  print(sum(beta>0))
}

plot(coefs_true, beta / X.n)


# 4. using a for loop with avg of previous channels as priors ####
# Initialize a matrix to store the fitted coefficients for each channel
beta.matrix <- matrix(0, nrow = p, ncol = m)

# Iterate through the channels
for (channel in 1:m) {
  # Calculate the average of previously fitted coefficients as the starting weights
  if (channel > 1) {
    start <- rowMeans(beta.matrix[, 1:(channel - 1), drop = FALSE])
  } else {
    start <- rep(0, p)
  }
  
  # Fit the model for the current channel
  beta.matrix[, channel] <- lmridge(
    X = sim$X,
    y = as.vector(sim$y[, channel]),
    lambda = 1,
    start = start,
    lower = rep(0, p),
    upper = rep(Inf, p),
    nonnegative = TRUE,
    solver = "osqp",  
    method = NULL,
    control.solver = NULL,
    verbose = TRUE
  )
}
cor(as.vector(beta.matrix), as.vector(true_coefs))



# Create beta matrix
beta.matrix <- matrix(0, nrow = p, ncol = m)
weights.matrix <- matrix(0, nrow = p, ncol = m)
prior.matrix <- matrix(0, nrow = p, ncol = m)
fitted.matrix <- matrix(0, nrow = p, ncol = m)
eta.matrix <- matrix(0, nrow = p, ncol = m)
converged.iwls <- vector("numeric", length = 0)
iter.iwls <- vector("list", m)
iter.l0 <- vector("list", m)
lambda.matrix <- matrix(0, nrow = p, ncol = m)
lambda.w.matrix <- matrix(0, nrow = p, ncol = m)
converged.l0 <- vector("numeric", length = 0)

# Order columns of y decreasingly based on sparsity (decreasing)
# colnames(Y) <- seq(ncol(Y))
# y_dec <- Y[, order(colSums(Y != 0), decreasing = TRUE)] # Order by decreasing level of sparsity
# y_dec <- Y[, order(colMeans(Y), decreasing = TRUE)] # Order by decreasing column average (doesn't change results)

# Iterate through the channels
system.time(for (channel in 1:m) {
  # Calculate the average of previously fitted coefficients as the starting weights
  if (channel > 1) {
    start <- rowMeans(beta.matrix[, 1:(channel - 1), drop = FALSE])
  } else {
    start <- NULL
  }
  # Fit the model for the current channel
  current.fit <-  L0glm.fit(X = sim$X,
                     y = as.vector(sim$y[, channel]),
                     lambda = 1,
                     start = start,
                     family = poisson(identity),
                     nonnegative = TRUE,
                     solver = "glmnet",
                     method = 1,
                     control.iwls = do.call("control.iwls.gen", list(maxit = 1)),
                     control.l0 = do.call("control.l0.gen", list(warn = FALSE, maxit = 100, delta = 1E-10, gamma = 2)))
  beta.matrix[, channel] <- current.fit$coefficients
  weights.matrix[, channel] <- current.fit$weights
  prior.matrix[, channel] <- current.fit$prior.weights
  fitted.matrix[, channel] <- current.fit$fitted.values
  eta.matrix[, channel] <- as.vector(current.fit$eta)
  converged.iwls[channel] <- c(converged.iwls, current.fit$converged_iwls)
  iter.iwls[[channel]] <- current.fit$iter.iwls
  names(iter.iwls)[channel] <- paste("channel", channel)
  iter.l0[[channel]] <- current.fit$iter.l0
  names(iter.l0)[channel] <- paste("channel", channel)
  lambda.matrix <- current.fit$lambda
  lambda.w.matrix <- current.fit$lambda.w
  converged.l0 <- c(converged.l0, current.fit$converged.l0)
})
# 37s osqp, 5.67s glmnet, 2.35s eigen.1,  2.11s eigen.2, 1.96s eigen.5-6, 1s eigen.7 for n=p=150
# 8.25 glmnet (vs. 21s wo start), 6.77 eigen.5-6 (vs 25s), 135.2s eigen.1 (vs. 242,1), 130.9s eigen.2 (vs 232.1) for n=p=500
# 127.45s for eigen.5 n=p=10000

# cor(as.vector(beta.matrix), as.vector(true_coefs)) #will only give correct results if y cols haven't been reordered

# Assigning same colnames to beta.matrix as y_dec and reordering the true_coefs as well in order to allow for comparison
colnames(beta.matrix) <- colnames(y_dec) 

true_dec <- true_coefs[, order(colSums(sim$y != 0), decreasing = TRUE)]
beta.matrix <- beta.matrix[, order(as.numeric(colnames(beta.matrix)))]

cor(as.vector(beta.matrix), as.vector(true_coefs))


image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=beta.matrix^0.1, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="Coefficients")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")

abline(v=(1:nrow(sim$X))[as.vector(rowMax(beta.matrix)!=0)], col="pink")


sum(rowMax(sim$beta_true)!=0)
sum(rowMax(beta.matrix)!=0)

# 5. Implemented multichannel
system.time(L0glm_osqp1 <- L0glm.fit(X=sim$X, y=sim$y, family = poisson(identity),
                                      solver = "eigen", method = 5, lambda = 1, nonnegative = TRUE,
                                      control.iwls = do.call("control.iwls.gen", list(maxit = 200,reltol = 3)),
                                      control.l0 = do.call("control.l0.gen", list(maxit = 100, delta = 1E-10, gamma = 2)),
                                      verbose = TRUE))

names(L0glm_osqp1)
cor(as.vector(L0glm_osqp1$coefficients), as.vector(true_coefs))

image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=L0glm_osqp1$coefficients^0.1, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="Coefficients (averaging previously fitted channels)")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(L0glm_osqp1$coefficients)!=0)], col="yellow")

true <- rowMax(sim$beta_true)!=0
res <- rowMax(L0glm_osqp1$coefficients)!=0
sum(true == res)

#------------

system.time(group.cpglm.fit <- group.cpglm.fit(X = sim$X, Y = sim$y,
                             family = poisson(identity),
                             lambda = 95,
                             solver = "osqp",
                             method = 8,
                             nonnegative = TRUE))

# For p=n=10000 and m = 168
# L2 group penalty => osqp= 80.4s, glmnet = 78.13s, chol = 59.94, QR = too long, CG = 
# L1 group penalty => slower and worse results
cor(as.vector(group.cpglm.fit$coefficients), as.vector(true_coefs))

image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=group.cpglm.fit$coefficients^0.1, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="Coefficients (L2 penalty, chol pivot = 1)")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(group.cpglm.fit$coefficients)!=0)], col="yellow")

sum(rowMax(sim$beta_true)!=0) 
sum(rowMax(group.cpglm.fit$coefficients)!=0) 

#------------------------
system.time(group.cpglm <- group.cpglm(X = sim$X, Y = sim$y, 
                             family = poisson(identity),
                             lambda = 1, # TODO fix error OSQP lambda = 0
                             solver = "eigen",
                             method = 5,
                             nonnegative = TRUE,
                             verbose = FALSE))
# For p=n=10000 and m = 168
# osqp=113.44s (100 peaks), 71.7s (50 peaks), glmnet = 150.17s (100 peaks) / 95s (50 peaks), chol = 56.39s / 36s, QR(pivot 1) = too slow, CG = 55.82s
names(group.cpglm)
cor(as.vector(group.cpglm$coefficients), as.vector(true_coefs))

image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=group.cpglm$coefficients^0.1, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="Coefficients (averaging previously fitted channels, glmnet)")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(group.cpglm$coefficients)!=0)], col="cyan")

sum(rowMax(sim$beta_true)!=0)
sum(rowMax(group.cpglm$coefficients)!=0)
