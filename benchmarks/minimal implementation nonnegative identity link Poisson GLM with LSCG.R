devtools::load_all(".")

# minimal implementation of nonnegative L0 penalized identity link Poisson GLM
nnL0poisglm = function(X, 
                       y, 
                       lambda = 1, 
                       delta = 1E-10, 
                       gamma = 2,
                       tol = 1E-10,
                       thresh = 1) {  
  # print(sum(sim$beta_true!=0))
  
  weights = 1/(as.vector(y)+0.1) # observation weights (approx 1/Poisson noise variance)
  yw = y*sqrt(weights) # weighted y
  n = nrow(X)
  p = ncol(X)
  ywaug = c(yw, rep(0, p)) # augmented weighted y
  penweights = rep(1, p) # ridge regression penalty weights on individual coefficients (sum to p)
  lambdas = lambda*penweights # ridge penalties on individual coefficients
  Xw = Diagonal(x=sqrt(weights)) %*% X # = X * sqrt(weights)
  Xwaug = rbind(Xw, Diagonal(x=sqrt(lambdas))) # row augmented covariate matrix
  require(microbenchmark)
  beta_hat <- solve_sparse_LSCG(A = Xwaug,
                               b = ywaug,
                               x0 = as.vector(crossprod(Xwaug, ywaug) / colSums(Xwaug^2)), # initialisation from diagonal/jacobi preconditioner
                               tol = tol,
                               maxit = 0L,
                               precond = 0L)$coefficients # 0.06s
  # system.time(beta_hat <- lmridge(X = Xwaug, y = ywaug, lambdas=rep(0, ncol(X)),
  #                     pen.weights = rep(1, ncol(X)),
  #                     weights = NULL,
  #                     lower = rep(0, ncol(X)),
  #                     upper = rep(.Machine$double.xmax, ncol(X)),
  #                     nonnegative = TRUE,
  #                     solver = "eigen",
  #                     method = 4L,
  #                     start = as.vector(crossprod(Xwaug, ywaug) / colSums(Xwaug^2)),
  #                     control.solver = NULL)) # 2.25s, TODO: figure out why this is so slow compared to solve_sparse_LSCG call?
  hist(beta_hat)
  beta_hat[beta_hat<thresh] = 0

  for (i in 1:20) {  
    nonzero = as.vector(beta_hat!=0) # nonzero coefficients
    # print(sum(nonzero))
    penweights[nonzero] = 1/(1+abs(beta_hat[nonzero])^gamma+delta^gamma) # update L2 norm ridge regression penalty weights
    lambdas[nonzero] = lambda*penweights[nonzero] # updated ridge regression penalties lambdas
    Xwaug = rbind(Xw[,nonzero,drop=F], Diagonal(x=sqrt(lambdas[nonzero])))
    ywaug = c(yw, rep(0, sum(nonzero)))
    beta_hat[nonzero] = solve_sparse_LSCG(A = Xwaug,
                                          b = matrix(ywaug),
                                          x0 = as.vector(crossprod(Xwaug, ywaug) / colSums(Xwaug^2)), # initialisation from diagonal/jacobi preconditioner, beta_hat[nonzero] also works
                                          tol = 1E-10,
                                          maxit = 0L,
                                          precond = 1L)$coefficients
    beta_hat[beta_hat<thresh] = 0
  }
  return(beta_hat)
}

# example
set.seed(1)
sim = simulate_spike_train(n=10000, p=10000, k=200, sparse=TRUE) # blurred spike train with 10000 timepoints & 200 blurred peaks
library(microbenchmark)
system.time(coefs <- nnL0poisglm(X = sim$X,
            y = sim$y, 
            lambda = 1, 
            delta = 1E-10, 
            gamma = 2,
            tol = 1E-10,
            thresh = 1E-4)) # 0.33s
plot(sim$beta_true, coefs)
mean(abs(sim$beta_true - coefs)) # 1.493188 - really accurate!!
length(intersect(which(sim$beta_true>0), which(as.vector(coefs)>0))) # 198 (should be 200)

cor(as.vector(sim$beta_true), as.vector(coefs))