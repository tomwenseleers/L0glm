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
  pen.weights = rep(1, p) # ridge regression penalty weights on individual coefficients (sum to p)
  lambdas = lambda*pen.weights # ridge penalties on individual coefficients
  Xw = Diagonal(x=sqrt(weights)) %*% X # = X * sqrt(weights)
  Xwaug = rbind(Xw, Diagonal(x=sqrt(lambdas))) # row augmented covariate matrix
  # require(microbenchmark)
  beta_hat <- lmridge.eigen(Xw = Xwaug, yw = ywaug,
                                        pen.weights = NULL, # PS remove this option?
                # lambdas = lambdas,
                lambdas = rep(0, ncol(X)),
                start = as.vector(crossprod(Xwaug, ywaug) / colSums(Xwaug^2)),
                lower = rep(0, ncol(X)), upper = rep(Inf, ncol(X)),
                control.solver = do.call("control.solver.gen",
                                         list(solver="eigen", method=4L))) # 0.06s
  # range(beta_hat)

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
    
    # beta_hat2=beta_hat
    beta_hat[nonzero] <- lmridge.eigen(Xw = Xw[,nonzero,drop=F], yw = yw,
                              pen.weights = NULL, # PS remove this option?
                              # lambdas = lambdas,
                              lambdas = lambdas[nonzero], # rep(0, ncol(X)),
                              start = NULL, # as.vector(crossprod(Xwaug, ywaug) / colSums(Xwaug^2)),
                              lower = rep(0, sum(nonzero)), upper = rep(Inf, sum(nonzero)),
                              control.solver = do.call("control.solver.gen",
                                                       list(solver="eigen", method=7L))) # 0.01s
    # beta_hat[nonzero] <- solve_sparse_LSCG(A = Xwaug,
    #                                       b = matrix(ywaug),
    #                                       x0 = as.vector(crossprod(Xwaug, ywaug) / colSums(Xwaug^2)), # initialisation from diagonal/jacobi preconditioner, beta_hat[nonzero] also works
    #                                       tol = 1E-10,
    #                                       maxit = 0L,
    #                                       precond = 1L)$coefficients
    # XwaugTXwaug = as(crossprod(Xwaug), "generalMatrix")
    # # X.n = sqrt(diag(XwaugTXwaug))
    # beta_hat[nonzero] <- solve_sparse_chol(A = XwaugTXwaug, # / X.n
    #                                        b = as.matrix(crossprod(Xwaug, ywaug)),
    #                                        pivot = 1L,
    #                                        ordering = 0L)$coefficients # / X.n
    beta_hat[beta_hat<thresh] = 0
    # print(range(beta_hat))
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
library(microbenchmark)
system.time(coefs <- nnL0poisglm(X = sim$X,
            y = sim$y, 
            lambda = 1, 
            delta = 1E-10, 
            gamma = 2,
            tol = 1E-10,
            thresh = 1)) # 0.2s
plot(sim$beta_true, coefs)
mean(abs(sim$beta_true - coefs)) # 37.53146 - really accurate!!

length(intersect(which(sim$beta_true>0), which(as.vector(coefs)>0))) # 859 (should be 1000)
sum(coefs>0) # 1079 (should be 1000)

