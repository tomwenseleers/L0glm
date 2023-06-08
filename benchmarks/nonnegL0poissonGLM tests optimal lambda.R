
nonnegL0poisglm = function(X = sim$X, 
                           y = sim$y, 
                           lambda = 1, 
                           delta = 1E-10, 
                           gamma = 2) {  
  print(sum(sim$beta_true!=0))
  
  weights = 1/(y+0.1) # observation weights (approx 1/Poisson noise variance)
  yw = y*sqrt(weights) # weighted y
  n = nrow(X)
  p = ncol(X)
  ywaug = c(yw, rep(0, p)) # augmented weighted y
  penweights = rep(1, p) # ridge regression penalty weights on individual coefficients (sum to p)
  lambdas = lambda*penweights # ridge penalties on individual coefficients
  Xw = X*Diagonal(x=sqrt(weights)) # weighted X
  Xwaug = rbind(Xw, Diagonal(x=sqrt(lambdas))) # row augmented covariate matrix
  P = Diagonal(x=1/colSums(Xwaug^2)) # Jacobi/diagonal type preconditioning matrix
  beta_hat = P %*% crossprod(Xwaug, ywaug) # coefficient initialisation based on preconditioning matrix
  beta_hat[beta_hat<0] = 0
  # these coefficients will be bounded between 0 and 1/lambda
  hist(as.vector(beta_hat))
  
  for (i in 1:10) {  
    nonzero = as.vector(beta_hat!=0) # nonzero coefficients
    print(sum(nonzero))
    penweights[nonzero] = 1/(1+abs(beta_hat[nonzero])^gamma+delta^gamma) # update L2 norm ridge regression penalty weights
    lambdas[nonzero] = lambda*penweights[nonzero] # updated ridge regression penalties lambdas
    # Xwaug = rbind(Xw[,nonzero,drop=F], Diagonal(x=sqrt(lambdas[nonzero])))  
    # y waug = c(yw, rep(0, sum(nonzero)))
    # P = Diagonal(x=1/colSums(Xwaug^2)) # Jacobi/diagonal type preconditioning matrix
    beta_hat[nonzero] = lmridge(X=Xw[,nonzero,drop=F],
                                y=yw,
                                lambdas=lambdas[nonzero],
                                nonnegative=TRUE,
                                solver="osqp") # P %*% crossprod(Xwaug, ywaug) # updated coefficients
    hist(as.vector(beta_hat))
  }
  
  plot(sim$beta_true, as.vector(beta_hat))  
  cor(as.vector(sim$beta_true), as.vector(beta_hat)) # 0.47   
  
}


