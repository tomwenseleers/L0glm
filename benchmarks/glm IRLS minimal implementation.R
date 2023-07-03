# minimal plain glm IRLS implementation

glm_irls = function(X, y, weights=rep(1,nrow(X)), family=binomial(logit), maxit=25, tol=1e-16) {
  if (!is(family, "family")) family = family()
  variance = family$variance
  linkinv = family$linkinv
  mu.eta = family$mu.eta
  etastart = NULL
  
  nobs = nrow(X)    # needed by the initialize expression below
  nvars = ncol(X)   # needed by the initialize expression below
  eval(family$initialize) # initializes n and fitted values mustart
  eta = family$linkfun(mustart) # we then initialize eta with this
  dev.resids = family$dev.resids
  dev = sum(dev.resids(y, linkinv(eta), weights))
  devold = 0
  beta_old = rep(1, nvars)
  
  for(j in 1:maxit)
  {
    mu = linkinv(eta) 
    varg = variance(mu)
    gprime = mu.eta(eta)
    z = eta + (y - mu) / gprime # potentially -offset if you would have an offset argument as well
    W = weights * as.vector(gprime^2 / varg)
    beta = solve(crossprod(X,W*X), crossprod(X,W*z), tol=2*.Machine$double.eps)
    eta = X %*% beta # potentially +offset if you would have an offset argument as well
    dev = sum(dev.resids(y, mu, weights))
    if (abs(dev - devold) / (0.1 + abs(dev)) < tol) break
    devold = dev
    beta_old = beta
  }
  list(coefficients=t(beta), iterations=j)
}