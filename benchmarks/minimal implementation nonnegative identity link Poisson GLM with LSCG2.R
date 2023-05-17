devtools::load_all(".")

# minimal implementation of nonnegative L0 penalized identity link Poisson GLM
L0glm.fit_new = function(X, 
                         y, 
                         weights = rep(1, ncol(X)),
                         family = poisson(identity),
                         offset = rep(0, nrow(X)),
                         start = NULL,
                         lower = rep(0, ncol(X)),
                         upper = rep(Inf, ncol(X)),
                         lambda = 1, 
                         delta = 1E-10, 
                         gamma = 2,
                         thresh = 1,
                         solver = "eigen", # or "osqp" or "glmnet"
                         method = 7L, # = Cholesky # 3(+4) = LSCG, 5(+6) = CG, 7(+8) = Cholesky  
                         maxit = 25,
                         tol = 1E-5,
                         warmupiters = 4) { # nr of iterations before obs weights are updated in addition to penalty weights in each iteration, 3 or 4 gave best performance
  
  if (!is(family, "family")) family = family()
  variance = family$variance # variance function for given family
  linkinv = family$linkinv # inverse link function
  if (!is.function(variance) || !is.function(linkinv) )
    stop("'family' argument seems not to be a valid family object", call. = FALSE)
  mu.eta = family$mu.eta # function of derivative of the expected values
  etastart = NULL
  
  nobs = nrow(X)    # needed by the initialize expression below
  nvars = ncol(X)   # needed by the initialize expression below
  
  # if (nvars == 0) { # In case a NULL model is fit, code taken from 'glm'
  #   return(list(coefficients = numeric(),
  #               weights = numeric(),
  #               prior.weights = weights,
  #               fitted.values = rep(0, n),
  #               family = family,
  #               converged.iwls = TRUE,
  #               iter.iwls = 0))
  # }
  
  eval(family$initialize) # initializes n and fitted values mustart based on y (mu=expected value on response scale)
  
  if (is.null(start)) {  beta <- beta_old <- rep(1, nvars) # initial model coefficients
                         eta = family$linkfun(mustart) } else {  # eta = expected value on link scale
                           beta <- beta_old <- start # when starting coefficient values were provided
                           eta = offset + as.vector(if (ncol(X) == 1L) X * start else X %*% start) # eta = expected value on link scale
                         }
  names(beta) <- colnames(X)
  
  mu = linkinv(eta) # expected value on response scale
  dev.resids = family$dev.resids # deviance residuals function
  dev = sum(dev.resids(y, mu, weights)) # deviance
  dev_old = 0
  nonzero = as.vector(beta!=0) # nonzero coefficients

  for (i in 1L:maxit) {
    print(i)
    varg = variance(mu) # variance
    gprime = mu.eta(eta) # derivative of the expected values

    good = (weights > 0) & (gprime != 0) # valid rows
    
    if (i==1|i>warmupiters) {
      if (family$link == "identity"){ # simplified formula and faster than formula in else statement
        W = as.vector(weights/variance(mu)) # 1/variance (* weights) with identity link, where variance = 1 with Gaussian noise, variance = mu with Poisson noise, variance = mu^2 with Gamma noise
        good = (W > 0) & (gprime != 0) # valid rows
        z = (y - offset)[good] # response
      } else {
        mu.prime <- family$mu.eta(eta) # derivative of the expected values
        good = (weights > 0) & (mu.prime != 0) # valid rows
        z <- (eta - offset)[good] + (y - mu)[good] / mu.prime[good] # working (adjusted) response
        W <- weights[good] * as.vector(mu.prime[good]^2 / variance(mu)[good]) # updated weights 
      }}
    
    if (all(!good)) {
      conv <- FALSE
      warning(gettextf("no observations informative at iteration %d",
                       i), domain = NA)
      break
    }
      
    pen.weights = AR.weights(beta = beta, delta = delta, gamma = gamma) # adaptive ridge pentalty weights
    lambdas = lambda*pen.weights # ridge penalties on individual coefficients
    
    # weighted least squares ridge regression
    # Xw = Diagonal(x=sqrt(W[good])) %*% X[good,,drop=F] # = X * sqrt(W)
    # Xwaug = rbind(Xw, Diagonal(x=sqrt(lambdas))) # row augmented covariate matrix
    # require(microbenchmark)
    # beta = solve(crossprod(Xw), crossprod(Xw,z*sqrt(W)), tol=2*.Machine$double.eps)
    beta_old = beta
    beta[nonzero] <- lmridge(X = X[good, nonzero, drop=F], y = z,
                   weights = W,
                   pen.weights = pen.weights[nonzero], # TODO: remove this argument
                   lambdas = lambdas[nonzero],
                   start = NULL, # as.vector(crossprod(Xwaug, ywaug) / colSums(Xwaug^2)),
                   lower = lower[nonzero], 
                   upper = upper[nonzero],
                   solver = solver,
                   method = method)
    beta[nonzero][abs(beta[nonzero])<thresh] = 0
    nonzero = as.vector(beta!=0) # nonzero coefficients
    if (any(!is.finite(beta))) {
      conv <- FALSE
      warning(gettextf("non-finite coefficients at iteration %d", i), domain = NA)
      break
    }
    ## calculate updated values of eta and mu with the new coefs:
    eta = X[,nonzero,drop=F] %*% beta[nonzero] + offset # expected value / predictions on link scale
    mu = linkinv(eta) # expected value / predictions on response scale
    dev_old = dev
    dev = sum(dev.resids(y, mu, weights))
    
    if (abs(dev - dev_old) / (0.1 + abs(dev)) < tol) break
    # or if (norm(beta-beta_old, "2") < tol) break
    dev_old = dev
    beta_old = beta
  }
  
  return(beta)
  
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


L0glm.fit_new(X = X, y = y)


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
