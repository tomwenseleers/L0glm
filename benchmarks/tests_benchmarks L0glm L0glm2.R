# Simulate data
set.seed(1)
n <- 100
p <- 100
s <- 0.1 # sparsity as proportion of p that have nonzero coefficients // p*s == #of peaks
k <- round(p*s) # nr of nonzero covariates
sim <- simulate_spike_train(n=n, p=p, k=k, mean_beta = 1000, sd_logbeta = 1, family="poisson")
# sim <- simulate_spike_train()
coefs_true <- as.vector(sim$beta_true)

Xsp <- sim$X
X.n <- apply(Xsp, 2, norm, type = "2")
Xsp <- sweep(Xsp, 2, X.n, "/")
X <- as.matrix(Xsp)
y <- as.vector(sim$y)
lambda <- 1
pen_weights <- rep(1, ncol(X))
lambdas <- lambda*pen_weights
lower <- rep(5, ncol(X))
upper <- rep(Inf, ncol(X))
nonnegative <- TRUE
start <- NULL
verbose <- TRUE
weights <- 1/(y+0.1) # observation weights

control.solver = do.call("control.solver.gen", list(solver = "eigen", tol=1e-3, method = 3))
lmridge.eigen(Xw, yw, pen.weights, lambdas, lower, upper, start, control.solver = control.solver)


# 1. OSQP
system.time(L0glm_osqp <- L0glm.fit(X=Xsp, 
                                     y=y, 
                                     family = poisson(identity),
                                     solver = "osqp", 
                                     lambda = 1, 
                                     nonnegative = TRUE, 
                                     control.iwls = do.call("control.iwls.gen", list(maxit = 100)),
                                     control.l0 = do.call("control.l0.gen", list(maxit = 100, delta = 1E-10, gamma = 2)),
                                     control.solver = do.call("control.solver.gen", list(solver = "osqp", rowaug = FALSE)),
                                     control.block = do.call("control.block.gen", list(maxit = 70, block.size = 100))
                                     )) # 35s for n=p=10000 lambda = 1


plot_L0glm_benchmark(x = sim$x, y = y,
                     coefs = L0glm_osqp$coefficients / X.n,
                     fitted = L0glm_osqp$fitted.values,
                     a.true = sim$beta_true,
                     main="Ground truth vs l0glm estimates")

cor(coefs_true, L0glm_osqp$coefficients / X.n) # 0.9990666

sum(coefs_true>0) # 10
sum(L0glm_osqp$coefficients>0) # 10 with lambda = 1

# 2. GLMNET
system.time(L0glm_glmnet <- L0glm.fit(X=Xsp, 
                                       y=y, 
                                       family = poisson(identity),
                                       solver = "glmnet",
                                       lambda = 1,
                                       nonnegative = T,
                                       control.iwls = do.call("control.iwls.gen", list(maxit = 100)),
                                       control.l0 = do.call("control.l0.gen", list(maxit = 100, delta = 1E-10, gamma = 2)),
                                       control.solver = do.call("control.solver.gen", list(solver = "glmnet", rowaug = TRUE)),
                                       #control.block = do.call("control.block.gen", list(maxit = 50, block.size = 15))
                                      )) # 46s for n=p=10000 lambda = 1


plot_L0glm_benchmark(x = sim$x, y = y,
                     coefs = L0glm_glmnet$coefficients / X.n,
                     fitted = L0glm_glmnet$fitted.values,
                     a.true = sim$beta_true,
                     main="Ground truth vs l0glm estimates")
cor(coefs_true, L0glm_glmnet$coefficients / X.n) # 0.99972
sum(coefs_true>0) # 10
sum(L0glm_glmnet$coefficients>0) # 21



# 3. Eigen - QR
system.time(L0glm_qr1 <- L0glm.fit(X=Xsp,
                                  y=y, 
                                  family = poisson(identity), 
                                  verbose = TRUE,
                                  solver = "eigen",
                                  method = 1,
                                  lambda = 1,
                                  start = rep(0.1, ncol(X)),
                                  nonnegative = T,
                                  control.iwls = do.call("control.iwls.gen", list(maxit = 100)),
                                  control.l0 = do.call("control.l0.gen", list(maxit = 100, delta = 1E-10, gamma = 1.5)),
                                  control.solver = do.call("control.solver.gen", list(solver = "eigen", method = 1, rowaug = F)),
                                  #control.block = do.call("control.block.gen", list(maxit = 50, block.size = 15)),
                                  )) # 46s for n=p=10000 lambda = 1 

plot_L0glm_benchmark(x = sim$x, y = y,
                     coefs = L0glm_qr1$coefficients / X.n,
                     fitted = L0glm_qr1$fitted.values,
                     a.true = sim$beta_true,
                     main="Ground truth vs l0glm estimates")
cor(coefs_true, L0glm_qr1$coefficients / X.n) # 0.9253677
sum(coefs_true>0) # 10
sum(L0glm_qr1$coefficients>0) # 10


#l0glm Eigen - LSCG

system.time(L0glm_lscg0 <- L0glm.fit(X=Xsp, y=y, family = poisson(identity),
                                    solver = "eigen", method = 3,
                                    lambda = 1, nonnegative = TRUE,
                                    verbose = TRUE,
                                    control.iwls = do.call("control.iwls.gen", list(maxit = 0)),
                                    control.l0 = do.call("control.l0.gen", list(maxit = 100, delta = 1E-10, gamma = 1.5)),
                                    control.solver = do.call("control.solver.gen", list(solver = "eigen", method = 3, rowaug = FALSE))))
# TODO gives bad results because it produces nonzero coefficients, CHECK ERROR replacement has length zero
plot_L0glm_benchmark(x = sim$x, y = y,
                     coefs = L0glm_lscg0$coefficients / X.n,
                     fitted = L0glm_lscg0$fitted.values,
                     a.true = sim$beta_true,
                     main="Ground truth vs l0glm estimates")
cor(coefs_true, L0glm_lscg0$coefficients / X.n) #
sum(coefs_true>0) # 10
sum(L0glm_lscg0$coefficients>0) #


system.time(L0glm_lscg1 <- L0glm.fit(X=Xsp, y=y, family = poisson(identity),
                                    solver = "eigen", method = 4,
                                    lambda = 1, nonnegative = TRUE,
                                    control.iwls = do.call("control.iwls.gen", list(maxit = 100)),
                                    control.l0 = do.call("control.l0.gen", list(maxit = 0, delta = 1E-10, gamma = 2)),
                                    control.solver = do.call("control.solver.gen", list(solver = "eigen", method = 4, rowaug = FALSE))))
# TODO FIX gives bad results bc produces nonzero coefficients ?
plot_L0glm_benchmark(x = sim$x, y = y,
                     coefs = L0glm_lscg1$coefficients / X.n,
                     fitted = L0glm_lscg1$fitted.values,
                     a.true = sim$beta_true,
                     main="Ground truth vs l0glm estimates")
cor(coefs_true, L0glm_lscg1$coefficients / X.n) # 0.7240002 NOT OK
sum(coefs_true>0) # 100
sum(L0glm_lscg1$coefficients>0) # 360



#l0glm Eigen - CG

system.time(L0glm_cg0 <- L0glm.fit(X=Xsp, y=y, family = poisson(identity),
                                      solver = "eigen", method = 5, # or 6?
                                      lambda = 1, nonnegative = TRUE,
                                      control.iwls = do.call("control.iwls.gen", list(maxit = 200)),
                                      control.l0 = do.call("control.l0.gen", list(maxit = 100, delta = 1E-10, gamma = 1.5))
)) 
plot_L0glm_benchmark(x = sim$x, y = y,
                     coefs = L0glm_cg0$coefficients / X.n,
                     fitted = L0glm_cg0$fitted.values,
                     a.true = sim$beta_true,
                     main="Ground truth vs l0glm estimates")
cor(coefs_true, L0glm_cg0$coefficients / X.n) #
sum(coefs_true>0) # 100
sum(L0glm_cg0$coefficients>0) #


system.time(L0glm_cg1 <- L0glm.fit(X=Xsp, y=y, family = poisson(identity),
                                      solver = "eigen", method = 6,
                                      lambda = 1, nonnegative = TRUE,
                                      control.iwls = do.call("control.iwls.gen", list(maxit = 200)),
                                      control.l0 = do.call("control.l0.gen", list(maxit = 1000, delta = 1E-10, gamma = 1.5))
)) 
plot_L0glm_benchmark(x = sim$x, y = y,
                     coefs = L0glm_cg1$coefficients / X.n,
                     fitted = L0glm_cg1$fitted.values,
                     a.true = sim$beta_true,
                     main="Ground truth vs l0glm estimates")
cor(coefs_true, L0glm_cg1$coefficients / X.n) # 0.7240002 NOT OK
sum(coefs_true>0) # 100
sum(L0glm_cg1$coefficients>0) # 360


#l0glm Eigen - CHOL - best results so far
system.time(L0glm_chol0 <- L0glm.fit(X=Xsp, y=y, family = poisson(identity),
                                    solver = "eigen", method = 8, upper = 2000,# or 8?
                                    lambda = 1, nonnegative = TRUE, 
                                    control.iwls = do.call("control.iwls.gen", list(maxit = 1000)),
                                    control.l0 = do.call("control.l0.gen", list(maxit = 100, delta = 1E-10, gamma = 2)),
                                    control.solver = do.call("control.solver.gen", list(solver = "eigen", method = 8, rowaug = FALSE))))

#19.14s for n=p=10000 (for sparse matrix)
plot_L0glm_benchmark(x = sim$x, y = sim$y,
                     coefs = L0glm_chol0$coefficients / X.n,
                     fitted = L0glm_chol0$fitted.values,
                     a.true = sim$beta_true,
                     main="Ground truth vs l0glm estimates")

cor(coefs_true, L0glm_chol0$coefficients) # 0.95 for n=p=10000
sum(coefs_true>0) # 1000
sum(L0glm_chol0$coefficients>0) # 1067

