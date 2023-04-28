devtools::load_all()

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
lower <- rep(0, ncol(X))
upper <- rep(Inf, ncol(X))
nonnegative <- TRUE
start <- NULL
verbose <- TRUE
weights <- 1/(y+0.1) # observation weights




# 1. L0glm ##
system.time(L0glm1.fit0 <- L0glm.fit(X=X, y=y, family = poisson(identity), weights = weights,
                                    lambda = 0, nonnegative = TRUE,
                                    control.iwls = list(maxit = 2, warn = FALSE),
                                    control.l0 = list(maxit = 1, rel.tol = 1E-4,
                                                      delta = 1E-10, gamma = 2, warn = FALSE)
)) # 948s for n=p=10000  # TODO make sure that function also works if you just give subset of options
# TODO make sure that it works for delta=0 (check nonzero subsetting)

# system.time(L0glm1 <- L0glm(formula = "y ~ 0 + .", data = data.frame(X, y = y), family = poisson(identity),
#                 lambda = 0, tune.meth = "none", nonnegative = TRUE, normalize = FALSE,
#                 control.iwls = list(maxit = 1),
#                 control.l0 = list(maxit = 10, delta = 1E-10, gamma = 2))) #Time = 7.2s

plot_L0glm_benchmark(x = sim$x, y = y,
                     coefs = L0glm1.fit0$coefficients / X.n,
                     fitted = L0glm1.fit0$fitted.values, a.true = sim$beta_true,
                     main="Ground truth vs L0glm estimates")

cor(coefs_true, L0glm1.fit0$coefficients / X.n) # 0.9634462

sum(coefs_true>0) # 1000
sum(L0glm1.fit0$coefficients>0) # 2022



system.time(L0glm1.fit1 <- L0glm.fit(X=X, y=y, family = poisson(identity),
                            lambda = 1, nonnegative = TRUE,
                            control.iwls = list(maxit = 1, warn = FALSE, thresh=1E-5),
                            control.l0 = list(maxit = 10, rel.tol = 1E-4,
                                              delta = 1E-30, gamma = 2, warn = FALSE)
                            )) #
# TO DO make sure that function also works if you just give subset of options
# system.time(L0glm1 <- L0glm(formula = "y ~ 0 + .", data = data.frame(X, y = y), family = poisson(identity),
#                 lambda = 1, tune.meth = "none", nonnegative = TRUE, normalize = FALSE,
#                 control.iwls = list(maxit = 1),
#                 control.l0 = list(maxit = 10, delta = 1E-10, gamma = 2))) #Time = 7.2s

plot_L0glm_benchmark(x = sim$x, y = y,
                     coefs = L0glm1.fit1$coefficients / X.n,
                     fitted = L0glm1.fit1$fitted.values, a.true = sim$beta_true,
                     main="Ground truth vs L0glm estimates")

cor(coefs_true, L0glm1.fit1$coefficients / X.n) # 0.9999708

sum(coefs_true>0) # 100
sum(L0glm1.fit1$coefficients>0) # 142



# 2. L0glm2 OSQP
system.time(L0glm_osqp0 <- L0glm.fit2(X=Xsp, y=y, family = poisson(identity),
                 solver = "osqp", lambda = 0, nonnegative = TRUE,
                 control.iwls = do.call("control.iwls.gen2", list(maxit = 100)),
                 control.l0 = do.call("control.l0.gen2", list(maxit = 100, delta = 1E-10, gamma = 2)),
                 control.solver = do.call("control.solver.gen", list(solver = "osqp", rowaug = FALSE))
            )) # 7.2s for n=p=10000
plot_L0glm_benchmark(x = sim$x, y = y,
                     coefs = L0glm_osqp0$coefficients / X.n,
                     fitted = L0glm_osqp0$fitted.values,
                     a.true = sim$beta_true,
                     main="Ground truth vs L0glm2 estimates")

cor(coefs_true, L0glm_osqp0$coefficients / X.n) # 0.9990666

sum(coefs_true>0) # 10
sum(L0glm_osqp0$coefficients>0) # 62


system.time(L0glm_osqp1 <- L0glm.fit2(X=Xsp, y=y, family = poisson(identity),
                                      solver = "osqp", lambda = 1, nonnegative = TRUE, verbose = TRUE,
                                      control.iwls = do.call("control.iwls.gen2", list(maxit = 10000, reltol  = 1E-3)),
                                      control.l0 = do.call("control.l0.gen2", list(maxit = 100, delta = 1E-10, gamma = 2)),
                                      control.solver = do.call("control.solver.gen", list(solver = "osqp", rowaug = TRUE))))

plot_L0glm_benchmark(x = sim$x, y = y,
                     coefs = L0glm_osqp1$coefficients / X.n,
                     fitted = L0glm_osqp1$fitted.values,
                     a.true = sim$beta_true,
                     main="Ground truth vs L0glm2 estimates")

cor(coefs_true, L0glm_osqp1$coefficients / X.n) # 0.938751
sum(coefs_true>0) # 10
sum(L0glm_osqp1$coefficients>0) # 18

cor(L0glm_osqp0$coefficients, L0glm_osqp1$coefficients)

#L0glm2 GLMNET
system.time(L0glm_glmnet0 <- L0glm.fit2(X=Xsp, y=y, family = poisson(identity), 
                                        solver = "glmnet", lambda = 0, nonnegative = TRUE,
                                        control.iwls = do.call("control.iwls.gen2", list(maxit = 1)),
                                        control.l0 = do.call("control.l0.gen2", list(maxit = 100, delta = 1E-10, gamma = 2))))



plot_L0glm_benchmark(x = sim$x, y = y,
                     coefs = L0glm_glmnet0$coefficients / X.n,
                     fitted = L0glm_glmnet0$fitted.values,
                     a.true = sim$beta_true,
                     main="Ground truth vs L0glm2 estimates")
cor(coefs_true, L0glm_glmnet0$coefficients / X.n) # 0.9990665
sum(coefs_true>0) # 10
sum(L0glm_glmnet0$coefficients>0) # 21


system.time(L0glm_glmnet1 <- L0glm.fit2(X=Xsp, y=y, family = poisson(identity), 
                                        solver = "glmnet", lambda = 1, nonnegative = TRUE,
                                        control.iwls = do.call("control.iwls.gen2", list(maxit = 100, warn = F)),
                                        control.l0 = do.call("control.l0.gen2", list(maxit = 100, delta = 1E-10, gamma = 2))
))

plot_L0glm_benchmark(x = sim$x, y = y,
                     coefs = L0glm_glmnet1$coefficients / X.n,
                     fitted = L0glm_glmnet1$fitted.values,
                     a.true = sim$beta_true,
                     main="Ground truth vs L0glm2 estimates")

cor(coefs_true, L0glm_glmnet1$coefficients / X.n) # 0.9990138
cor(L0glm_glmnet0$coefficients / X.n, L0glm_glmnet1$coefficients / X.n) # 0.9999587
sum(coefs_true>0) # 100
sum(L0glm_glmnet1$coefficients>0) # 111



#L0glm2 Eigen - QR

system.time(L0glm_qr1 <- L0glm.fit2(X=Xsp, y=y, family = poisson(identity), verbose = TRUE,
                                      solver = "eigen", lambda = 1, nonnegative = TRUE,
                                      control.iwls = do.call("control.iwls.gen2", list(maxit = 100)),
                                      control.l0 = do.call("control.l0.gen2", list(maxit = 100, delta = 1E-10, gamma = 1.5)),
                                      control.solver = do.call("control.solver.gen", list(solver = "eigen", method = 1, rowaug = TRUE)))) # 0.50s 

# TODO make sure that all solvers always give a result even for underdetermined cases and lambda=0 -> then set to smallest lambda value possible
plot_L0glm_benchmark(x = sim$x, y = y,
                     coefs = L0glm_qr0$coefficients / X.n,
                     fitted = L0glm_qr0$fitted.values,
                     a.true = sim$beta_true,
                     main="Ground truth vs L0glm2 estimates")
cor(coefs_true, L0glm_qr0$coefficients / X.n) # 0.9253677
sum(coefs_true>0) # 10
sum(L0glm_qr0$coefficients>0) # 15


system.time(L0glm_qr2 <- L0glm.fit2(X=Xsp, y=y, family = poisson(identity),
                                    solver = "eigen", lambda = 1, nonnegative = TRUE,
                                    control.iwls = do.call("control.iwls.gen2", list(maxit = 100)),
                                    control.l0 = do.call("control.l0.gen2", list(maxit = 100, delta = 1E-10, gamma = 1.5)),
                                    control.solver = do.call("control.solver.gen", list(solver = "eigen", method = 2, rowaug = TRUE))))

plot_L0glm_benchmark(x = sim$x, y = y,
                     coefs = L0glm_qr1$coefficients / X.n,
                     fitted = L0glm_qr1$fitted.values,
                     a.true = sim$beta_true,
                     main="Ground truth vs L0glm2 estimates")
cor(coefs_true, L0glm_qr1$coefficients / X.n) # 0.9879139
sum(coefs_true>0) # 10
sum(L0glm_qr1$coefficients>0) # 16


#L0glm2 Eigen - LSCG

system.time(L0glm_lscg0 <- L0glm.fit2(X=Xsp, y=y, family = poisson(identity),
                                    solver = "eigen", method = 3,
                                    lambda = 1, nonnegative = TRUE,
                                    verbose = TRUE,
                                    control.iwls = do.call("control.iwls.gen2", list(maxit = 0)),
                                    control.l0 = do.call("control.l0.gen2", list(maxit = 100, delta = 1E-10, gamma = 1.5)),
                                    control.solver = do.call("control.solver.gen", list(solver = "eigen", method = 3, rowaug = FALSE))))
# TODO gives bad results because it produces nonzero coefficients, CHECK ERROR replacement has length zzero
plot_L0glm_benchmark(x = sim$x, y = y,
                     coefs = L0glm_lscg0$coefficients / X.n,
                     fitted = L0glm_lscg0$fitted.values,
                     a.true = sim$beta_true,
                     main="Ground truth vs L0glm2 estimates")
cor(coefs_true, L0glm_lscg0$coefficients / X.n) #
sum(coefs_true>0) # 10
sum(L0glm_lscg0$coefficients>0) #


system.time(L0glm_lscg1 <- L0glm.fit2(X=Xsp, y=y, family = poisson(identity),
                                    solver = "eigen", method = 4,
                                    lambda = 1, nonnegative = TRUE,
                                    control.iwls = do.call("control.iwls.gen2", list(maxit = 100)),
                                    control.l0 = do.call("control.l0.gen2", list(maxit = 0, delta = 1E-10, gamma = 2)),
                                    control.solver = do.call("control.solver.gen", list(solver = "eigen", method = 4, rowaug = FALSE))))
# TODO FIX gives bad results bc produces nonzero coefficients ?
plot_L0glm_benchmark(x = sim$x, y = y,
                     coefs = L0glm_lscg1$coefficients / X.n,
                     fitted = L0glm_lscg1$fitted.values,
                     a.true = sim$beta_true,
                     main="Ground truth vs L0glm2 estimates")
cor(coefs_true, L0glm_lscg1$coefficients / X.n) # 0.7240002 NOT OK
sum(coefs_true>0) # 100
sum(L0glm_lscg1$coefficients>0) # 360



#L0glm2 Eigen - CG

system.time(L0glm_cg0 <- L0glm.fit2(X=Xsp, y=y, family = poisson(identity),
                                      solver = "eigen", method = 5, # or 6?
                                      lambda = 1, nonnegative = TRUE,
                                      control.iwls = do.call("control.iwls.gen2", list(maxit = 200)),
                                      control.l0 = do.call("control.l0.gen2", list(maxit = 100, delta = 1E-10, gamma = 1.5))
)) 
plot_L0glm_benchmark(x = sim$x, y = y,
                     coefs = L0glm_cg0$coefficients / X.n,
                     fitted = L0glm_cg0$fitted.values,
                     a.true = sim$beta_true,
                     main="Ground truth vs L0glm2 estimates")
cor(coefs_true, L0glm_cg0$coefficients / X.n) #
sum(coefs_true>0) # 100
sum(L0glm_cg0$coefficients>0) #


system.time(L0glm_cg1 <- L0glm.fit2(X=Xsp, y=y, family = poisson(identity),
                                      solver = "eigen", method = 6,
                                      lambda = 1, nonnegative = TRUE,
                                      control.iwls = do.call("control.iwls.gen2", list(maxit = 200)),
                                      control.l0 = do.call("control.l0.gen2", list(maxit = 1000, delta = 1E-10, gamma = 1.5))
)) 
plot_L0glm_benchmark(x = sim$x, y = y,
                     coefs = L0glm_cg1$coefficients / X.n,
                     fitted = L0glm_cg1$fitted.values,
                     a.true = sim$beta_true,
                     main="Ground truth vs L0glm2 estimates")
cor(coefs_true, L0glm_cg1$coefficients / X.n) # 0.7240002 NOT OK
sum(coefs_true>0) # 100
sum(L0glm_cg1$coefficients>0) # 360


#L0glm2 Eigen - CHOL - best results so far
system.time(L0glm_chol0 <- L0glm.fit2(X=Xsp, y=y, family = poisson(identity),
                                    solver = "eigen", method = 8, # or 8?
                                    lambda = 1, nonnegative = TRUE,
                                    control.iwls = do.call("control.iwls.gen2", list(maxit = 1000)),
                                    control.l0 = do.call("control.l0.gen2", list(maxit = 100, delta = 1E-10, gamma = 2)),
                                    control.solver = do.call("control.solver.gen", list(solver = "eigen", method = 8, rowaug = FALSE))))

#19.14s for n=p=10000 (for sparse matrix)
plot_L0glm_benchmark(x = sim$x, y = sim$y,
                     coefs = L0glm_chol0$coefficients / X.n,
                     fitted = L0glm_chol0$fitted.values,
                     a.true = sim$beta_true,
                     main="Ground truth vs L0glm2 estimates")

cor(coefs_true, L0glm_chol0$coefficients) # 0.95 for n=p=10000
sum(coefs_true>0) # 1000
sum(L0glm_chol0$coefficients>0) # 1067

