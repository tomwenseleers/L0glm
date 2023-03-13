# Simulate data
set.seed(1)
n <- 100
p <- 100
s <- 0.1 # sparsity as proportion of p that have nonzero coefficients // p*s == #of peaks
k <- round(p*s) # nr of nonzero covariates
sim <- simulate_spike_train(n=n, p=p, k=k, mean_beta = 1000, sd_logbeta = 1, family="poisson")
# sim <- simulate_spike_train()
coefs_true <- sim$beta_true

X <- sim$X
X.n <- apply(X, 2, norm, type = "2")
X <- sweep(X, 2, X.n, "/")
Xsp <- as(X, "dgCMatrix")
Xsp[Xsp<1E-6] <- 0
y <- sim$y
lambda <- 1
pen_weights <- rep(1, ncol(X))
lambdas <- lambda*pen_weights
lower <- rep(0, ncol(X))
upper <- rep(Inf, ncol(X))
nonnegative <- TRUE
start <- NULL
verbose <- TRUE
weights <- 1/(y+0.1) # observation weights




# L0glm
system.time(L0glm1.fit0 <- L0glm.fit(X=X, y=y, family = poisson(identity),
                                    lambda = 0, nonnegative = TRUE, 
                                    control.iwls = list(maxit = 100, warn = FALSE),
                                    control.l0 = list(maxit = 100, rel.tol = 1E-4, 
                                                      delta = 1E-10, gamma = 2, warn = FALSE) # 1.7s TODO make sure that it works for delta=0 (check nonzero subsetting)
)) #Time = 7.8s # TODO make sure that function also works if you just give subset of 
# system.time(L0glm1 <- L0glm(formula = "y ~ 0 + .", data = data.frame(X, y = y), family = poisson(identity),
#                 lambda = 0, tune.meth = "none", nonnegative = TRUE, normalize = FALSE,
#                 control.iwls = list(maxit = 1),
#                 control.l0 = list(maxit = 10, delta = 1E-10, gamma = 2))) #Time = 7.2s

plot_L0glm_benchmark(x = sim$x, y = y, 
                     coefs = L0glm1.fit0$coefficients / X.n, 
                     fitted = L0glm1.fit0$fitted.values, a.true = sim$beta_true,
                     main="Ground truth vs L0glm estimates")

cor(coefs_true, L0glm1.fit0$coefficients / X.n) # 0.9993351

sum(coefs_true>0) # 10
sum(L0glm1.fit0$coefficients>0) # 21



system.time(L0glm1.fit1 <- L0glm.fit(X=X, y=y, family = poisson(identity),
                            lambda = 1, nonnegative = TRUE, 
                            control.iwls = list(maxit = 1, warn = FALSE, thresh=1E-5),
                            control.l0 = list(maxit = 100, rel.tol = 1E-4, 
                                              delta = 1E-30, gamma = 2, warn = FALSE) # 11.06s
                            )) #Time = 7.8s # TO DO make sure that function also works if you just give subset of 
# system.time(L0glm1 <- L0glm(formula = "y ~ 0 + .", data = data.frame(X, y = y), family = poisson(identity),
#                 lambda = 1, tune.meth = "none", nonnegative = TRUE, normalize = FALSE,
#                 control.iwls = list(maxit = 1),
#                 control.l0 = list(maxit = 10, delta = 1E-10, gamma = 2))) #Time = 7.2s

plot_L0glm_benchmark(x = sim$x, y = y, 
                     coefs = L0glm1.fit1$coefficients / X.n, 
                     fitted = L0glm1.fit1$fitted.values, a.true = sim$beta_true,
                     main="Ground truth vs L0glm estimates")

cor(coefs_true, L0glm1.fit1$coefficients / X.n) # 0.9999708

sum(coefs_true>0) # 10
sum(L0glm1.fit1$coefficients>0) # 10






# test implementation of L0 nonnegative identity link Poisson
beta <- lmridge(X = X, 
                y = y, 
                lambda = rep(0, ncol(X)),
                weights = rep(1, nrow(X)),
                lower = lower,
                upper = upper,
                nonnegative = TRUE,
                solver = "nnls",
                method = NULL,
                start = NULL,
                control.solver = NULL,
                verbose = TRUE)
lambda <- 1
for (iter in 1:500) {
beta[beta<1E-20] <- 0
nonzero <- beta>0
lambdas <- lambda*1/beta^2
fitted <- X[,nonzero,drop=F] %*% beta[nonzero]
beta[nonzero] <- lmridge(X=X[,nonzero,drop=F], 
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
print(sum(beta>0))
}
                         
plot(coefs_true, beta / X.n)                         





#L0glm2 OSQP
system.time(L0glm_osqp0 <- L0glm.fit2(X=Xsp, y=y, family = poisson(identity), 
                 solver = "osqp", lambda = 0, nonnegative = TRUE,
                 control.iwls = do.call("control.iwls.gen2", list(maxit = 100)),
                 control.l0 = do.call("control.l0.gen2", list(maxit = 100, delta = 1E-10, gamma = 2))
            )) # 0.47s
plot_L0glm_benchmark(x = sim$x, y = y, 
                     coefs = L0glm_osqp0$coefficients / X.n, 
                     fitted = L0glm_osqp0$fitted.values, 
                     a.true = sim$beta_true,
                     main="Ground truth vs L0glm2 estimates")

cor(coefs_true, L0glm_osqp0$coefficients / X.n) # 0.9990666

sum(coefs_true>0) # 10
sum(L0glm_osqp0$coefficients>0) # 62


system.time(L0glm_osqp1 <- L0glm.fit2(X=Xsp, y=y, family = poisson(identity), 
                                      solver = "osqp", lambda = 1, nonnegative = TRUE,
                                      control.iwls = do.call("control.iwls.gen2", list(maxit = 100, thresh = 1E-5)),
                                      control.l0 = do.call("control.l0.gen2", 
                                                           list(maxit = 100, delta = 1E-10, gamma = 2))
                                      )) # 7.78s
# system.time(L0glm_osqp1 <- L0glm2(formula = "y ~ 0 + .", data = data.frame(X, y = y), 
#                                     family = poisson(identity), 
#                                     solver = "osqp", lambda = 1, nonnegative = TRUE, 
#                                     lower = lower, upper = upper, 
#                                     control.l0 = list(maxit = 100, delta = 1E-10, gamma = 2),
#                                     control.iwls = list(maxit = 100),
#                                     normalize = FALSE,
#                                     verbose = TRUE))

plot_L0glm_benchmark(x = sim$x, y = y, 
                     coefs = L0glm_osqp1$coefficients / X.n, 
                     fitted = L0glm_osqp1$fitted.values, 
                     a.true = sim$beta_true,
                     main="Ground truth vs L0glm2 estimates")

cor(coefs_true, L0glm_osqp1$coefficients / X.n) # 0.938751
sum(coefs_true>0) # 10
sum(L0glm_osqp1$coefficients>0) # 18

cor(L0glm_osqp0$coefficients / X.n, L0glm_osqp1$coefficients / X.n) 



#L0glm2 GLMNET
system.time(L0glm_glmnet0 <- L0glm.fit2(X=Xsp, y=y, family = poisson(identity), 
                                      solver = "glmnet", lambda = 0, nonnegative = TRUE,
                                      control.iwls = do.call("control.iwls.gen2", 
                                                             list(maxit = 100)),
                                      control.l0 = do.call("control.l0.gen2", 
                                                           list(maxit = 100, delta = 1E-10, gamma = 2))
)) # 0.53s

# system.time(L0glm_glmnet0 <- L0glm2(formula = "y ~ 0 + .", data = data.frame(X, y = y), family = poisson(identity), 
#                                   solver = "glmnet", lambda = 0, nonnegative = nonnegative, lower = lower,
#                                   control.l0 = list(maxit = 100, delta = 1E-10, gamma = 2),
#                                   control.iwls = list(maxit = 100),
#                                   upper = upper, verbose = TRUE))

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
                                        control.iwls = do.call("control.iwls.gen2", list(maxit = 100)),
                                        control.l0 = do.call("control.l0.gen2", list(maxit = 100, delta = 1E-10, gamma = 2))
)) # 3.14s

# system.time(L0glm_glmnet1 <- L0glm2(formula = "y ~ 0 + .", data = data.frame(X, y = y), family = poisson(identity), 
#                                     solver = "glmnet", lambda = 1, nonnegative = TRUE, lower = lower,
#                                     control.l0 = list(maxit = 100, delta = 1E-10, gamma = 2),
#                                     control.iwls = list(maxit = 100),
#                                     upper = upper, verbose = TRUE))

plot_L0glm_benchmark(x = sim$x, y = y, 
                     coefs = L0glm_glmnet1$coefficients / X.n, 
                     fitted = L0glm_glmnet1$fitted.values, 
                     a.true = sim$beta_true,
                     main="Ground truth vs L0glm2 estimates")

cor(coefs_true, L0glm_glmnet1$coefficients / X.n) # 0.9990138
cor(L0glm_glmnet0$coefficients / X.n, L0glm_glmnet1$coefficients / X.n) # 0.9999587
sum(coefs_true>0) # 10
sum(L0glm_glmnet1$coefficients>0) # 19



#L0glm2 Eigen - QR

system.time(L0glm_qr0 <- L0glm.fit2(X=Xsp, y=y, family = poisson(identity), 
                                      solver = "eigen", method = 1, # or 2?
                                      lambda = 1E-10, nonnegative = TRUE,
                                      control.iwls = do.call("control.iwls.gen2", list(maxit = 100)),
                                      control.l0 = do.call("control.l0.gen2", list(maxit = 100, delta = 1E-10, gamma = 2))
)) # 0.50s TODO make sure that all solvers always give a result even for underdetermined cases and lambda=0 -> then set to smallest lambda value possible
plot_L0glm_benchmark(x = sim$x, y = y, 
                     coefs = L0glm_qr0$coefficients / X.n, 
                     fitted = L0glm_qr0$fitted.values, 
                     a.true = sim$beta_true,
                     main="Ground truth vs L0glm2 estimates")
cor(coefs_true, L0glm_qr0$coefficients / X.n) # 0.9253677
sum(coefs_true>0) # 10
sum(L0glm_qr0$coefficients>0) # 15


system.time(L0glm_qr1 <- L0glm.fit2(X=Xsp, y=y, family = poisson(identity), 
                                    solver = "eigen", method = 1,  # or 2?
                                    lambda = 1, nonnegative = TRUE,
                                    control.iwls = do.call("control.iwls.gen2", list(maxit = 100, thresh=1E-5)),
                                    control.l0 = do.call("control.l0.gen2", list(maxit = 100, 
                                                                                 delta = 1E-10, 
                                                                                 gamma = 2))
)) # 48s
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
                                    solver = "eigen", method = 3, # or 4?
                                    lambda = 0, nonnegative = TRUE,
                                    control.iwls = do.call("control.iwls.gen2", list(maxit = 100)),
                                    control.l0 = do.call("control.l0.gen2", list(maxit = 100, delta = 1E-10, gamma = 2))
)) # TODO ERROR Not compatible with requested type: [type=NULL; target=double].
plot_L0glm_benchmark(x = sim$x, y = y, 
                     coefs = L0glm_lscg0$coefficients / X.n, 
                     fitted = L0glm_lscg0$fitted.values, 
                     a.true = sim$beta_true,
                     main="Ground truth vs L0glm2 estimates")
cor(coefs_true, L0glm_lscg0$coefficients / X.n) # 
sum(coefs_true>0) # 10
sum(L0glm_lscg0$coefficients>0) # 


system.time(L0glm_lscg1 <- L0glm.fit2(X=Xsp, y=y, family = poisson(identity), 
                                    solver = "eigen", method = 3, 
                                    lambda = 1, nonnegative = TRUE,
                                    control.iwls = do.call("control.iwls.gen2", list(maxit = 100)),
                                    control.l0 = do.call("control.l0.gen2", list(maxit = 100, delta = 1E-10, gamma = 2)
                                    # control.solver = list(solver="eigen", 
                                    #                       method = "LSCG", tol = 0, maxit = 0L, precond = 0L)
                                    )
)) # ERROR Not compatible with requested type: [type=NULL; target=double].
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
                                      lambda = 0, nonnegative = TRUE,
                                      control.iwls = do.call("control.iwls.gen2", list(maxit = 100)),
                                      control.l0 = do.call("control.l0.gen2", list(maxit = 100, delta = 1E-10, gamma = 2))
)) # TODO ERROR Not compatible with requested type: [type=NULL; target=double].
plot_L0glm_benchmark(x = sim$x, y = y, 
                     coefs = L0glm_cg0$coefficients / X.n, 
                     fitted = L0glm_cg0$fitted.values, 
                     a.true = sim$beta_true,
                     main="Ground truth vs L0glm2 estimates")
cor(coefs_true, L0glm_cg0$coefficients / X.n) # 
sum(coefs_true>0) # 100
sum(L0glm_cg0$coefficients>0) # 


system.time(L0glm_cg1 <- L0glm.fit2(X=Xsp, y=y, family = poisson(identity), 
                                      solver = "eigen", method = 5, 
                                      lambda = 1, nonnegative = TRUE,
                                      control.iwls = do.call("control.iwls.gen2", list(maxit = 100)),
                                      control.l0 = do.call("control.l0.gen2", list(maxit = 100, delta = 1E-10, gamma = 2))
)) # ERROR Not compatible with requested type: [type=NULL; target=double].
plot_L0glm_benchmark(x = sim$x, y = y, 
                     coefs = L0glm_cg1$coefficients / X.n, 
                     fitted = L0glm_cg1$fitted.values, 
                     a.true = sim$beta_true,
                     main="Ground truth vs L0glm2 estimates")
cor(coefs_true, L0glm_cg1$coefficients / X.n) # 0.7240002 NOT OK
sum(coefs_true>0) # 100
sum(L0glm_cg1$coefficients>0) # 360


