devtools::load_all()

# Simulate data
set.seed(1)
n <- 1000
p <- 1000
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

block.size <- 0.1 * p
block.ids <- ceiling((1:p)/block.size)

#OSQP
system.time(L0glm_osqp <- L0glm.fit2(X=Xsp, y=y, family = poisson(identity),
                                      solver = "osqp", lambda = 1, nonnegative = TRUE, block = FALSE, 
                                      control.iwls = do.call("control.iwls.gen2", list(maxit = 100)),
                                      control.l0 = do.call("control.l0.gen2", list(maxit = 100, delta = 1E-10, gamma = 2)),
                                     control.fit = do.call("control.fit.gen2", list(maxit = 100, block.size = block.size)),
                                      control.solver = do.call("control.solver.gen", list(solver = "osqp", rowaug = TRUE)))) #1.09s for n=p=1000

system.time(L0glm_block_osqp <- L0glm.fit2(X=Xsp, y=y, family = poisson(identity),
                                      solver = "osqp", lambda = 1, nonnegative = TRUE, block = TRUE,
                                      control.iwls = do.call("control.iwls.gen2", list(maxit = 100)),
                                      control.l0 = do.call("control.l0.gen2", list(maxit = 100, delta = 1E-10, gamma = 2)),
                                      control.fit = do.call("control.fit.gen2", list(maxit = 100, block.size = block.size)),
                                      control.solver = do.call("control.solver.gen", list(solver = "osqp", rowaug = TRUE)))) #3.93s for n=p=1000 and block.size = 300

plot_L0glm_benchmark(x = sim$x, y = y,
                     coefs = L0glm_block_osqp$coefficients / X.n,
                     fitted = L0glm_block_osqp$fitted.values,
                     a.true = sim$beta_true,
                     main="Ground truth vs L0glm2 estimates")

cor(coefs_true, L0glm_block_osqp$coefficients) # 0.93
cor(coefs_true, L0glm_osqp$coefficients)
cor(L0glm_osqp$coefficients, L0glm_block_osqp$coefficients)




#GLMNET
system.time(L0glm_glmnet <- L0glm.fit2(X=Xsp, y=y, family = poisson(identity),
                                     solver = "glmnet", lambda = 1, nonnegative = TRUE, block = FALSE, 
                                     control.iwls = do.call("control.iwls.gen2", list(maxit = 100)),
                                     control.l0 = do.call("control.l0.gen2", list(maxit = 100, delta = 1E-10, gamma = 2)),
                                     control.fit = do.call("control.fit.gen2", list(maxit = 100, block.size = block.size)),
                                     control.solver = do.call("control.solver.gen", list(solver = "glmnet")))) #34.81 for n=p=3000

system.time(L0glm_block_glmnet <- L0glm.fit2(X=Xsp, y=y, family = poisson(identity),
                                           solver = "glmnet", lambda = 1, nonnegative = TRUE, block = TRUE,
                                           control.iwls = do.call("control.iwls.gen2", list(maxit = 100)),
                                           control.l0 = do.call("control.l0.gen2", list(maxit = 100, delta = 1E-10, gamma = 2)),
                                           control.fit = do.call("control.fit.gen2", list(maxit = 100, block.size = block.size)),
                                           control.solver = do.call("control.solver.gen", list(solver = "glmnet"))))

plot_L0glm_benchmark(x = sim$x, y = y,
                     coefs = L0glm_block_glmnet$coefficients / X.n,
                     fitted = L0glm_block_glmnet$fitted.values,
                     a.true = sim$beta_true,
                     main="Ground truth vs L0glm2 estimates")

cor(coefs_true, L0glm_block_glmnet$coefficients) # 0.45
cor(coefs_true, L0glm_glmnet$coefficients) #0.997
cor(L0glm_glmnet$coefficients, L0glm_block_glmnet$coefficients)
