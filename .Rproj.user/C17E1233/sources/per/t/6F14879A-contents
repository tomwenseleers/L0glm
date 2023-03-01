####---- SCRIPT DESCRIPTION ----####


# This script shows how to use the package 'L0glm' and performs some test and
# benchmarks as comparison to other popular softwares and algorithms

####---- SETUP ENVIRONMENT ----####


setwd("D:/Documents/GitHub/L0glm/Paper")
# devtools::install_github("tomwenseleers/L0glm")
library(L0glm)
library(microbenchmark)
library(export)
library(ggplot2)
graph2ppt(file = "Github/graphs") # Initialize the ppt file


####---- SHOWCASE ON SIMULATED DATA ----####


# Simulate some data
sim <- simulate_spike_train(Plot = TRUE)
x <- sim$X
y <- sim$y

# Set up the parameters for controlling the algorithm
ctrl.fit <- control.fit.gen() # default
ctrl.iwls <- control.iwls.gen() # default
ctrl.l0 <- control.l0.gen() # default

# Fit a GLM with Poisson error structure and identity link, with nonnegativity
# constraints on the coefficients, and L0 penalty
L0glm.out <- L0glm(formula = y ~ 0 + ., data = data.frame(y = y, x),
                   family = poisson(identity), normalize = FALSE,
                   lambda = 1, tune.meth = "none", nonnegative = TRUE,
                   control.iwls = ctrl.iwls, control.l0 = ctrl.l0,
                   control.fit = ctrl.fit)

# Perform inference on the coefficients. The function will automatically choose
# the correct inference procedure (non parametric bootstrapping in this case):
L0glm.infer.out <- L0glm.inference(L0glm.out, level = 0.95, boot.repl = 1000,
                                   control.l0 = ctrl.l0, control.iwls = ctrl.iwls,
                                   control.fit = ctrl.fit)

# Plot the results
plot_L0glm_benchmark(x = sim$x, y = y, fit = L0glm.out, a.true = sim$a,
                     # inference = L0glm.infer.out,
                     main = "Ground truth vs L0 penalized L0glm estimates")


####---- BENCHMARK WITH GLM ----####


# From ?glm examples
# Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)

# Becnhmark
microbenchmark(
  # Glm fitting
  "glm" = {
    glm_fit <- glm(counts ~ outcome + treatment, family = poisson())
  },
  # L0glm fitting (using glm settings)
  "L0glm (glm settings)" = {
    L0glm_fit <- L0glm(counts ~ outcome + treatment,
                       family = poisson(),
                       lambda = 0, tune.meth = "none", nonnegative = FALSE,
                       control.iwls = list(maxit = 25, thresh = .Machine$double.eps),
                       control.l0 = list(maxit = 1),
                       control.fit = list(maxit = 1), verbose = FALSE)
  },
  times = 25
)

# Compare coefficients
df <- data.frame(coef.glm = coef(glm_fit),
                 coef.L0glm = coef(L0glm_fit))
abs(df$coef.glm - df$coef.L0glm)

# Plot coefficients
data <- data.frame(y = unlist(df),
                   x = rep(1:nrow(df), ncol(df)),
                   type = rep(c("glm", "L0glm"), each = nrow(df)))
pl <- ggplot(data = data, aes(x = x, y = y, color = type)) +
  geom_point() + geom_line() +
  ggtitle("Compare coefficients estimated using glm or L0glm") +
  ylab("Estimate") + xlab("Index") +
  scale_colour_manual(name = "Algorithm",
                      values = c(glm = "red3", L0glm = "green4"),
                      labels = c(glm = "glm", L0glm = "L0glm"))
graph2ppt(pl, file = "Github/graphs", scaling = 50, append = TRUE)

# Conclusion
# Both algorithms give almost exactly the same solution (up to 2E-15). The
# higher functionnality of the L0glm framework (possibility of nonnegative
# constraints and regularization) is a the cost of timing performance.


####---- BENCHMARK WITH GLMNET, RIDGE AND PENALIZED ----####


library(glmnet)
library(ridge)
library(penalized)

# Simulate data with Gaussian noise
set.seed(123)
n <- 100
p <- 20
x <- matrix(rnorm(n*p), nrow = n, ncol = p)
x <- scale(x)
beta <- runif(p, min = -1)
y0 <- x %*% beta
y <- y0 + rnorm(n, mean = 0, sd = 2.5)
y <- scale(y)
lam <- 1

microbenchmark(
  # Ridge regression using glmnet
  "glmnet" = {
    glmnet_fit <- glmnet(x = x, y = y, family = "gaussian", alpha = 0,
                         standardize = FALSE, thresh = .Machine$double.eps,
                         lambda = 10^seq(10,log10(lam)), intercept = TRUE)
  },
  # Ridge regression using ridge
  "ridge" = {
    ridge_fit <- linearRidge(y ~ 1 + ., data = data.frame(y = y, x),
                             lambda = lam*n, scaling = "none")
  },
  "penalized" = {
    penal_fit <- penalized(response = y, penalized = x,
                           lambda1 = 0, lambda2 = lam*n, positive = FALSE, model= "linear",
                           epsilon = .Machine$double.eps, maxiter = 25, trace = F)
  },
  # L0glm fitting (using glm settings)
  "L0glm (ridge settings)" = {
    L0glm_fit <- L0glm(y ~ 1 + ., data = data.frame(y = y, x),
                       family = gaussian(),
                       lambda = lam*n, tune.meth = "none", nonnegative = FALSE,
                       control.iwls = list(maxit = 25, thresh = .Machine$double.eps),
                       control.l0 = list(maxit = 1),
                       control.fit = list(maxit = 1, tol = .Machine$double.eps),
                       verbose = FALSE)
  },
  times = 25
)
# Note: all lambdas except for glmnet are scaled with 'n'. This is to allow
# comparing the same objective function as illustrated in https://stackoverflow.com/questions/39863367/ridge-regression-with-glmnet-gives-different-coefficients-than-what-i-compute

# Check results
df <- data.frame(coef.glmnet = as.vector(coef(glmnet_fit, s = 1)),
                 coef.ridge = coef(ridge_fit),
                 coef.penalized = coef(penal_fit),
                 coef.L0glm = coef(L0glm_fit),
                 coef.true = c(0, beta)) # no intercept is present
abs(df$coef.glmnet - df$coef.L0glm)
abs(df$coef.ridge - df$coef.L0glm)
abs(df$coef.penalized - df$coef.L0glm)

# Plot coefficients
data <- data.frame(y = unlist(df),
                   x = rep(1:nrow(df), ncol(df)),
                   type = rep(c("glmnet", "ridge", "penalized", "L0glm", "true"), each = nrow(df)))
pl <- ggplot(data = data, aes(x = x, y = y, color = type)) +
  geom_point() + geom_line() +
  ggtitle("Compare true coefficients with coefficients estimated \nusing glmnet, ridge, or L0glm") +
  ylab("Estimate") + xlab("Index") +
  scale_colour_manual(name = "Algorithm",
                      values = c(glmnet = "red3", ridge = "orange2", penalized = "purple", L0glm = "green4", true = "grey40"),
                      labels = c(glmnet = "glmnet", ridge = "ridge", penalized = "penalized", L0glm = "L0glm", true = "True"))
graph2ppt(pl, file = "Github/graphs", scaling = 50, append = TRUE)

# Conclusion
# Both algorithms exhibit coefficients following the same trend but the absolute
# coefficients values are quite different. However, L0glm seems to better
# estimate the true value of the coefficients

# TODO delete, optimizing lambda
test1 <- cv.glmnet(x = x, y = y, family = "gaussian", alpha = 0,
                   standardize = FALSE, thresh = .Machine$double.eps,
                   nfolds = 3,
                   lambda = 10^seq(10,-10), intercept = FALSE)
plot(test1$lambda, test1$cvm,log= "xy", type = 'l')
test1$lambda.min
test2 <- L0glm(y ~ 0 + ., data = data.frame(y = y, x), family = gaussian(),
               lambda = 10^seq(10,-10), tune.meth = "3-fold", tune.crit = "rss",
               nonnegative = FALSE, control.l0 = list(maxit = 1),
               control.iwls = list(maxit = 25, thresh = .Machine$double.eps),
               control.fit = list(maxit = 1), verbose = TRUE)
plot(test2$lambda.tune$lambdas, test2$lambda.tune$IC[, "rss"],log= "xy", type = 'l')
test2$lambda.tune$best.lam


####---- BENCHMARK WITH L0Learn, L0ara, bestsubset ----####


library(L0Learn)
library(l0ara)
library(bestsubset)

# From ?L0Learn.fit examples
# Generate synthetic data for this example
n <- 200
p <- 500
k <- 10
data <- GenSynthetic(n = n, p = p, k = k, seed = 123)
x <- data$X
# x <- scale(x)
# beta <-  c(rep(1, k), rep(0, p - k))
# y <- x %*% beta + rnorm(n)
# y <- scale(y)
y <- data$y

microbenchmark(
  # L0 penalized regression using L0Learn
  "L0Learn" = {
    L0Learn_fit <- L0Learn.fit(x = x, y = y, penalty="L0", maxSuppSize = ncol(X),
                               nGamma = 0, autoLambda = FALSE, lambdaGrid = list(1.56E-2),
                               tol = 1E-7)
  },
  # L0 penalized regression using L0ara
  "L0ara" = {
    L0ara_fit <- l0ara(x = x, y = y, family = "gaussian", lam = 2.5,
                       standardize = F, eps = 1E-7)
  },
  # Best subset regression using bestsubset
  "bestsubset" = {
    bs_fit <- bs(x = x, y = y, k = k, intercept = TRUE,
                 form = ifelse(nrow(x) < ncol(x), 2, 1), time.limit = 5, nruns = 50,
                 maxiter = 1000, tol = 1e-7, polish = TRUE, verbose = FALSE)
  },
  # L0 penalized regression using L0glm
  "L0glm" = {
    L0glm_fit <- L0glm(y ~ 1 + ., data = data.frame(y = y, x),
                       family = gaussian(),
                       lambda = 2.5, tune.meth = "none", nonnegative = FALSE,
                       control.iwls = list(maxit = 100, thresh = 1E-7),
                       control.l0 = list(maxit = 100, rel.tol = 1E-7),
                       control.fit = list(maxit = 1), verbose = FALSE)
  },
  times = 5
)
# Note that bestsubset is optimized using the true number of nonzero coefficient
# because tuning it was much to slow. The algorithm check solution using
# Gurobi's mixed integer program solver which is very slow for k = 10, so time
# limit was set to 5 s which dramatically overestimates to true time
# performance of bestsubset

# Check results
df <- data.frame(coef.L0Learn = as.numeric(L0Learn_fit$beta[[1]]),
                 coef.L0ara = L0ara_fit$beta,
                 coef.bestsubset = as.vector(bs_fit$beta),
                 coef.L0glm = coef(L0glm_fit)[-1],
                 coef.true = beta)
all(df[(k+1):p,] == 0)
# No false positives !
abs(df$coef.L0Learn - df$coef.bestsubset)[1:k]
abs(df$coef.L0glm - df$coef.L0ara)[1:k]
abs(df$coef.L0glm - df$coef.L0Learn)[1:k]

data <- data.frame(y = unlist(df[1:k,]),
                   x = rep(1:k, ncol(df[1:k,])),
                   type = rep(c("L0Learn", "L0ara", "bestsubset", "L0glm", "true"), each = k))
pl <- ggplot(data = data, aes(x = x, y = y, color = type)) +
  geom_point() + geom_line() +
  ggtitle("Compare true coefficients with coefficient estimated \nusing bestsubset, L0ara, L0glm, L0Learn") +
  ylab("Estimate") + xlab("Index") +
  scale_colour_manual(name = "Algorithm",
                      values = c(L0Learn = "red3", L0ara = "orange2",
                                 bestsubset = "purple", L0glm = "green4", true = "grey40"),
                      labels = c(L0Learn = "L0Learn", L0ara = "L0ara",
                                 bestsubset = "bestsubset", L0glm = "L0glm", true = "True"))
graph2ppt(pl, file = "Github/graphs", scaling = 50, append = TRUE)

# Conclusion
# All algorithms find the correct set of coefficients
# The solution with L0Learn is very similar to the solution found with bestsubset
# (up to 5E-4). The L0glm and L0ara find almost the same solution (up to 1E-8).
# However, there is a noticeable difference between L0glm (and hence L0ara) and
# bestsubset (and hence L0Learn) and this difference seems to be systematic across
# coefficients.



# TODO delete
# Select lambda
# L0Learn
test1 <- L0Learn.cvfit(x = x, y = y, loss = "SquaredError", penalty = "L0",
                       algorithm = "CD", maxSuppSize = 100, nLambda = 100, nGamma = 0,
                       maxIters = 200, tol = 1e-7, autoLambda = TRUE, nFolds = 3, seed = 123)
plot(test1$fit$lambda[[1]], test1$cvMeans[[1]], log = "xy", type = "l")
test1$fit$lambda[[1]][which.min(test1$cvMeans[[1]])]
# L0ara
test2 <- cv.l0ara(x = x, y = y, family = "gaussian", lam = 10^seq(-2,2, length.out = 21),
                  measure = "mse", nfolds = 3, eps = 1E-7)
plot(test2$lambda, test2$cv.error, type = "l", log = "xy")
test2$lam.min
# bestsubset
test3 <- bs(x = x, y = y, k = 0:15, intercept = TRUE,
            form = ifelse(nrow(x) < ncol(x), 2, 1), time.limit = 100, nruns = 50,
            maxiter = 1000, tol = 1e-07, polish = TRUE, verbose = T)
# This is super slow to compute ... Best k will be set to true k = 10
# L0glm
test4 <- L0glm(y ~ 1 + ., data = data.frame(y = y, x),
               family = gaussian(),
               lambda = 10^seq(-2,2, length.out = 21),
               tune.meth = "3-fold", tune.crit = "rss", nonnegative = FALSE,
               control.iwls = list(maxit = 100, thresh = 1E-7),
               control.l0 = list(maxit = 100, rel.tol = 1E-7),
               control.fit = list(maxit = 1), verbose = TRUE)
plot(test4$lambda.tune$lambdas, test4$lambda.tune$IC[,"rss"], type = "l", log = "xy")
test4$lambda.tune$best.lam


####---- COMPARE L0 PENALTY WITH LASSO, MCP, OR SCAD PENALTY ----####


library(ncvreg)

# From ?L0Learn.fit examples
# Generate synthetic data for this example
n <- 200
p <- 500
k <- 10
data <- GenSynthetic(n = n, p = p, k = k, seed = 123)
beta <-  c(rep(1, k), rep(0, p - k))
x <- data$X
y <- data$y

microbenchmark(
  "lasso" = {
    lasso_fit <- ncvreg(X = x, y = y, family = "gaussian", penalty = "lasso",
                        alpha = 1, lambda = 0.1278381, eps = 1e-7)
  },
  "MCP" = {
    mcp_fit <- ncvreg(X = x, y = y, family = "gaussian", penalty = "MCP",
                        alpha = 1, lambda = 0.1678555, eps = 1e-7)
  },
  "SCAD" = {
    scad_fit <- ncvreg(X = x, y = y, family = "gaussian", penalty = "SCAD",
                        alpha = 1, lambda = 0.1442869, eps = 1e-7)
  },
  "L0" = {
    L0glm_fit <- L0glm(y ~ 1 + ., data = data.frame(y = y, x),
                       family = gaussian(),
                       lambda = 2.5, tune.meth = "none", nonnegative = FALSE,
                       control.iwls = list(maxit = 100, thresh = 1E-7),
                       control.l0 = list(maxit = 100, rel.tol = 1E-7),
                       control.fit = list(maxit = 1), verbose = FALSE)
  },
  times = 5
)

# Check results
df <- data.frame(coef.lasso = coef(lasso_fit)[-1], # remove intercept
                 coef.mcp = coef(mcp_fit)[-1],
                 coef.scad = coef(scad_fit)[-1],
                 coef.L0 = coef(L0glm_fit)[-1],
                 coef.true = beta)
FP <- colSums(df[(k+1):p,] != 0)
TP <- colSums(df[1:k,] != 0)
FN <- colSums(df[1:k,] == 0)
TN <- colSums(df[(k+1):p,] == 0)
TP/(TP + FN) # sensitivity
TN/(TN + FP) # specificity

# Plot results
data <- data.frame(y = unlist(df),
                   x = rep(1:(p+1), ncol(df)),
                   type = rep(c("lasso", "mcp", "scad", "L0", "true"), each = p+1))
nz <- (1:(p+1)) %in% 2:(k+1)
pl <- ggplot(data = data[nz,], aes(x = x, y = y, color = type)) +
  geom_point() + geom_line() +
  ggtitle("Compare true nonzero coefficients with coefficient estimated \nusing lasso, MCP, SCAD, or L0 penalties") +
  ylab("Estimate") + xlab("Index") +
  scale_colour_manual(name = "Algorithm",
                      values = c(lasso = "red3", mcp = "orange2", scad = "purple", L0 = "green4", true = "grey40"),
                      labels = c(lasso = "lasso", mcp = "MCP", scad = "SCAD",  L0 = "L0glm", true = "True"))
graph2ppt(pl, file = "Github/graphs", scaling = 50, append = TRUE)


# TODO delete
# Opimize lambda
lasso_cv <- cv.ncvreg(X = x, y = y, family = "gaussian", penalty = "lasso",
                      alpha = 1, eps = 1e-9, nfolds = 3, seed = 123)
plot(x = lasso_cv$lambda, y = lasso_cv$cve, log = "xy", type = "l")
lasso_cv$lambda.min # 0.1278381
mcp_cv <- cv.ncvreg(X = x, y = y, family = "gaussian", penalty = "MCP",
                    alpha = 1, eps = 1e-9, nfolds = 3, seed = 123)
plot(x = mcp_cv$lambda, y = mcp_cv$cve, log = "xy", type = "l")
mcp_cv$lambda.min # 0.1678555
scad_cv <- cv.ncvreg(X = x, y = y, family = "gaussian", penalty = "SCAD",
                     alpha = 1, eps = 1e-9, nfolds = 3, seed = 123)
plot(x = scad_cv$lambda, y = scad_cv$cve, log = "xy", type = "l")
scad_cv$lambda.min # 0.1442869
L0glm_cv <- L0glm(y ~ 1 + ., data = data.frame(y = y, x),
                  family = gaussian(),
                  lambda = 10^seq(-1, 1, length.out = 100), tune.meth = "3-fold",
                  tune.crit = "rss", nonnegative = FALSE, seed = 123,
                  control.iwls = list(maxit = 100, thresh = 1E-4),
                  control.l0 = list(maxit = 100, rel.tol = 1E-7, warn = TRUE),
                  control.fit = list(maxit = 1), verbose = TRUE)
plot(x = L0glm_cv$lambda.tune$lambdas, y = L0glm_cv$lambda.tune$IC[,"rss"], log = "xy", type = "l")
L0glm_cv$lambda.tune$best.lam # 2.25702







#### EXTRA TESTS ####



#### TEST NEED FOR STANDARDIZATION ####


# Simulate some data
set.seed(1234)
n = 200
npeaks = 20
sim <- simulate_spike_train(Plot = TRUE, n = n)
x <- sim$X
scales <- runif(n, min = -1E2, max = 1E2)
x <- sweep(x, 2, scales, "*")
beta <- runif(n, min = -1, max = 1)
beta[sample(1:length(beta), size = length(beta)-25)] <- 0
y0 <- x %*% beta
y <- y0  + rnorm(n, sd = 2.5)
plot(y)
lines(y0)
abline(h = 0)

par(mfrow = c(3,1))
plot(y, main = "True model")
matlines(sweep(x, 2, beta, "*"), lty = 1, type = "l", col = "orange2")
# No normalization
L0glm.out <- L0glm(formula = y ~ 0 + ., data = data.frame(y = y, x),
                   family = gaussian(identity),
                   lambda = 1000, tune.meth = "none", nonnegative = FALSE,
                   control.iwls = control.iwls.gen(maxit = 1),
                   control.l0 = control.l0.gen(),
                   control.fit = control.fit.gen())
plot(y, main = "No normalization")
matlines(sweep(x, 2, L0glm.out$coefficients, "*"), lty = 1, type = "l", col = "orange2")
# With normalization
x <- sweep(x, 2, apply(x, 2, norm, type = "2"), "/")
L0glm.out <- L0glm(formula = y ~ 0 + ., data = data.frame(y = y, x),
                   family = gaussian(identity),
                   lambda = 1, tune.meth = "none", nonnegative = FALSE,
                   control.iwls = control.iwls.gen(maxit = 1),
                   control.l0 = control.l0.gen(),
                   control.fit = control.fit.gen())
plot(y, main = "With normalization")
matlines(sweep(x, 2, L0glm.out$coefficients, "*"), lty = 1, type = "l", col = "orange2")



#### TEST L0GLM AGAINST GLM ####

# Compare L0glm against stats::glm
X.sub <- X[, seq(1,ncol(X), by = 5)] # Subsampling needed otherwise
# gaussian: identity, log and inverse
# binomial: logit, probit, cauchit, log, and cloglog
# Gamma: inverse, identity and log
# poisson: log, identity, and sqrt
# inverse.gaussian: 1/mu^2, inverse, identity and log.
# TODO test the qusi families, quasi: logit, probit, cloglog, identity, inverse, log, 1/mu^2 and sqrt
fam <- poisson()
microbenchmark(cpglm <- pen.nnglm(X = X.sub, y = y, family = fam, lambda = 0, maxit.l0 = 1, maxit.iwnnls = 25, cv = "none", constr = "none"),
               glm <- glm.fit(x = X.sub, y = y, family = fam, control = list(maxit = 25, epsilon = 1E-8), intercept = FALSE),
               times = 50)
cbind(cpglm = cpglm$coefficients, glm = glm$coefficients, diff.abs = abs(cpglm$coefficients - glm$coefficients))

par(mfrow = c(2,1))
plot(cpglm$coefficients, type = "l", col = 4)
lines(glm$coefficients, col = 2)
plot(y, pch = 16, cex = 0.5, col = "grey")
lines(cpglm$fitted.values, col = 4)
lines(glm$fitted.values, col = 2)

par(mfrow = c(1,1))
plot(y, pch = 16, col = "grey")
lines(y_nonoise, col = 2)
lines(cpglm$fitted.values, col = "orange2")






# Simulate data
get.data.poisson.nn <- function(n = 200, npeaks = 20,seed = 123,
                                peakhrange = c(10,1E3), Plot = F,
                                envir = .GlobalEnv){
  set.seed(seed)
  x = 1:n
  # unkown peak locations
  u = sample(x, npeaks, replace=FALSE)
  # unknown peak heights
  h = 10^runif(npeaks, min=log10(min(peakhrange)), max=log10(max(peakhrange)))
  # locations of spikes of simulated spike train, which are assumed to be unknown
  # here, and which needs to be estimated from the measured total signal
  a = rep(0, n)
  a[u] = h
  # peak shape function
  gauspeak = function(x, u, w, h=1) h*exp(((x-u)^2)/(-2*(w^2)))
  # banded matrix with peak shape measured beforehand
  X = do.call(cbind, lapply(1:n, function (u) gauspeak(x, u=u, w=5, h=1) ))
  # X = Matrix(X, sparse = T)
  # colnames(X) <- paste("Var", 1:ncol(X))
  # noiseless simulated signal = linear convolution of spike train with peak shape function
  y_nonoise = as.vector(X %*% a)
  # simulated signal with random poisson noise
  y <- rpois(n, y_nonoise)

  # Plot the data
  if(Plot){
    par(mfrow=c(1,1))
    plot(y, type = "l", ylab = "Signal", xlab = "x",
         main = "Simulated spike train (red) to be estimated given known blur kernel & with Poisson noise")
    lines(a, type = "h", col = "red")
  }
  l <- mget(ls(environment()))
  invisible(list2env(x = l[!names(l) %in% c("envir", "Plot", "gauspeak", "seed")], envir = envir))
}

################################################################
#### COMPARE OUR PACKAGE VS STATS GLM FOR NON PENALIZED GLM ####
################################################################


sim <- simulate_spike_train(Plot = TRUE, n = n)
X <- sim$X
y <- sim$y
# Test our package vs glm
X.sub <- X[,seq(1,ncol(X), by = 5)]
# gaussian: identity, log and inverse
# binomial: logit, probit, cauchit, log, and cloglog
# Gamma: inverse, identity and log
# poisson: log, identity, and sqrt
# inverse.gaussian: 1/mu^2, inverse, identity and log.
# TODO test the qusi families, quasi: logit, probit, cloglog, identity, inverse, log, 1/mu^2 and sqrt
fam <- poisson()
microbenchmark(cpglm <- L0glm(y ~ 0 + ., data = data.frame(y = y, X.sub), family = fam, lambda = 0,
                              control.l0 = control.l0.gen(maxit = 1),
                              control.iwls = control.iwls.gen(maxit = 1),
                              nonnegative = FALSE, tune.meth = "none", verbose = FALSE),
               glm <- glm.fit(x = X.sub, y = y, family = fam, control = list(maxit = 25, epsilon = 1E-8), intercept = FALSE),
               times = 50)
cbind(cpglm = cpglm$coefficients, glm = glm$coefficients, diff.abs = abs(cpglm$coefficients - glm$coefficients))

par(mfrow = c(2,1))
plot(cpglm$coefficients, type = "l", col = 4)
lines(glm$coefficients, col = 2)
plot(y, pch = 16, cex = 0.5, col = "grey")
lines(cpglm$fitted.values, col = 4)
lines(glm$fitted.values, col = 2)

par(mfrow = c(1,1))
plot(y, pch = 16, col = "grey")
lines(y_nonoise, col = 2)
lines(cpglm$fitted.values, col = "orange2")


# Binomial case

set.seed(123)
n <- 200
ncov <- 25
beta.t <- rep(0, n)
nz <- sample(1:n, ncov)
beta.t[nz] <- runif(ncov, min = -10, max = 10)
eta.true <- X %*% beta.t
y.true <- binomial()$linkinv(eta.true)
y <- as.numeric(y.true >= 0.5)
plot(y, pch = 15)
lines(y.true, col = 2)

# Test our package vs glm
X.sub <- X[,seq(1,ncol(X), by = 5)]
# gaussian: identity, log and inverse
# binomial: logit, probit, cauchit, log, and cloglog
# Gamma: inverse, identity and log
# poisson: log, identity, and sqrt
# inverse.gaussian: 1/mu^2, inverse, identity and log.
# TODO test the qusi families, quasi: logit, probit, cloglog, identity, inverse, log, 1/mu^2 and sqrt
fam <- binomial()
microbenchmark(cpglm <- pen.nnglm(X = X.sub, y = y, family = fam, lambda = 0, maxit.l0 = 1, maxit.iwnnls = 25, cv = "none", constr = "none"),
               glm <- glm.fit(x = X.sub, y = y, family = fam, control = list(maxit = 25, epsilon = 1E-8), intercept = FALSE),
               times = 2)
cbind(cpglm = cpglm$coefficients, glm = glm$coefficients, diff.abs = abs(cpglm$coefficients - glm$coefficients))
# Same results
par(mfrow = c(2,1))
plot(cpglm$coefficients, type = "l", col = 4)
lines(glm$coefficients, col = 2)
plot(y, pch = 16, cex = 0.5, col = "grey")
lines(cpglm$fitted.values, col = 4)
lines(glm$fitted.values, col = 2)

# Test to extract only coefficients that are positively correlated with
# increased risk = nonegativity constrained
cpglm <- pen.nnglm(X = X, y = y, family = fam, lambda = 10^seq(-5,2), maxit.l0 = 25, maxit.iwnnls = 25, cv = "loocv", constr = "nonneg")
plot(eta.true, type = "l")
lines(beta.t, col = 2, type = "h")
lines(X %*% cpglm$coefficients, col = "green4")
lines(cpglm$coefficients, col = "green4", type = "h")


############################################################
#### COMPARE OUR PACKAGE VS OTHER PACKAGE FOR RIDGE GLM ####
############################################################


get.data.poisson.nn()
X.sub <- X[,seq(1,ncol(X), by = 5)]
cpglm <- L0glm(X = X.sub, y = y, family = gaussian(identity),
               control.l0 = control.l0(maxit = 1),
               control.fit = control.fit(constr = "none"),
               lambda = 1E-3, cv = "none")
library(ridge)
ridge.test <- linearRidge(y ~ 0 + ., data = data.frame(y = y, as.matrix(X.sub)),
                          lambda = 1E-3, nPCs = NULL, scaling = "none")
library(glmnet)
glmnet.fit <- glmnet(x = X.sub, y = y, family = "gaussian", lambda = 1E-3,
                     alpha = 0, standardize = F, intercept = F)
library(penalized)
pen.fit <- penalized(response = y, penalized = as.matrix(X.sub),
                     unpenalized = ~ 0, lambda1 = 0, lambda2 = 1E-3,
                     model = "linear", standardize = F)

plot(cpglm$coefficients, type = "l")
lines(ridge.test$coef, col = 2)
lines(glmnet.fit$beta, col = 4)
lines(pen.fit@penalized, col = "orange2")


cpglm <- pen.nnglm(X = X, y = y, family = gaussian(identity), constr = "none",
                   lambda = 200, maxit.l0 = 25, cv = "none")
glmnet.fit <- glmnet(x = X, y = y, family = "gaussian", lambda = 10, thresh = 1E-16,
                     alpha = 1, standardize = F, intercept = F)
par(mfrow = c(2,1))
plot.nnglm.benchmark(x = x, y = y, fit = nnglm_pois_L0, a.true = a,
                     main="Ground truth vs L0 penalized NNGLM estimates")
glmnet.fit$coefficients = as.vector(glmnet.fit$beta)
glmnet.fit$fitted.values = X %*% glmnet.fit$coefficients
plot_L0glm_benchmark(x = x, y = y, fit = glmnet.fit, a.true = a,
                     main="Ground truth vs L0 penalized NNGLM estimates")


################################################################
#### COMPARE OUR PACKAGE VS OTHER PACKAGE FOR L0 REGRESSION ####
################################################################


library(L0Learn)
# Simulate some data
get.data.poisson.nn()
# L0Learn
system.time(L0fit <- L0Learn.fit(x=X, y=y, penalty="L0", algorithm="CD", nLambda=50, maxSuppSize = 200))
# Note algortihm == "CDPSI"should give better results but this is not the case here
# L0glm
system.time(L0glmfit <- L0glm(X = X, y = y, family = gaussian(identity), intercept = F,
                              lambda = 10^seq(-1, 4, length.out = 51), cv = "none",
                              control.fit = control.fit(constr = "none"),
                              control.iwls = control.iwls(rel.tol = 1E-4, maxit = 100, thresh = 1E-1),
                              control.l0 = control.l0(rel.tol = 1E-4, maxit = 100, delta = 1E-3), verbose = F))
# check coefficient trace
matplot(x = unlist(L0fit$lambda), y = t(abs(as.matrix(L0fit$beta[[1]]))), type = "l", log = "xy")
matplot(x = colnames(L0glmfit$coefficients), y = t(abs(L0glmfit$coefficients)), type = "l", log = "xy")

# Check best solution
par(mfrow = c(3,1))
# True
plot(y, main = "True")
lines(X %*% a, col = 2)
matlines(sweep(X, 2, a, "*"), lty = 1, col = "green4")
# L0Learn
good <- 10
plot(y, main = "L0Learn")
lines(X %*% L0fit$beta[[1]][,good], col = 2)
matlines(sweep(X, 2, L0fit$beta[[1]][,good], "*"), lty = 1, col = "green4")
# L0glm
good <- 31
plot(y, main = "L0glm")
lines(L0glmfit$fitted.values[,good], col = 2)
matlines(sweep(X, 2, L0glmfit$coefficients[,good], "*"), lty = 1, col = "green4")


##################################
#### ISSUE WITH SPARSE MATRIX ####
##################################


X.test <- matrix(runif(1000), ncol = 10)
X.test[sample(1:1000, 500)] <- 0
X.test <- Matrix(X.test, sparse = T)
coef.test <- runif(10)

t1 <- X.test %*% coef.test
t2 <- as.matrix(X.test) %*% coef.test
sum(abs(t1-t2))
sum(abs(X.test-as.matrix(X.test)))


###########################################################
#### TRY MATCHING BEST LAMBDA TO THEORETICAL IC LAMBDA ####
###########################################################


get.data.poisson.nn()
system.time(nnglm_pois_L0 <- L0glm(X = X, y = y, family = poisson(identity), intercept = FALSE,
                                   nonnegative = TRUE,
                                   tune.crit = "bic", tune.meth = "IC",
                                   lambda = 10^seq(4,6, length.out = 21), # Use arbitrary sequence
                                   control.fit = list(),
                                   control.iwls = list(thresh = 1E-7, rel.tol = 1E-4),
                                   control.l0 = list(maxit = 50)))
matplot(nnglm_pois_L0$lambda.tune$coefficients.lam, type = "l", log = "", lty = 1)
par(mfrow = c(2,2))
for(i in 1:ncol(nnglm_pois_L0$lam.opt$IC)){
  sub <- 1:18
  matplot(x = nnglm_pois_L0$lam.opt$lambdas[sub], y = nnglm_pois_L0$lam.opt$IC[sub,i], type = "l", main = colnames(nnglm_pois_L0$lam.opt$IC)[i], log = "x")
  abline(v = nnglm_pois_L0$lam.opt$best.lam)
  abline(v = nnglm_pois_L0$lam.opt$lambdas[which.min(nnglm_pois_L0$lam.opt$IC[,i])], col = "red")

}
plot.nnglm.benchmark(x = x, y = y, fit = nnglm_pois_L0, a.true = a,
                     main="Ground truth vs L0 penalized NNGLM estimates (BIC on full data)")


###########################################
#### TEST BEST LAMBDA SELECTION METHOD ####
###########################################

# Common constants to all simulations
n = 200
x = 1:n
# unknown peak heights
peakhrange = c(10,1E3)
# peak shape function
gauspeak = function(x, u, w, h=1) h*exp(((x-u)^2)/(-2*(w^2)))
# banded matrix with peak shape measured beforehand
X = do.call(cbind, lapply(1:n, function (u) gauspeak(x, u=u, w=5, h=1) ))

# Fitting
fam <- poisson("identity")
lambdas <- 10^seq(-8,2)

# Range of nonzero coefficients for simulation
K.range <- 1:25
# Generate simulated data for every K
set.seed(123)
res <- list()
for(K in K.range){
  # Simulate data
  # Unkown peak coefficients
  u = sample(x, K, replace=FALSE)
  h = 10^runif(K, min=log10(min(peakhrange)), max=log10(max(peakhrange)))
  a = rep(0, n)
  a[u] = h
  # Response variable
  y_nonoise = as.vector(X %*% a)
  y <- rpois(n, y_nonoise)

  # Split data in training (90%) and validation (10%)
  train <- sort(sample(1:n, n*0.9))
  X.train <- X[train,]
  y.train <- y[train]
  X.test <- X[-train,]
  y.test <- y[-train]

  # Fit data after optimizing lambda for L0 penalty
  # BIC
  print("BIC")
  cpglm <- pen.nnglm(X = X.train, y = y.train, family = fam, constr = "none",
                     lambda = lambdas, maxit.l0 = 1, cv = "full.data", cv.crit = "bic", inference = TRUE)
  glm.fit(x = X.train, y = y.train, family = fam, intercept = F)

  y.hat.test <- X.test %*% cpglm$coefficients
  BIC.mse <- sum((y.test - y.hat.test)^2)/length(y.test)
  BIC.p0 <- sum(cpglm$coefficients != 0)
  # AIC
  print("AIC")
  cpglm <- pen.nnglm(X = X.train, y = y.train, family = fam, constr = "none",
                     lambda = lambdas, maxit.l0 = 1, cv = "full.data", cv.crit = "aic", inference = TRUE)
  y.hat.test <- X.test %*% cpglm$coefficients
  AIC.mse <- sum((y.test - y.hat.test)^2)/length(y.test)
  AIC.p0 <- sum(cpglm$coefficients != 0)
  # LOOCV
  print("LOOCV")
  cpglm <- pen.nnglm(X = X.train, y = y.train, family = fam, constr = "none",
                     lambda = lambdas, maxit.l0 = 1, cv = "loocv", inference = TRUE)
  y.hat.test <- X.test %*% cpglm$coefficients
  LOOCV.mse <- sum((y.test - y.hat.test)^2)/length(y.test)
  LOOCV.p0 <- sum(cpglm$coefficients != 0)
  # loglik
  print("loglik")
  cpglm <- pen.nnglm(X = X.train, y = y.train, family = fam, constr = "none",
                     lambda = lambdas, maxit.l0 = 1, cv = "full.data", cv.crit = "loglik", inference = TRUE)
  y.hat.test <- X.test %*% cpglm$coefficients
  loglik.mse <- sum((y.test - y.hat.test)^2)/length(y.test)
  loglik.p0 <- sum(cpglm$coefficients != 0)

  # Store results
  res[K] <- list(MSE = c(BIC = BIC.mse, AIC = AIC.mse, LOOCV = LOOCV.mse, loglik = loglik.mse),
                 p0 = c(BIC = BIC.p0, AIC = AIC.p0, LOOCV = LOOCV.p0, loglik = loglik.p0))

  print.progress(K, max(K.range))
  cat("\n")
}
#




########################################
#### FIND OPTIMAL LAMBDA WITHOUT CV ####
########################################

# See Theorem 1 in Frommlet et al 2016.

# 1. Generate orthogonal matrix
library(pracma)
set.seed(1234)
n <- 200
p <- 200
X <- randortho(n)
X <- X * sqrt(n)
X <- X[,1:p]

# 2. Simulate coefficients and generate data with gaussian error structure
beta <- rep(0, p)
s <- 5
beta[1:s] <- 1
y0 <- X %*% beta
sigma2 <- 10
y <- y0 + rnorm(mean = 0, sd = sigma2, n = n)

# 3. Tune lambda
# theoetically, BIC should be optimized when lambda = log(n)/4
lams <- 10^(log10(log(n)/4) + seq(-4, 2, length.out = 150))
L0glm_fit <- L0glm(X = X, y = y, family = gaussian(identity),
                   intercept = FALSE, lambda = lams, # Use arbitrary sequence
                   tune.crit = "bic", tune.meth = "IC", seed = 123,
                   nonnegative = F,
                   control.fit = list(),
                   control.iwls = list(rel.tol = 1E-4, maxit = 50, thres = 1E-6),
                   control.l0 = list(rel.tol = 1E-4, maxit = 100, delta = 1E-6))

# 4. Check results
k <- rowSums(L0glm_fit$lambda.tune$coefficients.lam!=0)
coef.dif <- apply(L0glm_fit$lambda.tune$coefficients.lam, 1, function(x) sum(abs(x - beta)) )
rss <- L0glm_fit$lambda.tune$IC[,"rss"]
min2ll <- -2*L0glm_fit$lambda.tune$IC[,"loglik"]
par(mfrow = c(3,2))
for(crit in c("aic", "bic")){
  # AIC and BIC based on log likelihood
  # plot(x = L0glm_fit$lambda.tune$lambdas, y = L0glm_fit$lambda.tune$IC[,crit],
  #      log = "x", type = "l", main = crit)

  # Custom IC as defined in Frommlet et Nuel where BIC = 1/sima2*RSS + k * log(n)
  # sigma2 is known from the data generation (not the case in real world problems)
  plot(x = L0glm_fit$lambda.tune$lambdas, y = rss/sigma2 + k * ifelse(crit == "aic", 2, log(n)),
       log = "x", type = "l", main = crit)

  abline(v = lams[which.min(L0glm_fit$lambda.tune$IC[,crit])], col = 2)
  abline(v = 1/4 * ifelse(crit == "aic", 2, log(n)))
}
sigs <- apply(L0glm_fit$lambda.tune$coefficients.lam,  1, function(beta) sum((y - X %*% beta)^2))
sigs <- sigs/(n-k)
plot(x = L0glm_fit$lambda.tune$lambdas, y = rss/sigs, log = "x", type = "l", col = 4, main = "RSS/sig")
plot(x = L0glm_fit$lambda.tune$lambdas, y = min2ll, log = "x", type = "l", col = 4, main = "-2logLik")
plot(x = L0glm_fit$lambda.tune$lambdas, y = k, main = "Number of non zero coefficients", type = "l", log = "x")
abline(h = s, lty=2)
plot(x = L0glm_fit$lambda.tune$lambdas, y = coef.dif, main = "Coefficient distance to true coefficients", type = "l", log = "x", ylim = c(0, max(coef.dif)))
abline(h = 0, lty = 2)
abline(v = L0glm_fit$lambda.tune$lambdas[which.min(coef.dif)], col = 2)

matplot(x = L0glm_fit$lambda.tune$lambdas, y = L0glm_fit$lambda.tune$coefficients.lam, type = "l", lty = 1, log = "x")

par(mfrow = c(3,2))
plot(x = L0glm_fit$lambda.tune$lambdas, y = min2ll, log = "x", type = "l", col = 4, main = "-2logLik")
plot(x = L0glm_fit$lambda.tune$lambdas, y = rss, log = "x", type = "l", col = 4, main = "RSS")
plot(x = L0glm_fit$lambda.tune$lambdas, y = min2ll+2*k, log = "x", type = "l", col = 4, main = "AIC.ll")
abline(v = 0.5)
abline(v = L0glm_fit$lambda.tune$lambdas[which.min(min2ll+2*k)], col = 2)
plot(x = L0glm_fit$lambda.tune$lambdas, y = min2ll+log(n)*k, log = "x", type = "l", col = 4, main = "BIC.ll")
abline(v = log(n)/4)
abline(v = L0glm_fit$lambda.tune$lambdas[which.min(min2ll+log(n)*k)], col = 2)
plot(x = L0glm_fit$lambda.tune$lambdas, y = rss+2*k, log = "x", type = "l", col = 4, main = "AIC.rss")
abline(v = 0.5)
abline(v = L0glm_fit$lambda.tune$lambdas[which.min(rss+2*k)], col = 2)
plot(x = L0glm_fit$lambda.tune$lambdas, y = rss+log(n)*k, log = "x", type = "l", col = 4, main = "BIC.rss")
abline(v = log(n)/4)
abline(v = L0glm_fit$lambda.tune$lambdas[which.min(rss+log(n)*k)], col = 2)


# Poisson family

get.data.poisson.nn()
lams <- 10^(log10(log(n)/4) + seq(-2, 2, length.out = 50))
L0glm.pois <- L0glm(X = X, y = y, family = poisson(identity),
                    lambda = lams, intercept = FALSE,
                    tune.crit = "bic", tune.meth = "IC", seed = 123,
                    control.fit = control.fit(constr = "nonneg"),
                    control.iwls = control.iwls(rel.tol = 1E-4, maxit = 50, thres = 1E-6),
                    control.l0 = control.l0(rel.tol = 1E-4, maxit = 100, delta = 1E-6))

k <- rowSums(L0glm.pois$lambda.tune$coefficients.lam!=0)
coef.dif <- apply(L0glm.pois$lambda.tune$coefficients.lam, 1, function(x) sum(abs(x - beta)) )
rss <- L0glm.pois$lambda.tune$IC[,"rss"]
min2ll <- -2*L0glm.pois$lambda.tune$IC[,"loglik"]
par(mfrow = c(3,2))
for(crit in c("aic", "bic")){
  plot(x = L0glm.pois$lambda.tune$lambdas, y = L0glm.pois$lambda.tune$IC[,crit],
       log = "x", type = "l", main = crit)
  abline(v = lams[which.min(L0glm.pois$lambda.tune$IC[,crit])], col = 2)
  abline(v = 1/4 * ifelse(crit == "aic", 2, log(n)))
}
plot(x = L0glm.pois$lambda.tune$lambdas, y = (rss), log = "x", type = "l", col = 4, main = "RSS")
plot(x = L0glm.pois$lambda.tune$lambdas, y = min2ll, log = "x", type = "l", col = 4, main = "-2logLik")
plot(x = L0glm.pois$lambda.tune$lambdas, y = k, main = "Number of non zero coefficients", type = "l", log = "x")
abline(h = npeaks, lty=2)
plot(x = L0glm.pois$lambda.tune$lambdas, y = coef.dif, main = "Coefficient distance to true coefficients", type = "l", log = "x", ylim = c(0, max(coef.dif)))
abline(h = 0, lty = 2)
abline(v = L0glm.pois$lambda.tune$lambdas[which.min(coef.dif)], col = 2)


# Find maximum lambda
# https://stats.stackexchange.com/questions/416144/minimum-and-maximum-regularization-in-l0-pseudonorm-penalized-regression/417080#417080

library(L0Learn)
# Simulate some data
data <- GenSynthetic(n=100,p=500,k=30,seed=123)
X = data$X
y = data$y
# get.data.poisson.nn()
# make L0 penalized fit:
lams <- list()
lams[[1]] <- 10^seq(-9,9 ,length.out = 101)
L0fit <- L0Learn.fit(x=X, y=y, penalty="L0", algorithm="CD", maxSuppSize = 100, autoLambda = F, lambdaGrid = lams)
# Does not work...
unlist(L0fit$lambda)[unlist(L0fit$suppSize)==0][1] # = 0.0618124
unlist(L0fit$lambda)[unlist(L0fit$suppSize)==max(unlist(L0fit$suppSize))][1] # = 6.5916e-09

L0glm.res <- L0glm(X = X, y = y, family = gaussian(identity), start = NULL,
                   no.pen = 0, intercept = F, lambda = lams[[1]], control.l0 = control.l0(maxit = 100, rel.tol = 1E-9),
                   control.fit = control.fit(constr = "none"), control.iwls = control.iwls(), tune.meth = "none")
k <- colSums(L0glm.res$coefficients != 0)

max.th <- 1/2 * max(1/diag(crossprod(X, X)) * crossprod(X, y)^2)
max.emp <- lams[[1]][min(which(k == 0))]
max.th
max.emp
max.th/max.emp

# CCL: there is not an exact match, especially for our shifted gaussian case where
# the fold difference can be up to +/-1000. However the theoretically value seems
# to be sysematically higher than the empirical max lambdda. So it can still be
# use as a threshold to avoid testing too large lambdas !

plot(L0glm.res$coefficients[,86])






############################################################
#### TEST IF AR AND IWLS CAN BE OPTIIZED SIMULTANEOUSLY ####
############################################################

# Simulated poisson data
get.data.poisson.nn(Plot = T)

# Nested update
system.time(L0glm.nu <- L0glm(X = X, y = y, family = poisson(identity),
                              intercept = FALSE, lambda = 1,
                              control.fit = control.fit(constr = "nonneg"),
                              control.iwls = control.iwls(rel.tol = 1E-4, maxit = 50, thresh = 1E-4),
                              control.l0 = control.l0(rel.tol = 1E-4, maxit = 100)))
plot_L0glm_benchmark(x = x, y = y, fit = L0glm.nu, a.true = a,
                     main="Nested updates")


# Simultaneous update
system.time(L0glm.su <- L0glm(X = X, y = y, family = poisson(identity),
                              intercept = FALSE, lambda = 1,
                              control.fit = control.fit(constr = "nonneg"),
                              control.iwls = control.iwls(rel.tol = 1E-4, maxit = 1, thresh = 1E-4),
                              control.l0 = control.l0(rel.tol = 1E-4, maxit = 100)))
plot_L0glm_benchmark(x = x, y = y, fit = L0glm.su, a.true = a,
                     main="Nested updates")

# CCL: updating both AR and IWLS simultaneously lead to similar fit (maybe evn better)
# than the nested version and is +/- 5x faster (+/-5 x less iterations)

#################################################################
#### CHECK STABILITY OF THE ESTIMATES OVER A SERIE OF LAMBDA ####
#################################################################

# Simulated poisson data
get.data.poisson.nn(Plot = T)

lams <- 10^seq(-2, 2, length.out = 101)
L0glm <- L0glm(X = X, y = y, family = poisson(identity),
               intercept = FALSE, lambda = lams,
               cv = "full.data", cv.crit = "bic",
               control.fit = control.fit(constr = "nonneg"),
               control.iwls = control.iwls(rel.tol = 1E-4, maxit = 1, thresh = 1E-4),
               control.l0 = control.l0(rel.tol = 1E-4, maxit = 100))
plot_L0glm_benchmark(x = x, y = y, fit = L0glm, a.true = a,
                     main="Nested updates")

matplot(L0glm$lambda.tune$coefficients.lam, type = "l", lty = 1)



############################
#### TEST ON LARGE DATA ####
############################


# Use a random GC-MS sample
source("D:/Documents/Dropbox/christophe/_tests deconvolution scripts/utils-0.0.2.R")
M <- read.cdf("D:/Documents/beer_chromatograms_miguel/20170129/54085015.cdf")
y <- rowSums(M)
n <- length(y)
peak.widths <- fit.average.peak.width(y = y,win.size = 300,
                                      peak.width.range = c(1,20))
peak.width.fun <- peak.widths$peak.width.fun
peak.widths <- peak.widths$peak.width.fun$fitted.values
peak.widths.preds <- predict(peak.width.fun, newdata = list(x = 1:length(y)))
X <- build.banded.gaus.mat(w = peak.widths.preds, tol = 1E-5)

deconv.fit <- L0glm(X = X, y = y, family = poisson(identity), intercept = F,
                    lambda = 1,
                    control.l0 = control.l0(maxit = 25),
                    control.iwls = control.iwls(thresh = 1, maxit = 1),
                    control.fit = control.fit(maxit = 25, block.size = 100, constr = "nonneg"))

xlims <- c(4600, 5100)
plot(x = 1:n, y = y, xlim = xlims, type = "l", log = "y")
lines(x = 1:n, y = deconv.fit$fitted.values, col = "orange2")
matlines(x = 1:n, y = sweep(X[,xlims[1]:xlims[2]], 2, deconv.fit$coefficients[xlims[1]:xlims[2]], "*"), lty = 1, col = "green4")

x <- 1:n
bas <- bs(x = 1:n, df = 15+3, degree = 3)
class(bas) <- "matrix"
X.full <- cbind(bas, X)

deconv.fit <- L0glm(X = X.full, y = y, family = poisson(identity), intercept = F,
                    lambda = 1, no.pen = 1:ncol(bas),
                    control.l0 = control.l0(maxit = 25),
                    control.iwls = control.iwls(thresh = 1, maxit = 1),
                    control.fit = control.fit(maxit = 25, block.size = 100, constr = "nonneg"))

xlims <- c(4600, 5100)
plot(x = 1:n, y = y, xlim = xlims, type = "l", log = "y", ylim = c(1, max(y)))
lines(x = 1:n, y = deconv.fit$fitted.values, col = "orange2")
matlines(x = 1:n, y = sweep(X.full[,xlims[1]:xlims[2]], 2, deconv.fit$coefficients[xlims[1]:xlims[2]], "*"), lty = 1, col = "green4")


####################
#### DEBUG ZONE ####
####################





unlist(L0fit$lambda)[unlist(L0fit$suppSize)==0][1] # = 0.0618124
unlist(L0fit$lambda)[unlist(L0fit$suppSize)==max(unlist(L0fit$suppSize))][1] # = 6.5916e-09





max(unlist(L0fit$suppSize)) # size of largest model = 100
max(diag(1/crossprod(X, X)) * (crossprod(X, y)^2)) # 677.1252


