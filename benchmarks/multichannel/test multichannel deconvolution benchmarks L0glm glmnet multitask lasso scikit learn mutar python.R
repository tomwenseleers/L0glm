source(if (file.exists("benchmarks/cv_glmnet_parallel.R")) "benchmarks/cv_glmnet_parallel.R" else "../benchmarks/cv_glmnet_parallel.R")
# TEST MULTITASK FEATURE SELECTION METHODS IN glmnet, scikitlearn (MultiTaskLasso=multivariate group LASSO) & mutar ####
# library(remotes)
# remotes::install_github("tomwenseleers/L0glm/L0glm")

devtools::load_all(".")

# simulate blurred multichannel spike train
set.seed(1)
s <- 0.1 # sparsity (% of timepoints where there is a peak)
p <- 500 # 500 variables
# simulate multichannel blurred spike train with Gaussian noise
sd_noise <- 1
sim <- simulate_spike_train(n=p, p=p, k=round(s*p), 
                            mean_beta = 10000,
                            sd_logbeta = 1,
                            family="gaussian",
                            sd_noise = sd_noise,
                            multichannel=TRUE, sparse=TRUE)
X <- sim$X # n x p matrix
Y <- sim$y # n x m matrix  
colnames(X) = paste0("x", 1:ncol(X)) # needed for abess
colnames(Y) = paste0("y", 1:ncol(Y))
true_coefs <- sim$beta_true
m <- ncol(Y) # nr of tasks
n <- nrow(X) # nr of observations
p <- ncol(X) # nr of independent variables (shifted copies of point spread functions)
W <- 1/(Y+0.1) # approx 1/variance Poisson observation weights with family="poisson", n x m matrix

# Create an array with dimensions m x n x p with weighted versions of X
library(abind)
X_weighted <- lapply(seq_len(m), function (channel) as.matrix(X * sqrt(W[, channel])))
X_weighted <- abind(X_weighted, along = 3)
X_weighted <- aperm(X_weighted, perm=c(3,1,2))
str(X_weighted) # m x n x p array

# Create weighted version of Y
Y_weighted <- Y * sqrt(W)

# to call scikitlearn & mutar
# install.packages("reticulate")
library(reticulate)
# use the correct version of python
use_python("C:/ProgramData/Anaconda3")

# import required python modules
np <- import("numpy")
sklearn <- import("sklearn.linear_model")
scipy_sparse <- import("scipy.sparse")
mutar <- import("mutar")

# convert  R objects to Python objects before using them with Python functions
scipy_sparse <- import("scipy.sparse")
X_py_sparse <- r_to_py(X) # scipy.sparse.csc.csc_matrix
X_py_dense <- r_to_py(as.matrix(X)) # X_py_sparse$toarray() # numpy.ndarray
Y_py <- r_to_py(as.matrix(Y))

# using relaxed group L0 GLM ####
library(parallel)
library(matrixStats)
library(qlcMatrix)
system.time(beta_relaxedgroupL0 <- L0glm.fit_new(X = sim$X, y = sim$y, 
                                     family = gaussian(identity), # poisson(identity),
                                     lambda = 1, # needs to be larger with larger sd_noise for family=gaussian, e.g. 0.5 with sd_noise=1 & 5 with sd_noise=5 
                                     lower = rep(0, ncol(X)),
                                     upper = rep(Inf, ncol(X)),
                                     nonnegative = TRUE,
                                     thresh = 0.01,
                                     warmupiters = 5,
                                     solver = "eigen",
                                     method = 1L, # using Eigen's sparse Cholesky solver
                                     parallellize = T, 
                                     relaxedGroupL0 = T)) # 10.4s
# TODO: osqp still gives errors if I set parallellize to T - then gives
# Error in checkForRemoteErrors(val) : 
# 8 nodes produced errors; first error: length(l) == m is not TRUE

beta_relaxedgroupL0 <- as.matrix(beta_relaxedgroupL0) # convert to dense matrix
# beta_relaxedgroupL0[beta_relaxedgroupL0>1E7] = 0 # TODO: check- osqp gives a couple of outlier coefficients

image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=as.matrix(sim$beta_true)^0.1, 
      col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="True coefficients")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim$X))[as.vector(rowMin(beta_relaxedgroupL0)!=0)], col="cyan")

image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=beta_relaxedgroupL0^0.1, 
      col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="nonnegative relaxed group L0")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim$X))[as.vector(rowMin(beta_relaxedgroupL0)!=0)], col="cyan")
sum(rowMax(sim$beta_true)!=0) # 50 true peaks
sum(rowMin(beta_relaxedgroupL0)!=0) # 54 peaks detected
length(intersect(which(rowMax(sim$beta_true)!=0), 
                 which(rowMin(beta_relaxedgroupL0)!=0))) # 47 out of 50 peaks correctly detected
length(which(sim$beta_true!=0)) # 3758 coefs nonzero
length(as.vector(sim$beta_true)) # 80500 coefs in total
length(intersect(which(sim$beta_true!=0), 
                 which(beta_relaxedgroupL0!=0))) # 2188 true pos coefs
length(intersect(which(sim$beta_true==0), 
                 which(beta_relaxedgroupL0!=0))) # 487 false pos coefs
plot(as.vector(sim$beta_true), as.vector(beta_relaxedgroupL0), pch=16)
cor(as.vector(sim$beta_true), as.vector(beta_relaxedgroupL0)) # 0.9839155


# MULTICHANNEL DECONVOLUTION USING nonnegative group LASSO with glmnet ####
library(glmnet)
system.time(cvfit <- cv_glmnet_parallel(X,
                               Y, 
                               # weights = 1/(rowSums(Y)+0.1), # does not allow weights to be passed as a matrix, i.e. 1/(Y+0.1)
                               alpha = 1,  # lasso
                               # gamma = 0, # set when relax = TRUE
                               # lambda = 10^seq(-3,2,length.out=20),
                               # dfmax = round(s*p),
                               # pmax = round(s*p),
                               family = "mgaussian", 
                               nlambda = 100,
                               # lambda.min.ratio = 1E-4,
                               # thresh = 1e-10,
                               nfolds = 5, 
                               standardize = FALSE,
                               standardize.response = FALSE, # set to TRUE solution is not sparse enough
                               intercept = FALSE, 
                               relax = FALSE,
                               lower.limits = rep(0, ncol(X)+1))) # 26.71s

plot(cvfit)

# Fit the glmnet model over a sequence of lambda values
system.time(fit_glmnet <- glmnet(X, 
                     Y,
                     # weights = 1/(rowSums(Y)+0.1), # does not allow weights to be passed as a matrix, i.e. 1/(Y+0.1)
                     alpha = 1, # lasso
                     # lambda = 10^seq(-3,2,length.out=20),
                     # dfmax = round(s*p),
                     # pmax = round(s*p),
                     family = "mgaussian", 
                     nlambda = 100,
                     # lambda.min.ratio = 1E-4,
                     # thresh = 1e-10,
                     standardize = FALSE,
                     standardize.response = FALSE, # set to TRUE solution is not sparse enough
                     intercept = FALSE, 
                     relax = FALSE, 
                     lower.limits = rep(0, ncol(X)+1))) # 8.55s

# Select the best lambda value (I am using cvfit$lambda.1se to get slightly sparser solution)
cvfit$lambda.1se # 2.954027
cvfit$lambda.min # 1.855215
best_lambda <- cvfit$lambda.1se

plot(fit_glmnet, xvar="lambda")

# get the coefficients for each task for the best lambda value
coefs <- coef(fit_glmnet, s = best_lambda)
beta_mgaussian_glmnet <- do.call(cbind, lapply(seq_len(m), function (channel) as.matrix(coefs[[channel]][-1,,drop=F])))
colnames(beta_mgaussian_glmnet) = 1:ncol(beta_mgaussian_glmnet)
beta_mgaussian_glmnet[abs(beta_mgaussian_glmnet)<0.01] = 0 # slight amount of thresholding of small coefficients

image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=beta_mgaussian_glmnet^0.01, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="nonnegative group Lasso glmnet (red=true support)")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")
# abline(v=(1:nrow(sim$X))[as.vector(rowMin(beta_mgaussian_glmnet)!=0)], col="cyan")
sum(rowMax(sim$beta_true)!=0) # 50 true peaks
sum(rowMin(beta_mgaussian_glmnet)!=0) # 71 peaks detected
length(intersect(which(rowMax(sim$beta_true)!=0), 
                 which(rowMin(beta_mgaussian_glmnet)!=0))) # 48 out of 50 peaks detected
length(which(sim$beta_true!=0)) # 3758 coefs nonzero
length(as.vector(sim$beta_true)) # 80500 coefs in total
length(intersect(which(sim$beta_true!=0), 
                 which(beta_mgaussian_glmnet!=0))) # 3542 true pos coefs
length(intersect(which(sim$beta_true==0), 
                 which(beta_mgaussian_glmnet!=0))) # 6429 false pos coefs
plot(as.vector(sim$beta_true), as.vector(beta_mgaussian_glmnet), pch=16)
cor(as.vector(sim$beta_true), as.vector(beta_mgaussian_glmnet)) # 0.9812194



# MULTITASK LEARNING USING abess PACKAGE ####
# see https://abess-team.github.io/abess/articles/v06-MultiTaskLearning.html
# https://abess-team.github.io/abess/reference/abess-package.html
# check my question here https://github.com/abess-team/abess/issues/510

library(abess)
# sim_dat <- generate.data(n = 1000, p = 1000, y.dim=100, support.size=10, family="mgaussian")
# abess(x = sim_dat[["x"]], y = sim_dat[["y"]], family="mgaussian")
# note: abess cannot take into account nonnegativity constraints
# because of that it is best to pre-select variables using nonnegative group LASSO (using glmnet)
# and then do best subset selection among those variables
# otherwise abess ends up with a model with a support size of 123 (using ebic)
preselvars = which(rowMin(beta_mgaussian_glmnet)!=0) # pre-select variables based on nonnegative mgaussian group LASSO
system.time(abess_fit <- abess(x = X[,preselvars,drop=F], 
                               y = Y, 
                               # weights = 1/(Y+0.1), 
                               family = "mgaussian",
                               # normalize = 3, # or 1?
                               tune.path = "sequence", # gsection", # or sequence
                               support.size = c(1:length(preselvars)),
                               # gs.range = c(1,200),
                               lambda=0, 
                               # init.active.set=,  # vector of integers
                               # splicing.type = 2,
                               max.splicing.iter = 20,
                               warm.start = TRUE,
                               tune.type = "ebic") # or cv or ebic or bic or aic
            ) # 0.29s
# abess_fit
extract(abess_fit)$support.size # 55
plot(abess_fit, type="tune")
            
# extract(abess_fit)
beta_abess = beta_mgaussian_glmnet
beta_abess[preselvars,] = as.matrix(extract(abess_fit)$beta) # coefficient matrix for best subset
beta_abess[beta_abess < 0] = 0
image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=beta_abess^0.1, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="abess multitask mgaussian (red=true peaks)")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim$X))[as.vector(rowMin(beta_abess)!=0)], col="cyan")
sum(rowMax(sim$beta_true)!=0) # 50 true peaks
sum(rowMin(beta_abess)!=0) # 55 peaks detected
length(intersect(which(rowMax(sim$beta_true)!=0), 
                 which(rowMin(beta_abess)!=0))) # 45 out of 50 peaks detected
length(which(sim$beta_true!=0)) # 3758 coefs nonzero
length(as.vector(sim$beta_true)) # 80500 coefs in total
length(intersect(which(sim$beta_true!=0), 
                 which(beta_abess!=0))) # 3216 true pos coefs
length(intersect(which(sim$beta_true==0), 
                 which(beta_abess!=0))) # 2867 false pos coefs
plot(as.vector(sim$beta_true), as.vector(beta_abess), pch=16)
cor(as.vector(sim$beta_true), as.vector(beta_abess)) # 0.9861599




# # MULTICHANNEL DECONVOLUTION USING group LASSO with grpreg ####
# # takes forever to run or crashes
# library(grpreg)
# system.time(cvfit_grpreg <- cv.grpreg(as.matrix(X), Y))
# # Select the best lambda value
# best_lambda <- cvfit_grpreg$lambda.min
# system.time(fit_grpreg <- grpreg(as.matrix(X), Y))
# beta_grpreg <- coef(fit_grpreg, lambda=best_lambda)
# 
# # also try https://search.r-project.org/CRAN/refmans/gglasso/html/gglasso.html
# and https://search.r-project.org/CRAN/refmans/RMTL/html/MTL.html



# MULTICHANNEL DECONVOLUTION USING MultiTaskLasso IN scikitlearn (using L1/L2 block norm penalty = multivariate group LASSO of https://projecteuclid.org/download/pdfview_1/euclid.aos/1291388368) ####

# define the MultiTaskLasso model
multitasklasso <- sklearn$MultiTaskLasso(alpha=10, max_iter=10000) # 1.8s
# NOTE alpha was chosen so that selected nr of variables matches with true nr
# check what IC would work to select this alpha

# specify nonnegativity constraints - not possible for the moment?
# multitasklasso$positive <- TRUE

# set observation weights - not possible for the moment?
# multitasklasso$sample_weights = W

# fit the model - note: currently does not support observation weights
system.time(multitasklasso$fit(X_py_dense, Y_py)) # 1.8s

# get coefficients
beta_MultiTaskLasso <- t(multitasklasso$coef_)

image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=beta_MultiTaskLasso^1E-10, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="scikitlearn multitask LASSO coefficients")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim$X))[as.vector(rowMin(beta_MultiTaskLasso)!=0)], col="cyan")
sum(rowMax(sim$beta_true)!=0) # 50 true peaks
sum(rowMin(beta_MultiTaskLasso)!=0) # 51 peaks detected with this alpha
length(intersect(which(rowMax(sim$beta_true)!=0), 
                 which(rowMin(beta_MultiTaskLasso)!=0))) # 36 peaks detected
length(which(sim$beta_true!=0)) # 3758 coefs nonzero
length(as.vector(sim$beta_true)) # 80500 coefs in total
length(intersect(which(sim$beta_true!=0), 
                 which(beta_MultiTaskLasso!=0))) # 2427 true pos coefs
length(intersect(which(sim$beta_true==0), 
                 which(beta_MultiTaskLasso!=0))) # 5784 false pos coefs
plot(as.vector(sim$beta_true), as.vector(beta_MultiTaskLasso), pch=16)
cor(as.vector(sim$beta_true), as.vector(beta_MultiTaskLasso)) # 0.9590041



# MULTICHANNEL DECONVOLUTION USING methods in mutar package https://hichamjanati.github.io/mutar/#, https://hichamjanati.github.io/mutar/auto_examples/index.html ####
# in Anaconda Powershell run pip install -U mutar

# define the different mutar models
indlasso <- mutar$IndLasso(alpha = rep(1, m), fit_intercept=FALSE, normalize=FALSE)
indrewlasso <- mutar$IndRewLasso(alpha = rep(1, m), fit_intercept=FALSE)
grouplasso <- mutar$GroupLasso(alpha = 1)
dirtymodel <- mutar$DirtyModel(alpha = 10, beta = 5)
multilevellasso <- mutar$MultiLevelLasso(alpha = 10, fit_intercept=FALSE)

# convert  R objects to Python objects before using them with Python functions
X_py_sparse <- r_to_py(X) # scipy.sparse.csc.csc_matrix
X_py_dense <- r_to_py(as.matrix(X)) # X_py_sparse$toarray() # numpy.ndarray
class(X_py_dense) # numpy.ndarray, 200x200
Y_py <- np$array(t(as.matrix(Y))) # needs to be transposed
dim(Y_py) # 142 200
class(Y_py) # matrix array


# Now you have to reshape X array to match the required shape.
# We use np.repeat function to create copies of your X matrix for each task/channel
X_py_expanded <- np$expand_dims(X_py_dense, 0L)
dim(X_py_expanded) # 1 200 200
class(X_py_expanded) # array
X_reshaped <- np$tile(X_py_expanded, c(m, 1L, 1L))
dim(X_reshaped) # 142 200 200
class(X_reshaped) # array

X_weighted_py <- np$array(X_weighted)
dim(X_weighted_py) # 142 200 200
class(X_weighted_py) # array
Y_weighted_py <- np$array(t(as.matrix(Y_weighted))) # needs to be transposed
dim(Y_weighted_py) # 142 200
class(Y_weighted_py) # matrix array


# fit the models - note: currently methods do not support observation weights
system.time(indlasso$fit(X_reshaped, Y_py)) # 1.97s
system.time(indrewlasso$fit(X_reshaped, Y_py)) # 29s 
system.time(grouplasso$fit(X_reshaped, Y_py)) # 123.95s
system.time(dirtymodel$fit(X_reshaped, Y_py)) # 26.97s
system.time(multilevellasso$fit(X_reshaped, Y_py)) # 491s

# indlasso$fit(X_weighted_py, Y_weighted_py)
# indrewlasso$fit(X_weighted_py, Y_weighted_py)
# grouplasso$fit(X_weighted_py, Y_weighted_py)
# dirtymodel$fit(X_weighted_py, Y_weighted_py)
# multilevellasso$fit(X_weighted_py, Y_weighted_py)

# get coefficients
beta_IndLasso <- indlasso$coef_
beta_IndRewLasso <- indrewlasso$coef_
beta_GroupLasso <- grouplasso$coef_
beta_DirtyModel <- dirtymodel$coef_
beta_MultiLevelLasso <- multilevellasso$coef_

# Individual Lasso
image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=beta_IndLasso^0.1, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="mutar Individual Lasso coefficients")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim$X))[as.vector(rowMin(beta_IndLasso)!=0)], col="cyan")
sum(rowMax(sim$beta_true)!=0) # 50 true peaks
sum(rowMin(beta_IndLasso)!=0) # 103 peaks detected - still too high
length(intersect(which(rowMax(sim$beta_true)!=0), 
                 which(rowMin(beta_IndLasso)!=0))) # 49 out of 50 peaks detected
length(which(sim$beta_true!=0)) # 3758 coefs nonzero
length(as.vector(sim$beta_true)) # 80500 coefs in total
length(intersect(which(sim$beta_true!=0), 
                 which(beta_IndLasso!=0))) # 973 true pos coefs
length(intersect(which(sim$beta_true==0), 
                 which(beta_IndLasso!=0))) # 434 false pos coefs


# Individual Reweighted Lasso
image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=beta_IndRewLasso^0.1, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="mutar Individual Reweighted Lasso coefficients")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim$X))[as.vector(rowMin(beta_IndRewLasso)!=0)], col="cyan")
sum(rowMax(sim$beta_true)!=0) # 50 true peaks
sum(rowMin(beta_IndRewLasso)!=0) # 92 peaks detected - still too high
length(intersect(which(rowMax(sim$beta_true)!=0), 
                 which(rowMin(beta_IndRewLasso)!=0))) # 49 peaks detected
length(which(sim$beta_true!=0)) # 3758 coefs nonzero
length(as.vector(sim$beta_true)) # 80500 coefs in total
length(intersect(which(sim$beta_true!=0), 
                 which(beta_IndRewLasso!=0))) # 913 true pos coefs
length(intersect(which(sim$beta_true==0), 
                 which(beta_IndRewLasso!=0))) # 200 false pos coefs


# Group Lasso
image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=beta_GroupLasso^0.1, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="mutar Group Lasso coefficients")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim$X))[as.vector(rowMin(beta_GroupLasso)!=0)], col="cyan")
sum(rowMax(sim$beta_true)!=0) # 50 true peaks
sum(rowMin(beta_GroupLasso)!=0) # 86 peaks detected - still too high
length(intersect(which(rowMax(sim$beta_true)!=0), 
                 which(rowMin(beta_GroupLasso)!=0))) # 50 peaks detected
length(which(sim$beta_true!=0)) # 3758 coefs nonzero
length(as.vector(sim$beta_true)) # 80500 coefs in total
length(intersect(which(sim$beta_true!=0), 
                 which(beta_GroupLasso!=0))) # 3758 true pos coefs
length(intersect(which(sim$beta_true==0), 
                 which(beta_GroupLasso!=0))) # 10088 false pos coefs

# Dirty Model
image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=beta_DirtyModel^0.1, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="mutar Dirty Model coefficients")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim$X))[as.vector(rowMin(beta_DirtyModel)!=0)], col="cyan")
sum(rowMax(sim$beta_true)!=0) # 50 true peaks
sum(rowMin(beta_DirtyModel)!=0) # 66 peaks detected - still too high
length(intersect(which(rowMax(sim$beta_true)!=0), 
                 which(rowMin(beta_DirtyModel)!=0))) # 36 peaks detected
length(which(sim$beta_true!=0)) # 3758 coefs nonzero
length(as.vector(sim$beta_true)) # 80500 coefs in total
length(intersect(which(sim$beta_true!=0), 
                 which(beta_DirtyModel!=0))) # 1839 true pos coefs
length(intersect(which(sim$beta_true==0), 
                 which(beta_DirtyModel!=0))) # 4633 false pos coefs

# MultiLevel Lasso
image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=beta_MultiLevelLasso^1E-10, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="L0 nonneg Poisson GLM coefficients")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim$X))[as.vector(rowMin(beta_MultiLevelLasso)!=0)], col="cyan")
sum(rowMax(sim$beta_true)!=0) # 50 true peaks
sum(rowMin(beta_MultiLevelLasso)!=0) # 34 peaks detected - too low
length(intersect(which(rowMax(sim$beta_true)!=0), 
                 which(rowMin(beta_MultiLevelLasso)!=0))) # 27 peaks detected
length(which(sim$beta_true!=0)) # 3758 coefs nonzero
length(as.vector(sim$beta_true)) # 80500 coefs in total
length(intersect(which(sim$beta_true!=0), 
                 which(beta_MultiLevelLasso!=0))) # 1406 true pos coefs
length(intersect(which(sim$beta_true==0), 
                 which(beta_MultiLevelLasso!=0))) # 654 false pos coefs



