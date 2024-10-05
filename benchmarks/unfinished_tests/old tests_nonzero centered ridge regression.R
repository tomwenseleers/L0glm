# Tests of L0glm package zero vs nonzero (beta_old) centered ridge regression ####

library(L0glm)
library(export)

n <- p <- 10000 # to check performance use e.g. n=p=10000, otherwise use e.g. n=p=100-300
s <- 0.1 # sparsity, nr of timepoints with a peak
k <- round(n*s)
mean_beta <- 2000
sd_logbeta <- 1
sd_noise <- 5 # SD of noise for Gaussian noise case
crit <- "bic" # choose lambda based on what would be expected to optimise this IC
# options in order of increasing regularisation: aic, bic, gic, ebic, mbic

# 1. SINGLE CHANNEL DECONVOLUTION ####
## 1.1 Single channel deconvolution assuming Poisson noise ####
set.seed(1)
sim <- simulate_spike_train(n = n, 
                            p = p,
                            k = k,
                            family = "poisson", 
                            mean_beta = mean_beta,
                            sd_logbeta = sd_logbeta)
X <- sim$X
y <- sim$y

# # using lambda 1
# L0glm.out <- L0glm.fit(X, y, 
#                            family = poisson(identity), 
#                            lambda = 1, 
#                            nonnegative = TRUE)
# # points(sim$x, L0glm.out, pch = 16, col= 'blue')

# using lambda expected to optimise given IC
lambda_preset = preset.lambda(IC = crit, # options aic, bic, gic, ebic, mbic 
                              y = y, 
                              X = X, 
                              family = poisson(identity)) 
lambda_preset # lambda = 1.15 would be expected to optimise BIC with n=p=100, lambda = 1.76 would optimise GIC

# couple of possible coefficient initialisations to use as starting values in L0glm.fit which could be used for start ####
weights = 1/(sim$y+0.1) # L0glm_gaus_presetlambda$weights
Xw = Matrix(Diagonal(x=sqrt(weights)) %*% sim$X, sparse=T) # = X * sqrt(weights)
Xwaug = Matrix(rbind(Xw, Diagonal(x=sqrt(rep(lambda_preset,ncol(sim$X))))), sparse=T) # row augmented covariate matrix
ywaug = as.vector(c(sim$y, rep(0, ncol(sim$X)))) # 

# orthogonal initialisation (assuming variables are orthogonal) ####
beta_init_orthogonal <- as.vector(crossprod(Xwaug, ywaug) / colSums(Xwaug^2)) 
cor(beta_init_orthogonal, as.vector(sim$beta_true)) 
plot(beta_init_orthogonal, as.vector(sim$beta_true), pch=16, col="blue")
# colSums(Xwaug^2) = diag(crossprod(Xwaug)))
# PS original variables in columns of X can be very far from orthogonal, 
# but the augmentation with a diagonal matrix with sqrt(lambdas) along the diagonal
# makes the assumption hold better

# diagonal initialisation (assuming system is dominated by diagonal of X) ####
# for signal deconvolution problem with n=p this gives the best initial solution
# only works if X is square though
# beta_init_diagonal <- y*sqrt(weights) / diag(Xw)


library(microbenchmark)
system.time(L0glm.out <- L0glm.fit(X, y, 
                           family = poisson(identity), 
                           lambda = lambda_preset,
                           nonzero.centered = FALSE, # gives slightly better solution based on IC, about the same speed or slightly faster
                           ridge.center = NULL, 
                           nonnegative = TRUE,
                           start = rep(1, ncol(X)),
                           relaxed = FALSE, 
                           thresh = 1E-10,                           
                           warmupiters = 1,
                           warmupiters.nzc = 10)) 
# 0.34s
points(sim$x, L0glm.out$coefficients, pch = 16, col= 'blue')
# graph2pdf(file="..//benchmarks/plots/1_1 single channel deconvolution Poisson.pdf", width=7, height=5)

# some metrics
do.call(cbind, get.metrics(beta = L0glm.out$coefficients, 
                           X = X, y = y, 
                           family = poisson(identity), 
                           beta_true = sim$beta_true))
plot(x=sim$beta_true, y=L0glm.out, 
     xlab="Estimated coefficients", ylab="True coefficients",
     main="L0glm with Poisson noise",
     pch = 16, col= 'blue', type="p")

# solution over a sequence of lambda values



## 1.2 Single channel deconvolution assuming Gaussian noise ####
set.seed(1)
sim <- simulate_spike_train(n = n, 
                            p = p,
                            k = k,
                            family = "gaussian", 
                            mean_beta = mean_beta,
                            sd_logbeta = sd_logbeta,
                            sd_noise = sd_noise)
X <- sim$X
y <- sim$y

# using lambda expected to optimise given IC
lambda_preset = preset.lambda(IC = crit, # options aic, bic, gic, ebic, mbic 
                              y = y, 
                              X = X, 
                              family = gaussian(identity)) 
lambda_preset # lambda that would be expected to optimise given IC

system.time(L0glm.out <- L0glm.fit(X, y, 
                           family = gaussian(identity), 
                           lambda = lambda_preset, 
                           nonzero.centered = TRUE, # gives slightly better solution based on IC, about the same speed or slightly faster
                           ridge.center = NULL, 
                           nonnegative = TRUE,
                           relaxed = FALSE, # TODO: CHECK: relaxed argument has no effect in single channel fit - is this as intended? Should one not be allowed to fit an L0glm and then refit all selected variables unpenalized with lambda=0? 
                           warmupiters = 1)) # should be set to 1 for gaussian
# 0.18s
points(sim$x, L0glm.out, pch = 16, col= 'blue')
# graph2pdf(file="..//benchmarks/plots/1_2 single channel deconvolution Gaussian.pdf", width=7, height=5)

# some metrics
do.call(cbind, get.metrics(beta = L0glm.out, 
                           X = X, y = y, 
                           family = poisson(identity), 
                           beta_true = sim$beta_true))
plot(x=sim$beta_true, y=L0glm.out, 
     xlab="Estimated coefficients", ylab="True coefficients",
     main="L0glm with Gaussian noise",
     pch = 16, col= 'blue', type="p")




# 2. MULTICHANNEL DECONVOLUTION ####
n = 100
p = 100
k = ifelse(round(n/10)<=100, round(n/10), 100) # nr of peaks

## 2.1 Multichannel deconvolution assuming Poisson noise ####
set.seed(1)
sim <- simulate_spike_train(n = n, 
                            p = p,
                            k = k,
                            family = "poisson", 
                            mean_beta = mean_beta,
                            sd_logbeta = sd_logbeta,
                            multichannel = TRUE)
X <- sim$X
y <- sim$y

# using lambda expected to optimise given IC (lambda = 1.15)
lambda_preset = preset.lambda(IC = crit, # options aic, bic, gic, ebic, mbic 
                              y = y, 
                              X = X, 
                              family = poisson(identity)) 
lambda_preset # lambda = 1.15 would be expected to optimise BIC, lambda = 1.76 would optimise GIC

system.time(L0glm.out <- L0glm.fit(X, y, 
                           family = poisson(identity), 
                           lambda = lambda_preset,
                           nonzero.centered = FALSE, # gives slightly better solution based on IC, about the same speed or slightly faster
                           ridge.center = NULL,
                           nonnegative = TRUE,
                           grouped = TRUE,
                           group.maxit = 20,
                           relaxed = FALSE,
                           warmupiters = 4))
# 30s for n=p=10000, 12.64s for n=p=100
abline(v=which(rowMin(L0glm.out)!=0), col = 'cyan')
# graph2pdf(file="..//benchmarks/plots/2_1 multichannel deconvolution Poisson.pdf", width=7, height=5)

# NOTE: nonzero.centered = FALSE gives better result here, though for single channel
# data result with nonzero.centered = TRUE was better

# some metrics
do.call(cbind, get.metrics(beta = L0glm.out, 
                           X = X, y = y, 
                           family = poisson(identity), 
                           beta_true = sim$beta_true))
sum(L0glm.out!=0) # 583
sum(sim$beta_true!=0) # 793
length(intersect(which(L0glm.out!=0), which(sim$beta_true!=0))) # 564
plot(x=sim$beta_true, y=L0glm.out, 
     xlab="Estimated coefficients", ylab="True coefficients",
     main="L0 group glm with Poisson noise",
     pch = 16, col= 'blue', type="p")



## 2.2 Multichannel deconvolution assuming Gaussian noise ####
sim <- simulate_spike_train(n = n, 
                            p = p,
                            k = k,
                            family = "gaussian", 
                            mean_beta = mean_beta,
                            sd_logbeta = sd_logbeta,
                            sd_noise = sd_noise, 
                            multichannel = TRUE)
X <- sim$X
y <- sim$y

# using lambda expected to optimise given IC
lambda_preset = preset.lambda(IC = crit, # options aic, bic, gic, ebic, mbic 
                              y = y, 
                              X = X,  
                              family = gaussian(identity)) 
lambda_preset # lambda that would be expected to optimise given IC

system.time(L0glm.out <- L0glm.fit(X, y, 
                           family = gaussian(identity), 
                           lambda = lambda_preset, 
                           nonzero.centered = FALSE, # gives slightly better solution based on IC, about the same speed or slightly faster
                           ridge.center = NULL,
                           nonnegative = TRUE,
                           grouped = TRUE,
                           group.maxit = 20,
                           relaxed = FALSE,
                           warmupiters = 1)) # TODO: FOR GAUSSIAN GIVES BEST SOLUTION
# XXs for n=p=10000, 25.87s for n=p=100
abline(v=which(rowMin(L0glm.out)!=0), col = 'cyan')
# graph2pdf(file="..//benchmarks/plots/2_2 multichannel deconvolution Gaussian.pdf", width=7, height=5)

# NOTE: nonzero.centered = FALSE gives better result here, though for single channel
# data result with nonzero.centered = TRUE was better

# some metrics
do.call(cbind, get.metrics(beta = L0glm.out, 
                           X = X, y = y, 
                           family = gaussian(identity), 
                           beta_true = sim$beta_true))
sum(L0glm.out!=0) # 604
sum(sim$beta_true!=0) # 821
length(intersect(which(L0glm.out!=0), which(sim$beta_true!=0))) # 491
plot(x=sim$beta_true, y=L0glm.out, 
     xlab="Estimated coefficients", ylab="True coefficients",
     main="L0 group glm with Gaussian noise",
     pch = 16, col= 'blue', type="p")
