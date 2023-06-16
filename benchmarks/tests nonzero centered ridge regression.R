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
# L0glm.out <- L0glm.fit_new(X, y, 
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

library(microbenchmark)
system.time(L0glm.out <- L0glm.fit_new(X, y, 
                           family = poisson(identity), 
                           lambda = lambda_preset,
                           nonzero.centered = TRUE, # gives slightly better solution based on IC, about the same speed or slightly faster
                           ridge.center = NULL, 
                           nonnegative = TRUE,
                           relaxed = FALSE, # TODO: CHECK: relaxed argument has no effect in single channel fit - is this as intended? Should one not be allowed to fit an L0glm and then refit all selected variables unpenalized with lambda=0? 
                           warmupiters = 4)) # should be zeto to 4 for nongaussian, eg poisson
# 0.25s with nonzero.centered = TRUE, 0.23s with nonzero.centered = FALSE
points(sim$x, L0glm.out, pch = 16, col= 'blue')
# graph2pdf(file="..//benchmarks/plots/1_1 single channel deconvolution Poisson.pdf", width=7, height=5)

# some metrics
do.call(cbind, get.metrics(beta = L0glm.out, 
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

system.time(L0glm.out <- L0glm.fit_new(X, y, 
                           family = gaussian(identity), 
                           lambda = lambda_preset, 
                           nonzero.centered = FALSE, # gives slightly better solution based on IC, about the same speed or slightly faster
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

system.time(L0glm.out <- L0glm.fit_new(X, y, 
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

system.time(L0glm.out <- L0glm.fit_new(X, y, 
                           family = gaussian(identity), 
                           lambda = lambda_preset, 
                           nonzero.centered = FALSE, # gives slightly better solution based on IC, about the same speed or slightly faster
                           ridge.center = NULL,
                           nonnegative = TRUE,
                           grouped = TRUE,
                           group.maxit = 20,
                           relaxed = FALSE,
                           warmupiters = 1))
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
