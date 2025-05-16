## BENCHMARK OF GC-MS DECONVOLUTION USING THE L0BestRankGNMF ALGORITHM & COMPARISON WITH MCR-ALS (ANNLS) #### 

## Below we deconvolute overlapping peaks in multichannel GC/MS data using L0BestRankGNMF.
# These data consist of a blurred mutichannel spike train with Poisson noise.
# We use n = 200 timepoints & a true rank k = 50,
# i.e. with 50 compounds being present. The aim is to simultaneously recover
# elution profiles and mass spectra of the true compounds that are present,
# without picking up too many false positives. We then compare the result
# with that obtained using the ALS package, which uses MCR-ALS (alternating
# nonnegative least squares with isotonic regression to impose unimodality
# constraints) and show that with that approach one ends up with a large number
# of false positives, as well as with worse prediction performance.

library(Matrix)
library(L0glm)
devtools::load_all(".")

n <- 200 # timepoints
k_true <- 50 # true rank (nr of compounds)
peak_width <- 2.5 # width / sigma of Gaussian peak shapes 
psf_fun <- function(x, u, w=peak_width) { exp(((x-u)^2)/(-2*(w^2))) } # gaussian point spread / peak shape function
psf <- psf_fun(x=0:20, u=10) # vector with PSF
plot(psf, type="o")
sim.pois.mc <- simulate_spike_train(n = n, k = k_true, family = "poisson", 
                                    seed = 123,
                                    multichannel = TRUE,
                                    mean_beta = 10000,
                                    sparse = TRUE,
                                    sd_logbeta = 1,
                                    psf = psf_fun)
D = sim.pois.mc$y # n x p with n timepoints and p mass channels
H_true = sim.pois.mc$beta_true[rowSums(sim.pois.mc$beta_true)!=0,,drop=F]
image(D)
image(H_true)

# initialise elution profiles (W matrix)
# we initialise W with an n x n banded matrix with 
# time-shifted copies of the point spread / peak shape function (PSF) 
# (this can be estimated by running standards or tuned based on the IC of different L0glm fits 
# with banded matrices with different peak widths)
Wstart = sim.pois.mc$X  

# TODO: might be worth adding another simulate_spike_train matrix to D
# with say 10 compounds and much wider & higher peaks, to see how the algo
# would then be able to resolve those

# this matrix can subsequently be pre-subsetted by regressing the row sums of D
# (=total ion chromatogram) on our banded covariate matrix Wstart using an
# unpenalised nonnegative identity link Poisson fit, fit using osqp and with
# lambda=0, after which we can retain only the columns of Wstart with nonzero
# coefficients and use this as initialisation for W
keep <- (L0glm.fit(X = Wstart, 
                        y = rowSums(D), # n x p data matrix, now regressing on row sums of D (=total ion chromatogram)
                        intercept = FALSE,
                        family = poisson(link=identity),
                        lambda = 0, 
                        nonnegative = TRUE,
                        solver = "osqp" # using unregularised nonnegative identity link Poisson fit using osqp and with lambda=0
           )$coefficients != 0)
sum(keep) # 77
Wstart_rowsums = as.matrix(Wstart[, keep, drop=F])

system.time(nmf <- L0BestRankGNMF(D, 
                                  Wstart = Wstart_rowsums,  # initially selected elution profiles
                                  grouped_Hupdate = TRUE, # here in the H=t(S) with S=mass spectra we use a grouped L0 penalty 
                                  IC = "gic", # information criterion to be optimised
                                  max_iter = 5, tol = 1e-4,
                                  method = 4L, parallelize = FALSE, 
                                  family = poisson(link="identity"))) # nonnegative identity link Poisson here
# 29s unparallelized
dim(nmf$W)
dim(nmf$H)
nmf$iter # run for 5 iterations here
nmf$k # 46 out of 50 components recovered
nmf$dev_trace # objective over subsequent iterations 
plot(nmf$dev_trace, type="o", ylab="Penalized Bregman (half-deviance) loss", xlab="Iteration")
loss_nmf <- neg_loglik(D, mu = as.matrix(nmf$W %*% nmf$H), fam = poisson(identity))
loss_nmf # unpenalised negative log likelihood 
# 11425.88
penloss_nmf <- nmf$dev # final penalised negative log-likelihood objective
penloss_nmf
# 12084.11

FP = max((nmf$k - k_true), 0) 
FP # 0 false positives
FN = max((k_true - nmf$k), 0) 
FN # 4 false negatives
nmf$k / k_true # 0.92, i.e. nearly all compounds recovered

matplot(nmf$W, type="l") # obtained elution profiles
image(H_true) # true spectra
image(nmf$H) # recovered spectra
# image(nmf$H-sim.pois.mc$beta_true[rowSums(sim.pois.mc$beta_true)!=0,, drop=F]) # residual against ground truth for spectra

image(as.matrix(nmf$W %*% nmf$H))    # sparse NMF reconstruction
image(as.matrix(sim.pois.mc$y))      # observed data matrix (with Poisson noise)
image(as.matrix(sim.pois.mc$y_true)) # ground truth (without Poisson noise)
image(as.matrix(nmf$W %*% nmf$H-sim.pois.mc$y_true)) # residual against ground truth
range(as.matrix(nmf$W %*% nmf$H-sim.pois.mc$y_true)) # -144.5284  129.1605
mean(abs(as.matrix(nmf$W %*% nmf$H-sim.pois.mc$y_true))) # MAE against ground truth = 2.811112
image(as.matrix(nmf$W %*% nmf$H-sim.pois.mc$y)) # residual against observed signal
range(as.matrix(nmf$W %*% nmf$H-sim.pois.mc$y)) # -150.4403  175.9993

# TODO: calculate maximum pairwise cosine similarity between each recovered spectrum (rows in nmf$H)
# and the best matching true spectra (rows in H_true) and report boxplot with average cosine similarities of best match 

# cf ALS::matchFactor(S[,1],test0$S[,1]) (dot product match factor)

# NOTE: there are still significant optimisations of our algorithm possible.
# E.g. we could pre-cluster the outcome mass channels, e.g. using XMeans clustering,
# aggregate the signal within each cluster, determine elution profiles and optimal
# rank within that reduced aggregate latent space, and finally refit on the original data
# In addition, it would be possible to use a hierarchical-ALS style (here
# hierarchical-maximum likelihood) algorithm with local component-wise updates, 
# instead of using global updates as in our L0BestRankGNMF algo. That would have
# the advantage that unimodality constraints could more easily and be more stably
# incorporated, e.g. by smoothing the obtained estimates with a log-concave spline
# fit obtained using the cobs package after each component-wise update.
# This would scale better and more stably to large datasets.
# Further extensions for tensor factorisation or simultaneous nonnegative matrix
# factorisation with automatic best rank estimation would also be possible, but 
# this remains an avenue for future research.


## comparison with solution provided by ALS package (which uses plain alternating nonnegative least squares, i.e. ignoring Poisson noise and using a fixed rank decomposition, but including optional unimodality constraints)
library(ALS)
system.time(nmf_als <- als(CList = list(Wstart_rowsums), # starting values for elution profiles
                           S = matrix(1, nrow=ncol(D), ncol=ncol(Wstart_rowsums)), # dummy matrix for S
                           PsiList = list(as.matrix(D)), # data matrix 
                           optS1st = TRUE,
                           uniC = TRUE, 
                           maxiter = 30,
                           normS = 1)) # 10.78s, very buggy though - often gives back NAs with higher maxitr
nmf_als$W = nmf_als$CList[[1]]
nmf_als$H = t(nmf_als$S)
nmf_als$k = nrow(nmf_als$H)
nmf_als$k # 77 compounds recovered, but 27 false positives
FP = (nmf_als$k - k_true) 
FP # 27 false positives
nmf_als$k / k_true # 1.54, i.e. estimated rank 54% too high
loss_nmf_als <- neg_loglik(D, mu = as.matrix(nmf_als$W %*% nmf_als$H), fam = poisson(identity))
loss_nmf_als # unpenalised negative log likelihood 
# 28886.19 vs 11425.88 for L0BestRankGNMF (so much worse)
penloss_nmf_als <- penalised_nll(W = nmf_als$W, H = nmf_als$H, D = D, fam = poisson(link="identity"), 
                                 nmf$lambdaW, nmf$lambdaH)
penloss_nmf_als
# 29988.01 vs 12084.11 for L0BestRankGNMF (so much worse)

matplot(nmf_als$W, type="l") # obtained elution profiles
image(sim.pois.mc$beta_true[rowSums(sim.pois.mc$beta_true)!=0,, drop=F]) # true spectra
image(Matrix(nmf_als$H, sparse=T)) # recovered spectra, with 36 false positives

image(as.matrix(nmf_als$W %*% nmf_als$H))    # sparse NMF reconstruction
image(as.matrix(sim.pois.mc$y))      # observed data matrix (with Poisson noise)
image(as.matrix(sim.pois.mc$y_true)) # ground truth (without Poisson noise)
image(as.matrix(nmf_als$W %*% nmf_als$H-sim.pois.mc$y_true)) # residual against ground truth
range(as.matrix(nmf_als$W %*% nmf_als$H-sim.pois.mc$y_true)) # -198.0029  216.6830
mean(abs(as.matrix(nmf_als$W %*% nmf_als$H-sim.pois.mc$y_true))) # MAE against ground truth = 4.589117
image(as.matrix(nmf_als$W %*% nmf_als$H-sim.pois.mc$y)) # residual against observed signal
range(as.matrix(nmf_als$W %*% nmf_als$H-sim.pois.mc$y)) # -181.44258   92.20318

