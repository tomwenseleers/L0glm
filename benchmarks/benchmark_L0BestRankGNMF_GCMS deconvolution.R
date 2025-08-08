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
# devtools::load_all(".")
library(L0glm)

n <- 500 # timepoints
k_true <- 100 # true rank (nr of compounds)
peak_width_range <- c(1, 5) # peak width range
sim <- simulate_GCMS(n = n, 
                     k = k_true,
                     peak_width_range = peak_width_range,
                     mean_beta   = 1e4,
                     sd_logbeta  = 1)
D = sim$D # n x m with n timepoints and m mass channels
H_true = sim$H_true
W_true = sim$W_true
matplot(W_true, type="l") # true elution profiles
plotmat(D, xlab="Elution time (min)", ylab="Mass (m/z)") # observed chromatogram
plotmat(H_true, xlab="Compound", ylab="Mass (m/z)") # true mass spectra

# initialise elution profiles (W matrix)
# we estimate a sensible starting estimate of the typical PSF using
# L0glm_blind_deco 
# (it would also be possible to estimate it from an alkane ladder standard)
L = 51 # length of PSF
g0 = convkernel("gaussian", L, list(sigma = 6))  # initial guess for PSF
system.time(fitpsf <- L0glm_blind_deco(y = rowSums(D), # we use TIC 
                          family = poisson(identity), 
                          baseline = NULL, # or "intercept" or "bspline"
                          L = L,                 # PSF kernel length 
                          kernel_est = g0, 
                          kernel_smoothfun = function(y) log_conc_smooth(y, nknots=round(L/5)),
                          perc = "mean",         # use centroid centering
                          IC = "gic",
                          plot = TRUE,
                          tol = 1e-3))
# 3s
Wstart = convmatrix(fitpsf$g_hat, nrow(D)) # n x n banded convolution matrix

# this matrix can subsequently be pre-subsetted by doing a best subset
# regression of the row sums of D
# (=total ion chromatogram) on our banded covariate matrix Wstart using L0glm
keep <- (L0glm.fit(X = Wstart, 
                        y = rowSums(D), # n x p data matrix, now regressing on row sums of D (=total ion chromatogram)
                        intercept = FALSE,
                        family = poisson(link=identity),
                        lambda = "gic",  
                        nonnegative = TRUE,
                        solver = "eigen",
                        method = 7L # using unregularised nonnegative identity link Poisson fit using osqp solver and with lambda=0
           )$coefficients != 0)
sum(keep) # 101
Wstart = as.matrix(Wstart[, which(keep), drop=F])
matplot(Wstart, type="l")

# for faster NMF performance we first cluster & pool correlated mass traces ####
Dagg <- glm_cluster(D,
            family   = poisson(identity),
            method   = "poissondist",
            rowbins  = 5)$Y_agg # aggregate together every bin of rowbins=5 subsequent timepoints for extra performance
matplot(Dagg^0.5, type="l", lty=1) # looks good

# now fit L0BestRankGNMF model on aggregated data Dagg
IC = "gic" # information criterion to be optimised

system.time(nmfDagg <- L0BestRankGNMF(Dagg, 
                                  Wstart = Wstart,  # initially selected elution profiles
                                  grouped_Hupdate = TRUE, # here in the H=t(S) with S=mass spectra we use a grouped L0 penalty 
                                  IC = IC, # information criterion to be optimised
                                  max_iter = 10, tol = 1e-4,
                                  method = 7L, parallelize = FALSE, 
                                  family = poisson(link="identity"))) # nonnegative identity link Poisson here
# 27s unparallelized
dim(nmfDagg$W) # 500 x 100
dim(nmfDagg$H) # 100 x 13
nmfDagg$iter # run for 7 iterations here
nmfDagg$k # 100 components retained
nmfDagg$dev_trace # objective over subsequent iterations 
plot(nmfDagg$dev_trace, type="o", log="y", 
     ylab="Penalized deviance loss", xlab="Iteration")
loss_nmf <- neg_loglik(Dagg, mu = as.matrix(nmfDagg$W %*% nmfDagg$H), fam = poisson(identity))
loss_nmf # unpenalised -2*negative log likelihood 
# 23647.15
penloss_nmf <- nmfDagg$dev # final penalised negative log-likelihood objective
penloss_nmf
# 9047.767

nmfDagg$k / k_true # 1 - approx correct?

matplot(W_true, type="l", lty=1, 
        xlab = "Time index",
        ylab = "Elution profile",
        main = "True (top) & estimated (bottom) elution profiles (W)",
        ylim=c(-max(W_true), max(W_true))) # true elution profiles
matplot(-nmfDagg$W, type="l", lty=1, ylim=c(-max(W_true), max(W_true)), add=TRUE) # obtained elution profiles
# for compounds with peaks that are wider than normal peaks would still be split
# these could still be regrouped based on the cosine similarity (in aggregated space) 
# in a post-clustering step

# finally refit selected features on original nonaggregated data to get mass spectra
# nmfD <- XXX FINISH XXX


plotmat(nmf$H) # recovered spectra
plotmat(H_true) # true spectra

plotmat(nmf$W %*% nmf$H)    # sparse NMF reconstruction in aggregated space

plotmat(sim$D)              # observed data matrix (with Poisson noise)
plotmat(sim$D_true)         # ground truth (without Poisson noise)

# plotmat(nmf$W %*% nmf$H-sim$D_true) # residual against ground truth
# range(as.matrix(nmf$W %*% nmf$H-sim$D_true)) # -239.2441  181.8252
# mean(abs(as.matrix(nmf$W %*% nmf$H-sim$D_true))) # MAE against ground truth = 2.743787
# plotmat(nmf$W %*% nmf$H-sim$D) # residual against observed signal
# range(as.matrix(nmf$W %*% nmf$H-sim$D)) # -299.5757  211.6712

# TODO: calculate maximum pairwise cosine similarity between each recovered spectrum (rows in nmf$H)
# and the best matching true spectra (rows in H_true) and report boxplot with average cosine similarities of best match 

# cf ALS::matchFactor(S[,1],test0$S[,1]) (dot product match factor)

# NOTE: there are still significant optimisations of our algorithm possible.
# E.g. it would be possible to use a hierarchical-ALS style (here
# hierarchical-maximum likelihood) algorithm with local component-wise updates, 
# instead of using global updates as in our L0BestRankGNMF algo. That would have
# the advantage that unimodality constraints could more easily and be more stably
# incorporated, e.g. by smoothing the obtained estimates with a log-concave spline
# fit (see ?log_conc_smooth).
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
                           normS = 1)) # 10.29s, very buggy though - often gives back NAs with higher maxitr & often diverging
nmf_als$W = nmf_als$CList[[1]]
nmf_als$H = t(nmf_als$S)
nmf_als$k = nrow(nmf_als$H)
nmf_als$k # 77 compounds recovered, but 27 false positives
FP = (nmf_als$k - k_true) 
FP # 27 false positives
nmf_als$k / k_true # 1.54, i.e. estimated rank 54% too high
loss_nmf_als <- neg_loglik(D, mu = as.matrix(nmf_als$W %*% nmf_als$H), fam = poisson(identity))
loss_nmf_als # unpenalised negative log likelihood 
# 57772.38 vs 24750.79 for L0BestRankGNMF (so much worse)
penloss_nmf_als <- penalised_nll(W = nmf_als$W, H = nmf_als$H, D = D, fam = poisson(link="identity"), 
                                 nmf$lambdaW, nmf$lambdaH, grouped_Hupdate=TRUE)
penloss_nmf_als
# 58572.09 vs 22251.9 for L0BestRankGNMF (so much worse)

matplot(nmf_als$W, type="l") # obtained elution profiles
image(sim$beta_true[rowSums(sim$beta_true)!=0,, drop=F]) # true spectra
image(Matrix(nmf_als$H, sparse=T)) # recovered spectra, with 36 false positives

image(as.matrix(nmf_als$W %*% nmf_als$H))    # sparse NMF reconstruction
image(as.matrix(sim$y))      # observed data matrix (with Poisson noise)
image(as.matrix(sim$y_true)) # ground truth (without Poisson noise)
image(as.matrix(nmf_als$W %*% nmf_als$H-sim$y_true)) # residual against ground truth
range(as.matrix(nmf_als$W %*% nmf_als$H-sim$y_true)) # -198.0029  216.6830
mean(abs(as.matrix(nmf_als$W %*% nmf_als$H-sim$y_true))) # MAE against ground truth = 4.589117
image(as.matrix(nmf_als$W %*% nmf_als$H-sim$y)) # residual against observed signal
range(as.matrix(nmf_als$W %*% nmf_als$H-sim$y)) # -181.44258   92.20318

# TODO: make proper graphs with nice colours etc
