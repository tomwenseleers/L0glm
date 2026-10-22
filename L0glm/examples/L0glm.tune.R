# Small versions of the L0glm.fit deconvolution and Gaussian ISIS examples.
# Five lambda values and three folds keep these examples short; increase
# nlambda and nfolds for a finer search. True coefficients are only used
# for plots and recovery metrics, never for fitting or choosing lambda.
old_par <- par(mgp=c(2, .7, 0))

# Single-channel Poisson deconvolution.
# Zoom around the GIC preset; use prop_range=1 for the full lambda range.
single <- simul(method=2, n=500, p=500, s=.1, seed=2)
single_gic <- L0glm.tune(single$X, single$y, family=poisson(identity),
  tune="gic", nlambda=5, prop_range=.25, coefs.true=single$beta_true,
  nonnegative=TRUE, relaxed=TRUE, parallelize=FALSE, plot=TRUE)

# Venetian-blind CV: train on two out of every three time points.
folds <- rep(1:3, length.out=nrow(single$X))
single_cv <- L0glm.tune(single$X, single$y, family=poisson(identity),
  tune="cv", nlambda=5, prop_range=.25, nfolds=3, foldid=folds,
  coefs.true=single$beta_true,
  nonnegative=TRUE, relaxed=TRUE, parallelize=FALSE, plot=TRUE)
c(GIC=single_gic$lambda_best, CV=single_cv$lambda_best,
  CV_1se=single_cv$lambda_1se, GIC_preset=single_cv$lambda_preset)
single_cv$metrics[c("lambda", "cv_deviance", "cv_se")]
L0glm.metrics(single_gic$fit_best, beta_true=single$beta_true)[
  c("gic", "TP", "FP", "FN", "TN", "F1_score", "MCC")]
L0glm.metrics(single_cv$fit_best, beta_true=single$beta_true)[
  c("gic", "TP", "FP", "FN", "TN", "F1_score", "MCC")]
# fit_1se is a full-data model, ready for prediction or coefficient extraction.
fitted(single_cv$fit_1se)[1:6, , drop=FALSE]

# Grouped multichannel Poisson deconvolution, keeping every channel.
multi <- simul(method=4, n=500, p=500, s=.05, seed=2)
multi_gic <- L0glm.tune(multi$X, multi$y, family=poisson(identity),
  tune="gic", nlambda=5, prop_range=.25, coefs.true=multi$beta_true,
  nonnegative=TRUE, grouped=TRUE, relaxed=TRUE, parallelize=FALSE, plot=TRUE)

# All channels of each time point share one fold assignment.
folds <- rep(1:3, length.out=nrow(multi$X))
multi_cv <- L0glm.tune(multi$X, multi$y, family=poisson(identity),
  tune="cv", nlambda=5, prop_range=.25, nfolds=3, foldid=folds,
  coefs.true=multi$beta_true,
  nonnegative=TRUE, grouped=TRUE, relaxed=TRUE, parallelize=FALSE, plot=TRUE)
c(GIC=multi_gic$lambda_best, CV=multi_cv$lambda_best,
  CV_1se=multi_cv$lambda_1se, GIC_preset=multi_cv$lambda_preset)
L0glm.metrics(multi_gic$fit_best, beta_true=multi$beta_true)[
  c("gic", "TP", "FP", "FN", "TN", "F1_score", "MCC",
    "within_row_F1_score", "within_row_MCC")]
L0glm.metrics(multi_cv$fit_best, beta_true=multi$beta_true)[
  c("gic", "TP", "FP", "FN", "TN", "F1_score", "MCC",
    "within_row_F1_score", "within_row_MCC")]
dim(coef(multi_cv$fit_1se))            # all original predictor rows and channels

# High-dimensional Gaussian regression with 60 positive/negative effects.
# X alone needs about 4 GB; set p <- 20000L for a smaller example.
\dontrun{
n <- 500L
p <- 1000000L
set.seed(123)
X_hd <- matrix(rnorm(n*p), n, p)
support <- sample.int(p, 60L)
beta_hd <- numeric(p)
beta_hd[support] <- runif(60L, -.8, .8)
mean_hd <- as.numeric(X_hd[, support, drop=FALSE] %*% beta_hd[support])
set.seed(1)
y_hd <- mean_hd + rnorm(n)

# screen="auto" applies p_eligible > 10*n independently at each lambda,
# with n equal to the training row count, including within each CV fold.
hd_gic <- L0glm.tune(X_hd, y_hd, family=gaussian(), tune="gic",
  nlambda=5, coefs.true=beta_hd, screen.nthreads=8, nthreads=1, plot=FALSE)
# Screening and preprocessing are recomputed within each training fold.
set.seed(42)
hd_cv <- L0glm.tune(X_hd, y_hd, family=gaussian(), tune="cv",
  nlambda=5, nfolds=3, coefs.true=beta_hd,
  screen.nthreads=8, nthreads=1, plot=FALSE)
c(GIC=hd_gic$lambda_best, CV=hd_cv$lambda_best,
  CV_1se=hd_cv$lambda_1se, GIC_preset=hd_cv$lambda_preset)
hd_cv$metrics[c("lambda", "cv_deviance", "cv_se")]
L0glm.metrics(hd_gic$fit_best, beta_true=beta_hd)[
  c("gic", "TP", "FP", "FN", "TN", "F1_score", "MCC")]
L0glm.metrics(hd_cv$fit_best, beta_true=beta_hd)[
  c("gic", "TP", "FP", "FN", "TN", "F1_score", "MCC")]
hd_cv$fit_best$screening$history
dim(coef(hd_cv$fit_1se))              # screened-out predictors remain zero rows
}

# CV plots show mean held-out deviance +/- one standard error. Vertical
# lines mark lambda_best, lambda_1se and the full-data GIC preset.
par(old_par)
