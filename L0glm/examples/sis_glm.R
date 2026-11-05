set.seed(1)
n <- 150; p <- 300
X <- matrix(rnorm(n*p),n,p)
beta <- c(1,-1,.8,-.8,rep(0,p-4))
y <- as.numeric(X%*%beta)+rnorm(n)
screen <- sis_glm(X,y,nsis=40,return_stats=FALSE)
screen$indices
threaded <- sis_glm(X,y,nsis=40,return_stats=FALSE,nthreads=4)
identical(threaded$indices,screen$indices)

# Independent channels are screened jointly without signed cancellation.
multi <- sis_glm(X,cbind(y,-y),nsis=40,return_stats=FALSE)
multi$indices

# Conditional screening includes the specified nuisance design.
conditional <- sis_glm(X,y,nuisance=cbind(1,X[,1]),nsis=40)
head(conditional$stats[conditional$ranking,])

# Offset-based outer ISIS; each candidate fit starts afresh.
iterated <- sis_glm(X,y,nsis=40,iterative=TRUE,niter=10,return_stats=FALSE)
iterated$history
coef(iterated$fit)                 # sparse, all original predictor rows
head(predict(iterated$fit,X[1:5,,drop=FALSE]))

# The same full-dimension model through the main fitting interface.
fit <- L0glm.fit(X,y,screen="ISIS",nsis=40,screen.maxit=10,screen.nthreads=4)
fit$screening$history
coef(fit)
# screen="auto" (the default) enables ISIS exactly when p_eligible > 10*n.
# Use screen="none" to fit the full design directly.

if (requireNamespace("mgcv",quietly=TRUE)) {
  counts <- rnbinom(n,mu=exp(.5+.4*X[,1]),size=2)
  nb <- sis_glm(X,counts,family=mgcv::nb(),nsis=20,return_stats=FALSE)
  nb$theta
  classes <- factor(sample(letters[1:3],n,TRUE))
  multinomial <- sis_glm(X,classes,family=mgcv::multinom(2),
                        nsis=20,return_stats=FALSE)
}
if (requireNamespace("mgcv",quietly=TRUE) && requireNamespace("survival",quietly=TRUE)) {
  time <- rexp(n,exp(.5*X[,1]))
  outcome <- survival::Surv(pmin(time,2),time<=2)
  cox <- sis_glm(Matrix::Matrix(X,sparse=TRUE),outcome,
                family=mgcv::cox.ph(),nsis=20,return_stats=FALSE)
}

\dontrun{
# One million predictors, 500 observations and 60 nonzero effects.
# The dense design occupies about 4 GB; data generation is outside the timings.
set.seed(123)
n <- 500L; p <- 1000000L; nsis <- 100L
X <- matrix(rnorm(as.double(n) * p), n, p)
colnames(X) <- paste0("x", seq_len(p))
active <- sort(sample.int(p, 60L))
beta <- sample(c(-1, 1), 60L, TRUE) * runif(60L, .5, 1.5)
y <- as.numeric(X[, active, drop = FALSE] %*% beta) + rnorm(n)

t_fast <- system.time(fast <- sis_glm(X, y, nsis = nsis))
t_one <- system.time(one <- sis_glm(X, y, nsis = nsis, nthreads = 1L))
t_SIS <- system.time(reference <- SIS::SIS(X, y, family = "gaussian",
  nsis = nsis, iter = FALSE, penalty = "lasso", tune = "bic",
  standardize = FALSE, parallel = FALSE, seed = 123L))

data.frame(method = c("sis_glm, default threads", "sis_glm, one thread", "SIS::SIS"),
  seconds = c(t_fast[["elapsed"]], t_one[["elapsed"]], t_SIS[["elapsed"]]))

# Compare the initial screening rank, including all one million predictors.
identical(fast$indices, as.integer(reference$sis.ix0))       # TRUE
identical(fast$ranking, as.integer(reference$ix_list[[1L]]))  # TRUE
identical(one$ranking, fast$ranking)                        # TRUE
cor(fast$stats$rank, order(reference$ix_list[[1L]]))          # 1

# Median elapsed times over three repetitions on an i9-11950H:
# sis_glm, default eight physical cores: 0.23 s
# sis_glm, one thread:                  0.82 s
# SIS::SIS:                            2.69 s
# R 4.5.0, SIS 1.5; identical full ranks in every repetition.
# SIS::SIS also fits a BIC-tuned LASSO on the 100 retained predictors.
# Gaussian correlation ranks are unchanged by standardize=FALSE, which avoids
# making a redundant standardized copy of this large design.
}
