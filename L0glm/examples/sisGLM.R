set.seed(1)
n <- 150; p <- 300
X <- matrix(rnorm(n*p),n,p)
beta <- c(1,-1,.8,-.8,rep(0,p-4))
y <- as.numeric(X%*%beta)+rnorm(n)
screen <- sisGLM(X,y,nsis=40,return_stats=FALSE)
screen$indices
threaded <- sisGLM(X,y,nsis=40,return_stats=FALSE,nthreads=4)
identical(threaded$indices,screen$indices)

# Independent channels are screened jointly without signed cancellation.
multi <- sisGLM(X,cbind(y,-y),nsis=40,return_stats=FALSE)
multi$indices

# Conditional screening includes the specified nuisance design.
conditional <- sisGLM(X,y,nuisance=cbind(1,X[,1]),nsis=40)
head(conditional$stats[conditional$ranking,])

# Offset-based outer ISIS; each candidate fit starts afresh.
iterated <- sisGLM(X,y,nsis=40,iterative=TRUE,niter=10,return_stats=FALSE)
iterated$history
coef(iterated$fit)                 # sparse, all original predictor rows
head(predict(iterated$fit,X[1:5,,drop=FALSE]))

# The same full-dimension model through the main fitting interface.
fit <- l0glmFit(X,y,screen="ISIS",nsis=40,screen.maxit=10,screen.nthreads=4)
fit$screening$history
coef(fit)
# screen="auto" (the default) enables ISIS exactly when p_eligible > 10*n.
# Use screen="none" to fit the full design directly.

if (requireNamespace("mgcv",quietly=TRUE)) {
  counts <- rnbinom(n,mu=exp(.5+.4*X[,1]),size=2)
  nb <- sisGLM(X,counts,family=mgcv::nb(),nsis=20,return_stats=FALSE)
  nb$theta
  classes <- factor(sample(letters[1:3],n,TRUE))
  multinomial <- sisGLM(X,classes,family=mgcv::multinom(2),
                        nsis=20,return_stats=FALSE)
}
if (requireNamespace("mgcv",quietly=TRUE) && requireNamespace("survival",quietly=TRUE)) {
  time <- rexp(n,exp(.5*X[,1]))
  outcome <- survival::Surv(pmin(time,2),time<=2)
  cox <- sisGLM(Matrix::Matrix(X,sparse=TRUE),outcome,
                family=mgcv::cox.ph(),nsis=20,return_stats=FALSE)
}

\dontrun{
# One million predictors, 1000 observations and 60 nonzero effects.
# The dense design occupies about 8 GB; data generation is outside the timings.
stopifnot(requireNamespace("SIS",quietly=TRUE))
set.seed(123)
n <- 1000L; p <- 1000000L; nsis <- 100L
X <- matrix(rnorm(as.double(n) * p), n, p)
active <- sample.int(p, 60L)
beta <- numeric(p)
beta[active] <- runif(60L, -.8, .8)
set.seed(1)
y <- as.numeric(X[, active, drop = FALSE] %*% beta[active]) + rnorm(n)

# The default fit uses the GIC preset and automatic ISIS screening.
t_default <- system.time(fit <- l0glmFit(X,y,family=gaussian()))
l0glmMetrics(fit,beta_true=beta)[c("gic","TP","FP","FN","TN","F1_score","MCC")]

t_fast <- system.time(fast <- sisGLM(X, y, nsis = nsis))
t_one <- system.time(one <- sisGLM(X, y, nsis = nsis, nthreads = 1L))
# Screening only: this helper runs the initial screen in SIS 1.5 without its
# subsequent LASSO fit. Its internal interface may vary between SIS versions.
sis_screen <- getFromNamespace("obtain.ix0","SIS")
t_screen <- system.time(reference <- sis_screen(X,y,s1=NULL,s2=NULL,
  family="gaussian",nsis=nsis,iter=FALSE,varISIS="vanilla",perm=FALSE,
  q=1,greedy=FALSE,greedy.size=1L,iterind=0L))

data.frame(method=c("sisGLM, default threads","sisGLM, one thread","SIS screen"),
  seconds=c(t_fast[["elapsed"]],t_one[["elapsed"]],t_screen[["elapsed"]]))

# Compare the initial screening rank, including all one million predictors.
identical(fast$indices,as.integer(reference$ix0))       # TRUE
identical(fast$ranking,as.integer(reference$ix0all))    # TRUE
identical(one$ranking,fast$ranking)                   # TRUE

# Complete models: both timings include screening and variable selection.
t_SIS <- system.time(sis_fit <- SIS::SIS(X,y,family="gaussian",nsis=nsis,
  iter=FALSE,penalty="lasso",tune="bic",standardize=FALSE,parallel=FALSE,seed=123L))
t_L0glm <- system.time(l0_fit <- l0glmFit(X,y,
  screen="SIS",nsis=nsis,lambda="gic",relaxed=TRUE))
data.frame(method=c("SIS + BIC LASSO","SIS + L0glm GIC preset"),
  seconds=c(t_SIS[["elapsed"]],t_L0glm[["elapsed"]]))
l0glmMetrics(l0_fit,beta_true=beta)[c("gic","TP","FP","F1_score","MCC")]

# Shared screening also makes end-to-end comparisons convenient for other
# methods. Change screen to "ISIS" for method-specific iterative screening.
comparison <- sparseGLM(X,y,method=c("L0glm","lasso","relaxed_lasso"),
  screen="SIS",nsis=nsis,tune="bic",nlambda=10,beta_true=beta,intercept=FALSE,
  lambda=list(L0glm=presetLambda("gic",p,n)),options=list(L0glm=list(relaxed=TRUE)))
l0glmMetrics(comparison)[c("method","screen","seconds","bic","TP","FP","F1_score","MCC")]
}
