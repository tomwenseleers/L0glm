set.seed(1)
n <- 150; p <- 300
X <- matrix(rnorm(n*p),n,p)
beta <- c(1,-1,.8,-.8,rep(0,p-4))
y <- as.numeric(X%*%beta)+rnorm(n)
screen <- sis_glm(X,y,nsis=40,return_stats=FALSE)
screen$indices

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
fit <- L0glm.fit(X,y,screen="ISIS",nsis=40,screen.maxit=10)
fit$screening$history
coef(fit)
# screen="auto" (the default) enables ISIS only for ultra-wide eligible designs.
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
