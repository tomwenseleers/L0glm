sis_thread_cases <- function() {
  set.seed(2710);n <- 120L;p <- 19L
  X <- matrix(rnorm(n*(p-3L)),n,p-3L);X[abs(X)<.7] <- 0
  X <- cbind(X,X[,1L],0,1)
  N <- cbind(1,runif(n,-.05,.05));w <- sample(0:3,n,TRUE);off <- runif(n,-.01,.01)
  families <- c(lapply(c("identity","log","inverse"),gaussian),
    lapply(c("logit","probit","cloglog","cauchit","identity","log"),binomial),
    lapply(c("log","identity","sqrt"),poisson),
    lapply(c("log","identity","inverse"),Gamma),
    lapply(c("log","identity","inverse","1/mu^2"),inverse.gaussian),
    lapply(c("logit","probit","cloglog","cauchit","identity","log"),quasibinomial),
    lapply(c("log","identity","sqrt"),quasipoisson),
    list(quasi("identity","constant"),quasi("log","mu"),quasi("log","mu^2"),
      quasi("log","mu^3"),quasi("logit","mu(1-mu)")),
    lapply(c("log","identity","sqrt"),function(link) MASS::negative.binomial(2,link)))
  if (requireNamespace("mgcv",quietly=TRUE)) for (theta in c(2,0))
    families <- c(families,lapply(c("log","identity","sqrt"),function(link)
      do.call(mgcv::nb,list(theta=theta,link=link))))
  cases <- list()
  for (i in seq_along(families)) {
    f <- families[[i]]
    binary <- f$family%in%c("binomial","quasibinomial") || f$family=="quasi" && f$link=="logit"
    eta <- if (binary) switch(f$link,identity=.4+off,log=-1+off,off) else
      if (f$link%in%c("identity","inverse","sqrt","1/mu^2")) 2+off else .7+off
    mu <- f$linkinv(eta)
    y <- if (binary) rbinom(n,1,mu) else if (f$family%in%c("poisson","quasipoisson")) rpois(n,mu) else
      if (grepl("negative binomial",f$family,ignore.case=TRUE)) rnbinom(n,mu=mu,size=2) else
        if (f$family=="gaussian") mu+rnorm(n,0,.03) else
          if (f$family=="inverse.gaussian") statmod::rinvgauss(n,mean=mu,dispersion=.02) else rgamma(n,40,scale=mu/40)
    label <- paste(i,f$family,f$link,sep=":")
    cases[[label]] <- list(X=X,y=y,family=f,weights=w,offset=off,nuisance=N)
  }
  for (q in 0:1) cases[[paste0("gaussian-null-",q)]] <- list(X=X,y=rnorm(n),family=gaussian(),
    weights=w,offset=off,nuisance=N[,seq_len(q),drop=FALSE])
  Y <- cbind(rnorm(n),rnorm(n));Y[,2] <- -Y[,1]
  cases$independent_gaussian <- list(X=X,y=Y,family=gaussian(),weights=cbind(w,rev(w)),
    offset=cbind(off,-off),nuisance=N)
  cases$independent_poisson <- list(X=X,y=matrix(rpois(n*2,2),n,2),family=poisson(identity),
    weights=cbind(w,rev(w)),offset=cbind(off,-off),nuisance=N)
  if (requireNamespace("mgcv",quietly=TRUE)) {
    cases$independent_nb <- list(X=X,y=matrix(rnbinom(n*2,mu=2,size=2),n,2),family=mgcv::nb(),
      weights=cbind(w,rev(w)),offset=cbind(off,-off),nuisance=N)
    for (q in c(0L,2L)) cases[[paste0("multinomial-counts-",q)]] <- list(X=X,
      y=matrix(rpois(n*3,2),n,3),family=mgcv::multinom(2),weights=w,
      offset=matrix(rnorm(n*3,0,.02),n,3),nuisance=N[,seq_len(q),drop=FALSE])
    cases$multinomial_factor <- list(X=X,y=factor(sample(letters[1:3],n,TRUE)),
      family=mgcv::multinom(2),weights=w,offset=matrix(0,n,3),nuisance=N)
    if (requireNamespace("survival",quietly=TRUE)) for (q in c(0L,2L))
      cases[[paste0("cox-ties-",q)]] <- list(X=X,y=survival::Surv(sample(1:30,n,TRUE),rbinom(n,1,.7)),
        family=mgcv::cox.ph(),weights=w,offset=off,nuisance=N[,seq_len(q),drop=FALSE])
  }
  cases
}
