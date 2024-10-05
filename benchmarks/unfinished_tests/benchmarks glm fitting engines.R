# see https://stackoverflow.com/questions/75341154/any-way-to-run-glm-or-similar-function-in-r-when-i-only-want-to-compute-mle-e

library(fastglm)
library(speedglm)
library(parglm)
library(oem)
library(glmnet)
library(microbenchmark)
library(ggplot2)
devtools::load_all()

set.seed(101)
n <- 1e6
dd <- data.frame(x = rnorm(n))
dd$y <- rbinom(n, size = 1, prob = plogis(2-3*dd$x))
X <- model.matrix(~x, data = dd)

b <- binomial()
m1 <- microbenchmark(
  glm = glm(y~x, dd, family = b),
  glm.fit = glm.fit(X, dd$y, family = b),
  fastglm = fastglm(X, dd$y, family = b),
  fpure1 = fastglmPure(X, dd$y, family = b, method = 1),
  fpure2 = fastglmPure(X, dd$y, family = b, method = 2),
  fpure3 = fastglmPure(X, dd$y, family = b, method = 3),
  speedglm.wfit = speedglm.wfit(dd$y, X, intercept = F, family = b, chunksize = 1000, sparse = F),
  parglm.fit = parglm.fit(X, dd$y, intercept = F, family = b, 
                      control = parglm.control(nthreads=8L, 
                             block_size=1000, 
                             method="FAST")), # this is using the arma::chol Cholesky LLT solver
  bigglm = bigGlm(X, dd$y,
                  alpha = 0, lambda = 0, intercept = FALSE, 
                  standardize = FALSE, 
                  standardize.response = FALSE,
                  family = b, path = FALSE),
  oem = oem(X, dd$y, family="binomial", penalty="ols", standardize=F, intercept=F),
  times  = 20)

autoplot(m1)


# Gamma GLM, tall problem here  ####

library(biglm)
library(speedglm)
library(glmnet)
library(parglm)
n<-10000000  # n<-1000000 # 
k<-70
y <- rgamma(n,1.5,1)
x <-round( matrix(rnorm(n*k),n,k),digits=3)
colnames(x) <-paste("s",1:k,sep = "")
da<- data.frame(y,x)
fo <- as.formula(paste("y~",paste(paste("s",1:k,sep=""),collapse="+")))   
X <- model.matrix(fo, data=da)

system.time(m1 <- glm(fo,data=da,family=Gamma(log))) 
# 14s with n=1000000, out of memory with n=10000000
system.time(m2 <- bigglm(fo,data=da,family=Gamma(log), maxit=100)) 
# 28s with n=1000000, 273s with n=10000000
system.time(m3 <- fastglmPure(X, y, family=Gamma(log), method = 2)) 
# 3.94s with n=1000000, 41s with n=10000000
system.time(m4 <- speedglm.wfit(y=y, X=X, intercept=F, data=da, family=Gamma(log))) 
# 3.8s with n=1000000, out of memory with n=10000000 - shglm should still work??
system.time(m5 <- parglm.fit(X, y, intercept = F, 
                             control = parglm.control(nthreads=8L, 
                                                      # maxit=1,
                                                      block_size=1000, 
                                                      method="FAST"), 
                             # this is using the arma::chol Cholesky LLT solver
                             data=da,family=Gamma(log))$coefficients) 
# 1s for n<-1000000, 9s for n<-10000000
system.time(m5b <- parglm:::parallelglm(X=X, Y=y, 
                                        family=paste0(Gamma(log)$family, "_", Gamma(log)$link), 
                                        start=rep(0,ncol(X)),
                                        weights=rep(1, nrow(X)), 
                                        offsets=rep(0,nrow(X)),
                                        tol=1E-8, nthreads=8L, 
                                        it_max=100L, trace=FALSE,
                                        method="FAST", block_size=1000, 
                                        use_start=FALSE)$coefficients) 
# 0.45s for n<-1000000, 4.8s for n<-10000000

system.time(m6 <- L0glm.fit(X, y, 
                            lambda=0, # "aic",
                            family=Gamma(log),
                            normalize=F,
                            solver="parglm",
                            method=3L)$coefficients)  
# 11s with parglm, 15s with eigen LSCG for n<-1000000
# 104s with parglm for n<-10000000
# NOTE: for dense n>>p problems parglm:::parallelglm solution could be used
# as starting values in L0glm
# TODO: fix L0glm.fit so that it doesn't run out of memory with n<-10000000 & lambda="aic"
# check if memory efficiency could be improved & where exactly things are 
# going wrong
# 
# prob still a lot of scope for further code optimisation, given that using parglm
# within my algo is so much slower than the full Rcpp coded parglm:::parallelglm
system.time(m7 <- bigGlm(x=X, y=y, 
                            family=Gamma(log),
                            alpha=0, 
                            standardize = F,
                            lambda=0, 
                            intercept = F, 
                            thresh = 1E-10)) # 18s, or Error: cannot allocate vector of size 2.7 Gb
