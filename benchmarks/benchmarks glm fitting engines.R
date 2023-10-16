# see https://stackoverflow.com/questions/75341154/any-way-to-run-glm-or-similar-function-in-r-when-i-only-want-to-compute-mle-e

library(fastglm)
library(speedglm)
library(parglm)
library(oem)
library(glmnet)
library(microbenchmark)
library(ggplot2)

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
                             method="FAST")),
  bigglm = bigGlm(X, dd$y,
                  alpha = 0, lambda = 0, intercept = FALSE, 
                  standardize = FALSE, 
                  standardize.response = FALSE,
                  family = b, path = FALSE),
  oem = oem(X, dd$y, family="binomial", penalty="ols", standardize=F, intercept=F),
  times  = 20)

autoplot(m1)
