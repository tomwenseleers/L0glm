\donttest{
## Synthetic Gaussian example (dense X): n = 500, p = 1e6
n <- 500; p <- 1e6; k <- 60
set.seed(1)
X <- matrix(rnorm(n * p), n, p)
beta <- numeric(p); beta[sample(p, k)] <- runif(k, -0.8, 0.8)
y <- as.numeric(X %*% beta + rnorm(n))

## Screen with classical SIS size nsis ~ n/log n
nsis <- round(n / log(n))

## Time our sis_glm
tm_sis_glm <- system.time(
   SIS_glm <- sis_glm(X, y,
                      family = gaussian(),
                      weights = NULL, offset = NULL,
                      nsis = nsis, return_stats = FALSE)
)
print(tm_sis_glm)

## Time vanilla SIS::SIS
if (requireNamespace("SIS", quietly = TRUE)) {
  tm_SIS <- system.time(
     SIS_pkg <- SIS::SIS(X, y, family = "gaussian", nsis = nsis, iter = FALSE)
  )
print(tm_SIS)
## On a typical workstation we observe >10x speedup of sis_glm over vanilla SIS::SIS,
## with near-identical selected sets (results vary by hardware and BLAS).
length(intersect(SIS_glm$indices, SIS_pkg$ix))
}
}
