# How to reproduce a ridge regression in glmnet
# see https://stats.stackexchange.com/questions/129179/why-is-glmnet-ridge-regression-giving-me-a-different-answer-than-manual-calculat

devtools::load_all(".")

set.seed(123)
n   <-  10000 # 2000
p   <-  10000 # 2000
xvals <- rnorm(n*p,0,1)
sparsity <- 0.9 # proportion of zeros
xvals[which(rbinom(n=n*p, size=1, prob=sparsity)==1)] <- 0
X   <- matrix(xvals,n,p)
library(Matrix)
X_sparse <- Matrix(X, sparse=T)
beta <- rnorm(p,0,1)
y    <- X%*%beta+rnorm(n,0,0.5)
sd_y <- sqrt(var(y)*(n-1)/n)[1,1]
lambda <- 1
prior.mean <- rep(0, p)
Xaug <- rbind(X, diag(x=rep(sqrt(lambda),p)))
Xaug_sp <- Matrix(Xaug, sparse=T)
yaug <- c(y, rep(prior.mean*sqrt(lambda)))
beta_init <- as.vector(crossprod(Xaug, yaug) / colSums(Xaug ^ 2))
tol <- 1E-6

# manual ridge regression ####

# Cholesky decomposition based solve, using normal equations
chol.solve = function (X, y, tol=tol) { # .Machine$double.eps^(2/3)) 
  ch <- base::chol(crossprod(X)) # solve using Cholesky decomposition - only fast with Intel MKL installed though
  backsolve(ch, forwardsolve(ch, crossprod(X, y), upper = TRUE, trans = TRUE))
}

system.time(beta_man1_dense <- chol.solve(rbind(X, 
                                                diag(x=rep(sqrt(lambda),ncol(X)))), 
                       c(y, prior.mean*sqrt(lambda)))) 
# for n=p=10000: 12.6s with Intel MKL
system.time(beta_man1_sparse <- chol.solve(rbind(X_sparse, 
                                                 Diagonal(x=rep(sqrt(lambda),ncol(X)))),
                       c(y, prior.mean*sqrt(lambda)))) 
# for n=p=10000: 22.6s with Intel MKL

system.time(beta_man1_dense <- solve(crossprod(X)+lambda*diag(x=1, nrow=p, ncol=p), 
                                     crossprod(X,y))[,1]) 
# for n=p=10000: 7.9s with Intel MKL
system.time(beta_man1_sparse <- solve(crossprod(X_sparse)+lambda*Diagonal(n=p), 
                                      crossprod(X_sparse,y))[,1]) 
# for n=p=10000: 30s with Intel MKL

# using Cholesky Eigen solver:
system.time(beta_eigen_chol_dense <- lmridge.eigen(Xw=X,
                                                   yw=y,
                                                   lambdas=rep(lambda,p),
                                                   start=rep(0,p),
                                                   control.solver = 
                                                     modifyList(do.call("control.solver.gen", 
                                                             list(solver = "eigen", method = 7L)),
                                                             list(rowaug = FALSE)),
                                                   prior.mean=rep(0,p)
)) 
# for n=p=10000: 70s with rowaug=T, 65s with rowaug=F
system.time(beta_eigen_chol_sparse <- lmridge.eigen(Xw=X_sparse,
              yw=y,
              lambdas=rep(lambda,p),
              start=rep(0,p),
              control.solver = modifyList(do.call("control.solver.gen", 
                                                  list(solver = "eigen", method = 7L)),
                                          list(rowaug = FALSE)),
              prior.mean=rep(0,p)
              )) 
# for n=p=10000: 143s with rowaug=T, 161s with rowaug=F


# using CG Eigen solver:
system.time(beta_eigen_cg_dense <- lmridge.eigen(Xw=X,
                                                   yw=y,
                                                   lambdas=rep(lambda,p),
                                                   start=rep(0,p),
                                                   control.solver = modifyList(do.call("control.solver.gen", 
                                                                                       list(solver = "eigen", method = 5L)),
                                                                               list(rowaug=T)),
                                                   prior.mean=rep(0,p)
)) 
# for n=p=10000: 12s
system.time(beta_eigen_cg_sparse <- lmridge.eigen(Xw=X_sparse,
                                                    yw=y,
                                                    lambdas=rep(lambda,p),
                                                    start=rep(0,p),
                                                    control.solver = 
                                                    modifyList(do.call("control.solver.gen", 
                                                                  list(solver = "eigen", method = 5L)),
                                                               list(rowaug=T)),
                                                    prior.mean=rep(0,p)
)) 
# for n=p=10000: 22s

# using LSCG Eigen solver: FASTEST OF ALL
system.time(beta_eigen_lscg_dense <- lmridge.eigen(Xw=X,
                                                 yw=y,
                                                 lambdas=rep(lambda,p),
                                                 start=rep(0,p),
                                                 control.solver = modifyList(do.call("control.solver.gen", 
                                                                                     list(solver = "eigen", method = 3L)),
                                                                             list(rowaug=T)),
                                                 prior.mean=rep(0,p)
)) 
# for n=p=10000: 10.5s 
system.time(beta_eigen_lscg_sparse <- lmridge.eigen(Xw=X_sparse,
                                                  yw=y,
                                                  lambdas=rep(lambda,p),
                                                  start=rep(0,p),
                                                  control.solver = modifyList(do.call("control.solver.gen", 
                                                                                      list(solver = "eigen", method = 3L)),
                                                                              list(rowaug=T)),
                                                  prior.mean=rep(0,p)
)) 
# for n=p=10000: 0.65s


# equivalent manual ridge regression:
# see https://stats.stackexchange.com/questions/76925/ridge-regression-in-r-where-coefficients-are-penalized-toward-numbers-other-than
system.time(beta_man2_dense <- lm.fit(x=rbind(X, diag(x=rep(sqrt(lambda), ncol(X)))),
       y=c(y, sqrt(lambda)*prior.mean) )$coefficients) 
# for n=p=10000: XXs with Intel MKL

# using glmnet - FASTEST OPTION?? - better than bigGlm? also scales better? ####
library(glmnet)
system.time({fit_glmnet1_dense <- glmnet(X, y, alpha=0, standardize = F,
                                         lambda=sd_y*lambda/n, nlambda=1,
                                         intercept = F, thresh = tol)
             beta_glmnet1_dense <- as.vector(coef(fit_glmnet1_dense, 
                               s = sd_y*lambda/n, exact = TRUE, x=X, y=y))[-1]
}) 
# for n=p=10000: 1.5s
system.time({fit_glmnet1_sparse <- glmnet(X_sparse, y, alpha=0, standardize = F,
                                         lambda=sd_y*lambda/n, nlambda=1,
                                         intercept = F, thresh = tol)
             beta_glmnet1_sparse <- as.vector(coef(fit_glmnet1_sparse, 
                                     s = sd_y*lambda/n, exact = TRUE, x=X_sparse, y=y))[-1]
}) 
# for n=p=10000: 0.28s


# # using glmnet.fit:
# library(glmnet)
# system.time({fit_glmnet2_dense <- glmnet:::glmnet.fit(X, y, 
#                                                       weights=rep(1,p),
#                                                       warm=list(beta=beta_init, a0=0),
#                                                       alpha=0, 
#                                          lambda=sd_y*lambda/n, 
#                                          intercept = F, 
#                                          thresh = tol)
#              beta_glmnet2_dense <- as.vector(fit_glmnet2_dense$beta)
# }) # 217s
# system.time({fit_glmnet2_sparse <- glmnet:::glmnet.fit(X_sparse, y, 
#                                                        weights=rep(1,p),
#                                                        warm=list(beta=beta_init, a0=0),
#                                                        alpha=0,
#                                           lambda=sd_y*lambda/n, 
#                                           intercept = F, 
#                                           thresh = tol)
#              beta_glmnet2_sparse <- as.vector(fit_glmnet2_sparse$beta)
# }) # XXs

# using glmnet & fit ridge regression via row augmentation:
system.time({fit_glmnet2_dense <- glmnet(rbind(X, Diagonal(x=rep(sqrt(lambda),ncol(X)))), 
                     c(y, rep(sqrt(lambda)*prior.mean)), lambda=0, nlambda=1, 
                     alpha=0, standardize = F, intercept = F, thresh = tol)
             beta_glmnet2_dense <- as.vector(coef(fit_glmnet2_dense, 
                               s = 0, exact = TRUE, x=X, y=y))[-1]
}) 
# for n=p=10000: 0.79s
system.time({fit_glmnet2_sparse <- glmnet(rbind(X_sparse, Diagonal(x=rep(sqrt(lambda),ncol(X)))), 
                                         c(y, rep(sqrt(lambda)*prior.mean)), lambda=0, nlambda=1, 
                                         alpha=0, standardize = F, intercept = F, thresh = tol)
             beta_glmnet2_sparse <- as.vector(coef(fit_glmnet2_sparse, 
                                     s = 0, exact = TRUE, x=X, y=y))[-1]
}) 
# for n=p=10000: 0.33s


# # using glmnet.fit & fit ridge regression via row augmentation:
# library(glmnet)
# system.time({fit_glmnet3_dense <- glmnet:::glmnet.fit(x=rbind(X, 
#                                                               Diagonal(x=rep(sqrt(lambda),ncol(X)))), 
#                                          y=c(y, rep(sqrt(lambda)*prior.mean)), 
#                                          weights=rep(1,p),
#                                          lambda=0, warm=list(beta=beta_init, a0=0),
#                                          alpha=0, intercept = F, thresh = tol) 
#              beta_glmnet3_dense <- as.vector(fit_glmnet3_dense$beta)
# }) 
# # for n=p=10000: XXs
# system.time({fit_glmnet3_sparse <- glmnet(rbind(X_sparse, 
#                                                 Diagonal(x=rep(sqrt(lambda),ncol(X)))), 
#                                           c(y, rep(sqrt(lambda)*prior.mean)), 
#                                           lambda=0,  warm=list(beta=beta_init, a0=0),
#                                           alpha=0, standardize = F, intercept = F, thresh = tol)
#              beta_glmnet3_sparse <- as.vector(fit_glmnet3_sparse$beta)
# }) 
# # for n=p=10000: too slow...


# using bigGlm:
library(glmnet)
system.time({fit_bigglm1_dense <- bigGlm(x=rbind(X, Diagonal(x=rep(sqrt(lambda),ncol(X)))), 
                                         y=c(y, rep(sqrt(lambda)*prior.mean)), 
                                         alpha=0, 
                                         standardize = F,
                                         lambda=0, 
                                         intercept = F, 
                                         thresh = tol)
beta_bigglm1_dense <- as.vector(fit_bigglm1_dense$beta)
}) 
# 0.97s for n=p=10000
system.time({fit_bigglm1_sparse <- bigGlm(x=rbind(X_sparse, Diagonal(x=rep(sqrt(lambda),ncol(X)))), 
                                          y=c(y, rep(sqrt(lambda)*prior.mean)), 
                                          alpha=0, 
                                          standardize=F,
                                          standardize.response = FALSE,
                                          lambda=0, 
                                          intercept = F, 
                                          thresh = tol,
                                          family = "gaussian", path = FALSE)
beta_bigglm1_sparse <- as.vector(fit_bigglm1_sparse$beta)
}) 
# 0.5s for n=p=10000


cbind(beta_man1_dense[1:10],
      beta_man1_sparse[1:10],
      # beta_man2_dense[1:10],
      beta_glmnet1_dense[1:10],
      beta_glmnet1_sparse[1:10],
      beta_glmnet2_dense[1:10],
      beta_glmnet2_sparse[1:10],
      # beta_glmnet3_dense[1:10],
      # beta_glmnet3_sparse[1:10],
      beta_bigglm1_dense[1:10],
      beta_bigglm1_sparse[1:10],
      beta_eigen_chol_dense[1:10],
      beta_eigen_chol_sparse[1:10],
      beta_eigen_cg_dense[1:10],
      beta_eigen_cg_sparse[1:10],
      beta_eigen_lscg_dense[1:10],
      beta_eigen_lscg_sparse[1:10] #
)



cbind(beta_man1_dense[1:10],
      beta_man1_sparse[1:10],
      beta_man2_dense[1:10],
      beta_glmnet1_dense[1:10],
      beta_glmnet1_sparse[1:10],
      beta_glmnet2_dense[1:10],
      beta_glmnet2_sparse[1:10],
      beta_glmnet3_dense[1:10],
      beta_glmnet3_sparse[1:10],
      beta_bigglm1_dense[1:10],
      beta_bigglm1_sparse[1:10],
      beta_eigen_chol_dense[1:10],
      beta_eigen_chol_sparse[1:10],
      beta_eigen_cg_dense[1:10],
      beta_eigen_cg_sparse[1:10],
      beta_eigen_lscg_dense[1:10],
      beta_eigen_lscg_sparse[1:10]
      )