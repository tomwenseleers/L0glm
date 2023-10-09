# Minimal tests of L0glm package plus comparison of solution quality to that returned by ordinis, abess & L0Learn ####

library(remotes)
# remotes::install_github("tomwenseleers/L0glm/L0glm")
library(export)
library(scales)
devtools::load_all(".")
# library(L0glm)
# remotes::install_github("jaredhuling/ordinis")
# library(ordinis)
# library(abess)
# remotes::install_github("hazimehh/L0Learn", ref="nng") # branch of L0Learn with nonnegativity constraints
# library(L0Learn)

# 1. SINGLE CHANNEL DECONVOLUTION ####

n <- p <- 500 # 10000 # to check performance use e.g. n=p=10000, otherwise use e.g. n=p=100-300
s <- 0.1 # sparsity, nr of timepoints with a peak 
k <- round(n*s)
mean_beta <- 2000
sd_logbeta <- 1
sd_noise <- 5 # SD of noise for Gaussian noise case
crit <- "bic" # "bic" # choose lambda based on what would be expected to optimise this IC
# options in order of increasing regularisation: aic, bic, gic, ebic, mbic

## 1.1 Single channel deconvolution assuming Poisson noise ####
### 1.1.1 L0glm solution ####
set.seed(1)
sim <- simulate_spike_train(n = n, 
                            p = p,
                            k = k,
                            family = "poisson", 
                            mean_beta = mean_beta,
                            sd_logbeta = sd_logbeta)
X <- sim$X
y <- sim$y
dispersion <- 1
family <- poisson(identity)
W <- 1/(y+1E-10) # approx 1/variance weights in case identity link Poisson is
                 # approximated using weighted least squares
beta_true <- sim$beta_true
L0overL1norm = sum(beta_true!=0)/sum(abs(beta_true))

# using lambda expected to optimise given IC
lambda_preset = preset.lambda(IC = crit, p=p, n=n) # options aic, bic, gic, ebic, mbic 
 
lambda_preset 

system.time(L0glm.out_preset <- L0glm.fit(X, y, 
                               family = family, 
                               lambda = lambda_preset, 
                               # ridge.center = start, # ridge.center
                               nonnegative = TRUE,
                               start = rep(50, ncol(sim$X)), # NULL, # rep(50, ncol(sim$X)), # 10-100 works best, coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                               solver = "eigen", # note: solver osqp & glmnet better than eigen for very small lambda values)
                               maxit = 200,
                               gamma = 2,
                               method = 7L, # 0.27s for n=p=10000 with eigen solver methods 7 & 8 (Cholesky), 0.61s with methods 5 & 6 (CG), 1.27s with methods 3, 4 (LSCG), 219s with methods 1 & 2 (QR), 4.8s for n=p=10000 with osqp solver, 2s for n=p=10000 with glmnet solver
                               thresh = 1E-5
)) # 0.02s for n=p=100, 0.04s for n=p=500, 0.25s for n=p=10000
L0glm.out_preset_coefs <- L0glm.out_preset$coefficients
sum(L0glm.out_preset_coefs!=0)  #55
points(sim$x, L0glm.out_preset_coefs, pch = 16, col= 'blue')
# graph2pdf(file="..//benchmarks/plots/1_1 single channel deconvolution Poisson.pdf", width=7, height=5)

mincoef <- min(abs(L0glm.out_preset_coefs)[L0glm.out_preset_coefs!=0]) # nonzero coefficient with smallest absolute value
do.call(cbind, get.metrics(beta = L0glm.out_preset_coefs,  # some metrics
                           X = X, y = y, 
                           family = family, 
                           beta_true = sim$beta_true))
# p      n   m   deviance deviance_true    loglik   MD_coefs MAD_coefs sMAPE_coefs bias_coefs
# [1,] 500 500 1 396.1279      103.3009 -1670.882 -0.2829917  39.18788  0.08039668   1.000992
#       R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true
# [1,] 0.9796187    105.6528      248.028      0.4763334      1.000992   0.9862615                50
#      support_size_estimated estimated_over_true_support_size TP FP  TN FN precision recall
# [1,]                     55                              1.1 43 12 438  7 0.7818182   0.86
#      specificity accuracy  F1_score      gic     ebic      bic      aic     mbic
# [1,]   0.9733333    0.962 0.8190476 3966.206 4025.371 3683.568 3451.764 4214.682

plot(x=L0glm.out_preset_coefs, y=sim$beta_true, pch=16, col='blue',
     main="L0 glm identity link Poisson with lambda preset to optimise BIC",
     xlab="Estimated coefficients", ylab="True coefficients")



# coefficients & metrics over a range of lambda values
# calculate max lambda where all coefs would be penalized to 0
maxlam <- maxlambda(X=X, weights=L0glm.out_preset$weights, 
                    z=y, dispersion=1, no.pen=NULL, coefs=L0glm.out_preset$coefficients,
                    family)
maxlam # 

# lambdas <- 10^seq(round(log10(lambda_preset))-1, log10(maxlam)+1, by=0.05) #
# length(lambdas) # 152 - maybe in benchmark we should use 20 values like for the other methods?
lambdas <- 10^seq(round(log10(lambda_preset))-1, log10(maxlam), length.out=100) #
length(lambdas) # 20 - like for other methods
coefmat <- Matrix(0, nrow=length(lambdas), ncol=ncol(sim$X))
warmstart <- 0.1 # % by which L0glm.out_preset will be used as start, here used lowest value that still gave stable coef paths

system.time(metrics <- data.frame(do.call(rbind, lapply(1:length(lambdas), function (i) { 
  lambda <- lambdas[[i]]
  # start <- L0glm.out_preset_coefs*warmstart+start_coefs(X, family)*(1-warmstart)
  # start[L0glm.out_preset_coefs==0] <- 0 # slightly faster & stable with this
  L0glm.out <- L0glm.fit(X, y, 
                             family = family, 
                             # family = gaussian(identity),
                             # weights = 1/(sim$y+0.1),
                             lambda = lambda, 
                             nonnegative = TRUE,
                             start = L0glm.out_preset_coefs, # start, # NULL, # L0glm.out_preset_coefs, # start, # 
                             # prior.mean = L0glm.out_preset_coefs*warmstart*0.2,
                             warm.start = warmstart,
                             solver = "eigen", # eigen, glmnet or osqp, # note: constrained solvers glmnet & osqp better than eigen for very small lambda values
                             method = 7L, 
                             # maxit = 200,
                             # gamma = 2, 
                             thresh = 1E-5,
                             relaxed = FALSE
                             ) 
  L0glm.out_coefs = L0glm.out$coefficients
  coefmat[i,] <<- as.vector(L0glm.out_coefs)
  return( unlist(c(lambda=lambda, get.metrics(beta = L0glm.out_coefs,
                                                       X = X, y = y,
                                                       family = family, 
                                                       beta_true = sim$beta_true) ))

          )
  } )))) # 0.81s for n=p=500, 11s for n=p=10000 & 20 lambda values with eigen solver methods 7 (Cholesky)
matplot(x=log10(metrics$lambda), y=as.matrix(coefmat), type="l", lty=1,
        xlab="Log10(lambda)", ylab="Coefficient estimates", 
        main="Regularisation path L0 Poisson GLM")
abline(v=log10(maxlam), col="red") # max lambda where all coefs would be penalized to 0
abline(h=as.vector(sim$beta_true[sim$beta_true!=0]), col=alpha(1:6,0.2))


# coefficients for solution that optimises given criterion & metrics
best = which.min(metrics[,crit])
L0glm.out_best_coefs = coefmat[best,]
plot(x=L0glm.out_best_coefs, y=sim$beta_true, pch=16, col='blue',
     main="L0 Poisson GLM with tuned lambda",
     xlab="Estimated coefficients", ylab="True coefficients")
abline(a=0, b=1, col="lightgrey")

do.call(cbind, get.metrics(beta = L0glm.out_best_coefs,  # some metrics
                           X = X, y = y, 
                           family = family, 
                           beta_true = sim$beta_true))
#        p   n m deviance deviance_true    loglik   MD_coefs MAD_coefs sMAPE_coefs bias_coefs
# [1,] 500 500 1 402.0993      97.53431 -1672.868 -0.3098657  37.86387  0.07571449    1.00061
#       R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true
# [1,] 0.9804216    107.3606     240.3617      0.4418009       1.00061   0.9870698                50
#      support_size_estimated estimated_over_true_support_size TP FP  TN FN precision recall
# [1,]                     54                             1.08 43 11 439  7 0.7962963   0.86
#      specificity accuracy  F1_score      gic     ebic      bic      aic     mbic
# [1,]   0.9755556    0.964 0.8269231 3958.824 4016.913 3681.325 3453.736 4202.782


plot(log10(metrics$lambda), metrics$deviance, type="l")
metrics$lambda[which.min(metrics$deviance)] # 1 for n=p=10000, 0.1 for n=p=100 # NOTE this is deviance relative to observed signal, not deviance relative to ground truth though
plot(log10(metrics$lambda), metrics$deviance_true, type="l")
metrics$lambda[which(metrics$deviance_true==min(metrics$deviance_true))] # XX for n=p=10000, 10 for n=p=100 # NOTE this is deviance relative to observed signal, not deviance relative to ground truth though
plot(log10(metrics$lambda), metrics$loglik, type="l")
metrics$lambda[which(metrics$loglik==max(metrics$loglik))] # 1 for n=p=10000, 0.1 for n=p=100
plot(log10(metrics$lambda), metrics$MD_coefs, type="l")
plot(log10(metrics$lambda), metrics$MAD_coefs, type="l")
min(metrics$lambda[which(round(metrics$MAD_coefs,0)==min(round(metrics$MAD_coefs,0)))]) # 2.238721 for n=p=10000, 0.5623413 for n=p=100 
plot(log10(metrics$lambda), metrics$sMAPE_coefs, type="l") # good one, optimal for lambda between 1 & 10
min(metrics$lambda[which(round(metrics$sMAPE_coefs,2)==min(round(metrics$sMAPE_coefs,2)))]) # 11.22018 for n=p=10000, 0.5623413 for n=p=100
plot(log10(metrics$lambda), metrics$bias_coefs, type="l")
metrics$lambda[which(abs(metrics$bias_coefs-1)==min(abs(metrics$bias_coefs-1)))] # 
plot(log10(metrics$lambda), metrics$R2_coefs, type="l")
plot(log10(metrics$lambda), metrics$estimated_over_true_support_size, type="l") # good one, optimal for lambda between 1 & 10
range(metrics$lambda[which(abs(metrics$estimated_over_true_support_size-1)==min(abs(metrics$estimated_over_true_support_size-1)))]) # 2.238721 for n=p=10000, 0.5623413-25.1188643 for n=p=100
plot(log10(metrics$lambda), metrics$TP, type="l")
plot(log10(metrics$lambda), metrics$FP, type="l")
plot(log10(metrics$lambda), metrics$TN, type="l")
plot(log10(metrics$lambda), metrics$FN, type="l")
plot(log10(metrics$lambda), metrics$precision, type="l")
plot(log10(metrics$lambda), metrics$recall, type="l")
plot(log10(metrics$lambda), metrics$specificity, type="l")
plot(log10(metrics$lambda), metrics$accuracy, type="l") # good one, optimal for lambda between 1 & 10
min(metrics$lambda[which(metrics$accuracy==max(metrics$accuracy))]) # 12.58925 for n=p=10000, 0.5623413 for n=p=100
range(metrics$lambda[which(metrics$accuracy==max(metrics$accuracy))]) # 12.58925 12.58925 for n=p=10000, 0.5623413-25.1188643 for n=p=100
plot(log10(metrics$lambda), metrics$F1_score, type="l")
metrics$lambda[which.max(metrics$F1_score)] # 12.58925 for n=p=10000, 0.5623413 for n=p=100
plot(log10(metrics$lambda), metrics$aic, type="l")
min(metrics$lambda[which(metrics$aic==min(metrics$aic))]) # 1 for n=p=10000, 0.5623413 for n=p=10000 (optimal at same value as where loglik is optimal)
plot(log10(metrics$lambda), metrics$bic, type="l")
max(metrics$lambda[which(metrics$bic==min(metrics$bic))]) # 1.584893 for n=p=10000, 0.5623413 for n=p=100
plot(log10(metrics$lambda), metrics$ebic, type="l")
max(metrics$lambda[which(metrics$ebic==min(metrics$ebic))]) # 1.584893 for n=p=10000, 0.5623413 for n=p=100
plot(log10(metrics$lambda), metrics$gic, type="l")
max(metrics$lambda[which(metrics$gic==min(metrics$gic))]) # 1.584893 for n=p=10000, 0.5623413 for n=p=100
abline(v=log10(max(metrics$lambda[which(metrics$gic==min(metrics$gic))])), col="red") 
plot(log10(metrics$lambda), metrics$mbic, type="l")
min(metrics$lambda[which(metrics$mbic==min(metrics$mbic))]) # 1.995262 for n=p=10000, 0.5623413 for n=p=100




### 1.1.2 Comparison of solution quality to that returned by ordinis (adaptive lasso, MCP, SCAD), glmnet (lasso), abess & L0Learn ####

#### nonnegative weighted lasso regression in glmnet ####
library(glmnet)
system.time(fit_glmnet <- glmnet(x = X, 
                                 y = y, 
                                 weights = W,
                                 family = "gaussian", # poisson(identity) should be supported but gave an error
                                 nlambda = 20L,
                                 intercept = F, 
                                 standardize = F, 
                                 # thresh = 1E-10,
                                 # maxit=100000,
                                 alpha = 1, # lasso (lower value, ie elastic net, made solution worse)
                                 lower.limits = 0)) # for nonnegativity constraints
# 0.01s for n=p=100, 0.01s for n=p=500, 0.05s for n=p=10000
plot(fit_glmnet, xvar='lambda')
abline(h=as.vector(sim$beta_true[sim$beta_true!=0]), col=alpha(1:6,0.2))

LLs <- logLik(fit_glmnet, weights=W, X=X, y=y)
edfs = attr(LLs, "df")
IC = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,crit]
which.min(IC) # 

bestlambda_lasso_glmnet = fit_glmnet$lambda[which.min(IC)]
coefs_lasso_glmnet = coef(fit_glmnet, s = bestlambda_lasso_glmnet)

sum(coefs_lasso_glmnet!=0) # not sparse enough
plot(x=coefs_lasso_glmnet[-1], y=sim$beta_true, pch=16, col='blue',
     main="glmnet nonnegative weighted LASSO",
     xlab="Estimated coefficients", ylab="True coefficients")


#### adaptive LASSO using glmnet ####

# since glmnet is so much faster than ordinis & supports sparse
# covariate matrices, let's use that to fit an adaptive LASSO model, cf
# https://rpubs.com/kaz_yos/alasso
# note that for the first step we use LASSO rather than ridge coefficients, i.e. we use
# coefs_lasso_glmnet and then in a second step use a glmnet LASSO fit with
# adaptive penalty weights penalty.factor = 1 / (abs(coefs_lasso_glmnet))+delta)
# with delta=1E-5
# NOTE: If this would be iterated further this would be a way to approximate L0 norm
# penalized regression (cf. Nuel 2017 talk)

system.time({
  # we start with a regular LASSO fit, with lambda tuned according to our chosen IC
  fit_lasso_glmnet <- glmnet(x = X, 
                                 y = y, 
                                 weights = W,
                                 family = "gaussian", # poisson(identity) should be supported but gave an error
                                 nlambda = 20L,
                                 intercept = F, 
                                 standardize = F, 
                                 # thresh = 1E-10,
                                 # maxit=100000,
                                 alpha = 1, # lasso (lower value, ie elastic net, made solution worse)
                                 lower.limits = 0) # for nonnegativity constraints
LLs <- logLik(fit_lasso_glmnet, weights=W, X=X, y=y)
edfs = attr(LLs, "df")
IC = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,crit]
bestlambda_lasso_glmnet = fit_lasso_glmnet$lambda[which.min(IC)]
coefs_lasso_glmnet = coef(fit_lasso_glmnet, s = bestlambda_lasso_glmnet)

# in second step we do another LASSO fit with adaptive penalty weights
delta = 1E-5
fit_alasso_glmnet <- glmnet(x = X, 
                           y = y, 
                           weights = W,
                           family = "gaussian", # poisson(identity) should be supported but gave an error
                           nlambda = 20L,
                           penalty.factor = 1 / (abs(coefs_lasso_glmnet[-1])+delta), # adaptive penalty weights
                           intercept = F, 
                           standardize = F, 
                           alpha = 1, # lasso (lower value, ie elastic net, made solution worse)
                           lower.limits = 0) # for nonnegativity constraints
LLs <- logLik(fit_alasso_glmnet, weights=W, X=X, y=y)
edfs = attr(LLs, "df")
IC = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,crit]
bestlambda_alasso_glmnet = fit_alasso_glmnet$lambda[which.min(IC)]
coefs_alasso_glmnet = coef(fit_alasso_glmnet, s = bestlambda_lasso_glmnet)
})
# 0.02s for n=p=500, 0.11s for n=p=10000

plot(fit_alasso_glmnet, xvar='lambda')
abline(h=as.vector(sim$beta_true[sim$beta_true!=0]), col=alpha(1:6,0.2))

sum(coefs_alasso_glmnet!=0) # 1092 for n=p=10000, pretty good
plot(x=coefs_alasso_glmnet[-1], y=sim$beta_true, pch=16, col='blue',
     main="glmnet nonnegative weighted adaptive LASSO",
     xlab="Estimated coefficients", ylab="True coefficients")





#### nonnegative L0Learn best subset fit after glmnet nonnegative LASSO ####
library(L0Learn)
selcoefs = which(coefs_lasso_glmnet[-1]!=0)
p_selcoefs = length(selcoefs)
system.time(fit_nnL0Learn_selcoefs <- L0Learn.fit(x=as.matrix(Diagonal(x=sqrt(W)) %*% X[,selcoefs,drop=F]),  # note: no support for sparse matrices
                                         y=as.vector(y*sqrt(W)),
                                         loss = "SquaredError",
                                         penalty="L0", algorithm="CD",
                                         maxSuppSize=(p_selcoefs-1), 
                                         intercept=FALSE, 
                                         nLambda=100, # NOTE: 100 gave much better results than 20
                                         # activeSetNum = 3, maxSwaps = 100,
                                         autoLambda=TRUE
                                         # lambdaGrid=list(10^seq(-5,5,length.out=nlam)/ # DECREASING sequence recommended, but INCREASING works better under nonnegativity, see https://cran.r-project.org/web/packages/L0Learn/vignettes/L0Learn-vignette.html
                                         #                   (norm(sim$y*sqrt(1/(sim$y+1)),"2")^2)) 
))  # lambda of 1 is scaled by 1/(norm(weighted y,"2")) because L0Learn scales y by L2 norm of y
# 0.29s for n=p=500
plot(fit_nnL0Learn_selcoefs)



LLs = logLik(fit_nnL0Learn_selcoefs, weights=W, X[,selcoefs,drop=F], y)
edfs = attr(LLs, "df")
LLs = as.vector(LLs)
IC = get.IC(loglik = LLs, n = n, p = p_selcoefs, edf = edfs)[,crit]
coefs_nnL0Learn_selcoefs = coef(fit_nnL0Learn_selcoefs,
                       lambda = fit_nnL0Learn_selcoefs$lambda[[1]][which.min(IC)],
                       gamma = 0)
coefs_nnL0Learn_after_lasso = Matrix(0, p)
coefs_nnL0Learn_after_lasso[selcoefs] = coefs_nnL0Learn_selcoefs
sum(coefs_nnL0Learn_after_lasso!=0) # 52

plot(x=coefs_nnL0Learn_after_lasso, y=sim$beta_true, pch=16, col='blue',
     main="weighted nnL0Learn after nonnegative LASSO",
     xlab="Estimated coefficients", ylab="True coefficients")


#### # full nonnegative L0Learn best subset fit ####
# # quite good, but slow for n=p=10000
system.time(fit_nnL0Learn <- L0Learn.fit(x=as.matrix(Diagonal(x=sqrt(W)) %*% X),  # note: no support for sparse matrices
                                         y=as.vector(y*sqrt(W)),
                                         loss = "SquaredError",
                                         penalty="L0", algorithm="CD",
                                         maxSuppSize=round(p/2),
                                         intercept=FALSE,
                                         nLambda=100, # NOTE: to make timings fair we should probably use 100 or 20 lambda values everywhere?
                                         # activeSetNum = 3, maxSwaps = 100,
                                         autoLambda=TRUE
                                         # lambdaGrid=list(10^seq(-5,5,length.out=nlam)/ # DECREASING sequence recommended, but INCREASING works better under nonnegativity, see https://cran.r-project.org/web/packages/L0Learn/vignettes/L0Learn-vignette.html
                                         #                   (norm(sim$y*sqrt(1/(sim$y+1)),"2")^2))
))  # lambda of 1 is scaled by 1/(norm(weighted y,"2")) because L0Learn scales y by L2 norm of y
# 0.86s for n=p=500
plot(fit_nnL0Learn)

LLs = logLik(fit_nnL0Learn, weights=W, X, y)
edfs = attr(LLs, "df")
LLs = as.vector(LLs)
IC = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,crit]
coefs_nnL0Learn = coef(fit_nnL0Learn,
                       lambda = fit_nnL0Learn$lambda[[1]][which.min(IC)],
                       gamma = 0)
sum(coefs_nnL0Learn!=0) # very close
plot(x=coefs_nnL0Learn, y=sim$beta_true, pch=16, col='blue',
     main="nnL0Learn weighted least squares best subset",
     xlab="Estimated coefficients", ylab="True coefficients")



#### abess best subset fit after nonnegative LASSO ####

library(abess)
# very good solution, and fast for small problems, but
# becomes very slow for large problems, e.g. n=p=10000
selcoefs = which(coefs_lasso_glmnet[-1]!=0)
p_selcoefs = length(selcoefs)
system.time(fit_abessafterlasso_selcoefs <- abess(x = X[,selcoefs,drop=F], 
                                                  y = y, 
                                                  weights = W,  # approx 1/variance weights to approximate family="poisson", which is not supported
                                                  family = "gaussian",
                                                  tune.path = "sequence", 
                                                  support.size = c(1:(p_selcoefs-1)),
                                                  lambda = 0, 
                                                  warm.start = TRUE, # FALSE worked better here, but much slower
                                                  tune.type = crit))  # options bic, ebic, aic or cv
# 0.06s for n=p=500 
# note: abess currently has no support for nonnegativity constraints, so 
# abess after nonnegative LASSO works much better
plot(fit_abessafterlasso_selcoefs, type="tune")  # optimal support size 
coefs_abessafterlasso_selcoefs = as.matrix(extract(fit_abessafterlasso_selcoefs)$beta) # coefficient matrix for best subset
coefs_abessafterlasso = Matrix(rep(0, p), sparse=T)
coefs_abessafterlasso[selcoefs] = coefs_abessafterlasso_selcoefs
coefs_abessafterlasso[coefs_abessafterlasso<0] = 0
sum(coefs_abessafterlasso!=0) # pretty good
plot(x=coefs_abessafterlasso, y=sim$beta_true, pch=16, col='blue',
     main="abess best subset weighted least squares\nafter nonnegative LASSO",
     xlab="Estimated coefficients", ylab="True coefficients")


#### # full abess best subset fit (with EBIC criterion) ####
# # note: currently no support for nonnegativity constraints
# # solution not great, so marked out for now & slow for n=p=10000
# system.time(fit_abess <- abess(x = X, 
#                                y = y, 
#                                weights = W,  # approx 1/variance weights to approximate family="poisson", which is not supported
#                                family = "gaussian",
#                                tune.path = "sequence", 
#                                support.size = c(1:round(p*0.5)),
#                                lambda=0, 
#                                warm.start = TRUE, # TRUE worked best here
#                                tune.type = "ebic")) # had to use ebic here for optimal performance, XXs for n=p=XX, options bic, ebic, aic or cv
# # fast for small problems, but becomes slow for n=p=10000
# plot(fit_abess, type="tune")  # optimal support size 
# coefs_abess = as.matrix(extract(fit_abess)$beta) # coefficient matrix for best subset
# coefs_abess[coefs_abess<0] = 0
# sum(coefs_abess!=0) # 11, quite OK
# plot(x=coefs_abess, y=sim$beta_true, pch=16, col='blue',
#      main="abess best subset weighted least squares",
#      xlab="Estimated coefficients", ylab="True coefficients")



#### nonnegative weighted MCP penalized regression in ordinis ####

library(ordinis)
# NOTE: ncvreg also has MCP penalized models, but not with nonnegativity constraints
system.time(fit_MCP <- ordinis(X, y, 
                           weights=W, # approx 1/variance Poisson weights
                           nlambda = 100L, # 1000L gives much better solution than 100L, 20L doesn't work & coefficient path looks very different
                           penalty = "mcp", # using MCP penalty
                           intercept = FALSE,
                           lower.limits = rep(0, p), # force all coefficients to be positive
                           standardize = FALSE,
                           # maxit = 1000,
                           # tol = 1E-100,
                           # maxit.irls = 1000,
                           # tol.irls = 1E-100,
                           # gamma = 1.1, # ifelse(penalty == "scad", 3.7, 1.4),
                           alpha = 1))  # ie no elastic net
# 12.25s for n=p=500, XXs for n=p=10000
plot(fit_MCP)
abline(h=as.vector(sim$beta_true[sim$beta_true!=0]), col=alpha(1:6,0.2)) # note : ordinis uses different colours, so cols don't match now


LLs = logLik2(fit_MCP, weights=W, X, y)
edfs = attr(LLs, "df")
LLs = as.vector(LLs)
IC = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,crit]
coefs_MCP_ordinis = fit_MCP$beta[,which.min(IC)]
sum(coefs_MCP_ordinis!=0) # not sparse enough
plot(x=coefs_MCP_ordinis[-1], y=sim$beta_true, pch=16, col='blue',
     main="ordinis nonnegative MCP",
     xlab="Estimated coefficients", ylab="True coefficients")


#### nonnegative weighted SCAD penalized regression in ordinis ####
# NOTE: ncvreg also has SCAD penalized models, but not with nonnegativity constraints
system.time(fit_SCAD <- ordinis(X, 
                                y, 
                               weights=W, # approx 1/variance Poisson weights
                               nlambda = 20L, # 1000L gives much better solution than 100L or 20L & coefficient path looks very different
                               penalty = "scad", # using SCAD penalty
                               intercept = FALSE,
                               lower.limits = rep(0, p), # force all coefficients to be positive
                               standardize = FALSE,
                               # maxit = 1000,
                               # tol = 1E-100,
                               # maxit.irls = 1000,
                               # tol.irls = 1E-100,
                               # gamma = 100, # ifelse(penalty == "scad", 3.7, 1.4),
                               alpha = 1))  # ie no elastic net
# 4.75s for n=p=500
plot(fit_SCAD)
abline(h=as.vector(sim$beta_true[sim$beta_true!=0]), col=alpha(1:6,0.2)) # note : ordinis uses different colours, so cols don't match now
LLs = logLik2(fit_SCAD, weights=W, X, y)
edfs = attr(LLs, "df")
LLs = as.vector(LLs)
IC = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,crit]
coefs_SCAD_ordinis = fit_SCAD$beta[,which.min(IC)]
sum(coefs_SCAD_ordinis!=0) # not sparse enough
plot(x=coefs_SCAD_ordinis[-1], y=sim$beta_true, pch=16, col='blue',
     main="ordinis nonnegative SCAD",
     xlab="Estimated coefficients", ylab="True coefficients")


#### # nonnegative weighted lasso regression in ordinis ####
# system.time(fit_lasso_ordinis <- ordinis(X, 
#                                          y, 
#                                  weights=W, # approx 1/variance Poisson weights
#                                  nlambda = 20L, # 1000L gives much better solution than 100L or 20L & coefficient path looks very different
#                                  lambda.min.ratio = 0.0001,
#                                  penalty = "lasso", # using LASSO penalty
#                                  intercept = FALSE,
#                                  lower.limits = rep(0, p), # force all coefficients to be positive
#                                  standardize = FALSE,
#                                  maxit = 1000,
#                                  # tol = 1E-100,
#                                  # maxit.irls = 1000,
#                                  # tol.irls = 1E-100,
#                                  # gamma = 100, # ifelse(penalty == "scad", 3.7, 1.4),
#                                  alpha = 1))  # ie no elastic net
# # XXs for n=p=100, 522s for n=p=10000
# plot(fit_lasso_ordinis)
# abline(h=as.vector(sim$beta_true[sim$beta_true!=0]), col=alpha(1:6,0.2)) # note : ordinis uses different colours, so cols don't match now

# LLs = logLik2(object=fit_lasso_ordinis, weights=W, X=X, y=y)
# edfs = attr(LLs, "df")
# LLs = as.vector(LLs)
# IC = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,crit]
# which.min(IC) # 
# bestlambda_lasso_ordinis = fit_lasso_ordinis$lambda[which.min(IC)]
# bestlambda_lasso_ordinis # 
# coefs_lasso_ordinis = fit_lasso_ordinis$beta[,which.min(IC)]
# sum(coefs_lasso_ordinis!=0) # not sparse enough
# plot(x=coefs_lasso_ordinis[-1], y=sim$beta_true, pch=16, col='blue',
#      main="ordinis nonnegative weighted LASSO",
#      xlab="Estimated coefficients", ylab="True coefficients")


#### # nonnegative weighted adaptive lasso regression in ordinis ####
# system.time(fit_alasso_ordinis <- ordinis(X, y, 
#                                   weights=W, # approx 1/variance Poisson weights
#                                   nlambda = 20L, # 1000L gives much better solution than 100L or 20L & coefficient path looks very different
#                                   penalty = "alasso", # using adaptive LASSO penalty
#                                   intercept = FALSE,
#                                   lower.limits = rep(0, p), # force all coefficients to be positive
#                                   standardize = FALSE,
#                                   # maxit = 1000,
#                                   # tol = 1E-100,
#                                   # maxit.irls = 1000,
#                                   # tol.irls = 1E-100,
#                                   # gamma = 100, # ifelse(penalty == "scad", 3.7, 1.4),
#                                   alpha = 1))  # ie no elastic net
# # XXs for n=p=100, XXs for n=p=10000
# plot(fit_alasso_ordinis)
# abline(h=as.vector(sim$beta_true[sim$beta_true!=0]), col=alpha(1:6,0.2)) # note : ordinis uses different colours, so cols don't match now
# 
# LLs = logLik2(fit_alasso_ordinis, weights=W, X, y)
# edfs = attr(LLs, "df")
# LLs = as.vector(LLs)
# IC = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,crit]
# coefs_alasso_ordinis = fit_alasso_ordinis$beta[,which.min(IC)]
# sum(coefs_alasso_ordinis!=0) # 9 - slightly too sparse
# plot(x=coefs_alasso_ordinis[-1], y=sim$beta_true, pch=16, col='blue',
#      main="ordinis nonnegative weighted adaptive LASSO",
#      xlab="Estimated coefficients", ylab="True coefficients")



#### SUMMARY TABLE WITH BENCHMARKS ####
benchmarks_1_1 = rbind(
      data.frame(method="L0glm_tuned_lambda", do.call(cbind, get.metrics(beta = L0glm.out_best_coefs,  # some metrics
                                                                     X = X, y = y, 
                                                                     family = family, 
                                                                     beta_true = sim$beta_true))),
  
      data.frame(method="L0glm_preset_lambda", do.call(cbind, get.metrics(beta = L0glm.out_preset_coefs,  # some metrics
                           X = X, y = y, 
                           family = family, 
                           beta_true = sim$beta_true))),
      data.frame(method="nnalasso_glmnet", do.call(cbind, get.metrics(beta = coefs_alasso_glmnet,  # some metrics
                                                                      X = cbind(1, X), y = y, 
                                                                      weights = W,
                                                                      family = gaussian(identity), 
                                                                      beta_true = c(0, as.vector(sim$beta_true))) )),
      data.frame(method="nnL0Learn_after_glmnet_nnlasso", do.call(cbind, get.metrics(beta = coefs_nnL0Learn_after_lasso,  # some metrics
                                                                                     X = X, y = y, 
                                                                                     weights = W,
                                                                                     family = gaussian(identity),
                                                                                     beta_true = as.vector(sim$beta_true)) )),
      data.frame(method="abess_after_glmnet_nnlasso", do.call(cbind, get.metrics(beta = coefs_abessafterlasso,  # some metrics
                                                                                 X = X, y = y, 
                                                                                 weights = W,
                                                                                 family = gaussian(identity),
                                                                                 beta_true = as.vector(sim$beta_true)) )),
      data.frame(method="nnL0Learn", do.call(cbind, get.metrics(beta = coefs_nnL0Learn,  # some metrics
                                                                      X = X, y = y, 
                                                                      weights = W,
                                                                      family = gaussian(identity), 
                                                                      beta_true = as.vector(sim$beta_true)) )),
      data.frame(method="nnlasso_glmnet", do.call(cbind, get.metrics(beta = coefs_lasso_glmnet,  # some metrics
                                                                     X = cbind(1, X), y = y, 
                                                                     weights = W,
                                                                     family = gaussian(identity), 
                                                                     beta_true = c(0, as.vector(sim$beta_true))) )),
      data.frame(method="nnMCP_ordinis", do.call(cbind, get.metrics(beta = coefs_MCP_ordinis,  # some metrics
                                                                    X = cbind(1, X), y = y, 
                                                                    weights = W,
                                                                    family = gaussian(identity),
                                                                    beta_true = c(0, as.vector(sim$beta_true))) )),
      data.frame(method="nnSCAD_ordinis", do.call(cbind, get.metrics(beta = coefs_SCAD_ordinis,  # some metrics
                                                             X = cbind(1, X), y = y, 
                                                             weights = W,
                                                             family = gaussian(identity), 
                                                             beta_true = c(0, as.vector(sim$beta_true))) ))
      # data.frame(method="nnlasso_ordinis", do.call(cbind, get.metrics(beta = coefs_lasso_ordinis,  # some metrics
      #                                                                 X = cbind(1, X), y = y, 
      #                                                                 weights = W,
      #                                                                 family = gaussian(identity),
      #                                                                 beta_true = c(0, as.vector(sim$beta_true))) )),
      # data.frame(method="nnalasso_ordinis", do.call(cbind, get.metrics(beta = coefs_alasso_ordinis,  # some metrics
      #                                                       X = cbind(1, X), y = y, 
      #                                                       weights = W,
      #                                                       family = gaussian(identity), 
      #                                                       beta_true = c(0, as.vector(sim$beta_true))) )),
      # data.frame(method="abess", do.call(cbind, get.metrics(beta = coefs_abess,  # some metrics
      #                                                                X = X, y = y, 
      #                                                                weights = W,
      #                                                                family = gaussian(identity), 
      #                                                                beta_true = as.vector(sim$beta_true)) ))
      )
benchmarks_1_1 = benchmarks_1_1[order(benchmarks_1_1[,crit]),] # sort by IC
benchmarks_1_1
# CONCLUSION: L0glm with tuned lambda best followed by L0glm with preset lambda,
# followed by L0Learn after lasso, alasso & L0Learn; 
# all other methods much worse

# For n=p=500

#                           method   p   n m     deviance deviance_true    loglik    MD_coefs
# 1             L0glm_tuned_lambda 500 500 1 3.816367e+02  8.174847e+01 -1664.194  -0.2732649
# 2            L0glm_preset_lambda 500 500 1 3.901860e+02  9.186823e+01 -1669.468  -0.2586797
# 4 nnL0Learn_after_glmnet_nnlasso 501 500 1 1.911327e+04  5.047348e+10 -2217.407  -6.2095115
# 3                nnalasso_glmnet 502 500 1 1.887068e+04  5.047348e+10 -2216.214  -6.5259568
# 6                      nnL0Learn 501 500 1 2.056724e+04  5.047348e+10 -2234.737  -6.6674122
# 7                 nnlasso_glmnet 502 500 1 1.880836e+04  5.047332e+10 -2272.387  -6.4909290
# 8                  nnMCP_ordinis 502 500 1 2.218307e+04  5.047349e+10 -2320.644  -8.1817751
# 9                 nnSCAD_ordinis 502 500 1 3.644955e+04  5.046555e+10 -2462.794 -11.5835933
# 5     abess_after_glmnet_nnlasso 501 500 1 2.541601e+09  3.600837e+10 -5168.888  25.7496713
#   MAD_coefs sMAPE_coefs bias_coefs  R2_coefs support_size_true support_size_estimated
# 1  41.71338  0.08770819  0.9980317 0.9785675                50                     55
# 2  44.89718  0.09978243  0.9988826 0.9778247                50                     56
# 4 117.90380  0.17375413  0.9651048 0.8999311                50                     52
# 3 115.79339  0.15415604  0.9810797 0.8915773                50                     55
# 6 205.63036  0.23467638  0.8126596 0.6770243                50                     51
# 7 144.37886  0.36529360  0.9994436 0.8652952                50                    112
# 8 333.05557  0.41517219  1.1096016 0.4792685                50                    119
# 9 366.28820  0.47824117  1.1203732 0.3898567                50                    137
# 5 215.22493  0.22115776  0.7853456 0.6375838                50                     54
#   estimated_over_true_support_size TP FP  TN FN precision recall specificity  accuracy  F1_score
# 1                             1.10 42 13 437  8 0.7636364   0.84   0.9711111 0.9580000 0.8000000
# 2                             1.12 41 15 435  9 0.7321429   0.82   0.9666667 0.9520000 0.7735849
# 4                             1.04 30 22 428 20 0.5769231   0.60   0.9511111 0.9160000 0.5882353
# 3                             1.10 35 20 431 15 0.6363636   0.70   0.9556541 0.9301397 0.6666667
# 6                             1.02 22 29 421 28 0.4313725   0.44   0.9355556 0.8860000 0.4356436
# 7                             2.24 39 73 378 11 0.3482143   0.78   0.8381375 0.8323353 0.4814815
# 8                             2.38 39 80 371 11 0.3277311   0.78   0.8226164 0.8183633 0.4615385
# 9                             2.74 42 95 356  8 0.3065693   0.84   0.7893570 0.7944112 0.4491979
# 5                             1.08 25 29 421 25 0.4629630   0.50   0.9355556 0.8920000 0.4807692
#         gic      ebic       bic       aic      mbic
# 1  3952.829  4011.994  3670.191  3438.388  4201.305
# 2  3974.732  4034.973  3686.955  3450.937  4227.726
# 4  5036.743  5093.775  4764.189  4540.815  5276.202
# 3  5068.632  5128.911  4780.446  4544.428  5321.664
# 6  5060.044  5116.000  4792.633  4573.473  5294.985
# 7  5828.542  5950.178  5247.025  4770.774  6339.126
# 8  6004.581  6133.752  5387.041  4881.288  6546.794
# 9  6493.375  6641.922  5783.204  5201.588  7116.920
# 5 10962.418 11021.602 10679.579 10447.776 11210.913

# TODO: probably best to package some of the above in an individual sparse learning
# function sparseglm with identical arguments (e.g. method=lasso, alasso, mcp, scad, L0learn,
# L0learn_after_lasso, abess, abess_after_lasso, etc with arguments X=, y=, 
# family=identity link gaussian or poisson, nonnegative=T/F)
# and coefficients as output & to then make a benchmark function with the same arguments
# plus fitfunc, which could be a list of fit functions you would like to try plus a
# benchmark_sim function that also runs the simulation function for given
# sim_options simulation option parameters (the latter also best packaged in
# a single simulate() function with sim=1 e.g. spike train with Gaussian noise,
# sim=2 spike train with Poisson noise, sim=3 orthogonal problem etc).
# That benchmark function can then also record the timings & output that
# as a column of the dataframe above.
# TODO: do the same for Gaussian case below (just set W=rep(1,p)) & for multichannel
# and also run similar benchmarks for orthogonal case & high dimensional cases
# (I had a try in test_high dimensional, but didn't get very far - there
# ordinis gave a better result right now, but maybe I was doing something wrong).


## 1.2 Single channel deconvolution assuming Gaussian noise ####
### 1.2.1 L0glm solution ####
set.seed(1)
sim <- simulate_spike_train(n = n, 
                            p = p,
                            k = k,
                            family = "gaussian", 
                            mean_beta = mean_beta,
                            sd_logbeta = sd_logbeta,
                            sd_noise = sd_noise)
X <- sim$X
y <- sim$y
family <- gaussian(identity)
dispersion <- sd_noise^2
W <- rep(1,p) # with Gaussian equal observation weights

# using lambda expected to optimise given IC
lambda_preset <- preset.lambda(IC = crit, p=p, n=n) # options aic, bic, gic, ebic, mbic 
lambda_preset # lambda that would be expected to optimise given IC

system.time(L0glm.out_preset <- L0glm.fit(X, y, 
                               family = family, 
                               lambda = lambda_preset, 
                               nonnegative = TRUE,
                               start = NULL, # rep(1, ncol(X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                               solver = "eigen", # note: solver osqp better than eigen for very small lambda values))
                               method = 7L # Cholesky
                              )) 
# 0.03s
L0glm.out_preset_coefs = L0glm.out_preset$coefficients
points(sim$x, L0glm.out_preset_coefs, pch = 16, col= 'blue')
sum(L0glm.out_preset_coefs!=0) # 54

# some metrics
do.call(cbind, get.metrics(beta = L0glm.out_preset_coefs, 
                           X = X, y = y, 
                           family = family, 
                           beta_true = sim$beta_true))
# graph2pdf(file="..//benchmarks/plots/1_2 single channel deconvolution Gaussian.pdf", width=7, height=5)

#        p   n m deviance deviance_true    loglik    MD_coefs MAD_coefs sMAPE_coefs bias_coefs
# [1,] 501 500 1 10189.95      1390.415 -1017.107 -0.02009234  6.407866  0.03245251   1.001921
#      R2_coefs support_size_true support_size_estimated estimated_over_true_support_size TP FP  TN FN
# [1,] 0.998688                50                     54                             1.08 48  6 444  2
#      precision recall specificity accuracy  F1_score      gic    ebic      bic      aic     mbic
# [1,] 0.8888889   0.96   0.9866667    0.984 0.9230769 2658.855 2718.04 2376.017 2144.213 2907.351

# coefficients & metrics over a range of lambda values
# calculate max lambda where all coefs would be penalized to 0
maxlam <- maxlambda(X=X, weights=L0glm.out_preset$weights, 
                    z=y, dispersion=L0glm.out_preset$dispersion, 
                    no.pen=NULL, coefs=L0glm.out_preset$coefficients,
                    family=family)
maxlam 

# lambdas <- 10^seq(round(log10(lambda_preset))-1, log10(maxlam)+1, by=0.05) #
# length(lambdas) # 152 - maybe in benchmark we should use 20 values like for the other methods?
lambdas <- 10^seq(round(log10(lambda_preset))-2, log10(maxlam), length.out=100) #
length(lambdas) # 20 - like for other methods
coefmat <- Matrix(0, nrow=length(lambdas), ncol=ncol(sim$X))
warmstart <- 0.1 # % by which L0glm.out_preset will be used as start, this value gave best solution quality
gamma <- 2
system.time(metrics <- data.frame(do.call(rbind, lapply(1:length(lambdas), function (i) { 
  lambda <- lambdas[[i]]
  
  # gamma = 2
  start = L0glm.out_preset_coefs # *scaling # TODO incorporate scaling if normalize=T in L0glm.fit.internal
  # start = warmstart*start + (1-warmstart)*rep(sd(start)^(1/gamma), ncol(X))
  
  start <- warmstart*start+(1-warmstart)*start_coefs(X, family)*100
  # start[L0glm.out_preset_coefs==0] <- 0
  L0glm.out <- L0glm.fit(X, y, 
                             family = gaussian(identity), 
                             lambda = lambda, 
                             nonnegative = TRUE,
                             dispersion = L0glm.out_preset$dispersion,
                             start = start, # L0glm.out_preset_coefs, # start, # rep(0.1, ncol(sim$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                             # warm.start = warmstart, 
                             solver = "eigen", # note: solver osqp better than eigen for very small lambda values))
                             method = 7L,
                             relaxed = FALSE,
                             thresh = 1E-5,
)
  coefmat[i,] <<- as.vector(L0glm.out$coefficients)
  return( unlist(c(lambda=lambda, get.metrics(beta = L0glm.out$coefficients, 
                                              X = X, y = y, 
                                              family = gaussian(identity), 
                                              beta_true = sim$beta_true))) )
} )))) # 0.87s for n=p=500 and using Eigen Cholesky solver nr 7
matplot(x=log10(lambdas), y=as.matrix(coefmat), type="l", lty=1,
        xlab="Log10(lambda)", ylab="Coefficient estimates", 
        main="Regularisation path L0 Gaussian GLM", col=1:6)
abline(v=log10(maxlam), col="red") # max lambda where all coefs would be penalized to 0
library(scales)
abline(h=as.vector(sim$beta_true[sim$beta_true!=0]), col=alpha(1:6,0.2))

# coefficients for solution that optimises given criterion & metrics
best = which.min(metrics[,crit])
L0glm.out_best_coefs = coefmat[best,]
plot(x=L0glm.out_best_coefs, y=sim$beta_true, pch=16, col='blue',
     main="L0 Gaussian GLM with tuned lambda",
     xlab="Estimated coefficients", ylab="True coefficients")


do.call(cbind, get.metrics(beta = L0glm.out_best_coefs,  # some metrics
                           X = X, y = y, 
                           family = family, 
                           beta_true = sim$beta_true))
#        p   n m deviance deviance_true  loglik    MD_coefs MAD_coefs  sMAPE_coefs bias_coefs
# [1,] 501 500 1 10250.99      1029.794 -1014.6 -0.03125375 0.2979296 0.0002203519   1.000016
# R2_coefs support_size_true support_size_estimated estimated_over_true_support_size TP FP  TN
# [1,] 0.9999993                50                     50                                1 50  0 450
#      FN precision recall specificity accuracy F1_score      gic     ebic      bic      aic     mbic
# [1,]  0         1      1           1        1        1 2608.413 2663.293 2346.145 2131.199 2838.836
plot(log10(metrics$lambda), metrics$deviance, type="l")
metrics$lambda[which.min(metrics$deviance)] # 1 # NOTE this is deviance relative to observed signal, not deviance relative to ground truth though
plot(log10(metrics$lambda), metrics$deviance_true, type="l")
metrics$lambda[which(metrics$deviance_true==min(metrics$deviance_true))] # 100 # deviance relative to ground truth
plot(log10(metrics$lambda), metrics$loglik, type="l")
abline(v=log10(maxlam), col="red") # max lambda where all coefs would be penalized to 0
metrics$lambda[which.max(metrics$loglik)] # 1.985 - this is AIC optimal value
plot(log10(metrics$lambda), metrics$MD_coefs, type="l")
plot(log10(metrics$lambda), metrics$MAD_coefs, type="l")
metrics$lambda[which.min(metrics$MAD_coefs)] # 4.466836
plot(log10(metrics$lambda), metrics$sMAPE_coefs, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(metrics$sMAPE_coefs)] # 4.466836
plot(log10(metrics$lambda), metrics$bias_coefs, type="l")
metrics$lambda[which.min(abs(metrics$bias_coefs-1))] # 1.778279
plot(log10(metrics$lambda), metrics$R2_coefs, type="l")
plot(log10(metrics$lambda), metrics$estimated_over_true_support_size, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.min(abs(metrics$estimated_over_true_support_size-1))] # 4.466836
plot(log10(metrics$lambda), metrics$TP, type="l")
plot(log10(metrics$lambda), metrics$FP, type="l")
plot(log10(metrics$lambda), metrics$TN, type="l")
plot(log10(metrics$lambda), metrics$FN, type="l")
plot(log10(metrics$lambda), metrics$precision, type="l")
plot(log10(metrics$lambda), metrics$recall, type="l")
plot(log10(metrics$lambda), metrics$specificity, type="l")
plot(log10(metrics$lambda), metrics$accuracy, type="l") # good one, optimal for lambda between 1 & 10
metrics$lambda[which.max(metrics$accuracy)] # 4.466836
plot(log10(metrics$lambda), metrics$F1_score, type="l")
metrics$lambda[which.max(metrics$F1_score)] # 
plot(log10(metrics$lambda), metrics$aic, type="l")
metrics$lambda[which.min(metrics$aic)] # 3.04
plot(log10(metrics$lambda), metrics$bic, type="l")
metrics$lambda[which.min(metrics$bic)] # 4.66
plot(log10(metrics$lambda), metrics$ebic, type="l")
metrics$lambda[which.min(metrics$ebic)] # 4.66
plot(log10(metrics$lambda), metrics$gic, type="l")
metrics$lambda[which.min(metrics$gic)] # 4.66
abline(v=log10(metrics$lambda[which.min(metrics$gic)]), col="blue") 
plot(log10(metrics$lambda), metrics$mbic, type="l")
metrics$lambda[which(metrics$mbic==min(metrics$mbic))] # 4.66


### 1.2.2 Comparison of solution quality to that returned by ordinis (adaptive lasso, MCP, SCAD), glmnet (lasso), abess & L0Learn ####

#### nonnegative lasso regression in glmnet ####
system.time(fit_glmnet <- glmnet(x = X, 
                                 y = y, 
                                 weights = W,
                                 family = "gaussian", # poisson(identity) should be supported but gave an error
                                 nlambda = 20L,
                                 intercept = F, 
                                 standardize = F, 
                                 # thresh = 1E-10,
                                 # maxit=100000,
                                 alpha = 1, # lasso (lower value, ie elastic net, made solution worse)
                                 lower.limits = 0)) # for nonnegativity constraints
# 0.01s for n=p=100, 0.01s for n=p=500, 0.05s for n=p=10000
plot(fit_glmnet, xvar='lambda')
abline(h=as.vector(sim$beta_true[sim$beta_true!=0]), col=alpha(1:6,0.2))


LLs <- logLik(fit_glmnet, weights=W, X=X, y=y)
edfs = attr(LLs, "df")
IC = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,crit]
which.min(IC) # 

bestlambda_lasso_glmnet = fit_glmnet$lambda[which.min(IC)]
coefs_lasso_glmnet = coef(fit_glmnet, s = bestlambda_lasso_glmnet)

sum(coefs_lasso_glmnet!=0) # not sparse enough
plot(x=coefs_lasso_glmnet[-1], y=sim$beta_true, pch=16, col='blue',
     main="glmnet nonnegative LASSO",
     xlab="Estimated coefficients", ylab="True coefficients")


#### adaptive LASSO using glmnet ####

# since glmnet is so much faster than ordinis & supports sparse
# covariate matrices, let's use that to fit an adaptive LASSO model, cf
# https://rpubs.com/kaz_yos/alasso
# note that for the first step we use LASSO rather than ridge coefficients, i.e. we use
# coefs_lasso_glmnet and then in a second step use a glmnet LASSO fit with
# adaptive penalty weights penalty.factor = 1 / (abs(coefs_lasso_glmnet))+delta)
# with delta=1E-5
# NOTE: If this would be iterated further this would be a way to approximate L0 norm
# penalized regression (cf. Nuel 2017 talk)

system.time({
  # we start with a regular LASSO fit, with lambda tuned according to our chosen IC
  fit_lasso_glmnet <- glmnet(x = X, 
                             y = y, 
                             weights = W,
                             family = "gaussian", # poisson(identity) should be supported but gave an error
                             nlambda = 20L,
                             intercept = F, 
                             standardize = F, 
                             # thresh = 1E-10,
                             # maxit=100000,
                             alpha = 1, # lasso (lower value, ie elastic net, made solution worse)
                             lower.limits = 0) # for nonnegativity constraints
  LLs <- logLik(fit_lasso_glmnet, weights=W, X=X, y=y)
  edfs = attr(LLs, "df")
  IC = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,crit]
  bestlambda_lasso_glmnet = fit_lasso_glmnet$lambda[which.min(IC)]
  coefs_lasso_glmnet = coef(fit_lasso_glmnet, s = bestlambda_lasso_glmnet)
  
  # in second step we do another LASSO fit with adaptive penalty weights
  delta = 1E-5
  fit_alasso_glmnet <- glmnet(x = X, 
                              y = y, 
                              weights = W,
                              family = "gaussian", # poisson(identity) should be supported but gave an error
                              nlambda = 20L,
                              penalty.factor = 1 / (abs(coefs_lasso_glmnet[-1])+delta), # adaptive penalty weights
                              intercept = F, 
                              standardize = F, 
                              # thresh = 1E-10,
                              # maxit=100000,
                              alpha = 1, # lasso (lower value, ie elastic net, made solution worse)
                              lower.limits = 0) # for nonnegativity constraints
  LLs <- logLik(fit_alasso_glmnet, weights=W, X=X, y=y)
  edfs = attr(LLs, "df")
  IC = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,crit]
  bestlambda_alasso_glmnet = fit_alasso_glmnet$lambda[which.min(IC)]
  coefs_alasso_glmnet = coef(fit_alasso_glmnet, s = bestlambda_lasso_glmnet)
})
# 0.02s for n=p=500, 0.11s for n=p=10000

plot(fit_alasso_glmnet, xvar='lambda')
abline(h=as.vector(sim$beta_true[sim$beta_true!=0]), col=alpha(1:6,0.2))

sum(coefs_alasso_glmnet!=0) # pretty good
plot(x=coefs_alasso_glmnet[-1], y=sim$beta_true, pch=16, col='blue',
     main="glmnet nonnegative adaptive LASSO",
     xlab="Estimated coefficients", ylab="True coefficients")





#### nonnegative L0Learn best subset fit after glmnet nonnegative LASSO ####
selcoefs = which(coefs_lasso_glmnet[-1]!=0)
p_selcoefs = length(selcoefs)
system.time(fit_nnL0Learn_selcoefs <- L0Learn.fit(x=as.matrix(X[,selcoefs,drop=F]),  # note: no support for sparse matrices
                                                  y=as.vector(y),
                                                  loss = "SquaredError",
                                                  penalty="L0", algorithm="CD",
                                                  maxSuppSize=(p_selcoefs-1), 
                                                  intercept=FALSE, 
                                                  nLambda=100, # NOTE: 100 gave much better results than 20
                                                  # activeSetNum = 3, maxSwaps = 100,
                                                  autoLambda=TRUE
                                                  # lambdaGrid=list(10^seq(-5,5,length.out=nlam)/ # DECREASING sequence recommended, but INCREASING works better under nonnegativity, see https://cran.r-project.org/web/packages/L0Learn/vignettes/L0Learn-vignette.html
                                                  #                   (norm(sim$y*sqrt(1/(sim$y+1)),"2")^2)) 
))  # lambda of 1 is scaled by 1/(norm(weighted y,"2")) because L0Learn scales y by L2 norm of y
# 0.61s for n=p=500
plot(fit_nnL0Learn_selcoefs)


LLs = logLik(fit_nnL0Learn_selcoefs, weights=W, X[,selcoefs,drop=F], y)
edfs = attr(LLs, "df")
LLs = as.vector(LLs)
IC = get.IC(loglik = LLs, n = n, p = p_selcoefs, edf = edfs)[,crit]
coefs_nnL0Learn_selcoefs = coef(fit_nnL0Learn_selcoefs,
                                lambda = fit_nnL0Learn_selcoefs$lambda[[1]][which.min(IC)],
                                gamma = 0)
coefs_nnL0Learn_after_lasso = Matrix(0, p)
coefs_nnL0Learn_after_lasso[selcoefs] = coefs_nnL0Learn_selcoefs
sum(coefs_nnL0Learn_after_lasso!=0) # 

plot(x=coefs_nnL0Learn_after_lasso, y=sim$beta_true, pch=16, col='blue',
     main="nnL0Learn after nonnegative LASSO",
     xlab="Estimated coefficients", ylab="True coefficients")


#### # full nonnegative L0Learn best subset fit ####
# quite good, but slow for n=p=10000
system.time(fit_nnL0Learn <- L0Learn.fit(x=as.matrix(X),  # note: no support for sparse matrices
                                         y=as.vector(y),
                                         loss = "SquaredError",
                                         penalty="L0", algorithm="CD",
                                         maxSuppSize=round(p/2),
                                         intercept=FALSE,
                                         nLambda=100, # NOTE: to make timings fair we should probably use 100 or 20 lambda values everywhere?
                                         # activeSetNum = 3, maxSwaps = 100,
                                         autoLambda=TRUE
                                         # lambdaGrid=list(10^seq(-5,5,length.out=nlam)/ # DECREASING sequence recommended, but INCREASING works better under nonnegativity, see https://cran.r-project.org/web/packages/L0Learn/vignettes/L0Learn-vignette.html
                                         #                   (norm(sim$y*sqrt(1/(sim$y+1)),"2")^2))
))  # lambda of 1 is scaled by 1/(norm(weighted y,"2")) because L0Learn scales y by L2 norm of y
# 1.13s for n=p=500
plot(fit_nnL0Learn)


LLs = logLik(fit_nnL0Learn, weights=W, X, y)
edfs = attr(LLs, "df")
LLs = as.vector(LLs)
IC = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,crit]
coefs_nnL0Learn = coef(fit_nnL0Learn,
                       lambda = fit_nnL0Learn$lambda[[1]][which.min(IC)],
                       gamma = 0)
sum(coefs_nnL0Learn!=0) # not sparse enough
plot(x=coefs_nnL0Learn, y=sim$beta_true, pch=16, col='blue',
     main="nnL0Learn least squares best subset",
     xlab="Estimated coefficients", ylab="True coefficients")




#### abess best subset fit after nonnegative LASSO ####
# very good solution, and fast for small problems, but
# becomes very slow for large problems, e.g. n=p=10000
selcoefs = which(coefs_lasso_glmnet[-1]!=0)
p_selcoefs = length(selcoefs)
system.time(fit_abessafterlasso_selcoefs <- abess(x = X[,selcoefs,drop=F], 
                                                  y = y, 
                                                  # weights = W,  # approx 1/variance weights to approximate family="poisson", which is not supported
                                                  family = "gaussian",
                                                  tune.path = "sequence", 
                                                  support.size = c(1:(p_selcoefs-1)),
                                                  lambda = 0, 
                                                  warm.start = TRUE, # FALSE worked better here, but much slower
                                                  tune.type = crit))  # options bic, ebic, aic or cv
# 0.08s for n=p=500 
# note: abess currently has no support for nonnegativity constraints, so 
# abess after nonnegative LASSO works much better
plot(fit_abessafterlasso_selcoefs, type="tune")  # optimal support size 
coefs_abessafterlasso_selcoefs = as.matrix(extract(fit_abessafterlasso_selcoefs)$beta) # coefficient matrix for best subset
coefs_abessafterlasso = Matrix(rep(0, p), sparse=T)
coefs_abessafterlasso[selcoefs] = coefs_abessafterlasso_selcoefs
coefs_abessafterlasso[coefs_abessafterlasso<0] = 0
sum(coefs_abessafterlasso!=0) # pretty good
plot(x=coefs_abessafterlasso, y=sim$beta_true, pch=16, col='blue',
     main="abess best subset weighted least squares\nafter nonnegative LASSO",
     xlab="Estimated coefficients", ylab="True coefficients")



#### # full abess best subset fit (with EBIC criterion) ####
# # note: currently no support for nonnegativity constraints
# # solution not great, so marked out for now & slow for n=p=10000
# system.time(fit_abess <- abess(x = X, 
#                                y = y, 
#                                weights = W,  # approx 1/variance weights to approximate family="poisson", which is not supported
#                                family = "gaussian",
#                                tune.path = "sequence", 
#                                support.size = c(1:round(p*0.5)),
#                                lambda=0, 
#                                warm.start = TRUE, # TRUE worked best here
#                                tune.type = "ebic")) # had to use ebic here for optimal performance, XXs for n=p=XX, options bic, ebic, aic or cv
# # fast for small problems, but becomes slow for n=p=10000
# plot(fit_abess, type="tune")  # optimal support size 
# coefs_abess = as.matrix(extract(fit_abess)$beta) # coefficient matrix for best subset
# coefs_abess[coefs_abess<0] = 0
# sum(coefs_abess!=0) # 11, quite OK
# plot(x=coefs_abess, y=sim$beta_true, pch=16, col='blue',
#      main="abess best subset weighted least squares",
#      xlab="Estimated coefficients", ylab="True coefficients")



#### nonnegative weighted MCP penalized regression in ordinis ####
# NOTE: ncvreg also has MCP penalized models, but not with nonnegativity constraints
system.time(fit_MCP <- ordinis(X, y, 
                               # weights=W, # approx 1/variance Poisson weights
                               nlambda = 100L, # 1000L gives much better solution than 100L, 20L doesn't work & coefficient path looks very different
                               penalty = "mcp", # using MCP penalty
                               intercept = FALSE,
                               lower.limits = rep(0, p), # force all coefficients to be positive
                               standardize = FALSE,
                               # maxit = 1000,
                               # tol = 1E-100,
                               # maxit.irls = 1000,
                               # tol.irls = 1E-100,
                               # gamma = 100, # ifelse(penalty == "scad", 3.7, 1.4),
                               alpha = 1))  # ie no elastic net
# 0.09s for n=p=500, XXs for n=p=10000
plot(fit_MCP)
abline(h=as.vector(sim$beta_true[sim$beta_true!=0]), col=alpha(1:6,0.2)) # note : ordinis uses different colours, so cols don't match now


LLs = logLik2(fit_MCP, weights=W, X, y)
edfs = attr(LLs, "df")
LLs = as.vector(LLs)
IC = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,crit]
coefs_MCP_ordinis = fit_MCP$beta[,which.min(IC)]
sum(coefs_MCP_ordinis!=0) # not sparse enough
plot(x=coefs_MCP_ordinis[-1], y=sim$beta_true, pch=16, col='blue',
     main="ordinis nonnegative MCP",
     xlab="Estimated coefficients", ylab="True coefficients")


#### nonnegative weighted SCAD penalized regression in ordinis ####
# NOTE: ncvreg also has SCAD penalized models, but not with nonnegativity constraints
system.time(fit_SCAD <- ordinis(X, 
                                y, 
                                # weights=W, # approx 1/variance Poisson weights
                                nlambda = 100L, # 1000L gives much better solution than 100L or 20L & coefficient path looks very different
                                penalty = "scad", # using SCAD penalty
                                intercept = FALSE,
                                lower.limits = rep(0, p), # force all coefficients to be positive
                                standardize = FALSE,
                                # maxit = 1000,
                                # tol = 1E-100,
                                # maxit.irls = 1000,
                                # tol.irls = 1E-100,
                                # gamma = 100, # ifelse(penalty == "scad", 3.7, 1.4),
                                alpha = 1))  # ie no elastic net
# 0.14s for n=p=500
plot(fit_SCAD)
abline(h=as.vector(sim$beta_true[sim$beta_true!=0]), col=alpha(1:6,0.2)) # note : ordinis uses different colours, so cols don't match now
LLs = logLik2(fit_SCAD, weights=W, X, y)
edfs = attr(LLs, "df")
LLs = as.vector(LLs)
IC = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,crit]
coefs_SCAD_ordinis = fit_SCAD$beta[,which.min(IC)]
sum(coefs_SCAD_ordinis!=0) # not sparse enough
plot(x=coefs_SCAD_ordinis[-1], y=sim$beta_true, pch=16, col='blue',
     main="ordinis nonnegative SCAD",
     xlab="Estimated coefficients", ylab="True coefficients")


#### # nonnegative weighted lasso regression in ordinis ####
# system.time(fit_lasso_ordinis <- ordinis(X, 
#                                          y, 
#                                  weights=W, # approx 1/variance Poisson weights
#                                  nlambda = 20L, # 1000L gives much better solution than 100L or 20L & coefficient path looks very different
#                                  lambda.min.ratio = 0.0001,
#                                  penalty = "lasso", # using LASSO penalty
#                                  intercept = FALSE,
#                                  lower.limits = rep(0, p), # force all coefficients to be positive
#                                  standardize = FALSE,
#                                  maxit = 1000,
#                                  # tol = 1E-100,
#                                  # maxit.irls = 1000,
#                                  # tol.irls = 1E-100,
#                                  # gamma = 100, # ifelse(penalty == "scad", 3.7, 1.4),
#                                  alpha = 1))  # ie no elastic net
# # XXs for n=p=100, 522s for n=p=10000
# plot(fit_lasso_ordinis)
# abline(h=as.vector(sim$beta_true[sim$beta_true!=0]), col=alpha(1:6,0.2)) # note : ordinis uses different colours, so cols don't match now

# LLs = logLik2(object=fit_lasso_ordinis, weights=W, X=X, y=y)
# edfs = attr(LLs, "df")
# LLs = as.vector(LLs)
# IC = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,crit]
# which.min(IC) # 
# bestlambda_lasso_ordinis = fit_lasso_ordinis$lambda[which.min(IC)]
# bestlambda_lasso_ordinis # 
# coefs_lasso_ordinis = fit_lasso_ordinis$beta[,which.min(IC)]
# sum(coefs_lasso_ordinis!=0) # not sparse enough
# plot(x=coefs_lasso_ordinis[-1], y=sim$beta_true, pch=16, col='blue',
#      main="ordinis nonnegative weighted LASSO",
#      xlab="Estimated coefficients", ylab="True coefficients")


#### # nonnegative weighted adaptive lasso regression in ordinis ####
# system.time(fit_alasso_ordinis <- ordinis(X, y, 
#                                   weights=W, # approx 1/variance Poisson weights
#                                   nlambda = 20L, # 1000L gives much better solution than 100L or 20L & coefficient path looks very different
#                                   penalty = "alasso", # using adaptive LASSO penalty
#                                   intercept = FALSE,
#                                   lower.limits = rep(0, p), # force all coefficients to be positive
#                                   standardize = FALSE,
#                                   # maxit = 1000,
#                                   # tol = 1E-100,
#                                   # maxit.irls = 1000,
#                                   # tol.irls = 1E-100,
#                                   # gamma = 100, # ifelse(penalty == "scad", 3.7, 1.4),
#                                   alpha = 1))  # ie no elastic net
# # XXs for n=p=100, XXs for n=p=10000
# plot(fit_alasso_ordinis)
# abline(h=as.vector(sim$beta_true[sim$beta_true!=0]), col=alpha(1:6,0.2)) # note : ordinis uses different colours, so cols don't match now
# 
# LLs = logLik2(fit_alasso_ordinis, weights=W, X, y)
# edfs = attr(LLs, "df")
# LLs = as.vector(LLs)
# IC = get.IC(loglik = LLs, n = n, p = p, edf = edfs)[,crit]
# coefs_alasso_ordinis = fit_alasso_ordinis$beta[,which.min(IC)]
# sum(coefs_alasso_ordinis!=0) # 9 - slightly too sparse
# plot(x=coefs_alasso_ordinis[-1], y=sim$beta_true, pch=16, col='blue',
#      main="ordinis nonnegative weighted adaptive LASSO",
#      xlab="Estimated coefficients", ylab="True coefficients")



#### SUMMARY TABLE WITH BENCHMARKS ####
benchmarks_1_2 = rbind(
  data.frame(method="L0glm_tuned_lambda", do.call(cbind, get.metrics(beta = L0glm.out_best_coefs,  # some metrics
                                                                     X = X, y = y, 
                                                                     family = family, 
                                                                     beta_true = sim$beta_true))),
  
  data.frame(method="L0glm_preset_lambda", do.call(cbind, get.metrics(beta = L0glm.out_preset_coefs,  # some metrics
                                                                      X = X, y = y, 
                                                                      family = family, 
                                                                      beta_true = sim$beta_true))),
  data.frame(method="nnalasso_glmnet", do.call(cbind, get.metrics(beta = coefs_alasso_glmnet,  # some metrics
                                                                  X = cbind(1, X), y = y, 
                                                                  weights = W,
                                                                  family = gaussian(identity), 
                                                                  beta_true = c(0, as.vector(sim$beta_true))) )),
  data.frame(method="nnL0Learn_after_glmnet_nnlasso", do.call(cbind, get.metrics(beta = coefs_nnL0Learn_after_lasso,  # some metrics
                                                                                 X = X, y = y, 
                                                                                 weights = W,
                                                                                 family = gaussian(identity),
                                                                                 beta_true = as.vector(sim$beta_true)) )),
  data.frame(method="abess_after_glmnet_nnlasso", do.call(cbind, get.metrics(beta = coefs_abessafterlasso,  # some metrics
                                                                             X = X, y = y, 
                                                                             weights = W,
                                                                             family = gaussian(identity),
                                                                             beta_true = as.vector(sim$beta_true)) )),
  data.frame(method="nnL0Learn", do.call(cbind, get.metrics(beta = coefs_nnL0Learn,  # some metrics
                                                            X = X, y = y, 
                                                            weights = W,
                                                            family = gaussian(identity), 
                                                            beta_true = as.vector(sim$beta_true)) )),
  data.frame(method="nnlasso_glmnet", do.call(cbind, get.metrics(beta = coefs_lasso_glmnet,  # some metrics
                                                                 X = cbind(1, X), y = y, 
                                                                 weights = W,
                                                                 family = gaussian(identity), 
                                                                 beta_true = c(0, as.vector(sim$beta_true))) )),
  data.frame(method="nnMCP_ordinis", do.call(cbind, get.metrics(beta = coefs_MCP_ordinis,  # some metrics
                                                                X = cbind(1, X), y = y, 
                                                                weights = W,
                                                                family = gaussian(identity),
                                                                beta_true = c(0, as.vector(sim$beta_true))) )),
  data.frame(method="nnSCAD_ordinis", do.call(cbind, get.metrics(beta = coefs_SCAD_ordinis,  # some metrics
                                                                 X = cbind(1, X), y = y, 
                                                                 weights = W,
                                                                 family = gaussian(identity), 
                                                                 beta_true = c(0, as.vector(sim$beta_true))) ))
  # data.frame(method="nnlasso_ordinis", do.call(cbind, get.metrics(beta = coefs_lasso_ordinis,  # some metrics
  #                                                                 X = cbind(1, X), y = y, 
  #                                                                 weights = W,
  #                                                                 family = gaussian(identity),
  #                                                                 beta_true = c(0, as.vector(sim$beta_true))) )),
  # data.frame(method="nnalasso_ordinis", do.call(cbind, get.metrics(beta = coefs_alasso_ordinis,  # some metrics
  #                                                       X = cbind(1, X), y = y, 
  #                                                       weights = W,
  #                                                       family = gaussian(identity), 
  #                                                       beta_true = c(0, as.vector(sim$beta_true))) )),
  # data.frame(method="abess", do.call(cbind, get.metrics(beta = coefs_abess,  # some metrics
  #                                                                X = X, y = y, 
  #                                                                weights = W,
  #                                                                family = gaussian(identity), 
  #                                                                beta_true = as.vector(sim$beta_true)) ))
)
benchmarks_1_2 = benchmarks_1_2[order(benchmarks_1_2[,crit]),] # sort by IC
benchmarks_1_2
# CONCLUSION: L0glm with tuned lambda best followed by L0glm with preset lambda,
# followed by L0Learn after lasso; 
# all other methods much worse

# For n=p=500
#                           method   p   n m     deviance deviance_true    loglik     MD_coefs
# 1             L0glm_tuned_lambda 501 500 1     10250.99  1.029794e+03 -1014.600  -0.03125375
# 2            L0glm_preset_lambda 501 500 1     10189.95  1.390415e+03 -1017.107  -0.02009234
# 4 nnL0Learn_after_glmnet_nnlasso 501 500 1    128245.96  1.212676e+05 -1664.244  -0.30244913
# 9                 nnSCAD_ordinis 502 500 1    334347.66  3.305220e+05 -1892.801   0.25641215
# 8                  nnMCP_ordinis 502 500 1   2131544.51  2.114008e+06 -2362.907   2.27259794
# 6                      nnL0Learn 501 500 1   2125871.34  2.115412e+06 -2372.240   2.36233714
# 3                nnalasso_glmnet 502 500 1   3712409.86  3.720694e+06 -2483.615  -5.47144854
# 7                 nnlasso_glmnet 502 500 1   4078061.04  4.098858e+06 -2566.100 -11.28178925
# 5     abess_after_glmnet_nnlasso 501 500 1 335541753.31  3.356881e+08 -3620.631  31.01370750
#     MAD_coefs  sMAPE_coefs bias_coefs  R2_coefs support_size_true support_size_estimated
# 1   0.2979296 0.0002203519  1.0000161 0.9999993                50                     50
# 2   6.4078660 0.0324525129  1.0019211 0.9986880                50                     54
# 4  88.5999833 0.2076315006  0.9911148 0.9470501                50                     68
# 9  82.4210909 0.1628268383  0.9596194 0.9259290                50                     58
# 8 239.6124410 0.2738004153  0.8680802 0.6620018                50                     65
# 6 275.4875187 0.3351748124  0.8355767 0.6399233                50                     74
# 3 101.8715539 0.1389233611  0.9769416 0.9138261                50                     47
# 7 117.1300908 0.3250827351  1.0364852 0.9213728                50                    106
# 5 156.2071531 0.1567333836  0.8144065 0.7004816                50                     57
#   estimated_over_true_support_size TP FP  TN FN precision recall specificity  accuracy  F1_score
# 1                             1.00 50  0 450  0 1.0000000   1.00   1.0000000 1.0000000 1.0000000
# 2                             1.08 48  6 444  2 0.8888889   0.96   0.9866667 0.9840000 0.9230769
# 4                             1.36 34 34 416 16 0.5000000   0.68   0.9244444 0.9000000 0.5762712
# 9                             1.16 34 24 427 16 0.5862069   0.68   0.9467849 0.9201597 0.6296296
# 8                             1.30 24 41 410 26 0.3692308   0.48   0.9090909 0.8662675 0.4173913
# 6                             1.48 21 53 397 29 0.2837838   0.42   0.8822222 0.8360000 0.3387097
# 3                             0.94 33 14 437 17 0.7021277   0.66   0.9689579 0.9381238 0.6804124
# 7                             2.12 43 63 388  7 0.4056604   0.86   0.8603104 0.8602794 0.5512821
# 5                             1.14 35 22 428 15 0.6140351   0.70   0.9511111 0.9260000 0.6542056
#        gic     ebic      bic      aic     mbic
# 1 2608.413 2663.293 2346.145 2131.199 2838.836
# 2 2658.855 2718.040 2376.017 2144.213 2907.351
# 4 4112.129 4186.379 3757.295 3466.487 4423.878
# 9 4455.889 4519.398 4152.265 3903.603 4722.477
# 8 5475.624 5546.668 5135.977 4857.813 5773.842
# 6 5596.266 5676.972 5210.576 4894.481 5935.123
# 3 5512.548 5564.216 5265.532 5063.230 5729.433
# 7 6347.804 6462.981 5797.164 5346.201 6831.277
# 5 7899.977 7962.389 7601.710 7357.263 8162.026


# 2. MULTICHANNEL DECONVOLUTION ####

n <- p <- 200 # 500 # 10000 # to check performance use e.g. n=p=10000, otherwise use e.g. n=p=500
s <- 0.1 # sparsity, nr of timepoints with a peak 
k <- round(n*s)
mean_beta <- 10000
sd_logbeta <- 1
sd_noise <- 5 # SD of noise for Gaussian noise case

# lambda expected to optimise given IC
crit <- "bic" # choose lambda based on what would be expected to optimise this IC
# options in order of increasing regularisation: aic, bic, gic, ebic, mbic
lambda_preset = preset.lambda(IC = crit, p=p, n=n) # options aic, bic, gic, ebic, mbic  
lambda_preset 


## 2.1 Multichannel deconvolution assuming Poisson noise ####
### 2.1.1 L0glm solution ####
set.seed(1)
sim <- simulate_spike_train(n = n, 
                            p = p,
                            k = ifelse(round(n/10)<=100, round(n/10), 100),
                            family="poisson", 
                            mean_beta = mean_beta,
                            sd_logbeta = sd_logbeta,
                            multichannel=TRUE)
X <- sim$X
y <- sim$y
family <- poisson(identity)
m <- ncol(y) # nr of tasks/channels

# relaxed L0 group fit, with preset lambda
system.time(L0glm.out_preset <- L0glm.fit(X, y, 
                               family = family, 
                               lambda = lambda_preset, 
                               grouped = TRUE,
                               relaxed = TRUE, # NOTE: CHECK setting relaxed to TRUE does not seem to improve solution quality
                               nonnegative = TRUE,
                               parallelize = TRUE,
                               maxit = 10, # default 25
                               group.maxit = 4, # default 6
                               dispersion = 1,
                               start = rep(1, ncol(X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                               solver = "eigen", # note: solver osqp better than eigen for very small lambda values))
                               method = 7L,
                               thresh = 1E-5
)) # 10.46s for n=p=100, 10s for n=p=500, 30s for n=p=10000

# TODO: get rid of warnings
# Warning messages:
# 1: In (function (..., deparse.level = 1)  :
#        number of rows of result is not a multiple of vector length (arg 8)
       
L0glm.out_preset_coefs = L0glm.out_preset$coefficients
abline(v=which(rowMin(L0glm.out_preset_coefs)!=0), col = 'cyan')
graph2pdf(file="..//benchmarks/plots/2_1 multichannel deconvolution Poisson.pdf", width=7, height=5)

# some metrics
metrics_L0glm = do.call(cbind, get.metrics(beta = L0glm.out_preset_coefs, 
                           X = X, y = y, 
                           family = poisson(identity), 
                           beta_true = sim$beta_true))
metrics_L0glm
#        p   n   m deviance deviance_true    loglik     MD_coefs MAD_coefs sMAPE_coefs bias_coefs
# [1,] 500 500 159 37570.71      4.48e+12 -105213.2 -0.005010161 0.2475819  0.02308121   1.000824
#       R2_coefs MD_selcoefs MAD_selcoefs sMAPE_selcoefs bias_selcoefs R2_selcoefs support_size_true
# [1,] 0.9998294 0.007370507      4.77074      0.3679277      1.000824   0.9998312                50
#      support_size_estimated estimated_over_true_support_size TP FP  TN FN precision recall
# [1,]                     50                                1 50  0 450  0         1      1
#      specificity accuracy F1_score    gic     ebic      bic      aic     mbic
# [1,]           1        1        1 256249 260590.7 235508.5 218498.4 274482.7

# coefficients & metrics over a range of lambda values
# calculate max lambda where all coefs would be penalized to 0
maxlam <- maxlambda(X=X, weights=L0glm.out_preset$weights, 
                    z=y, dispersion=1, no.pen=NULL, coefs=L0glm.out_preset$coefficients,
                    family=family)
maxlam # 

lambdas <- 10^seq(round(log10(lambda_preset))-1, log10(maxlam)+1, by=1) #
coefmat <- Matrix(0, nrow=length(lambdas), ncol=ncol(sim$X)) # matrix with avg coefficients over tasks/channels
system.time(metrics <- data.frame(do.call(rbind, lapply(1:length(lambdas), function (i) { 
  lambda <- lambdas[[i]] 
  L0glm.out <- L0glm.fit(X, y, 
                             family = family, 
                             lambda = lambda, 
                             grouped = TRUE,
                             relaxed = TRUE, 
                             nonnegative = TRUE,
                             dispersion = L0glm.out_preset$dispersion,
                             start = rep(1, ncol(sim$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                             solver = "eigen", # note: solver osqp better than eigen for very small lambda values))
                             method = 7L,
                             parallelize = TRUE,
                             maxit = 10, # default 25
                             group.maxit = 4, # default 6
                             thresh = 1E-5
  )$coefficients
  coefmat[i,] <<- rowMeans(L0glm.out) # avg coefficients across tasks/channels
  return( unlist(c(lambda=lambda, get.metrics(beta = L0glm.out, 
                                              X = X, y = y, 
                                              family = poisson(identity), 
                                              beta_true = sim$beta_true) )) )
} ) ))) # 86s for n=p=100 unparallelized, 64s parallelized
# TODO fix bug
# Error in `[<-`(`*tmp*`, variable.indices, , value = Matrix(do.call(cbind,  : 
# subscript out of bounds
# In addition: There were 36 warnings (use warnings() to see them)
# GET rid of warning
# <sparse>[ <logic> ]: .M.sub.i.logical() maybe inefficient

matplot(x=log10(metrics$lambda), y=as.matrix(coefmat), type="l", lty=1,
        xlab="Log10(lambda)", ylab="Coefficient estimates", main="Regularisation path Group L0 Poisson GLM")
# abline(v=log10(maxlam), col="red") # max lambda where all coefs would be penalized to 0

plot(log10(metrics$lambda), metrics$deviance, type="l")
metrics$lambda[which.min(metrics$deviance)] # 0.1 # NOTE this is deviance relative to observed signal, not deviance relative to ground truth though
plot(log10(metrics$lambda), metrics$deviance_true, type="l")
metrics$lambda[which(metrics$deviance_true==min(metrics$deviance_true))] # 0.2511886 0.3981072 1.0000000 # deviance relative to ground truth
plot(log10(metrics$lambda), metrics$loglik, type="l")
metrics$lambda[which.max(metrics$loglik)] # 0.6309573
plot(log10(metrics$lambda), metrics$MD_coefs, type="l")
plot(log10(metrics$lambda), metrics$MAD_coefs, type="l")
metrics$lambda[which.min(metrics$MAD_coefs)] # 1.584893
plot(log10(metrics$lambda), metrics$sMAPE_coefs, type="l") # good one
metrics$lambda[which.min(metrics$sMAPE_coefs)] # 1.584893
plot(log10(metrics$lambda), metrics$bias_coefs, type="l")
metrics$lambda[which.min(abs(metrics$bias_coefs-1))] # 10
plot(log10(metrics$lambda), metrics$R2_coefs, type="l")
plot(log10(metrics$lambda), metrics$estimated_over_true_support_size, type="l") # good one
range(metrics$lambda[which(abs(metrics$estimated_over_true_support_size-1)==min(abs(metrics$estimated_over_true_support_size-1)))]) # 1.584893 10.000000
plot(log10(metrics$lambda), metrics$TP, type="l")
plot(log10(metrics$lambda), metrics$FP, type="l")
plot(log10(metrics$lambda), metrics$TN, type="l")
plot(log10(metrics$lambda), metrics$FN, type="l")
plot(log10(metrics$lambda), metrics$precision, type="l")
plot(log10(metrics$lambda), metrics$recall, type="l")
plot(log10(metrics$lambda), metrics$specificity, type="l")
plot(log10(metrics$lambda), metrics$accuracy, type="l") # good one
metrics$lambda[which.max(metrics$accuracy)] # 1.584893
plot(log10(metrics$lambda), metrics$F1_score, type="l")
metrics$lambda[which.max(metrics$F1_score)] # 1.584893
plot(log10(metrics$lambda), metrics$aic, type="l")
metrics$lambda[which(metrics$aic==min(metrics$aic))] # 0.6309573
plot(log10(metrics$lambda), metrics$bic, type="l")
metrics$lambda[which.min(metrics$bic)] # 0.6309573
plot(log10(metrics$lambda), metrics$ebic, type="l")
metrics$lambda[which.min(metrics$ebic)] # 0.6309573
plot(log10(metrics$lambda), metrics$gic, type="l")
metrics$lambda[which.min(metrics$gic)] # 0.6309573
abline(v=log10(metrics$lambda[which.min(metrics$gic)]), col="blue") 
plot(log10(metrics$lambda), metrics$mbic, type="l")
metrics$lambda[which(metrics$mbic==min(metrics$mbic))] # 0.6309573
# TODO fix NaNs or Infs



### 2.1.2 Comparison with abess weighted mgaussian best subset glm #### 

# This was our L0 glm metrics:
metrics_L0glm = do.call(cbind, get.metrics(beta = L0glm.out_preset_coefs, 
                                           X = X, y = y, 
                                           family = family, 
                                           beta_true = sim$beta_true))
metrics_L0glm
# p   n   m deviance deviance_true    loglik    MD_coefs MAD_coefs sMAPE_coefs bias_coefs
# [1,] 100 100 133 7624.287      1.75e+12 -20813.86 -0.01060908 0.3244337  0.02987475   0.999666
# R2_coefs support_size_true support_size_estimated estimated_over_true_support_size TP FP TN
# [1,] 0.9996625                10                     10                                1 10  0 90
# FN precision recall specificity accuracy F1_score      gic     ebic      bic      aic     mbic
# [1,]  0         1      1           1        1        1 47225.92 48959.14 45293.43 43219.71 50417.88

# abess mgaussian best subset fit plus metrics:
# best subset multitask learning using family="mgaussian" using abess
library(abess)
system.time(abess_fit <- abess(x = X, 
                               y = y, 
                               weights = 1/(y+0.1),  # approx 1/variance weights to approximate family="poisson", which is not supported
                               family = "mgaussian",
                               tune.path = "sequence", 
                               support.size = c(1:(p-1)),
                               lambda=0, 
                               warm.start = TRUE,
                               tune.type = "gic")) # 0.28s for n=p=100, 20.99s for n=p=500, or cv or ebic or bic or aic
# note: currently no support for nonnegativity constraints
plot(abess_fit, type="tune")  # optimal support size 
beta_abess = as.matrix(extract(abess_fit)$beta) # coefficient matrix for best subset
beta_abess[abs(beta_abess)<1E-5] = 0 # apply same thresholding as in L0glm
library(qlcMatrix)
image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=beta_abess^0.1, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="abess multitask weighted mgaussian (red=true peaks)")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim$X))[as.vector(rowMin(beta_abess)!=0)], col="cyan")
sum(rowMax(sim$beta_true)!=0) # 10 true peaks
sum(rowMin(beta_abess)!=0) # 10 peaks detected

# note: for larger p and lower sparsity, solution quality quickly becomes a lot worse
# see https://github.com/abess-team/abess/issues/510

metrics_abess = do.call(cbind, get.metrics(beta = beta_abess, 
                                           X = X, y = y, 
                                           family = poisson(identity), 
                                           beta_true = sim$beta_true))
metrics_abess
# TODO fix warnings
# Warning messages:
# 1: In log(y/mu) : NaNs produced
# 2: In log(y/mu) : NaNs produced
# 3: In dpois(y, mu, log = TRUE) : NaNs produced




## 2.2 Multichannel deconvolution assuming Gaussian noise ####

n <- p <- 200

### 2.2.1 L0glm solution ####
set.seed(1)
sim <- simulate_spike_train(n = n, 
                            p = p,
                            k = ifelse(round(n/10)<=100, round(n/10), 100),
                            family = "gaussian", 
                            mean_beta = mean_beta,
                            sd_logbeta = sd_logbeta,
                            sd_noise = sd_noise, 
                            multichannel = TRUE)
X <- sim$X
y <- sim$y
family <- gaussian(identity)
dispersion <- sd_noise^2
m <- ncol(y) # nr of tasks/channels

# fit with preset lambda
system.time(L0glm.out_preset <- L0glm.fit(X, y, 
                               family = family, 
                               lambda = "ebic", # lambda_preset, 
                               grouped = TRUE,
                               relaxed = TRUE,
                               nonnegative = TRUE,
                               dispersion = NULL,
                               # start = rep(1, ncol(sim$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                               solver = "eigen", # note: solver osqp better than eigen for very small lambda values))
                               method = 7L,
                               parallelize = TRUE, 
                               maxit = 20, # default 25
                               group.maxit = 4, # default 6
                               thresh = 1E-5
)) # 9.79s for n=p=100, 12s for n=p=500, 33s for n=p=10000
# TODO fix bug L0glm.out_preset$dispersion comes out as NULL

L0glm.out_preset_coefs = L0glm.out_preset$coefficients
abline(v=which(rowMin(L0glm.out_preset_coefs)!=0), col = 'cyan')
graph2pdf(file="..//benchmarks/plots/2_2 multichannel deconvolution Gaussian.pdf", width=7, height=5)

# some metrics
metrics_L0glm = do.call(cbind, get.metrics(beta = as.matrix(L0glm.out_preset_coefs), 
                           X = X, y = y, 
                           family = family, 
                           beta_true = sim$beta_true))
metrics_L0glm
#        p   n   m deviance deviance_true   loglik MD_coefs  MAD_coefs sMAPE_coefs bias_coefs  R2_coefs
# [1,] 501 500 161 48273.95      8403.259 -99995.1 0.018487 0.05023418  0.08479029   1.000029 0.9999976
#      support_size_true support_size_estimated estimated_over_true_support_size TP FP  TN FN precision recall
# [1,]                50                     50                                1 50  0 450  0         1      1
#      specificity accuracy F1_score    gic     ebic      bic      aic   mbic
# [1,]           1        1        1 277832 285207.4 242585.1 213698.2 308799

# coefficients & metrics over a range of lambda values
# calculate max lambda where all coefs would be penalized to 0
maxlam <- maxlambda(X=X, weights=L0glm.out_preset$weights, 
                    z=y, dispersion=1, # L0glm.out_preset$dispersion, 
                    no.pen=NULL, coefs=L0glm.out_preset$coefficients,
                    family=family)
maxlam # 

lambdas <- 10^seq(round(log10(lambda_preset))-1, log10(maxlam)+1, by=1) #
coefmat <- Matrix(0, nrow=length(lambdas), ncol=ncol(sim$X)) # matrix with avg coefficients over tasks/channels
system.time(metrics <- data.frame(do.call(rbind, lapply(1:length(lambdas), function (i) { 
  lambda <- lambdas[[i]]
  L0glm.out <- L0glm.fit(X, y, 
                             family = gaussian(identity), 
                             lambda = lambda, 
                             grouped = TRUE,
                             relaxed = TRUE, 
                             nonnegative = TRUE,
                             dispersion = 1,
                             start = rep(1, ncol(sim$X)), # coefficient starting values - high value in combination with nonzero centere ridge regression enforces nonnegativity even when using unconstrained eigen solver
                             solver = "eigen", # note: solver osqp better than eigen for very small lambda values))
                             maxit = 10, # default 25
                             group.maxit = 4, # default 6
                             method = 7L,
                             parallelize = TRUE, 
                             thresh = 1E-5
  )$coefficients
  coefmat[i,] <<- rowMeans(L0glm.out) # avg coefficients across tasks/channels
  return( unlist(c(lambda=lambda, get.metrics(beta = L0glm.out, 
                                              X = X, y = y, 
                                              family = family, 
                                              beta_true = sim$beta_true))) )
} )))) # 85s for n=p=100

matplot(x=log10(metrics$lambda), y=as.matrix(coefmat), type="l", lty=1,
        xlab="Log10(lambda)", ylab="Coefficient estimates", main="Regularisation path Group L0 Gaussian GLM")
# abline(v=log10(maxlam), col="red") # max lambda where all coefs would be penalized to 0


plot(log10(metrics$lambda), metrics$deviance, type="l")
metrics$lambda[which.min(metrics$deviance)] # 1# NOTE this is deviance relative to observed signal, not deviance relative to ground truth though
plot(log10(metrics$lambda), metrics$deviance_true, type="l")
range(metrics$lambda[which(metrics$deviance_true==min(metrics$deviance_true))]) # 1.0000 251.1886 # deviance relative to ground truth
plot(log10(metrics$lambda), metrics$loglik, type="l")
metrics$lambda[which.max(metrics$loglik)] # 1
plot(log10(metrics$lambda), metrics$MD_coefs, type="l")
plot(log10(metrics$lambda), metrics$MAD_coefs, type="l")
metrics$lambda[which.min(metrics$MAD_coefs)] # 1
plot(log10(metrics$lambda), metrics$sMAPE_coefs, type="l") # good one
metrics$lambda[which.min(metrics$sMAPE_coefs)] # 1
plot(log10(metrics$lambda), metrics$bias_coefs, type="l")
metrics$lambda[which.min(abs(metrics$bias_coefs-1))] # 1
plot(log10(metrics$lambda), metrics$R2_coefs, type="l")
plot(log10(metrics$lambda), metrics$estimated_over_true_support_size, type="l") # good one
range(metrics$lambda[which(abs(metrics$estimated_over_true_support_size-1)==min(abs(metrics$estimated_over_true_support_size-1)))]) # 1.000 1584.893
plot(log10(metrics$lambda), metrics$TP, type="l")
plot(log10(metrics$lambda), metrics$FP, type="l")
plot(log10(metrics$lambda), metrics$TN, type="l")
plot(log10(metrics$lambda), metrics$FN, type="l")
plot(log10(metrics$lambda), metrics$precision, type="l")
plot(log10(metrics$lambda), metrics$recall, type="l")
plot(log10(metrics$lambda), metrics$specificity, type="l")
plot(log10(metrics$lambda), metrics$accuracy, type="l") # good one
metrics$lambda[which.max(metrics$accuracy)] # 1
plot(log10(metrics$lambda), metrics$F1_score, type="l")
metrics$lambda[which.max(metrics$F1_score)] # 1
plot(log10(metrics$lambda), metrics$aic, type="l")
range(metrics$lambda[which(metrics$aic==min(metrics$aic))]) # 1.0000 251.1886
plot(log10(metrics$lambda), metrics$bic, type="l")
metrics$lambda[which.min(metrics$bic)] # 1
plot(log10(metrics$lambda), metrics$ebic, type="l")
metrics$lambda[which.min(metrics$ebic)] # 1
plot(log10(metrics$lambda), metrics$gic, type="l")
metrics$lambda[which.min(metrics$gic)] # 1
abline(v=log10(metrics$lambda[which.min(metrics$gic)]), col="blue") 
plot(log10(metrics$lambda), metrics$mbic, type="l")
range(metrics$lambda[which(metrics$mbic==min(metrics$mbic))]) # 1.0000 251.1886



### 2.2.2 Comparison with abess mgaussian best subset glm #### 

# This was our L0 glm metrics:
metrics_L0glm = do.call(cbind, get.metrics(beta = L0glm.out_preset_coefs, 
                             X = X, y = y, 
                             family = family, 
                             beta_true = sim$beta_true))
metrics_L0glm
#        p   n   m deviance deviance_true   loglik MD_coefs  MAD_coefs sMAPE_coefs bias_coefs  R2_coefs
# [1,] 501 500 161 48273.95      8403.259 -99995.1 0.018487 0.05023418  0.08479029   1.000029 0.9999976
#      support_size_true support_size_estimated estimated_over_true_support_size TP FP  TN FN precision recall
# [1,]                50                     50                                1 50  0 450  0         1      1
#      specificity accuracy F1_score    gic     ebic      bic      aic   mbic
# [1,]           1        1        1 277832 285207.4 242585.1 213698.2 308799

# abess mgaussian best subset fit plus metrics:
# best subset multitask learning using family="mgaussian" using abess
library(abess)
system.time(abess_fit <- abess(x = X, 
                   y = y, 
                   # weights = 1/(y+0.1),  # NOTE: can be used as approx 1/variance weights instead of using family="poisson"
                   family = "mgaussian",
                   tune.path = "sequence", 
                   support.size = c(1:(p-1)),
                   lambda=0, 
                   warm.start = TRUE,
                   tune.type = "gic")) # 0.31s for n=p=100, 21s for np=500, or cv or ebic or bic or aic
# note: currently no support for nonnegativity constraints
plot(abess_fit, type="tune")  # optimal support size 
beta_abess = as.matrix(extract(abess_fit)$beta) # coefficient matrix for best subset
library(qlcMatrix)
image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=beta_abess^0.1, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="abess multitask mgaussian (red=true peaks)")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim$X))[as.vector(rowMin(beta_abess)!=0)], col="cyan")
sum(rowMax(sim$beta_true)!=0) # 10 true peaks
sum(rowMin(beta_abess)!=0) # 29 peaks detected - no good, mainly because abess currently does not support nonnegativity constraints
# note: for larger p and lower sparsity, solution quality quickly becomes a lot worse
# see https://github.com/abess-team/abess/issues/510

metrics_abess = do.call(cbind, get.metrics(beta = beta_abess, 
                                           X = X, y = y, 
                                           family = family, 
                                           beta_true = sim$beta_true))
metrics_abess
#        p   n   m deviance deviance_true  loglik    MD_coefs MAD_coefs sMAPE_coefs bias_coefs   R2_coefs
# [1,] 501 500 161 107224.6      82175.66 -150540 -0.02288781  62.32001   0.6122153  0.0966366 0.08835559
#      support_size_true support_size_estimated estimated_over_true_support_size TP  FP  TN FN precision recall
# [1,]                50                    157                             3.14 34 123 327 16 0.2165605   0.68
#      specificity accuracy  F1_score      gic     ebic      bic      aic   mbic
# [1,]   0.7266667    0.722 0.3285024 588165.6 615366.7 458172.8 351635.9 702374


### 2.2.3 Comparison with glmnet nonnegative mgaussian multivariate group LASSO fit #### 

library(glmnet)
system.time(cvfit <- cv.glmnet(X, 
                               y, 
                               alpha = 1,  # lasso
                               family = "mgaussian", 
                               nlambda = 100,
                               nfolds = 5, 
                               standardize = FALSE,
                               standardize.response = FALSE,
                               intercept = FALSE, 
                               relax = FALSE,
                               lower.limits = rep(0, ncol(X)+1))) # for nonnegativity constraints
# 6.32s for n=p=100, 27s for n=p=500

plot(cvfit)

# Fit the nonnegative group LASSO glmnet model over a sequence of lambda values
system.time(fit_glmnet <- glmnet(X, 
                                 y,
                                 alpha = 1, # lasso
                                 family = "mgaussian", 
                                 nlambda = 100,
                                 standardize = FALSE,
                                 standardize.response = FALSE,
                                 intercept = FALSE, 
                                 relax = FALSE, 
                                 lower.limits = rep(0, ncol(X)+1))) # for nonnegativity constraints
# 3.11s for n=p=100, 3.45s for n=p=500

# Select the best lambda value (I am using cvfit$lambda.1se to get slightly sparser solution)
best_lambda <- cvfit$lambda.1se

plot(fit_glmnet, xvar="lambda")

# get the coefficients for each task for the best lambda value
coefs <- coef(fit_glmnet, s = best_lambda)
beta_mgaussian_glmnet <- do.call(cbind, lapply(seq_len(m), function (channel) as.matrix(coefs[[channel]][-1,,drop=F])))
colnames(beta_mgaussian_glmnet) = 1:ncol(beta_mgaussian_glmnet)
beta_mgaussian_glmnet[abs(beta_mgaussian_glmnet)<1E-5] = 0 # slight amount of thresholding of small coefficients using same thresh as used in L0glm

image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=beta_mgaussian_glmnet^0.1, col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="nonnegative group Lasso glmnet (red=true support)")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim$X))[as.vector(rowMin(beta_mgaussian_glmnet)!=0)], col="cyan")
sum(rowMax(sim$beta_true)!=0) # 10 true peaks
sum(rowMin(beta_mgaussian_glmnet)!=0) # 13 peaks detected

metrics_glmnet = do.call(cbind, get.metrics(beta = beta_mgaussian_glmnet, 
                                           X = X, y = y, 
                                           family = family, 
                                           beta_true = sim$beta_true))
metrics_glmnet
#        p   n   m deviance deviance_true    loglik   MD_coefs MAD_coefs sMAPE_coefs bias_coefs  R2_coefs
# [1,] 501 500 161 15972207      15920873 -336848.1 -0.3822619  2.590571   0.1882265   1.050774 0.9650188
#      support_size_true support_size_estimated estimated_over_true_support_size TP FP  TN FN precision recall
# [1,]                50                     72                             1.44 48 24 426  2 0.6666667   0.96
#      specificity accuracy  F1_score      gic     ebic      bic      aic   mbic
# [1,]   0.9466667    0.948 0.7868852 789402.6 800365.7 737010.6 694072.1 835433

# note: for larger p solution quality quickly becomes a lot worse
# see https://github.com/abess-team/abess/issues/510

