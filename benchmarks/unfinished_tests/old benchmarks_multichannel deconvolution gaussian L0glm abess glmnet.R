# library(remotes)
# remotes::install_github("tomwenseleers/L0glm/L0glm")
# library(L0glm)
library(export)
devtools::load_all()

# EXAMPLE MULTICHANNEL DECONVOLUTION ####
# simulate blurred multichannel spike train with Gaussian noise ####
set.seed(1)
s <- 0.1 # sparsity (% of timepoints where there is a peak / variables that are nonnegative)
p <- 200 # nr of variables
# simulate multichannel blurred spike train with Gaussian noise
sd_noise <- 10
sim <- simulate_spike_train(n=p, p=p, k=round(s*p), 
                            mean_beta = 10000,
                            sd_logbeta = 1,
                            family="gaussian", # alt option "poisson" for nonnegative identity link poisson
                            sd_noise = sd_noise,
                            multichannel=TRUE, 
                            sparse=TRUE)
# graph2png(file="../benchmarks/plots/1 simulated multichannel spike train.png", width=7, height=5)
X <- sim$X # n x p matrix
Y <- sim$y # n x m matrix  
true_coefs <- sim$beta_true
true_support_set <- which(rowMax(sim$beta_true)!=0)
m <- ncol(Y) # nr of tasks
n <- nrow(X) # nr of observations
p <- ncol(X) # nr of independent variables (shifted copies of point spread functions)
W <- 1/(Y+0.1) # approx 1/variance Poisson observation weights with family="poisson", n x m matrix
family <- gaussian(identity)

# relaxed group L0 penalized GLM performance metrics over lambda sequence  ####

# lambda values that in principle would be expected to maximize aic, bic, gic, ebic or mbic 
# can be calculated using
# (NOTE: with high dispersion/variance these come out way too high though,  except if sd_noise=1 - check this)

# if dispersion/variance is known
lambda_preset_aic_knowndisp = preset.lambda(IC="aic", y=Y, X=X, dispersion=sd_noise^2)  
lambda_preset_aic_knowndisp 
lambda_preset_bic_knowndisp = preset.lambda(IC="bic", y=Y, X=X, dispersion=sd_noise^2)
lambda_preset_bic_knowndisp
lambda_preset_gic_knowndisp = preset.lambda(IC="gic", y=Y, X=X, dispersion=sd_noise^2)  
lambda_preset_gic_knowndisp 
lambda_preset_ebic_knowndisp = preset.lambda(IC="ebic", y=Y, X=X, dispersion=sd_noise^2) 
lambda_preset_ebic_knowndisp
lambda_preset_ebic_knowndisp = preset.lambda(IC="mbic", y=Y, X=X, dispersion=sd_noise^2) 
lambda_preset_ebic_knowndisp
# NOTE: if dispersion/variance is not known it can first be estimated
lambda_preset_aic = preset.lambda(IC="aic", y=Y, X=X, 
              family=gaussian(identity))  
lambda_preset_aic 
lambda_preset_bic = preset.lambda(IC="bic", y=Y, X=X, 
              family=gaussian(identity))
lambda_preset_bic 
lambda_preset_gic = preset.lambda(IC="gic", y=Y, X=X, 
              family=gaussian(identity))
lambda_preset_gic 
lambda_preset_ebic = preset.lambda(IC="ebic", y=Y, X=X, 
              family=gaussian(identity))
lambda_preset_ebic 
lambda_preset_mbic = preset.lambda(IC="mbic", y=Y, X=X, 
              family=gaussian(identity))  
lambda_preset_mbic 

# below we fit over sequence of lambda values to see if these IC are indeed optimal at these values
# NOTE: we could also use golden section search
# using https://search.r-project.org/CRAN/refmans/cmna/html/goldsect.html
# as in abess

# for the lambda sequence we will below use lambda_preset_aic 
# (round(log10(lambda_preset_aic))-3) as the lower bound and the theoretical
# lambda that would penalize all variables to zero with an L0 penalty as maximum
# see https://stats.stackexchange.com/questions/416144/minimum-and-maximum-regularization-in-l0-pseudonorm-penalized-regression
maxlambda = maxlambda(X=X, weights=W, z=Y, family=family) 
# lambda that would penalize all variables to 0 in L0 regression 

lambdas = 10^seq((round(log10(lambda_preset_aic))-3), round(log10(maxlambda))+1, by=1)

library(parallel)
library(qlcMatrix)
system.time(warmup <- L0glm.fit(X = X, y = Y,
                               family = gaussian(identity), # poisson(identity),
                               lambda = min(lambdas), #  # or 0.1 or 1
                               start = rep(1, ncol(X)),
                               nonnegative = TRUE,
                               thresh = 1,
                               grouped = T, # use group L0 penalty
                               group.maxit = 5,
                               relaxed = T, # return unpenalized fit for selected variables at the end, similar to concept of relaxed LASSO
                               warmupiters = 5,
                               solver = "eigen",
                               method = 1L # using Eigen's sparse Cholesky solver
                               )) # 5s
sum(rowMin(abs(warmup$coefficients))!=0) # 151, nr of selected variables at lowest lambda used

warm_start = warmup$coefficients
no_warm_start = matrix(1, nrow=nrow(warm_start), ncol=ncol(warm_start))

w = 0.9 # 0 = no warm starts, 1 = full warm starts
system.time(coefficients <- lapply(lambdas, function (lambda) {
                                new_coefs <- L0glm.fit(X = sim$X, y = sim$y, 
                                                                            family = gaussian(identity), # poisson(identity),
                                                                            start = w*warm_start+(1-w)*no_warm_start*1, # use warm starts
                                                                            lambda = lambda, # needs to be larger with larger sd_noise for family=gaussian, e.g. 0.5 with sd_noise=1 & 5 with sd_noise=5 
                                                                            nonnegative = TRUE,
                                                                            thresh = 1, # coefficients below this threshold are set to zero
                                                                            grouped = T, # use group L0 penalty
                                                                            group.maxit = 6,
                                                                            relaxed = T, # return unpenalized fit for selected variables at the end, similar to concept of relaxed LASSO
                                                                            warmupiters = 5,
                                                                            solver = "eigen",
                                                                            method = 1L, # using Eigen's sparse Cholesky solver
                                                                            parallellize = T 
                                                                            )$coefficients
                                warm_startt <<- new_coefs
                                print(length(new_coefs))
                                avgcoefs <- as.vector(rowMeans(new_coefs))
                                print(length(avgcoefs))
                                return(avgcoefs)
}
)) # 83s
# metrics relaxed group L0 over sequence of lambda values
library(qlcMatrix)
metrics <- data.frame(do.call(rbind, lapply(coefficients, function (coefs) get.metrics(beta=coefs, X=X, y=Y, 
                                                                                       family=gaussian(identity), 
                                                                                       beta_true=sim$beta_true)))) 


system.time(metrics2 <- L0glm.tune(X = sim$X, y = sim$y, family = gaussian(identity), nonnegative = T, coefs.true = sim$beta_true, IC = "aic", group.maxit = 6, plot = T, grouped = T, relaxed = T, parallellize=T))


# best empirical lambda based on calculated GIC over lambda sequence:
metric = "gic"
bestlambda = lambdas[which.min(metrics2$metrics[,metric])]
bestlambda # 100

bestlambda = lambdas[which.min(abs(lambdas-bestlambda))] 
beta_relaxedgroupL0 = as.matrix(coefficients[[which(lambdas==bestlambda)]]) # estimated coefficients at best lambda
# beta_relaxedgroupL0[beta_relaxedgroupL0<10] = 0

image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=as.matrix(sim$beta_true)^0.1, 
      col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="True coefficients\n(red=true support)")
abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")
# abline(v=(1:nrow(sim$X))[as.vector(rowMin(beta_relaxedgroupL0)!=0)], col="magenta")
# graph2png(file="../benchmarks/plots/2 true coefficients true support.png", width=7, height=5)

image(x=1:nrow(sim$y), y=1:ncol(sim$y), z=beta_relaxedgroupL0^0.1, 
      col = topo.colors(255),
      useRaster=TRUE,
      xlab="Time", ylab="Channel", main="nonnegative relaxed group L0\n(red=estimated support)")
# abline(v=(1:nrow(sim$X))[as.vector(rowMax(sim$beta_true)!=0)], col="red")
abline(v=(1:nrow(sim$X))[as.vector(rowMin(beta_relaxedgroupL0)!=0)], col="red")
# graph2png(file="../benchmarks/plots/3 estimated coefficients and support nonneg relaxed group L0.png", width=7, height=5)
sum(rowMax(sim$beta_true)!=0) # true support size of 20
sum(rowMin(beta_relaxedgroupL0)!=0) # estimate support size correct: also 20
length(intersect(which(rowMax(sim$beta_true)!=0), 
                 which(rowMin(beta_relaxedgroupL0)!=0))) # 20 out of 20 variables correctly selected, no false pos & no false negs
length(which(sim$beta_true!=0)) # 1606 coefs nonzero
length(as.vector(sim$beta_true)) # 29800 coefs in total
length(intersect(which(sim$beta_true!=0), 
                 which(beta_relaxedgroupL0!=0))) # 1464 true pos coefs
length(intersect(which(sim$beta_true==0), 
                 which(beta_relaxedgroupL0!=0))) # 968 false pos coefs
plot(as.vector(beta_relaxedgroupL0), as.vector(sim$beta_true), pch=16, 
     xlab="estimated coefficients nonnegative relaxed L0 group GLM",
     ylab="true coefficients", col="steelblue",
     main="True vs estimated coefficients\nnonnegative relaxed L0 group GLM\nR2=0.9988")
# graph2png(file="../benchmarks/plots/4 true vs estimated coefficients relaxed group L0.png", width=7, height=5)
cor(as.vector(sim$beta_true), as.vector(beta_relaxedgroupL0))^2 # R2=0.9987671


plot(x=log10(lambdas), y=metrics$support_size_estimated,
     ylab="estimated support size",
     main="Estimated support size in function of regularisation\nusing nonnegative relaxed group L0 GLM\n(red=estimate at best GIC, blue=true support size)",
     pch=16, type="o")
abline(v=log10(bestlambda), col="red")
abline(h=metrics$support_size_true, col="blue")
# graph2png(file="../benchmarks/plots/5 support size vs lambda relaxed group L0.png", width=7, height=5)

plot(x=log10(lambdas), y=metrics$deviance, pch=16, type="o", 
     ylim=c(0,max(unlist(metrics$deviance))),
     main = "Theoretical deviance relative to true underlying signal without noise\n(red=estimate at best GIC, blue=best deviance" )  # theoretical deviance relative to true underlying signal without noise
lambda_bestdeviance = lambdas[which.min(metrics$deviance)]
lambda_bestdeviance
abline(v=log10(bestlambda), col="red")
abline(v=log10(lambda_bestdeviance), col="blue")
par(oma=c(0,2,0,0))

plot(x=log10(lambdas), y=metrics$loglik,
     main="Log likelihood in function of regularisation\nusing nonnegative relaxed group L0 GLM\n(red=estimate at best GIC)",
     ylab="log likelihood",
     pch=16, type="o") # log likelihood
lambda_bestloglik = lambdas[which.max(metrics$loglik)]
lambda_bestloglik
abline(v=log10(bestlambda), col="red")
# abline(v=log10(lambda_bestloglik), col="blue")
# graph2png(file="../benchmarks/plots/6 log likelihood vs lambda relaxed group L0.png", width=7, height=5)

plot(x=log10(lambdas), y=metrics$MD_coefs, pch=16, type="o")
abline(v=log10(bestlambda), col="red")


plot(x=log10(lambdas), y=metrics$MAD_coefs, 
     main="Mean absolute deviation of coefficients vs true coefs\nusing nonnegative relaxed group L0 GLM\n(red=estimate at best GIC)",
     ylab="Mean absolute deviation coefficients",
     pch=16, type="o")
lambdas[which.min(metrics$MAD_coefs)] 
abline(v=log10(bestlambda), col="red")
# graph2png(file="../benchmarks/plots/7 mean abs deviation true coefficients vs lambda relaxed group L0.png", width=7, height=5)

plot(x=log10(lambdas), y=metrics$sMAPE_coefs, pch=16, type="o")
lambdas[which.min(metrics$sMAPE_coefs)]
abline(v=log10(bestlambda), col="red")

plot(x=log10(lambdas), y=metrics$bias_coefs, pch=16, type="o", ylim=c(0.9,1.1)) # closest to 1 best
lambda_minimal_bias = min(lambdas[which(abs(unlist(metrics$bias_coefs)-1)<=min(abs(unlist(metrics$bias_coefs)-1), na.rm=T)*1.3)]) # 1
lambda_minimal_bias 
abline(v=log10(bestlambda), col="red")
abline(h=1, col="blue")

plot(x=log10(lambdas), y=metrics$R2_coefs, pch=16, type="o", ylim=c(0, 1)) # highest lambda before sudden drop in R2 optimal
abline(v=log10(bestlambda), col="red") # best would be lambda of 1

plot(x=log10(lambdas), y=metrics$estimated_over_true_support_size, pch=16, type="o") # lambda 1 close to optimal
lambda_correct_support_size = min(lambdas[which.min(abs(unlist(metrics$estimated_over_true_support_size)-1))]) #
lambda_correct_support_size 
abline(v=log10(bestlambda), col="red") 
abline(v=log10(lambda_correct_support_size), col="blue")

plot(x=log10(lambdas), y=metrics$precision, pch=16, type="o")
abline(v=log10(bestlambda), col="red")

plot(x=log10(lambdas), y=metrics$recall, pch=16, type="o")
abline(v=log10(bestlambda), col="red")

plot(x=log10(lambdas), y=metrics$specificity, pch=16, type="o")
abline(v=log10(bestlambda), col="red")

plot(x=log10(lambdas), y=metrics$accuracy,
     main="Accuracy to recover true support\nusing nonnegative relaxed group L0 GLM\n(red=estimate at best GIC)",
     ylab="accuracy",
     pch=16, type="o") # GOOD FOR THEORETICAL BEST SUPPORT RECOVERY
lambda_best_accuracy = lambdas[which.max(metrics$accuracy)] # 
lambda_best_accuracy
abline(v=log10(bestlambda), col="red")
# abline(v=log10(lambda_best_accuracy), col="blue")
# graph2png(file="../benchmarks/plots/8 accuracy support recovery vs lambda relaxed group L0.png", width=7, height=5)

plot(x=log10(lambdas), y=metrics$F1_score, pch=16, type="o") # BEST FOR THEORETICAL BEST SUPPORT RECOVERY, has max at lambda=1
lambda_best_F1score = lambdas[which.max(metrics$F1_score)] # 
lambda_best_F1score # 2.8
abline(v=log10(bestlambda), col="red")
abline(v=log10(lambdas[which.max(metrics$F1_score)]), col="blue")

plot(x=log10(lambdas), y=metrics$ebic, pch=16, type="o") # keep, cf https://www.researchgate.net/publication/46468889_Consistency_of_the_Empirical_Bayes_Information_Criterion_for_Selecting_Variables_in_Linear_Mixed_Models#:~:text=Compared%20to%20mAIC%2C%20cAIC%20and%20BIC%2C%20the%20EBIC,of%20covariates%20and%20the%20number%20of%20clusters%20increase.
lambda_best_EBIC = lambdas[which.min(metrics$ebic)] #
lambda_best_EBIC # 1
abline(v=log10(bestlambda), col="red")
abline(v=log10(lambda_best_EBIC), col="blue") # corresponds to lambda giving optimal support recovery

plot(x=log10(lambdas), y=metrics$gic, pch=16, type="o") # keep
lambda_best_gic = lambdas[which.min(metrics$gic)] # 1
abline(v=log10(bestlambda), col="red")
abline(v=log10(lambda_best_gic), col="blue") # corresponds to lambda giving optimal support recovery

plot(x=log10(lambdas), y=metrics$aic, pch=16, type="o") # keep - no good though
lambda_best_aic = lambdas[which.min(metrics$aic)] # 0.088 - no good
abline(v=log10(bestlambda), col="red")
abline(v=log10(lambda_best_aic), col="blue") #

plot(x=log10(lambdas), y=metrics$bic, pch=16, type="o") # keep
lambda_best_bic = lambdas[which.min(metrics$bic)] # 1
abline(v=log10(bestlambda), col="red")
abline(v=log10(lambda_best_aic), col="blue") # keep


plot(x=log10(lambdas), y=metrics$mbic, pch=16, type="o") # keep
lambda_best_mbic = lambdas[which.min(metrics$mbic)] # 1
abline(v=log10(bestlambda), col="red")
abline(v=log10(lambda_best_mbic), col="blue")

# determine union of selected variables by L0GLM fit over several bootstrap replicates (USVB)
# this set could be used to select variables & do inference using the base R GLM function
system.time({
  B <- 100 # nr of bootstrap replicates, 100 here for the sake of the argument
  beta_relaxedgroupL0_boot <- vector("list", B) # list of coefficient matrices for each bootstrap replicate
  set.seed(1)
  for (i in 1:B) { 
    sampled = sample(sim$x, replace=T) # sampled observations
    beta_relaxedgroupL0_boot[[i]] = L0glm.fit(X = sim$X[sampled,,drop=F], y = sim$y[sampled,,drop=F], 
                                  family = gaussian(identity), # poisson(identity),
                                  # start = w*warm_start+(1-w)*no_warm_start, # to use warm starts
                                  lambda = bestlambda, 
                                  nonnegative = TRUE,
                                  thresh = 1, # coefficients below this threshold are set to zero
                                  delta = 0, # TODO: if we have the thresh argument we can get rid of delta & set delta to zero, no need to have both
                                  group = T, # use group L0 penalty
                                  group.maxit = 6,
                                  relaxedGroupL0 = T, # return unpenalized fit for selected variables at the end, similar to concept of relaxed LASSO
                                  warmupiters = 5,
                                  solver = "eigen",
                                  method = 1L, # using Eigen's sparse Cholesky solver
                                  parallellize = T 
  )
}}
) # 910s

# do inference by calculating the
# union of selected variables across L0glm best subset fits on bootstrap replicates
# and then fitting a (constrained) (g)lm on that set of selected variables
boot_selectedvars_relaxedGroupL0 = lapply(beta_relaxedgroupL0_boot, function (coefs) which(rowMin(coefs)!=0)) # set of selected variables for each fit
USVBset_relaxedGroupL0 = unique(unlist(boot_selectedvars_relaxedGroupL0)) # union of selected variables
length(USVBset_relaxedGroupL0) # 35 variables selected out of 200 (true support size=20)

USVBset_relaxedGroupL0_setsize_vs_B = sapply(1:B, function (nr_bootstrap_repls) { 
                                                  length(unique(unlist(boot_selectedvars_relaxedGroupL0[1:nr_bootstrap_repls]))) }) 
plot(1:B, USVBset_relaxedGroupL0_setsize_vs_B, type="o", pch=16, 
     main="Size of union of variables selected by nonnegative\nrelaxed group L0 GLM after bootstrapping\n(blue=true support size)",
     xlab="Number of bootstrap replicates",
     ylab="Size of union of selected variables")
abline(h=length(true_support_set), col="blue")
# graph2png(file="../benchmarks/plots/9 size of union of selected variables after bootstrap relaxed group L0.png", width=7, height=5)

# gaussian generalised linear model fits on union of selected variables for each outcome variable ####
# since we have nonnegativity constraints I use glm.cons here
library(zetadiv)
library(broom)

signvariables = lapply(1:m, # variables that are significant for each of the m channels (after Bonferroni correction)
       function (var) {
         glmconsfit_on_subset = glm.cons(sim$y[,var] ~ 0+as.matrix(sim$X[,USVBset_relaxedGroupL0,drop=F]), 
                                         family=stats::gaussian(identity),
                                         cons=1, # for nonnegativity constraints
                                         cons.inter=1)
         which(tidy(glmconsfit_on_subset)$p.value<(0.05/m)) }) # /m for Bonferroni correction
set_signvariables = unique(unlist(signvariables))
length(set_signvariables) # 17 variables with coefficient sign >0 in at least one of the channels
intersect(set_signvariables, true_support_set) # 15 out of 17 are in the true support set, so only 2 false positives


# a regular unconstrained linear model fit on all variables obviously 
# would give really crap results, given that p=n & that such a fit
# does not take into account nonnegativity constraints
system.time(beta_lm_onallvars <- lm.fit(x=as.matrix(sim$X), y=as.matrix(sim$y))$coefficients)
range(beta_lm_onallvars) # -52381709  52281909
range(sim$beta_true) # 
plot(beta_lm_onallvars, sim$beta_true, pch=16, 
     xlab="estimated coefficients linear model",
     ylab="true coefficients", col="steelblue",
     main="True vs estimated coefficients\nleast square regression fit on all variables")
# graph2png(file="../benchmarks/plots/10 true vs estimated coefficients least square regr all variables.png", width=7, height=5)



# Alternatively, similar to the bolasso we could also identify the variables that 
# are selected by the L0glm fit at least threshold % of the time, similar to what 
# the https://cran.r-project.org/web/packages/bolasso/bolasso.pdf
# does. If we set threshold to 0.3, this seems to converge very quickly to approx the true support.
# (PROBABLY THE BEST SOLUTION)
threshold = 0.3 # value much lower than in bolasso, where threshold to achieve consistency=0.9
# variables selected at least threshold % of the time:
BOL0GLMset = USVBset_relaxedGroupL0[colMeans(do.call(rbind,lapply(1:length(beta_relaxedgroupL0_boot), function (i) USVBset_relaxedGroupL0 %in% boot_selectedvars_relaxedGroupL0[[i]]))) >= threshold]
length(BOL0GLMset) # 20 variables selected
length(intersect(BOL0GLMset, true_support_set)) # 19/20 selected variables overlap with true support, 1 false positive but highly collinear with the correct variable

BOL0GLMset_setsize_vs_B = sapply(1:B, function (nr_bootstrap_repls) { length(USVBset_relaxedGroupL0[colMeans(do.call(rbind,lapply(1:nr_bootstrap_repls, function (i) USVBset_relaxedGroupL0 %in% boot_selectedvars_relaxedGroupL0[[i]]))) >= threshold]) })
plot(1:B, 
     BOL0GLMset_setsize_vs_B, type="o", pch=16,
     main="Size of set of variables selected 30% of the time by\nnonnegative relaxed group L0 GLM after\nbootstrapping (blue=true support size)",
     xlab="Number of bootstrap replicates",
     ylab="Size of set of selected variables")
abline(h=length(true_support_set), col="blue")
# graph2png(file="../benchmarks/plots/9b size of set of variables selected 30 percent of the time after bootstrap relaxed group L0.png", width=7, height=5)


# we could then do inference on these selected variables using glm.cons or regular glm



# NOTE: if you use nonnegative group LASSO one also has to use this approach, as in https://cran.r-project.org/web/packages/bolasso/bolasso.pdf ####
# but with threshold = 0.9 to achieve variable selection consistency

# determine optimal lambda in nonnegative group LASSO
library(glmnet)
system.time(cvfit <- cv.glmnet(X, 
                               Y, 
                               alpha = 1,  # lasso
                               family = "mgaussian",  # for group LASSO
                               nlambda = 100,
                               # lambda.min.ratio = 1E-4,
                               # thresh = 1e-10,
                               nfolds = 5, 
                               standardize = FALSE,
                               standardize.response = FALSE,
                               intercept = FALSE, 
                               relax = FALSE,
                               lower.limits = rep(0, ncol(X)+1))) # 15s
plot(cvfit)

# Select the best lambda value (I am using cvfit$lambda.1se to get slightly sparser solution)
cvfit$lambda.1se
cvfit$lambda.min 
best_lambda_mgaussian_glmnet <- cvfit$lambda.1se

# Fit the glmnet model over a sequence of lambda values
system.time(fit_glmnet <- glmnet(
  sim$X, 
  sim$y,
  alpha = 1,
  # lasso
  family = "mgaussian",
  nlambda = 100,
  # lambda.min.ratio = 1E-4,
  # thresh = 1e-10,
  standardize = FALSE,
  standardize.response = FALSE,
  intercept = FALSE,
  relax = FALSE,
  lower.limits = rep(0, ncol(X) + 1)  # for nonnegativity constraints
)) # 3s


# get the coefficients for each task for the best lambda value
coefs <- coef(fit_glmnet, s = best_lambda_mgaussian_glmnet)
beta_mgaussian_glmnet <- predict(fit_glmnet, s = best_lambda_mgaussian_glmnet, newx=sim$X)
plot(as.vector(beta_mgaussian_glmnet), as.vector(sim$beta_true), pch=16, 
     xlab="estimated coefficients nonnegative group LASSO glmnet",
     ylab="true coefficients", col="steelblue",
     main="True vs estimated coefficients\nnonnegative group LASSO glmnet\nR2=0.24")
# graph2png(file="../benchmarks/plots/11 true vs estimated coefficients relaxed group L0.png", width=7, height=5)
cor(as.vector(sim$beta_true), as.vector(beta_mgaussian_glmnet))^2 # R2=0.2364216


# determine union of  variables selected by group LASSO glmnet fit over several bootstrap replicates (USVB)
system.time({
  B <- 100 # nr of bootstrap replicates, 100 here for the sake of the argument
  beta_mgaussian_glmnet_boot <- vector("list", B) # list of coefficient matrices for each bootstrap replicate
  set.seed(1)
  for (i in 1:B) { 
    sampled = sample(sim$x, replace=T) # sampled observations
    
    # Fit the glmnet model over a sequence of lambda values
    fit_glmnet <- glmnet(
      sim$X[sampled,,drop=F], 
      sim$y[sampled,,drop=F],
      alpha = 1,
      # lasso
      family = "mgaussian",
      nlambda = 100,
      # lambda.min.ratio = 1E-4,
      # thresh = 1e-10,
      standardize = FALSE,
      standardize.response = FALSE,
      intercept = FALSE,
      relax = FALSE,
      lower.limits = rep(0, ncol(X) + 1)  # for nonnegativity constraints
    )
    
    # get the coefficients for each task for the best lambda value
    coefs <- coef(fit_glmnet, s = best_lambda_mgaussian_glmnet)
    beta_mgaussian_glmnet <- do.call(cbind, lapply(seq_len(m), function (channel) as.matrix(coefs[[channel]][-1,,drop=F])))
    colnames(beta_mgaussian_glmnet) = 1:ncol(beta_mgaussian_glmnet)
    beta_mgaussian_glmnet[abs(beta_mgaussian_glmnet)<1] = 0 # same thresholding as used in L0 glm fits
    
    # store coefficients in list
    beta_mgaussian_glmnet_boot[[i]] = beta_mgaussian_glmnet 
  }}
) # 427s

# union of selected variables across L0glm best subset fits on bootstrap replicates
# and then fitting a (constrained) (g)lm on that set of selected variables
# this doesn't work as it would select too large a set & set would keep on growing
boot_selectedvars_mgaussian_glmnet = lapply(beta_mgaussian_glmnet_boot, function (coefs) which(rowMin(coefs)!=0)) # set of selected variables for each fit
USVBset_mgaussian_glmnet = unique(unlist(boot_selectedvars_mgaussian_glmnet)) # union of selected variables
length(USVBset_mgaussian_glmnet) # 181 variables selected out of 200 (true support size=20)
# so this doesn't work!

USVBset_mgaussian_glmnet_setsize_vs_B = sapply(1:B, function (nr_bootstrap_repls) { 
  length(unique(unlist(boot_selectedvars_mgaussian_glmnet[1:nr_bootstrap_repls]))) }) 
plot(1:B, USVBset_mgaussian_glmnet_setsize_vs_B, type="o", pch=16, 
     main="Size of union of variables selected by nonnegative\ngroup LASSO glmnet fit after bootstrapping\n(blue=true support size)",
     xlab="Number of bootstrap replicates",
     ylab="Size of union of selected variables",
     ylim=c(0,200))
abline(h=length(true_support_set), col="blue")
# graph2png(file="../benchmarks/plots/12 size of union of selected variables after bootstrap nonnegative group LASSO.png", width=7, height=5)

# set of variables selected 70% of the time, as in https://arxiv.org/abs/0804.1302, 
# does work & achieves variable selection consistency
threshold = 0.7 # optimal value here to achieve variable selection consistency, slightly lower than suggested 0.9
# variables selected at least threshold % of the time:
bolassoset = USVBset_mgaussian_glmnet[colMeans(do.call(rbind,lapply(1:length(beta_mgaussian_glmnet_boot), function (i) USVBset_mgaussian_glmnet %in% boot_selectedvars_mgaussian_glmnet[[i]]))) >= threshold]
length(bolassoset) # 17 variables selected
length(intersect(bolassoset, true_support_set)) # 17 selected variables all overlap with true support, no false positives

bolassoset_setsize_vs_B = sapply(1:B, function (nr_bootstrap_repls) { length(USVBset_mgaussian_glmnet[colMeans(do.call(rbind,lapply(1:nr_bootstrap_repls, function (i) USVBset_mgaussian_glmnet %in% boot_selectedvars_mgaussian_glmnet[[i]]))) >= threshold]) })
plot(1:B, 
     bolassoset_setsize_vs_B, type="o", pch=16,
     main="Size of set of variables selected 70% of the time by\nnonnegative group LASSO after\nbootstrapping (blue=true support size)",
     xlab="Number of bootstrap replicates",
     ylab="Size of set of selected variables")
abline(h=length(true_support_set), col="blue")
# graph2png(file="../benchmarks/plots/12b size of set of variables selected 70 percent of the time after bootstrap nonneg group LASSO.png", width=7, height=5)





# abess MULTITASK mgaussian FIT OVER SEQUENCE OF SUPPORT SIZES ####
# note: this also gives poor results, probably mainly because right now it
# doesn't yet support nonnegativity constraints; only gives a good
# solution if you pre-select variables using nonnegative group LASSO glmnet fit
# but not ideal (as there would also be uncertainty on this pre-selection)

library(abess)
allvars = 1:ncol(X)
preselvars = allvars # which(rowMin(beta_mgaussian_glmnet)!=0) # pre-select variables based on nonnegative mgaussian group LASSO
system.time(abess_fit <- abess(x = X[,preselvars,drop=F], 
                               y = Y, 
                               # weights = 1/(Y+0.1), # 1/variance Poisson weights do work in abess if needed
                               family = "mgaussian",
                               tune.path = "sequence", # or "gsection"
                               support.size = c(1:150),
                               # gs.range = c(1,200),
                               lambda=0, 
                               max.splicing.iter = 20,
                               warm.start = TRUE,
                               tune.type = "ebic") # ebic gave sparsest solution; alt options cv or bic or aic or gic
) # 1.09s


# abess_fit
extract(abess_fit)$support.size # 31
sum(rowMin(abs(sim$beta_true))!=0) # true support 20
plot(abess_fit) # quite unstable coefficient paths
plot(abess_fit, type="tune")

# metrics multitask mgaussian abess fit
support.size <- abess_fit$support.size
abess_fit$edf
coefficients <- abess_fit$beta # lapply(abess_fit$beta, function (coefs) { out = Matrix(0, nrow=p, ncol=m, sparse=T) } )
library(qlcMatrix)
metrics <- data.frame(do.call(rbind, lapply(1:length(coefficients), function (i) { coefs=coefficients[[i]]
                                                                             # edf=sum(rowMin(abs(coefs))!=0) # *ncol(coefs)
                                                                             coefs[abs(coefs)<1] = 0; # same thresholding as in L0glm
                                                                             get.metrics(beta=coefs, X=sim$X, y=sim$y, 
                                                                                       weights=matrix(1, nrow=nrow(sim$y), ncol=ncol(sim$y)), 
                                                                                       offset=Matrix(0, nrow=nrow(sim$X), ncol=ncol(sim$y), sparse=T), 
                                                                                       family=gaussian(identity), 
                                                                                       beta_true=sim$beta_true,
                                                                                       edf=sum(coefs!=0), 
                                                                                       lower=rep(0, ncol(sim$X)),
                                                                                       upper=rep(Inf,ncol(sim$X))
                                                                                       ) } )))
sel=1:100

# TODO: I need to check the calculation of these metrics, since the number I get for my
# information criteria & the curves look different than the ones returned by abess 
plot(x=support.size[sel], y=metrics$deviance[sel], pch=16, type="o", ylim=c(0,max(unlist(metrics$deviance[sel])))) # highest lambda before sudden jump in deviance corresponds with optimal lambda
abline(v=metrics$support_size_true, col="red")
plot(x=support.size[sel], y=metrics$loglik[sel], pch=16, type="o")
abline(v=metrics$support_size_true, col="red")
plot(x=support.size[sel], y=metrics$MD_coefs[sel], pch=16, type="o")
abline(v=metrics$support_size_true, col="red")
plot(x=support.size[sel], y=metrics$MAD_coefs[sel], pch=16, type="o") # optimal at support=11
support.size[which.min(metrics$MAD_coefs)] # 12
abline(v=metrics$support_size_true, col="red")
plot(x=support.size[sel], y=metrics$sMAPE_coefs[sel], pch=16, type="o") # optimal at support=11
support.size[which.min(metrics$sMAPE_coefs)] # 12
abline(v=metrics$support_size_true, col="red")
plot(x=support.size[sel], y=metrics$bias_coefs[sel], pch=16, type="o", ylim=c(0.1,2)) #
support.size[which.min(metrics$bias_coefs)]
abline(v=metrics$support_size_true, col="red")
abline(h=1, col="blue")
plot(x=support.size[sel], y=metrics$R2_coefs[sel], pch=16, type="o", ylim=c(0, 1)) #
abline(v=metrics$support_size_true, col="red")
plot(x=support.size[sel], y=metrics$precision[sel], pch=16, type="o")
abline(v=metrics$support_size_true, col="red")
plot(x=support.size[sel], y=metrics$recall[sel], pch=16, type="o")
abline(v=metrics$support_size_true, col="red")
plot(x=support.size[sel], y=metrics$specificity[sel], pch=16, type="o")
abline(v=metrics$support_size_true, col="red")
plot(x=support.size[sel], y=metrics$accuracy[sel], pch=16, type="o")
support.size[which.max(metrics$accuracy)] # 12
abline(v=metrics$support_size_true, col="red")
abline(v=support.size[sel][which.max(as.vector(metrics$accuracy[sel]))], col="blue")
plot(x=support.size[sel], y=metrics$F1_score[sel], pch=16, type="o") 
support.size[which.max(metrics$F1_score)] # 12
abline(v=metrics$support_size_true, col="red")
abline(v=support.size[sel][which.max(as.vector(metrics$F1_score[sel]))], col="blue")
plot(x=support.size[sel], y=metrics$gic[sel], pch=16, type="o") # keep
support.size[which.min(metrics$gic)] # 12
abline(v=metrics$support_size_true, col="red")
plot(x=support.size[sel], y=metrics$aic[sel], pch=16, type="o") # keep
support.size[which.min(metrics$aic)] # 12
abline(v=metrics$support_size_true, col="red")
plot(x=support.size[sel], y=metrics$bic[sel], pch=16, type="o") # keep
support.size[which.min(metrics$bic)] # 12
abline(v=metrics$support_size_true, col="red")
plot(x=support.size[sel], y=metrics$ebic[sel], pch=16, type="o") # keep
support.size[which.min(metrics$ebic)] # 12
abline(v=metrics$support_size_true, col="red")