options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx12000m")) # increase JAVA memory to 12 Gb
# install.packages("pacman")
library(pacman)
# this command will load all required packages & will install them first if
# missing - you can place this command at the start of all your scripts
# this is more elegant than a lot of separate install.packages() and library() statements
# If you are asked "Would you like to use a personal library instead?" answer YES
suppressWarnings(pacman::p_load(rJava, glmnet, caret, mgcv, gam,
                                nnet, e1071, randomForest, xgboost,
                                Metrics, dplyr, effects, ggplot2, ggthemes, 
                                officer, flextable, export, abess,
                                update=FALSE))


## 1. Load covariate matrix with fingerprints and chemical properties calculated from molecular structures to predict outcome variable (aquatic toxicity) ####

# file with outcome variable (toxicity) and molecular structures (INCHI & SMILES)
DF = read.csv("../benchmarks/2_benchmark_vs_abess_relaxed_LASSO_toxicityQSAR/aquatictox_molecules_withsmiles.csv")
head(DF)
nrow(DF) # 322 molecules, toxicity=outcome variable, inchi & smiles encodes chemical structure

y = DF$toxicity # dependent variable = aquatic toxicity
names(y) =  paste0("Mol_", 1:length(y))

# file with calculated molecular fingerprints and chemical properties calculated
# from these molecular structures to use to predict toxicity and be able to
# predict the toxicity of other similar chemicals (see Part 1 tutorial to see
# how this was done) i.e. this is our predictor / covariate matrix
X = read.csv("../benchmarks/2_benchmark_vs_abess_relaxed_LASSO_toxicityQSAR/aquatictox_covariates.csv")


## 2. Covariate matrix preprocessing ####

# center and/or L2 norm normalize columns, to center all variables & get them
# all on the same scale 

# function to L2 norm normalize the columns of the design  matrix 
# (the L2 norm is the square root of the sum of the squared values)
# (i.e. to have them on the same scale, same aim as standardization) :
L2norm = function (X) apply(X, 2, function (col) col/norm(as.matrix(col),"2"))

# function to center & L2 norm normalize the columns of the design matrix :
centL2norm = function (X) apply(X, 2, function (col) (col-mean(col))/norm((col-mean(col)),"2"))

# here we will just L2 norm normalize all columns of X to get them all on the
# same scale (same function as standardisation)
X = L2norm(X=as.matrix(X)) # L2 norm normalize columns of X
Xc = centL2norm(X=as.matrix(X)) # centered & L2 norm normalized version of X

# for more help on preprocessing using caret, e.g. to include 
# one-hot encoded dummy-coded factors see 
# https://topepo.github.io/caret/pre-processing.html


## 3. Partition dataset into training & test set & filter out some variables ####

# we will use 75% of our data (p=0.75) for training the models and use the
# remaining 25% to test the performance
set.seed(1) # set random number generator seed to make script reproducible
data_index = createDataPartition(y, p = .75, list = FALSE)

dim(X) # 322 obs x 4397 variables

X_train = X[data_index, ] # training set covariate values
X_test = X[-data_index, ] # validation test set covariate values
Xc_train = Xc[data_index, ] # training set covariate values (also centering variables first)
Xc_test = Xc[-data_index, ] # validation test set covariate values (also centering variables first)


# we will now remove variables that are near-invariant in the training set :
X = X[,-nearZeroVar(X_train, freqCut = 99/1, uniqueCut = 1)]
dim(X) # 322 obs x 2776 variables
Xc = Xc[,-nearZeroVar(Xc_train, freqCut = 99/1, uniqueCut = 1)]
dim(Xc) # 322 obs x 2776 variables

# sometimes linearly dependent variables are removed using
X_nolindep = X[, -findLinearCombos(X_train)$remove]
dim(X_nolindep) # 322 240

# and sometimes highly collinear variables are also removed using
X_nocol = X_nolindep[,-findCorrelation(cor(X_nolindep), cutoff = 0.99, verbose = TRUE)]
dim(X_nocol) # 322 195

### here we will only remove the near-invariant variables as above, 
### so we continue working with X (2776 variables)

X_train = X[data_index, ] # training set covariate values
X_test = X[-data_index, ] # validation test set covariate values
Xc_train = Xc[data_index, ] # training set covariate values (also centering variables first)
Xc_test = Xc[-data_index, ] # validation test set covariate values (also centering variables first)

y_train = y[data_index] # training set outcome variable
y_test = y[-data_index] # validation test set outcome variable

n=nrow(Xc_train)
p=ncol(Xc_train)



## 1. RELAXED LASSO REGRESSION USING glmnet PACKAGE (L1 norm penalized regression) ####

# The L1 penalty in LASSO regression causes coefficients to have a downward bias.
# To reduce this bias, glmnet also offers the possibility to fit relaxed LASSO
# models, whereby the LASSO selected variables are refit unpenalised using a 
# regular least square regression or glm (if gamma is set to 0). With mixing
# hyperparameter gamma set to a value between 0 and 1 it will give a weighted
# average of the LASSO and relaxed LASSO solution.

# The optimal level of regularisation in relaxed LASSO comes out higher
# than in LASSO, which also reduces the nr of false positives.

# first determine optimal lambda regularization parameter using cross validation
set.seed(1)
system.time(relaxedlasso_cv <- cv.glmnet(x = X, y = y, 
                            family = "gaussian", # normally distributed noise, can also be "binomial", "poisson", "multinomial" or "cox" (for survival data)
                            alpha = 1, # to get LASSO penalty (ie with L1 norm penalty, 0<alpha<1 gives elastic net, with weighted L1 and L2 norm penalty, which can be better with many collinear variables)
                            type.measure = "deviance", # loss to use for cross-validation, "deviance" more general and also works for binomial & Poisson data
                            standardize = TRUE,
                            nfolds = 5, # we do 5-fold cross validation
                            gamma = c(0, 0.5), relax = TRUE, path = TRUE, # for relaxed LASSO, here using 2 different values of gamma
                            intercept = TRUE)) 
# 9.3s 
optlambda = relaxedlasso_cv$relaxed$lambda.min # optimal lambda that gives minimum out-of-sample MSE
optlambda # 0.07920818 - this is always higher than for regular LASSO, so model will result in fewer false positives
optgamma = relaxedlasso_cv$relaxed$gamma.min # optimal gamma value for relaxed LASSO (gamma=0=full relaxed LASSO, gamma=1=LASSO, in between: weighted avg)
optgamma # 0.5

plot(relaxedlasso_cv) # MSE ifo lambda sequence
# graph2ppt(file="relaxed LASSO_cv_lambda_MSE.pptx", width=4, height=4)

# now again fit relaxed LASSO model over sequence of lambda values (these are automatically set)
system.time(relaxedlassomodel <- glmnet(X_train, y_train, family = "gaussian", 
                           alpha = 1, # to get LASSO penalty
                           standardize = TRUE,
                           gamma = c(0, 0.5), relax = TRUE, path = TRUE, # for relaxed LASSO, here using 2 different values of gamma
                           intercept = TRUE)) 
# 1.3s
plot(relaxedlassomodel, xvar="lambda", gamma=optgamma) # coefficients ifo lambda for best gamma value
abline(v=log(optlambda), col="green4", lty=3, lwd=3) # optimal lambda
# graph2ppt(file="3_relaxed LASSO_coefficient paths.pptx", width=8, height=6)

beta = coef(relaxedlassomodel, s=optlambda, gamma=optgamma) # model coefficients at optimal lambda & gamma value determined using 5-fold cross validation
sum(beta!=0) # nr nonzero coefficients: 40 (39+intercept)
coefsrelaxedlasso = as.vector(beta)

# show results & validate on test set:
pred_test_relaxedlasso <- yhat_test <- predict(relaxedlassomodel, newx=X_test, s=optlambda, gamma=optgamma) # predictions test set
pred_train_relaxedlasso <- predict(relaxedlassomodel, newx=X_train, s=optlambda, gamma=optgamma) # predictions training set
relaxedlassoset = which(beta[-1]!=0)-1 # set of variables included in the model
rownames(beta)[relaxedlassoset] # variables included in the model

# benchmarks calculated on test set
benchm_relaxedlasso = c(RMSE = RMSE(yhat_test, y_test), # root mean squared error, sqrt(mean((pred - obs)^2))
                        Rsq = 1 - (sum((yhat_test - y_test) ^ 2))/(sum((y_test - mean(y_test)) ^ 2)), # R squared 
                        MAE = MAE(yhat_test, y_test), # mean absolute error
                        MdAE = median(abs(yhat_test-y_test)), # median absolute error
                        maxdev = max(abs(yhat_test-y_test))) # maximum deviation
benchm_relaxedlasso
#      RMSE  Rsquared       MAE      MdAE    maxdev 
# 0.4483345 0.8510729 0.3527127 0.2958129 1.4820455
plot(pred_train_relaxedlasso, y_train, 
     xlab="Predicted toxicity", ylab="Observed toxicity", pch=16, col="grey", xlim=c(0,7), ylim=c(0,7),
     main=paste0("relaxed LASSO regression (test R2=", as.character(round(benchm_relaxedlasso["Rsq"]*100,1)), "%)") )
points(pred_test_relaxedlasso, y_test, pch=16, col="blue" )
title("(grey=training set, blue=validation set, 2776 features, 39 features selected)", line=0.8, cex.main=0.9)
# graph2ppt(file="relaxed LASSO_obs vs pred.pptx", width=6, height=6)



## 2. BEST SUBSET SELECTION VIA abess PACKAGE ####

# preselected_vars = relaxedlassoset # to work with relaxed LASSO preselected variables
preselected_vars = 1:ncol(X_train) # to use all variables 
system.time(bestsubsetmodel <- abess(x = X_train[,preselected_vars], 
                        y = y_train, 
                        family = "gaussian",
                        tune.type = "gic", # use "gic" for best subset selection
                        fit.intercept = TRUE,
                        support.size = 0:(n-1), # min(n, round(n/(log(log(n))*log(p)))) = 18,
                        warm.start = FALSE # slower but gives better solution
)) 
# 0.06s if max support size is in(n, round(n/(log(log(n))*log(p)))) = 18
# 1s with max support size of n-1
plot(bestsubsetmodel, type="tune") # information criterion in function of support size (nr of selected variables / model size)
# graph2ppt(file="bestsubset_abess_cv_error_in_function_of_modelsize.pptx", width=4, height=4)

plot(bestsubsetmodel, type="coef") # coefficient path in functions of support size (nr of selected variables / model size)
# graph2ppt(file="bestsubset_abess_coefficient paths.pptx", width=8, height=6)

bestsubsetcoefs = as.vector(coef(bestsubsetmodel, 
                                 support.size = bestsubsetmodel[["best.size"]]))
sum(bestsubsetcoefs!=0) # 13=12+intercept

# show results & validate on test set:
coefs = bestsubsetcoefs 
pred_test_bestsubset <- yhat <- cbind(1, X_test[,preselected_vars]) %*% coefs 
pred_train_bestsubset <- cbind(1, X_train[,preselected_vars])  %*% coefs
abess_bestsubset_set = which(coefs[-1]!=0) # set of variables included in the model
colnames(X_train[,preselected_vars])[abess_bestsubset_set] # variables included in the model

# selected variables with coefficients
bestsubset_set_coefs = bestsubsetcoefs[-1][abess_bestsubset_set]
names(bestsubset_set_coefs) = colnames(X_train[,preselected_vars])[abess_bestsubset_set]
bestsubset_set_coefs

# benchmarks calculated on test set
benchm_bestsubset <- c(RMSE = RMSE(yhat, y_test), # root mean squared error, sqrt(mean((pred - obs)^2)) 
                       Rsquared = 1 - (sum((yhat - y_test) ^ 2))/(sum((y_test - mean(y_test)) ^ 2)), # R squared 
                       MAE = MAE(yhat, y_test), # mean absolute error
                       MdAE = median(abs(yhat-y_test)), # median absolute error
                       maxdev = max(abs(yhat-y_test))) # maximum deviation
benchm_bestsubset
#      RMSE  Rsquared       MAE      MdAE    maxdev 
# 0.5517107 0.7744763 0.4123050 0.3338744 1.9931948  
plot(pred_train_bestsubset, y_train, 
     xlab="Predicted toxicity", ylab="Observed toxicity", pch=16, col="grey", xlim=c(0,7), ylim=c(0,7),
     main=paste0("abess best subset (R2=", as.character(round(benchm_bestsubset["Rsquared"]*100,1)), "%)") )
points(pred_test_bestsubset, y_test, pch=16, col="blue" )
title("(grey=training set, blue=validation set, 2776 features, 12 features selected)", line=0.8, cex.main=0.9)
# graph2ppt(file="bestsubset_abess_obs vs pred.pptx", width=6, height=6)


# NOTE: abess currently does not report confidence intervals or p values
# for the selected variables; in principle those could be calculated using
# bootstrapping (resampling the original data with replacement & fitting the model
# on each resampled dataset, after which the 95% confidence intervals on the
# fitted coefficients would be given by the 2.5% and 97.5% percentiles of the
# coefficients obtained over all fits done on all resampled datasets).




## 3. APPROXIMATE BEST SUBSET SELECTION VIA L0glm PACKAGE ####

# NOTE: you can also approximately fit best subset glm models (L0 penalized GLM models)
# using my experimental L0glm package. This package uses an adapted iteratively reweighted least squares
# algorithm (IRLS) with adaptive ridge penalties to approximate best subset (cf. Frommlet & Nuel 2016).
# The level of penalisation (lambda) can be chosen a priori to approximately optimise an information
# criterion like AIC, BIC or GIC. The final version
# will appear on https://github.com/tomwenseleers/L0glm
# install.packages("L0glm_0.2.0.tar.gz")
library(L0glm)
preselected_vars_all = 1:ncol(X_train) # to use all variables
preselected_vars_sis = sis_glm( # select top n variables using sure independence screening
  X = Xc_train,
  y = y_train,
  family = gaussian(),
  nsis = n*2)$indices
preselected_vars_relaxedlasso = relaxedlassoset # relaxed LASSO preselected variables
preselected_vars_abess = abess_bestsubset_set # abess preselected variables

length(preselected_vars_sis) # n=484
length(preselected_vars_abess) # 12
length(preselected_vars_relaxedlasso) # 39

length(intersect(preselected_vars_abess, preselected_vars_relaxedlasso)) # 2
length(intersect(preselected_vars_abess, preselected_vars_sis)) # 6


# choose which set to use
# preselected_vars <- preselected_vars_abess # warm starting from the abess solution works, but that's defeating the purpose
# preselected_vars <- preselected_vars_sis # warm starting from the sis preselect set also still OKish
# preselected_vars <- preselected_vars_relaxedlasso # works badly
preselected_vars <- preselected_vars_all # using all variables


# # # CV based lambda tuning
# system.time(lambdatune <- L0glm.tune(X = Xc_train[,preselected_vars], # coding design matrix as sparse & also adding an intercept
#                         y = y_train,
#                         intercept = TRUE,
#                         nonnegative = FALSE,
#                         family = gaussian(identity),
#                         tune = "cv",
#                         nfolds = 5,
#                         nlambda = 10,
#                         lambda_min = 0.1,
#                         solver = "eigen",
#                         method = 1L
#                         )) 
# bestlambda_cv = lambdatune$lambda_best
# bestlambda_cv 

# GIC based lambda tuning
system.time(lambdatune <- L0glm.tune(X = Xc_train[,preselected_vars], # coding design matrix as sparse & also adding an intercept
                        y = y_train, 
                        intercept = TRUE,
                        nonnegative = FALSE,
                        family = gaussian(identity),
                        tune = "gic",
                        nlambda = 10,
                        lambda_min = 2,
                        solver = "eigen", 
                        method = 1L
)) # 188s
bestlambda_gic = lambdatune$lambda_best
bestlambda_gic # 13.5012
lambdatune$fit_best$dispersion # 0.3236964

system.time(bestsubsetcoefs2 <- L0glm.fit(X = Xc_train[,preselected_vars], # coding design matrix as sparse & also adding an intercept
                                                             y = y_train, 
                                                             intercept = TRUE,
                                                             family = gaussian(identity),
                                                             lambda = "gic", # bestlambda_gic also works, "aic" doesn't
                                                             nonnegative = FALSE,
                                                             start = 1, # this one works!, start=NULL does not
                                                             solver = "eigen" # using the glmnet::bigGlm solver internally
                                )$coefficients) # 17s
sum(bestsubsetcoefs2!=0) # 15
bestsubset_set2 = which(bestsubsetcoefs2[-1]!=0) # set of variables included in the model
bestsubsetcoefs2[-1][bestsubset_set2]
colnames(X_train[,preselected_vars])[bestsubset_set2] # 14 variables
# selected variables with coefficients
bestsubset_set_coefs2 = bestsubsetcoefs2[-1][bestsubset_set2]
names(bestsubset_set_coefs2) = colnames(X_train[,preselected_vars])[bestsubset_set2]
bestsubset_set_coefs2


# show results & validate on test set:
coefs = bestsubsetcoefs2 
pred_test_bestsubset <- yhat <- cbind(1, Xc_test[,preselected_vars])  %*% coefs
pred_train_bestsubset <- yhat_train <- cbind(1, Xc_train[,preselected_vars])  %*% coefs
bestsubset_set2 = which(coefs[-1]!=0) # set of variables included in the model
colnames(X_train[,preselected_vars])[bestsubset_set2] # variables included in the model

# benchmarks calculated on test set
benchm_bestsubset <- c(RMSE = RMSE(as.vector(yhat), y_test), # root mean squared error, sqrt(mean((pred - obs)^2)) 
                       Rsquared = 1 - (sum((as.vector(yhat) - y_test) ^ 2))/(sum((y_test - mean(y_test)) ^ 2)), # R squared 
                       sqcor = cor(as.vector(yhat), y_test)^2, # squared correlation
                       MAE = MAE(as.vector(yhat), y_test), # mean absolute error
                       MdAE = median(abs(yhat-y_test)), # median absolute error
                       maxdev = max(abs(yhat-y_test))) # maximum deviation
benchm_bestsubset
#      RMSE  Rsquared     sqcor       MAE      MdAE    maxdev 
# 0.5010755 0.8139731 0.8354221 0.3568975 0.2884061 1.9826291  
plot(pred_train_bestsubset, y_train, 
     xlab="Predicted toxicity", ylab="Observed toxicity", pch=16, col="grey", xlim=c(0,7), ylim=c(0,7),
     main=paste0("L0glm best subset (R2=", as.character(round(benchm_bestsubset["Rsquared"]*100,1)), "%)") )
points(pred_test_bestsubset, y_test, pch=16, col="steelblue" )
title("(grey=training set, blue=validation set, 2757 features, 14 features selected)", line=0.8, cex.main=0.9)
# graph2ppt(file="bestsubset_L0glm_obs vs pred.pptx", width=6, height=6)
