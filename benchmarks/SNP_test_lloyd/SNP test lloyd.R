load(file.path(path, "envSmall2 1.RData"))
path = "~/Github/L0glm/benchmarks/SNP_test_lloyd"

library(L0glm)
library(ordinis)
library(ncvreg)
library(glmnet)
library(caret)
library(ggplot2)
library(ggthemes)
library(scales)
library(export)

# A. SET UP SNP MATRIX ####

remove = which(min(y)==y) # remove 1 observation with y=0.001 (original obs was 0 lag?)
Y = y[-remove]
X = reduced_dataset[-remove,]

# we recode the reference strain to be the strain with the lowest lag time,
# so that we can use nonnegativity constraints in our fitted model
new_refstrain = which.min(Y) # 14th row is strain with lowest lag time (min(Y)=0.4128311) (not sure which strain this is, there is no rownames in reduced_dataset)
# we code this strain with the lowest trait value (lag time) as our reference strain, 
# so that all SNPs are defined relative to this & coefficients in the model can
# be made to be all positive. Nonnegativity constraints will help to regularize
# the problem & will increase the statistical power to pick up causal SNPs.
# We reverse the coding for all SNPs coded as 2 in that strain, so that all
# SNPs will be defined relative to that.
columns_to_be_flipped = which(X[new_refstrain,]==2)
X[,columns_to_be_flipped] = 2-X[,columns_to_be_flipped]

# Y = y
# X = reduced_dataset
n = nrow(X)
n # 99
p = ncol(X)
p # 67397

# log (y+a) transform outcome variable
a = 0.01 # NOTE: a should be chosen so that variance around predicted mean is approx constant (i.e. that variance does not increase in function of predicted mean, with 0.01 that was the case)
logY = log(Y+a)
hist(logY, breaks=20) # with this value of a distribution looks fairly symmetrical
minv = log(min(Y)+a)
maxv = log(max(Y)+a)

# split data in training & test set using 60%:40% split using stratified sampling using deciles of logY
set.seed(1)
prop_train = 0.75
train_indices = createDataPartition(y=logY, p=prop_train, groups=10, list=F)
X_train = X[train_indices,]
Y_train = Y[train_indices]
X_test = X[-train_indices,]
Y_test = Y[-train_indices]
n_train = nrow(X_train)
n_train # 81
n_test = nrow(X_test)
n_test # 18

# Log-transformed variables for training and test sets
logY_train = log(Y_train + a)
hist(logY_train, breaks=20)
length(logY_train) # 81
logY_test = log(Y_test + a)
length(logY_test) # 18
hist(logY_test, breaks=20)

# sum(colSums(reduced_dataset==1)) # 1742


# B. TAKE INTO ACCOUNT GENETIC BACKGROUND / ANCESTRY ####

# To take into account ancestry, there is 2 main options:

## PHYLOGENETIC MIXED MODELS (pglmm) OR PHYLOGENETIC GENERALISED LEAST SQUARES MODELS (pgls) ####

# if we would happen to have a nexus file of the phylogeny, you could use
# phylo <- ape::read.nexus("XXX.nex")
# and calculate a strain phylogenetic variance-covariance matrix
# Z <- ape::vcv.phylo(phylo)
# which you could then take into account in a phylogenetic mixed model (pglmmm) or
# phylogenetic generalised least squares (pgls) fit, after which you could
# work with the residuals of that model.
# E.g. using
# - pglmm in phyr
# - Almer  in evolvability
# - lme4 through this hack: 
# Brain dump: phylogenetic comparative methods in the lme4-verse (bbolker.github.io)
# - brms (Bayesian): Estimating Phylogenetic Multilevel Models with brms (r-project.org)
# https://cran.r-project.org/web/packages/brms/vignettes/brms_phylogenetics.html
# - MCMCglmm (Bayesian): Chapter 11 – 1) A simple model in MCMCglmm - Book MPCM Evolution 2014 (mpcm-evolution.com)
# Online Practical Materials Chapter 11 – General Quantitative Genetic Methods for Comparative Biology 1) A simple model in MCMCglmm ← back to the list of Chapter 11 OPMs
# www.mpcm-evolution.com
# - lme, I think syntax is:
#  group <- ... # Grouping factor corresponding to rows/cols of Z
#
## Create a correlation structure using the phylogenetic vcov matrix Z
#phyloCorr <- corSymm(value = Z, form = ~ 1 | group)
#
# Fit the lme model
# model <- lme(fixed = response ~ fixedEffect,
#             data = data,
#             random = list(group = pdBlocked(list(pdIdent(~ 1 | group)))),
#             correlation = phyloCorr)
#
#- work with residuals of phylogenetic generalised least squares model (pgls model): Phylogenetic generalized linear models in caper


## ADD PRINCIPAL COMPONENTS OF SNP MATRIX AS UNPENALISED FIXED EFFECTS ####

# Simpler alternative is to use some leading principal components of the SNP
# matrix for Z and include these as unpenalised fixed effects in the model.
# (via penalty.factor in glmnet & ordinis, always.include in abess & no.pen in L0glm).

# Here I'll go for this simpler option:
# Calculate PCA on X_train
pca_res <- prcomp(X_train, scale. = TRUE)
screeplot(pca_res, type="lines")
num_pcs <- 2 # 10 # Number of PCA factors to use
# Get PCA factors for training data
pca_factors_train <- pca_res$x[, 1:num_pcs]
# Project test data & full X matrix data using training PCA
pca_factors_test <- predict(pca_res, newdata = X_test)[, 1:num_pcs]
pca_factors_X <- predict(pca_res, newdata = X)[, 1:num_pcs]

# add PCA factors to training & test data
X2 = cbind(pca_factors_X, X)
X_train2 = cbind(pca_factors_train, X_train)
X_test2 = cbind(pca_factors_test, X_test)

# check predictive performance of pure principal component regression
pcreg = lm(logY_train ~ ., data=data.frame(logY_train, pca_factors_train))
summary(pcreg)
coefs_pcreg = coef(pcreg)

preds_pcreg_train = exp(cbind(1, pca_factors_train) %*% coefs_pcreg) - a
preds_pcreg_test = exp(cbind(1, pca_factors_test) %*% coefs_pcreg) - a

plot(log(preds_pcreg_train+a), log(Y_train+a), 
     pch=16, col="blue",
     xlim=c(minv,maxv), ylim=c(minv,maxv))
points(log(preds_pcreg_test+a), log(Y_test+a), 
       pch=16, col="red")
plot(preds_pcreg_train, Y_train, pch=16, col="blue")
points(preds_pcreg_test, Y_test, pch=16, col="red")

# predictive performance on test set
postResample(pred = log(Y_test+a), 
             obs = log(preds_pcreg_test+a))

# predictive performance on training set
postResample(pred = log(Y_train+a), 
             obs = log(preds_pcreg_train+a))


# Get the loadings
loadings <- pca_res$rotation

# Determine the number of top SNPs and components you're interested in
num_top_snps <- 10

# Initialize a list to store SNP names
top_snps_by_component <- list()

# Loop over the first 'num_components' principal components
for (i in 1:num_pcs) {
  # Get indices of top 'num_top_snps' SNPs based on loading absolute values
  top_snps_indices <- order(-abs(loadings[, i]))[1:num_top_snps]
  
  # Retrieve the SNP names using the indices
  top_snps_names <- colnames(X_train)[top_snps_indices]
  
  # Store the names in the list
  top_snps_by_component[[paste("PC", i, sep = "")]] <- top_snps_names
}

# The list 'top_snps_by_component' now contains the SNP names for each component
top_snps_by_component

# NOTE: only downside of including these PCA components is that we can then
# no longer use nonnegativity constraints on the fitted coefficients
# If you would like to retain those then we would have to use e.g.
# a generalised low rank approximation of our SNP matrix with nonnegativity
# constraints & assuming Poisson noise, e.g. using
# https://rdrr.io/cran/h2o/man/h2o.glrm.html
# or an NMF, e.g. using RcppML::nmf
# Haven't tried this yet.
# Especially L0glm tends to perform much better with nonnegativity constraints.
# (a nonnegative double SVD, https://rdrr.io/github/zdebruine/LSMF/man/nndsvd.html
# could be another option, though a generalised low rank model is probably more
# principled - not sure if distribution of SNPs would best be considered as
# Poisson distributed or as binomial distributed though)



# C. SET UP ADAPTIVE PENALTIES ####

# The ROR1 SNPs that you already known to be involved in determining lag.

ROR1_SNP_names = c('chromosome12_366378','chromosome12_366845',
                   'chromosome12_367241','chromosome12_367293',
                   'chromosome12_367545')
ROR1_SNP = (colnames(X_train) %in% ROR1_SNP_names) # leave these SNPs unpenalised in the models
ROR1_SNP_IDs = which(ROR1_SNP)
colnames(X_train)[ROR1_SNP_IDs]
# "chromosome12_366378" "chromosome12_367241" "chromosome12_367293" "chromosome12_367545"
length(ROR1_SNP_IDs) # 4
# NOTE: chromosome12_366845 was collinear with chromosome12_366378 & filtered out beforehand 

# These ROR1 SNPs and the PCA factors we leave unpenalised in the models.
penalty.factor = c(rep(0,num_pcs),(!ROR1_SNP)*1) # to leave PCA factors & ROR1 SNPs unpenalised
penalty.factor2 = c(rep(1,num_pcs),(!ROR1_SNP)*1) # to only leave ROR1 SNPs unpenalised



# 1.1 LASSO fit using glmnet (GOOD) ####

library(glmnet)

set.seed(1)
# determine optimal lambda
alpha = 1 # 1, ie LASSO, gave best performance 
lassocvfit = cv.glmnet(x=X2, 
                       y=logY, 
                       intercept = T,
                       # lower.limits = c(rep(-Inf,num_pcs), rep(0,ncol(X))), # for nonnegativity constraints (possible given that we coded strain with lowest trait value as reference strain)
                       # upper.limits = c(rep(Inf,num_pcs), rep(Inf,ncol(X))), # for nonnegativity constraints (possible given that we coded strain with lowest trait value as reference strain)
                       penalty.factor = penalty.factor, # do not penalize PCA factors & ROR1 SNPs
                       nlambda = 100L, 
                       lambda.min.ratio = 0.01,
                       dfmax = n,
                       alpha = alpha,
                       nfolds = 10)
bestlambda_lasso = lassocvfit$lambda.min
plot(lassocvfit)

lassofit = glmnet(x=X_train2, 
                  y=logY_train, 
                  intercept = T,
                  # lower.limits = c(rep(-Inf,num_pcs), rep(0,ncol(X))), # for nonnegativity constraints (possible given that we coded strain with lowest trait value as reference strain)
                  # upper.limits = c(rep(Inf,num_pcs), rep(Inf,ncol(X))), # for nonnegativity constraints (possible given that we coded strain with lowest trait value as reference strain)
                  penalty.factor = penalty.factor, # do not penalize PCA factors & RO1 SNPs
                  nlambda = 100L, 
                  lambda.min.ratio = 0.01,
                  dfmax = n,
                  alpha = alpha)
plot(lassofit, xvar="lambda")
lassocoefs = as.vector(coef(lassofit, s=bestlambda_lasso))
names(lassocoefs) = c("Intercept", colnames(X2))
lassoset = which(lassocoefs[-1]!=0)
range(lassocoefs[-c(1:num_pcs)])
sum(lassocoefs!=0) # 18

lassocoefs[which(lassocoefs!=0)]
#    Intercept                 PC1                 PC2  chromosome2_647503 
# 7.090706e-01       -6.887220e-05        1.852857e-03        2.252114e-02 
# chromosome3_98952   chromosome3_98969  chromosome4_589643 chromosome10_275952 
#      3.660608e-02        7.963131e-05        9.171777e-02        3.138695e-02 
# chromosome12_366378 chromosome12_367241 chromosome12_367293 chromosome12_367545 
#        1.916885e-01       -3.234361e-02       -3.086797e-01        2.089371e-01 
# chromosome15_82213  chromosome15_91150 chromosome15_533398 chromosome15_533423 
#       3.230805e-04        1.204716e-02        5.299003e-02        9.289180e-16 
# chromosome16_500828 chromosome16_610743 
#        1.096702e-01        4.378103e-02

lassopreds_train = exp(cbind(1,X_train2) %*% lassocoefs)-a
lassopreds_test = exp(cbind(1,X_test2) %*% lassocoefs)-a

# predictive performance on test set
predperf_test = data.frame(t(postResample(pred = log(Y_test+a), 
                                          obs = log(lassopreds_test+a))))
predperf_test
#      RMSE  Rsquared       MAE 
# 0.605987 0.4769911 0.4553376

# predictive performance on training set
predperf_train = data.frame(t(postResample(pred = log(Y_train+a), 
                                           obs = log(lassopreds_train+a))))
predperf_train
#      RMSE  Rsquared       MAE 
# 0.5160314 0.4587799 0.3932061

model="LASSO model"
test_R2_rounded <- round(predperf_test$Rsquared * 100, 0)
test_RMSE_rounded <- round(predperf_test$RMSE, 2)
train_R2_rounded <- round(predperf_train$Rsquared * 100, 0)
train_RMSE_rounded <- round(predperf_train$RMSE, 2)
subtitle_test <- bquote("test " ~ R^2 ~ " = " ~ .(test_R2_rounded) ~ "%, test RMSE = " ~ .(test_RMSE_rounded) ~ "  ")
subtitle_train <- bquote("training " ~ R^2 ~ " = " ~ .(train_R2_rounded) ~ "%, training RMSE = " ~ .(train_RMSE_rounded) ~ "  ")
data_plot <- data.frame(
  Predicted = log(c(lassopreds_train, lassopreds_test)  + a),
  Observed = log(c(Y_train, Y_test) + a),
  Group = relevel(factor(c(rep("training set", length(Y_train)),
                   rep("test set", length(Y_test)))), ref="training set")
)
ggplot(data_plot, aes(x = Predicted, y = Observed, color = Group)) +
  geom_point(shape = 16, size=2) +
  scale_x_continuous("Predicted lag time",
                     breaks = log(c(1:10, 12, 15, 20) + a),
                     labels = c(1:10, 12, 15, 20),
                     limits = c(log(0.5 + a), log(25 + a))) +
  scale_y_continuous("Observed lag time",
                     breaks = log(c(1:10, 12, 15, 20) + a),
                     labels = c(1:10, 12, 15, 20),
                     limits = c(log(0.5 + a), log(25 + a))) +
  scale_color_manual(values=c("blue", "red")) +
  labs(title=model) +
  annotate("text", x = Inf, y = Inf, label = subtitle_test, hjust = 1, vjust = 1.5, size = 3.5) +
  annotate("text", x = Inf, y = Inf, label = subtitle_train, hjust = 1, vjust = 2.5, size = 3.5) +
  theme_few(base_size=14) +
  theme(legend.position="none")
graph2ppt(file=file.path(path, "1_1_LASSO model.pptx"), width=5, height=4)


# 1.2 relaxed LASSO fit using glmnet (GOOD, 2nd best test RMSE) ####

library(glmnet)
set.seed(1)
# determine optimal lambda
relaxedlassocvfit = cv.glmnet(x=X2,
                              y=logY,
                              intercept = T,
                              # lower.limits = c(rep(-Inf,num_pcs), rep(0,ncol(X))), # for nonnegativity constraints (possible given that we coded strain with lowest trait value as reference strain)
                              # upper.limits = c(rep(Inf,num_pcs), rep(Inf,ncol(X))), # for nonnegativity constraints (possible given that we coded strain with lowest trait value as reference strain)
                              penalty.factor = penalty.factor, # do not penalize PCA factors & ROR1 SNPs
                              nlambda = 100L,
                              lambda.min.ratio = 0.01,
                              dfmax = n,
                              alpha = 1, # LASSO
                              relax = TRUE, path = TRUE,
                              gamma = c(0, 0.25, 0.5, 0.75, 1), # for relaxed LASSO
                              nfolds = 10)
bestlambda_relaxedlasso = relaxedlassocvfit$relaxed$lambda.min
bestlambda_relaxedlasso
bestgamma_relaxedlasso = relaxedlassocvfit$relaxed$gamma.min
bestgamma_relaxedlasso # 0.5 -> = ie intermediate between pure LASSO & fully relaxed LASSO

plot(relaxedlassocvfit)

relaxedlassofit = glmnet(x=X_train2,
                         y=logY_train,
                         intercept = T,
                         # lower.limits = c(rep(-Inf,num_pcs), rep(0,ncol(X))), # for nonnegativity constraints (possible given that we coded strain with lowest trait value as reference strain)
                         # upper.limits = c(rep(Inf,num_pcs), rep(Inf,ncol(X))), # for nonnegativity constraints (possible given that we coded strain with lowest trait value as reference strain)
                         penalty.factor = penalty.factor, # do not penalize PCA factors & ROR1 SNPs
                         nlambda = 100L,
                         lambda.min.ratio = 0.01,
                         dfmax = n,
                         alpha = 1, # LASSO
                         relax = TRUE, path = TRUE,
                         gamma = c(0, 0.25, 0.5, 0.75, 1))
plot(relaxedlassofit, xvar="lambda", gamma=bestgamma_relaxedlasso)
relaxedlassocoefs = as.vector(coef(relaxedlassofit,
                         s=bestlambda_relaxedlasso,
                         gamma=bestgamma_relaxedlasso))
names(relaxedlassocoefs) = c("Intercept", colnames(X2))

relaxedlassoset = which(relaxedlassocoefs[-1]!=0)
sum(relaxedlassocoefs!=0) # 10
range(relaxedlassocoefs)

relaxedlassocoefs[which(relaxedlassocoefs!=0)]
#   Intercept                 PC1                 PC2  chromosome4_589643 
# 0.775078092         0.000111360         0.001428462         0.219140428 
# chromosome12_366378 chromosome12_367241 chromosome12_367293 chromosome12_367545 
#         0.236316883        -0.071868877        -0.298411668         0.217066606 
# chromosome16_500828 chromosome16_610743 
#         0.215016126         0.108752702 

relaxedlassopreds_train = exp(cbind(1,X_train2) %*% relaxedlassocoefs)-a
relaxedlassopreds_test = exp(cbind(1,X_test2) %*% relaxedlassocoefs)-a

# predictive performance on test set
predperf_test = data.frame(t(postResample(pred = log(Y_test+a), 
                                          obs = log(relaxedlassopreds_test+a))))
predperf_test
#      RMSE  Rsquared       MAE 
# 0.5755353 0.5036339 0.4457884

# predictive performance on training set
predperf_train = data.frame(t(postResample(pred = log(Y_train+a), 
                                           obs = log(relaxedlassopreds_train+a))))
predperf_train
#      RMSE  Rsquared       MAE 
# 0.5046236 0.4300272 0.3912162

model="relaxed LASSO model"
test_R2_rounded <- round(predperf_test$Rsquared * 100, 0)
test_RMSE_rounded <- round(predperf_test$RMSE, 2)
train_R2_rounded <- round(predperf_train$Rsquared * 100, 0)
train_RMSE_rounded <- round(predperf_train$RMSE, 2)
subtitle_test <- bquote("test " ~ R^2 ~ " = " ~ .(test_R2_rounded) ~ "%, test RMSE = " ~ .(test_RMSE_rounded) ~ "  ")
subtitle_train <- bquote("training " ~ R^2 ~ " = " ~ .(train_R2_rounded) ~ "%, training RMSE = " ~ .(train_RMSE_rounded) ~ "  ")
data_plot <- data.frame(
  Predicted = log(c(relaxedlassopreds_train, relaxedlassopreds_test)  + a),
  Observed = log(c(Y_train, Y_test) + a),
  Group = relevel(factor(c(rep("training set", length(Y_train)),
                           rep("test set", length(Y_test)))), ref="training set")
)
ggplot(data_plot, aes(x = Predicted, y = Observed, color = Group)) +
  geom_point(shape = 16, size=2) +
  scale_x_continuous("Predicted lag time",
                     breaks = log(c(1:10, 12, 15, 20) + a),
                     labels = c(1:10, 12, 15, 20),
                     limits = c(log(0.5 + a), log(25 + a))) +
  scale_y_continuous("Observed lag time",
                     breaks = log(c(1:10, 12, 15, 20) + a),
                     labels = c(1:10, 12, 15, 20),
                     limits = c(log(0.5 + a), log(25 + a))) +
  scale_color_manual(values=c("blue", "red")) +
  labs(title=model) +
  annotate("text", x = Inf, y = Inf, label = subtitle_test, hjust = 1, vjust = 1.5, size = 3.5) +
  annotate("text", x = Inf, y = Inf, label = subtitle_train, hjust = 1, vjust = 2.5, size = 3.5) +
  theme_few(base_size=14) +
  theme(legend.position="none")
graph2ppt(file=file.path(path, "1_2_relaxed LASSO model.pptx"), width=5, height=4)


# 1.3 SCAD fit using ordinis (GOOD) ####

library(ordinis)

set.seed(1)
# determine optimal lambda
cvfit = cv.ordinis(x=X2,
                   y=logY,
                   intercept = T,
                   penalty = "scad",  # SCAD worked better than MCP
                   # gamma = 3.7, # 3.7 for SCAD, 3 for MCP default
                   # lower.limits = c(rep(-Inf,num_pcs), rep(0,ncol(X))), # for nonnegativity constraints (possible given that we coded strain with lowest trait value as reference strain)
                   # upper.limits = c(rep(Inf,num_pcs), rep(Inf,ncol(X))), # for nonnegativity constraints (possible given that we coded strain with lowest trait value as reference strain)
                   penalty.factor = penalty.factor, # do not penalize PCA factors & ROR1 SNPs
                   nlambda = 100L,
                   # lambda.min.ratio = 0.01,
                   dfmax = n,
                   alpha = 1,
                   nfolds = 10)
plot(cvfit)
bestlambda = cvfit$lambda.min

ordfit = ordinis(x=X_train2,
                 y=logY_train,
                 intercept = T,
                 penalty = "scad",
                 # gamma = 3.7 for SCAD, 3 for MCP default
                 # lower.limits = c(rep(-Inf,num_pcs), rep(0,ncol(X))), # for nonnegativity constraints (possible given that we coded strain with lowest trait value as reference strain)
                 # upper.limits = c(rep(Inf,num_pcs), rep(Inf,ncol(X))), # for nonnegativity constraints (possible given that we coded strain with lowest trait value as reference strain)
                 penalty.factor = penalty.factor, # do not penalize PCA factors & ROR1 SNPs
                 nlambda = 100L,
                 # lambda.min.ratio = 0.001,
                 dfmax = n,
                 alpha = 1)
plot(ordfit)
ordcoefs = ordfit$beta[,which.min(abs(bestlambda-ordfit$lambda))]
ordset = which(ordcoefs[-1]!=0)
sum(ordcoefs!=0) # 20

ordcoefs[which(ordcoefs!=0)]
#  (Intercept)                 PC1                 PC2  chromosome2_647503 
# 7.557267e-01       -8.805448e-05        1.792761e-03        1.322597e-02 
# chromosome3_98952   chromosome3_98969   chromosome3_99037   chromosome3_99197 
#      2.226103e-02        3.339154e-03        5.008731e-04        7.513097e-05 
# chromosome3_99244  chromosome4_589643 chromosome10_275952 chromosome12_366378 
#      1.126965e-05        8.372214e-02        2.539921e-02        1.968820e-01 
# chromosome12_367241 chromosome12_367293 chromosome12_367545  chromosome15_91150 
#       -2.969407e-02       -3.134395e-01        2.067862e-01        8.116644e-03 
# chromosome15_533398 chromosome15_533423 chromosome16_500828 chromosome16_610743 
#        3.652886e-02        5.479329e-03        9.924741e-02        4.217771e-02

range(ordcoefs)
ordpreds_train = exp(cbind(1,X_train2) %*% ordcoefs)-a
ordpreds_test = exp(cbind(1,X_test2) %*% ordcoefs)-a

# predictive performance on test set
predperf_test = data.frame(t(postResample(pred = log(Y_test+a), 
                                          obs = log(ordpreds_test+a))))
predperf_test
#      RMSE  Rsquared       MAE 
# 0.6088008 0.4911842 0.4572762

# predictive performance on training set
predperf_train = data.frame(t(postResample(pred = log(Y_train+a), 
                                           obs = log(ordpreds_train+a))))
predperf_train
#      RMSE  Rsquared       MAE 
# 0.5310102 0.407099 0.403658

model="SCAD model"
test_R2_rounded <- round(predperf_test$Rsquared * 100, 0)
test_RMSE_rounded <- round(predperf_test$RMSE, 2)
train_R2_rounded <- round(predperf_train$Rsquared * 100, 0)
train_RMSE_rounded <- round(predperf_train$RMSE, 2)
subtitle_test <- bquote("test " ~ R^2 ~ " = " ~ .(test_R2_rounded) ~ "%, test RMSE = " ~ .(test_RMSE_rounded) ~ "  ")
subtitle_train <- bquote("training " ~ R^2 ~ " = " ~ .(train_R2_rounded) ~ "%, training RMSE = " ~ .(train_RMSE_rounded) ~ "  ")
data_plot <- data.frame(
  Predicted = log(c(ordpreds_train, ordpreds_test)  + a),
  Observed = log(c(Y_train, Y_test) + a),
  Group = relevel(factor(c(rep("training set", length(Y_train)),
                           rep("test set", length(Y_test)))), ref="training set")
)
ggplot(data_plot, aes(x = Predicted, y = Observed, color = Group)) +
  geom_point(shape = 16, size=2) +
  scale_x_continuous("Predicted lag time",
                     breaks = log(c(1:10, 12, 15, 20) + a),
                     labels = c(1:10, 12, 15, 20),
                     limits = c(log(0.5 + a), log(25 + a))) +
  scale_y_continuous("Observed lag time",
                     breaks = log(c(1:10, 12, 15, 20) + a),
                     labels = c(1:10, 12, 15, 20),
                     limits = c(log(0.5 + a), log(25 + a))) +
  scale_color_manual(values=c("blue", "red")) +
  labs(title=model) +
  annotate("text", x = Inf, y = Inf, label = subtitle_test, hjust = 1, vjust = 1.5, size = 3.5) +
  annotate("text", x = Inf, y = Inf, label = subtitle_train, hjust = 1, vjust = 2.5, size = 3.5) +
  theme_few(base_size=14) +
  theme(legend.position="none")
graph2ppt(file=file.path(path, "1_3_SCAD model.pptx"), width=5, height=4)



# # 1.3B SCAD or MCP fit using ncvreg ####
# 
# library(ncvreg)
# 
# # determine optimal lambda
# cvfit = cv.ncvreg(X=X, 
#                   y=logY, 
#                   intercept = T,
#                   penalty = "SCAD", 
#                   # gamma = 3, # 3.7 for SCAD, 3 for MCP default
#                   nlambda = 100L, 
#                   lambda.min = 0.05,
#                   dfmax = n,
#                   alpha = 1,
#                   nfolds = 5)
# bestlambda = cvfit$lambda.min
# plot(cvfit)



# 1.4 abess fit (GOOD, best test set RMSE) ####

library(abess)
preselected_vars = relaxedlassoset # ordset # lassoset # using relaxed LASSO preselection of variables, lassoset also works but fit worse
# preselected_vars = 1:ncol(X2) # to use all variables
X_train_subs2 = X_train2[,preselected_vars] # only main effects
X_test_subs2 = X_test2[,preselected_vars] # only main effects
# also tried adding all 1st order interactions, but fit worse
# X_train_subs2 = model.matrix(~0+.^2, # all main effects+1st order interaction effects
#                             data=data.frame(X_train[,preselected_vars])) 
# X_test_subs2 = model.matrix(~0+.^2, # all main effects+1st order interaction effects
#                            data=data.frame(X_test[,preselected_vars])) 

abessfit = abess(x=X_train_subs2, 
                 y=logY_train,
                 # beta.low = 0, # for nonnegativity constraints (possible given that we coded strain with lowest trait value as reference strain)
                 always.include = which(penalty.factor[preselected_vars]==0), # always include ROR1 SNPs (if you also always include PCA factors model was worse)
                 fit.intercept = T, 
                 tune.path="sequence",
                 tune.type="aic", # you can use aic if you use a preselected set, otherwise use gic or bic
                 normalize=F,
                 # nfolds=10,
                 max.splicing.iter = 100,
                 warm.start = T #,
                 #support.size = 1:round(n/2)
                 )
bestsupportsize = abessfit[["best.size"]] 
bestsupportsize # 9

abesscoefs = as.vector(coef(abessfit, 
                            support.size = bestsupportsize))
names(abesscoefs) = c("Intercept", colnames(X_train_subs2))
abessset = which(abesscoefs[-1]!=0)
abesscoefs[which(abesscoefs!=0)]
#    Intercept                 PC1                 PC2  chromosome4_589643 
# 0.6550205477        0.0003655571        0.0012095673        0.4074851997 
# chromosome12_366378 chromosome12_367241 chromosome12_367293 chromosome12_367545 
#        0.2900666797       -0.1604696426       -0.2551567153        0.2380030727 
# chromosome16_500828 chromosome16_610743 
#        0.3857281160        0.1932357216

range(abesscoefs) 
# abesspreds_train = exp(cbind(1,X_train) %*% abesscoefs)-a
# abesspreds_test = exp(cbind(1,X_test) %*% abesscoefs)-a
abesspreds_train = exp(cbind(1,X_train_subs2) %*% abesscoefs)-a
abesspreds_test = exp(cbind(1,X_test_subs2) %*% abesscoefs)-a

# predictive performance on test set
predperf_test = data.frame(t(postResample(pred = log(Y_test+a), 
                                          obs = log(abesspreds_test+a))))
predperf_test
#      RMSE  Rsquared       MAE 
# 0.5638767 0.4204861 0.4421836

# predictive performance on training set
predperf_train = data.frame(t(postResample(pred = log(Y_train+a), 
                                           obs = log(abesspreds_train+a))))
predperf_train
#      RMSE  Rsquared       MAE 
# 0.4722008 0.4620488 0.3624877

model="abess model on relaxed LASSO preselected variables"
test_R2_rounded <- round(predperf_test$Rsquared * 100, 0)
test_RMSE_rounded <- round(predperf_test$RMSE, 2)
train_R2_rounded <- round(predperf_train$Rsquared * 100, 0)
train_RMSE_rounded <- round(predperf_train$RMSE, 2)
subtitle_test <- bquote("test " ~ R^2 ~ " = " ~ .(test_R2_rounded) ~ "%, test RMSE = " ~ .(test_RMSE_rounded) ~ "  ")
subtitle_train <- bquote("training " ~ R^2 ~ " = " ~ .(train_R2_rounded) ~ "%, training RMSE = " ~ .(train_RMSE_rounded) ~ "  ")
data_plot <- data.frame(
  Predicted = log(c(abesspreds_train, abesspreds_test)  + a),
  Observed = log(c(Y_train, Y_test) + a),
  Group = relevel(factor(c(rep("training set", length(Y_train)),
                           rep("test set", length(Y_test)))), ref="training set")
)
ggplot(data_plot, aes(x = Predicted, y = Observed, color = Group)) +
  geom_point(shape = 16, size=2) +
  scale_x_continuous("Predicted lag time",
                     breaks = log(c(1:10, 12, 15, 20) + a),
                     labels = c(1:10, 12, 15, 20),
                     limits = c(log(0.5 + a), log(25 + a))) +
  scale_y_continuous("Observed lag time",
                     breaks = log(c(1:10, 12, 15, 20) + a),
                     labels = c(1:10, 12, 15, 20),
                     limits = c(log(0.5 + a), log(25 + a))) +
  scale_color_manual(values=c("blue", "red")) +
  labs(title=model) +
  annotate("text", x = Inf, y = Inf, label = subtitle_test, hjust = 1, vjust = 1.5, size = 3.5) +
  annotate("text", x = Inf, y = Inf, label = subtitle_train, hjust = 1, vjust = 2.5, size = 3.5) +
  theme_few(base_size=14) +
  theme(legend.position="none")
graph2ppt(file=file.path(path, "1_4_abess model on relaxed LASSO preselected variables.pptx"), width=5, height=4)



# 1.5 L0glm fit (same solution as abess) ####

library(L0glm)
preselected_vars = relaxedlassoset # using relaxed LASSO preselection of variables
# preselected_vars = 1:ncol(X_train2) # to use all variables
X_train_subs2 = X_train2[,preselected_vars] # only main effects
X_test_subs2 = X_test2[,preselected_vars] # only main effects
# also tried adding all 1st order interactions, but fit worse
# X_train_subs2 = model.matrix(~0+.^2, # all main effects+1st order interaction effects
#                                             data=data.frame(X_train[,preselected_vars])) 
# X_test_subs2 = model.matrix(~0+.^2, # all main effects+1st order interaction effects
#                                             data=data.frame(X_test[,preselected_vars])) 
dim(as.matrix(X_train_subs2))
colnames(X_train_subs2)
L0glmfit = L0glm.fit(X=Matrix(cbind(1,X_train_subs2), sparse=T), 
                     y=logY_train,
                     no.pen=c(1, which(penalty.factor[preselected_vars]==0)+1), # do not penalize intercept & PCA factors & ROR1 SNPs
                     # lower = c(rep(-Inf,num_pcs+1), rep(0,ncol(X_train_subs2)-num_pcs)), # for nonnegativity constraints on coefs (would be best to use strain with lowest lag time as reference though & code all SNPs relative to that)
                     # upper = c(rep(Inf,num_pcs+1), rep(Inf,ncol(X_train_subs2)-num_pcs)),
                     # start = relaxedlassocoefs[c(1,preselected_vars+1)], # using SCAD coefficient estimates as starting estimates
                     # warm.start = 0.9,
                     lambda=preset.lambda("aic", 
                                          p=p, 
                                          n=n_train),
                     solver="eigen", # eigen or glmnet or osqp all work OK after variable preselection, without preselection use glmnet
                     method=7L)
# L0glmcoefs = rep(0, ncol(X)+1) 
# L0glmcoefs[c(1,preselected_vars+1)] = as.vector(L0glmfit$coefficients)
L0glmcoefs = as.vector(L0glmfit$coefficients)
names(L0glmcoefs) = c("Intercept", colnames(X_train_subs2))
L0glmset = preselected_vars[which(L0glmcoefs[-1]!=0)]
sum(L0glmcoefs!=0) # 10
L0glmcoefs[L0glmcoefs!=0]
#    Intercept                 PC1                 PC2  chromosome4_589643 
# 0.6550205477        0.0003655571        0.0012095673        0.4074851997 
# chromosome12_366378 chromosome12_367241 chromosome12_367293 chromosome12_367545 
#        0.2900666797       -0.1604696426       -0.2551567153        0.2380030727 
# chromosome16_500828 chromosome16_610743 
#        0.3857281160        0.1932357216  

# L0glmpreds_train = exp(cbind(1,X_train) %*% L0glmcoefs)-a
# L0glmpreds_test = exp(cbind(1,X_test) %*% L0glmcoefs)-a
L0glmpreds_train = exp(cbind(1,X_train_subs2) %*% L0glmcoefs)-a
L0glmpreds_test = exp(cbind(1,X_test_subs2) %*% L0glmcoefs)-a

# predictive performance on test set
predperf_test = data.frame(t(postResample(pred = log(Y_test+a), 
                                          obs = log(L0glmpreds_test+a))))
predperf_test
#      RMSE  Rsquared       MAE 
# 0.5638767 0.4204861 0.4421836

# predictive performance on training set
predperf_train = data.frame(t(postResample(pred = log(Y_train+a), 
                                           obs = log(L0glmpreds_train+a))))
predperf_train
#      RMSE  Rsquared       MAE 
# 0.4722008 0.4620488 0.3624877

model="L0glm fit on relaxed LASSO preselected variables"
test_R2_rounded <- round(predperf_test$Rsquared * 100, 0)
test_RMSE_rounded <- round(predperf_test$RMSE, 2)
train_R2_rounded <- round(predperf_train$Rsquared * 100, 0)
train_RMSE_rounded <- round(predperf_train$RMSE, 2)
subtitle_test <- bquote("test " ~ R^2 ~ " = " ~ .(test_R2_rounded) ~ "%, test RMSE = " ~ .(test_RMSE_rounded) ~ "  ")
subtitle_train <- bquote("training " ~ R^2 ~ " = " ~ .(train_R2_rounded) ~ "%, training RMSE = " ~ .(train_RMSE_rounded) ~ "  ")
data_plot <- data.frame(
  Predicted = log(c(L0glmpreds_train, L0glmpreds_test)  + a),
  Observed = log(c(Y_train, Y_test) + a),
  Group = relevel(factor(c(rep("training set", length(Y_train)),
                           rep("test set", length(Y_test)))), ref="training set")
)
ggplot(data_plot, aes(x = Predicted, y = Observed, color = Group)) +
  geom_point(shape = 16, size=2) +
  scale_x_continuous("Predicted lag time",
                     breaks = log(c(1:10, 12, 15, 20) + a),
                     labels = c(1:10, 12, 15, 20),
                     limits = c(log(0.5 + a), log(25 + a))) +
  scale_y_continuous("Observed lag time",
                     breaks = log(c(1:10, 12, 15, 20) + a),
                     labels = c(1:10, 12, 15, 20),
                     limits = c(log(0.5 + a), log(25 + a))) +
  scale_color_manual(values=c("blue", "red")) +
  labs(title=model) +
  annotate("text", x = Inf, y = Inf, label = subtitle_test, hjust = 1, vjust = 1.5, size = 3.5) +
  annotate("text", x = Inf, y = Inf, label = subtitle_train, hjust = 1, vjust = 2.5, size = 3.5) +
  theme_few(base_size=14) +
  theme(legend.position="none")
graph2ppt(file=file.path(path, "1_5_L0glm fit on relaxed LASSO preselected variables.pptx"), width=5, height=4)

