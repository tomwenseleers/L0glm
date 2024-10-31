#### SIMLULATE DATA FOR RANDOM INTERCEPT & RANDOM INTERCEPT+SLOPE MODEL ####

library(lme4)
library(Matrix)
library(afex)
library(ggplot2)
library(ggthemes)
library(openxlsx)


# Set random seed for reproducibility
set.seed(12345)
set_treatment_contrasts()

# Parameters
n_clones = 1000
n_timepoints = 5
n_regimes = 3

# Create factors for clones & regimes & covariate for times
clones = factor(1:n_clones)
times = seq(0, 4, length.out = n_timepoints)
regimes = factor(1:n_regimes, labels = c("control", "heavy_metals", "elevated_temp"))

# Expand data frame
data = expand.grid(CLONE = clones, TIME = times, REGIME = regimes)
dim(data) # 15 000 x 3

# Fixed effects design matrix X (includes intercept) ####
X = model.matrix(~ TIME * REGIME, data)
dim(X) # 15000 x 6

# Random intercept and slope matrices
Z_intercept = model.matrix(~ 0 + CLONE, data) # used for random intercept data daphnia.xlsx
dim(Z_intercept) # 15 000 x 1000
Z_slope = model.matrix(~ 0 + CLONE:TIME, data)
dim(Z_slope) # 15 000 x 1000
# Combined random effect design matrix Z used for random slope data daphnia2.xlsx
Z = cbind(Z_intercept, Z_slope)
dim(Z)
# 15 000 x 2 000

# also try coding random intercepts & slopes using sum/effect coding
set_sum_contrasts()
Z_intercept_sumcoded = model.matrix(~ 0 + CLONE, data) # used for random intercept data daphnia.xlsx
Z_slope_sumcoded = model.matrix(~ 0 + CLONE:TIME, data)

# Simulate random effects ####
sd_intercept = 3 # random intercept SD
sd_slope = 2 # random slope SD
u_intercept = rnorm(n_clones, mean=0, sd=sd_intercept)  # Normally distributed random intercepts
u_slope = rnorm(n_clones, mean=0, sd=sd_slope)      # Normally distributed random slopes in function of time
# Combined random effects
u = c(u_intercept, u_slope)
length(u)
# 2000

# Normally distributed errors ####
sd_error = 1
e = rnorm(nrow(data), mean=0, sd=sd_error)

# Example fixed effect coefficients ####
beta = c(intercept = 0, 
          b_TIME = 1, 
          b_REGIMEheavy_metals = 0, b_REGIMEelevated_temp = 0, 
          b_TIME_REGIMEheavy_metals = 1, b_TIME_REGIME_elevated_temp = 2)  
# Intercept, slope for time, intercepts for each regime & slopes for each regime


# Simulate response for random intercept model ####
# a regular linear model had structure Y = X %*% beta + e
# a random intercept model will instead have the structure
Y = X %*% beta + Z_intercept %*% u_intercept + e
# where Z_intercept is a model matrix
# model.matrix(~ 0 + CLONE, data), i.e. a
# matrix with 10 columns for each clone & a one if a datapoint is from that particular clone
head(Z_intercept)
# and u_intercept are the random effect intercepts for each clone (assumed to be normally distributed)
data$LOG2_EXPR_HSP_VS_ACTIN = as.vector(Y) # update dataframe (outcome Y=log2 transformed expression of HSP gene relative to housekeeping actin gene)

# Simulate response for random slope model ####
# a regular linear model had structure Y = X %*% beta + e
# a random intercept+slope model will instead have the structure
Y2 = X %*% beta + Z %*% u + e
# where Z are two column bound model matrices
# Z_intercept = model.matrix(~ 0 + CLONE, data)
# and Z_slope = model.matrix(~ 0 + CLONE:TIME, data)
# i.e. Z = cbind(Z_intercept, Z_slope)
# and u are the contatenated random effect intercepts & random effect slopes
# (both assumed to be normally distributed)
data2 = data # make new data frame for this dataset
data2$LOG2_EXPR_HSP_VS_ACTIN = as.vector(Y2) # update dataframe (outcome Y=log2 transformed expression of HSP gene relative to housekeeping actin gene)

# plot & export random intercept data
qplot(data=data, x=TIME, y=LOG2_EXPR_HSP_VS_ACTIN, group=CLONE, geom="line", colour=CLONE, alpha=I(0.01)) +
  facet_wrap(~ REGIME, ncol=3) + theme_few() + theme(legend.position="none")
# write.xlsx(data, file="daphnia.xlsx")

# plot & export random slope data
qplot(data=data2, x=TIME, y=LOG2_EXPR_HSP_VS_ACTIN, group=CLONE, geom="line", colour=CLONE, alpha=I(0.1)) +
  facet_wrap(~ REGIME, ncol=3) + theme_few() + theme(legend.position="none")
# write.xlsx(data2, file="daphnia2.xlsx")



# comparison lmer fit & L0glm fit for random intercept model ####

system.time(fit1_lmer <- lmer(LOG2_EXPR_HSP_VS_ACTIN ~ REGIME*TIME + (1|CLONE), data=data)) # 0.22s
fit1_lmer_fixed = fixef(fit1_lmer)
fit1_lmer_fixed
# (Intercept)         TIME               REGIMEheavy_metals      REGIMEelevated_temp          REGIMEheavy_metals:TIME 
# 0.12784781    1.00070889                       0.01832759               0.03206825                       0.99492947 
# REGIMEelevated_temp:TIME 
# 1.98938563 
# true values: 0 1 0 0 1 2
fit1_lmer_u = ranef(fit1_lmer)$CLONE[,1]
mean(fit1_lmer_u) # -4.921331e-12, i.e. ca. 0
sd(fit1_lmer_u) # 2.995728
cor(fit1_lmer_u, u_intercept) # 0.9962977
hist(fit1_lmer_u)
hist(u_intercept)


# fit using L0glm
library(L0glm)
system.time(fit1_L0glm <- L0glm.fit(X=Matrix(cbind(X, Z_intercept), sparse=T),  # fastest if coded as sparse
                                    y=data$LOG2_EXPR_HSP_VS_ACTIN, no.pen=1, lambda="aic")) # 0.14s
fit1_L0glm_fixed = head(coef(fit1_L0glm),-length(levels(data$CLONE)))
fit1_L0glm_u = tail(coef(fit1_L0glm),length(levels(data$CLONE)))
mean(fit1_L0glm_u) # 3.938169, should be 0
# center u at zero and add shift to fixed effect intercept
fit1_L0glm_fixed[1] = fit1_L0glm_fixed[1]+mean(fit1_L0glm_u)
fit1_L0glm_fixed 
# 0.1446464 0.9951094 0.0000000 0.0000000 1.0010387 2.0000750 
# closely match true values: 0 1 0 0 1 2
fit1_L0glm_u = fit1_L0glm_u-mean(fit1_L0glm_u)
mean(fit1_L0glm_u) # 1.132214e-16
sd(fit1_L0glm_u) # 3.015034
cor(fit1_L0glm_u, u_intercept) # 0.9954688
hist(fit1_L0glm_u) # artificial 2nd peak at -mean(fit1_L0glm_u1) (coefficients that were originally shrunk to zero) with lambda="gic", but not with lambda="aic"
# TODO this artificial 2nd peak could be resolved if random effect coefficients would be zero centered 
# internally in each IRLS iteration, TODO allow extra argument Z to pass list of random effect design matrices
hist(u_intercept)



# comparison lmer fit & L0glm fit for random slope model ####

system.time(fit2_lmer <- lmer(LOG2_EXPR_HSP_VS_ACTIN ~ REGIME*TIME + (1+TIME|CLONE), data=data2)) # 0.44s
summary(fit2_lmer)
# Random effects:
# Groups   Name        Variance Std.Dev. Corr
# CLONE    (Intercept) 8.9154   2.9859       
# TIME        4.0820   2.0204   0.04
# Residual             0.9904   0.9952       
# Number of obs: 15000, groups:  CLONE, 1000
# 
# Fixed effects:
#                           Estimate Std. Error        df t value Pr(>|t|)    
# (Intercept)              1.278e-01  9.752e-02 1.088e+03   1.311    0.190    
# REGIMEheavy_metals       1.833e-02  3.447e-02 1.300e+04   0.532    0.595    
# REGIMEelevated_temp      3.207e-02  3.447e-02 1.300e+04   0.930    0.352    
# TIME                     9.399e-01  6.466e-02 1.031e+03  14.536   <2e-16 ***
# REGIMEheavy_metals:TIME  9.949e-01  1.407e-02 1.300e+04  70.691   <2e-16 ***
# REGIMEelevated_temp:TIME 1.989e+00  1.407e-02 1.300e+04 141.349   <2e-16 ***

fit2_lmer_fixed = fixef(fit2_lmer)
fit2_lmer_fixed
# (Intercept)         TIME               REGIMEheavy_metals      REGIMEelevated_temp          REGIMEheavy_metals:TIME 
# 0.12784781    0.93992393                       0.01832759               0.03206825                       0.99492947 
# REGIMEelevated_temp:TIME 
# 1.98938563 
# true values: 0 1 0 0 1 2
fit2_lmer_u1 = ranef(fit2_lmer)$CLONE[,1]
fit2_lmer_u2 = ranef(fit2_lmer)$CLONE[,2]
mean(fit2_lmer_u1) # 4.908764e-12, i.e. ca. 0
sd(fit2_lmer_u1) # 2.953434
cor(fit2_lmer_u1, u_intercept) # 0.989619
hist(fit2_lmer_u1)
hist(u_intercept)
mean(fit2_lmer_u2) # 3.544503e-12, i.e. ca. 0
sd(fit2_lmer_u2) # 2.0124
cor(fit2_lmer_u2, u_slope) # 0.9961354
hist(fit2_lmer_u2)
hist(u_slope)



# fit using L0glm
library(L0glm)
system.time(fit2_L0glm <- L0glm.fit(X=Matrix(cbind(X, Z_intercept, Z_slope), sparse=T),  # fastest if coded as sparse
                                    y=data2$LOG2_EXPR_HSP_VS_ACTIN, no.pen=c(1,2), lambda="aic")) # 0.20s
# NOTE/TODO solution with "aic" best - with "gic" some artefacts in random effect coefficients (some shrunk to zero)
# maybe prior centers should be updated after each IRLS iteration for random effect coefficients to
# make sure random effect coefficients are non-sparse
# Maybe use lambda="aic" for random effect coefficients only?
# After each IRLS iteration random effect coefficients should also be centered
fit2_L0glm_fixed = head(coef(fit2_L0glm),-2*length(levels(data$CLONE)))
fit2_L0glm_u = tail(coef(fit2_L0glm),2*length(levels(data$CLONE)))
fit2_L0glm_u1 = head(fit2_L0glm_u,length(levels(data$CLONE)))
fit2_L0glm_u2 = tail(fit2_L0glm_u,length(levels(data$CLONE)))
mean(fit2_L0glm_u1) # 1.313248, should be 0
# center u1 at zero and add shift to fixed effect intercept
fit2_L0glm_fixed[1] = fit2_L0glm_fixed[1]+mean(fit2_L0glm_u1)
mean(fit2_L0glm_u2) # 1.665777, should be 0
# center u2 at zero and add shift to fixed effect of TIME
fit2_L0glm_fixed[2] = fit2_L0glm_fixed[2]+mean(fit2_L0glm_u2)
fit2_L0glm_fixed 
# 0.1446464 0.9343244 0.0000000 0.0000000 1.0010387 2.0000750
# closely match true values: 0 1 0 0 1 2
# center u1 & u2
fit2_L0glm_u1 = fit2_L0glm_u1-mean(fit2_L0glm_u1)
fit2_L0glm_u2 = fit2_L0glm_u2-mean(fit2_L0glm_u2)
mean(fit2_L0glm_u1) # -1.198724e-16
sd(fit2_L0glm_u1) # 2.970462
cor(fit2_L0glm_u1, u_intercept) # 0.976875
mean(fit2_L0glm_u2) # -6.654746e-17
sd(fit2_L0glm_u2) # 2.018681
cor(fit2_L0glm_u2, u_slope) # 0.9918059


hist(fit2_L0glm_u1) # artificial 2nd peak at -mean(fit2_L0glm_u1) (coefficients that were originally shrunk to zero) with lambda="gic", but not with lambda="aic"
# this artificial 2nd peak could be resolved if random effect coefficients were zero centered 
# internally in each IRLS iteration, TODO allow extra argument Z to pass list of random effect design matrices
# & add this zero centering in each IRLS iteration of random effect coefficients
hist(u_intercept)
hist(fit2_L0glm_u2) # artificial 2nd peak at -mean(fit2_L0glm_u2) (coefficients that were originally shrunk to zero) with lambda="gic", but not with lambda="aic"
# this artificial 2nd peak could be resolved if random effect coefficients were zero centered 
# internally in each IRLS iteration, TODO allow extra argument Z to pass list of random effect design matrices
# & add this zero centering in each IRLS iteration of random effect coefficients
hist(u_slope)
