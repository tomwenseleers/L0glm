# Sparse Gaussian regression.
set.seed(1)
X <- matrix(rnorm(360),60,6)
truth <- c(2,-1,rep(0,4))
y <- as.numeric(X %*% truth + rnorm(60))
comparison <- sparseGLM(X,y,method=c("lasso","relaxed_lasso"),
  tune="gic",nlambda=5,beta_true=truth)
l0glmMetrics(comparison)[c("method","seconds","gic","adj_R2_weighted","F1_score","MCC")]
plot(comparison,which=c(1,2))

# The same SIS/ISIS controls apply to every supported method. Each ISIS round
# uses that method's own selected model to screen for additional predictors.
set.seed(2)
wide <- matrix(rnorm(100*2000),100,2000)
truth_wide <- c(3,-2,1.5,rep(0,1997))
y_wide <- as.numeric(wide%*%truth_wide)+rnorm(100)
screened <- sparseGLM(wide,y_wide,method=c("L0glm","lasso","relaxed_lasso"),
  screen="ISIS",nsis=10,nlambda=5,beta_true=truth_wide,
  lambda=list(L0glm=presetLambda("gic",ncol(wide),nrow(wide))),
  options=list(L0glm=list(relaxed=TRUE)))
l0glmMetrics(screened)[c("method","screen","seconds","gic","TP","FP","F1_score","MCC")]
# With tune="cv", all screening and adaptive pilots use training rows only.
screened_cv <- sparseGLM(wide,y_wide,method="lasso",screen="SIS",nsis=20,
  tune="cv",nfolds=3,nlambda=5,beta_true=truth_wide)
l0glmMetrics(screened_cv)[c("method","screen","seconds","TP","FP","F1_score","MCC")]

# Blurred Poisson spike trains: smaller versions of the single/multichannel
# deconvolution benchmarks. MCP and SCAD are omitted because they were slow
# on the large deconvolution examples. All times below are measured locally.
comparison_methods <- c("L0glm","abess","L0Learn","lasso","relaxed_lasso",
                        "adaptive_lasso","adaptive_relaxed_lasso")
fields <- c("method","seconds","gic","deviance","adj_R2_weighted","F1_score","MCC")

# Single channel: 50 spikes among 500 candidate locations.
single <- simul(method=2,n=500,p=500,s=.1,seed=2)
lambda_gic <- log(ncol(single$X))*log(log(nrow(single$X)))
comparison_single <- sparseGLM(single,family=poisson(identity),
  method=comparison_methods,tune="gic",nlambda=10,
  nonnegative=TRUE,intercept=FALSE,approximation="wls",
  support.size=1:100,lambda=list(L0glm=lambda_gic),
  options=list(L0glm=list(relaxed=FALSE),abess=list(warm.start=TRUE)))
l0glmMetrics(comparison_single)[fields]
plot(comparison_single,which=c(1,2,4,5))

# Multichannel: use 100 time points and six abundant mass channels for speed.
multi <- simul(method=4,n=100,p=100,s=.1,seed=2)
keep <- head(order(Matrix::colSums(multi$y),decreasing=TRUE),6)
for (name in c("y","y_true","beta_true")) multi[[name]] <- multi[[name]][,keep,drop=FALSE]
lambda_gic <- log(ncol(multi$X))*log(log(nrow(multi$X)))
comparison_multi <- sparseGLM(multi,family=poisson(identity),
  method=comparison_methods,tune="gic",nlambda=10,
  nonnegative=TRUE,intercept=FALSE,approximation="wls_grouped",
  support.size=1:20,lambda=list(L0glm=lambda_gic),
  options=list(L0glm=list(relaxed=TRUE),abess=list(warm.start=TRUE)))
# ABESS and glmnet select rows jointly using Gaussian shared-row-weight WLS;
# L0glm fits the exact grouped Poisson model. L0Learn fits channels independently.
# GIC, deviance and prediction R-squared always use the original Poisson model.
l0glmMetrics(comparison_multi)[c(fields,"within_row_F1_score","within_row_MCC")]

plot(comparison_multi,which=c(1,4),type="coefs",channels="all")
plot(comparison_multi,which=c(1,2,4,5),channels="mean")
# One ground-truth/estimate pair; heatmaps default to which=1.
plot(comparison_multi,type="heatmap")
# Two pages, each with its own ground-truth/estimate pair.
plot(comparison_multi,type="heatmap",which=c(1,2))
