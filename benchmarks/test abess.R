library(abess)

# MULTICHANNEL EXAMPLE ####

n = 500 # obs
p = 1000 # features
y.dim = 10 # nr of tasks/channels
k = 100 # features with nonzero coefs

coef_real = matrix(0, nrow=p, ncol=y.dim)
nonzero_features = sample(1:p, k, replace = FALSE)
coef_real[nonzero_features,] = pmax(rnorm(k*y.dim, mean=1, sd=1), 0)

data = generate.data(n, p, support.size = k, 
                      family = "mgaussian", # here simulating multitask case with Gaussian noise, but also interested in multitask case with Poisson noise (& identity link)
                      beta = coef_real,
                      y.dim = y.dim)

abess_fit = abess(
  x=data$x, 
  y=data$y,
  # weight=1/(data$y+0.1), # I'd like to use this to approximate identity link Poisson case - now gives error "Error in weight.glm(para) : is.vector(para$weight) is not TRUE"
  family="mgaussian",
  beta.low=0, # nonnegativity constraints
  fit.intercept=FALSE,
  tune.path="sequence",
  tune.type="aic"
)

coef_abess = coef(abess_fit, 
                  support.size = abess_fit[["best.size"]], 
                  sparse = FALSE)[[1]][-1,]
# NOTE: extract(abess_fit)
# gives error with fit.intercept=FALSE
# Error in h(simpleError(msg, call)) : 
# error in evaluating the argument 'x' in selecting a method for function 'as.matrix': incorrect number of dimensions
# this looks like a bug

plot(x=rowMeans(coef_abess), y=rowMeans(coef_real), pch=16, col='blue')

image(x=1:p, y=1:y.dim, z=coef_real, col = topo.colors(255), useRaster=TRUE) 
abline(v=(1:p)[as.vector(rowMax(coef_abess)!=0)], col="red")

table(rowMeans(coef_abess)!=0, 
      rowMeans(coef_real)!=0,
      dnn=c("estimated coefs nonzero", "true coefs nonzero"))
#                         true coefs nonzero
# estimated coefs nonzero FALSE TRUE
#                   FALSE   900   60
#                    TRUE     0   40


# SINGLE CHANNEL EXAMPLE ####

n = 500 # obs
p = 1000 # features
k = 100 # features with nonzero coefs

coef_real = rep(0, p)
nonzero_features = sample(1:p, k, replace = FALSE)
coef_real[nonzero_features] = pmax(rnorm(k, mean=1, sd=1), 0)

data = generate.data(n, p, support.size = k, 
                     family = "gaussian", # here simulating with Gaussian noise, but also interested in identity link Poisson case
                     beta = coef_real)

abess_fit = abess(
  x=data$x, 
  y=data$y,
  weight=as.vector(1/(pmax(data$y,0)+0.1)), 
  family="gaussian",
  beta.low=0, # nonnegativity constraints
  fit.intercept=FALSE,
  tune.path="sequence",
  tune.type="aic"
)

coef_abess = coef(abess_fit, 
                  support.size = abess_fit[["best.size"]], 
                  sparse = FALSE)[-1,]
# NOTE: extract(abess_fit)$beta also works

plot(x=coef_abess, y=coef_real, pch=16, col='blue')

table(coef_abess!=0, 
      coef_real!=0,
      dnn=c("estimated coefs nonzero", "true coefs nonzero"))
#                         true coefs nonzero
# estimated coefs nonzero FALSE TRUE
#                   FALSE   912   49
#                    TRUE     5   34


