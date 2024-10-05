# TEST USE OF SPARSE MODEL MATRICES VIA MatrixModels & model.Matrix(*, sparse=TRUE) ####

library(MatrixModels)
library(microbenchmark)
library(emmeans)
library(effects)

# see ?model.Matrix & ?

data(CO2, package="datasets")

## dense linear model
fit_dense = glm4(uptake ~ 0 + Type*Treatment, data=CO2, doFit = TRUE)
# class(fit_dense) = "glm"
emmeans(fit_dense, ~Treatment) # not supported
effects(fit_dense) # not supported
microbenchmark(glm4(uptake ~ 0 + Type*Treatment, data=CO2, doFit = TRUE)) # 2.3ms

str(glm4(uptake ~ 0 + Type*Treatment, data=CO2, doFit = TRUE), 4)

## sparse linear model
fit_sparse = glm4(uptake ~ 0 + Type*Treatment, data=CO2, doFit = TRUE,
     sparse = TRUE)
emmeans(fit_sparse, ~Treatment) # not supported
effects(fit_sparse) # not supported
microbenchmark(glm4(uptake ~ 0 + Type*Treatment, data=CO2, doFit = TRUE,
                    sparse = TRUE)) # 4.9ms
microbenchmark(MatrixModels:::lm.fit.sparse(x=model.Matrix(uptake ~ 0 + Type*Treatment,
                                              data=CO2, sparse=TRUE),
                                            y=CO2$uptake, method="cholesky")) # 3.6ms
str(glm4(uptake ~ 0 + Type*Treatment, data=CO2, doFit = TRUE,
         sparse = TRUE), 4)
