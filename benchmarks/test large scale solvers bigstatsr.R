# test some large-scale solvers ####


# bigstatsr::big_spLinReg ####

library(bigstatsr)
n <- 10000
p <- 10000
X <- FBM(nrow = n, ncol = p, 
         init = rnorm(n*p, sd = 1),
         type="double", backingfile = "big_matrix2")
y <- rnorm(n)

system.time(model <- big_spLinReg(X, y, ind.train = rows_along(X),
                                  alphas = 1E-4, nlambda = 10,
                                  pf.X = runif(ncol(X)), ncores=1)) # 17s
plot(model)


# biglm & biganalytics ####

library(biglm)
library(biganalytics)
library(bigmemory)
X <- big.matrix(init=rnorm(n*p, sd = 1), 
                nrow=n, ncol=p,
                   backingfile = "big_matrix4.bin",
                   descriptorfile = "big_matrix4.desc",
                   dimnames = list(NULL, paste0("V", 1:p)))
head(X)
# Generate the new column vector y
y <- rnorm(n)

# Create a new big.matrix to store the combined data (original + new column)
combined_matrix <- filebacked.big.matrix(nrow = n, ncol = p + 1,
                                         backingfile = "combined_matrix4.bin",
                                         descriptorfile = "combined_matrix4.desc",
                                         dimnames = list(NULL, c(paste0("V", 1:p), "y")))
# Copy the original matrix into the new big.matrix
combined_matrix[, 1:p] <- X[, ]

# Add the new column y to the new big.matrix
combined_matrix[, p + 1] <- y

# or use read.big.matrix to read in file

# Create a formula for the regression
predictor_names <- paste0("V", 1:p)
formula_string <- paste("y ~", paste(predictor_names, collapse = " + "))
options(expressions=p+1)
options(max_ppsize=p) # --max-ppsize 50000
formula <- as.formula(formula_string)

biglm <- bigglm.big.matrix(formula=formula_string, data=combined_matrix, family = gaussian())
coef(biglm.2)

# also
bigalgebra

bigFastlm

biglasso





# speedglm ####

library(speedglm)



# torch ####
library(torch)



# gpuMatrix ####



# sgd ####



# spark MLLib using sparkR ####
# see test sparkR.R
