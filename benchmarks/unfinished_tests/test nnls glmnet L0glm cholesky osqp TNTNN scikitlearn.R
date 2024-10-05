# BENCHMARK nnls, glmnet & L0glm to solve nonnegative least square problem ####

# see pytntnn tntnn_x_sklearn_simple_comparison.py
# for comparison with sci-kit learn & TNT-NN for this example
# also see Chou et al. 2022
# Non-negative Least Squares via Overparametrization

options(repos = c(CRAN="https://cloud.r-project.org/"))
library(ggplot2)
library(dplyr)
library(Matrix)
# library(RcppML)
devtools::load_all(".")
# library(L0glm)

library(nnls)
library(glmnet)
library(keras)

set.seed(1)

# Generate some random data
n_samples <- 20000
n_features <- 10000
X <- matrix(rnorm(n_samples * n_features), n_samples, n_features)
true_coef <- 3 * rnorm(n_features)

# Threshold coefficients to render them non-negative
true_coef[true_coef < 0] <- 0
true_nonzero = sum(true_coef!=0)
print(paste0("True nonzero coefs: ", true_nonzero))
# "True nonzero coefs: 5004"

y <- X %*% true_coef

# Add some noise
y <- y + 5 * rnorm(n_samples)

# Split the data in train set and test set
idx_train <- sample(1:n_samples, n_samples / 2)
X_train <- X[idx_train, ]
y_train <- y[idx_train]
X_train_sp <- Matrix(X_train, sparse=T)
X_test <- X[-idx_train, ]
y_test <- y[-idx_train]


# Solve using nnls from the nnls package (classic Fortran Lawson & Hanson approach) ####
t <- system.time( fit_nnls <- nnls::nnls(A = X_train, b = y_train)$x )

print(paste0("NNLS execution time: ", t[[3]]))
predicted_nnls <- X_test %*% fit_nnls
r2_score_nnls <- 1 - sum((y_test - predicted_nnls)^2) / sum((y_test - mean(y_test))^2)
print(paste0("R2 Score NNLS: ", r2_score_nnls))
nonzero_nnls <- sum(fit_nnls!=0)
print(paste0("Nonzero coefs NNLS: ", nonzero_nnls))

# "NNLS execution time:  902.740000000005"
# "R2 Score NNLS: 0.997749375716445"
# "Nonzero coefs NNLS: 7471"


# Solve using glmnet (analogous to scikit-learn's coordinate descent) ####
t <- system.time(fit_cd <- glmnet(X_train, y_train, 
                                  alpha = 0, 
                                  lambda = 0,
                                  lower.limits = 0))

print(paste0("glmnet Coordinate Descent execution time: ", t[[3]]))
coefs_cd <- coef(fit_cd, s=0)
predicted_cd <- predict(fit_cd, s = 0, newx = X_test)
r2_score_cd <- 1 - sum((y_test - predicted_cd)^2) / sum((y_test - mean(y_test))^2)
print(paste0("R2 Score glmnet Coordinate Descent: ", r2_score_cd))
nonzero_cd <- sum(coefs_cd!=0)
print(paste0("Nonzero coefs glmnet Coordinate Descent: ", nonzero_cd))

# "glmnet Coordinate Descent execution time: 2.71999999998661"
# "R2 Score glmnet Coordinate Descent: 0.995265957065691"
# "Nonzero coefs glmnet Coordinate Descent: 7797"
# TODO: THIS ONE COULD BE A GOOD INITIALISATION FOR THE DENSE CASE WHERE THERE 
# ARE NONNEGATIVITY CONSTRAINTS


# Solve using glmnet::bigGlm ####
t <- system.time(fit_cd <- glmnet::bigGlm(X_train, y_train, 
                                  lower.limits = 0))

print(paste0("glmnet::bigGlm execution time: ", t[[3]]))
coefs_cd <- coef(fit_cd, s=0)
predicted_cd <- predict(fit_cd, s = 0, newx = X_test)
r2_score_cd <- 1 - sum((y_test - predicted_cd)^2) / sum((y_test - mean(y_test))^2)
print(paste0("R2 Score glmnet::bigGlm: ", r2_score_cd))
nonzero_cd <- sum(coefs_cd!=0)
print(paste0("Nonzero coefs glmnet::bigGlm: ", nonzero_cd))

# "glmnet Coordinate Descent execution time: 2.66999999999825"
# "R2 Score glmnet Coordinate Descent: 0.995265957065691"
# "Nonzero coefs glmnet Coordinate Descent: 7797"


# # Solve using FAST-nnls (Fast Approximate Solution Trajectory-NNLS) ####
# see https://stackoverflow.com/questions/58006606/rcpp-implementation-of-fast-nonnegative-least-squares#:~:text=4.%20EDIT%20Jan%202022:%20Use%20the%20nnls%20function
t <- system.time(fit_fastnnls <- RcppML::nnls(a = crossprod(X_train),
                                              b = crossprod(X_train, y_train),
                                              cd_maxit = 0L, # set to 0L to use FAST-NNLS without coordinate descent step
                                              fast_nnls=TRUE
                                              )) 
# # t <- system.time( fit_fastnnls <- fastnnls_dense(A = X_train, b = y_train))
# # t <- system.time( fit_fastnnls <- fastnnls_sparse(A = X_train_sp, b = y_train))

print(paste0("FAST-NNLS execution time: ", t[[3]]))
predicted_fastnnls <- X_test %*% fit_fastnnls
r2_score_fastnnls <- 1 - sum((y_test - predicted_fastnnls)^2) / sum((y_test - mean(y_test))^2)
print(paste0("R2 Score FAST-NNLS: ", r2_score_fastnnls))
nonzero_fastnnls <- sum(fit_fastnnls!=0)
print(paste0("Nonzero coefs FAST-NNLS: ", nonzero_fastnnls))

"FAST-NNLS execution time:  486.139999999999"
"R2 Score FAST-NNLS: 0.391182631236245" # realy poor
"Nonzero coefs FAST-NNLS: 4244"



# Solve using osqp ####
t <- system.time( fit_osqp <- L0glm.fit(
  X = X_train_sp, y = y_train,
  start = rep(1, ncol(X_train)),
  warm.start = 1,
  family = gaussian(link="identity"), 
  lambda = 0, 
  nonnegative = TRUE,
  solver = "osqp")$coefficients )

print(paste0("osqp execution time: ", t[[3]]))
predicted_osqp <- as.vector(X_test %*% fit_osqp)
r2_score_osqp <- 1 - sum((y_test - predicted_osqp)^2) / sum((y_test - mean(y_test))^2)
print(paste0("R2 Score osqp: ", r2_score_osqp))
nonzero_osqp <- sum(fit_osqp!=0)
print(paste0("Nonzero coefs osqp: ", nonzero_osqp))





# Solve using L0glm via LSCG ####
t <- system.time( fit_L0glm_LSCG <- L0glm.fit(
  X = Matrix(X_train, sparse=T), # note: due to recent change in Matrix package error when using dense X Error in as(Matrix(rbind(Xw, Diagonal(x = sqrt(lambdas))), sparse = T),  : no method or default for coercing “dgRMatrix” to “dgCMatrix” 
  y = y_train,
  start = rep(1, ncol(X_train)),
  warm.start = 1,
  family = gaussian(link="identity"), 
  lambda = "gic", 
  nonnegative = TRUE,
  solver = "eigen", method = 3L)$coefficients )

print(paste0("L0glm LSCG execution time: ", t[[3]]))
predicted_L0glm_LSCG <- as.vector(X_test %*% fit_L0glm_LSCG)
r2_score_L0glm_LSCG <- 1 - sum((y_test - predicted_L0glm_LSCG)^2) / sum((y_test - mean(y_test))^2)
print(paste0("R2 Score L0glm LSCG: ", r2_score_L0glm_LSCG))
nonzero_L0glm_LSCG <- sum(fit_L0glm_LSCG!=0)
print(paste0("Nonzero coefs L0glm LSCG: ", nonzero_L0glm_LSCG))

# "L0glm LSCG execution time: 107.760000000009"
# "R2 Score L0glm LSCG: 0.995738817265223"
# "Nonzero coefs L0glm LSCG: 4264"




# Solve using L0glm via Cholesky LLT ####
# t <- system.time( fit_L0glm_chol <- L0glm.fit(
#   X = Matrix(X_train, sparse=T), y = y_train,
#   start = rep(1, ncol(X_train)),
#   warm.start = 1,
#   family = gaussian(link="identity"),
#   lambda = "bic",
#   nonnegative = TRUE,
#   solver = "eigen", method = 7L)$coefficients )
# 
# print(paste0("L0glm chol execution time: ", t[[3]]))
# predicted_L0glm_chol <- as.vector(X_test %*% fit_L0glm_chol)
# r2_score_L0glm_chol <- 1 - sum((y_test - predicted_L0glm_chol)^2) / sum((y_test - mean(y_test))^2)
# print(paste0("R2 Score L0glm chol: ", r2_score_L0glm_chol))
# nonzero_L0glm_chol <- sum(fit_L0glm_chol!=0)
# print(paste0("Nonzero coefs L0glm chol: ", nonzero_L0glm_chol))
# 
# # "L0glm chol execution time: XXX"
# # "R2 Score L0glm chol: XXX"
# # "Nonzero coefs L0glm chol: XXX"
# THIS ONE TAKES FOREVER


# Solve using L0glm via glmnet ####
t <- system.time( fit_L0glm_glmnet <- L0glm.fit(
  X = X_train, y = y_train,
  start = rep(1, ncol(X_train)), # slightly faster with coefs_cd[-1]
  warm.start = 1,
  family = gaussian(link="identity"), 
  lambda = "bic", 
  nonnegative = TRUE,
  solver = "glmnet")$coefficients )

print(paste0("L0glm glmnet execution time: ", t[[3]]))
predicted_L0glm_glmnet <- as.vector(X_test %*% fit_L0glm_glmnet)
r2_score_L0glm_glmnet <- 1 - sum((y_test - predicted_L0glm_glmnet)^2) / sum((y_test - mean(y_test))^2)
print(paste0("R2 Score L0glm glmnet: ", r2_score_L0glm_glmnet))
nonzero_L0glm_glmnet <- sum(fit_L0glm_glmnet!=0)
print(paste0("Nonzero coefs L0glm glmnet: ", nonzero_L0glm_glmnet))

# "L0glm glmnet execution time: 140.770000000004"
# "R2 Score L0glm glmnet: 0.99866311501799"
# "Nonzero coefs L0glm glmnet: 4690"




# Solve using L0glm via osqp ####
t <- system.time( fit_L0glm_osqp <- L0glm.fit(
  X = X_train_sp, y = y_train,
  start = rep(1, ncol(X_train)),
  warm.start = 1,
  family = gaussian(link="identity"), 
  lambda = "bic", 
  nonnegative = TRUE,
  solver = "osqp")$coefficients )

print(paste0("L0glm osqp execution time: ", t[[3]]))
predicted_L0glm_osqp <- as.vector(X_test %*% fit_L0glm_osqp)
r2_score_L0glm_osqp <- 1 - sum((y_test - predicted_L0glm_osqp)^2) / sum((y_test - mean(y_test))^2)
print(paste0("R2 Score L0glm osqp: ", r2_score_L0glm_osqp))
nonzero_L0glm_osqp <- sum(fit_L0glm_osqp!=0)
print(paste0("Nonzero coefs L0glm chol: ", nonzero_L0glm_osqp))

# "L0glm osqp execution time: XXs"
# "R2 Score L0glm osqp: XX"
# "Nonzero coefs L0glm osqp: XX"
# THIS ONE TAKES FOREVER




# TEST ON MNIST number 3 PROBLEM ####
# see Chou et al. 2022
# Non-negative Least Squares via Overparametrization

# Load necessary libraries
# install.packages("keras")
library(keras)

# Load MNIST data
mnist <- dataset_mnist()
train_images <- mnist$train$x

# Select an image of the digit "3"
digit_three_images <- train_images[mnist$train$y == 3,,]
selected_image <- digit_three_images[1,,]
image(1:28, 1:28, t(selected_image^3), 
      col=grey.colors(256), useRaster=T,
      axes=FALSE, main="MNIST Digit '3'")

# Flatten and normalize the image to have values between 0 and 1
x_plus <- as.vector(selected_image) / 255

# Set matrix A dimensions
M <- 300
N <- 784
A <- matrix(rnorm(M * N), nrow=M)
image(t(A^3), 
      col=grey.colors(256), useRaster=T,
      axes=FALSE)

# Define function to corrupt the signal based on q value
corrupt_signal <- function(x_plus, q) {
  x_minus <- rnorm(length(x_plus))
  x <- x_plus + q * x_minus
  return(list(x = x, x_minus = x_minus))
}
# define function to normalize to 2 norm
normalize_to_2norm <- function(x, desired_norm) {
  current_norm <- sqrt(sum(x^2))
  return((desired_norm / current_norm) * x)
}

# Example: For q = 0.05
q <- 0.05
corrupted_data <- corrupt_signal(x_plus, q)
x <- corrupted_data$x
image(matrix(x^3, 
             ncol=sqrt(length(x)),
             nrow=sqrt(length(x)), byrow=T), 
      col=grey.colors(256), useRaster=T,
      axes=FALSE)
x_minus <- corrupted_data$x_minus
image(matrix(x_minus^3, 
             ncol=sqrt(length(x)),
             nrow=sqrt(length(x)), byrow=T), 
      col=grey.colors(256), useRaster=T,
      axes=FALSE)

# Rescale the ground truth image to have 2-norm (1 - q)
x_plus_rescaled <- normalize_to_2norm(x_plus, 1 - q)
image(matrix(x_plus_rescaled^3, 
                  ncol=sqrt(length(x_plus_rescaled)),
                  nrow=sqrt(length(x_plus_rescaled)), byrow=T), 
      col=grey.colors(256), useRaster=T,
      axes=FALSE)

# Measurements
y <- A %*% matrix(x, ncol=1)

# Use the nnls function to solve for x
system.time(reconstruction_nnls <- nnls(A, y)$x) # 0.08s
image(matrix(reconstruction_nnls^3, 
             ncol=sqrt(length(reconstruction_nnls)),
             nrow=sqrt(length(reconstruction_nnls)), byrow=T), 
      col=grey.colors(256), useRaster=T,
      axes=FALSE)

# Use L0glm to solve for x
system.time(reconstruction_L0glm_osqp <- L0glm.fit(
  X = A, y = y,
  start = rep(1, ncol(A)),
  warm.start = 1,
  gamma = 1, # 1 worked better than 2 here
  family = gaussian(link="identity"), 
  lambda = preset.lambda(IC="gic", p=ncol(A), n=nrow(A)),
  nonnegative = TRUE,
  solver = "osqp", # 0.22s with LSCG, 3.41s with osqp, 12s with glmnet
  method = 3L)$coefficients) 
image(matrix(reconstruction_L0glm_osqp^3, 
             ncol=sqrt(length(reconstruction_L0glm_osqp)),
             nrow=sqrt(length(reconstruction_L0glm_osqp)), byrow=T), 
      col=grey.colors(256), useRaster=T,
      axes=FALSE)

system.time(reconstruction_L0glm_LSCG <- L0glm.fit(
  X = A, y = y,
  start = rep(1, ncol(A)),
  warm.start = 1,
  gamma = 1, # 1 worked better than 2 here
  family = gaussian(link="identity"), 
  lambda = preset.lambda(IC="gic", p=ncol(A), n=nrow(A)),
  nonnegative = TRUE,
  solver = "eigen", # 0.22s with LSCG, 3.41s with osqp, 12s with glmnet
  method = 3L)$coefficients)  # = LSCG
image(matrix(reconstruction_L0glm_LSCG^3, 
             ncol=sqrt(length(reconstruction_L0glm_LSCG)),
             nrow=sqrt(length(reconstruction_L0glm_LSCG)), byrow=T), 
      col=grey.colors(256), useRaster=T,
      axes=FALSE)
# superior performance compared to nnls, though slower for this small example