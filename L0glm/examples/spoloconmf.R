set.seed(12)
n <- 300L
psf <- convkernel("gaussian", 21L, list(sigma = 2))
X <- convmatrix(psf, n)
H <- cbind(X[, c(45, 95, 160, 235)],
           pmax(0, 1 - abs(seq_len(n) - 185) / 22))
H <- H %*% Matrix::Diagonal(x = 1 / Matrix::colSums(H))
S <- matrix(rexp(5 * 8, rate = 1 / 3000), 5, 8)
B <- splines::bs(seq_len(n), df = 4, intercept = TRUE)
mu <- as.matrix(H %*% S) + B %*% matrix(20, 4, 8)
Y <- matrix(rpois(length(mu), mu), n, 8)
fit <- spoloconmf(X, Y, background = B, ncores = 2L)
fit
matplot(fit$H, type = "l", lty = 1, xlab = "Scan", ylab = "Elution profile")
