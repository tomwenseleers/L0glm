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
fit <- splogconNMF(X, Y, background = B, family = poisson("identity"),
  precluster = 4L, correlation = "none", ncores = 2L)
fit
c(compounds = fit$rank, GIC = GIC(fit), seconds = fit$seconds[["total"]])
matplot(fit$H, type = "l", lty = 1, xlab = "Scan", ylab = "Elution profile")

# For overdispersed intensities, estimate one common noise scale from temporal
# differences. Selection and fusion use four groups; the result has all channels.
Yq <- 5 * matrix(rpois(length(mu), mu / 5), n, 8)
fitq <- splogconNMF(X, Yq, background = B, family = quasipoisson("identity"),
  precluster = 4L, correlation = "none", ncores = 2L)
GIC(fitq)

# Other identity-link choices include gaussian(), mgcv::nb(), Gamma(),
# nb1(dispersion=5, link="identity"), and
# Gamma_scale(scale=5, link="identity", censor=0.5).
