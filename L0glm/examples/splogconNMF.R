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
# Broad profiles with nearly identical shapes can share complementary ion sets.
# A conservative fusion pass runs on all original channels after preclustering.
fit$profile_fusion$merges
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

\dontrun{
# Six co-eluting compounds, retaining all 671 mass channels.
data(GCMS_data1)
Y <- GCMS_data1
X <- convmatrix(convkernel("gaussian", 65L, list(sigma = 8)), nrow(Y))
f1 <- splogconNMF(X, Y, background = matrix(1, nrow(Y), 1L),
  family = quasipoisson("identity"), precluster = 16L,
  correlation = "none", ncores = 8L)
c(compounds = f1$rank, GIC = GIC(f1), seconds = f1$seconds[["total"]])
matplot(as.numeric(rownames(Y)), as.matrix(f1$H), type = "l", lty = 1,
  xlab = "Retention time (min)", ylab = "Unit-area elution profile")

# Beer volatiles with a broad acid peak and a low-df spline background.
data(GCMS_data2)
Y <- GCMS_data2
time <- as.numeric(rownames(Y))
X <- convmatrix(convkernel("gaussian", 65L, list(sigma = 6)), nrow(Y))
B <- splines::bs(time, df = 4L, intercept = TRUE)
f2 <- splogconNMF(X, Y, background = B,
  family = quasipoisson("identity"), precluster = 32L,
  correlation = "none", ncores = 8L)
c(compounds = f2$rank, GIC = GIC(f2), seconds = f2$seconds[["total"]])
plotmat(fitted(f2), xlab = "Retention time (min)", ylab = "Mass (m/z)",
  main = "Fitted GC/MS chromatogram")
# GCMSdeco() adds CDF input, MSP/Excel export, and traces/heatmaps/profile plots.
}
