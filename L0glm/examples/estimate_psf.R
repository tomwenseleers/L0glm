set.seed(12)
n <- 300L
time <- seq_len(n) / 60
kernel <- dnorm(-15:15, sd = 2.5)
kernel <- kernel / sum(kernel)
A <- matrix(0, n, 3L)
A[c(40, 120, 210, 270), ] <- 30000
mu <- as.matrix(convmatrix(kernel, n) %*% A) + 10
Y <- matrix(rpois(length(mu), mu), n, 3L)
rownames(Y) <- time

# Isolated standard peaks provide a PSF on the sample's retention-time grid.
calibration <- estimate_psf(Y, standard = list(Y = Y, time = time),
  method = "standard", nwindows = 3L, ncores = 2L)
dim(calibration$X)
Matrix::colSums(calibration$X)[1:5]
matplot(calibration$psf_offsets, t(calibration$psf_kernels),
  type = "l", lty = 1, xlab = "Relative retention time (min)", ylab = "PSF")

# Blind refinement can use the original channels or their grouped sums.
blind <- estimate_psf(Y, standard = NULL, method = "blind", precluster = 2L,
  ncores = 2L, control = list(kernel_est = kernel, maxiter = 2L))
blind$width_scans

\dontrun{
data(GCMS_data2)
Y <- GCMS_data2
# The alkane TIC and its time axis are taken from attr(Y, "alkanes_TIC").
psf <- estimate_psf(Y, nwindows = 3L, precluster = 32L)
B <- splines::bs(as.numeric(rownames(Y)), df = 4L, intercept = TRUE)
fit <- splogconNMF(psf$X, Y, background = B,
  family = quasipoisson("identity"), precluster = 32L)

# Calibration alone, or a single common PSF, are available through the same API.
calibrated <- estimate_psf(Y, method = "standard", nwindows = 3L)
constant <- estimate_psf(Y, nwindows = 1L)
}
