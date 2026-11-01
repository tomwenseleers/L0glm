\dontrun{
library(L0glm)
data(GCMS_data1)
data(GCMS_data2)
out <- file.path(tempdir(), "L0glm-GCMS")
dir.create(out, showWarnings = FALSE)
write_GCMS_cdf(GCMS_data1, file.path(out, "GCMS_data1.cdf"))
write_GCMS_cdf(GCMS_data2, file.path(out, "GCMS_data2.cdf"))

fit1 <- GCMSdeco(file.path(out, "GCMS_data1.cdf"),
  file.path(out, "GCMS_data1.msp"), nwindows = 1L, precluster = 16L,
  baseline_df = 1L,
  correlation = "none", ncores = 8L,
  psf_control = list(kernel_est = convkernel("gaussian", 65L,
    list(sigma = 8)), maxiter = 0L))
fit2 <- GCMSdeco(file.path(out, "GCMS_data2.cdf"),
  file.path(out, "GCMS_data2.msp"), nwindows = 1L, precluster = 32L,
  correlation = "none", ncores = 8L,
  psf_control = list(kernel_est = convkernel("gaussian", 65L,
    list(sigma = 6)), maxiter = 0L))
fit1$compounds
GIC(fit2)
coef(fit2, normalize = "l1")

plot(fit1, type = "lines", nrows = 1L,
  minutes_per_panel = diff(range(fit1$time)) + .01, labels = TRUE)
plot(fit1, type = "heatmap")
plot(fit1, type = "profiles", nrows = 3L)
plot(fit2, type = "lines", file = file.path(out, "GCMS_data2_traces.pdf"),
  nrows = 1L, minutes_per_panel = diff(range(fit2$time)) + .01)
plot(fit2, type = "heatmap", file = file.path(out, "GCMS_data2_heatmaps.pdf"))
plot(fit2, type = "profiles", file = file.path(out, "GCMS_data2_profiles_and_spectra.pdf"))

# For a full chromatogram with a separate calibration standard:
beer <- GCMSdeco("beer.cdf", "beer.msp", alkanes_file = "alkanes.cdf",
  xlsx_file = "beer.xlsx", nwindows = 10L, precluster = 32L)
}
