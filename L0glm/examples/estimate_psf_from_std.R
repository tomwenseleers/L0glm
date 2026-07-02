data(GCMS_data1)
tic <- attr(GCMS_data1, "alkanes_TIC") # total ion chromatogram of C9 to C40 alkane ladder

plot(x=as.numeric(names(attr(GCMS_data1, "alkanes_TIC"))),  
           y=attr(GCMS_data1, "alkanes_TIC"), 
           xlab="Retention time (min)", 
           ylab="Total ion intensity", 
           main="C9 to C40 alkane ladder standard", 
           col = "blue",
           type="l")

# estimate time-varying PSF from alkane ladder standard:
fit <- estimate_psf_from_std(
  y = tic,
  family = poisson(identity),
  baseline_df = 10,
  output_timepoints = as.numeric(names(tic)),
  n_peaks = 32,
  return_details = TRUE
)

# fitted peaks:
matplot(fit$psf_offsets, t(fit$psf_kernels[1:32, ]),
        type = "l", lty = 1,
        xlab = "Offset from centroid (min)",
        ylab = "PSF density",
        main = "Fitted peak shapes")

# estimated time-varying PSF:
matplot(fit$X[seq(1, nrow(fit$X), by = 1), seq(1, ncol(fit$X), by = 200)], 
        type = "l", col = "blue", lty = 1,
        xlab = "Time index", ylab = "Intensity",
        main = "Estimated time-varying PSF")
