data(simulated_spectra)
sim <- simulate_GCMS(n = 10000,
                     k = 500,
                     S = simulated_spectra,
                     peak_width_range = c(3, 5),
                     mean_beta  = 1e4,
                     sd_logbeta = 1,
                     width_height_coef = 3,
                     rt_height_coef    = 0.02,
                     sample_nr = 1,
                     rt_drift_coef = 0.001,
                     bg_scale = 100,
                     family = "poisson")
cols <- colorRampPalette(c("black", "black","blue","green","yellow","orange", "red"))(1000)
Dvis <- sim$D[seq(1, 10000, by = 3), , drop = FALSE]^0.12
image(x = as.numeric(rownames(Dvis)),
      y = as.numeric(colnames(Dvis)),
      z = Dvis,
      col    = cols,
      useRaster = TRUE,
      xlab   = "Time (min)",
      ylab   = "Mass (m/z)",
      main   = "Simulated GC/MS chromatogram")
