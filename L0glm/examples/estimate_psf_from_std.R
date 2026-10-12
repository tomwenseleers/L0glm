data(GCMS_data1)
# total ion chromatogram of C9 to C40 alkane ladder
tic <- attr(GCMS_data1, "alkanes_TIC")

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
matplot(fit$X[seq(1, nrow(fit$X), by = 1), seq(1, ncol(fit$X), by = 300)], 
        type = "l", col = "blue", lty = 1,
        xlab = "Time index", ylab = "Intensity",
        main = "Estimated time-varying PSF")

# deconvolute GC/MS chromatogram using estimated time-varying PSF:
X <- estimate_psf_from_std(
  y = tic,
  family = poisson(identity),
  baseline_df = 10,
  output_timepoints = as.numeric(rownames(GCMS_data1)),
  n_peaks = 32,
  return_details = FALSE
)
matplot(X, 
        type = "l", col = "blue", lty = 1,
        xlab = "Time index", ylab = "Intensity",
        main = "Estimated time-varying PSF\nfor chromatogram subsection")
system.time(fit.mc <- L0nnpois(
  X,
  GCMS_data1,
  lambda = "gic",
  grouped = TRUE,
  relaxed = TRUE,
  parallelize = TRUE,
  ncores = round(parallel::detectCores(logical = FALSE))
)) # 1.4s

plotmat(GCMS_data1, 
        xlab = "Elution time (min)",
        ylab = "Mass (m/z)",
        main = "Section of Cataglyphis iberica\nant queen cuticular GC/MS profile")
sum(Matrix::rowSums(fit.mc$coefficients)!=0) # 63
rows = Matrix::rowSums(fit.mc$coefficients)[which(Matrix::rowSums(fit.mc$coefficients)!=0)]
nrows = rows/max(rows)
abline(v = which(Matrix::rowSums(fit.mc$coefficients)!=0),
       col = sapply(nrows, function (alpha) adjustcolor("cyan", alpha = alpha)))
# TODO: cluster together compounds with similar mass spectra when joined elution profiles
# are log concave

hist(Matrix::rowSums(fit.mc$coefficients)[which(Matrix::rowSums(fit.mc$coefficients)!=0)])

for (i in 1:length(which(Matrix::rowSums(fit.mc$coefficients)!=0)) ) {
  plot(x = as.numeric(colnames(GCMS_data1)),
       y = fit.mc$coefficients[which(Matrix::rowSums(fit.mc$coefficients)!=0)[i],],
       main = paste0("Spectrum deconvoluted compound ", as.character(i),
                     "\nRT ", rownames(GCMS_data1)[which(Matrix::rowSums(fit.mc$coefficients)!=0)[1]], " min"),
       ylab = "Intensity", xlab = "Mass",
       col = "blue", type = "h")
}

