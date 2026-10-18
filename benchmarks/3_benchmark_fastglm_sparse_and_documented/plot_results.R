# Run after the two benchmark scripts. No additional plotting packages required.
repo <- if (file.exists("L0glm/DESCRIPTION")) "." else ".."
root <- file.path(repo, "benchmarks/3_benchmark_fastglm_sparse_and_documented")
grid <- read.csv(file.path(root, "sparse_truth/timings_metrics.csv"))
docs <- read.csv(file.path(root, "documented_examples/timings_metrics.csv"))
cases <- read.csv(file.path(root, "sparse_truth/cases.csv"))
families <- c("gaussian_identity", "binomial_logit", "binomial_probit", "binomial_cloglog",
              "poisson_log", "Gamma_log", "inverse_gaussian_log", "nb_log")
family_labels <- c("Gaussian\nidentity", "Binomial\nlogit", "Binomial\nprobit", "Binomial\ncloglog",
                   "Poisson\nlog", "Gamma\nlog", "Inverse\nGaussian\nlog", "Negative\nbinomial\nlog")
matrix_by_case <- function(data, value) {
  mat <- matrix(NA_real_, length(families), nrow(cases))
  for (i in seq_len(nrow(data))) mat[match(data$family[i], families), match(data$case[i], cases$case)] <- value[i]
  mat
}
selected <- grid[grid$method == "L0_gic", ]
reference <- grid[grid$method == "fastglmPure", ]
paired <- merge(selected, reference, by = c("case", "family"), suffixes = c("_L0", "_fastglm"))
ratio <- paired$seconds_L0/paired$seconds_fastglm
ratio[paired$p_L0 + 1 >= paired$n_L0 | !is.na(paired$error_fastglm)] <- NA_real_
cost <- matrix_by_case(paired, log2(ratio))
quality <- matrix_by_case(selected, selected$F1)
heatmap <- function(mat, title, colors, limits, labels) {
  par(mar = c(6, 13, 4, 1))
  plot(0, type = "n", xlim = c(.5, length(families)+.5), ylim = c(.5, nrow(cases)+.5),
       axes = FALSE, xlab = "", ylab = "", main = title)
  for (j in seq_len(ncol(mat))) for (i in seq_len(nrow(mat))) {
    value <- mat[i,j]
    index <- if (is.na(value)) NA_integer_ else 1L + round((length(colors)-1) *
      (min(limits[2], max(limits[1], value))-limits[1])/diff(limits))
    rect(i-.5, j-.5, i+.5, j+.5, col = if (is.na(index)) "grey90" else colors[index], border = "white")
    if (is.finite(value)) text(i, j, labels(value), cex = .68)
  }
  axis(1, seq_along(families), labels = FALSE, tick = FALSE)
  mtext(family_labels, side = 1, at = seq_along(families), line = 2, cex = .7)
  axis(2, seq_len(nrow(cases)), cases$case, las = 1, tick = FALSE, cex.axis = .8)
  box()
}
render <- function(path, draw, width = 14, height = 9) {
  png(file.path(root, paste0(path, ".png")), width = width, height = height, units = "in", res = 180)
  draw(); dev.off()
  pdf(file.path(root, paste0(path, ".pdf")), width = width, height = height)
  draw(); dev.off()
}
render("sparse_truth_comparison", function() {
  par(mfrow = c(1,2), oma = c(3,0,0,0))
  heatmap(cost, "L0 GIC / fastglmPure elapsed time", colorRampPalette(c("#2166ac","white","#b2182b"))(101),
    c(-4,4), function(x) sprintf("%.1fx", 2^x))
  heatmap(quality, "L0 GIC support recovery (F1)", colorRampPalette(c("#fff7bc","#238443"))(101),
    c(0,1), function(x) sprintf("%.2f", x))
  mtext("Blue: L0 faster; red: fastglmPure faster. Grey: no comparable unique full-design reference or failed fit.",
    side = 1, outer = TRUE, line = 1, cex = .8)
}, width = 18, height = 10)
methods <- c("L0_unpenalized", "L0_default", "L0_gic", "abess_default")
labels <- c("L0 lambda=0", "L0 lambda=1", "L0 GIC", "abess defaults")
colors <- c("#377eb8", "#984ea3", "#168a33", "#e41a1c")
dense_cases <- c("readme_dense_logit_5000", "readme_dense_logit_50000", "vignette_factor_dense",
  "large_data_factor_dense", "additional_dense_gaussian_421", "vignette_nb_challenging")
dense_names <- c("README\nn=5000,p=30", "Scaled README\nn=50000,p=300", "Factor interactions\nn=50000,p=199",
  "Factor interactions\nn=5000,p=99", "Additional Gaussian\nn=400,p=100", "Negative binomial\nn=200,p=8")
fractions <- matrix(NA_real_, length(methods), length(dense_cases))
for (j in seq_along(dense_cases)) for (i in seq_along(methods)) {
  row <- docs[docs$case == dense_cases[j] & docs$method == methods[i], ]
  if (nrow(row)) fractions[i,j] <- row$support_fraction[1]
}
render("dense_truth_support", function() {
  par(mar = c(7,5,4,1))
  barplot(fractions, beside = TRUE, col = colors, border = NA, names.arg = dense_names,
    ylim = c(0,1.13), ylab = "Fraction of slopes retained", main = "Dense truth: selected support depends on signal and penalty", cex.names = .8)
  abline(h = 1, lty = 3, col = "grey40")
  legend("top", legend = labels, fill = colors, horiz = TRUE, bty = "n", cex = .85)
  mtext("Blank bars indicate unsupported, omitted or timed-out comparisons; see the CSV error column. Intercepts are excluded.",
    side = 1, line = 5, cex = .8)
}, width = 13, height = 7)
cat("Benchmark figures written to", root, "\n")
