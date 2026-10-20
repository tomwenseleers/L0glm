# Run from the repository root: Rscript L0glm/pkgdown/build.R
package <- if (file.exists("L0glm/DESCRIPTION")) "L0glm" else "."
assets <- file.path(package, "pkgdown/assets")
dir.create(assets, recursive = TRUE, showWarnings = FALSE)
file.copy(file.path(package, "man/figures/logo.png"), file.path(assets, "logo.png"), overwrite = TRUE)
figures <- c("examples-single-mirrored-1.png", "examples-multichannel-heatmaps-1.png",
  "examples-multichannel-heatmaps-2.png", "examples-multichannel-coefficients-1.png")
names <- c("deconvolution-single.png", "deconvolution-multi-l0glm.png",
  "deconvolution-multi-lasso.png", "deconvolution-coefficients.png")
for (i in seq_along(figures)) file.copy(file.path(package, "../benchmarks/plots/readme", figures[i]),
  file.path(assets, names[i]), overwrite = TRUE)
pkgdown::build_site(package, examples = FALSE, install = FALSE,
  new_process = FALSE, preview = FALSE, quiet = FALSE)
file.copy(file.path(package, "pkgdown/extra.css"),
  file.path(pkgdown::as_pkgdown(package)$dst_path, "extra.css"), overwrite = TRUE)
