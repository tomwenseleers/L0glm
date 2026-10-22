# Run from the repository root: Rscript L0glm/pkgdown/build.R
package <- if (file.exists("L0glm/DESCRIPTION")) "L0glm" else "."
assets <- file.path(package, "pkgdown/assets")
dir.create(assets, recursive = TRUE, showWarnings = FALSE)
logos <- c("logo.png", "logo-dark.svg", "logo-light.svg", "logo-dark.png", "logo-light.png")
file.copy(file.path(package, "man/figures", logos), assets, overwrite = TRUE)
figures <- c("examples-single-mirrored-1.png", "examples-multichannel-heatmaps-1.png",
  "examples-multichannel-heatmaps-2.png", "examples-multichannel-coefficients-1.png")
names <- c("deconvolution-single.png", "deconvolution-multi-l0glm.png",
  "deconvolution-multi-lasso.png", "deconvolution-coefficients.png")
for (i in seq_along(figures)) file.copy(file.path(package, "../benchmarks/plots/readme", figures[i]),
  file.path(assets, names[i]), overwrite = TRUE)
website_figures <- file.path(package, "../benchmarks/plots/website")
file.copy(list.files(website_figures, pattern="\\.png$", full.names=TRUE),
  assets, overwrite=TRUE)
pkgdown::build_site(package, examples = FALSE, install = FALSE,
  new_process = FALSE, preview = FALSE, quiet = FALSE)
destination <- pkgdown::as_pkgdown(package)$dst_path
# Remove the automatic home header after logo insertion, without freeing nodes.
home <- file.path(destination, "index.html")
html <- xml2::read_html(home)
xml2::xml_remove(xml2::xml_find_all(html, ".//main//div[@class='page-header']"), free = FALSE)
xml2::write_html(html, home)
pkgdown::build_search(package)
file.copy(file.path(package, "pkgdown/extra.css"),
  file.path(destination, "extra.css"), overwrite = TRUE)
