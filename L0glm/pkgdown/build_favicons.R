# Run from the repository root after updating man/figures/logo-dark.png.
package <- if (file.exists("L0glm/DESCRIPTION")) "L0glm" else "."
destination <- file.path(package, "pkgdown/favicon")
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
logo <- magick::image_read(file.path(package, "man/figures/logo-dark.png"))
icon <- function(size) magick::image_extent(
  magick::image_resize(logo, sprintf("%dx%d", size, size)),
  sprintf("%dx%d", size, size), gravity = "center", color = "none")
sizes <- c("favicon-96x96.png" = 96L, "apple-touch-icon.png" = 180L,
  "web-app-manifest-192x192.png" = 192L, "web-app-manifest-512x512.png" = 512L)
for (name in names(sizes)) magick::image_write(icon(sizes[[name]]),
  file.path(destination, name), format = "png")
magick::image_write(magick::image_join(lapply(c(16L, 32L, 48L, 64L, 128L, 256L), icon)),
  file.path(destination, "favicon.ico"), format = "ico")
file.copy(file.path(package, "man/figures/logo-dark.svg"),
  file.path(destination, "favicon.svg"), overwrite = TRUE)
