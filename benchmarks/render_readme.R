# Regenerate README.md and its figures from README.Rmd with the installed package.
# Run from the repository root or the L0glm/ package directory.
render_readme <- function() {
  root <- if (file.exists("README.Rmd")) "." else ".."
  if (!file.exists(file.path(root, "README.Rmd")))
    stop("Run from the repository root or L0glm/ package directory.")
  previous <- setwd(root)
  on.exit(setwd(previous))
  if (!requireNamespace("knitr", quietly=TRUE)) stop("Install knitr to render the README.")
  library(L0glm)
  # Consistent BLAS settings; screening explicitly requests eight threads.
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
  knitr::knit("README.Rmd", output="README.md", envir=new.env(parent=globalenv()))
  rendered <- sub("[ \t]+$", "", readLines("README.md"))
  writeLines(gsub("\n{3,}", "\n\n", paste(rendered, collapse="\n")), "README.md")
}
render_readme()
