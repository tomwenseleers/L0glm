# Regenerate README.md with the website's published benchmark tables and figures.
# Run from the repository root or the L0glm/ package directory.
render_readme <- function() {
  root <- if (file.exists("README.Rmd")) "." else ".."
  if (!file.exists(file.path(root, "README.Rmd")))
    stop("Run from the repository root or L0glm/ package directory.")
  previous <- setwd(root)
  on.exit(setwd(previous))
  if (!requireNamespace("knitr", quietly=TRUE)) stop("Install knitr to render the README.")
  status <- system2("python", "benchmarks/update_website_tables.py")
  if (status != 0L) stop("Could not synchronize the published benchmark tables.")
  knitr::knit("README.Rmd", output="README.md", envir=new.env(parent=globalenv()))
  rendered <- sub("[ \t]+$", "", readLines("README.md"))
  writeLines(gsub("\n{3,}", "\n\n", paste(rendered, collapse="\n")), "README.md")
}
render_readme()
