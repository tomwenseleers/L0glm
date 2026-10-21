# Building the documentation website

From the repository root, run:

```r
source("L0glm/pkgdown/build.R")
```

The site is generated in `docs/`. Reference examples are displayed as source;
the tutorials evaluate their compact runnable examples. The million-predictor
benchmark is documented with its complete code and recorded output, so a site
build does not allocate that large design.

* `_pkgdown.yml` defines navigation and the reference index.
* `pkgdown/index.md` provides the home page.
* `pkgdown/extra.css` contains the visual theme.
* `vignettes/` contains tutorials and the algorithm/theory article.
* `man/figures/` contains the primary logo and its light and dark variants.
* `pkgdown/build_favicons.R` regenerates browser and app icons from the logo.

The website is published at <https://tomwenseleers.github.io/L0glm/>.
The GitHub Actions workflow builds and publishes it on pushes to `master`.
GitHub Pages uses the Actions deployment source, with the repository variable
`PAGES_ENABLED=true`. KaTeX renders mathematical notation in tutorials and help
pages.
