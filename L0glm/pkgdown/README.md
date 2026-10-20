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
* `man/figures/logo.png` is the source hex logo.

The GitHub Actions workflow builds the site on pushes to `master`. Its deployment
job is enabled with the repository variable `PAGES_ENABLED=true` after GitHub
Pages has been enabled for the chosen publication repository. Pages currently
requires a public repository or a GitHub plan supporting private-repository Pages.
