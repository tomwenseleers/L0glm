# L0glm visual identity

The hexagon contains a sparse coefficient matrix with selected predictor rows
highlighted across response channels. Turquoise and amber rows represent best
subset selection for individual predictors and groups of predictors.

The artwork comes from the supplied `L0glm_hex_logo_assets.zip` design set.

| Asset | Use |
|---|---|
| `man/figures/logo.png` | Primary package logo: supplied dark PNG, 480 × 556 |
| `man/figures/logo-dark.svg` | Dark backgrounds and the website banner |
| `man/figures/logo-light.svg` | White and light backgrounds |
| `man/figures/logo-dark.png`, `logo-light.png` | Transparent PNGs, 1200 × 1400 |
| `man/figures/logo-dark-editable.svg`, `logo-light-editable.svg` | Editable text versions |

The display SVGs have lettering converted to vector outlines. The README
selects the light or dark SVG using the viewer's colour scheme. Browser and
app icons use the primary dark design.

Regenerate icons from the repository root with
`source("L0glm/pkgdown/build_favicons.R")`, then rebuild the site with
`source("L0glm/pkgdown/build.R")`. Commit the source assets and generated `docs/`
alongside website changes.

| Role | Colour |
|---|---|
| Navy | `#0A182D` |
| Turquoise | `#66E8D7` |
| Teal | `#14B6AA` |
| Amber | `#F3A64D` |
| Light mint | `#DCF7F2` |

Website body text uses system fonts; mathematical notation and code remain
selectable. The layout supports keyboard navigation, responsive tables and
reduced-motion preferences.
