"""Publish the same measured benchmark tables in the website and README."""
import csv
import re
from statistics import median
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "benchmarks/plots/website"
LABELS = {"L0glm": "L0glm", "lasso_bic": "LASSO", "relaxed_lasso_bic": "Relaxed LASSO",
          "lasso_gic": "LASSO", "relaxed_lasso_gic": "Relaxed LASSO",
          "lasso_cv": "LASSO", "relaxed_lasso_cv": "Relaxed LASSO",
          "abess": "ABESS", "lasso": "LASSO", "relaxed_lasso": "Relaxed LASSO",
          "adaptive_lasso": "Adaptive LASSO", "adaptive_relaxed_lasso": "Adaptive relaxed LASSO"}

def rows(name):
    with (DATA / name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))

def value(row, field):
    raw = row[field]
    if field == "method":
        return LABELS.get(raw, raw)
    if field == "tune":
        return raw
    x = float(raw)
    if field == "seconds":
        return f"{x:.2f}"
    if field in {"F1_score", "MCC", "within_row_F1_score", "within_row_MCC"}:
        return f"{x:.3f}"
    return f"{x:.0f}"

def table(data, fields, bold_first=False):
    names = {"method": "Method", "tune": "Selection", "seconds": "Seconds", "gic": "GIC",
             "F1_score": "F1", "within_row_F1_score": "Within-row F1",
             "within_row_MCC": "Within-row MCC"}
    lines = ["| " + " | ".join(names.get(f, f) for f in fields) + " |",
             "| " + " | ".join("---" if f in {"method", "tune"} else "---:" for f in fields) + " |"]
    for i, row in enumerate(data):
        values = [value(row, f) for f in fields]
        if bold_first and i == 0:
            values = [f"**{v}**" for v in values]
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)

million = rows("million-metrics.csv")[0]
million_lasso = rows("million-lasso-metrics.csv")
million_display = [{**million, "method": "L0glm + ISIS", "tune": "GIC preset"}] + million_lasso
million_cv = {row["method"]: row for row in million_lasso if row["tune"] == "5-fold CV"}
single = rows("single-metrics.csv")
multi = rows("multi-metrics.csv")
selection = rows("selection-comparison.csv")
multi_by_method = {row["method"]: row for row in multi}
multi_relaxed_ratio = float(multi_by_method["relaxed_lasso_gic"]["seconds"]) / float(multi_by_method["L0glm"]["seconds"])
fields = ["method", "tune", "seconds", "gic", "TP", "FP", "FN", "TN", "F1_score", "MCC"]
ratio = float(million_cv["lasso_cv"]["seconds"]) / float(million["seconds"])
def group_label(raw):
    label = LABELS[raw].replace("LASSO", "group LASSO")
    return label[0].upper() + label[1:]


multi_display = [{**r, "method": group_label(r["method"])} for r in multi]


def deconvolution_speed(data, grouped=False):
    fitted = data[0]
    by_method = {row["method"]: row for row in data}
    elapsed = float(fitted["seconds"])
    plain_ratio = float(by_method["lasso_gic"]["seconds"]) / elapsed
    relaxed_ratio = float(by_method["relaxed_lasso_cv"]["seconds"]) / elapsed
    name = "group LASSO" if grouped else "LASSO"
    return (f'L0glm recovers **{int(fitted["TP"]):,} of {int(fitted["support_size_true"]):,} true peaks '
            f'with {int(fitted["FP"]):,} false positives** in **{elapsed:.2f} seconds**. '
            f'In this benchmark it is **{plain_ratio:.1f}× faster than GIC-tuned {name}** and '
            f'**{relaxed_ratio:.0f}× faster than relaxed {name} tuned by five-fold CV**.')


def abess_table(model):
    methods = ("L0glm", "abess", "L0Learn", "lasso", "relaxed_lasso",
               "adaptive_lasso", "adaptive_relaxed_lasso", "MCP", "SCAD")
    data = { (r["rho"], r["method"]): r for r in rows("abess-simulations.csv")
             if r["model"] == model }
    prediction = "Test RMSE" if model == "gaussian" else "Test AUC"
    lines = [f"| Correlation | Method | {prediction} | Sensitivity | FPR | MCC | Seconds |",
             "|:--|:--|--:|--:|--:|--:|--:|"]
    for rho in ("0.1", "0.7"):
        for method in methods:
            r = data[rho, method]
            lines.append(f'| {rho} | {LABELS.get(method, method)} | '
                         f'{float(r["prediction"]):.3f} | {100*float(r["TPR"]):.0f}% | '
                         f'{100*float(r["FPR"]):.2f}% | {float(r["MCC"]):.3f} | '
                         f'{float(r["seconds"]):.2f} |')
    return "\n".join(lines)


def local_search_table():
    data = { (r["model"], r["rho"], r["variant"]): r for r in rows("local-search.csv") }
    lines = ["| Family | Correlation | Local search | GIC | TP | FP | MCC | Seconds |",
             "|:--|--:|:--|--:|--:|--:|--:|--:|"]
    for model, label in (("gaussian", "Gaussian"), ("binomial", "Logistic")):
        for rho in ("0.1", "0.7"):
            for variant in ("disabled", "auto"):
                r = data[model, rho, variant]
                lines.append(f'| {label} | {rho} | {"FALSE" if variant == "disabled" else "auto"} | '
                             f'{float(r["GIC"]):.1f} | {float(r["TP"]):.1f} | '
                             f'{float(r["FP"]):.1f} | {float(r["MCC"]):.3f} | '
                             f'{float(r["seconds"]):.2f} |')
    return "\n".join(lines)


def difficulty_table(model):
    methods = (("L0glm", "gic"), ("abess", "gic"), ("L0Learn", "gic"),
               ("lasso", "gic"), ("relaxed_lasso", "gic"),
               ("adaptive_lasso", "gic"), ("adaptive_relaxed_lasso", "gic"),
               ("MCP", "gic"), ("SCAD", "gic"), ("lasso", "cv"), ("relaxed_lasso", "cv"))
    data = {(int(r["p"]), r["method"], r["tune"]): r
            for r in rows("support-difficulty.csv") if r["model"] == model}
    prediction = "Test RMSE" if model == "gaussian" else "Test AUC"
    lines = [f"| Predictors | Method | Selection | {prediction} | TP | FP | F1 | MCC | Seconds |",
             "|--:|:--|:--|--:|--:|--:|--:|--:|--:|"]
    for p in sorted({key[0] for key in data}):
        for method, tune in methods:
            r = data[p, method, tune]
            selection = "5-fold CV" if tune == "cv" else "GIC preset" if method == "L0glm" else "GIC path"
            lines.append(f'| {p:,} | {LABELS.get(method, method)} | {selection} | '
                         f'{float(r["prediction"]):.3f} | {float(r["TP"]):.1f} | '
                         f'{float(r["FP"]):.1f} | {float(r["F1"]):.3f} | '
                         f'{float(r["MCC"]):.3f} | {float(r["seconds"]):.2f} |')
    return "\n".join(lines)


def local_coordinate_table():
    data = rows("local-search-million-runs.csv")
    lines = ["| Search coordinates | Total seconds | Local-search seconds | GIC | TP | FP | F1 | MCC |",
             "|:--|--:|--:|--:|--:|--:|--:|--:|"]
    for coordinates, label in (("original", "Original predictor space"),
                               ("screened", "Retained screening pool")):
        runs = [r for r in data if r["coordinates"] == coordinates]
        r = runs[0]
        lines.append(f'| {label} | {median(float(x["seconds"]) for x in runs):.2f} | '
                     f'{median(float(x["local_seconds"]) for x in runs):.2f} | '
                     f'{float(r["gic"]):.1f} | {r["TP"]} | {r["FP"]} | '
                     f'{float(r["F1_score"]):.3f} | {float(r["MCC"]):.3f} |')
    return "\n".join(lines)


content = {
    "million-intro": (
        "This Gaussian simulation has **one million predictors**, 500 observations and "
        f'{million["support_size_true"]} nonzero effects. One default L0glm call takes '
        f'**{float(million["seconds"]):.2f} seconds** and recovers '
        f'**{million["TP"]} true effects with {million["FP"]} false positives**.'
    ),
    "million-screening-summary": (
        f'The default L0glm fit with ISIS recovers **{million["TP"]} true effects with '
        f'{million["FP"]} false positives** in **{float(million["seconds"]):.2f} s** '
        'using eight physical cores. Five-fold-CV LASSO on the unscreened design takes '
        f'**{float(million_cv["lasso_cv"]["seconds"]):.2f} s** and recovers '
        f'{million_cv["lasso_cv"]["TP"]} true effects with {million_cv["lasso_cv"]["FP"]} false positives. '
        'The [benchmark article](performance.html#one-million-predictor-example) gives the tuning settings.'
    ),
    "million-comparison": table(million_display, fields, bold_first=True),
    "million-speed": f'L0glm is **{ratio:.0f}× faster than five-fold-CV LASSO** in this comparison, with higher sensitivity and fewer false positives.',
    "million-metrics": table([million], fields[2:]),
    "single-metrics": table(single, fields),
    "multi-metrics": table(multi_display, fields + ["within_row_F1_score", "within_row_MCC"]),
    "multi-home": table(multi_display, ["method","tune","seconds","gic","TP","FP","F1_score","MCC"]),
    "multi-spectral": (
        'Among correctly located peaks, the recovered mass spectra have a within-row F1 of '
        f'**{float(multi[0]["within_row_F1_score"]):.3f}** and MCC of '
        f'**{float(multi[0]["within_row_MCC"]):.3f}**.'
    ),
    "single-speed": deconvolution_speed(single),
    "multi-speed": deconvolution_speed(multi, grouped=True),
    "selection-metrics": table(selection, ["method", "gic", "TP", "FP", "FN", "TN", "F1_score", "MCC"]),
    "speed-stat": f'<a href="articles/performance.html#one-million-predictor-example"><strong>{ratio:.0f}× faster</strong><span>than five-fold-CV LASSO in that example</span></a>',
    "multi-stat": f'<a href="articles/performance.html#multichannel-deconvolution"><strong>{multi[0]["TP"]} / {multi[0]["support_size_true"]}</strong><span>peaks recovered in simulated GC/MS deconvolution and <b>{multi_relaxed_ratio:.0f}× faster</b> than GIC-tuned relaxed group LASSO</span></a>',
    "abess-gaussian": abess_table("gaussian"),
    "abess-binomial": abess_table("binomial"),
    "local-search-comparison": local_search_table(),
    "difficulty-gaussian": difficulty_table("gaussian"),
    "difficulty-binomial": difficulty_table("binomial"),
    "local-coordinate-comparison": local_coordinate_table(),
}
paths = [ROOT/"README.Rmd", ROOT/"README.md", ROOT/"L0glm/pkgdown/index.md", ROOT/"L0glm/vignettes/L0glm.Rmd"]
paths += list((ROOT/"L0glm/vignettes/articles").glob("*.Rmd"))
for path in paths:
    text = original = path.read_text(encoding="utf-8")
    for key, replacement in content.items():
        pattern = re.compile(r"(<!-- BEGIN " + re.escape(key) + r" -->).*?(<!-- END " + re.escape(key) + r" -->)", re.S)
        text = pattern.sub(lambda m: m[1] + "\n" + replacement + "\n" + m[2], text)
    if path.name == "performance.Rmd":
        compared = [million, million_cv["lasso_cv"], million_cv["relaxed_lasso_cv"]]
        text = re.sub(r"seconds <- c\([^\n]+", 'seconds <- c(' + ', '.join(f'{float(row["seconds"]):.2f}' for row in compared) + ')', text)
        for field in ("TP", "FP"):
            text = re.sub(rf'{field.lower()} <- c\([^\n]+', field.lower() + ' <- c(' + ', '.join(row[field] for row in compared) + ')', text)
        text = re.sub(r"\*\*\d+ times faster than five-fold-CV LASSO\*\*", f'**{ratio:.0f} times faster than five-fold-CV LASSO**', text)
    if text != original:
        path.write_text(text, encoding="utf-8")
print("Updated measured tables and homepage statistics.")
