"""Publish the same measured benchmark tables in the website and README."""
import csv
import re
from math import log
from statistics import median
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "benchmarks/plots/website"
LABELS = {"L0glm": "L0glm", "biglasso": "biglasso LASSO", "lasso_bic": "LASSO", "relaxed_lasso_bic": "Relaxed LASSO",
          "lasso_gic": "LASSO", "relaxed_lasso_gic": "Relaxed LASSO",
          "lasso_cv": "LASSO", "relaxed_lasso_cv": "Relaxed LASSO",
          "abess": "ABESS", "lasso": "LASSO", "relaxed_lasso": "Relaxed LASSO",
          "adaptive_lasso": "Adaptive LASSO", "adaptive_relaxed_lasso": "Adaptive relaxed LASSO"}

def rows(name):
    return read_csv(DATA / name)


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))

def value(row, field):
    raw = row[field]
    if field == "method":
        return LABELS.get(raw, raw)
    if field in {"tune", "screen", "model"}:
        return raw
    if field == "interaction_recovered":
        return "Yes" if raw == "TRUE" else "No"
    if raw in {"NA", "NaN", ""}:
        return "—"
    x = float(raw)
    if field == "seconds":
        return f"{x:.2f}"
    if field in {"F1_score", "MCC", "within_row_F1_score", "within_row_MCC"}:
        return f"{x:.3f}"
    if field in {"test_RMSE", "test_AUC", "test_logloss"}:
        return f"{x:.3f}"
    return f"{x:.0f}"

def table(data, fields, bold_first=False):
    names = {"method": "Method", "tune": "Selection", "seconds": "Seconds", "gic": "GIC",
             "F1_score": "F1", "within_row_F1_score": "Within-row F1",
             "within_row_MCC": "Within-row MCC", "model": "Complete model",
             "screen": "Screen", "p": "Predictors", "test_RMSE": "Test RMSE",
             "test_AUC": "Test AUC", "test_logloss": "Test log loss",
             "interaction_recovered": "Interaction recovered"}
    lines = ["| " + " | ".join(names.get(f, f) for f in fields) + " |",
             "| " + " | ".join("---" if f in {"method", "tune", "screen", "model"} else "---:" for f in fields) + " |"]
    for i, row in enumerate(data):
        values = [value(row, f) for f in fields]
        if bold_first and i == 0:
            values = [f"**{v}**" for v in values]
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)

million = rows("million-metrics.csv")[0]
million_lasso = rows("million-lasso-metrics.csv")
assert int(million["n"]) == 1000 and int(million["support_size_true"]) == 60
assert all(int(r["n"]) == int(million["n"]) and int(r["p"]) == int(million["p"])
           and int(r["support_size_true"]) == 60 and r["status"] == "ok" for r in million_lasso)
million_display = [{**million, "method": "L0glm + ISIS", "tune": "GIC preset"}] + million_lasso
million_cv = {row["method"]: row for row in million_lasso if row["tune"] == "5-fold CV"}
single = rows("single-metrics.csv")
multi = rows("multi-metrics.csv")
selection = rows("selection-comparison.csv")
single_by_method = {row["method"]: row for row in single}
multi_by_method = {row["method"]: row for row in multi}
multi_relaxed_ratio = float(multi_by_method["relaxed_lasso_gic"]["seconds"]) / float(multi_by_method["L0glm"]["seconds"])
fields = ["method", "tune", "seconds", "gic", "TP", "FP", "FN", "TN", "F1_score", "MCC"]
highlight_fields = [f for f in fields if f not in {"FN", "TN"}]
ratio = float(million_cv["lasso_cv"]["seconds"]) / float(million["seconds"])
tuning = rows("million-tuning.csv")[0]
million_false = ("no false positives" if int(million["FP"]) == 0
                 else f'{million["FP"]} false positives')


def screening_tables():
    directory = ROOT / "benchmarks/3_benchmarks_sis_vs_SIS_results"
    screening = read_csv(directory / "timings.csv")
    labels = {"sis_glm_auto": "`sisGLM()`, eight physical cores",
              "sis_glm_1": "`sisGLM(nthreads=1)`", "SIS_screen": "SIS Gaussian screen (`obtain.ix0`)"}
    speeds = {method: median(float(r["seconds"]) for r in screening if r["method"] == method)
              for method in labels}
    screen_table = table([{"method": label, "seconds": speeds[method]}
                          for method, label in labels.items()], ["method", "seconds"])
    models = read_csv(directory / "models.csv")
    assert all(int(r["n"]) == 1000 for r in models)
    model_labels = {"SIS_package": "SIS::SIS + BIC LASSO", "l0glmFit": "l0glmFit, GIC preset",
                    "sparseglm_lasso": "sparseGLM LASSO, BIC grid",
                    "sparseglm_relaxed_lasso": "sparseGLM relaxed LASSO, BIC grid"}
    complete = []
    for key in dict.fromkeys((r["p"], r["screen"], r["method"]) for r in models):
        group = [r for r in models if (r["p"], r["screen"], r["method"]) == key]
        complete.append({**group[0], "model": model_labels[key[2]],
                         "seconds": median(float(r["seconds"]) for r in group)})
    complete_table = table(complete,["p", "screen", "model", "seconds", "gic", "TP", "FP", "F1_score", "MCC"])
    adapters = read_csv(directory / "all_methods.csv")
    assert all(int(r["n"]) == 1000 and r["status"] == "ok" for r in adapters)
    adapter_table = table([r for r in adapters if r["screen"] == "ISIS"],
                          ["method", "seconds", "gic", "TP", "FP", "F1_score", "MCC"])
    (directory / "README.md").write_text(
        "# SIS screening and complete-model comparisons\n\n"
        "The Gaussian simulation has 1,000 observations, one million standard-normal predictors "
        "and 60 signed nonzero effects, matching the homepage exactly. Seeds are 123 for the "
        "design/effects and 1 for the noise. All timing medians use three repetitions and exclude "
        "simulation. Hardware: HP ZBook Fury G8, Intel i9-11950H, eight physical cores, "
        "16 logical processors, Windows, R 4.5.0 and one BLAS thread.\n\n"
        "## Screening only\n\n" + screen_table + "\n\n"
        "All complete predictor rankings agree. SIS's `obtain.ix0` helper excludes its "
        "subsequent LASSO fit; Gaussian rankings are invariant to column standardization.\n\n"
        "## Complete selection pipelines\n\n" + complete_table + "\n\n"
        "One-pass SIS retains 100 predictors. Million-predictor ISIS uses the default "
        "GIC batch and ten rounds; the SIS-package ISIS comparison uses 5,000 predictors "
        "and three rounds. `nsis=100` is a batch size for L0glm and a total support cap "
        "for SIS. These comparisons set `local.search=FALSE` for standalone L0glm, "
        "isolating screening and adaptive ridge. SIS fits a free intercept; L0glm "
        "and glmnet use the simulation's zero intercept.\n\n"
        "## All supported adapters with ISIS\n\n" + adapter_table + "\n\n"
        "L0glm uses the GIC preset; competitors select BIC paths on their own screening "
        "pools. The adapter table reports one complete fitting/scoring run per method. "
        "The CSV files retain every measured run, criterion and setting. Reproduce with "
        "[3_benchmarks_sis_vs_SIS.R](../3_benchmarks_sis_vs_SIS.R) and "
        "[3_benchmarks_sparseglm_screening.R](../3_benchmarks_sparseglm_screening.R).\n",
        encoding="utf-8")
    return {"screening-times": screen_table, "screening-models": complete_table,
            "screening-adapters": adapter_table,
            "screening-speed": f'The OpenMP screen is **{speeds["SIS_screen"] / speeds["sis_glm_auto"]:.1f} times faster** on this machine.'}
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
    data = {(r["method"], r["tune"]): r for r in rows("support-difficulty.csv")
            if r["model"] == model and int(r["p"]) == 20000}
    prediction = "Test RMSE" if model == "gaussian" else "Test AUC"
    lines = [f"| Method | Selection | {prediction} | TP | FP | F1 | MCC | Seconds |",
             "|:--|:--|--:|--:|--:|--:|--:|--:|"]
    for method, tune in methods:
        r = data[method, tune]
        selection = "5-fold CV" if tune == "cv" else "GIC preset" if method == "L0glm" else "GIC path"
        lines.append(f'| {LABELS.get(method, method)} | {selection} | '
                     f'{float(r["prediction"]):.3f} | {float(r["TP"]):.1f} | '
                     f'{float(r["FP"]):.1f} | {float(r["F1"]):.3f} | '
                     f'{float(r["MCC"]):.3f} | {float(r["seconds"]):.2f} |')
    return "\n".join(lines)


def abess_highlight():
    data = {(r["model"], r["method"]): r for r in rows("support-difficulty.csv")
            if int(r["p"]) == 20000 and r["tune"] == "gic"}
    gaussian = data["gaussian", "L0glm"]
    abess = data["gaussian", "abess"]
    logistic = data["binomial", "L0glm"]
    false_positives = ("no false positives" if float(logistic["FP"]) == 0
                       else f'{float(logistic["FP"]):.1f} false positives')
    return (f'With **20,000 predictors**, L0glm matches ABESS\'s Gaussian GIC and support recovery '
            f'in **{float(gaussian["seconds"]):.2f} s versus {float(abess["seconds"]):.2f} s**. '
            f'Its logistic fit takes **{float(logistic["seconds"]):.2f} s**, recovering '
            f'**{float(logistic["TP"]):.1f} of {int(logistic["k"])} true effects** with '
            f'**{false_positives}** on average across three replications.')


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
        f'This Gaussian simulation has **one million predictors**, {int(million["n"]):,} observations and '
        f'{million["support_size_true"]} nonzero effects. One default L0glm call takes '
        f'**{float(million["seconds"]):.2f} seconds** and recovers '
        f'**{million["TP"]} true effects with {million_false}**.'
    ),
    "million-screening-summary": (
        f'The default L0glm fit with ISIS recovers **{million["TP"]} true effects with '
        f'{million_false}** in **{float(million["seconds"]):.2f} s** '
        'using eight physical cores. Five-fold-CV LASSO on the unscreened design takes '
        f'**{float(million_cv["lasso_cv"]["seconds"]):.2f} s** and recovers '
        f'{million_cv["lasso_cv"]["TP"]} true effects with {million_cv["lasso_cv"]["FP"]} false positives. '
        'The [benchmark article](performance.html#one-million-predictor-example) gives the tuning settings.'
    ),
    "million-comparison": table(million_display, fields, bold_first=True),
    "readme-million-comparison": table(
        [million_display[0], million_cv["lasso_cv"], million_cv["relaxed_lasso_cv"]],
        highlight_fields, bold_first=True),
    "million-speed": f'L0glm is **{ratio:.0f}× faster than five-fold-CV LASSO** in this comparison, recovering **{million["TP"]} true effects with {million_false}**.',
    "million-path-summary": f'The best grid value is **{float(tuning["lambda_best"]):.2f}**, with a GIC preset of **{float(tuning["lambda_preset"]):.2f}**.',
    "million-metrics": table([million], fields[2:]),
    "single-metrics": table(single, fields),
    "multi-metrics": table(multi_display, fields + ["within_row_F1_score", "within_row_MCC"]),
    "multi-home": table(multi_display, ["method","tune","seconds","gic","TP","FP","F1_score","MCC"]),
    "readme-multi-comparison": table(
        [{**multi_by_method[method], "method": group_label(method)}
         for method in ("L0glm", "relaxed_lasso_gic", "relaxed_lasso_cv")],
        highlight_fields, bold_first=True),
    "readme-single-highlight": (
        f'L0glm achieves **F1 {float(single[0]["F1_score"]):.3f} in {float(single[0]["seconds"]):.2f} s** '
        'with its GIC preset, compared with '
        f'**F1 {float(single_by_method["relaxed_lasso_cv"]["F1_score"]):.3f} '
        f'in {float(single_by_method["relaxed_lasso_cv"]["seconds"]):.2f} s** '
        'for fully relaxed LASSO tuned by five-fold CV.'
    ),
    "readme-abess-highlight": abess_highlight(),
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
content.update(screening_tables())
thread_directory = ROOT / "benchmarks/3_benchmarks_sis_openmp_results"
thread_summary = read_csv(thread_directory / "summary.csv")
thread_runs = read_csv(thread_directory / "timings.csv")
thread_table = ["| Screening threads | SIS seconds (median of 5) | SIS speedup | ISIS + L0glm seconds (median of 3) | ISIS speedup |",
                "|---:|---:|---:|---:|---:|"]
for threads in (1, 2, 4, 8):
    pair = {r["stage"]: r for r in thread_summary if int(r["threads"]) == threads}
    a, b = pair["SIS"], pair["ISIS + L0glm"]
    thread_table.append(f'| {threads} | {float(a["median_seconds"]):.2f} | {float(a["speedup"]):.2f} | '
                        f'{float(b["median_seconds"]):.2f} | {float(b["speedup"]):.2f} |')
thread_fit = next(r for r in thread_runs if r["stage"] == "ISIS + L0glm")
assert abs(float(thread_fit["gic"]) - float(million["gic"])) < 1e-8
assert abs(float(thread_fit["MCC"]) - float(million["MCC"])) < 1e-12
(thread_directory / "README.md").write_text(
    "# OpenMP SIS/ISIS benchmark\n\n"
    "The design has 1,000 observations, one million iid standard-normal predictors (8 GB), "
    "60 signed effects in [-0.8, 0.8], problem seed 123 and noise seed 1. "
    "HP ZBook Fury G8, Intel i9-11950H, eight physical cores, 16 logical cores and 64 GB RAM; "
    "Windows and R 4.5.0. BLAS and fitting use one thread; screening concurrency varies.\n\n"
    f'SIS retains {2*int(1000/log(1000))} candidates. Complete ISIS uses the GIC-based '
    "batch, up to ten rounds, relaxed fitting and automatic local refinement. Simulation, "
    "explicit GC, equality checks and final reporting are outside timed calls.\n\n"
    + "\n".join(thread_table) + "\n\n"
    "Every tested score, sign, rank, coefficient and ISIS history agrees exactly across "
    "thread counts. The original-dimension sparse coefficients are retained. All fits give "
    f'GIC {float(thread_fit["gic"]):.2f}, F1 {float(thread_fit["F1_score"]):.3f} and '
    f'MCC {float(thread_fit["MCC"]):.3f}.\n\n'
    "Workers retain private top-k heaps ordered by score and original column index. "
    "Each predictor keeps its summation order, and merging the heaps preserves exact "
    "serial selection, including ties. No floating-point reduction spans workers.\n\n"
    "Reproduce with [3_benchmarks_sis_openmp_threads.R](../3_benchmarks_sis_openmp_threads.R). "
    "`timings.csv`, `summary.csv` and `isis_history.csv` retain the full results.\n",
    encoding="utf-8")
normal_path = ROOT / "benchmarks/plots/file-backed/normal-10000-1000000-metrics.csv"
lasso_path = ROOT / "benchmarks/plots/file-backed/normal-10000-1000000-lasso.csv"
gwas_path = ROOT / "benchmarks/plots/file-backed/gwas-4000-1000000-metrics.csv"
if normal_path.exists():
    normal = read_csv(normal_path)
    assert int(normal[0]["observations"]) == 10000 and int(normal[0]["p"]) == 1000000
    content["file-normal-metrics"] = table(normal,["seconds", "gic", "TP", "FP", "FN", "TN", "F1_score", "MCC", "test_RMSE"]) + (
        f'\n\nGenerating the {float(normal[0]["file_GB"]):.0f} GB file took '
        f'**{float(normal[0]["generation_seconds"]):.2f} seconds**, excluded from the fitting time.')
if gwas_path.exists():
    gwas = read_csv(gwas_path)
    content["file-gwas-metrics"] = table(gwas,["method", "seconds", "gic", "TP", "FP", "FN", "F1_score", "MCC", "interaction_recovered", "test_AUC", "test_logloss"])
if lasso_path.exists():
    lasso = read_csv(lasso_path)
    if len(lasso) == 2:
        assert all(int(r["observations"]) == 10000 and r["intercept"] == "TRUE" for r in lasso)
        assert next(r for r in lasso if r["method"] == "biglasso")["grid_interior"] == "TRUE"
        content["file-lasso-metrics"] = table(lasso,["method", "tune", "seconds", "gic", "TP", "FP", "FN", "F1_score", "MCC", "test_RMSE"])
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
