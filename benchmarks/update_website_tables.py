"""Publish measured website tables after benchmark_website_examples.R finishes."""
import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "benchmarks/plots/website"
LABELS = {"L0glm": "L0glm", "lasso_gic": "LASSO", "relaxed_lasso_gic": "Relaxed LASSO",
          "lasso_cv": "LASSO", "relaxed_lasso_cv": "Relaxed LASSO"}

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

def table(data, fields):
    names = {"method": "Method", "tune": "Selection", "seconds": "Seconds", "gic": "GIC",
             "F1_score": "F1", "within_row_F1_score": "Within-row F1",
             "within_row_MCC": "Within-row MCC"}
    lines = ["| " + " | ".join(names.get(f, f) for f in fields) + " |",
             "| " + " | ".join("---" if f in {"method", "tune"} else "---:" for f in fields) + " |"]
    lines += ["| " + " | ".join(value(row, f) for f in fields) + " |" for row in data]
    return "\n".join(lines)

million = rows("million-metrics.csv")[0]
single = rows("single-metrics.csv")
multi = rows("multi-metrics.csv")
selection = rows("selection-comparison.csv")
fields = ["method", "tune", "seconds", "gic", "TP", "FP", "FN", "TN", "F1_score", "MCC"]
ratio = 123.13 / float(million["seconds"])
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
            f'with just {int(fitted["FP"]):,} false positives** in **{elapsed:.2f} seconds**: '
            f'**{plain_ratio:.0f}× faster than GIC-tuned {name}** and '
            f'**{relaxed_ratio:.0f}× faster than relaxed {name} tuned by five-fold CV** in this comparison.')


content = {
    "million-metrics": table([million], fields[2:]),
    "million-confusion": "| Ground truth | Selected | Excluded |\n|---|---:|---:|\n"
        f"| Nonzero effects | {million['TP']} TP | {million['FN']} FN |\n"
        f"| Zero effects | {million['FP']} FP | {million['TN']} TN |",
    "single-metrics": table(single, fields),
    "multi-metrics": table(multi_display, fields + ["within_row_F1_score", "within_row_MCC"]),
    "multi-home": table(multi_display, ["method","tune","seconds","TP","FP","F1_score","MCC"]),
    "single-speed": deconvolution_speed(single),
    "multi-speed": deconvolution_speed(multi, grouped=True),
    "selection-metrics": table(selection, ["method", "gic", "TP", "FP", "FN", "TN", "F1_score", "MCC"]),
    "speed-stat": f'<a href="articles/performance.html#one-million-predictor-example"><strong>{ratio:.0f}× faster</strong><span>than five-fold-CV LASSO in that example</span></a>',
    "multi-stat": f'<a href="articles/performance.html#multichannel-deconvolution"><strong>{multi[0]["TP"]} / {multi[0]["support_size_true"]}</strong><span>true peaks recovered in multichannel GC/MS deconvolution</span></a>',
}
paths = [ROOT/"L0glm/pkgdown/index.md", ROOT/"L0glm/vignettes/L0glm.Rmd"]
paths += list((ROOT/"L0glm/vignettes/articles").glob("*.Rmd"))
for path in paths:
    text = original = path.read_text(encoding="utf-8")
    for key, replacement in content.items():
        pattern = re.compile(r"(<!-- BEGIN " + re.escape(key) + r" -->).*?(<!-- END " + re.escape(key) + r" -->)", re.S)
        text = pattern.sub(lambda m: m[1] + "\n" + replacement + "\n" + m[2], text)
    if path.name == "performance.Rmd":
        line = (f'| **L0glm + ISIS** | GIC preset | **{float(million["seconds"]):.2f}** | **{float(million["gic"]):.0f}** | '
                f'**{million["TP"]}** | **{million["FP"]}** | {million["FN"]} | {million["TN"]} | '
                f'**{float(million["F1_score"]):.3f}** | **{float(million["MCC"]):.3f}** |')
        text = re.sub(r"^\| \*\*L0glm \+ ISIS.*$", lambda _: line, text, flags=re.M)
        text = re.sub(r"seconds <- c\([^\n]+", f'seconds <- c({float(million["seconds"]):.2f}, 123.13, 91.76)', text)
        text = re.sub(r"\*\*\d+ times faster than five-fold-CV LASSO\*\*", f'**{ratio:.0f} times faster than five-fold-CV LASSO**', text)
    if text != original:
        path.write_text(text, encoding="utf-8")
print("Updated measured tables and homepage statistics.")
