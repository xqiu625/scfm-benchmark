"""Programmatic API for scfm-benchmark.

Functions here are import-friendly so notebooks can call them directly,
e.g. `from scfm_benchmark.api import build_report`.
"""

from __future__ import annotations

from pathlib import Path

SECTIONS = [
    ("Primary benchmark", "primary_benchmark*.csv"),
    ("Calibration per cohort", "calibration*.csv"),
    ("Adaptation curve", "adaptation_curve*.csv"),
    ("GSEA per FM dimension", "gsea*.csv"),
    ("SHAP attribution", "*shap*.csv"),
    ("Cross-disease transfer", "cross_disease*.csv"),
]


def build_report(results_dir: Path, out_path: Path) -> int:
    """Render results_dir/*.{csv,json,png} into a self-contained HTML report.

    Looks for canonical artifact filenames produced by the Option C scripts.
    Each section is optional; missing files are silently skipped.
    """
    if not results_dir.exists():
        print(f"results_dir does not exist: {results_dir}")
        return 1

    out_path.parent.mkdir(parents=True, exist_ok=True)

    body: list[str] = [
        "<h1>scFM Benchmark Report</h1>",
        f"<p><em>Source: <code>{results_dir}</code></em></p>",
    ]
    for label, glob in SECTIONS:
        files = sorted(results_dir.rglob(glob))
        if not files:
            continue
        body.append(f"<h2>{label}</h2><ul>")
        for f in files:
            body.append(f"<li><code>{f.relative_to(results_dir)}</code></li>")
        body.append("</ul>")

    figures = sorted(results_dir.rglob("*.png"))
    if figures:
        body.append("<h2>Figures</h2>")
        for fig in figures:
            rel = fig.relative_to(results_dir)
            body.append(f"<figure><img src='{rel}' alt='{rel}'>")
            body.append(f"<figcaption><code>{rel}</code></figcaption></figure>")

    html = (
        "<!doctype html><html><head><meta charset='utf-8'>"
        "<title>scFM Benchmark Report</title>"
        "<style>"
        "body{font-family:system-ui,sans-serif;max-width:960px;"
        "margin:2em auto;padding:1em;line-height:1.5}"
        "code{background:#f4f4f4;padding:1px 4px;border-radius:3px}"
        "figure{margin:1em 0}img{max-width:100%}"
        "</style></head><body>"
        + "".join(body)
        + "</body></html>"
    )
    out_path.write_text(html)
    print(f"Wrote {out_path}")
    return 0
