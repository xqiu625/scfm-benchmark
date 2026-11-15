"""Shared style constants + helpers for the v2 revision figure set.

Fonts, colors, and mismatch-flagging come from
reports/manuscript/revision/figure_generation_spec.md's "全局风格" section.
"""
import os

import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
matplotlib.rcParams["pdf.fonttype"] = 42  # embed as real text, not paths
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["axes.spines.top"] = False
matplotlib.rcParams["axes.spines.right"] = False
matplotlib.rcParams["figure.dpi"] = 100

DISEASE_COLOR = {"COVID-19": "#C0392B", "Sepsis": "#E67E22", "Healthy": "#7F8C8D"}

MODEL_COLOR = {
    "STATE": "#2C3E50",
    "UCE": "#8E44AD",
    "Geneformer": "#2980B9",
    "scPRINT-2": "#27AE60",
    "scGPT": "#D35400",
    "PCA-200": "#7F8C8D",
    "PCA-200+LR": "#7F8C8D",
    "null": "#BDC3C7",
}

MODEL_ORDER = ["STATE", "UCE", "Geneformer", "scPRINT-2", "scGPT"]

OUT_DIR = "reports/manuscript/figures/v2"
DPI = 300

_MISMATCHES = []


def flag_mismatch(figure, quantity, spec_value, csv_value, csv_path, note=""):
    """Record a spec-anchor vs. CSV-recomputed mismatch for the final report."""
    msg = (f"[{figure}] {quantity}: spec/P0 anchor = {spec_value}, "
           f"recomputed from {csv_path} = {csv_value}. {note}")
    _MISMATCHES.append(msg)
    print("MISMATCH: " + msg, flush=True)
    return msg


def get_mismatches():
    return list(_MISMATCHES)


def check_close(figure, quantity, spec_value, csv_value, csv_path, tol=0.0005, note=""):
    """Assert-and-record: verify csv_value is within tol of spec_value; flag if not."""
    if abs(spec_value - csv_value) > tol:
        flag_mismatch(figure, quantity, spec_value, csv_value, csv_path, note)
        return False
    return True


def save_fig(fig, name):
    os.makedirs(OUT_DIR, exist_ok=True)
    pdf_path = os.path.join(OUT_DIR, f"{name}.pdf")
    png_path = os.path.join(OUT_DIR, f"{name}.png")
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"SAVED {pdf_path}, {png_path}", flush=True)
    return pdf_path, png_path


def sig_stars(q):
    if q < 0.001:
        return "***"
    if q < 0.01:
        return "**"
    if q < 0.05:
        return "*"
    return "n.s."
