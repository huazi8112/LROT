# -*- coding: utf-8 -*-
"""Build Supplementary Figure S14 from machine-readable matched-mask outputs.

Both panels are regenerated from the same 3,000-spot section pair and the
fixed HUMAN_ALL_LR expression-feature mask. Panel B uses Δ = (reduction after
family removal) − (full-library reduction), in percentage points. The ±0.1 pp
band is a descriptive display threshold, not a statistical test.
"""
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
from PIL import Image, ImageDraw, ImageFont

_HERE = Path(__file__).resolve()
ROOT = next(p for p in _HERE.parents if (p / "lrot_core.py").is_file())
OUT_DIR = Path(os.environ.get("LROT_OUTPUT_DIR", ROOT / "lrot_output"))
FIG_DIR = Path(os.environ.get("LROT_FIGURE_OUTPUT_DIR", ROOT.parent / "outputs" / "generated_figures" / "S14"))
if FIG_DIR.resolve() == ROOT or ROOT in FIG_DIR.resolve().parents:
    raise ValueError("Use LROT_FIGURE_OUTPUT_DIR outside the code package")
OUT_FINAL = FIG_DIR / "S14.png"
ERBB_JSON = OUT_DIR / "lrot_cancer_erbb2_ablation_truecoord_results.json"
FAMILY_JSON = OUT_DIR / "lrot_cancer_ablation_scan_truecoord_results.json"

W, H = 4980, 3027
H_A, CUT_A, SEAM = 1483, 1483, 105
H_B = H - CUT_A - SEAM
FRAME = {"left": 373, "right": 4931, "top": 66, "bottom": 1303}
NAMES = [
    "ERBB2/HER", "Cell Adhesion", "Semaphorin", "Neurotrophin/Growth",
    "Ephrin", "Notch", "Slit/Robo", "Hedgehog", "Wnt", "BMP", "Reelin",
    "Chemokine/Cytokine",
]
DISPLAY_NAMES = [
    "ERBB2/HER", "Cell Adhesion\n", "Semaphorin", "Neurotrophin/\nGrowth",
    "Ephrin", "Notch", "Slit/Robo", "Hedgehog", "Wnt", "BMP", "Reelin",
    "Chemokine/\nCytokine",
]
BLUE, RED, GRAY = "#2166AC", "#D6604D", "#888888"
BAND_PP = 0.1

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Microsoft YaHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["savefig.facecolor"] = "white"


def load_results():
    with ERBB_JSON.open(encoding="utf-8") as f:
        erbb = json.load(f)
    with FAMILY_JSON.open(encoding="utf-8") as f:
        family = json.load(f)
    # The two analysis scripts use slightly different explanatory suffixes in
    # their JSON labels. Validate the actual frozen feature set and its size.
    masks = (erbb.get("expression_mask", ""), family.get("expression_mask", ""))
    if any("HUMAN_ALL_LR" not in mask for mask in masks):
        raise ValueError("Both S14 panels must use the HUMAN_ALL_LR expression mask.")
    if erbb.get("expression_mask_symbols") != family.get("expression_mask_symbols"):
        raise ValueError("The S14 panels used different expression-mask sizes.")
    if erbb.get("expression_mask_symbols") != 115:
        raise ValueError("Unexpected HUMAN_ALL_LR expression-mask size.")
    by_family = {x["family"]: x for x in family["family_ablations"]}
    missing = set(NAMES) - set(by_family)
    if missing:
        raise ValueError(f"Missing family results: {sorted(missing)}")
    return erbb, family, by_family


def render_panel(path, height, draw_panel, title):
    fig = plt.figure(figsize=(W / 300.0, height / 300.0), dpi=300)
    ax = fig.add_axes([
        FRAME["left"] / float(W),
        1.0 - FRAME["bottom"] / float(height),
        (FRAME["right"] - FRAME["left"]) / float(W),
        (FRAME["bottom"] - FRAME["top"]) / float(height),
    ])
    draw_panel(ax)
    ax.set_title(title, fontsize=11.5, pad=2)
    ax.grid(True, axis="y", alpha=0.3, zorder=0)
    ax.set_axisbelow(True)
    fig.savefig(path, dpi=300)
    plt.close(fig)
    size = Image.open(path).size
    if size != (W, height):
        raise RuntimeError(f"Panel size mismatch: {size}, expected {(W, height)}")


def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    work_dir = FIG_DIR / "panels"
    work_dir.mkdir(parents=True, exist_ok=True)
    erbb, family, by_family = load_results()
    full = erbb["variants"]["full"]["reduction_percent_vs_matched_fgw"]
    erbb_values = [
        erbb["variants"][k]["reduction_percent_vs_matched_fgw"]
        for k in ("full", "no_erbb2", "erbb2_only")
    ]
    family_deltas = [by_family[name]["delta_percentage_points"] for name in NAMES]

    panel_a = work_dir / "_S14_panelA_new.png"
    panel_b = work_dir / "_S14_panelB_new.png"

    def draw_a(ax):
        labels = ["Full LR library", "Remove ERBB2/HER2", "ERBB2/HER2 only"]
        colors = [BLUE, RED, "#4DAF4A"]
        bars = ax.bar(np.arange(3), erbb_values, color=colors, alpha=0.88,
                      width=0.55, zorder=3)
        top = max(erbb_values)
        ax.set_xticks(np.arange(3), labels, fontsize=11)
        ax.set_ylabel("Transport-entropy reduction vs matched-mask FGW (%)", fontsize=12)
        ax.set_ylim(0, max(0.1, top * 1.32))
        for bar, value in zip(bars, erbb_values):
            ax.text(bar.get_x() + bar.get_width()/2, value + top*0.025,
                    f"{value:.2f}%", ha="center", va="bottom",
                    fontsize=12, fontweight="bold", color="#222222", zorder=4)

    def draw_b(ax):
        x = np.arange(len(NAMES))
        colors = [
            BLUE if d < -BAND_PP else (RED if d > BAND_PP else GRAY)
            for d in family_deltas
        ]
        ax.axhspan(-BAND_PP, BAND_PP, color="#cccccc", alpha=0.4, zorder=0)
        bars = ax.bar(x, family_deltas, color=colors, width=0.62,
                      edgecolor="#333333", linewidth=0.6, zorder=3)
        ax.axhline(0, color="black", linewidth=1.0, zorder=4)
        span = max(BAND_PP * 1.5, max(abs(d) for d in family_deltas) * 1.25)
        span = np.ceil(span * 10) / 10
        ax.set_ylim(-span, span)
        ax.set_yticks(np.arange(-span, span + 0.05, 0.1))
        ax.set_xticks(x, DISPLAY_NAMES, fontsize=10.5)
        ax.set_xlim(-0.75, len(NAMES) - 0.25)
        ax.set_ylabel("Change in entropy reduction (percentage points)", fontsize=12)
        ax.text(0.99, 0.96, f"Full-library reduction: {full:.3f}%",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=10.5, color="#333333")
        for i, (bar, delta) in enumerate(zip(bars, family_deltas)):
            # Label every saved result; small values must not appear to be zero.
            offset = span * 0.035 * (1 if delta >= 0 else -1)
            label = f"{delta:+.2e}" if 0 < abs(delta) < 0.001 else f"{delta:+.3g}"
            ax.text(bar.get_x() + bar.get_width()/2, delta + offset,
                    label, ha="center",
                    va="bottom" if delta >= 0 else "top",
                    fontsize=11, fontweight="bold" if abs(delta) > BAND_PP else "normal",
                    color="#222222", zorder=5)
        ax.legend(handles=[
            Patch(facecolor=BLUE, edgecolor="#333333",
                label="Removal lowers reduction (Δ < −0.1 pp)"),
            Patch(facecolor=GRAY, edgecolor="#333333",
                label="Within descriptive band (|Δ| ≤ 0.1 pp)"),
            Patch(facecolor=RED, edgecolor="#333333",
                label="Removal raises reduction (Δ > +0.1 pp)"),
        ], fontsize=10, loc="upper center", ncol=3, frameon=False,
           bbox_to_anchor=(0.5, 0.98))
    render_panel(panel_a, H_A, draw_a,
                 "Breast-cancer ERBB2/HER2 ablation (fixed expression-feature mask)")
    render_panel(panel_b, H_B, draw_b,
                 "Breast-cancer family-level ablation (12 LR families; descriptive)")

    canvas = Image.new("RGB", (W, H), "white")
    canvas.paste(Image.open(panel_a).convert("RGB"), (0, 0))
    canvas.paste(Image.open(panel_b).convert("RGB"), (0, CUT_A + SEAM))
    try:
        font = ImageFont.truetype("arialbd.ttf", int(H * 0.032))
    except OSError:
        font = ImageFont.load_default()
    draw = ImageDraw.Draw(canvas)
    for letter, y in (("(A)", 42), ("(B)", CUT_A + SEAM + 12)):
        draw.text((60, y), letter, font=font, fill=(20, 30, 40))
    canvas.save(OUT_FINAL, dpi=(300, 300))
    if Image.open(OUT_FINAL).size != (W, H):
        raise RuntimeError("S14 composite size mismatch")
    print(f"S14 saved: {OUT_FINAL} ({W}×{H})")
    print(f"Matched-mask full reduction={full:.6f}%; family Δ pp={family_deltas}")


if __name__ == "__main__":
    main()
