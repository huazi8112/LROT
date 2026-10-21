# -*- coding: utf-8 -*-
"""Matched-mask ERBB2/HER2 ablation for breast-cancer Visium sections.

All variants use the same HUMAN_ALL_LR exclusion mask for expression cost.
Only the LR guidance library changes. The tissue has no spot-level ground-truth
labels, so this analysis reports entropy and transport cost, not mapping accuracy.
"""
import json
import os
import sys
import warnings
warnings.filterwarnings("ignore")

_R = os.path.dirname(os.path.abspath(__file__))
while _R != os.path.dirname(_R) and not os.path.exists(os.path.join(_R, "lrot_core.py")):
    _R = os.path.dirname(_R)
sys.path.insert(0, _R)
sys.path.insert(0, os.path.join(_R, "experiments"))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from run_cancer_visium_truecoord import (
    load_breast_slice, apply_misalignment, compute_lr_strength_human,
)
from run_cancer_visium import HUMAN_LR_DB, HUMAN_ALL_LR
from lrot_core import fgw_lr_solver, compute_ot_cost

OUT_DIR = os.environ.get("LROT_OUTPUT_DIR", os.path.join(_R, "lrot_output"))
N = 3000
GAMMA = 0.1
SEED_A, SEED_B = 42, 99
ERBB2_PAIRS = {
    ("ERBB2", "ERBB2"), ("ERBB2", "ERBB3"), ("NRG1", "ERBB3"),
    ("NRG1", "ERBB4"), ("EGF", "EGFR"), ("EREG", "EGFR"),
    ("BTC", "EGFR"), ("HBEGF", "EGFR"), ("TGFA", "EGFR"),
}
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False


def mean_entropy(P):
    Pn = P / np.maximum(P.sum(1, keepdims=True), 1e-10)
    return float(-np.sum(Pn * np.log(np.maximum(Pn, 1e-10)), axis=1).mean())


def build_db(variant):
    if variant == "full":
        return HUMAN_LR_DB
    if variant == "no_erbb2":
        return {k: v for k, v in HUMAN_LR_DB.items() if k not in ERBB2_PAIRS}
    if variant == "erbb2_only":
        return {k: v for k, v in HUMAN_LR_DB.items() if k in ERBB2_PAIRS}
    raise ValueError(variant)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    print("=" * 72)
    print("Breast cancer ERBB2/HER2 ablation (matched expression mask)")
    print(f"Fixed mask: HUMAN_ALL_LR ({len(HUMAN_ALL_LR)} symbols); "
          "same mask in FGW and every LROT variant")
    A = load_breast_slice("V1_Breast_Cancer_Block_A_Section_1", N, SEED_A)
    B = load_breast_slice("V1_Breast_Cancer_Block_A_Section_2", N, SEED_B)
    B["coords"] = apply_misalignment(B["coords"], rotation_deg=2.0, translation=0.5)
    print(f"Sections: {A['n_spots']} × {A['n_genes']}; "
          f"Section 1 -> 2; seeds {SEED_A}/{SEED_B}; target misalignment 2°/0.5")

    zero = np.zeros((A["n_spots"], B["n_spots"]), dtype=np.float64)
    P0, loss0 = fgw_lr_solver(
        A, B, zero, gamma=0.0, beta=0.3, alpha=0.5,
        max_iter=200, verbose=False, lr_gene_list=HUMAN_ALL_LR)
    ent0 = mean_entropy(P0)
    cost0 = float(compute_ot_cost(P0, A["expr"], B["expr"]))
    print(f"Matched-mask FGW: entropy={ent0:.6f}; transport cost={cost0:.6f}; "
          f"iterations={len(loss0)}")

    variants = [
        ("full", "Full LR library"),
        ("no_erbb2", "Remove ERBB2/HER2"),
        ("erbb2_only", "ERBB2/HER2 only"),
    ]
    results = {}
    genes_A = {str(g).upper() for g in A["gene_names"]}
    genes_B = {str(g).upper() for g in B["gene_names"]}
    print("\nVariant                 pairs  entropy   reduction(%)  Δ vs full(pp)  cost")
    print("-" * 79)
    for variant, label in variants:
        db = build_db(variant)
        S, top = compute_lr_strength_human(A, B, db)
        pair_matches = sum(1 for lig, rec in db
                           if lig.upper() in genes_A and rec.upper() in genes_B)
        P, loss = fgw_lr_solver(
            A, B, S, gamma=GAMMA, beta=0.3, alpha=0.5,
            max_iter=200, verbose=False, lr_gene_list=HUMAN_ALL_LR)
        ent = mean_entropy(P)
        cost = float(compute_ot_cost(P, A["expr"], B["expr"]))
        reduction = (ent0 - ent) / ent0 * 100
        results[variant] = {
            "label": label, "lr_pairs_in_database": len(db),
            "lr_pairs_with_symbols_in_both_matrices": pair_matches,
            "entropy": ent, "reduction_percent_vs_matched_fgw": reduction,
            "cost": cost, "iterations": len(loss),
            "top_pairs": [[str(l), str(r), float(w)] for l, r, w in top[:5]],
        }
        print(f"{label:<23} {pair_matches:>3}/68  {ent:>8.4f}  {reduction:>12.4f}  "
              f"{reduction:>+13.4f}  {cost:>7.4f}")

    full_reduction = results["full"]["reduction_percent_vs_matched_fgw"]
    for key in ("no_erbb2", "erbb2_only"):
        results[key]["delta_reduction_percentage_points"] = (
            results[key]["reduction_percent_vs_matched_fgw"] - full_reduction)
    print("\nDescriptive contrasts relative to the full library:")
    for key in ("no_erbb2", "erbb2_only"):
        print(f"  {results[key]['label']}: "
              f"Δ={results[key]['delta_reduction_percentage_points']:+.4f} percentage points")
    print("No spot-level accuracy is reported: these breast sections lack true spot labels.")
    print("These single section-pair ablations are descriptive and do not establish causal "
          "or mechanistic pathway effects.")

    payload = {
        "analysis": "matched-mask ERBB2/HER2 ablation",
        "expression_mask": "HUMAN_ALL_LR (frozen for matched FGW and all LROT variants)",
        "expression_mask_symbols": len(HUMAN_ALL_LR),
        "sample_size_per_section": N,
        "seeds": {"A": SEED_A, "B": SEED_B},
        "gamma": GAMMA,
        "fgw_entropy": ent0,
        "fgw_transport_cost": cost0,
        "fgw_iterations": len(loss0),
        "delta_definition": "variant reduction minus full-library reduction, percentage points",
        "variants": results,
        "interpretation_limit": "one section pair; descriptive only; no spot-level ground-truth labels",
    }
    json_path = os.path.join(OUT_DIR, "lrot_cancer_erbb2_ablation_truecoord_results.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)
        f.write("\n")
    txt_path = os.path.join(OUT_DIR, "lrot_cancer_erbb2_ablation_truecoord_results.txt")
    with open(txt_path, "w", encoding="utf-8") as f:
        f.write("Matched-mask ERBB2/HER2 ablation; same HUMAN_ALL_LR expression mask in every arm.\n")
        f.write(f"FGW entropy={ent0:.6f}; cost={cost0:.6f}; iterations={len(loss0)}\n")
        for key, label in variants:
            r = results[key]
            delta = r.get("delta_reduction_percentage_points", 0.0)
            f.write(f"{label}: pairs={r['lr_pairs_with_symbols_in_both_matrices']}/68; "
                    f"entropy={r['entropy']:.6f}; reduction={r['reduction_percent_vs_matched_fgw']:.6f}%; "
                    f"delta_vs_full={delta:+.6f} pp; cost={r['cost']:.6f}; "
                    f"iterations={r['iterations']}\n")
            f.write("  top pairs: " + ", ".join(f"{x[0]}-{x[1]}" for x in r["top_pairs"]) + "\n")
        f.write("No spot-level accuracy: the breast sections lack true spot labels. "
                "Descriptive single-pair ablation; no causal or mechanistic inference.\n")

    reductions = [results[k]["reduction_percent_vs_matched_fgw"]
                  for k, _ in variants]
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    labels = [label for _, label in variants]
    colors = ["#2166AC", "#D6604D", "#4DAF4A"]
    bars = ax.bar(labels, reductions, color=colors, alpha=0.88, width=0.6)
    for bar, value in zip(bars, reductions):
        ax.text(bar.get_x() + bar.get_width()/2, value + max(reductions)*0.025,
                f"{value:.2f}%", ha="center", fontsize=10, fontweight="bold")
    ax.set_ylabel("Transport-entropy reduction vs matched-mask FGW (%)")
    ax.set_title("ERBB2/HER2 ablation with a fixed expression-feature mask")
    ax.grid(True, axis="y", alpha=0.3)
    ax.set_ylim(min(0.0, min(reductions)*1.1), max(reductions)*1.18)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "lrot_cancer_erbb2_ablation_truecoord.png"),
                dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Saved:", txt_path, json_path,
          os.path.join(OUT_DIR, "lrot_cancer_erbb2_ablation_truecoord.png"), flush=True)


if __name__ == "__main__":
    main()
