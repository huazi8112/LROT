"""Reproduce main-text Figure 8: descriptive DLPFC score-ranking analyses.

The annotation-derived score uses target-section layer labels and is a
retrospective diagnostic. It should not be interpreted as a label-free
confidence estimate or as probability calibration.
"""
from __future__ import annotations

import argparse
import gc
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import h5py
import numpy as np
from scipy.sparse import csc_matrix
from scipy.stats import spearmanr


DEFAULT_ROOT = Path(__file__).resolve().parents[1]
SECTION_SEEDS = {"151507": 42, "151508": 99, "151509": 99}
SECTION_PAIRS = (("151507", "151508"), ("151508", "151509"))
LAYER_ORDER = ("L1", "L2", "L3", "L4", "L5", "L6", "WM")
N_SPOTS = 1000
ROTATION_DEG = 2.0
TRANSLATION = 0.5
ALPHA = 0.5
BETA = 0.3
GAMMA_LROT = 0.1
MAX_ITER = 200


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=DEFAULT_ROOT,
                        help="LROT repository root (defaults to the parent of figures/).")
    parser.add_argument("--data-dir", type=Path,
                        help="DLPFC raw data directory (defaults to <repo-root>/dlpfc_data).")
    parser.add_argument("--output-dir", type=Path,
                        help="Output directory (defaults to <repo-root>/lrot_output/figures_png).")
    return parser.parse_args()


def unique_gene_names(names):
    """Match the manuscript loader's duplicate-gene-name handling."""
    seen = set()
    next_suffix = {}
    unique = []
    for raw in names:
        name = str(raw)
        if name not in seen:
            seen.add(name)
            unique.append(name)
            continue
        suffix = next_suffix.get(name, 1)
        candidate = f"{name}-{suffix}"
        while candidate in seen:
            suffix += 1
            candidate = f"{name}-{suffix}"
        next_suffix[name] = suffix + 1
        seen.add(candidate)
        unique.append(candidate)
    return unique


def read_10x_h5(sample_id: str, allowed_barcodes: set[str], data_dir: Path):
    """Read selected gene-expression barcodes from a standard 10x HDF5 file."""
    path = data_dir / f"{sample_id}_filtered_feature_bc_matrix.h5"
    with h5py.File(path, "r") as handle:
        matrix = handle["matrix"]
        decode = lambda value: value.decode("utf-8") if isinstance(value, (bytes, np.bytes_)) else str(value)
        all_barcodes = [decode(value) for value in matrix["barcodes"][:]]
        features = matrix["features"]
        names_all = [decode(value) for value in features["name"][:]]
        feature_types = [decode(value) for value in features["feature_type"][:]]
        gene_indices = np.array(
            [index for index, feature_type in enumerate(feature_types)
             if feature_type == "Gene Expression"], dtype=int
        )
        gene_names = unique_gene_names([names_all[index] for index in gene_indices])
        shape = tuple(int(value) for value in matrix["shape"][:])
        expression = csc_matrix(
            (matrix["data"][:], matrix["indices"][:], matrix["indptr"][:]),
            shape=shape,
        )

    barcode_to_index = {barcode: index for index, barcode in enumerate(all_barcodes)}
    selected = [barcode for barcode in all_barcodes if barcode in allowed_barcodes]
    cell_indices = np.array([barcode_to_index[barcode] for barcode in selected], dtype=int)
    selected_expression = expression[gene_indices, :][:, cell_indices].T.toarray()
    return selected, gene_names, selected_expression.astype(np.float32, copy=False)


def load_sections(data_dir: Path):
    required = [
        data_dir / f"{sid}_filtered_feature_bc_matrix.h5" for sid in SECTION_SEEDS
    ] + [
        data_dir / f"{sid}_tissue_positions_list.txt" for sid in SECTION_SEEDS
    ] + [data_dir / "barcode_level_layer_map.tsv"]
    missing = [path for path in required if not path.is_file()]
    if missing:
        names = "\n".join(f"  - {path.name}" for path in missing)
        raise FileNotFoundError(
            "Missing DLPFC inputs:\n" + names +
            "\nRun `python experiments/supp/download_dlpfc_data.py` from the repository root."
        )

    label_map = {label: index for index, label in enumerate(LAYER_ORDER)}
    sections = {}
    for sid, seed in SECTION_SEEDS.items():
        positions = {}
        with (data_dir / f"{sid}_tissue_positions_list.txt").open(encoding="utf-8") as handle:
            for line in handle:
                fields = line.rstrip("\n").split(",")
                if len(fields) == 6:
                    positions[fields[0]] = (float(fields[4]), float(fields[5]))

        layers = {}
        with (data_dir / "barcode_level_layer_map.tsv").open(encoding="utf-8") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) == 3 and fields[1] == sid:
                    layers[fields[0]] = fields[2]

        eligible = set(positions).intersection(layers)
        barcodes, genes, expression = read_10x_h5(sid, eligible, data_dir)
        labels = np.array([label_map[layers[barcode]] for barcode in barcodes], dtype=int)
        coordinates = np.array([positions[barcode] for barcode in barcodes], dtype=np.float64)
        coordinates -= coordinates.mean(axis=0)
        rng = np.random.RandomState(seed)
        indices = rng.choice(len(barcodes), min(N_SPOTS, len(barcodes)), replace=False)
        sections[sid] = {
            "coords": coordinates[indices], "expr": expression[indices],
            "gene_names": genes, "region_labels": labels[indices],
            "n_spots": len(indices), "n_genes": len(genes), "sample_id": sid,
        }
    expected = set(range(len(LAYER_ORDER)))
    for sid, section in sections.items():
        labels = set(np.unique(section["region_labels"]).tolist())
        if not labels.issubset(expected):
            raise ValueError(f"Unexpected layer labels in section {sid}: {sorted(labels)}")
    return sections


def compute_ranking_data(repo_root: Path, sections):
    sys.path.insert(0, str(repo_root))
    sys.path.insert(0, str(repo_root / "experiments"))
    from lrot_core import fgw_lr_solver
    from run_cancer_visium import (
        HUMAN_ALL_LR,
        HUMAN_LR_DB,
        apply_misalignment,
        compute_lr_strength_human,
    )

    pooled = {
        name: {key: [] for key in ("score", "correct", "row_entropy")}
        for name in ("FGW", "LROT")
    }
    metrics = []

    for sid_a, sid_b in SECTION_PAIRS:
        source = sections[sid_a]
        target = dict(sections[sid_b])
        target["coords"] = apply_misalignment(
            target["coords"], rotation_deg=ROTATION_DEG, translation=TRANSLATION
        )
        strengths, _ = compute_lr_strength_human(source, target, HUMAN_LR_DB)

        for name, gamma, lr_strength in (
            ("FGW", 0.0, np.zeros_like(strengths)),
            ("LROT", GAMMA_LROT, strengths),
        ):
            plan, _ = fgw_lr_solver(
                source, target, lr_strength, gamma=gamma, beta=BETA, alpha=ALPHA,
                max_iter=MAX_ITER, verbose=False, lr_gene_list=HUMAN_ALL_LR,
            )
            row_plan = plan / np.maximum(plan.sum(axis=1, keepdims=True), 1e-12)
            target_one_hot = np.eye(len(LAYER_ORDER))[target["region_labels"]]
            layer_votes = row_plan @ target_one_hot
            score = layer_votes.max(axis=1)
            correct = (layer_votes.argmax(axis=1) == source["region_labels"]).astype(float)
            row_entropy = (-row_plan * np.log(np.maximum(row_plan, 1e-12))).sum(axis=1)

            pooled[name]["score"].append(score)
            pooled[name]["correct"].append(correct)
            pooled[name]["row_entropy"].append(row_entropy)
            top20 = np.argsort(score)[-int(round(len(score) * 0.2)):]
            metrics.append(
                f"PAIR\t{sid_a}->{sid_b}\t{name}\tn={len(score)}\t"
                f"spearman_rho={spearmanr(score, correct).statistic:.3f}\t"
                f"accuracy={correct.mean():.3f}\t"
                f"top20_accuracy={correct[top20].mean():.3f}"
            )
            del plan, row_plan, target_one_hot, layer_votes, score, correct, row_entropy
            gc.collect()
        del target, strengths
        gc.collect()

    for method in pooled:
        for metric in pooled[method]:
            pooled[method][metric] = np.concatenate(pooled[method][metric])
    return pooled, metrics


def make_figure(pooled, metrics, output_dir: Path):
    colors = {"LROT": "#D9823C", "FGW": "#8FA9C9"}
    fill_colors = {"LROT": "#DE9459", "FGW": "#A0B6D1"}
    plt.rcParams.update({
        "font.size": 11, "axes.titlesize": 13, "axes.labelsize": 12,
        "legend.fontsize": 10, "xtick.labelsize": 10, "ytick.labelsize": 10,
        "font.family": "DejaVu Sans", "axes.grid": True,
        "grid.color": "#E7E7E7", "grid.linewidth": 0.8,
    })
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 8.25), constrained_layout=True)

    # A: observed correctness by score decile; score is annotation-derived.
    ax = axes[0, 0]
    for name in ("FGW", "LROT"):
        score, correct = pooled[name]["score"], pooled[name]["correct"]
        bins = np.array_split(np.argsort(score), 10)
        x = np.array([score[indices].mean() for indices in bins])
        y = np.array([correct[indices].mean() for indices in bins])
        rho = spearmanr(score, correct).statistic
        ax.plot(x, y, marker="o", linewidth=2.0, markersize=6,
                color=colors[name], label=f"{name} (ρ={rho:.3f})")
    ax.set(title="(A) Annotation-derived score and observed correctness",
           xlabel="Mean score within decile", ylabel="Observed layer accuracy",
           ylim=(0.12, 1.0))
    ax.legend(frameon=False, loc="upper left")
    ax.grid(alpha=0.8)

    # B: cumulative retained accuracy after ranking spots by score.
    ax = axes[0, 1]
    for name in ("FGW", "LROT"):
        score, correct = pooled[name]["score"], pooled[name]["correct"]
        order = np.argsort(score)[::-1]
        fractions = np.linspace(0.05, 1.0, 20)
        accuracies = [
            correct[order[:max(1, int(round(len(order) * fraction))) ]].mean()
            for fraction in fractions
        ]
        gamma = GAMMA_LROT if name == "LROT" else 0
        ax.plot(fractions, accuracies, marker="o", markersize=5, linewidth=2.0,
                color=colors[name], label=f"{name} (γ={gamma:g})")
        top20 = order[:int(round(len(order) * 0.2))]
        metrics.append(
            f"POOLED\t{name}\tn={len(score)}\t"
            f"spearman_rho={spearmanr(score, correct).statistic:.3f}\t"
            f"accuracy={correct.mean():.3f}\t"
            f"top20_accuracy={correct[top20].mean():.3f}"
        )
    for name in ("LROT", "FGW"):
        ax.axhline(pooled[name]["correct"].mean(), color=colors[name],
                   linestyle="--", linewidth=1.2, alpha=0.8,
                   label=f"{name} full-set accuracy")
    ax.axhline(1 / len(LAYER_ORDER), color="#777777", linestyle=":",
               linewidth=1.2, label="Uniform 7-class guess (1/7)")
    ax.set(title="(B) Retention-accuracy ranking curve",
           xlabel="Retained fraction (descending score)",
           ylabel="Observed layer accuracy", xlim=(0.03, 1.02), ylim=(0.1, 1.0))
    ax.legend(frameon=True, framealpha=0.95, loc="upper right", fontsize=8.5)
    ax.grid(alpha=0.8)

    # C: coverage at score thresholds, with retained accuracy in parentheses.
    ax = axes[1, 0]
    thresholds = np.array([0.4, 0.5, 0.6, 0.7])
    x = np.arange(len(thresholds))
    width = 0.36
    for offset, name in ((-width / 2, "LROT"), (width / 2, "FGW")):
        score, correct = pooled[name]["score"], pooled[name]["correct"]
        coverage, accuracy = [], []
        for threshold in thresholds:
            keep = score >= threshold
            coverage.append(keep.mean() * 100)
            accuracy.append(correct[keep].mean() if keep.any() else np.nan)
        bars = ax.bar(x + offset, coverage, width, color=fill_colors[name],
                      edgecolor="white", linewidth=0.7,
                      label=f"{name} (retained accuracy)")
        for bar, pct, acc, threshold in zip(bars, coverage, accuracy, thresholds):
            ax.annotate(f"{pct:.1f}%\n({acc:.2f})",
                        (bar.get_x() + bar.get_width() / 2, bar.get_height()),
                        xytext=(0, 4), textcoords="offset points", ha="center",
                        va="bottom", fontsize=9)
            metrics.append(
                f"THRESHOLD\t{name}\tthreshold={threshold:.1f}\t"
                f"coverage={pct:.1f}%\tretained_accuracy={acc:.3f}"
            )
    ax.set_xticks(x, [f"score≥{threshold:.1f}" for threshold in thresholds])
    ax.set(title="(C) High-score coverage at fixed thresholds",
           xlabel="Annotation-derived score threshold",
           ylabel="Retained spots (%)", ylim=(0, 82))
    ax.legend(frameon=True, framealpha=0.95, loc="upper right")
    ax.grid(axis="y", alpha=0.8)

    # D: row entropy is computed from the transport plan and does not use labels.
    ax = axes[1, 1]
    for name in ("FGW", "LROT"):
        values = np.sort(pooled[name]["row_entropy"])
        ecdf = np.arange(1, len(values) + 1) / len(values)
        gamma = GAMMA_LROT if name == "LROT" else 0
        ax.plot(values, ecdf, linewidth=2.0, color=colors[name],
                label=f"{name} (γ={gamma:g})")
    ax.set(title="(D) Row entropy and overall plan ambiguity",
           xlabel="Transport row entropy", ylabel="Cumulative fraction",
           ylim=(0, 1.02))
    ax.legend(frameon=False, loc="lower right")
    ax.grid(alpha=0.8)

    output_dir.mkdir(parents=True, exist_ok=True)
    figure_path = output_dir / "lrot_figure8_dlpfc_score_ranking.png"
    fig.savefig(figure_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    metadata = [
        "Figure 8: descriptive DLPFC score-ranking analysis.",
        "Dataset: Human DLPFC Visium, donor Br8325; source data: spatialLIBD / HumanPilot.",
        "Sections: 151507(seed=42), 151508(seed=99), 151509(seed=99); up to 1,000 labelled spots per section.",
        "Pairs: 151507->151508 and 151508->151509; the pairs share section 151508 and are not independent replicates.",
        f"Target-coordinate perturbation: rotation={ROTATION_DEG:g} degrees, translation={TRANSLATION:g} pixels.",
        f"Solver: alpha={ALPHA:g}, beta={BETA:g}, max_iter={MAX_ITER}; FGW gamma=0, LROT gamma={GAMMA_LROT:g}.",
        "Score: maximum layer vote after row-normalizing P and aggregating by target-section layer labels.",
        "Panels A-C use target labels to evaluate a retrospective annotation-derived ranking score; they are not label-free confidence or calibration analyses.",
        "Panel D reports row entropy from P and is label-free.",
        "Pooled summaries are descriptive; no spot-level inferential test or confidence interval is reported.",
        "",
        *metrics,
    ]
    metrics_path = output_dir.parent / "figure8_score_ranking_results.txt"
    metrics_path.write_text("\n".join(metadata) + "\n", encoding="utf-8", newline="\n")
    return figure_path, metrics_path


def main():
    args = parse_args()
    repo_root = args.repo_root.resolve()
    data_dir = (args.data_dir or repo_root / "dlpfc_data").resolve()
    output_dir = (args.output_dir or repo_root / "lrot_output" / "figures_png").resolve()
    sections = load_sections(data_dir)
    pooled, metrics = compute_ranking_data(repo_root, sections)
    figure_path, metrics_path = make_figure(pooled, metrics, output_dir)
    print(f"Figure saved: {figure_path}")
    print(f"Metrics saved: {metrics_path}")


if __name__ == "__main__":
    main()
