from __future__ import annotations

import argparse
import ast
import csv
import json
import sys
import time
from datetime import date
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy import stats


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = ROOT / "lrot_output" / "s23_support"
MANIFEST = ROOT / "experiments" / "synthetic_design_manifest.json"
SOURCE = ROOT / "experiments" / "lrot_prototype.py"
ARCHIVE = ROOT / "lrot_output" / "arrays" / "lrot_data_synthetic.npz"
N_REPLICATES = 30
SEED_PAIRS = [(1000 + i, 2000 + i) for i in range(N_REPLICATES)]
PARAMETERS = {
    "alpha": 0.5,
    "beta": 0.3,
    "gamma_fgw": 0.0,
    "gamma_lrot": 0.1,
    "reg_gw": 0.01,
    "reg_ot": 0.01,
    "max_iter": 120,
    "tol": 1e-5,
    "spots_per_slice": 1000,
    "genes": 200,
    "regions": 4,
}
OLD_TO_CURRENT = {"Artemin": "Artn", "Neurturin": "Nrtn"}

sys.path.insert(0, str(ROOT))
core = None  # Selected after argument parsing; solver identities must not mix.


def get_lr_db() -> dict:
    tree = ast.parse(SOURCE.read_text(encoding="utf-8"))
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
            isinstance(target, ast.Name) and target.id == "LR_DB" for target in node.targets
        ):
            return ast.literal_eval(node.value)
    raise RuntimeError("Could not read LR_DB from experiments/lrot_prototype.py")


def configure_core(db: dict) -> None:
    core.LR_DB = db
    core.LIGANDS = sorted({ligand for ligand, _ in db})
    core.RECEPTORS = sorted({receptor for _, receptor in db})
    core.ALL_LR_GENES = sorted(set(core.LIGANDS + core.RECEPTORS))
    core.N_REGIONS = PARAMETERS["regions"]


def simulate(seed: int, feature_order: list[str], groups: dict, weights: dict) -> dict:
    rng = np.random.RandomState(seed)
    n_spots = PARAMETERS["spots_per_slice"]
    n_genes = PARAMETERS["genes"]
    side = int(np.ceil(np.sqrt(n_spots)))
    xs = np.linspace(0, 20, side)
    ys = np.linspace(0, 20, side)
    xx, yy = np.meshgrid(xs, ys)
    coords = np.column_stack([xx.ravel(), yy.ravel()])
    coords = coords[:n_spots] + rng.randn(n_spots, 2) * 0.3
    center = np.array([10.0, 10.0])
    angles = np.arctan2(coords[:, 1] - center[1], coords[:, 0] - center[0])
    labels = np.floor((angles + np.pi) / (2 * np.pi / PARAMETERS["regions"])).astype(int)
    labels = np.clip(labels, 0, PARAMETERS["regions"] - 1)

    gene_names = (list(feature_order) + [f"Gene_{i}" for i in range(max(0, n_genes - len(feature_order)))])[:n_genes]
    expr = rng.gamma(1, 1, (n_spots, n_genes)).astype(np.float32) * 0.5
    gene_index = {gene: i for i, gene in enumerate(gene_names)}
    region_weights = {}
    for region in range(PARAMETERS["regions"]):
        selected = {}
        for gene in groups[str(region)]["ligands"]:
            if gene in gene_index:
                selected[gene_index[gene]] = float(weights["ligand_enrichment"])
        for gene in groups[str(region)]["receptors"]:
            if gene in gene_index:
                selected[gene_index[gene]] = float(weights["receptor_enrichment"])
        region_weights[region] = selected
    for spot, region in enumerate(labels):
        for gene_i, multiplier in region_weights[int(region)].items():
            expr[spot, gene_i] *= multiplier
    expr += rng.gamma(0.5, 0.5, expr.shape) * 0.3
    return {
        "coords": coords,
        "expr": expr,
        "region_labels": labels,
        "gene_names": gene_names,
        "n_spots": n_spots,
        "n_genes": n_genes,
    }


def make_pair(feature_order: list[str], groups: dict, weights: dict, seed_a: int, seed_b: int):
    A = simulate(seed_a, feature_order, groups, weights)
    B = simulate(seed_b, feature_order, groups, weights)
    theta = 0.15
    rotation = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    center = B["coords"].mean(0)
    B["coords"] = (B["coords"] - center) @ rotation.T + center + np.array([0.5, 0.3])
    return A, B


def verify_archive(db: dict, feature_order: list[str], groups: dict, weights: dict) -> dict:
    if not ARCHIVE.is_file():
        raise FileNotFoundError(f"Historical benchmark archive is missing: {ARCHIVE}")
    archived = np.load(ARCHIVE, allow_pickle=True)
    A, B = make_pair(feature_order, groups, weights, 42, 123)
    checks = {}
    for label, generated in (("A", A), ("B", B)):
        archived_names = [OLD_TO_CURRENT.get(str(name), str(name)) for name in archived[f"{label}_gene_names"]]
        checks[label] = {
            "coords_exact": bool(np.array_equal(generated["coords"], archived[f"{label}_coords"])),
            "labels_exact": bool(np.array_equal(generated["region_labels"], archived[f"{label}_region_labels"])),
            "gene_names_exact_after_alias_canonicalization": generated["gene_names"] == archived_names,
            "expr_equal_after_archival_float16": bool(
                np.array_equal(generated["expr"].astype(np.float16), archived[f"{label}_expr"])
            ),
        }
    if not all(all(value for value in item.values()) for item in checks.values()):
        raise RuntimeError("The frozen generator does not reproduce the archived 42/123 input; stopping.")
    return checks


def load_design():
    db = get_lr_db()
    design = json.loads(MANIFEST.read_text(encoding="utf-8"))
    configure_core(db)
    feature_order = design["feature_order"]
    groups = design["region_enrichment"]
    archive_check = verify_archive(db, feature_order, groups, design["weights"])
    return db, feature_order, groups, design["weights"], archive_check


def solve(A: dict, B: dict, db: dict, gamma: float) -> dict:
    S_lr, _ = core.compute_lr_strength_matrix(A, B, db, scale=10.0)
    started = time.time()
    P, _ = core.fgw_lr_solver(
        A,
        B,
        S_lr,
        gamma=gamma,
        beta=PARAMETERS["beta"],
        alpha=PARAMETERS["alpha"],
        reg_gw=PARAMETERS["reg_gw"],
        reg_ot=PARAMETERS["reg_ot"],
        max_iter=PARAMETERS["max_iter"],
        tol=PARAMETERS["tol"],
        verbose=False,
    )
    accuracy = float(core.compute_mapping_consistency(P, A["region_labels"], B["region_labels"]))
    return {"accuracy": accuracy, "elapsed_seconds": round(time.time() - started, 2)}


def read_jsonl(path: Path) -> dict[tuple[int, int], dict]:
    completed = {}
    if not path.exists():
        return completed
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        item = json.loads(line)
        pair = tuple(item["seed_pair"])
        if pair in completed:
            raise RuntimeError(f"Duplicate seed pair in checkpoint: {pair}")
        completed[pair] = item
    return completed


def summarize(rows: list[dict]) -> dict:
    fgw = np.asarray([row["fgw_accuracy"] for row in rows], dtype=float)
    lrot = np.asarray([row["lrot_accuracy"] for row in rows], dtype=float)
    delta = np.asarray([row["delta_lrot_minus_fgw"] for row in rows], dtype=float)
    n = len(delta)
    mean_delta = float(np.mean(delta))
    sd_delta = float(np.std(delta, ddof=1))
    half_width = float(stats.t.ppf(0.975, n - 1) * sd_delta / np.sqrt(n))
    draws = np.random.default_rng(20260930).choice(delta, size=(100000, n), replace=True).mean(axis=1)
    return {
        "n": n,
        "fgw_mean": float(np.mean(fgw)),
        "fgw_sd_sample": float(np.std(fgw, ddof=1)),
        "lrot_mean": float(np.mean(lrot)),
        "lrot_sd_sample": float(np.std(lrot, ddof=1)),
        "delta_mean": mean_delta,
        "delta_sd_sample": sd_delta,
        "delta_median": float(np.median(delta)),
        "delta_min": float(np.min(delta)),
        "delta_max": float(np.max(delta)),
        "delta_ci95_t": [mean_delta - half_width, mean_delta + half_width],
        "delta_ci95_bootstrap_percentile_sensitivity": [float(v) for v in np.quantile(draws, [0.025, 0.975])],
        "bootstrap_sensitivity_seed": 20260930,
        "positive_delta_n": int(np.sum(delta > 0)),
        "zero_delta_n": int(np.sum(delta == 0)),
        "negative_delta_n": int(np.sum(delta < 0)),
        "all_primary_runs_are_new_prespecified_seed_pairs": True,
        "interpretation_scope": "frozen historical synthetic grouping for technical robustness only; not biological pathway validation",
        "software": {
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
        },
    }


def write_derived(rows: list[dict], archive_check: dict, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    delta_rows = np.asarray([row["delta_lrot_minus_fgw"] for row in rows], dtype=float)
    summary = summarize(rows)
    (out_dir / "historical_archive_reconstruction_check.json").write_text(
        json.dumps(archive_check, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    (out_dir / "synthetic_seed_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    fields = ["replicate", "seed_A", "seed_B", "fgw_accuracy", "lrot_accuracy", "delta_lrot_minus_fgw", "fgw_elapsed_seconds", "lrot_elapsed_seconds"]
    with (out_dir / "Supplementary_File_S2_synthetic_seed_pairs.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({
                "replicate": row["replicate"],
                "seed_A": row["seed_pair"][0],
                "seed_B": row["seed_pair"][1],
                "fgw_accuracy": row["fgw_accuracy"],
                "lrot_accuracy": row["lrot_accuracy"],
                "delta_lrot_minus_fgw": row["delta_lrot_minus_fgw"],
                "fgw_elapsed_seconds": row["fgw_elapsed_seconds"],
                "lrot_elapsed_seconds": row["lrot_elapsed_seconds"],
            })
    table = [
        "| Outcome | FGW (γ=0) | LROT (γ=0.1) | Paired difference (LROT−FGW) |",
        "|---|---:|---:|---:|",
        f"| Mean (n={summary['n']}) | {summary['fgw_mean']:.3f} | {summary['lrot_mean']:.3f} | {summary['delta_mean']:+.3f} |",
        f"| Sample SD | {summary['fgw_sd_sample']:.3f} | {summary['lrot_sd_sample']:.3f} | {summary['delta_sd_sample']:.3f} |",
        f"| 95% paired t confidence interval | — | — | [{summary['delta_ci95_t'][0]:+.3f}, {summary['delta_ci95_t'][1]:+.3f}] |",
        f"| Median paired difference | — | — | {summary['delta_median']:+.3f} |",
        f"| Percentile-bootstrap 95% interval (sensitivity) | — | — | [{summary['delta_ci95_bootstrap_percentile_sensitivity'][0]:+.3f}, {summary['delta_ci95_bootstrap_percentile_sensitivity'][1]:+.3f}] |",
        f"| Range of paired differences | — | — | [{summary['delta_min']:+.3f}, {summary['delta_max']:+.3f}] |",
        f"| Positive / negative pairs | — | — | {summary['positive_delta_n']}/{summary['negative_delta_n']} |",
        "",
        "The historical 42/123 realization and three earlier exploratory pairs are excluded. The frozen synthetic subsets support technical reproduction and do not validate biological pathway labels.",
    ]
    (out_dir / "synthetic_seed_table.md").write_text("\n".join(table) + "\n", encoding="utf-8")
    fig, ax = plt.subplots(figsize=(7.2, 3.8), constrained_layout=True)
    mean = summary["delta_mean"]
    lo, hi = summary["delta_ci95_t"]
    ax.axhline(0, color="#777777", linewidth=1.0, linestyle="--", zorder=0)
    ax.axhspan(lo, hi, color="#D55E00", alpha=0.14, linewidth=0, label="95% CI of mean difference")
    ax.axhline(mean, color="#D55E00", linewidth=1.7, label=f"Mean difference = {mean:+.3f}")
    ax.scatter(np.arange(1, len(rows) + 1), delta_rows, s=34, color="#0072B2", edgecolor="white", linewidth=0.45, zorder=3)
    ax.set_xlabel("Prespecified replicate")
    ax.set_ylabel("Mapping-accuracy difference (LROT − FGW)")
    ax.set_xlim(0.5, len(rows) + 0.5)
    ax.set_xticks([1, 5, 10, 15, 20, 25, 30])
    low = min(float(np.min(delta_rows)), lo, 0.0)
    high = max(float(np.max(delta_rows)), hi, 0.0)
    pad = max(0.02, (high - low) * 0.10)
    ax.set_ylim(low - pad, high + pad)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#DDDDDD", linewidth=0.6, alpha=0.8)
    ax.legend(frameon=False, loc="best", fontsize=8)
    fig.savefig(out_dir / "supplementary_paired_difference.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    registry = {
        "n_primary_replicates": N_REPLICATES,
        "seed_pairs_A_B": [list(pair) for pair in SEED_PAIRS],
        "pairing": "FGW and LROT use the same generated A/B pair within each replicate",
        "excluded_historical_seed_pair": [42, 123],
        "excluded_earlier_exploratory_pairs": 3,
        "seed_pair_selection": "prespecified consecutive paired ranges; frozen before solving",
    }
    (out_dir / "synthetic_seed_registry.json").write_text(
        json.dumps(registry, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


def verify_only(out_dir: Path) -> None:
    """Verify current CSV/F3 values and archived JSONL in their own scopes.

    The public CSV intentionally has no timing fields. Archived observations are
    checked against their archived summary, not silently equated to revalidation.
    """
    _db, _feature_order, _groups, _weights, archive_check = load_design()
    registry = json.loads((out_dir / "synthetic_seed_registry.json").read_text(encoding="utf-8"))
    if registry["seed_pairs_A_B"] != [list(p) for p in SEED_PAIRS] or registry["n_primary_replicates"] != N_REPLICATES:
        raise RuntimeError("The registry does not match the fixed 30 seed pairs")
    stored = json.loads((out_dir / "historical_archive_reconstruction_check.json").read_text(encoding="utf-8"))
    if any(stored.get(label, {}).get(key) != value for label, checks in archive_check.items() for key, value in checks.items()):
        raise RuntimeError("Archived input reconstruction differs from the frozen design")
    with (out_dir / "Supplementary_File_S2_synthetic_seed_pairs.csv").open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"replicate", "seed_A", "seed_B", "fgw_accuracy", "lrot_accuracy", "delta_lrot_minus_fgw"}
        if not required.issubset(reader.fieldnames or []):
            raise RuntimeError("Current CSV is missing scientific columns")
        table_rows = list(reader)
    if len(table_rows) != N_REPLICATES:
        raise RuntimeError("Current CSV must contain exactly 30 rows")
    current = json.loads((ROOT / "figures/current_release/inputs/F3.json").read_text(encoding="utf-8"))
    pairs = {(int(p["seed_A"]), int(p["seed_B"])): p for p in current["pairs"]}
    if set(pairs) != set(SEED_PAIRS):
        raise RuntimeError("Current F3 payload has different seed pairs")
    rows = []
    for i, item in enumerate(table_rows, start=1):
        pair = (int(item["seed_A"]), int(item["seed_B"]))
        if pair != SEED_PAIRS[i - 1] or int(item["replicate"]) != i:
            raise RuntimeError(f"Unexpected seed/replicate order in row {i}")
        row = {"replicate": i, "seed_pair": list(pair), **{k: float(item[k]) for k in ["fgw_accuracy", "lrot_accuracy", "delta_lrot_minus_fgw"]}}
        if not all(np.isfinite(row[k]) for k in ["fgw_accuracy", "lrot_accuracy", "delta_lrot_minus_fgw"]):
            raise RuntimeError(f"Nonfinite value in row {i}")
        if not (0 <= row["fgw_accuracy"] <= 1 and 0 <= row["lrot_accuracy"] <= 1):
            raise RuntimeError(f"Accuracy outside [0,1] in row {i}")
        if not np.isclose(row["delta_lrot_minus_fgw"], row["lrot_accuracy"] - row["fgw_accuracy"], atol=1e-12, rtol=0):
            raise RuntimeError(f"Difference does not match accuracies in row {i}")
        for field, ref in [("fgw_accuracy", "FGW"), ("lrot_accuracy", "LROT")]:
            if not np.isclose(row[field], pairs[pair][ref], atol=1e-12, rtol=0):
                raise RuntimeError(f"Current CSV and F3 differ for {field} in row {i}")
        # Optional historical timing columns are checked only when supplied.
        for field in ["fgw_elapsed_seconds", "lrot_elapsed_seconds"]:
            if item.get(field) not in [None, ""] and (not np.isfinite(float(item[field])) or float(item[field]) < 0):
                raise RuntimeError(f"Invalid optional timing field in row {i}")
        rows.append(row)
    observed = summarize(rows)
    if not np.allclose([r["delta_lrot_minus_fgw"] for r in rows], current["paired_differences"], atol=1e-12, rtol=0):
        raise RuntimeError("Current differences disagree with F3")
    cs = current["summary"]
    for key, target in [("delta_mean", cs["mean_difference"]), ("delta_sd_sample", cs["sd_difference_ddof1"]), ("positive_delta_n", cs["positive"])]:
        if not np.isclose(observed[key], target, atol=1e-12, rtol=0):
            raise RuntimeError("Current F3 summary mismatch: " + key)
    if not np.allclose(observed["delta_ci95_t"], cs["ci95_t"], atol=1e-12, rtol=0):
        raise RuntimeError("Current F3 interval differs from CSV")
    raw = read_jsonl(out_dir / "synthetic_seed_results.jsonl")
    if set(raw) != set(SEED_PAIRS):
        raise RuntimeError("Archived checkpoint is incomplete or has different seeds")
    for pair, item in raw.items():
        if not np.isclose(item["delta_lrot_minus_fgw"], item["lrot_accuracy"] - item["fgw_accuracy"], atol=1e-12, rtol=0):
            raise RuntimeError("Archived checkpoint difference is inconsistent")
    archived = summarize([raw[p] for p in SEED_PAIRS])
    archived_summary = json.loads((out_dir / "synthetic_seed_summary.json").read_text(encoding="utf-8"))
    fields = ["n", "positive_delta_n", "zero_delta_n", "negative_delta_n", "fgw_mean", "fgw_sd_sample", "lrot_mean", "lrot_sd_sample", "delta_mean", "delta_sd_sample", "delta_median", "delta_min", "delta_max", "delta_ci95_t", "delta_ci95_bootstrap_percentile_sensitivity"]
    for key in fields:
        if not np.allclose(archived[key], archived_summary[key], atol=1e-12, rtol=0):
            raise RuntimeError("Archived summary mismatch: " + key)
    changed = sum(any(not np.isclose(raw[p][field], rows[i][field], atol=1e-12, rtol=0) for field in ["fgw_accuracy", "lrot_accuracy"]) for i, p in enumerate(SEED_PAIRS))
    print(f"Verified current CSV/F3: n={len(rows)}, mean difference {observed['delta_mean']:+.12f}, 95% CI {observed['delta_ci95_t']}; positive {observed['positive_delta_n']}, negative {observed['negative_delta_n']}.")
    print(f"Verified archived JSONL/summary separately; {changed} pair(s) differ from current revalidation. No solver executed; no experimental values modified.")


def run(resume: bool, out_dir: Path) -> None:
    db, feature_order, groups, weights, archive_check = load_design()
    out_dir.mkdir(parents=True, exist_ok=True)
    check_path = out_dir / "historical_archive_reconstruction_check.json"
    check_path.write_text(json.dumps(archive_check, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    result_path = out_dir / "synthetic_seed_results.jsonl"
    if result_path.exists() and not resume:
        raise RuntimeError(
            f"Checkpoint already exists at {result_path}; use --resume to continue it or choose a fresh --output-dir."
        )
    completed = read_jsonl(result_path) if resume else {}
    unknown = set(completed) - set(SEED_PAIRS)
    if unknown:
        raise RuntimeError(f"Checkpoint contains unregistered seed pairs: {sorted(unknown)}")
    for index, pair in enumerate(SEED_PAIRS, start=1):
        if pair in completed:
            print(f"[{index:02d}/{N_REPLICATES}] retaining existing checkpoint for {pair}", flush=True)
            continue
        A, B = make_pair(feature_order, groups, weights, *pair)
        fgw = solve(A, B, db, 0.0)
        lrot = solve(A, B, db, 0.1)
        item = {
            "replicate": index,
            "seed_pair": list(pair),
            "fgw_accuracy": fgw["accuracy"],
            "fgw_elapsed_seconds": fgw["elapsed_seconds"],
            "lrot_accuracy": lrot["accuracy"],
            "lrot_elapsed_seconds": lrot["elapsed_seconds"],
            "delta_lrot_minus_fgw": lrot["accuracy"] - fgw["accuracy"],
        }
        with result_path.open("a", encoding="utf-8", newline="") as handle:
            handle.write(json.dumps(item, ensure_ascii=False) + "\n")
            handle.flush()
        completed[pair] = item
        print(f"[{index:02d}/{N_REPLICATES}] A/B={pair[0]}/{pair[1]} "
              f"FGW={fgw['accuracy']:.3f} LROT={lrot['accuracy']:.3f} "
              f"difference={item['delta_lrot_minus_fgw']:+.3f}", flush=True)
    if set(completed) != set(SEED_PAIRS):
        raise RuntimeError("The prespecified run set is incomplete")
    rows = [completed[pair] for pair in SEED_PAIRS]
    write_derived(rows, archive_check, out_dir)
    print(f"Updated paired results, summary, table and figure under {out_dir}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Reproduce the Supplementary S23 30-pair seed-robustness analysis.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--verify-only", action="store_true", help="Validate the stored seed pairs and summary without solving.")
    group.add_argument("--run", action="store_true", help="Run missing registered pairs and regenerate all summary artifacts.")
    parser.add_argument("--resume", action="store_true", help="With --run, keep completed JSONL seed pairs and run only missing pairs.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT,
                        help="Output/checkpoint directory (default: lrot_output/s23_support).")
    parser.add_argument("--solver", choices=["revalidated", "legacy"], default="revalidated", help="Current bound solver by default; legacy is historical only.")
    args = parser.parse_args()
    out_dir = args.output_dir if args.output_dir.is_absolute() else ROOT / args.output_dir
    global core
    if args.verify_only:
        import importlib
        core = importlib.import_module("lrot_core")
        verify_only(out_dir)
    else:
        out_dir = out_dir.resolve()
        if out_dir == ROOT or ROOT in out_dir.parents:
            parser.error("New run output must be outside the code package; retained results are read-only")
        if args.solver == "revalidated":
            if args.resume or out_dir.exists():
                parser.error("Current revalidation requires a fresh output directory and no checkpoint reuse")
            from paper_protocol.binding import bind
            wrapper = bind(out_dir / "solver_binding")
            core = wrapper.core
            core.fgw_lr_solver = wrapper.solver
        else:
            import importlib
            core = importlib.import_module("lrot_core")
            print("Historical legacy solver selected; do not substitute these results for current revalidation.")
        run(resume=args.resume, out_dir=out_dir)


if __name__ == "__main__":
    main()
