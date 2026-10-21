from __future__ import annotations
import ast, hashlib, json, platform, sys, time, argparse
from pathlib import Path
import numpy as np
import scipy
from scipy.spatial.distance import cdist
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'experiments' / 'lrot_prototype.py'
DESIGN = ROOT / 'experiments' / 'synthetic_design_manifest.json'
ARCHIVE = next((ROOT / 'lrot_output').rglob('lrot_data_synthetic.npz'))
APPROVED = ROOT / 'figures' / 'current_release' / 'assets' / 'F2.png'
parser = argparse.ArgumentParser(description='Recompute the adopted four 80-update diagnostic histories with fixed original parameters.')
parser.add_argument('--output-dir', type=Path, required=True)
args = parser.parse_args()
OUT = args.output_dir.resolve()
if OUT.exists() or OUT == ROOT or ROOT in OUT.parents:
    parser.error('Output must be a new directory outside the code package')
OUT.mkdir(parents=True, exist_ok=False)
src_bytes = SOURCE.read_bytes()
src_sha = hashlib.sha256(src_bytes).hexdigest()

def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()

def arr_sha(a) -> str:
    a = np.ascontiguousarray(a)
    return sha(a.tobytes(order="C"))

tree = ast.parse(src_bytes.decode("utf-8-sig"), filename=str(SOURCE))
keep_functions = {
    "lr_gene_list", "simulate_slice", "gene_to_idx",
    "compute_lr_strength_matrix", "cosine_distance",
    "euclidean_distance_matrix", "gromov_cost_matrix",
    "sinkhorn_stabilized", "fgw_lr_solver",
}
keep_globals = {"SEED", "LR_DB", "N_REGIONS", "REGION_NAMES"}
ns = {"np": np, "Path": Path, "cdist": cdist, "__file__": str(SOURCE), "__name__": "_f2_current_source_subset"}
base_nodes = []
solver_nodes = []
def targets(node):
    out = []
    for t in getattr(node, "targets", []):
        if isinstance(t, ast.Name):
            out.append(t.id)
        elif isinstance(t, (ast.Tuple, ast.List)):
            out.extend(x.id for x in t.elts if isinstance(x, ast.Name))
    return out

for node in tree.body:
    if isinstance(node, ast.FunctionDef) and node.name == "lr_gene_list":
        base_nodes.append(node)
    elif isinstance(node, (ast.Assign, ast.AnnAssign)) and set(targets(node)) & keep_globals:
        base_nodes.append(node)
    elif isinstance(node, ast.FunctionDef) and node.name in keep_functions - {"lr_gene_list"}:
        solver_nodes.append(node)
for node in base_nodes:
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(SOURCE), "exec"), ns)

design_bytes = DESIGN.read_bytes()
design = json.loads(design_bytes.decode("utf-8"))
ns["SYNTHETIC_DESIGN"] = design
ligands, receptors = ns["lr_gene_list"](ns["LR_DB"])
historical_count = len(set(ligands + receptors))
ns["ALL_LR_GENES"] = list(design["feature_order"][:historical_count])
ns["SYNTHETIC_REGION_ENRICHMENT"] = design["region_enrichment"]
if len(ns["ALL_LR_GENES"]) != historical_count or set(ns["ALL_LR_GENES"]) != set(ligands + receptors):
    raise RuntimeError("Current feature order does not match the LR gene set")
for region, sets in ns["SYNTHETIC_REGION_ENRICHMENT"].items():
    if (set(sets["ligands"]) | set(sets["receptors"])) - set(ns["ALL_LR_GENES"]):
        raise RuntimeError(f"Unknown LR gene in synthetic design region {region}")
for node in solver_nodes:
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(SOURCE), "exec"), ns)

slice_A = ns["simulate_slice"](1000, 200, ns["ALL_LR_GENES"], seed=42)
slice_B = ns["simulate_slice"](1000, 200, ns["ALL_LR_GENES"], seed=123)
theta = 0.15
c, s = np.cos(theta), np.sin(theta)
R = np.array([[c, -s], [s, c]])
center = slice_B["coords"].mean(0)
slice_B["coords"] = (slice_B["coords"] - center) @ R.T + center + np.array([0.5, 0.3])
S_lr, top_lr_pairs = ns["compute_lr_strength_matrix"](slice_A, slice_B, ns["LR_DB"])

input_arrays = {
    "A_coords": slice_A["coords"], "A_expr_full_precision": slice_A["expr"],
    "A_region_labels": slice_A["region_labels"], "A_gene_names": np.asarray(slice_A["gene_names"]),
    "B_coords": slice_B["coords"], "B_expr_full_precision": slice_B["expr"],
    "B_region_labels": slice_B["region_labels"], "B_gene_names": np.asarray(slice_B["gene_names"]),
    "S_lr_reconstructed": S_lr,
}
archive_checks = {}
alias_map = design.get("gene_symbol_aliases", {})
with np.load(ARCHIVE, allow_pickle=False) as z:
    for saved in ["A_coords", "A_expr", "A_region_labels", "A_gene_names",
                  "B_coords", "B_expr", "B_region_labels", "B_gene_names"]:
        if saved in ("A_expr", "B_expr"):
            current = input_arrays[saved + "_full_precision"].astype(z[saved].dtype)
            saved_values = z[saved]
            alias_changes = []
        elif saved in ("A_gene_names", "B_gene_names"):
            current_values = input_arrays[saved].astype(str).tolist()
            saved_values_list = z[saved].astype(str).tolist()
            canonical_saved = [alias_map.get(x, x) for x in saved_values_list]
            current = np.asarray(current_values)
            saved_values = np.asarray(canonical_saved)
            alias_changes = [{"index": i, "saved": old, "current": new}
                             for i,(old,new) in enumerate(zip(saved_values_list, canonical_saved)) if old != new]
        else:
            current = input_arrays[saved]
            saved_values = z[saved]
            alias_changes = []
        archive_checks[saved] = {
            "shape_match": list(current.shape) == list(saved_values.shape),
            "equal_at_saved_precision_or_manifest_aliases": bool(np.array_equal(current, saved_values)),
            "documented_alias_replacements": alias_changes,
            "saved_raw_array_sha256": arr_sha(z[saved]),
            "current_comparison_array_sha256": arr_sha(current),
            "saved_dtype": str(z[saved].dtype),
        }
    saved_main_loss = np.array(z["loss_hist"], copy=True)
    main_length = int(saved_main_loss.size)
failed = [k for k,v in archive_checks.items()
          if not v["shape_match"] or not v["equal_at_saved_precision_or_manifest_aliases"]]
if failed:
    raise RuntimeError("Input identity check failed: " + ", ".join(failed))
if main_length != 123:
    raise RuntimeError(f"Expected saved main loss history length 123, observed {main_length}")

gammas = [0.0, 0.05, 0.1, 0.2]
curves = {}
run_rows = []
for gamma in gammas:
    started = time.perf_counter()
    plan, history = ns["fgw_lr_solver"](
        slice_A, slice_B, S_lr,
        gamma=gamma, beta=0.3, alpha=0.5,
        max_iter=80, verbose=False
    )
    elapsed = time.perf_counter() - started
    curve = np.asarray(history, dtype=np.float64)
    key = f"gamma_{gamma:.2f}"
    curves[key] = curve
    run_rows.append({
        "gamma": gamma, "alpha": 0.5, "beta": 0.3,
        "reg_gw": 0.01, "reg_ot": 0.01, "max_iter": 80,
        "tol": 1e-5, "sinkhorn_max_iter": 200, "sinkhorn_tol": 1e-4,
        "outer_updates": int(curve.size), "first_loss": float(curve[0]),
        "last_loss": float(curve[-1]), "plan_sha256": arr_sha(plan),
        "loss_history_sha256": arr_sha(curve), "elapsed_seconds": elapsed,
    })

npz_out = OUT / "F2_four_gamma_histories_current_code.npz"
np.savez_compressed(npz_out, saved_main_gamma_0_10=saved_main_loss, **curves)
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(saved_main_loss, "b-", linewidth=1.5, label=r"LROT ($\gamma=0.1$)")
colors = ["gray", "orange", "green", "red"]
for i, gamma in enumerate(gammas):
    ax.plot(curves[f"gamma_{gamma:.2f}"], color=colors[i], linewidth=1.0, alpha=0.7,
            label=rf"$\gamma={gamma:.2f}$")
ax.set_title("F. Loss Curves (Ablation)", fontsize=11)
ax.set_xlabel("Iteration")
ax.set_ylabel("Loss")
ax.legend(fontsize=8)
fig.tight_layout()
png_out = OUT / "F2_panel_F_current_code_reconstruction.png"
fig.savefig(png_out, dpi=300, bbox_inches="tight")
plt.close(fig)

def file_sha(path):
    return sha(Path(path).read_bytes())
report = {
    "status": "completed_current_code_reconstruction",
    "interpretation": "Four histories were newly computed from current source. The archive confirms positional inputs at its stored precision; two gene symbols differ only by the exact aliases declared in the design manifest. This does not authenticate the absent historical source/environment snapshot.",
    "source": {"path": str(SOURCE), "sha256": src_sha},
    "design": {"path": str(DESIGN), "sha256": file_sha(DESIGN), "seed_pair": design["seed_pair"],
               "aliases_used_for_identity_check": alias_map},
    "saved_archive": {
        "path": str(ARCHIVE), "sha256": file_sha(ARCHIVE), "main_loss_length": main_length,
        "contains_four_ablation_loss_histories": False,
        "archive_input_checks": archive_checks,
        "expression_caveat": "Archived expression arrays are float16; equality is checked after casting current full-precision generated expression to that stored precision."
    },
    "approved_current_figure": {"path": str(APPROVED), "sha256": file_sha(APPROVED)},
    "source_protocol": {
        "synthetic_shape_each_slice": [1000, 200], "slice_seeds": [42, 123],
        "B_rotation_radians": 0.15, "B_translation": [0.5, 0.3],
        "gamma_values": gammas, "alpha": 0.5, "beta": 0.3,
        "reg_gw": 0.01, "reg_ot": 0.01, "outer_max_iter": 80, "outer_tol": 1e-5,
        "sinkhorn_max_iter": 200, "sinkhorn_tol": 1e-4,
        "stopping_rule": "relative plan change < tol and outer iteration index > 5",
    },
    "environment": {
        "python": sys.version, "platform": platform.platform(),
        "numpy": np.__version__, "scipy": scipy.__version__, "matplotlib": matplotlib.__version__,
    },
    "input_array_sha256": {k: arr_sha(v) for k,v in input_arrays.items()},
    "runs": run_rows,
    "generated": {
        "histories_npz": {"path": npz_out.name, "sha256": file_sha(npz_out)},
        "panel_png": {"path": png_out.name, "sha256": file_sha(png_out)},
        "plot_scope": "standalone panel F only; not a byte-identical replacement for the six-panel approved PNG",
    },
}
json_out = OUT / "F2_reconstruction_record.json"
json_out.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
md = [
    "# F2 four-curve reconstruction from current code",
    "",
    "Status: completed as a new current-code calculation; not authenticated as the missing historical run.",
    "",
    f"- Source SHA-256: {src_sha}",
    f"- Synthetic design SHA-256: {file_sha(DESIGN)}",
    f"- Existing archive SHA-256: {file_sha(ARCHIVE)} (main loss history: {main_length} points; no four ablation histories)",
    "- Input identity: coordinates, labels, feature positions, and expression values matched the archive at stored precision; the two archived symbols Artemin/Neurturin map to the current names Artn/Nrtn using the manifest's explicit aliases.",
    "- Current code calls: gamma=0/0.05/0.1/0.2; alpha=.5, beta=.3, max_iter=80, tol=1e-5; Sinkhorn max_iter=200, tol=1e-4.",
    "- Output is a new standalone F panel. Original approved image, manuscript, and experiment archive were not changed.",
    "- It does not prove the historical plotted curves came from this exact source/environment; original per-run logs/curve arrays and historical source/environment hashes remain absent.",
    "",
    "| gamma | updates | first loss | final loss | curve SHA-256 |",
    "|---:|---:|---:|---:|---|",
]
for row in run_rows:
    md.append(f"| {row['gamma']:.2f} | {row['outer_updates']} | {row['first_loss']:.12g} | {row['last_loss']:.12g} | {row['loss_history_sha256']} |")
(OUT / "F2_reconstruction_notes.md").write_text("\n".join(md) + "\n", encoding="utf-8")
print(json.dumps({"report": str(json_out), "npz": str(npz_out), "png": str(png_out), "runs": run_rows}, ensure_ascii=False))