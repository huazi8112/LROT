# -*- coding: utf-8 -*-
"""
run_cancer_visium_truecoord.py - Breast cancer LROT rerun with REAL 10x Visium coordinates
Plan A: replace the synthetic hexagonal-grid coordinates (old run_cancer_visium.py) with
the real tissue pixel coordinates from tissue_positions_list.csv for
  Section 1 -> Section 2  (V1_Breast_Cancer_Block_A_Section_1/2)

Pipeline mirrors run_dlpfc_visium.py real-coordinate conventions:
  - coordinates centered (mean subtracted)
  - Section 2 receives the standard 2-degree rotation + 0.5 translation misalignment
  - 3000 tissue spots per slice, sampled with seeds 42 (A) and 99 (B)

Outputs (new names; do not overwrite the old synthetic-coordinate results):
  lrot_output/lrot_cancer_truecoord_results.txt
  lrot_output/lrot_data_breast_cancer_truecoord.npz
"""
import os
import csv
import sys
import time
import warnings
warnings.filterwarnings('ignore')

import numpy as np
# 自定位仓库根：代码包解压到任意路径均可运行
_R = os.path.dirname(os.path.abspath(__file__))
while _R != os.path.dirname(_R) and not os.path.exists(os.path.join(_R, 'lrot_core.py')):
    _R = os.path.dirname(_R)


sys.path.insert(0, _R)
sys.path.insert(0, os.path.join(_R, 'experiments'))

from lrot_core import fgw_lr_solver, compute_ot_cost, compute_alignment
from run_cancer_visium import (HUMAN_LR_DB, HUMAN_ALL_LR,
                               apply_misalignment, compute_lr_strength_human)

DATA_DIR = os.environ.get('LROT_DATA_DIR', os.path.join(_R, 'cancer_data'))
OUT_DIR = os.environ.get('LROT_OUTPUT_DIR', os.path.join(_R, 'lrot_output'))
N_SAMPLE = 3000   # 按论文要求以 3,000 spots/slice 出图
GAMMA = 0.1
SEED_A = 42
SEED_B = 99
SAMPLES = ["V1_Breast_Cancer_Block_A_Section_1",
           "V1_Breast_Cancer_Block_A_Section_2"]

def _decode_h5_string(value):
    return value.decode("utf-8") if isinstance(value, (bytes, np.bytes_)) else str(value)

def _make_unique_feature_names(names):
    """Match AnnData's first-occurrence names plus numeric duplicate suffixes."""
    seen = set()
    next_suffix = {}
    unique = []
    for raw in names:
        name = str(raw)
        if name not in seen:
            seen.add(name)
            unique.append(name)
            continue
        n = next_suffix.get(name, 1)
        candidate = f"{name}-{n}"
        while candidate in seen:
            n += 1
            candidate = f"{name}-{n}"
        next_suffix[name] = n + 1
        seen.add(candidate)
        unique.append(candidate)
    return unique

def load_breast_slice(sample_id, target_spots=1000, seed=42):
    """Read a 10x filtered H5 and real tissue coordinates with barcode matching.

    The H5 matrix is decoded directly with h5py/scipy. This preserves the
    scanpy 10x semantics (gene-expression features, symbol names, barcode order)
    while avoiding a heavyweight scanpy import for this single input step.
    """
    import h5py
    from scipy.sparse import csc_matrix

    h5_path = os.path.join(DATA_DIR, sample_id + ".h5")
    pos_path = os.path.join(DATA_DIR, sample_id, "spatial",
                            "tissue_positions_list.csv")

    t0 = time.time()
    with h5py.File(h5_path, "r") as h5:
        matrix = h5["matrix"]
        all_barcodes = [_decode_h5_string(x) for x in matrix["barcodes"][:]]
        features = matrix["features"]
        all_names = [_decode_h5_string(x) for x in features["name"][:]]
        feature_types = [_decode_h5_string(x) for x in features["feature_type"][:]]
        gene_ix = np.asarray(
            [i for i, feature_type in enumerate(feature_types)
             if feature_type == "Gene Expression"], dtype=int)
        gene_names = _make_unique_feature_names([all_names[i] for i in gene_ix])
        shape = tuple(int(x) for x in matrix["shape"][:])
        counts = csc_matrix(
            (matrix["data"][:], matrix["indices"][:], matrix["indptr"][:]),
            shape=shape)
        # 10x stores features × barcodes in CSC order; transpose to cells × genes.
        expression_by_cell = counts[gene_ix, :].T.tocsr()
    del counts

    pos = {}
    with open(pos_path, newline="", encoding="utf-8") as f:
        rd = csv.reader(f)
        next(rd)
        for p in rd:
            if len(p) >= 6 and p[1] == "1":
                pos[p[0]] = (float(p[4]), float(p[5]))

    barcodes = [bc for bc in all_barcodes if bc in pos]
    if len(barcodes) == 0:
        raise RuntimeError(sample_id + ": no barcode matched to tissue positions")

    rng = np.random.RandomState(seed)
    n = min(target_spots, len(barcodes))
    idx = rng.choice(len(barcodes), n, replace=False)
    selected_barcodes = [barcodes[i] for i in idx]
    coords = np.array([pos[bc] for bc in barcodes], dtype=np.float64)
    coords = coords - coords.mean(axis=0)
    barcode_to_ix = {bc: i for i, bc in enumerate(all_barcodes)}
    selected_ix = np.asarray([barcode_to_ix[bc] for bc in selected_barcodes], dtype=int)

    print(f"  {sample_id}: H5 loaded and {len(barcodes)} barcodes matched "
          f"({time.time()-t0:.1f}s); densifying {n} sampled rows", flush=True)
    expr = expression_by_cell[selected_ix].toarray().astype(np.float32, copy=False)
    del expression_by_cell
    print(f"  {sample_id}: h5={len(all_barcodes)}, in-tissue matched={len(barcodes)}, "
          f"sampled={n} (seed={seed})", flush=True)
    print(f"    coord range x=[{coords[idx,0].min():.0f},{coords[idx,0].max():.0f}] "
          f"y=[{coords[idx,1].min():.0f},{coords[idx,1].max():.0f}]", flush=True)

    return {
        "coords": coords[idx],
        "expr": expr,
        "gene_names": gene_names,
        "region_labels": np.zeros(n, dtype=int),
        "n_spots": n,
        "n_genes": expr.shape[1],
        "sample_id": sample_id,
        "barcodes": selected_barcodes,
    }

def solve(slice_A, slice_B, S, gamma, lr_gene_list=None):
    return fgw_lr_solver(slice_A, slice_B, S, gamma=gamma, beta=0.3,
                         alpha=0.5, max_iter=200, verbose=False,
                         lr_gene_list=lr_gene_list)

def row_norm(P):
    return P / np.maximum(P.sum(1, keepdims=True), 1e-10)

def plan_stats(P):
    Pn = row_norm(P)
    ent_row = -np.sum(Pn * np.log(np.maximum(Pn, 1e-10)), axis=1)
    eff = np.exp(ent_row)
    return {
        "entropy": float(ent_row.mean()),
        "eff_median": float(np.median(eff)),
        "eff_mean": float(eff.mean()),
        "top1_mean": float(Pn.max(axis=1).mean()),
    }

def main():
    print("=" * 72)
    print("  LROT on Human Breast Cancer Visium - TRUE tissue coordinates")
    print(f"  {SAMPLES[0]}  ->  {SAMPLES[1]}")
    print(f"  {N_SAMPLE} spots/slice | gamma={GAMMA} | seeds {SEED_A}/{SEED_B}")
    print("=" * 72, flush=True)

    t0 = time.time()
    slice_A = load_breast_slice(SAMPLES[0], N_SAMPLE, SEED_A)
    slice_B = load_breast_slice(SAMPLES[1], N_SAMPLE, SEED_B)
    slice_B["coords"] = apply_misalignment(slice_B["coords"],
                                           rotation_deg=2.0,
                                           translation=0.5)
    print(f"  data load: {time.time()-t0:.1f}s", flush=True)

    S_lr, top_pairs = compute_lr_strength_human(slice_A, slice_B, HUMAN_LR_DB)
    S_zero = np.zeros_like(S_lr)
    genes_A = {str(g).upper() for g in slice_A['gene_names']}
    genes_B = {str(g).upper() for g in slice_B['gene_names']}
    matched_lr_pairs = sum(
        1 for lig, rec in HUMAN_LR_DB
        if lig.upper() in genes_A and rec.upper() in genes_B
    )
    print(f"  S_lr {S_lr.shape}  matched gene-symbol pairs "
          f"{matched_lr_pairs}/{len(HUMAN_LR_DB)} | top contributors: "
          + ", ".join(f"{l}-{r}" for l, r, _ in top_pairs[:6]), flush=True)

    t0 = time.time()
    # Matched-mask comparison: both methods exclude the same human LR genes from C_expr.
    # This isolates the effect of gamma/LR guidance while keeping expression features fixed.
    P_fgw, _ = solve(slice_A, slice_B, S_zero, 0.0,
                     lr_gene_list=HUMAN_ALL_LR)
    t_fgw = time.time() - t0
    m_fgw = plan_stats(P_fgw)
    m_fgw["cost"] = float(compute_ot_cost(P_fgw, slice_A["expr"], slice_B["expr"]))
    m_fgw["time"] = t_fgw
    e0 = m_fgw["entropy"]; c0 = m_fgw["cost"]
    em0 = m_fgw["eff_median"]; tp0 = m_fgw["top1_mean"]
    print(f"  FGW(gamma=0)  : {t_fgw:.1f}s  entropy={e0:.4f}  "
          f"cost={c0:.4f}  eff_targets(med)={em0:.0f}  top1={tp0:.3f}", flush=True)

    t0 = time.time()
    P_lrot, loss_hist = solve(slice_A, slice_B, S_lr, GAMMA,
                              lr_gene_list=HUMAN_ALL_LR)
    t_lrot = time.time() - t0
    m_lrot = plan_stats(P_lrot)
    m_lrot["cost"] = float(compute_ot_cost(P_lrot, slice_A["expr"], slice_B["expr"]))
    m_lrot["time"] = t_lrot
    e1 = m_lrot["entropy"]; c1 = m_lrot["cost"]
    em1 = m_lrot["eff_median"]; tp1 = m_lrot["top1_mean"]
    print(f"  LROT(gamma=0.1): {t_lrot:.1f}s  entropy={e1:.4f}  "
          f"cost={c1:.4f}  eff_targets(med)={em1:.0f}  top1={tp1:.3f}", flush=True)

    ent_imp = (e0 - e1) / e0 * 100
    cost_imp = (c0 - c1) / c0 * 100
    print(f"  Entropy reduction {ent_imp:.2f}% | Cost change {cost_imp:+.2f}%",
          flush=True)

    os.makedirs(OUT_DIR, exist_ok=True)
    txt_path = os.path.join(OUT_DIR, "lrot_cancer_truecoord_results.txt")
    with open(txt_path, "w", encoding="utf-8") as f:
        f.write("=" * 72 + "\n")
        f.write("  LROT on Human Breast Cancer Visium - TRUE tissue coordinates\n")
        f.write(f"  {SAMPLES[0]}  ->  {SAMPLES[1]}\n")
        f.write("  Real 10x tissue pixel coords (centered), matched expression mask "
                "(HUMAN_ALL_LR excluded in both arms), B misaligned "
                "(2deg + 0.5), " + str(N_SAMPLE) + " spots\n")
        f.write("=" * 72 + "\n\n")
        f.write(f"Slice A: {SAMPLES[0]} ({slice_A['n_spots']} spots, "
                f"{slice_A['n_genes']} genes, seed {SEED_A})\n")
        f.write(f"Slice B: {SAMPLES[1]} ({slice_B['n_spots']} spots, "
                f"{slice_B['n_genes']} genes, seed {SEED_B})\n\n")
        f.write(f"LR database pairs with gene symbols present in both matrices: "
                f"{matched_lr_pairs}/{len(HUMAN_LR_DB)}\n")
        f.write("Top LR pairs:\n")
        for i, (l, r, w) in enumerate(top_pairs[:8]):
            f.write(f"  {i+1}. {l}-{r} (weight={w:.2f})\n")
        f.write("\n")
        for tag, m in [("FGW (gamma=0)", m_fgw), ("LROT (gamma=0.1)", m_lrot)]:
            ee = m["entropy"]; cc = m["cost"]; tt = m["time"]
            emm = m["eff_median"]; emn = m["eff_mean"]; tpp = m["top1_mean"]
            f.write(tag + ":\n")
            f.write(f"  Entropy={ee:.4f}  Cost={cc:.4f}  Time={tt:.1f}s\n")
            f.write(f"  Effective targets/row: median={emm:.1f}  "
                    f"mean={emn:.1f}  top1={tpp:.3f}\n")
        f.write(f"\nEntropy reduction: {ent_imp:.2f}%\n")
        f.write(f"Cost change: {cost_imp:+.2f}%\n")
        f.write(f"LROT iterations: {len(loss_hist)}\n")
    print(f"  results -> {txt_path}", flush=True)

    aligned_coords = compute_alignment(P_lrot, slice_A["coords"], slice_B["coords"])
    npz_path = os.path.join(OUT_DIR, "lrot_data_breast_cancer_truecoord.npz")
    np.savez_compressed(
        npz_path,
        P=P_lrot, P_fgw=P_fgw,
        A_coords=slice_A["coords"], A_expr=slice_A["expr"].astype(np.float16),
        A_region_labels=slice_A["region_labels"],
        A_gene_names=np.array(slice_A["gene_names"], dtype=object),
        B_coords=slice_B["coords"], B_expr=slice_B["expr"].astype(np.float16),
        B_region_labels=slice_B["region_labels"],
        B_gene_names=np.array(slice_B["gene_names"], dtype=object),
        loss_hist=np.array(loss_hist), aligned_coords=aligned_coords,
        top_pairs=np.array(top_pairs, dtype=object),
        A_barcodes=np.array(slice_A["barcodes"], dtype=object),
        B_barcodes=np.array(slice_B["barcodes"], dtype=object),
    )
    print(f"  data -> {npz_path} ({os.path.getsize(npz_path)/1e6:.1f} MB)",
          flush=True)

    print("\n" + "=" * 72)
    print("  DONE (true-coordinate run). Figure/text update deferred until review.")
    print("=" * 72)

if __name__ == "__main__":
    main()
