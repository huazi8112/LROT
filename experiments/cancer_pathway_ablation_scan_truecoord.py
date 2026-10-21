# -*- coding: utf-8 -*-
"""
cancer_pathway_ablation_scan_truecoord.py — 乳腺癌通路级系统消融扫描（leave-one-out）
设计: 全部条件固定为 HUMAN_ALL_LR 非LR表达特征掩膜，仅改变 LR 引导库；与 Figure 6 使用同一表达特征口径。

做法: 固定 HUMAN_ALL_LR 表达特征掩膜和其余求解参数，只改变 LR 引导库；
      Δ 为移除家族后的熵降减去完整库熵降（百分点）。结果仅作描述性敏感性分析。

指标: 熵降%(相对FGW基线)。full 为完整库; 逐个移除后对比。
输出: lrot_output/lrot_cancer_ablation_scan_truecoord_results.txt, lrot_cancer_ablation_scan_truecoord.png
"""
import json
import os
import sys
import warnings
warnings.filterwarnings('ignore')
# 自定位仓库根：代码包解压到任意路径均可运行
_R = os.path.dirname(os.path.abspath(__file__))
while _R != os.path.dirname(_R) and not os.path.exists(os.path.join(_R, 'lrot_core.py')):
    _R = os.path.dirname(_R)

sys.path.insert(0, _R)
sys.path.insert(0, os.path.join(_R, 'experiments'))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from run_cancer_visium_truecoord import (load_breast_slice as download_and_prepare, apply_misalignment,
                               compute_lr_strength_human, HUMAN_LR_DB, HUMAN_ALL_LR)
from lrot_core import fgw_lr_solver, compute_ot_cost

OUT_DIR = os.environ.get('LROT_OUTPUT_DIR', os.path.join(_R, 'lrot_output'))
N = 3000
GAMMA = 0.1
SEED_A, SEED_B = 42, 99

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei']
plt.rcParams['axes.unicode_minus'] = False

buf = []
def log(s):
    print(s, flush=True)
    buf.append(s)

def mean_entropy(P):
    Pn = P / np.maximum(P.sum(1, keepdims=True), 1e-10)
    return -np.sum(Pn * np.log(np.maximum(Pn, 1e-10)), axis=1).mean()

def compute_ent0(A, B):
    """Fit the shared matched-mask FGW baseline used by every variant."""
    S0 = np.zeros((A['n_spots'], B['n_spots']))
    P0, _ = fgw_lr_solver(A, B, S0, gamma=0.0, beta=0.3, alpha=0.5,
                          max_iter=200, verbose=False, lr_gene_list=HUMAN_ALL_LR)
    return mean_entropy(P0)

def run_lrot(A, B, db):
    S, top = compute_lr_strength_human(A, B, db)
    # Keep expression features fixed to the full-library HUMAN_ALL_LR mask
    # for the full library and every leave-one-family-out variant.
    Pl, _ = fgw_lr_solver(A, B, S, gamma=GAMMA, beta=0.3, alpha=0.5,
                          max_iter=200, verbose=False, lr_gene_list=HUMAN_ALL_LR)
    return mean_entropy(Pl), top

def run_imp(ent0, entl):
    return (ent0 - entl) / ent0 * 100

# 互斥的信号家族定义（基于 HUMAN_LR_DB 的注释分组；EGF-EGFR/TGFA-EGFR 归入 ERBB2/HER 家族）
FAMILIES = {
    'ERBB2/HER': {("ERBB2","ERBB2"), ("ERBB2","ERBB3"), ("EREG","EGFR"), ("BTC","EGFR"),
                  ("HBEGF","EGFR"), ("NRG1","ERBB3"), ("NRG1","ERBB4"), ("EGF","EGFR"), ("TGFA","EGFR")},
    'Neurotrophin/Growth': {("NTNG1","NTRK2"), ("BDNF","NTRK2"), ("NTF3","NTRK3"), ("NGF","NTRK1"),
                            ("GDNF","GFRA1"), ("ARTN","GFRA3"), ("NRTN","GFRA2"), ("FGF8","FGFR1"),
                            ("FGF9","FGFR2"), ("FGF10","FGFR2"), ("FGF2","FGFR1"), ("HGF","MET"),
                            ("VEGFA","KDR"), ("VEGFB","FLT1"), ("PDGFA","PDGFRA"), ("PDGFB","PDGFRB"),
                            ("IGF1","IGF1R"), ("IGF2","IGF1R")},
    'Wnt': {("WNT3A","FZD1"), ("WNT5A","FZD5"), ("WNT7A","FZD10"), ("WNT1","FZD1")},
    'Notch': {("DLL1","NOTCH1"), ("DLL4","NOTCH4"), ("JAG1","NOTCH1"), ("JAG2","NOTCH2")},
    'Ephrin': {("EFNB1","EPHB2"), ("EFNA1","EPHA4"), ("EFNB2","EPHB4"), ("EFNA5","EPHA3")},
    'Semaphorin': {("SEMA3A","NRP1"), ("SEMA3F","NRP2"), ("SEMA4D","PLXNB1"), ("SEMA6A","PLXNA2")},
    'Slit/Robo': {("SLIT1","ROBO1"), ("SLIT2","ROBO2"), ("SLIT3","ROBO2")},
    'Chemokine/Cytokine': {("CXCL12","CXCR4"), ("CXCL13","CXCR5"), ("CCL2","CCR2"), ("CCL5","CCR5"),
                           ("IL1B","IL1R1"), ("TNF","TNFRSF1A"), ("TGFB1","TGFBR1"), ("TGFB2","TGFBR2")},
    'BMP': {("BMP4","BMPR1A"), ("BMP7","BMPR2"), ("BMP2","ACVR1")},
    'Cell Adhesion': {("CDH1","CDH1"), ("CDH2","CDH2"), ("NCAM1","NCAM1"), ("L1CAM","L1CAM"),
                      ("COL1A1","ITGA2"), ("VTN","ITGAV")},
    'Hedgehog': {("SHH","PTCH1"), ("IHH","PTCH2"), ("DHH","PTCH1")},
    'Reelin': {("RELN","LRP8"), ("RELN","VLDLR")},
}

_flat_family_pairs = [pair for family_pairs in FAMILIES.values() for pair in family_pairs]
if len(_flat_family_pairs) != len(set(_flat_family_pairs)) or set(_flat_family_pairs) != set(HUMAN_LR_DB):
    raise ValueError("Family definitions must be a disjoint partition of HUMAN_LR_DB.")

def main():
    log("=" * 72)
    log("  乳腺癌通路级系统消融扫描（leave-one-out）")
    log("=" * 72)
    A = download_and_prepare("V1_Breast_Cancer_Block_A_Section_1", N, SEED_A)
    B = download_and_prepare("V1_Breast_Cancer_Block_A_Section_2", N, SEED_B)
    B['coords'] = apply_misalignment(B['coords'], rotation_deg=2.0, translation=0.5)

    # Reuse the exact matched-mask FGW/full-library checkpoint just produced
    # by the ERBB2 ablation when the fixed dataset, seed and solver settings match.
    reference_path = os.environ.get(
        'LROT_MATCHED_REFERENCE_JSON',
        os.path.join(OUT_DIR, 'lrot_cancer_erbb2_ablation_truecoord_results.json'))
    reference = None
    if os.path.isfile(reference_path):
        with open(reference_path, encoding='utf-8') as f:
            reference = json.load(f)
        if (reference.get('expression_mask') !=
                'HUMAN_ALL_LR (frozen for matched FGW and all LROT variants)'
                or reference.get('sample_size_per_section') != N
                or reference.get('seeds') != {'A': SEED_A, 'B': SEED_B}
                or abs(float(reference.get('gamma', -1)) - GAMMA) > 1e-12):
            raise ValueError('Matched-mask reference does not match this scan configuration.')
        ent0 = float(reference['fgw_entropy'])
        entl_full = float(reference['variants']['full']['entropy'])
        _, full_top = compute_lr_strength_human(A, B, HUMAN_LR_DB)
        # Report a portable source path without exposing a personal checkout location.
        from pathlib import Path
        try:
            reference_display = Path(reference_path).resolve().relative_to(Path(_R).resolve()).as_posix()
        except ValueError:
            reference_display = f"external checkpoint: {Path(reference_path).name}"
        log(f"  Reused matched-mask FGW/full-library checkpoints from {reference_display}")
    else:
        ent0 = compute_ent0(A, B)
        entl_full, full_top = run_lrot(A, B, HUMAN_LR_DB)
    log(f"  Matched-mask FGW entropy: {ent0:.6f}")
    full_imp = run_imp(ent0, entl_full)
    log(f"  完整LR库({len(HUMAN_LR_DB)}对): 熵降 {full_imp:.4f}%")
    log(f"  Top LR 信号: {[f'{l}-{r}' for l,r,w in full_top[:5]]}")

    # 家族级消融: 移除每个信号家族
    log("\n  家族级消融: 移除每个信号家族后的熵降% (Δ=熵降-完整库)")
    log(f"  {'移除的家族':<24} {'对数':<6} {'熵降%':<10} {'Δ':<10} {'判定':<10}")
    log(f"  {'-'*62}")
    rows = []
    family_details = []
    for fname, pairs in FAMILIES.items():
        db = {k: v for k, v in HUMAN_LR_DB.items() if k not in pairs}
        n_removed = sum(1 for k in HUMAN_LR_DB if k in pairs)
        entl, _ = run_lrot(A, B, db)
        imp = run_imp(ent0, entl)
        d = imp - full_imp
        verdict = ('removal lowers' if d < -0.1 else
                   ('removal raises' if d > 0.1 else 'within ±0.1 pp'))
        rows.append((fname, n_removed, imp, d, verdict))
        family_details.append({
            "family": fname, "pairs_removed": n_removed, "entropy": entl,
            "reduction_percent": imp, "delta_percentage_points": d,
            "descriptive_band": verdict,
        })
        log(f"  {fname:<24} {n_removed:<6} {imp:<10.3f} {d:<+10.3f} {verdict}")

    # 排序
    rows_sorted = sorted(rows, key=lambda x: x[3])
    log("\n" + "=" * 72)
    log("  消融谱（Δ=移除家族后的熵降变化，单位为百分点；仅作描述）")
    log("=" * 72)
    for fname, npair, imp, d, v in rows_sorted:
        log(f"  {fname:<24} 熵降{imp:.3f}% Δ={d:+.3f} pp {v}")
    lowers = [r for r in rows_sorted if r[3] < -0.1]
    raises = [r for r in rows_sorted if r[3] > 0.1]
    within = len(rows) - len(lowers) - len(raises)
    log(f"\n  Δ < -0.1 pp: {len(lowers)} families; Δ > +0.1 pp: {len(raises)}; within band: {within}")
    log("  This is a descriptive sensitivity analysis; the ±0.1 percentage-point band is not a")
    log("  statistical test and does not establish a causal or mechanistic pathway contribution.")

    # 保存
    os.makedirs(OUT_DIR, exist_ok=True)
    result_txt = os.path.join(OUT_DIR, 'lrot_cancer_ablation_scan_truecoord_results.txt')
    with open(result_txt, 'w', encoding='utf-8') as f:
        f.write("\n".join(buf) + "\n")
    result_json = {
        "analysis": "matched-mask leave-one-family-out sensitivity",
        "expression_mask": "HUMAN_ALL_LR (frozen for FGW and all LROT variants)",
        "expression_mask_symbols": len(HUMAN_ALL_LR),
        "sample_size_per_section": N, "seeds": {"A": SEED_A, "B": SEED_B},
        "gamma": GAMMA, "fgw_entropy": ent0, "full_library_entropy": entl_full,
        "full_library_reduction_percent": full_imp,
        "delta_definition": "reduction after family removal minus full-library reduction, percentage points",
        "descriptive_band_pp": 0.1,
        "family_ablations": family_details,
    }
    with open(os.path.join(OUT_DIR, 'lrot_cancer_ablation_scan_truecoord_results.json'), 'w', encoding='utf-8') as f:
        json.dump(result_json, f, ensure_ascii=False, indent=2)
        f.write("\n")

    # 图
    fig, ax = plt.subplots(figsize=(9, 4.6))
    labels = [r[0] for r in rows_sorted]
    ds = [r[3] for r in rows_sorted]
    colors = ['#2166AC' if d < -0.1 else ('#D6604D' if d > 0.1 else '#888888') for d in ds]
    bars = ax.bar(labels, ds, color=colors, alpha=0.85, width=0.6)
    for b, d in zip(bars, ds):
        ax.text(b.get_x() + b.get_width()/2, d + (0.03 if d >= 0 else -0.06), f'{d:+.2f}',
                ha='center', fontsize=8)
    ax.axhline(0, color='black', lw=0.8)
    ax.axhspan(-0.1, 0.1, color='#cccccc', alpha=0.35, zorder=0)
    ax.text(0.99, 0.97, f'Matched-mask full-library reduction: {full_imp:.3f}%',
            transform=ax.transAxes, fontsize=8, color='gray', ha='right', va='top')
    ax.set_xlabel('Removed signaling family (leave-one-family-out)')
    ax.set_ylabel('Change in entropy reduction (percentage points)')
    ax.set_title('(B) Breast-cancer family ablation (descriptive; fixed expression mask)')
    plt.xticks(rotation=45, ha='right')
    ax.grid(True, axis='y', alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, 'lrot_cancer_ablation_scan_truecoord.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)
    log(f"\n  图已保存: {os.path.join(OUT_DIR, 'lrot_cancer_ablation_scan_truecoord.png')}")

if __name__ == '__main__':
    main()
