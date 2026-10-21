"""
run_real_visium.py — 在真实 10x Visium 数据上运行 LROT 对齐
读取 h5 表达矩阵 + 采样 spot + 生成网格坐标
"""
import numpy as np
import time
import sys
import os
# 自定位仓库根：代码包解压到任意路径均可运行
_R = os.path.dirname(os.path.abspath(__file__))
while _R != os.path.dirname(_R) and not os.path.exists(os.path.join(_R, 'lrot_core.py')):
    _R = os.path.dirname(_R)

sys.path.insert(0, _R)

import argparse
import hashlib
import json
import platform
from pathlib import Path

parser = argparse.ArgumentParser(description='Fixed mouse-expression perturbation experiment; generated coordinates, no region annotations.')
parser.add_argument('--input-h5', type=Path, default=Path(_R) / 'cancer_data' / 'mouse_brain.h5')
parser.add_argument('--output-dir', type=Path, required=True, help='New output directory outside the code package; never overwrites retained results.')
parser.add_argument('--solver', choices=['legacy', 'revalidated'], default='revalidated', help='Revalidated uses the paper source binding and archived float16 expression precision; legacy preserves the original script computation.')
args = parser.parse_args()
H5_PATH = args.input_h5.resolve()
out_dir = args.output_dir.resolve()
if out_dir == Path(_R) or Path(_R) in out_dir.parents:
    parser.error('output directory must be outside the code package')
if out_dir.exists():
    parser.error('output directory must not already exist')
if not H5_PATH.is_file():
    parser.error('input H5 not found; supply --input-h5')
def file_sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1048576), b''):
            h.update(block)
    return h.hexdigest()
expected_hash = 'fdae7180e8a2727737429e492053cd1f1900eece3a31807516f0ccdc4146430e'
input_hash = file_sha256(H5_PATH)
if input_hash != expected_hash:
    parser.error('input H5 SHA-256 differs from the frozen mouse protocol')
out_dir.mkdir(parents=True, exist_ok=False)
if args.solver == 'revalidated':
    from paper_protocol.binding import bind
    wrapper = bind(out_dir / 'solver_binding')
    from lrot_core import *
    fgw_lr_solver = wrapper.solver
else:
    from lrot_core import *
import lrot_core as selected_core
import scanpy as sc

IS_MOUSE = True

print("=" * 60)
print("  LROT on 10x Visium Mouse Brain (Real Data)")
print("=" * 60)

print("\n[1/4] 读取真实 Visium 表达矩阵...")
adata = sc.read_10x_h5(H5_PATH)
if not adata.var_names.is_unique:
    adata.var_names_make_unique()
expr = adata.X.toarray().astype(np.float32)
gene_names = list(adata.var_names)
n_spots = len(expr)
print(f"  原始: {n_spots} spots, {expr.shape[1]} genes")

# 固定协议：seed 42，采样 1000 spots
N_SAMPLE = 1000
rng = np.random.RandomState(42)
idx = rng.choice(n_spots, N_SAMPLE, replace=False)
expr = expr[idx]
import csv
with open(out_dir / 'ordered_barcodes.csv', 'w', newline='', encoding='utf-8') as handle:
    writer = csv.writer(handle)
    writer.writerow(['row', 'raw_index', 'barcode'])
    writer.writerows((i, int(j), str(adata.obs_names[j])) for i, j in enumerate(idx))
print(f"  采样后: {N_SAMPLE} spots")

print("\n[2/4] 生成配对切片...")
side = int(np.ceil(np.sqrt(N_SAMPLE)))
coords = []
for i in range(side):
    for j in range(side):
        x = i * np.sqrt(3)
        y = j * 1.5 + (0.75 if i % 2 == 1 else 0.0)
        coords.append([x, y])
coords = np.array(coords[:N_SAMPLE]) + rng.randn(N_SAMPLE, 2) * 0.1

slice_A = dict(coords=coords.astype(float), expr=expr, gene_names=gene_names[:],
               region_labels=np.zeros(N_SAMPLE, dtype=int),
               n_spots=N_SAMPLE, n_genes=expr.shape[1])

center = coords.mean(0)
theta = np.deg2rad(12.0)
c, s = np.cos(theta), np.sin(theta)
M = np.array([[c, -s], [s, c]]) @ np.diag([0.97, 0.97])
coords_B = (coords - center) @ M.T + center
expr_B = expr.copy().astype(float) + rng.gamma(0.5, 0.5, expr.shape) * 0.2
expr_B *= rng.uniform(0.9, 1.1, N_SAMPLE)[:, None]
slice_B = dict(coords=coords_B.astype(float), expr=expr_B.astype(np.float32),
               gene_names=gene_names[:], region_labels=np.zeros(N_SAMPLE, dtype=int),
               n_spots=N_SAMPLE, n_genes=expr.shape[1])

# The paper revalidation used the retained float16 expression arrays promoted to float32.
if args.solver == 'revalidated':
    for data in (slice_A, slice_B):
        data['expr'] = data['expr'].astype(np.float16).astype(np.float32)
print("\n[3/4] 计算 LR 强度矩阵...")
t0 = time.time()
# 根据数据源选择 LR 数据库
if IS_MOUSE:
    lr_db_to_use = LR_DB
    lr_genes_to_use = ALL_LR_GENES
else:
    # 人类数据使用大写 LR 数据库
    HUMAN_LR_DB = {(k[0].upper(), k[1].upper()): v for k, v in LR_DB.items()}
    # 额外添加乳腺癌相关通路
    HUMAN_LR_DB.update({
        ("ERBB2", "ERBB2"): 0.98,
        ("ERBB2", "ERBB3"): 0.92,
        ("EGF", "EGFR"): 0.95,
        ("TGFA", "EGFR"): 0.88,
        ("HGF", "MET"): 0.90,
        ("VEGFA", "KDR"): 0.92,
        ("NRG1", "ERBB3"): 0.90,
        ("CCL2", "CCR2"): 0.85,
        ("TGFB1", "TGFBR1"): 0.92,
        ("COL1A1", "ITGA2"): 0.85,
    })
    lr_db_to_use = HUMAN_LR_DB
    lr_genes_to_use = sorted(set([k[0] for k in HUMAN_LR_DB] + [k[1] for k in HUMAN_LR_DB]))  # sorted: 避免 PYTHONHASHSEED 跨进程顺序不同
S_lr, top_pairs = compute_lr_strength_matrix(slice_A, slice_B, lr_db_to_use)
print(f"  LR矩阵: {S_lr.shape}, 范围 [{S_lr.min():.4f}, {S_lr.max():.4f}]")
print(f"  Top LR: {[(l, r) for l, r, w in top_pairs[:6]]}")
print(f"  耗时: {time.time()-t0:.2f}s")

print(f"\n[4a/4] LROT (gamma=0.1)...")
t0 = time.time()
if not IS_MOUSE:
    P, loss_hist = fgw_lr_solver(slice_A, slice_B, S_lr, gamma=0.1,
        beta=0.3, alpha=0.5, max_iter=200, verbose=True,
        lr_gene_list=lr_genes_to_use)
else:
    P, loss_hist = fgw_lr_solver(slice_A, slice_B, S_lr, gamma=0.1,
        beta=0.3, alpha=0.5, max_iter=200, verbose=True)
print(f"  耗时: {time.time()-t0:.2f}s")

print(f"\n[4b/4] FGW baseline (gamma=0)...")
P_fgw, _ = fgw_lr_solver(slice_A, slice_B, np.zeros_like(S_lr), gamma=0.0,
    beta=0.3, alpha=0.5, max_iter=200, verbose=False)

lrot_cost = compute_ot_cost(P, slice_A['expr'], slice_B['expr'])
Pn = P / np.maximum(P.sum(1, keepdims=True), 1e-10)
lrot_ent = -np.sum(Pn * np.log(np.maximum(Pn, 1e-10)), axis=1).mean()
fgw_cost = compute_ot_cost(P_fgw, slice_A['expr'], slice_B['expr'])
P0n = P_fgw / np.maximum(P_fgw.sum(1, keepdims=True), 1e-10)
fgw_ent = -np.sum(P0n * np.log(np.maximum(P0n, 1e-10)), axis=1).mean()

print("\n" + "=" * 60)
print("  Results Comparison (Real Visium Data)")
print("=" * 60)
print("{:<22} {:<12} {:<12}".format("Method", "Cost", "Entropy"))
print("-" * 46)
print("{:<22} {:<12.4f} {:<12.4f}".format("FGW (gamma=0)", fgw_cost, fgw_ent))
print("{:<22} {:<12.4f} {:<12.4f}".format("LROT (gamma=0.1)", lrot_cost, lrot_ent))
print("=" * 60)

print("\nGenerating visualization...")
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
# Output directory was validated before input preparation.

# ===== 保存真实 Visium 结果 =====
with open(out_dir / 'lrot_real_visium_results.txt', 'w', encoding='utf-8') as f:
    f.write("=" * 72 + "\n")
    f.write("  LROT on Real 10x Visium Mouse Brain - RESULTS\n")
    f.write("=" * 72 + "\n")
    f.write(f"Data: {N_SAMPLE} sampled spots (from {n_spots}), {expr.shape[1]} genes\n")
    f.write("Transform: rotation=12.0 deg, scale=0.97, batch noise=10%\n\n")
    f.write("  Method            Cost     Entropy\n")
    f.write("  " + "-" * 40 + "\n")
    f.write(f"  FGW (gamma=0)   {fgw_cost:.4f}   {fgw_ent:.4f}\n")
    f.write(f"  LROT (gamma=0.1){lrot_cost:.4f}   {lrot_ent:.4f}\n\n")
    # 正文引用的是相对基线的降幅；4 位小数显示会掩盖真实比值（如 0.0070/0.0073
    # 的比值给出 4.7% 而非 4.1%），故此处同时给出 6 位小数的精确值。
    f.write(f"  (exact) FGW cost={fgw_cost:.6f} entropy={fgw_ent:.6f}\n")
    f.write(f"  (exact) LROT cost={lrot_cost:.6f} entropy={lrot_ent:.6f}\n\n")
    f.write(f"  Cost improvement: {(1 - lrot_cost/fgw_cost)*100:.2f}%\n")
    f.write(f"  Entropy improvement: {(1 - lrot_ent/fgw_ent)*100:.2f}%\n")
    f.write(f"  LROT loss iterations: {len(loss_hist)}\n")
    f.write(f"  Top LR pairs: {[(l, r) for l, r, w in top_pairs[:6]]}\n")
print(f"  [SAVE] 结果已保存: {out_dir / 'lrot_real_visium_results.txt'}")

coords_aligned = compute_alignment(P, slice_A['coords'], slice_B['coords'])

import fig_align  # 同目录共享工具(增强 Before/After: 高对比+放大插图+位移量化)
fig, axes = plt.subplots(1, 3, figsize=(17, 5.2))

d_b, d_a = fig_align.draw_before_after(axes, 0, 1,
                                       slice_A['coords'], slice_B['coords'], P,
                                       A_name='Slice A', B_name='Slice B',
                                       title_before='Before Alignment\n(raw, misaligned)',
                                       title_after='After LROT Alignment',
                                       sA=6, sB=6)

axes[2].plot(loss_hist, 'b-', linewidth=1.5)
axes[2].set_xlabel('Iteration', fontsize=10)
axes[2].set_ylabel('Loss', fontsize=10)
axes[2].set_title('Loss Curve ({} iter)'.format(len(loss_hist)), fontsize=11)
axes[2].grid(True, alpha=0.3)

plt.tight_layout()
save_path = out_dir / 'lrot_real_visium.png'
fig.savefig(str(save_path), dpi=300, bbox_inches='tight')
plt.close()
print("  Saved:", str(save_path))
# Preserve computation precision; do not overwrite the archival float16 NPZ.
np.savez_compressed(out_dir / 'mouse_transport.npz', P=P, P_fgw=P_fgw,
    A_coords=coords, B_coords=coords_B, A_expr=slice_A['expr'], B_expr=slice_B['expr'],
    gene_names=np.asarray(gene_names, dtype=str), raw_indices=idx,
    barcodes=np.asarray(adata.obs_names[idx], dtype=str), loss_hist=np.asarray(loss_hist))

def plan_metrics(plan):
    q = plan / np.maximum(plan.sum(1, keepdims=True), 1e-300)
    projected = (plan / np.maximum(plan.sum(0, keepdims=True), 1e-300)).T @ coords
    return {'row_entropy': float(-(q * np.log(np.maximum(q, 1e-300))).sum(1).mean()),
        'expression_cost': float(compute_ot_cost(plan, slice_A['expr'], slice_B['expr'])),
        'median_projection_displacement': float(np.median(np.linalg.norm(projected - coords, axis=1))),
        'same_index_argmax_fraction': float((plan.argmax(1) == np.arange(N_SAMPLE)).mean()),
        'total_mass': float(plan.sum()), 'finite': bool(np.isfinite(plan).all()),
        'minimum': float(plan.min()),
        'row_uniform_max_error': float(np.max(np.abs(plan.sum(1) - 1/N_SAMPLE))),
        'column_uniform_max_error': float(np.max(np.abs(plan.sum(0) - 1/N_SAMPLE)))}
metadata = {'input_sha256': input_hash, 'solver': ('paper_protocol source-bound revalidated solver' if args.solver == 'revalidated' else 'root lrot_core.fgw_lr_solver legacy protocol'),
    'solver_file_sha256': file_sha256(Path(selected_core.__file__)),
    'runner_sha256': file_sha256(Path(__file__)),
    'seed': 42, 'n_sample': N_SAMPLE, 'raw_spots': n_spots, 'genes': len(gene_names),
    'coordinates': 'generated staggered grid plus Gaussian noise; not measured tissue coordinates',
    'rotation_degrees': 12.0, 'scale': 0.97,
    'expression_perturbation': 'Gamma(shape=.5,scale=.5)*.2 additive; per-row Uniform(.9,1.1) multiplicative',
    'solver_parameters': {'alpha': .5, 'beta': .3, 'gamma': [.1,0.0], 'reg_gw': .01, 'reg_ot': .01, 'max_iter': 200, 'tol': 1e-5},
    'input_precision': ('float16 archival expression promoted to float32; float64 geometry' if args.solver == 'revalidated' else 'float32 expression; float64 geometry'), 'output_precision': 'full computation precision; archival expression was float16',
    'region_accuracy': None, 'region_accuracy_status': 'not applicable: no region annotations; zero labels are placeholders',
    'matching_scope': 'known same-index correspondence for a synthetic perturbation of measured expression; not biological ground truth',
    'LROT': plan_metrics(P), 'FGW': plan_metrics(P_fgw), 'loss_iterations': len(loss_hist),
    'environment': {'python': platform.python_version(), 'numpy': np.__version__, 'scanpy': __import__('importlib.metadata', fromlist=['version']).version('scanpy')},
    'artifacts': {p.name: file_sha256(p) for p in out_dir.iterdir() if p.is_file()}}
(out_dir / 'run_evidence.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')
print('Evidence saved:', out_dir / 'run_evidence.json')
print('Done!')
