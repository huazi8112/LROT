# -*- coding: utf-8 -*-
"""verify_env.py — 一键环境自检 + 最小复现验证。

用途：解压后**先跑这一条**，即可确认环境是否齐备、能否复现图件。

    python -X utf8 tools/verify_env.py

检查项
------
1. 解释器与关键包版本（与补充表 S14 / requirements.txt 对比）
2. 仓库自定位是否指向本包（而非包外的其他副本）
3. 当前正式清单指定的批准资产与绘图输入是否存在且哈希匹配
4. 可选：分别执行正文保存结果绘图与补图资产导出，记录图像哈希；不认证实验复现

退出码 0 = 全部通过；非 0 = 有缺项（会列出缺什么、怎么补）。
"""
import hashlib
import importlib.metadata
import importlib.util
import io
import os
import subprocess
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import lrot_paths as P                                                        # noqa: E402

FAIL = []

# ---------- 1. 解释器与依赖 ----------
print('=' * 78)
print('1. 解释器与依赖')
print('=' * 78)
print('  解释器      %s' % sys.executable)
print('  Python      %s' % sys.version.split()[0])
if sys.version_info[:2] != (3, 10):
    print('  ⚠ 本包在本机实测于 Python 3.10.11；其他版本可能因 matplotlib/scanpy '
          'API 变化导致个别脚本失败（尤其不要用 3.12）。')

NEED = {
    'numpy': '2.2.6', 'scipy': '1.15.3', 'matplotlib': '3.10.7', 'pandas': '2.3.3',
    'sklearn': '1.7.2', 'PIL': '12.0.0', 'h5py': '3.16.0', 'ot': '0.9.6.post1',
    'anndata': '0.11.4', 'scanpy': '1.11.5', 'networkx': '3.4.2', 'statsmodels': '0.14.6',
    'igraph': '1.0.0', 'leidenalg': '0.12.0',
}
DIST = {'sklearn': 'scikit-learn', 'PIL': 'pillow', 'ot': 'POT'}
miss, drift = [], []
for mod, want in NEED.items():
    if importlib.util.find_spec(mod) is None:
        miss.append(mod)
        print('  ✘ %-14s 未安装（期望 %s）' % (mod, want))
        continue
    try:
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            __import__(mod)
            got = importlib.metadata.version(DIST.get(mod, mod))
    except Exception as e:                                     # noqa: BLE001
        print('  ✘ %-14s 导入失败: %s' % (mod, e))
        miss.append(mod)
        continue
    flag = '✔' if got == want else '⚠'
    if got != want:
        drift.append((mod, want, got))
    print('  %s %-14s %-12s（期望 %s）' % (flag, mod, got, want))
if miss:
    FAIL.append('缺包: %s —— 请先 `pip install -r requirements.txt`' % ', '.join(miss))
if drift:
    print('  ⚠ 版本与实测环境不同：%s。多数情况可运行，若出图不一致请按本行装回原版本。'
          % ', '.join('%s %s→%s' % (m, w, g) for m, w, g in drift))
try:
    louvain_version = importlib.metadata.version('louvain')
except importlib.metadata.PackageNotFoundError:
    pass
else:
    print('  ⚠ 检测到可选包 louvain %s；仓库不使用它，且该包要求 igraph<0.12，'
          '与本项目 igraph 1.0.0 不兼容。不要安装 scanpy[louvain] extra。' % louvain_version)

# ---------- 2. 自定位 ----------
print('')
print('=' * 78)
print('2. 仓库自定位')
print('=' * 78)
print('  ROOT        %s' % P.ROOT)
print('  OUT         %s' % P.OUT)
if not (P.ROOT / 'lrot_core.py').exists():
    FAIL.append('ROOT 下没有 lrot_core.py，包可能不完整')
if not (P.ROOT / 'lrot_output').is_dir():
    FAIL.append('ROOT 下没有 lrot_output/，结果目录缺失')
print('  ✔ ROOT 内已含 lrot_core.py 与 lrot_output/，解压路径无关')

# ---------- 3. Tier 1 输入 ----------
print('')
print('=' * 78)
print('3. 当前清单中的批准资产与保存结果绘图输入')
print('=' * 78)
import json
import csv
from pathlib import Path
main_manifest = json.loads((P.ROOT / 'figures/current_release/figure_manifest.json').read_text(encoding='utf-8'))
supp_manifest = json.loads((P.ROOT / 'figures/supplement/figure_manifest.json').read_text(encoding='utf-8'))
required = []
for row in main_manifest['figures']:
    if row['id'].startswith('F') and row['render_status'] == 'ready':
        required.extend(row.get('inputs', []))
with (P.ROOT / 'figures/current_release/FINAL_PAPER_IMAGE_INDEX.csv').open(encoding='utf-8-sig', newline='') as handle:
    for row in csv.DictReader(handle):
        required.append({'path': row['selected_repository_path'], 'sha256': row['selected_file_sha256']})
for row in supp_manifest['figures']:
    required.append(row['asset'])
    if row.get('record_render_supported'):
        required.extend([row['template'], row['record_config'], row['record_renderer'], *row['records']])
checked = set()
for spec in required:
    key = (spec['path'], spec['sha256'])
    if key in checked:
        continue
    checked.add(key)
    path = (P.ROOT / spec['path']).resolve()
    if P.ROOT.resolve() not in path.parents or not path.is_file():
        FAIL.append('Missing or out-of-package pinned drawing input: ' + spec['path'])
    elif hashlib.sha256(path.read_bytes()).hexdigest() != spec['sha256']:
        FAIL.append('Drawing input hash mismatch: ' + spec['path'])
print('  Checked %d unique pinned assets/drawing inputs from current manifests.' % len(checked))
print('  F2 adopted four-curve reconstruction has a separate execution entrypoint; asset export alone does not verify it.')

# ---------- 4. 最小复现 ----------
print('')
print('=' * 78)
print('4. 可选绘图检查（生成状态与图像身份分别记录）')
print('=' * 78)
import argparse
from pathlib import Path
_parser = argparse.ArgumentParser()
_parser.add_argument('--render-output-dir', type=Path)
_args = _parser.parse_args()
if _args.render_output_dir is None:
    print('  SKIP figure rendering; supply --render-output-dir with a new project-external directory.')
    print('  Environment checks do not establish current manuscript or full experiment reproduction.')
else:
    _out = _args.render_output_dir.resolve()
    if _out.exists() or _out == P.ROOT.resolve() or P.ROOT.resolve() in _out.parents:
        _parser.error('render output must be new and outside the code project')
    import json
    _out.mkdir(parents=True, exist_ok=False)
    checks = []
    for figure, mode, subdir in [('F1', 'records', 'main'), ('S23', 'asset', 'supplement')]:
        result = subprocess.run(
            [sys.executable, '-X', 'utf8', str(P.ROOT / 'figures/render_current.py'),
             '--figure', figure, '--mode', mode, '--output-dir', str(_out / subdir)],
            cwd=str(P.ROOT), capture_output=True, text=True, encoding='utf-8')
        (_out / (figure + '.log')).write_text(result.stdout + result.stderr, encoding='utf-8')
        check = {'figure': figure, 'mode': mode, 'exit_code': result.returncode,
                 'generated': result.returncode == 0}
        if result.returncode:
            FAIL.append(figure + ': figure check failed; see isolated output log')
        else:
            report_name = 'render_run.json' if figure.startswith('F') else 'validation.json'
            report = json.loads((_out / subdir / report_name).read_text(encoding='utf-8'))
            output = report['figures' if figure.startswith('F') else 'outputs'][0]
            check['sha256'] = output['sha256']
            if figure.startswith('F'):
                manifest = json.loads((P.ROOT / 'figures/current_release/figure_manifest.json').read_text(encoding='utf-8'))
                expected = next(item['reference_sha256'] for item in manifest['figures'] if item['id'] == figure)
            else:
                expected = output['embedded_sha256']
            check['approved_sha256'] = expected
            check['byte_equal_approved'] = output['sha256'] == expected
            if mode == 'asset' and not check['byte_equal_approved']:
                FAIL.append(figure + ': exported asset differs from approved image')
            print('  %s %s: generated; byte_equal_approved=%s' % (figure, mode, check['byte_equal_approved']))
        checks.append(check)
    (_out / 'verification.json').write_text(json.dumps(
        {'checks': checks, 'solver_executed': False, 'full_experiment_reproduction_certified': False,
         'scope': 'Rendering success and image identity are separate; a schematic redraw can differ from the approved asset.'},
        ensure_ascii=False, indent=2), encoding='utf-8')

print('')
print('=' * 78)
if FAIL:
    print('结论：%d 项需要处理' % len(FAIL))
    for x in FAIL:
        print('  - %s' % x)
    print('=' * 78)
    sys.exit(1)
print('结论：环境与包自检全部通过，可以开始复现（见 RUN_REPRO.md）。')
print('=' * 78)
