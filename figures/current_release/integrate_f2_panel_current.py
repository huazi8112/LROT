from __future__ import annotations
import hashlib, json, argparse
from pathlib import Path
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'figures' / 'current_release'
PROJECT = ROOT
BASE = PROJECT / 'figures' / 'current_release' / 'assets' / 'F2.png'
HIST = PROJECT / 'figures' / 'current_release' / 'inputs' / 'F2_four_gamma_histories_current_code.npz'
parser = argparse.ArgumentParser(description='Render adopted F2 panel F from frozen histories; retain the other raster panels.')
parser.add_argument('--input-npz', type=Path)
parser.add_argument('--output-png', type=Path, required=True)
args = parser.parse_args()
if args.input_npz is not None:
    HIST = args.input_npz.resolve()
if args.output_png is not None:
    output = args.output_png.resolve()
    if output.exists() or output == ROOT or ROOT in output.parents:
        parser.error('Output must be a new PNG outside the code package')
    OUT = output.parent
(OUT / 'evidence').mkdir(parents=True, exist_ok=True)

SOURCE = PROJECT / 'experiments' / 'lrot_prototype.py'
CROP = (3540, 1500, 5220, 2973)  # left, top, right, bottom; exact native-pixel panel patch
AXIS_BOUNDS_NATIVE = (3699, 1603, 5081, 2806)  # detected black plot-frame boundary

sha = lambda b: hashlib.sha256(b).hexdigest()
base = Image.open(BASE).convert('RGB')
if base.size != (5220, 2973):
    raise RuntimeError(f'Unexpected baseline raster size: {base.size}')
with np.load(HIST, allow_pickle=False) as z:
    main = np.asarray(z['saved_main_gamma_0_10'], dtype=np.float64)
    curves = {g: np.asarray(z[f'gamma_{g:.2f}'], dtype=np.float64) for g in (0.0, 0.05, 0.1, 0.2)}
if main.size != 123 or any(v.size != 80 for v in curves.values()):
    raise RuntimeError('Saved histories have unexpected lengths')

left, top, right, bottom = CROP
width, height = right-left, bottom-top
x0, y0, x1, y1 = AXIS_BOUNDS_NATIVE
fig = plt.figure(figsize=(width/300, height/300), dpi=300, facecolor='white')
ax = fig.add_axes([(x0-left)/width, (bottom-y1)/height, (x1-x0)/width, (y1-y0)/height])
ax.plot(main, 'b-', linewidth=1.5, label=r'LROT ($\gamma=0.1$)')
colors = ['gray', 'orange', 'green', 'red']
for i, gamma in enumerate((0.0, 0.05, 0.1, 0.2)):
    ax.plot(curves[gamma], color=colors[i], linewidth=1.0, alpha=0.7,
            label=rf'$\gamma={gamma:.2f}$')
ax.set_title('F. Loss Curves (Ablation)', fontsize=11)
ax.set_xlabel('Iteration')
ax.set_ylabel('Loss')
ax.legend(fontsize=8)
panel_path = OUT / 'evidence' / 'F2_panel_F_current_code_native_layout.png'
fig.savefig(panel_path, dpi=300, facecolor='white', edgecolor='white', bbox_inches=None, pad_inches=0)
plt.close(fig)
panel = Image.open(panel_path).convert('RGB')
if panel.size != (width, height):
    raise RuntimeError(f'Panel render has unexpected size {panel.size}, expected {(width,height)}')
base.paste(panel, (left, top))
composite_path = args.output_png.resolve()
base.save(composite_path, format='PNG', dpi=(550.0116, 550.0116), optimize=False, compress_level=9)

# Verify all pixels outside the designated F-panel patch are unchanged.
new = np.asarray(Image.open(composite_path).convert('RGB'))
old = np.asarray(Image.open(BASE).convert('RGB'))
mask = np.ones(old.shape[:2], dtype=bool)
mask[top:bottom, left:right] = False
outside_changed = int(np.any(new != old, axis=2)[mask].sum())
if outside_changed:
    raise RuntimeError(f'{outside_changed} pixels changed outside panel F patch')
record = {
  'status':'adopted_panel_F_saved_result_redraw',
  'operation':'render_only_from_previously_saved_current_code_histories; no solver called in this integration step',
  'baseline_image':{'path':str(BASE),'sha256':sha(BASE.read_bytes()),'size_px':list(old.shape[1::-1])},
  'source_code':{'path':str(SOURCE),'sha256':sha(SOURCE.read_bytes())},
  'saved_histories':{'path':str(HIST),'sha256':sha(HIST.read_bytes()),'main_length':int(main.size),'ablation_lengths':{str(k):int(v.size) for k,v in curves.items()}},
  'panel_patch_native_px':{'left':left,'top':top,'right':right,'bottom':bottom,'width':width,'height':height},
  'axis_frame_native_px':{'left':x0,'top':y0,'right':x1,'bottom':y1},
  'layout':'panel F redrawn at target native pixel dimensions; no resizing or non-uniform stretching; all pixels outside F patch are unchanged',
  'gamma_01_prefix_matches_saved_main':bool(np.array_equal(curves[0.1],main[:80])),
  'composite_output':{'path':str(composite_path),'sha256':sha(composite_path.read_bytes()),'size_px':list(new.shape[1::-1])},
  'panel_output':{'path':str(panel_path),'sha256':sha(panel_path.read_bytes()),'size_px':list(panel.size)},
  'pixels_changed_outside_F_patch':outside_changed,
  'historical_status':'The four original per-run ablation histories were not retained; this candidate uses newly reconstructed current-code histories and does not claim historical run recovery.'
}
(OUT/'evidence'/'F2_composition_record.json').write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps(record,ensure_ascii=False,indent=2))


