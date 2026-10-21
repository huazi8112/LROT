"""Render the approved F1 canvas with its revised LR-prior label only.

The approved image supplies retained pixels. The declared label band is
replaced completely; this is a local label render, not a full schematic redraw.
No experimental values or equations are redrawn.
"""
from pathlib import Path
from io import BytesIO
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
LABEL = 'Species-specific manually specified LR prior'
BAND = (140, 770, 1600, 850)

def render(output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    canvas = Image.open(ROOT / 'figures/current_release/assets/F1.png').convert('RGB')
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
    fig = plt.figure(figsize=(4095 / 650, 5077 / 650), dpi=650, facecolor='white')
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 124)
    ax.axis('off')
    ax.text(3.6, 104.2, LABEL, fontsize=6.4, fontweight='bold', color='#20272E', ha='left', va='center')
    buf = BytesIO()
    fig.savefig(buf, dpi=650, facecolor='white')
    plt.close(fig)
    buf.seek(0)
    layer = Image.open(buf).convert('RGB')
    canvas.paste(layer.crop(BAND), BAND[:2])
    target = output_dir / 'lrot_method_flow.png'
    canvas.save(target, dpi=(650, 650))
    return target

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    print(render(args.output_dir))
