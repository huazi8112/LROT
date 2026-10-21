"""Compatibility entry for current S2/S4/S11/S15 asset export.

The historical stacked-component builder is retired: its components predate
the adopted figures. --verify checks current approved asset hashes only;
export copies current assets and does not reconstruct numerical results.
Use render_supplement.py --mode records for supported saved-record drawing.
"""
from pathlib import Path
import argparse, json, hashlib
from render_supplement import main as export
ROOT=Path(__file__).resolve().parents[1]
IDS=('S2','S4','S11','S15')
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify',action='store_true')
    parser.add_argument('--output-dir',type=Path)
    args=parser.parse_args()
    if args.verify:
        if args.output_dir:parser.error('--verify does not export files')
        manifest=json.loads((ROOT/'figures/supplement/figure_manifest.json').read_text(encoding='utf-8'))
        rows={r['id']:r for r in manifest['figures']}
        for fid in IDS:
            spec=rows[fid]['asset'];path=(ROOT/spec['path']).resolve()
            if ROOT not in path.parents or not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=spec['sha256']:
                raise ValueError(fid+': current approved asset missing or changed')
            print(fid, 'current asset hash verified; no composite regeneration')
        return 0
    if args.output_dir is None:parser.error('--output-dir is mandatory')
    argv=['--mode','asset','--output-dir',str(args.output_dir)]
    for fid in IDS:argv+=['--figure',fid]
    return export(argv)
if __name__=='__main__':raise SystemExit(main())
