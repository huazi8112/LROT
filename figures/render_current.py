"""Render current paper figures from pinned saved results, never from a solver.

The output directory is mandatory, must be new, and must be outside the code
project. Missing inputs fail closed; historical charts are never a fallback.
"""
from pathlib import Path
import argparse, hashlib, json, os, sys, subprocess, csv, shutil

ROOT = Path(__file__).resolve().parents[1]
BUNDLE = ROOT / 'figures' / 'current_release'
sys.path.insert(0, str(ROOT / 'figures'))

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--figure', action='append', help='F1, F3, S6 etc.; repeat to select')
    parser.add_argument('--list', action='store_true')
    parser.add_argument('--mode', choices=['records', 'asset'],
                        help='records: saved-result drawing; asset: exact approved-image export, no numerical reconstruction')
    parser.add_argument('--output-dir', type=Path)
    args = parser.parse_args()
    # Supplementary assets follow the dated panel decisions, separately from main figures.
    if args.figure and any(i.startswith('S') for i in args.figure):
        if any(not i.startswith('S') for i in args.figure):
            parser.error('Render main and supplementary figures separately.')
        from render_supplement import main as render_supplement
        selected=[]
        for i in args.figure:selected.extend(['--figure',i])
        if args.output_dir is not None:selected.extend(['--output-dir',str(args.output_dir)])
        if args.list:selected.append('--list')
        if args.mode is not None:selected.extend(['--mode', args.mode])
        return render_supplement(selected)
    index = json.loads((BUNDLE / 'figure_manifest.json').read_text(encoding='utf8'))
    rows = {r['id']: r for r in index['figures']}
    if args.list:
        for r in rows.values():
            if r['id'].startswith('S'):continue
            print(r['id'], r['render_status'], r.get('blocked_reason') or r.get('provenance_warning', ''))
        from render_supplement import main as render_supplement
        render_supplement(['--list'])
        return
    if args.mode == 'asset':
        ids = args.figure or [i for i in rows if i.startswith('F')]
        if not ids or len(ids) != len(set(ids)) or any(i not in rows or not i.startswith('F') for i in ids):
            parser.error('Select unique known main figure IDs for asset export.')
        with (BUNDLE / 'FINAL_PAPER_IMAGE_INDEX.csv').open(encoding='utf-8-sig', newline='') as handle:
            assets = {item['figure']: item for item in csv.DictReader(handle) if item['document'] == 'main'}
        pinned = {}
        for i in ids:
            item = assets[i]
            source = (ROOT / item['selected_repository_path']).resolve()
            if ROOT not in source.parents or not source.is_file():
                parser.error(i + ': approved asset absent or outside repository')
            expected = rows[i]['reference_sha256']
            if sha(source) != expected or item['selected_file_sha256'] != expected or item['embedded_sha256'] != expected:
                parser.error(i + ': approved asset hash mismatch')
            pinned[i] = source
        if args.output_dir is None:
            parser.error('--output-dir is mandatory')
        out = args.output_dir.resolve()
        if out == ROOT or ROOT in out.parents or out.exists():
            parser.error('Output must be a new directory outside the code project.')
        out.mkdir(parents=True, exist_ok=False)
        outputs = []
        for i in ids:
            dest = out / (i + '.png')
            shutil.copyfile(pinned[i], dest)
            outputs.append({'figure': i, 'file': str(dest), 'sha256': sha(dest),
                            'source': str(pinned[i].relative_to(ROOT)),
                            'render_method': 'approved_asset_export', 'byte_equal_approved': sha(dest) == rows[i]['reference_sha256']})
            print(i, dest, flush=True)
        report = {'mode': 'asset', 'figures': outputs, 'solver_executed': False,
                  'verification_scope': 'Approved image bytes only; no saved-array reconstruction or experimental reproduction.'}
        (out / 'render_run.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
        return
    ids = args.figure or [r['id'] for r in rows.values() if r['id'].startswith('F') and r['render_status'] == 'ready']
    if not ids or len(ids) != len(set(ids)) or any(i not in rows for i in ids):
        parser.error('Select unique, known figure IDs; use --list.')
    for i in ids:
        r = rows[i]
        if r['render_status'] != 'ready':
            parser.error(f'{i}: {r.get("blocked_reason", "saved drawing inputs unavailable")}')
        for inp in r.get('inputs', []):
            p = ROOT / inp['path']
            if not p.is_file() or sha(p) != inp['sha256']:
                parser.error(f'{i}: pinned input absent or changed: {p}')
    if args.output_dir is None:
        parser.error('--output-dir is mandatory; paper assets are never overwritten.')
    out = args.output_dir.resolve()
    if out == ROOT or ROOT in out.parents or out.exists():
        parser.error('Output must be a new directory outside the code project.')
    os.environ.update(MPLBACKEND='Agg', MPLCONFIGDIR=str(out / 'mpl_cache'),
                      OPENBLAS_NUM_THREADS='2', OMP_NUM_THREADS='2')
    out.mkdir(parents=True, exist_ok=False)
    outputs=[]
    for i in ids:
        r=rows[i]
        renderer_sha256=None
        render_method='pinned_saved_payload'
        if i == 'F1':
            render_method='approved_canvas_label_render'
            run=subprocess.run([sys.executable, '-X', 'utf8', str(ROOT/'figures/render_f1_approved.py'),
                                '--output-dir', str(out/'F1')],capture_output=True,text=True,encoding='utf8')
            (out/'F1.log').write_text(run.stdout+run.stderr,encoding='utf8')
            if run.returncode: raise RuntimeError('F1 rendering failed; see F1.log')
            png=out/'F1'/'lrot_method_flow.png'
        elif r.get('saved_renderer'):
            spec=r['saved_renderer']
            script=ROOT/spec['script_path']
            png=out/i/spec['output_name']
            cmd=[sys.executable,'-X','utf8',str(script),
                 '--input-npz',str(ROOT/spec['input_npz']),
                 '--output-png',str(png)]
            if spec.get('results_txt'):
                cmd.extend(['--results-txt',str(ROOT/spec['results_txt'])])
            run=subprocess.run(cmd,capture_output=True,text=True,encoding='utf8')
            (out/f'{i}.log').write_text(run.stdout+run.stderr,encoding='utf8')
            if run.returncode: raise RuntimeError(f'{i} saved-result rendering failed; see {i}.log')
            if not png.is_file(): raise RuntimeError(f'{i} renderer returned without creating {png}')
            renderer_sha256=sha(script)
            render_method='standalone_saved_result_script'
        else:
            from render_current_panels import render
            payload=json.loads((ROOT/r['inputs'][0]['path']).read_text(encoding='utf8'))
            png=render(i,payload,out)
        outputs.append({'figure':i,'file':str(png),'sha256':sha(png),
                        'inputs':r.get('inputs',[]),'source_sha256':sha(Path(__file__)),
                        'renderer_sha256':renderer_sha256,'render_method':render_method,
                        'interpretation':r['interpretation']})
        print(i,png,flush=True)
    import platform, importlib.metadata
    env={}
    for p in ['numpy','matplotlib','pillow']:
        try: env[p]=importlib.metadata.version(p)
        except importlib.metadata.PackageNotFoundError: env[p]=None
    report={'python':sys.version,'executable':sys.executable,'platform':platform.platform(),
            'packages':env,'figures':outputs,'solver_executed':False,
            'verification_scope':'Saved-results rendering only; no new experiment or ranking.'}
    (out/'render_run.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf8')

if __name__ == '__main__':
    main()
