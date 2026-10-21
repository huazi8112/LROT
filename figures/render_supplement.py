"""Export locked supplementary assets or apply approved local edits, without experiments.

The current manifest is authoritative for S1-S23. Pending whole figures fail
closed. Output must be a new directory outside the repository. Asset equality
does not certify experiment reproduction or trigger manuscript replacement.
"""
from pathlib import Path
import argparse,hashlib,json,shutil,sys,platform

ROOT=Path(__file__).resolve().parents[1]
BUNDLE=ROOT/'figures/supplement'

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def pinned(item):
    path=(ROOT/item['path']).resolve()
    if ROOT not in path.parents or not path.is_file() or sha(path)!=item['sha256']:
        raise ValueError('Pinned asset or result absent/changed: '+item['path'])
    return path

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--list',action='store_true')
    parser.add_argument('--figure',action='append',help='Repeat for selected S figures')
    parser.add_argument('--mode',choices=['asset','patch','panels','panel-patch','records'],default='asset')
    parser.add_argument('--panel',action='append',help='Repeat for selected approved panels')
    parser.add_argument('--output-dir',type=Path)
    parser.add_argument('--font-regular',type=Path)
    parser.add_argument('--font-bold',type=Path)
    args=parser.parse_args(argv)
    index=json.loads((BUNDLE/'figure_manifest.json').read_text(encoding='utf-8'))
    rows={r['id']:r for r in index['figures']}
    if args.list:
        for r in rows.values():
            print(r['id'],r['status'],r['treatment'],'pending='+','.join(r['pending_panels']))
        return 0
    ids=args.figure or [r['id'] for r in rows.values() if (r['status']=='approved_asset' if args.mode=='asset' else r.get('record_render_supported',False) if args.mode=='records' else r['patch_supported'] if args.mode=='patch' else any(p.get('patch_supported') for p in r['confirmed_panels']) if args.mode=='panel-patch' else r['confirmed_panels'])]
    if not ids or len(ids)!=len(set(ids)) or any(f not in rows for f in ids):
        parser.error('Choose unique known S figures; see --list.')
    if args.panel and (args.mode not in ['panels','panel-patch'] or len(args.panel)!=len(set(args.panel))):
        parser.error('--panel requires panel mode and unique panel identifiers')
    selected={}
    for f in ids:
        r=rows[f]
        if args.mode in ['panels','panel-patch']:
            selected[f]=[p for p in r['confirmed_panels'] if (not args.panel or p['panel'] in args.panel) and (args.mode!='panel-patch' or p.get('patch_supported'))]
            if not selected[f] or (args.panel and set(args.panel)!={p['panel'] for p in selected[f]}):
                parser.error(f+': requested panel is not independently approved for this operation')
            if not r['confirmed_panels']:parser.error(f+': no independently approved panel asset')
            for p in selected[f]:
                pinned(p['source'] if 'crop_box' in p else p)
                if args.mode=='panel-patch':
                    pinned(p['template'])
                    for record in p['records']:pinned(record)
        else:
            if r['status']!='approved_asset':parser.error(f+': whole figure pending; no old/new value selection or redraw permitted. Pending: '+','.join(r['pending_panels']))
            pinned(r['asset'])
            if args.mode=='records':
                if not r.get('record_render_supported'):parser.error(f+': no saved-record drawing implementation; use asset export')
                pinned(r['template']);pinned(r['record_config']);pinned(r['record_renderer'])
                for record in r['records']:pinned(record)
            if args.mode=='patch':
                if not r['patch_supported']:parser.error(f+': retained asset, not result regeneration; use --mode asset')
                pinned(r['template'])
                for p in r['records']:pinned(p)
    if args.mode in ['patch','panel-patch']:
        for f,name in [(args.font_regular,'arial.ttf'),(args.font_bold,'arialbd.ttf')]:
            if f is None or not f.is_file() or sha(f)!=index['font_sha256'][name]:
                parser.error('Local patch validation requires the exact pinned '+name+' via --font-regular/--font-bold. Fonts are not bundled.')
    if args.output_dir is None:parser.error('--output-dir is mandatory')
    out=args.output_dir.resolve()
    if out==ROOT or ROOT in out.parents or out.exists():
        parser.error('Output must be new and outside the code repository; never overwrite manuscript assets.')
    out.mkdir(parents=True,exist_ok=False)
    outputs=[]
    for f in ids:
        r=rows[f]
        if args.mode in ['panels','panel-patch']:
            for p in selected[f]:
                dest=out/(f+'_'+p['panel']+'.png');proofs=[]
                if args.mode=='panels':
                    if 'crop_box' in p:
                        from PIL import Image
                        Image.open(pinned(p['source'])).convert('RGB').crop(p['crop_box']).save(dest)
                        if sha(dest)!=p['sha256']:raise ValueError('Current asset crop hash mismatch')
                    else:shutil.copyfile(pinned(p),dest)
                else:
                    if f!='S2' or p['panel']!='D':raise ValueError('No approved local patch implementation for this panel')
                    from patch_s2d_labels import apply
                    proofs=[apply(pinned(p['template']),pinned(p['records'][0]),dest,args.font_regular)]
                    if sha(dest)!=p['sha256']:raise ValueError('Approved panel patch did not match locked asset')
                outputs.append({'figure':f,'panel':p['panel'],'sha256':sha(dest),'file':dest.name,'verification':p['meaning'],'byte_equal_locked_panel':sha(dest)==p['sha256'],'pixel_proofs':proofs})
        else:
            dest=out/(f+'.png');proofs=[]
            if args.mode=='asset':shutil.copyfile(pinned(r['asset']),dest)
            elif args.mode=='records':
                if f=='S19':
                    from render_s19_posthoc import render_one
                elif f=='S4':
                    from render_s4_saved import render_one as render_s4
                    render_one = lambda fid, bundle, out: render_s4(bundle, out)
                else:
                    from render_selected_panels import render_one
                proof=render_one(f,BUNDLE,out);proofs=[proof]
                if sha(dest)!=r['embedded_sha256']:raise ValueError(f+': saved-record output differs from locked DOCX asset; inspect runtime/fonts, never replace manuscript automatically')
            else:
                from apply_supplement_revisions import apply
                _,proofs=apply(int(f[1:]),BUNDLE,dest,args.font_regular,args.font_bold)
                if any(v is False for p in proofs for k,v in p.items() if k.endswith('unchanged') or k=='pixel_identical'):
                    raise ValueError(f+': protected original pixels changed')
            outputs.append({'figure':f,'sha256':sha(dest),'embedded_sha256':r['embedded_sha256'],'byte_equal_embedded':sha(dest)==r['embedded_sha256'],'file':dest.name,'verification':r['treatment'],'records':r['records'],'pixel_proofs':proofs})
        print(f,args.mode,flush=True)
    prohibited={'lrot_core','lrot_prototype','scanpy','torch','ot'} & set(sys.modules)
    if prohibited:raise RuntimeError('Unexpected experiment dependency imported: '+','.join(sorted(prohibited)))
    report={'mode':args.mode,'python':sys.version,'executable':sys.executable,'platform':platform.platform(),'manifest_sha256':sha(BUNDLE/'figure_manifest.json'),'entry_sha256':sha(Path(__file__)),'patch_script_sha256':sha(ROOT/'figures/apply_supplement_revisions.py'),'outputs':outputs,'experiment_or_solver_run':False,'manuscript_overwritten':False,'scope':'Asset/approved local-edit validation only; no end-to-end experiment certification.'}
    (out/'validation.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    return 0

if __name__=='__main__':
    raise SystemExit(main())
