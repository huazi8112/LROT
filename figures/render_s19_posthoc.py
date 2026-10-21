"""Rebuild S19B from pinned saved P and post-hoc bound ordered inputs. No solver."""
from pathlib import Path
import argparse,hashlib,json,io,sys

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return h.hexdigest()

def render_one(figure,bundle,out):
    if figure!='S19':raise ValueError('S19 renderer only')
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from PIL import Image,ImageChops
    def ah(a):return hashlib.sha256(str(a.shape).encode()+str(a.dtype).encode()+np.ascontiguousarray(a).tobytes()).hexdigest()
    folder=Path(bundle)/'s19_posthoc';out=Path(out);out.mkdir(parents=True,exist_ok=True)
    cfg=json.loads((folder/'render_config.json').read_text(encoding='utf-8'))
    bp=folder/cfg['binding']
    if sha(bp)!=cfg['binding_sha256']:raise ValueError('Binding file changed')
    rec=json.loads(bp.read_text(encoding='utf-8'))
    template=folder/cfg['template'];inp=folder/rec['ordered_inputs']['file']
    if sha(template)!=cfg['template_sha256'] or sha(inp)!=rec['ordered_inputs']['sha256']:raise ValueError('Pinned template or ordered inputs changed')
    with np.load(inp,allow_pickle=False) as z:arr={k:z[k] for k in z.files}
    for c in rec['ordered_input_checks']:
        k=c['side']+('_coords' if c['field']=='coords' else '_region_labels')
        if ah(arr[k])!=c['expected']:raise ValueError('Ordered input hash mismatch: '+k)
    chunks=[]
    for part in rec['plan']['parts']:
        p=folder/part['file']
        if p.stat().st_size!=part['size_bytes'] or sha(p)!=part['sha256']:raise ValueError('Saved plan chunk changed: '+part['file'])
        chunks.append(p.read_bytes())
    raw=b''.join(chunks);del chunks
    if hashlib.sha256(raw).hexdigest()!=rec['plan']['original_file_sha256']:raise ValueError('Reassembled NPZ does not match original saved plan file')
    with np.load(io.BytesIO(raw),allow_pickle=False) as z:P=z['P']
    del raw
    if list(P.shape)!=rec['plan']['shape'] or ah(P)!=rec['plan']['P_array_sha256']:raise ValueError('Saved plan identity mismatch')
    if P.shape!=(len(arr['A_coords']),len(arr['B_coords'])):raise ValueError('Plan row/column counts do not match ordered inputs')
    col=P.sum(axis=0)
    if not np.isfinite(P).all() or (P<0).any() or (col<=0).any():raise ValueError('Invalid saved plan for mapping')
    mapped=(P/col[None,:]).T@arr['A_coords']
    if ah(mapped)!=rec['mapping']['mapped_B_array_sha256']:raise ValueError('Mapped coordinate identity differs: inspect runtime; never reorder inputs')
    if not np.allclose(mapped,arr['mapped_B'],rtol=1e-12,atol=1e-10):raise ValueError('Saved mapping and independent recomputation differ')
    if len(arr['B_barcodes'])!=len(mapped) or len(arr['B_region_labels'])!=len(mapped):raise ValueError('Target barcode/label mismatch')
    A=arr['A_coords'];B=arr['B_coords'];low=np.minimum(np.minimum(A.min(axis=0),B.min(axis=0)),mapped.min(axis=0));high=np.maximum(np.maximum(A.max(axis=0),B.max(axis=0)),mapped.max(axis=0));center=(low+high)/2;half=max(high-low)/2*cfg['style']['spatial_margin_factor']
    box=cfg['B_interior_box'];width,height=box[2]-box[0],box[3]-box[1];dpi=cfg['style']['dpi']
    # Agg truncates a floating-point canvas size; an epsilon avoids a one-pixel roundoff.
    fig=Figure(figsize=((width+0.01)/dpi,(height+0.01)/dpi),dpi=dpi,facecolor='white');canvas=FigureCanvasAgg(fig);ax=fig.add_axes([0,0,1,1]);ax.set_xlim(center[0]-half,center[0]+half);ax.set_ylim(center[1]-half,center[1]+half);ax.set_aspect('auto');ax.set_axis_off()
    ax.scatter(A[:,0],A[:,1],c=cfg['style']['reference_gray'],s=cfg['style']['reference_size_pt2'],alpha=cfg['style']['reference_alpha'],linewidths=0,zorder=1)
    for i,color in enumerate(cfg['style']['colors']):
        mark=arr['B_region_labels']==i
        ax.scatter(mapped[mark,0],mapped[mark,1],s=cfg['style']['mapped_size_pt2'],c=[color],alpha=cfg['style']['mapped_alpha'],linewidths=0,zorder=2)
    canvas.draw();patch=Image.fromarray(np.asarray(canvas.buffer_rgba()).copy()).convert('RGB')
    if patch.size!=(width,height):raise ValueError('Unexpected raster dimensions')
    before=Image.open(template).convert('RGB');after=before.copy();after.paste(patch,(box[0],box[1]))
    delta=ImageChops.difference(before,after);delta.paste((0,0,0),box)
    if delta.getbbox() is not None:raise ValueError('Pixels outside approved B rectangle changed')
    dest=out/'S19.png';after.save(dest,dpi=(300,300))
    if cfg.get('expected_output_sha256') and sha(dest)!=cfg['expected_output_sha256']:raise ValueError('Output raster differs from pinned review figure; inspect plotting runtime')
    proof={'figure':'S19','output':str(dest),'sha256':sha(dest),'template_sha256':sha(template),'B_interior_box':box,'A_and_all_outside_B_pixels_unchanged':True,'coordinate_window':{'x':[float(center[0]-half),float(center[0]+half)],'y':[float(center[1]-half),float(center[1]+half)]},'P_array_sha256':ah(P),'mapped_B_array_sha256':ah(mapped),'binding_sha256':sha(bp),'ordered_coordinate_and_label_hashes_match':True,'posthoc_reconstruction':True,'old_historical_image_binding_restored':False,'solver_or_training_run':False,'environment':{'python':sys.version,'numpy':np.__version__,'matplotlib':matplotlib.__version__},'renderer_sha256':sha(Path(__file__))}
    (out/'S19_render_validation.json').write_text(json.dumps(proof,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    return proof

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--bundle',type=Path,default=Path(__file__).resolve().parent/'supplement');p.add_argument('--output-dir',required=True,type=Path);args=p.parse_args()
    root=Path(__file__).resolve().parents[1];out=args.output_dir.resolve()
    if out.exists() or out==root or root in out.parents:p.error('Use a NEW output directory outside the code repository')
    render_one('S19',args.bundle,out)
