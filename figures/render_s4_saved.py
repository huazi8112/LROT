"""Rebuild S4A from retained four-decimal summaries; preserve B/C pixels.

No solver, training or raw tissue input is used. Output is outside the package.
"""
from pathlib import Path
import argparse,hashlib,json
import numpy as np
from PIL import Image,ImageChops
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def render_one(bundle,output):
    bundle=Path(bundle);output=Path(output);output.mkdir(parents=True,exist_ok=True)
    manifest=json.loads((bundle/'figure_manifest.json').read_text(encoding='utf-8'))
    row=next(x for x in manifest['figures'] if x['id']=='S4')
    root=bundle.parents[1]
    template=root/row['template']['path']
    record=root/row['panel_A_summary']['path']
    for p,expected in [(template,row['template']['sha256']),(record,row['panel_A_summary']['sha256'])]:
        if sha(p)!=expected:raise ValueError('S4 pinned source changed: '+p.name)
    data=json.loads(record.read_text(encoding='utf-8'));v=np.array(data['rows'],dtype=float)
    assert v.shape==(6,6)
    base=Image.open(template).convert('RGB');assert base.size==(3561,3200)
    with plt.rc_context({'font.family':'DejaVu Sans','font.size':13}):
        fig,ax=plt.subplots(figsize=(3561/300,2016/300),dpi=300)
        fig.subplots_adjust(left=.12,right=.85,bottom=.15,top=.89)
        styles=[('o-','#2166AC','Neurotrophin/Growth'),('s--','#D6604D','Ephrin/Semaphorin'),('^-','#4DAF4A','Chemokine/Cytokine'),('d--','#984EA3','Non-LR control')]
        for col,(sty,color,label) in enumerate(styles,2):ax.plot(v[:,0],v[:,col],sty,color=color,lw=2,ms=6,label=label)
        ax.set_ylim(0,1);ax.set_xlabel('LR guidance weight γ',fontsize=15);ax.set_ylabel('Pathway activity correlation',fontsize=14);ax.grid(alpha=.2);ax.set_title('(A) Simulation: pathway agreement and entropy vs γ',fontsize=14)
        a2=ax.twinx();a2.plot(v[:,0],v[:,1],'*-',color='#b2182b',lw=2.4,ms=9,label='Mean row entropy');a2.set_ylabel('Mean row entropy (nats)',color='#b2182b',fontsize=14);a2.tick_params(axis='y',labelcolor='#b2182b');a2.set_ylim(3.4,4.05)
        h,l=ax.get_legend_handles_labels();h2,l2=a2.get_legend_handles_labels();ax.legend(h+h2,l+l2,loc='upper right',fontsize=11,framealpha=.95)
        a2.text(.20,3.98,'Entropy\n3.99 → 3.45\n(↓13%)',fontsize=13,color='#b2182b',va='top')
        panel_path=output/'S4_A.png';fig.savefig(panel_path,dpi=300);plt.close(fig)
    panel=Image.open(panel_path).convert('RGB').resize((3561,2016),Image.Resampling.LANCZOS);panel.save(panel_path)
    full=base.copy();full.paste(panel,(0,0));assert ImageChops.difference(base.crop((0,2016,3561,3200)),full.crop((0,2016,3561,3200))).getbbox() is None
    dest=output/'S4.png';full.save(dest,dpi=(300,300))
    proof={'figure':'S4','template_sha256':sha(template),'summary_sha256':sha(record),'output_sha256':sha(dest),'rows':data['rows'],'precision':data['precision'],'panels_B_C_pixels_unchanged':True,'experiment_executed':False,'verification_level':'saved-summary plotting only','environment':{'matplotlib':matplotlib.__version__,'numpy':np.__version__}}
    (output/'S4_drawing_validation.json').write_text(json.dumps(proof,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    return proof

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--bundle',type=Path,default=Path(__file__).parent/'supplement');ap.add_argument('--output',required=True,type=Path);a=ap.parse_args();root=Path(__file__).resolve().parents[1];out=a.output.resolve()
    if out==root or root in out.parents or out.exists():ap.error('Output must be new and outside the code package')
    print(json.dumps(render_one(a.bundle,out),ensure_ascii=False))
