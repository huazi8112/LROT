"""Draw explicitly selected saved-result panels on immutable Oct04 canvases.

This is saved-record plotting, not experimental reproduction. No solver, raw
array, training dependency or external file is read. Unselected pixels are
copied exactly from the locked canvas. S19 is deliberately unsupported.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np
from PIL import Image, ImageChops, ImageOps

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def new_panel(size, box, dpi, font=9):
    w,h=size; x0,y0,x1,y1=box
    fig=plt.figure(figsize=((w+0.001)/dpi,(h+0.001)/dpi),dpi=dpi,facecolor='white')
    ax=fig.add_axes([x0/w,(h-y1)/h,(x1-x0)/w,(y1-y0)/h])
    ax.tick_params(labelsize=font)
    return fig,ax

def finish(fig, path, dpi):
    fig.savefig(path,dpi=dpi,facecolor='white',metadata={'Software':'LROT selected saved-record panel export'})
    plt.close(fig)

def line_style(ax):
    ax.grid(True,color='#e5e5e5',linewidth=.6,alpha=.7)
    ax.set_axisbelow(True)

def draw_panel(fid,panel,data,size,axes,dpi,out):
    font=10 if fid in ('S2','S4','S7','S15') else 9
    with plt.rc_context({'font.family':'DejaVu Sans','font.size':font,
                         'axes.titlesize':font+1,'axes.labelsize':font,
                         'legend.fontsize':font-1,'lines.linewidth':1.5}):
        fig,ax=new_panel(size,axes,dpi,font)
        report={}
        if fid=='S2':
            rr=data['resource_rows'] if isinstance(data,dict) and 'resource_rows' in data else []; colors=['#427db8','#dc7868','#68bb65']; x=np.arange(3)
            if panel in ('A','C'):
                g=.1 if panel=='A' else 1.
                rows=[next(r for r in rr if r['resource']==rname and r['gamma']==g)
                      for rname in ['Full LR Database','Shuffled LR','Adhesion Only']]
                vals=[r['accuracy'] if panel=='A' else r['entropy'] for r in rows]
                ax.bar(x,vals,color=colors,width=.3,edgecolor='gray',linewidth=.5)
                ax.set_xticks(x,['Full','Shuffled','Adhesion'])
                ax.set_ylim(0,1 if panel=='A' else 5.4)
                ax.set_ylabel('Region accuracy' if panel=='A' else 'Mean row entropy')
                ax.set_title('A: Accuracy (γ=0.1)' if panel=='A' else '(C) Row entropy: expression-free LR + geometry')
                for xx,v in zip(x,vals): ax.text(xx,v+.025 if panel=='A' else v+.06,f'{v:.3f}',ha='center',fontsize=font-1,color='black' if panel=='A' else '#C93D3D',fontweight='bold')
                if panel=='A': ax.axhline(1/6,color='#999999',ls='--',lw=.8,label='Random 1/6');ax.legend(loc='upper right')
                else:
                    for xx,r in zip(x[1:],rows[1:]): ax.text(xx,r['entropy']-.18,f"{(r['entropy']/rows[0]['entropy']-1)*100:+.0f}%",ha='center',va='top',fontsize=font-2,color='white')
                report={'rows':rows,'plotted_values':vals,'gamma':g}
            elif panel=='E':
                vals=[]
                for g,c,off in [(.1,'#5fa3cc',-.15),(1.,'#dc7868',.15)]:
                    rows=[next(r for r in rr if r['resource']==name and r['gamma']==g) for name in ['Full LR Database','Shuffled LR','Adhesion Only']]
                    y=[r['fraction_P_gt_1e_minus10'] for r in rows];ax.bar(x+off,y,width=.3,color=c,label=f'γ={g:g}');vals.append(y)
                ax.set_xticks(x,['Full','Shuffled','Adhesion']);ax.set_ylim(0,.68);ax.set_ylabel('Fraction of entries > 10⁻¹⁰');ax.set_title('E: Transport plan support');ax.legend(loc='upper right')
                report={'gamma_order':[.1,1.],'plotted_values':vals,'threshold':1e-10,'independent_recount':False}
            elif panel=='D':
                # S2D is drawn from the six accepted diagnostics already saved
                # with this package. This changes presentation only; no solver
                # is imported or called here.
                rows=data
                resources=['Full LR Database','Shuffled LR','Adhesion Only']
                expected={(name,gamma) for name in resources for gamma in (.1,1.)}
                observed={}
                for row in rows:
                    gamma=float(row['gamma'])
                    alpha=float(row['alpha'])
                    params=row['parameters']
                    if row.get('accepted') is not True:
                        raise ValueError('S2D requires accepted saved diagnostics only')
                    if int(row['outer_iterations']) < 1:
                        raise ValueError('S2D outer_iterations must be positive')
                    if float(params['gamma']) != gamma or float(params['alpha']) != alpha:
                        raise ValueError('S2D parameter fields do not match the saved row')
                    if float(params['beta']) != .3:
                        raise ValueError('Unexpected S2D beta in saved diagnostics')
                    key=(row['resource'],gamma)
                    if key in observed:
                        raise ValueError('Duplicate S2D resource/gamma row')
                    observed[key]=row
                if set(observed) != expected or len(rows) != 6:
                    raise ValueError('S2D requires exactly the six recorded resource/gamma rows')
                x=np.arange(len(resources));width=.34;plotted=[]
                for gamma,color,offset,label in [(.1,'#5fa3cc',-.18,'γ = 0.1'),(1.,'#dc7868',.18,'γ = 1')]:
                    values=[int(observed[(name,gamma)]['outer_iterations']) for name in resources]
                    ax.bar(x+offset,values,width=width,color=color,label=label,
                           edgecolor='white',linewidth=.45)
                    plotted.append({'gamma':gamma,'values':values})
                all_values=[int(row['outer_iterations']) for row in rows]
                upper=(max(all_values)//20+1)*20
                ax.set_ylim(0,upper)
                ax.set_yticks(np.arange(0,upper+1,40))
                ax.set_xticks(x,['Full','Shuffled','Adhesion'])
                ax.set_xlabel('LR resource')
                ax.set_ylabel('Accepted outer iterations')
                ax.set_title('D: Accepted outer iterations')
                ax.legend(title='LR weight γ',loc='upper left',fontsize=font-1,
                          title_fontsize=font-1,frameon=False)
                line_style(ax)
                report={'resource_order':resources,'plotted_values':plotted,
                        'diagnostic':'accepted outer_iterations; not elapsed time',
                        'source_row_count':len(rows)}
            else:
                vals=[]
                for name,color,ls,label in [('LR_DB (61对)','#2166ac','-','Curated DB (61 pairs, 22 effective)'),('CellChatDB (2059对)','#d6604d','--','CellChatDB (2059 pairs, 38 effective)')]:
                    rows=[r for r in data['cellchat_rows'] if r[0]==name];xx=[r[1] for r in rows];yy=[r[2] for r in rows]
                    ax.plot(xx,yy,marker='o' if ls=='-' else 's',color=color,ls=ls,label=label);vals.append({'resource':name,'x':xx,'accuracy':yy})
                ax.set_xlim(.095,.205);ax.set_ylim(.705,.715);ax.set_xlabel('LR weight γ');ax.set_ylabel('Alignment accuracy');ax.set_title('(H) Alignment Accuracy',fontsize=14,fontweight='bold');ax.set_xlabel('LR weight γ',fontsize=12);ax.set_ylabel('Alignment accuracy',fontsize=12);ax.tick_params(labelsize=11);ax.legend(loc='upper left',fontsize=10);line_style(ax);report={'points':vals}
        elif fid=='S4':
            rows=data['mean_rows'];x=[r['gamma'] for r in rows]
            if panel=='B':
                for key,label,c,m in [('Neurotrophin / Growth Factor','Neurotrophin','#2166ac','o'),('Ephrin / Semaphorin','Ephrin','#d6604d','s'),('Chemokine / Cytokine','Chemokine','#4daf4a','^'),('Non-LR Control','Non-LR control','#984ea3','d')]: ax.plot(x,[r[key] for r in rows],color=c,marker=m,ms=4,label=label)
                ax.set_ylim(0,.5);ax.set_ylabel('Pathway activity correlation');ax.legend(loc='lower right',fontsize=7.5);ax.axvline(.1,c='#aaaaaa',ls=':',lw=.7);ax.text(.108,.03,'γ=0.1',c='#777777',fontsize=7)
                ax.set_title('(B) DLPFC: pathway-activity fidelity vs γ')
            else:
                ax.plot(x,[r['same_mass'] for r in rows],'-o',color='#2c3e50',ms=3,label='same-layer mass');ax.plot(x,[r['soft_acc'] for r in rows],'--s',color='#d6604d',ms=3,label='soft accuracy')
                ax.set_ylim(.325,.54);ax.set_ylabel('Layer concordance');a2=ax.twinx();a2.plot(x,[r['entropy'] for r in rows],'-o',color='#b2182b',ms=4,label='transport entropy');a2.set_ylim(3.8,6.3);a2.set_ylabel('Entropy (nats)',color='#b2182b');a2.tick_params(labelsize=font-1,labelcolor='#b2182b')
                h,l=ax.get_legend_handles_labels();h2,l2=a2.get_legend_handles_labels();ax.legend(h+h2,l+l2,loc='upper right',fontsize=7.5);ax.set_title('(C) DLPFC: layer concordance & entropy vs γ')
            ax.set_xlabel('LR Guidance Weight γ');line_style(ax);report={'points':rows,'aggregation':data['aggregation']}
        elif fid=='S5':
            aa=sorted({r['final_parameters']['alpha'] for r in data['points']});bb=sorted({r['final_parameters']['beta'] for r in data['points']})
            metric='entropy' if panel=='A' else 'soft_acc';z=np.array([[next(r['metrics'][metric] for r in data['points'] if r['final_parameters']['alpha']==a and r['final_parameters']['beta']==b) for b in bb] for a in aa])
            im=ax.imshow(z,aspect='auto',cmap='YlGnBu_r' if panel=='A' else 'YlGn',interpolation='nearest')
            ax.set_xticks(range(len(bb)),[f'β={b:g}' for b in bb]);ax.set_yticks(range(len(aa)),[f'α={a:g}' for a in aa]);ax.set_title('(A) LROT mean row entropy' if panel=='A' else '(B) LROT accuracy')
            for i in range(len(aa)):
                for j in range(len(bb)): ax.text(j,i,f'{z[i,j]:.3f}',ha='center',va='center',fontsize=font)
            w,h=size;x0,y0,x1,y1=axes;cax=fig.add_axes([(x1+55)/w,(h-y1)/h,45/w,(y1-y0)/h]);fig.colorbar(im,cax=cax).ax.tick_params(labelsize=font-1);report={'alphas':aa,'betas':bb,'metric':metric,'values':z.tolist()}
        elif fid=='S7':
            lines=[]
            for method,c,label in [('Expr OT only','#2ca02c','Expr OT only'),('FGW (γ=0)','#1f77b4','FGW (γ=0)'),('LROT (γ=0.1)','#d62728','LROT (γ=0.1)')]:
                rows=sorted([r for r in data['points'] if r['method']==method],key=lambda r:r['size']);x=[r['size'] for r in rows];y=[r['solve_time'] if panel=='A' else r['acc'] for r in rows];ax.plot(x,y,'-o',color=c,label=label,ms=4);lines.append({'method':method,'size':x,'values':y,'field':'solve_time' if panel=='A' else 'acc'})
            ax.set_xscale('log');ax.set_xlabel('Number of spots (n)');ax.set_ylabel('Outer solve time (s)' if panel=='A' else 'Mapping accuracy');ax.set_title('(A) Runtime scaling (log-log)' if panel=='A' else '(B) Accuracy vs. size');ax.legend(loc='upper left' if panel=='A' else 'upper right');line_style(ax)
            if panel=='A': ax.set_yscale('log')
            else: ax.set_ylim(.3,1)
            report={'lines':lines,'time_scope':data['time_scope'],'LR_preparation':'separate table column, not included in solve_time'}
        elif fid=='S8':
            order=data['method_order'];ys=[data['methods'][m]['values'] for m in order]
            bp=ax.boxplot(ys,patch_artist=True,widths=.52,showfliers=False,medianprops={'color':'#34495e','lw':1.4})
            for patch,c in zip(bp['boxes'],data['baseline_box_fill_colors']):patch.set_facecolor(c);patch.set_edgecolor('#5e6a71')
            rng=np.random.RandomState(0)
            for i,(yy,c) in enumerate(zip(ys,data['baseline_box_fill_colors']),1):
                ax.scatter(i+rng.uniform(-.18,.18,len(yy)),yy,s=16,c='white',edgecolors='#8b969c',linewidths=.65,zorder=3)
                for obj in bp['whiskers'][2*(i-1):2*i]+bp['caps'][2*(i-1):2*i]: obj.set_color(c)
            labels=['Soft FGW','LROT','PASTE2 (KL)','PASTE2 (Euclidean)'];ax.set_xticks(range(1,5),[f"{lab}\n{data['methods'][m]['mean']:.3f} ± {data['methods'][m]['sample_SD']:.3f}" for lab,m in zip(labels,order)])
            ax.tick_params(axis='x',length=0,pad=9);ax.set_ylim(0,1);ax.set_yticks(np.arange(0,1.01,.1));ax.set_xlim(.5,4.5);ax.grid(axis='y',c='#e3e7e9',lw=.7);ax.set_axisbelow(True);ax.spines[['top','right']].set_visible(False);ax.axhline(1/6,c='#999999',lw=.8);ax.text(4.49,1/6+.007,'Random reference = 1/6',ha='right',fontsize=7,c='#777777')
            fig.text(.5,.918,'All observations shown: mean ± sample SD',ha='center',fontsize=8)
            fig.text(.5,.088,'LROT > both tested PASTE2 configurations in all 20 paired seeds.',ha='center',fontsize=8)
            fig.text(.5,.055,'Descriptive comparison; no new experiments or spot-independent inference.',ha='center',fontsize=7,c='#666666')
            report={'method_order':order,'seed_order':data['seed_order'],'values':ys,'mean_and_sample_SD':{m:{'mean':float(np.mean(data['methods'][m]['values'])),'sample_SD':float(np.std(data['methods'][m]['values'],ddof=1))} for m in order},'external_parameters':data['external_parameters']}
        elif fid=='S9':
            y=[r['outer_iterations'] for r in data['points']];ax.bar(np.array([0,1])-.18,[y[0],y[2]],width=.36,color='#e89595',label='FGW (γ=0)');ax.bar(np.array([0,1])+.18,[y[1],y[3]],width=.36,color='#439668',label='LROT (γ=0.1)')
            for xx,v in zip([-.18,.18,.82,1.18],y):ax.text(xx,v+.45,str(v),ha='center',fontsize=8)
            ax.set_ylim(0,24);ax.set_xticks([0,1],['Sampled\n(1,000 spots)','Full-scale\n(4,221/4,381)']);ax.set_ylabel('Accepted outer iterations');ax.set_title('(C) Stopping diagnostics');ax.spines[['top','right']].set_visible(False);ax.legend(loc='upper right',fontsize=7);report={'points':data['points'],'values':y,'metric':data['metric']}
        elif fid=='S15':
            if panel=='B':
                for metric,c,label in [('C_expr','#c0392b','C_expr'),('C_lr','#2c3e50','C_lr'),('C_gw','#95A5A6','C_gw')]:
                    for model,ls,model_label in [('gene','--','lognormal'),('geneadd','-','log1p additive')]:
                        rr=sorted([r for r in data['cost_change_means'] if r['model']==model],key=lambda r:r['sigma']);ax.plot([r['sigma'] for r in rr],[r['mean_absolute_changes'][metric]['mean'] for r in rr],ls+'o',color=c,ms=3,label=f'{label} ({model_label})')
                ax.set_xlim(.29,.51);ax.set_ylim(0,.11);ax.set_ylabel('Mean |C(sigma) − C(0)|');ax.set_xlabel('Batch strength sigma');ax.set_title('(B) Cost-matrix perturbation entering solver');ax.legend(loc='upper left',ncol=2,fontsize=7);report={'means':data['cost_change_means'],'errorbars':False,'Cg_definition':data['Cg_definition']}
            else:
                x=sorted(map(float,data['stress_grid']));factor=math.sqrt(data['n_unique_seeds']/(data['n_unique_seeds']-1));isacc=panel=='C';idx=0 if isacc else 2;lines=[]
                for a,c in [(.3,'#2C3E50'),(.5,'#C0392B'),(.7,'#16A085'),(.9,'#8E44AD')]:
                    rows=[data['stress_grid'][str(s)][str(a)] for s in x];means=[r[idx] for r in rows];sd=[r[idx+1]*factor for r in rows];ax.errorbar(x,means,yerr=sd,fmt='-o',ms=3,capsize=2,color=c,label=f'alpha={a:g} (LR)');lines.append({'alpha':a,'x':x,'mean':means,'sample_SD':sd})
                rows=[data['fgw_reference'][str(s)] for s in x];fi=2 if isacc else 6;means=[r[fi] for r in rows];sd=[r[fi+1]*factor for r in rows];ax.errorbar(x,means,yerr=sd,fmt='--s',ms=3,capsize=2,c='#7F8C8D',label='FGW alpha=0.5');lines.append({'method':'FGW','x':x,'mean':means,'sample_SD':sd})
                ax.set_xlim(-.075,1.575);ax.set_ylim(.275,.83) if isacc else ax.set_ylim(2.3,6.3)
                if isacc:ax.axvspan(0,.5,color='#EAF5EC',alpha=.65,zorder=-1);ax.text(.25,.28,'fixed-weight scan range',ha='center',fontsize=6,c='#61A56B')
                ax.set_title('(C) Accuracy: where fixed weights fail' if isacc else '(D) Mapping sharpness under batch stress');ax.set_ylabel('Region mapping accuracy' if isacc else 'Transport entropy');ax.set_xlabel('Batch strength sigma (log1p additive)');ax.legend(loc='center left' if isacc else 'lower right',fontsize=7);report={'lines':lines,'SD_definition':data['SD_definition_CD'],'n':data['n_unique_seeds']}
            line_style(ax)
        elif fid=='S17':
            vals=[]
            for g,c,label in [(0.,'#4E79A7','FGW (γ=0)'),(.1,'#F28E2B','LROT (γ=0.1)')]:
                rows=sorted([r for r in data['rows'] if r['gamma']==g],key=lambda r:r['eps_mult']);ax.plot([r['eps_mult'] for r in rows],[r['entropy'] for r in rows],'-o',c=c,ms=4,label=label)
                for r in rows:ax.annotate(f"acc={r['acc']:.3f}",(r['eps_mult'],r['entropy']),xytext=(3,3 if g==0 else -10),textcoords='offset points',fontsize=7,color='#444444')
                vals.extend(rows)
            ax.set_xlim(.425,2.2);ax.set_ylim(2,5.65);ax.set_xticks([.5,.75,1,1.25,1.5,1.75,2]);ax.set_xlabel('Regularization multiplier (both reg_ot and reg_gw)');ax.set_ylabel('Mean row entropy');ax.set_title('(C) Regularization scan (realistic, seed42)',fontsize=10);ax.legend(loc='upper left',fontsize=8);ax.spines[['top','right']].set_visible(False)
            red=[]
            for m in [.5,1.,2.]:
                f=next(r['entropy'] for r in data['rows'] if r['gamma']==0 and r['eps_mult']==m);l=next(r['entropy'] for r in data['rows'] if r['gamma']==.1 and r['eps_mult']==m);red.append(100*(f-l)/f)
            ax.text(.975,2.05,'Entropy reduction '+ '/'.join(f'{v:.2f}%' for v in red),fontsize=6,color='#777777');report={'points':vals,'entropy_reduction_percent':red,'multiplier_definition':data['multiplier_definition']}
        else: raise ValueError(fid)
        finish(fig,out,dpi);return report


def apply_title_overlays(image, rules, dpi):
    """Replace declared title bands only, preserving every other pixel."""
    from PIL import ImageDraw, ImageFont
    from matplotlib import font_manager
    draw = ImageDraw.Draw(image)
    font_path = font_manager.findfont(font_manager.FontProperties(family='DejaVu Sans', weight='normal'))
    font = ImageFont.truetype(font_path, round(11 * dpi / 72))
    for rule in rules:
        x0, y0, x1, y1 = rule['box']
        draw.rectangle((x0, y0, x1-1, y1-1), fill='white')
        center = ((x0+x1)/2, (y0+y1)/2)
        bounds = draw.textbbox(center, rule['text'], font=font, anchor='mm')
        if not (x0 <= bounds[0] and y0 <= bounds[1] and bounds[2] < x1 and bounds[3] < y1):
            raise ValueError('Title does not fit its declared band')
        draw.text(center, rule['text'], font=font, fill='black', anchor='mm')
    return image

def render_one(fid,bundle,output):
    if fid=='S4':
        from render_s4_saved import render_one as render_s4
        return render_s4(bundle,output)
    bundle=Path(bundle);output=Path(output);output.mkdir(parents=True,exist_ok=True)
    configs=json.loads((bundle/'selected_panel_config.json').read_text(encoding='utf-8'));cfg=next(c for c in configs if c['id']==fid)
    template=bundle/cfg['template_path']
    if sha(template)!=cfg['template_sha256']:raise ValueError('Baseline template hash mismatch')
    input_path=bundle/'inputs'/f'{fid}.json';data=json.loads(input_path.read_text(encoding='utf-8'))
    base=Image.open(template).convert('RGB');full=base.copy();mask=Image.new('L',base.size,0);logs=[]
    for p in cfg['panels']:
        b=cfg['panel_boxes'][p];sz=(b[2]-b[0],b[3]-b[1]);dest=output/f'{fid}_{p}.png'
        if fid=='S2' and p=='H':
            # All four selected values match the baseline display. Keep pixels;
            # a new driver/context binding alone is not a reason to redraw.
            base.crop(b).save(dest)
            logs.append({'panel':'H','operation':'retained baseline pixels; selected saved values unchanged',
                         'png_sha256':sha(dest),'selected_values':data['cellchat_rows'],
                         'new_drawing':False})
            continue
        plotted=draw_panel(fid,p,data,sz,cfg['axes_boxes'][p],cfg['dpi'],dest)
        image=Image.open(dest).convert('RGB');assert image.size==sz,(image.size,sz)
        full.paste(image,(b[0],b[1]));mask.paste(255,b)
        logs.append({'panel':p,'box':b,'png_sha256':sha(dest),'plotted':plotted})
    if fid=='S2':
        d_cfg=cfg['panel_data']['D']
        d_input=bundle/d_cfg['path']
        if sha(d_input)!=d_cfg['sha256']:
            raise ValueError('S2D saved-diagnostic input hash mismatch')
        d_data=json.loads(d_input.read_text(encoding='utf-8'))
        b=cfg['panel_boxes']['D'];sz=(b[2]-b[0],b[3]-b[1])
        dest=output/'S2_D.png'
        plotted=draw_panel('S2','D',d_data,sz,cfg['axes_boxes']['D'],cfg['dpi'],dest)
        image=Image.open(dest).convert('RGB');assert image.size==sz,(image.size,sz)
        full.paste(image,(b[0],b[1]));mask.paste(255,b)
        logs.append({'panel':'D','box':b,'png_sha256':sha(dest),
                     'input_path':d_cfg['path'],'input_json_sha256':sha(d_input),
                     'operation':'redrawn from six accepted saved diagnostics on standard panel geometry',
                     'new_drawing':True,'plotted':plotted})
    title_rules=cfg.get('title_overlays', [])
    apply_title_overlays(full, title_rules, cfg['dpi'])
    for rule in title_rules: mask.paste(255, tuple(rule['box']))
    if fid=='S15':
        from repair_supplement_annotations import repair_s15
        full, annotation_boxes = repair_s15(full)
        for box in annotation_boxes: mask.paste(255, box)
        logs.append({'operation':'S15 annotation correction only','boxes':annotation_boxes,'scientific_values_changed':False})
    delta=ImageChops.difference(base,full);outside=ImageChops.multiply(delta,ImageOps.invert(mask).convert('RGB'));assert outside.getbbox() is None,'Protected pixels changed'
    dest=output/f'{fid}.png';full.save(dest,dpi=(cfg['dpi'],cfg['dpi']),software='LROT selected saved-record panel export')
    prohibited=cfg['scientific_imports_prohibited'];assert not any(m in sys.modules for m in prohibited)
    proof={'figure':fid,'template_sha256':sha(template),'input_json_sha256':sha(input_path),'output_sha256':sha(dest),'size':list(full.size),'protected_pixels_unchanged':True,'title_overlays':title_rules,'panels':logs,'verification_level':'saved-record plotting; no raw-array verification or experiment execution','environment':{'python':sys.version,'matplotlib':matplotlib.__version__,'numpy':np.__version__}}
    (output/f'{fid}_drawing_validation.json').write_text(json.dumps(proof,ensure_ascii=False,indent=2),encoding='utf-8');return proof

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--bundle',type=Path,default=Path(__file__).parent/'supplement');ap.add_argument('--output',type=Path,required=True);ap.add_argument('--figure',choices=['S2','S4','S5','S7','S8','S9','S15','S17','all'],default='all');a=ap.parse_args()
    ids=[a.figure] if a.figure!='all' else ['S2','S4','S5','S7','S8','S9','S15','S17'];proof=[render_one(fid,a.bundle,a.output) for fid in ids];(a.output/'selected_drawing_validation.json').write_text(json.dumps(proof,ensure_ascii=False,indent=2),encoding='utf-8')
if __name__=='__main__':main()
