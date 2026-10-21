"""Redraw S18 from the published Table S18 values, without running alignment."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
ROOT=Path(__file__).resolve().parents[2]
def render(source,output):
 data=json.loads(Path(source).read_text(encoding='utf-8'))
 pairs=list(dict.fromkeys(r['pair'] for r in data['rows']))
 rows={pair:sorted([{k:float(v) if k!='pair' else v for k,v in r.items()} for r in data['rows'] if r['pair']==pair],key=lambda r:r['gamma']) for pair in pairs}
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11})
 fig,axes=plt.subplots(3,2,figsize=(11,14.4));x=np.arange(4);labels=['D1-P1\n151507–151508','D1-P2\n151508–151509','D2-P1\n151669–151670','D2-P2\n151670–151671'];styles=[('#4e79a7','-','o'),('#4e79a7','-','s'),('#f28e2b','--','o'),('#f28e2b','--','s')]
 ax=axes[0,0];fg=[rows[q][0]['entropy'] for q in pairs];lr=[next(r['entropy'] for r in rows[q] if r['gamma']==.1) for q in pairs]
 for offset,vals,color,label in [(-.19,fg,'#e89595','FGW (γ=0)'),(.19,lr,'#439668','LROT (γ=0.1)')]:
  bars=ax.bar(x+offset,vals,.38,color=color,label=label);ax.bar_label(bars,fmt='%.2f',padding=3,fontsize=8)
 for j,(u,v) in enumerate(zip(fg,lr)):ax.text(j,max(u,v)+.45,f'−{100*(u-v)/u:.1f}%',ha='center',fontsize=8)
 ax.set_ylim(0,8);ax.set_ylabel('Mean row entropy (nats)');ax.set_title('(A) Entropy at γ=0.1');ax.set_xticks(x,labels,fontsize=8);ax.legend(fontsize=8,loc='upper left')
 for ax,metric,title,ylabel in [(axes[0,1],'entropy','(B) Normalized entropy','Entropy (% of γ=0)'),(axes[1,0],'same_mass','(C) Same-layer mass','Same-layer transported mass'),(axes[1,1],'soft_ari','(D) Soft ARI','Soft-vote ARI')]:
  for j,q in enumerate(pairs):
   rr=rows[q];ys=[r[metric] for r in rr]
   if metric=='entropy':ys=[100*v/ys[0] for v in ys]
   color,ls,marker=styles[j];ax.plot([r['gamma'] for r in rr],ys,color=color,ls=ls,marker=marker,ms=4,label=labels[j].replace('\n',' '))
  ax.set_title(title);ax.set_xlabel('LR guidance weight γ');ax.set_ylabel(ylabel);ax.set_xticks([0,.02,.05,.1,.15,.2],['0','0.02','0.05','0.1','0.15','0.2'])
  if metric=='entropy':ax.legend(fontsize=7,loc='best')
 deltas=[]
 for ax,g,panel in [(axes[2,0],.1,'E'),(axes[2,1],.05,'F')]:
  dm=[];da=[]
  for q in pairs:
   base=rows[q][0];rr=next(r for r in rows[q] if r['gamma']==g);dm.append(base['same_mass']-rr['same_mass']);da.append(base['soft_ari']-rr['soft_ari'])
  for offset,ys,col,fmt in [(-.19,dm,'#4e79a7','%+.4f'),(.19,da,'#f28e2b','%+.3f')]:
   bars=ax.bar(x+offset,ys,.38,color=col);annotations=ax.bar_label(bars,fmt=fmt,padding=5,fontsize=7.5)
   if panel=='E' and offset>0:annotations[3].set_position((0,18))
  lo=min(0,*dm,*da);hi=max(0,*dm,*da);margin=max(.012,(hi-lo)*.28);ax.set_ylim(lo-margin,hi+margin);ax.axhline(0,color='#888',lw=.9);ax.set_xticks(x,labels,fontsize=8);ax.set_ylabel('FGW (γ=0) − LROT');ax.set_title(f'({panel}) Differences at γ={g:g}')
  deltas.append({'gamma':g,'FGW_minus_LROT_same_mass':dm,'FGW_minus_LROT_soft_ARI':da})
 for ax in axes.flat:ax.grid(axis='y',color='#ddd',alpha=.65,lw=.6);ax.set_axisbelow(True);ax.spines[['top','right']].set_visible(False)
 fig.legend(handles=[Patch(color='#4e79a7',label='FGW−LROT same-layer mass (E/F)'),Patch(color='#f28e2b',label='FGW−LROT soft ARI (E/F)')],loc='lower center',ncol=2,frameon=False,fontsize=9,bbox_to_anchor=(.5,.018))
 fig.subplots_adjust(left=.09,right=.98,top=.97,bottom=.09,hspace=.48,wspace=.30)
 fig.savefig(output,dpi=260);plt.close(fig)
 return {'input_sha256':hashlib.sha256(Path(source).read_bytes()).hexdigest(),'output_sha256':hashlib.sha256(Path(output).read_bytes()).hexdigest(),'scope':data['scope'],'deltas':deltas,'data_precision':data['precision']}
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--input',type=Path,default=ROOT/'figures/supplement/inputs/S18_table_values.json');ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();dest=a.output.resolve()
 if ROOT==dest or ROOT in dest.parents:ap.error('Choose an output outside the repository')
 dest.parent.mkdir(parents=True,exist_ok=True);proof=render(a.input,dest);dest.with_suffix('.verification.json').write_text(json.dumps(proof,indent=2),encoding='utf-8')
if __name__=='__main__':main()
