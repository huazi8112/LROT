"""Build an F7 drawing input from one accepted prospective FGW/LROT run.

Reads only saved arrays and diagnostics. Does not call a solver or replace a paper asset.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np

FINGERPRINT='2a6841164acdc771ebc8bcddb24ac766f17adea723852478dd8db6599f5cb132'

def sha(p:Path)->str:
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
 return h.hexdigest()

def entropy(P):
 if P.shape!=(3000,3000) or not np.isfinite(P).all() or P.min()<0:
  raise ValueError('Saved plan shape, finite or nonnegative check failed')
 q=P/np.maximum(P.sum(1,keepdims=True),1e-300)
 return float(-(q*np.log(np.maximum(q,1e-300))).sum(1).mean())

def main():
 p=argparse.ArgumentParser()
 p.add_argument('--prepared-dir',type=Path,required=True)
 p.add_argument('--solver-dir',type=Path,required=True)
 p.add_argument('--output-dir',type=Path,required=True)
 a=p.parse_args();prep=a.prepared_dir.resolve();solve=a.solver_dir.resolve();out=a.output_dir.resolve()
 if out.exists():raise ValueError('Output directory must be new')
 meta=json.loads((prep/'input_identity.json').read_text(encoding='utf-8'))
 prepared=prep/'prepared_f7_raw_input.npz'
 if sha(prepared)!=meta['prepared_npz_sha256']:raise ValueError('Prepared input SHA changed')
 run=solve/'runs/paper_protocol'
 calls=json.loads((run/'calls.json').read_text(encoding='utf-8'))
 if len(calls)!=2:raise ValueError('Expected exactly two locked solver calls')
 plans={};losses={};evidence={}
 for name,gamma,call in zip(('FGW','LROT'),(0.0,0.1),calls):
  j=run/'evidence'/(call['key']+'.json');n=run/'evidence'/(call['key']+'.npz')
  record=json.loads(j.read_text(encoding='utf-8'))
  ident=record['input_identity']
  if not record['accepted'] or not record['diagnostics']['success'] or len(record['attempts'])!=1:
   raise ValueError(name+' was not an accepted single-attempt run')
  if ident['solver']!=FINGERPRINT or float(ident['params']['gamma'])!=gamma:
   raise ValueError(name+' solver or gamma identity changed')
  for section in ('A','B'):
   for field in ('coords','expr','region_labels'):
    if ident[section][field]!=meta['ordered_input_arrayhash'][section][field]:
     raise ValueError(name+' input identity changed: '+section+'/'+field)
  if name=='LROT' and ident['lr']!=meta['ordered_input_arrayhash']['S_lr']:
   raise ValueError('LROT LR matrix identity changed')
  with np.load(n,allow_pickle=False) as z:
   plans[name]=z['P'];losses[name]=z['loss']
  evidence[name]={'record_sha256':sha(j),'plan_npz_sha256':sha(n),'entropy':entropy(plans[name]),
                  'outer_steps':int(len(losses[name]))}
 with np.load(prepared,allow_pickle=False) as z:
  A=z['A_coords'];B=z['B_coords']
 top=np.array(meta['top_pairs_preview'],dtype=object)
 if top.shape!=(10,3):raise ValueError('Expected ten saved top-pair records')
 out.mkdir(parents=True,exist_ok=False)
 plot=out/'F7_same_run_plot_input.npz'
 np.savez_compressed(plot,P=plans['LROT'],P_fgw=plans['FGW'],A_coords=A,B_coords=B,
                     loss_hist=losses['LROT'],top_pairs=top)
 fg=evidence['FGW']['entropy'];lr=evidence['LROT']['entropy']
 m={'FGW_entropy':fg,'LROT_entropy':lr,'reduction_percent':100*(fg-lr)/fg,
    'FGW_evidence_npz_sha256':evidence['FGW']['plan_npz_sha256'],
    'LROT_evidence_npz_sha256':evidence['LROT']['plan_npz_sha256'],
    'plot_input_sha256':sha(plot),'solver_fingerprint':FINGERPRINT,
    'scope':'One prospective accepted raw-input run; does not prove historical same-run provenance'}
 (out/'F7_same_run_metrics.json').write_text(json.dumps(m,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
 (out/'F7_same_run_evidence.json').write_text(json.dumps({'input_sha256':sha(prepared),
    'methods':evidence,'plot_input_sha256':sha(plot),'metrics_json_sha256':sha(out/'F7_same_run_metrics.json')},indent=2)+'\n',encoding='utf-8')
 print('plot_input_sha256',sha(plot));print('entropy',fg,lr);print('no_solver_called',True)

if __name__=='__main__':main()
