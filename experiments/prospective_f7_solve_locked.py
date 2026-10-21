"""Prospective F7 matched-mask comparison using the hash-bound archived solver.

Requires a completed prospective_f7_prepare.py input. This script has no legacy
output fallback, no parameter search, and does not render or replace paper figures.
"""
from __future__ import annotations
import argparse, hashlib, importlib.metadata, inspect, json, os, sys, gc
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
FINGERPRINT='2a6841164acdc771ebc8bcddb24ac766f17adea723852478dd8db6599f5cb132'
REFERENCE={'FGW':5.849528636550965,'LROT':5.817591851227303}
TOL_ENTROPY=0.001  # prespecified descriptive comparison; no image replacement

def sha(path:Path)->str:
 h=hashlib.sha256()
 with path.open('rb') as f:
  for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
 return h.hexdigest()

def section(z,prefix):
 return {'coords':z[prefix+'_coords'],'expr':z[prefix+'_expr'],
         'region_labels':z[prefix+'_region_labels'],
         'gene_names':z[prefix+'_gene_names'].tolist(),
         'barcodes':z[prefix+'_barcodes'].tolist(),
         'n_spots':int(z[prefix+'_coords'].shape[0]),
         'n_genes':int(z[prefix+'_expr'].shape[1])}

def identity(w,A,B,S,parameters):
 bound=w.signature.bind(A,B,S,**parameters);bound.apply_defaults()
 params=dict(list(bound.arguments.items())[3:])
 params.update(strict=False,return_diagnostics=True,verbose=False,inner_max_iter=10000,
               marginal_abs_tol=1e-8,marginal_rel_tol=1e-5,tol=1e-5)
 return {'A':{k:w.arrayhash(A[k]) for k in ['coords','expr','region_labels']},
         'B':{k:w.arrayhash(B[k]) for k in ['coords','expr','region_labels']},
         'gene_names_A':list(A['gene_names']),'gene_names_B':list(B['gene_names']),
         'lr':w.arrayhash(S),'params':w.clean(params),'solver':w.fingerprint}

def score(P):
 if P.shape!=(3000,3000) or not np.isfinite(P).all() or P.min() < -1e-12:
  raise RuntimeError('Plan shape/finite/nonnegative check failed')
 a=P.sum(axis=1);b=P.sum(axis=0)
 q=P/np.maximum(a[:,None],1e-300)
 return {'entropy':float((-q*np.log(np.maximum(q,1e-300))).sum(axis=1).mean()),
         'row_marginal_max_abs':float(np.max(np.abs(a-1/3000))),
         'col_marginal_max_abs':float(np.max(np.abs(b-1/3000))),
         'min_P':float(P.min()),'max_P':float(P.max()),'sum_P':float(P.sum())}

def main():
 pa=argparse.ArgumentParser();pa.add_argument('--prepared-dir',type=Path,required=True)
 pa.add_argument('--output-dir',type=Path,required=True);args=pa.parse_args()
 prep=args.prepared_dir.resolve();out=args.output_dir.resolve()
 if out.exists() or out==ROOT or ROOT in out.parents:
  raise RuntimeError('Solver output must be a fresh directory outside the repository')
 meta=json.loads((prep/'input_identity.json').read_text(encoding='utf-8'))
 input_path=prep/'prepared_f7_raw_input.npz'
 if sha(input_path)!=meta['prepared_npz_sha256']:
  raise RuntimeError('Prepared input SHA-256 changed')
 if sha(ROOT/'experiments/prospective_f7_prepare.py')!=meta['preparation_script_sha256']:
  raise RuntimeError('Preparation code changed')
 from paper_protocol.binding import bind,solve_locked,verify
 manifest=verify()
 if manifest['solver_fingerprint']!=FINGERPRINT:raise RuntimeError('Solver fingerprint changed')
 with np.load(input_path,allow_pickle=False) as z:
  A=section(z,'A');B=section(z,'B');S=z['S_lr']
 from paper_protocol.archived_wrapper import arrayhash
 for name,sec in [('A',A),('B',B)]:
  for field in ['coords','expr','region_labels']:
   if arrayhash(sec[field])!=meta['ordered_input_arrayhash'][name][field]:
    raise RuntimeError('Prepared array identity differs: '+name+'/'+field)
 if arrayhash(S)!=meta['ordered_input_arrayhash']['S_lr']:
  raise RuntimeError('Prepared LR matrix identity differs')
 w=bind(out)
 if w.fingerprint!=FINGERPRINT:raise RuntimeError('Wrong bound solver')
 sys.path.insert(0,str(ROOT/'experiments'))
 from run_cancer_visium import HUMAN_ALL_LR
 if len(HUMAN_ALL_LR)!=115:raise RuntimeError('Unexpected HUMAN_ALL_LR symbol count')
 methods=[('FGW',np.zeros_like(S),0.0),('LROT',S,0.1)]
 locked=[]
 for name,matrix,gamma in methods:
  params={'gamma':gamma,'beta':0.3,'alpha':0.5,'max_iter':200,
          'lr_gene_list':HUMAN_ALL_LR,'verbose':False}
  expected=identity(w,A,B,matrix,params)
  locked.append({'method':name,'gamma':gamma,'expected':expected})
  (prep/f'locked_identity_{name}.json').write_text(json.dumps(expected,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
 plan={'scope':'One prospective fixed-input breast comparison; not historical run recovery',
       'input_sha256':sha(input_path),'source_files':manifest['source_files'],
       'solver_fingerprint':FINGERPRINT,
       'parameters':{'FGW':locked[0]['expected']['params'],'LROT':locked[1]['expected']['params']},
       'reference_entropy':REFERENCE,'entropy_abs_tolerance':TOL_ENTROPY,
       'required':'both solver diagnostics accepted; 3000x3000 finite nonnegative P; marginal tolerances from bound solver',
       'no_paper_image_replacement':True,
       'environment':{x:importlib.metadata.version(x) for x in ['numpy','scipy','scikit-learn','threadpoolctl','h5py']},
       'runner_sha256':sha(Path(__file__))}
 (prep/'prespecified_solver_plan.json').write_text(json.dumps(plan,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
 print('PREFLIGHT PASSED: input, source, solver and both method identities locked',flush=True)
 results={}
 for name,matrix,gamma in methods:
  entry=next(x for x in locked if x['method']==name)
  print('SOLVE_START',name,'gamma',gamma,flush=True)
  P,loss=solve_locked(w,A,B,matrix,entry['expected'],
    gamma=gamma,beta=0.3,alpha=0.5,max_iter=200,
    lr_gene_list=HUMAN_ALL_LR,verbose=False)
  metric=score(P);metric['outer_iterations']=len(loss)
  metric['reference_entropy']=REFERENCE[name]
  metric['abs_difference']=abs(metric['entropy']-REFERENCE[name])
  metric['within_prespecified_tolerance']=metric['abs_difference']<=TOL_ENTROPY
  results[name]=metric
  (prep/f'prospective_metric_{name}.json').write_text(json.dumps(metric,indent=2)+'\n',encoding='utf-8')
  print('SOLVE_DONE',name,json.dumps(metric,sort_keys=True),flush=True)
  del P,loss,matrix;gc.collect()
 outcome={'results':results,'numeric_comparison_passed':all(x['within_prespecified_tolerance'] for x in results.values()),
          'solver_run_root':str(out),'input_sha256':plan['input_sha256'],
          'scope':'prospective source-bound numerical run only; figure image and manuscript unchanged'}
 (prep/'prospective_outcome.json').write_text(json.dumps(outcome,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
 print('OUTCOME',json.dumps(outcome,ensure_ascii=False),flush=True)

if __name__=='__main__':main()
