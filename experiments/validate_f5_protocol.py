from pathlib import Path
import argparse,hashlib,json,sys,importlib.util,csv
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--repo-root',type=Path,required=True);p.add_argument('--data-dir',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args()
r=a.repo_root.resolve();sys.path.insert(0,str(r))
from paper_protocol.binding import bind
ref=json.loads((r/'figures/current_release/inputs/F5_input_binding.json').read_text(encoding='utf-8'))
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for c in iter(lambda:f.read(1048576),b''):h.update(c)
 return h.hexdigest()
def ah(x):return hashlib.sha256(np.ascontiguousarray(x).tobytes()).hexdigest()
for n,h in ref['input_file_sha256'].items():
 if sha(a.data_dir/n)!=h:raise RuntimeError('Raw hash mismatch: '+n)
source=r/ref['source_script']
if sha(source)!=ref['source_script_sha256']:raise RuntimeError('Frozen F5 loader source changed')
w=bind(a.output_dir);out=a.output_dir
spec=importlib.util.spec_from_file_location('f5_fixed',source);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
sections=m.load_sections(a.data_dir);checks={}
import h5py
orders=[]
for sid,sec in sections.items():
 identity=ref['ordered_section_identities'][sid]
 for field,key in [('coords','coords_sha256'),('expr','expression_sha256'),('region_labels','labels_sha256')]:
  if ah(sec[field])!=identity[key]:raise RuntimeError(sid+' '+field+' differs')
 if hashlib.sha256('\n'.join(sec['gene_names']).encode()).hexdigest()!=identity['gene_names_sha256']:raise RuntimeError('gene order differs')
 positions={line.split(',')[0] for line in (a.data_dir/(sid+'_tissue_positions_list.txt')).read_text().splitlines() if len(line.split(','))==6}
 layers={v[0] for line in (a.data_dir/'barcode_level_layer_map.tsv').read_text().splitlines() if len(v:=line.split('\t'))==3 and v[1]==sid}
 with h5py.File(a.data_dir/(sid+'_filtered_feature_bc_matrix.h5'),'r') as f:
  bars=[v.decode() for v in f['matrix/barcodes'][:] if v.decode() in positions & layers]
 idx=np.random.RandomState(m.SECTION_SEEDS[sid]).choice(len(bars),1000,replace=False);selected=[bars[i] for i in idx]
 if ah(idx)!=identity['index_sha256'] or hashlib.sha256('\n'.join(selected).encode()).hexdigest()!=identity['barcode_sha256']:raise RuntimeError('barcode order differs')
 orders.extend((sid,i,int(j),bars[j]) for i,j in enumerate(idx));checks[sid]=True
with open(out/'ordered_barcodes.csv','w',newline='',encoding='utf-8') as f:
 c=csv.writer(f);c.writerow(['section','row','eligible_index','barcode']);c.writerows(orders)
(out/'raw_input_checks.json').write_text(json.dumps(checks,indent=2));print('F5 RAW INPUT CHECKS PASSED',flush=True)
w.core.fgw_lr_solver=w.solver
pooled,metrics=m.compute_ranking_data(r,sections)
np.savez_compressed(out/'per_spot.npz',**{name+'_'+field:values for name,fields in pooled.items() for field,values in fields.items()})
fig,txt=m.make_figure(pooled,metrics,out/'figures')
result={}
from scipy.stats import spearmanr
for name,fields in pooled.items():
 score=fields['score'];correct=fields['correct'];top=np.argsort(score)[::-1][:400]
 result[name]={'rho':float(spearmanr(score,correct).statistic),'top20_accuracy':float(correct[top].mean()),'accuracy':float(correct.mean())}
(out/'outcome.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
