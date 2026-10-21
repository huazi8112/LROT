"""Prepare or rerun one explicitly recorded PASTE/PASTE2 comparison, offline."""
from pathlib import Path
import argparse,ast,hashlib,inspect,json,os,socket,subprocess,sys,time
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--record',required=True,help='Record ID from provenance/comparator_audit/experiments.json')
p.add_argument('--mode',choices=['inputs','run'],default='inputs')
p.add_argument('--data-root',type=Path,help='Required for a measured DLPFC record')
p.add_argument('--output-dir',required=True,type=Path)
a=p.parse_args();ROOT=Path(__file__).resolve().parents[1];OUT=a.output_dir.resolve()
if OUT.exists() or OUT==ROOT or ROOT in OUT.parents:p.error('Output must be new and outside the package.')
if os.environ.get('_LROT_EXTERNAL_CHILD')!='1':
 env=os.environ.copy();env.update(_LROT_EXTERNAL_CHILD='1',PYTHONHASHSEED='0',PYTHONDONTWRITEBYTECODE='1',OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',MKL_NUM_THREADS='2',NUMBA_CACHE_DIR=str(OUT/'numba_cache'),MPLCONFIGDIR=str(OUT/'mpl_cache'))
 raise SystemExit(subprocess.call([sys.executable,'-B',__file__,*sys.argv[1:]],env=env))
OUT.mkdir(parents=True);sys.dont_write_bytecode=True
def deny(*a,**k):raise RuntimeError('Offline run; networking disabled')
socket.create_connection=deny;socket.socket.connect=deny;socket.socket.connect_ex=deny
import numpy as np
from threadpoolctl import threadpool_limits
threadpool_limits(2)
base=ROOT/'provenance/comparator_audit'
records=json.loads((base/'experiments.json').read_text('utf-8'))
ref=next((x for x in records if x['record']==a.record),None)
if ref is None:p.error('Unknown record ID')
dataset=ref['dataset'];identity=json.loads((base/'datasets.json').read_text('utf-8'))[dataset]
def sha(f):return hashlib.sha256(Path(f).read_bytes()).hexdigest()
def h(x):
 x=np.ascontiguousarray(x);return hashlib.sha256(str(x.shape).encode()+str(x.dtype).encode()+x.tobytes()).hexdigest()
def js(n,x):(OUT/n).write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n','utf-8')
if dataset.startswith('realistic-seed'):
 snap=ROOT/'provenance/numerical_validation_snapshot/code';assert sha(snap/'real_st_loader.py')==identity['generator_source_sha256']
 sys.path.insert(0,str(snap));from real_st_loader import RealisticSTGenerator
 seed=int(dataset.replace('realistic-seed',''));pair=RealisticSTGenerator(n_spots_A=1000,n_spots_B=1000,n_genes=150,seed=seed).generate_paired_slices(rotation=5.,batch_effect=.15)
 A,B=pair.slice_A.to_dict(),pair.slice_B.to_dict()
else:
 if a.data_root is None:p.error('--data-root is required for DLPFC')
 from dlpfc_input_loader import load_dlpfc_input
 source,target=dataset.split('-');A,ai=load_dlpfc_input(a.data_root,source,42);B,bi=load_dlpfc_input(a.data_root,target,99)
 f=ROOT/'experiments/run_cancer_visium.py';tree=ast.parse(f.read_text('utf-8-sig'));fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='apply_misalignment')
 exec(compile(ast.Module(body=[fn],type_ignores=[]),str(f),'exec'),globals());B['coords']=apply_misalignment(B['coords'],rotation_deg=2.,translation=.5)
actual=[]
for sd in [A,B]:actual.append({'expression_hash':h(sd['expr'].astype(np.float64)),'spatial_hash':h(sd['coords'].astype(np.float64)),'gene_order_sha256':h(np.array(sd['gene_names'],dtype=str)),'shape':[sd['n_spots'],len(sd['gene_names'])]})
assert actual==identity['package_input_hashes'],'Reconstructed ordered inputs differ; output retained and no optimiser called.'
js('input_identity.json',{'record':a.record,'dataset':dataset,'BLAS_threads':2,'actual_package_inputs':actual,'matched_fixed_identity':True,'parameters':ref['parameters'],'labels_withheld_from_optimisation':True})
if a.mode=='inputs':print('INPUTS VERIFIED',a.record);raise SystemExit(0)
import anndata as an
from sklearn.metrics import adjusted_rand_score
def to_ad(sd):
 ad=an.AnnData(X=sd['expr'].astype(np.float64));ad.var_names=list(sd['gene_names']);ad.obs_names=[f'spot{i}' for i in range(sd['n_spots'])];ad.obsm['spatial']=sd['coords'].astype(np.float64);return ad
env=json.loads((base/'environment.json').read_text('utf-8'));params=ref['parameters']
if ref['method']=='PASTE':
 from paste_math import pairwise_align as fn
 pin=env['implementations']['paste_math/PASTE.py']
 kwargs={'alpha':params['alpha'],'dissimilarity':params['dissimilarity'],'numItermax':params.get('numItermax',200),'verbose':False}
else:
 from paste2.PASTE2 import partial_pairwise_align as fn
 pin=env['implementations']['paste2/PASTE2.py']
 kwargs={'s':params['s'],'alpha':params['alpha'],'dissimilarity':params['dissimilarity'],'a_distribution':np.ones(1000)/1000,'b_distribution':np.ones(1000)/1000,'verbose':False}
assert sha(inspect.getsourcefile(fn))==pin,'External implementation differs; no silent package substitution.'
start=time.perf_counter();plan=fn(to_ad(A),to_ad(B),**kwargs);seconds=time.perf_counter()-start
np.savez_compressed(OUT/'plan.npz',P=plan)
classes=np.eye(int(max(A['region_labels'].max(),B['region_labels'].max()))+1)[B['region_labels']];pred=(plan@classes).argmax(1)
expected=np.load(base/'plans.npz',allow_pickle=False)[ref['plan_key']]
result={'record':a.record,'accuracy':float((pred==A['region_labels']).mean()),'ARI':float(adjusted_rand_score(A['region_labels'],pred)),'seconds':seconds,'max_abs_saved_plan_difference':float(np.max(np.abs(plan-expected))),'source_file_sha256':pin,'timing_is_new_benchmark':False}
js('result.json',result);assert np.allclose(plan,expected,rtol=0,atol=1e-12),'Fresh result differs; retained, not adopted automatically.'
print(json.dumps(result))
