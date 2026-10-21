"""Prepare the two paper DLPFC cases; stop before any model training."""
from pathlib import Path
import argparse,ast,csv,hashlib,json,os,socket,subprocess,sys,time,types
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--data-root',required=True,type=Path)
p.add_argument('--output-dir',required=True,type=Path)
a=p.parse_args();ROOT=Path(__file__).resolve().parents[1];OUT=a.output_dir.resolve()
if OUT.exists() or OUT==ROOT or ROOT in OUT.parents:p.error('Output must be new and outside the package.')
if os.environ.get('_LROT_DLPFC_PREP_CHILD')!='1':
 env=os.environ.copy();env.update(_LROT_DLPFC_PREP_CHILD='1',PYTHONHASHSEED='0',PYTHONDONTWRITEBYTECODE='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMBA_CACHE_DIR=str(OUT/'numba_cache'),MPLCONFIGDIR=str(OUT/'mpl_cache'))
 raise SystemExit(subprocess.call([sys.executable,'-B',__file__,*sys.argv[1:]],env=env))
OUT.mkdir(parents=True);sys.dont_write_bytecode=True
def deny(*a,**k):raise RuntimeError('Offline preparation')
socket.create_connection=deny;socket.socket.connect=deny;socket.socket.connect_ex=deny
import numpy as np
import scanpy as sc
from scipy import sparse
from sklearn.neighbors import NearestNeighbors
from threadpoolctl import threadpool_limits
threadpool_limits(1)
from dlpfc_input_loader import load_dlpfc_input,arrayhash,sha
def js(d,n,x):(d/n).write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n','utf-8')
class Prepared(Exception):pass
def capture(ad,*args,**kw):
 assert ad.shape==(2000,50) and set(ad.obs.columns)=={'batch_name'}
 np.savez_compressed(current/'actual_training_input.npz',X=ad.X.toarray(),edgeList=ad.uns['edgeList'],obs_names=np.array(ad.obs_names,dtype=str),batch=np.array(ad.obs['batch_name'],dtype=str))
 js(current,'preparation.json',{'complete':True,'training_executed':False,'actual_training_input_sha256':sha(current/'actual_training_input.npz'),'ordered_spots_sha256':sha(current/'ordered_spots.csv'),'recipe_run_pair_sha256':hashlib.sha256(recipe.encode()).hexdigest()})
 raise Prepared()
sta=types.SimpleNamespace(train_STAligner=capture)
case=ROOT/'provenance/f3_staligner/scientific_case.py';text=case.read_text('utf-8-sig');tree=ast.parse(text)
fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='run_pair');recipe=ast.get_source_segment(text,fn)
exec(compile(ast.Module(body=[fn],type_ignores=[]),str(case),'exec'),globals())
f=ROOT/'experiments/run_cancer_visium.py';tree=ast.parse(f.read_text('utf-8-sig'));fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='apply_misalignment')
exec(compile(ast.Module(body=[fn],type_ignores=[]),str(f),'exec'),globals())
for source,target in [('151507','151508'),('151508','151509')]:
 current=OUT/(source+'-'+target);current.mkdir()
 A,ai=load_dlpfc_input(a.data_root,source,42);B,bi=load_dlpfc_input(a.data_root,target,99);B['coords']=apply_misalignment(B['coords'],rotation_deg=2.,translation=.5)
 bi['arrays_after_target_transform']={k:arrayhash(B[k]) for k in ['expr','coords','region_labels']}
 js(current,'input_identity.json',{'A':ai,'B':bi})
 with (current/'ordered_spots.csv').open('w',encoding='utf-8',newline='') as f:
  w=csv.writer(f);w.writerow(['batch','section','row','barcode','label'])
  for batch,sd in [('A',A),('B',B)]:w.writerows((batch,sd['sample_id'],i,b,int(sd['region_labels'][i])) for i,b in enumerate(sd['barcodes']))
 np.savez_compressed(current/'ordered_labels.npz',A_labels=A['region_labels'],B_labels=B['region_labels'])
 try:run_pair(source+'-'+target,A,B)
 except Prepared:print('PREPARED',source+'-'+target,flush=True)
js(OUT,'completed.json',{'complete':True,'pairs':2,'training_steps':0})
