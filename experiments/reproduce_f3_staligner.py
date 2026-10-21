from pathlib import Path
import os,sys,json,time,hashlib,ast,inspect,importlib,shutil,socket
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'provenance/f3_staligner'
import argparse,subprocess
parser=argparse.ArgumentParser(description='Verify or execute the pinned deterministic F3 STAligner case; no network or installation.')
parser.add_argument('--mode',choices=['verify','train'],default='verify')
parser.add_argument('--output-dir',required=True,type=Path)
args=parser.parse_args()
OUT=args.output_dir.resolve()
if OUT==ROOT or ROOT in OUT.parents or OUT.exists():parser.error('Output must be a new directory outside the code package.')
fixed={k:'1' for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS','NUMEXPR_NUM_THREADS']};fixed.update(PYTHONHASHSEED='0',PYTHONDONTWRITEBYTECODE='1',PYTHONIOENCODING='utf-8')
if os.environ.get('_LROT_F3_DETERMINISTIC_CHILD')!='1':
 child=os.environ.copy();child.update(fixed);child['_LROT_F3_DETERMINISTIC_CHILD']='1'
 child['NUMBA_CACHE_DIR']=str(OUT/'numba_cache');child['MPLCONFIGDIR']=str(OUT/'mpl_cache')
 raise SystemExit(subprocess.call([sys.executable,'-u',__file__,*sys.argv[1:]],env=child))
OUT.mkdir(parents=True,exist_ok=False)
if args.mode=='verify':
 import numpy as np
 from scipy.spatial.distance import cdist
 from sklearn.metrics import adjusted_rand_score
 z=np.load(BASE/'ordered_inputs.npz',allow_pickle=False)
 expected=json.loads((BASE/'verification.json').read_text('utf-8'))
 records=[]
 for name in ['a','b']:
  e=np.load(BASE/f'embedding_{name}.npy',allow_pickle=False)
  pred=z['B_region_labels'][cdist(e[:1000],e[1000:]).argmin(1)]
  record={'repeat':name,'accuracy':float(np.mean(pred==z['A_region_labels'])),'ARI':float(adjusted_rand_score(z['A_region_labels'],pred))}
  assert record['accuracy']==expected['first_metrics']['accuracy'] and record['ARI']==expected['first_metrics']['ARI']
  records.append(record)
 assert np.array_equal(np.load(BASE/'embedding_a.npy'),np.load(BASE/'embedding_b.npy'))
 (OUT/'verification.json').write_text(json.dumps({'scoring':records,'embedding_equal':True,'training_executed':False},indent=2)+'\n','utf-8')
 print(json.dumps(records));raise SystemExit(0)

def deny(*a,**k):raise RuntimeError('Network forbidden')
socket.create_connection=deny;socket.socket.connect=deny;socket.socket.connect_ex=deny
import numpy as np
from scipy import sparse
from scipy.spatial.distance import cdist
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import adjusted_rand_score
import scanpy as sc
import torch
import STAligner as sta
from threadpoolctl import threadpool_limits,threadpool_info
import importlib.metadata as md
threadpool_limits(1);torch.set_num_threads(1);torch.set_num_interop_threads(1)
torch.use_deterministic_algorithms(True,warn_only=False)
torch.backends.cudnn.deterministic=True;torch.backends.cudnn.benchmark=False
torch.backends.mkldnn.deterministic=True
assert torch.are_deterministic_algorithms_enabled() and torch.get_num_threads()==1 and torch.get_num_interop_threads()==1
assert os.environ.get('PYTHONHASHSEED')=='0'
print('DETERMINISM',torch.are_deterministic_algorithms_enabled(),'threads',torch.get_num_threads(),flush=True)
sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()
def js(name,obj): (OUT/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n','utf-8')
source=BASE/'scientific_case.py'
shutil.copy2(source,OUT/'original_staligner_default_epochs.py')
identity=json.loads((BASE/'input_identity.json').read_text('utf-8'))
assert sha(BASE/'ordered_inputs.npz')==identity['ordered_input_npz_sha256']
shutil.copy2(BASE/'ordered_inputs.npz',OUT/'ordered_inputs.npz')
shutil.copy2(BASE/'ordered_spots.csv',OUT/'ordered_spots.csv')
z=np.load(OUT/'ordered_inputs.npz',allow_pickle=False)
print('INPUT_KEYS',z.files,flush=True)
A={k:z['A_'+k] for k in ['coords','expr','region_labels']};B={k:z['B_'+k] for k in ['coords','expr','region_labels']}
for prefix,d in [('A',A),('B',B)]:d.update(n_spots=1000,gene_names=z[prefix+'_genes'].tolist())
module=importlib.import_module('STAligner.train_STAligner')
func=sta.train_STAligner
snap=OUT/'source_snapshot';snap.mkdir()
sourcehash={}
for name in ['train_STAligner.py','STALIGNER.py','gat_conv.py','mnn_utils.py','ST_utils.py','__init__.py']:
 f=Path(module.__file__).parent/name;shutil.copy2(f,snap/name);sourcehash[name]=sha(f)
versions={k:md.version(k) for k in ['numpy','scipy','torch','torch-geometric','STAligner','scanpy','anndata','pandas','scikit-learn']}
js('environment.json',{'python':sys.version,'packages':versions,'installed_source_sha256':sourcehash,'threads':1,'interop_threads':torch.get_num_interop_threads(),'deterministic_algorithms':torch.are_deterministic_algorithms_enabled(),'warn_only':torch.is_deterministic_algorithms_warn_only_enabled(),'cudnn_benchmark':torch.backends.cudnn.benchmark,'mkldnn_deterministic':torch.backends.mkldnn.deterministic,'PYTHONHASHSEED':os.environ.get('PYTHONHASHSEED'),'threadpools':threadpool_info(),'old_STAligner_source_hash_available':False,'nearest_neighbor_backend':'installed mnn_utils optional exact sklearn fallback; source preserved'})
js('protocol.json',{'scope':'One fresh run using original realistic-seed42 run_pair function; no DLPFC runs or parameter search','original_script_sha256':sha(source),'ordered_inputs_sha256':sha(OUT/'ordered_inputs.npz'),'preprocessing':'original function: concat -> normalize_total 1e4 -> log1p -> scanpy PCA50, unchanged defaults','graph':'original coords KNN7 query drop first then six neighbors; concatenate forward/reverse edges unchanged including duplicates','parameters':{'n_epochs':1000,'lr':.001,'verbose':False,'device':'cpu','random_seed':666,'knn_neigh':50},'default_parameters':str(inspect.signature(func)),'evaluation':'A-to-B Euclidean latent nearest neighbour; labels withheld from training','expected_old_accuracy':.934,'selection':'single prescribed seed; output retained regardless of score'})
logs=[];models=[]
def capture(epoch,loss,model,local):
 logs.append({'epoch':int(epoch),'phase':'pretrain' if epoch<500 else 'joint','loss':float(loss.detach().cpu())});models[:]=[model]
 if epoch>=500 and epoch%100==0 and 'anchor_ind' in local:
  np.savez_compressed(OUT/f'triplets_epoch_{epoch}.npz',anchor=np.asarray(local['anchor_ind']),positive=np.asarray(local['positive_ind']),negative=np.asarray(local['negative_ind']))
 if (epoch+1)%50==0:
  js('training_losses.json',logs);print('PROGRESS',epoch+1,'/1000',logs[-1]['loss'],flush=True)
tree=ast.parse(inspect.getsource(func))
class Instrument(ast.NodeTransformer):
 def visit_Subscript(self,node):
  self.generic_visit(node)
  if isinstance(node.slice,ast.Call) and isinstance(node.slice.func,ast.Attribute) and node.slice.func.attr=='keys':node.slice=ast.Call(func=ast.Name(id='list',ctx=ast.Load()),args=[node.slice],keywords=[])
  return node
 def visit_Expr(self,node):
  self.generic_visit(node)
  if isinstance(node.value,ast.Call) and isinstance(node.value.func,ast.Attribute) and node.value.func.attr=='step' and isinstance(node.value.func.value,ast.Name) and node.value.func.value.id=='optimizer':return [node,ast.parse('_capture_training(epoch,loss,model,locals())').body[0]]
  return node
tree=Instrument().visit(tree);ast.fix_missing_locations(tree)
(snap/'instrumented_train.py').write_text(ast.unparse(tree)+'\n','utf-8')
ns=dict(module.__dict__);ns.update(_capture_training=capture,tqdm=lambda iterable:iterable)
exec(compile(tree,'<recording_only_training>','exec'),ns)
trainfn=ns['train_STAligner']
def train_record(ad,*args,**kw):
 assert ad.shape==(2000,50)
 assert 'region_labels' not in ad.obs
 np.savez_compressed(OUT/'actual_training_input.npz',X=ad.X.toarray(),edgeList=ad.uns['edgeList'],obs_names=np.array(ad.obs_names,dtype=str),batch=np.array(ad.obs['batch_name'],dtype=str))
 print('TRAINING_INPUT',ad.shape,ad.uns['edgeList'].shape,flush=True)
 result=trainfn(ad,*args,**kw)
 assert len(logs)==1000 and sum(v['phase']=='joint' for v in logs)==500
 np.save(OUT/'embedding.npy',ad.obsm['STAligner']);torch.save(models[0].state_dict(),OUT/'final_model.pt');js('training_losses.json',logs)
 return result
sta.train_STAligner=train_record
expected_environment=json.loads((BASE/'environment.json').read_text('utf-8'))
assert versions==expected_environment['packages'], 'Package versions differ from pinned environment'
assert sourcehash==expected_environment['installed_source_sha256'], 'STAligner installed sources differ; do not silently substitute an upstream build'
neighbor_module=importlib.import_module('STAligner.mnn_utils')
assert not getattr(neighbor_module,'_has_hnswlib',False) and not getattr(neighbor_module,'_has_annoy',False), 'Pinned execution requires the recorded exact-neighbour fallback'
oldtree=ast.parse(source.read_text('utf-8-sig'))
node=next(x for x in oldtree.body if isinstance(x,ast.FunctionDef) and x.name=='run_pair')
exec(compile(ast.Module(body=[node],type_ignores=[]),str(source),'exec'),globals())
t=time.perf_counter()
row=run_pair('realistic-seed42',A,B)
row.update(completed_epochs=len(logs),actual_training_input_sha256=sha(OUT/'actual_training_input.npz'),embedding_file_sha256=sha(OUT/'embedding.npy'),training_scope='fresh complete single prescribed run; package-returned embedding before final weight update')
js('result.json',row);js('completed.json',{'success':True,'training_steps':1000,'seconds':time.perf_counter()-t});print('COMPLETED',json.dumps(row),flush=True)

expected_embedding=np.load(BASE/'embedding_a.npy',allow_pickle=False)
actual_embedding=np.load(OUT/'embedding.npy',allow_pickle=False)
assert np.array_equal(actual_embedding,expected_embedding), 'Fresh embedding differs from pinned deterministic result; preserve output and investigate'
js('regression_verification.json',{'passed':True,'fresh_embedding_exact_equal':True,'accuracy':row['accuracy_nn'],'ARI':row['ari_nn'],'scope':'fixed case only'})
