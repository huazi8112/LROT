"""Train the pinned DLPFC model using hash-bound, separately prepared PCA inputs."""
from pathlib import Path
import argparse,ast,hashlib,importlib,inspect,json,os,shutil,socket,sys,time
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--input-dir',required=True,type=Path)
p.add_argument('--source-dir',required=True,type=Path)
p.add_argument('--code-root',type=Path,default=Path(__file__).resolve().parents[1])
p.add_argument('--output-dir',required=True,type=Path)
a=p.parse_args();OUT=a.output_dir.resolve()
if OUT.exists():p.error('Output must be new.')
if OUT==a.code_root.resolve() or a.code_root.resolve() in OUT.parents:p.error('Output must be outside the code package.')
if os.environ.get('_LROT_DLPFC_TRAIN_CHILD')!='1':
 env=os.environ.copy();env.update(_LROT_DLPFC_TRAIN_CHILD='1',PYTHONHASHSEED='0',PYTHONDONTWRITEBYTECODE='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMBA_CACHE_DIR=str(OUT/'numba_cache'),MPLCONFIGDIR=str(OUT/'mpl_cache'))
 raise SystemExit(__import__('subprocess').call([sys.executable,'-B',__file__,*sys.argv[1:]],env=env))
OUT.mkdir(parents=True);sys.dont_write_bytecode=True;sys.path.insert(0,str(a.source_dir.resolve()))
def deny(*a,**k):raise RuntimeError('Offline training; network disabled')
socket.create_connection=deny;socket.socket.connect=deny;socket.socket.connect_ex=deny
import numpy as np
from scipy import sparse
from scipy.spatial.distance import cdist
from sklearn.metrics import adjusted_rand_score
from threadpoolctl import threadpool_limits
import anndata as an
import torch
import STAligner as sta
threadpool_limits(1);torch.set_num_threads(1);torch.set_num_interop_threads(1)
torch.use_deterministic_algorithms(True,warn_only=False)
torch.backends.cudnn.deterministic=True;torch.backends.cudnn.benchmark=False;torch.backends.mkldnn.deterministic=True
assert os.environ.get('PYTHONHASHSEED')=='0'
def sha(f):return hashlib.sha256(Path(f).read_bytes()).hexdigest()
def js(d,n,v):(d/n).write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n','utf-8')
module=importlib.import_module('STAligner.train_STAligner')
sourcehash={n:sha(Path(module.__file__).parent/n) for n in ['train_STAligner.py','STALIGNER.py','gat_conv.py','mnn_utils.py','ST_utils.py','__init__.py']}
expected=json.loads((a.code_root/'provenance/f3_staligner/environment.json').read_text('utf-8'))
assert sourcehash==expected['installed_source_sha256']
import importlib.metadata as md
versions={n:md.version(n) for n in expected['packages']};assert versions==expected['packages']
js(OUT,'environment.json',{'packages':versions,'source_sha256':sourcehash,'CPU_threads':1,'deterministic_algorithms':True,'PYTHONHASHSEED':'0','scope':'Fresh independently installed dependencies, pinned reconstructed source; no changes to original environment.'})
logs=[];models=[];current=None
def capture(epoch,loss,model):
 logs.append({'epoch':int(epoch),'phase':'pretrain' if epoch<500 else 'joint','loss':float(loss.detach().cpu())});models[:]=[model]
 if (epoch+1)%100==0:js(current,'training_losses.json',logs);print('PROGRESS',current.name,epoch+1,'/1000',flush=True)
class Instrument(ast.NodeTransformer):
 def visit_Subscript(self,node):
  self.generic_visit(node)
  if isinstance(node.slice,ast.Call) and isinstance(node.slice.func,ast.Attribute) and node.slice.func.attr=='keys':node.slice=ast.Call(func=ast.Name(id='list',ctx=ast.Load()),args=[node.slice],keywords=[])
  return node
 def visit_Expr(self,node):
  self.generic_visit(node)
  if isinstance(node.value,ast.Call) and isinstance(node.value.func,ast.Attribute) and node.value.func.attr=='step' and isinstance(node.value.func.value,ast.Name) and node.value.func.value.id=='optimizer':return [node,ast.parse('_capture(epoch,loss,model)').body[0]]
  return node
tree=Instrument().visit(ast.parse(inspect.getsource(sta.train_STAligner)));ast.fix_missing_locations(tree)
(OUT/'instrumented_train.py').write_text(ast.unparse(tree)+'\n','utf-8')
ns=dict(module.__dict__);ns.update(_capture=capture,tqdm=lambda x:x);exec(compile(tree,'<recording-only>','exec'),ns)
results=[]
for pair in ['151507-151508','151508-151509']:
 inp=a.input_dir/pair;current=OUT/pair;current.mkdir();logs.clear();models.clear()
 bind=json.loads((inp/'preparation.json').read_text('utf-8'))
 assert sha(inp/'actual_training_input.npz')==bind['actual_training_input_sha256']
 assert sha(inp/'ordered_spots.csv')==bind['ordered_spots_sha256']
 for name in ['actual_training_input.npz','ordered_spots.csv','ordered_labels.npz','input_identity.json','preparation.json']:shutil.copy2(inp/name,current/name)
 z=np.load(inp/'actual_training_input.npz',allow_pickle=False);labs=np.load(inp/'ordered_labels.npz',allow_pickle=False)
 ad=an.AnnData(X=sparse.csr_matrix(z['X']));ad.obs_names=z['obs_names'].tolist();ad.obs['batch_name']=__import__('pandas').Categorical(z['batch']);ad.uns['edgeList']=z['edgeList']
 assert ad.shape==(2000,50) and ad.uns['edgeList'].shape==(2,24000)
 start=time.perf_counter();ns['train_STAligner'](ad,n_epochs=1000,lr=.001,verbose=False,device='cpu',random_seed=666,knn_neigh=50)
 assert len(logs)==1000 and sum(r['phase']=='pretrain' for r in logs)==500
 e=ad.obsm['STAligner'];assert np.isfinite(e).all();np.save(current/'embedding.npy',e);torch.save(models[0].state_dict(),current/'final_model.pt');js(current,'training_losses.json',logs)
 nearest=cdist(e[:1000],e[1000:]).argmin(1);pred=labs['B_labels'][nearest];np.save(current/'predictions.npy',pred)
 r={'pair':pair,'accuracy':float((pred==labs['A_labels']).mean()),'ARI':float(adjusted_rand_score(labs['A_labels'],pred)),'training_steps':len(logs),'pretrain_steps':500,'joint_steps':500,'seconds':time.perf_counter()-start,'actual_training_input_sha256':sha(current/'actual_training_input.npz'),'embedding_sha256':sha(current/'embedding.npy'),'labels_withheld_from_training':True,'scope':'One prescribed deterministic single-thread fresh run; not a claim of exact historical training repeatability.'}
 js(current,'result.json',r);results.append(r);js(OUT,'results.json',results);print('COMPLETED',json.dumps(r),flush=True)
js(OUT,'completed.json',{'complete':True,'total_steps':2000,'pairs':2})
