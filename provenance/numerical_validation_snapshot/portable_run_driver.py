"""Portable derivative of the archived driver; not the original source bytes.
Requires explicit LROT_DATA_ROOT. The formal API extracts unchanged definitions
and does not execute this historical experiment driver.
"""
from pathlib import Path
import os,sys,json,hashlib,ast,inspect,time,socket,traceback,runpy,warnings,faulthandler
# Portable derivative: historical original bytes are retained in a private archive.
# The computational definitions below are unchanged; no default personal path.
_data_root = os.environ.get('LROT_DATA_ROOT')
if not _data_root:
    raise RuntimeError('Set LROT_DATA_ROOT explicitly to the authorized data project root; no default data search is permitted.')
DATA = Path(_data_root)
if not DATA.is_absolute():
    raise ValueError('LROT_DATA_ROOT must be an absolute authorized project directory.')
DATA = DATA.resolve(strict=True)
if not DATA.is_dir():
    raise ValueError('LROT_DATA_ROOT must identify a directory.')
faulthandler.dump_traceback_later(120,repeat=True)
print('BOOTSTRAP portable numerical driver',flush=True)
B=Path(__file__).resolve().parent;C=B/'code'
driver=sys.argv[1];stage=sys.argv[2] if len(sys.argv)>2 else ''
group=Path(driver).stem+('_'+stage if stage else '')
O=B/'runs'/group;O.mkdir(parents=True,exist_ok=True);E=O/'evidence';E.mkdir(exist_ok=True)
lock=O/'running.lock'
while True:
    try:
        with lock.open('x') as f:f.write(str(os.getpid()))
        break
    except FileExistsError:time.sleep(2)
os.environ.update(LROT_OUTPUT_DIR=str(O),LROT_OUT=str(O),LROT_DATA_DIR=str(DATA/'cancer_data'),LROT_DLPFC_DATA=str(DATA/'dlpfc_data'),MPLBACKEND='Agg',MPLCONFIGDIR=str(B/'mpl'),OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',MKL_NUM_THREADS='2')
(B/'numba_cache').mkdir(exist_ok=True);(B/'temp').mkdir(exist_ok=True)
os.environ.update(NUMBA_CACHE_DIR=str(B/'numba_cache'),TEMP=str(B/'temp'),TMP=str(B/'temp'))
sys.dont_write_bytecode=True;sys.path[:0]=[str(C),str(C/'experiments'),str(C/'experiments/supp')]
def deny(*a,**k):raise RuntimeError('Private-paper workflow: network forbidden')
socket.create_connection=deny;socket.socket.connect=deny;socket.socket.connect_ex=deny
import numpy as np
from threadpoolctl import threadpool_limits
threadpool_limits(2)
import lrot_core as core
import lrot_numerics_certified as cert
from lrot_numerics import plan_diagnostics
core.balanced_sinkhorn=cert.balanced_sinkhorn
actual=core.fgw_lr_solver;signature=inspect.signature(actual)
fingerprint=hashlib.sha256(b''.join((C/n).read_bytes() for n in ['lrot_core.py','lrot_numerics.py','lrot_numerics_certified.py','real_st_loader.py'])).hexdigest()
calls=[];hashes={}
class RejectedPlan(RuntimeError):pass
def arrayhash(x):
    x=np.ascontiguousarray(x);h=hashlib.sha256();h.update(str(x.shape).encode());h.update(str(x.dtype).encode());h.update(x.tobytes());return h.hexdigest()
def clean(x):
    if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [clean(v) for v in x]
    if isinstance(x,np.generic):return clean(x.item())
    if isinstance(x,np.ndarray):return clean(x.tolist()) if x.size<=5000 else {'shape':list(x.shape),'dtype':str(x.dtype),'sha256':arrayhash(x)}
    if isinstance(x,float) and not np.isfinite(x):return None
    if x is None or isinstance(x,(str,int,float,bool)):return x
    return str(x)
def context():
    out={}
    for f in inspect.stack()[2:6]:
        if str(C) not in f.filename:continue
        v={k:clean(x) for k,x in f.frame.f_locals.items() if k in ['seed','ds','s','sigma','sig','mode','name','variant','gamma','g','alpha','beta','i','sample_A','sample_B','lr_program','eps_mult','method','rotation','n_spots'] and isinstance(x,(str,int,float,bool,np.generic))}
        out[Path(f.filename).name+':'+f.function]=v
    return out
def evaluate(P,A,D):
    Q=P/P.sum(1,keepdims=True);out={'entropy':float(-(Q*np.log(np.maximum(Q,1e-300))).sum(1).mean())}
    if 'region_labels' in A and 'region_labels' in D and max(np.max(A['region_labels']),np.max(D['region_labels']))>0:
        labs=np.asarray(D['region_labels'],int);la=np.asarray(A['region_labels'],int);n=int(max(labs.max(),la.max()))+1
        pred=(Q@np.eye(n)[labs]).argmax(1)
        from sklearn.metrics import adjusted_rand_score
        out.update(soft_acc=float((pred==la).mean()),soft_ari=float(adjusted_rand_score(la,pred)),same_mass=float(P[la[:,None]==labs[None,:]].sum()/P.sum()))
    else:out['annotation_metrics_available']=False
    return out
def solver(*a,**kw):
    bound=signature.bind(*a,**kw);bound.apply_defaults();v=bound.arguments
    A=v.pop('slice_A',None);D=v.pop('slice_B',None);S=v.pop('S_lr',None)
    # Resolve actual argument names from signature rather than assuming a name.
    if A is None:
        values=list(bound.arguments.values());A,D,S=values[:3];v={k:x for k,x in list(bound.arguments.items())[3:]}
    params=dict(v);params.update(strict=False,return_diagnostics=True,verbose=False,inner_max_iter=10000,marginal_abs_tol=1e-8,marginal_rel_tol=1e-5,tol=1e-5)
    identities={'A':{k:arrayhash(A[k]) for k in ['coords','expr','region_labels']},'B':{k:arrayhash(D[k]) for k in ['coords','expr','region_labels']},'gene_names_A':list(A['gene_names']),'gene_names_B':list(D['gene_names']),'lr':arrayhash(S),'params':clean(params),'solver':fingerprint}
    key=hashlib.sha256(json.dumps(identities,sort_keys=True).encode()).hexdigest()
    if driver=='serial_runtime_audit.py':key+='__trial_'+str(len(calls)+1)
    j=E/(key+'.json');p=E/(key+'.npz')
    ctx=context();calls.append({'key':key,'context':ctx,'call':len(calls)+1});(O/'calls.json').write_text(json.dumps(calls,ensure_ascii=False,indent=2),encoding='utf-8')
    if j.exists():
        record=json.loads(j.read_text(encoding='utf-8'))
        if record['accepted']:
            with np.load(p) as z:return z['P'],z['loss'].tolist()
        raise RejectedPlan('Previously unaccepted result; no silent replacement '+key)
    # Reuse only an exactly identical, already accepted configuration in this
    # isolated revalidation. Never reuse plans for timing/scaling measurements.
    if driver not in ['lrot_scaling.py','run_runtime_benchmark.py','serial_runtime_audit.py']:
        existing=next((x for x in (B/'runs').glob('*/evidence/'+key+'.json') if x!=j),None)
        if existing is not None:
            record=json.loads(existing.read_text(encoding='utf-8'))
            if record.get('accepted') and record['input_identity']==identities:
                import shutil
                record['exact_checkpoint_reuse']=str(existing);record['driver']=driver;record['context']=ctx
                j.write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf-8');shutil.copy2(existing.with_suffix('.npz'),p)
                with np.load(p) as z:return z['P'],z['loss'].tolist()
    cert.reset_dual_cache();t=time.time();attempts=[]
    P,loss,d=actual(A,D,S,**params);attempts.append(d)
    if not d['success'] and d['stop_reason']=='outer_iteration_cap':
        params['max_iter']=1000;cert.reset_dual_cache();P,loss,d=actual(A,D,S,**params);attempts.append(d)
    core._last_compute_seconds=time.time()-t
    record={'key':key,'driver':driver,'stage':stage,'accepted':bool(d['success']),'input_identity':identities,'final_parameters':clean(params),'context':ctx,'runtime_seconds':core._last_compute_seconds,'diagnostics':d,'attempts':attempts,'metrics':evaluate(P,A,D)}
    j.write_text(json.dumps(clean(record),ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8');np.savez_compressed(p,P=P,loss=np.array(loss))
    print('AUDIT',group,len(calls),d['success'],d['stop_reason'],round(record['runtime_seconds'],2),record['metrics'],flush=True)
    if not d['success']:raise RejectedPlan('Unaccepted plan cannot enter figure/report: '+key)
    return P,loss
core.fgw_lr_solver=solver
# Input resolution stays within authorized original experiment directories.
import lrot_paths
original_find=lrot_paths.find_input
def find_input(name,*extras):
    choices=[Path(x)/name for x in extras]+[O/name]
    roots=[DATA/'LROT_code/lrot_output',DATA/'lrot_output']
    for root in roots:
        choices.extend([root/name,root/'数据_npz'/name,root/'figures_png'/name])
        choices.extend(root.rglob(name))
    for p in choices:
        if p.is_file():
            hashes[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest();return p
    raise FileNotFoundError('Authorized local input unavailable: '+name)
lrot_paths.find_input=find_input
# External algorithms retain their own implementation; capture finite feasibility and convergence logs.
pending_external_evidence=[]
def external_hook(module,name,tag):
    fn=getattr(module,name);sig=inspect.signature(fn)
    def wrapped(*a,**kw):
        info={'method':tag,'parameters':clean(kw),'context':context()};logged='return_log' in sig.parameters
        info['input_identity']=[{'expression_hash':arrayhash(ad.X),'spatial_hash':arrayhash(ad.obsm['spatial']),'gene_order':list(ad.var_names),'shape':list(ad.shape)} for ad in a[:2]]
        if logged:kw['return_log']=True
        elif 'log' in sig.parameters:kw['log']=True;logged=True
        t=time.time();result=fn(*a,**kw);compute_seconds=time.time()-t
        module._last_compute_seconds=compute_seconds
        for pending_key,pending_info,pending_P,pending_loss in pending_external_evidence:
            (E/(pending_key+'.json')).write_text(json.dumps(pending_info,indent=2),encoding='utf-8')
            np.savez_compressed(E/(pending_key+'.npz'),P=pending_P,loss=pending_loss)
        pending_external_evidence.clear()
        P,log=result if logged and isinstance(result,tuple) else (result,{})
        P=np.asarray(P);d=plan_diagnostics(P,np.ones(P.shape[0])/P.shape[0],np.ones(P.shape[1])/P.shape[1])
        info.update(seconds=compute_seconds,diagnostics=clean(d),log=clean(log),finite=bool(np.isfinite(P).all()),nonnegative=bool((P>=0).all()),mass=float(P.sum()))
        info['feasible']=info['finite'] and info['nonnegative'] and abs(info['mass']-1)<1e-6 and max(d['row_max_relative'],d['col_max_relative'])<=1e-5
        key=tag+'_'+str(len(list(E.glob(tag+'*.json')))+1)
        (E/(key+'.json')).write_text(json.dumps(info,ensure_ascii=False,indent=2),encoding='utf-8');np.savez_compressed(E/(key+'.npz'),P=P)
        if not info['feasible']:raise RuntimeError('External plan is infeasible: '+key)
        return P
    setattr(module,name,wrapped)
if driver in ['run_dlpfc_baselines.py','run_official_paste2.py','run_paste2_multiseed.py','run_paste2_multiseed_euclid.py','run_runtime_benchmark.py','external_convergence_audit.py','synthetic_staligner_reconstruction.py','serial_runtime_audit.py','staligner_full_protocol.py','staligner_default_epochs.py']:
    import paste_math,paste_math.PASTE as pm,paste2.PASTE2 as p2
    import importlib.metadata as md
    env={'implementations':{str(Path(m.__file__)):hashlib.sha256(Path(m.__file__).read_bytes()).hexdigest() for m in [pm,p2]}}
    for pkg in ['paste-math','paste2','STAligner','POT','torch','torch-geometric','anndata','scanpy','numpy','scipy']:
        try:env[pkg]=md.version(pkg)
        except md.PackageNotFoundError:env[pkg]='local import; distribution metadata unavailable'
    (O/'external_environment.json').write_text(json.dumps(env,indent=2),encoding='utf-8')
    external_hook(paste_math,'pairwise_align','PASTE');external_hook(p2,'partial_pairwise_align','PASTE2')
    low_level=pm.my_fused_gromov_wasserstein
    def low_level_record(*a,**kw):
        kw['log']=True
        result=low_level(*a,**kw);P,log=result
        loss=np.asarray(log.get('loss',[]),float);cap=int(kw.get('numItermax',200))
        key='PASTE_CG_'+str(len(list(E.glob('PASTE_CG*.json')))+1)
        rel=float(abs(loss[-1]-loss[-2])/max(abs(loss[-2]),1e-300)) if len(loss)>1 else None
        info={'parameters':clean(kw),'loss':clean(loss.tolist()),'iterations':max(0,len(loss)-1),'iteration_cap':cap,'cap_reached':len(loss)-1>=cap,'last_relative_change':rel,'last_absolute_change':float(abs(loss[-1]-loss[-2])) if len(loss)>1 else None,'objective':clean(log.get('fgw_dist')),'plan_hash':arrayhash(P),'status':'conditional gradient stopping diagnostic; not global optimality'}
        if driver=='serial_runtime_audit.py':pending_external_evidence.append((key,info,np.asarray(P),loss))
        else:
            (E/(key+'.json')).write_text(json.dumps(info,indent=2),encoding='utf-8')
            np.savez_compressed(E/(key+'.npz'),P=np.asarray(P),loss=loss)
        return result
    pm.my_fused_gromov_wasserstein=low_level_record
    p2_inner=p2.partial_fused_gromov_wasserstein;p2_signature=inspect.signature(p2_inner)
    def p2_inner_record(*a,**kw):
        kw['log']=True;bound=p2_signature.bind(*a,**kw);bound.apply_defaults();values=bound.arguments
        final={};old_profile=sys.getprofile()
        def profiler(frame,event,arg):
            if event=='return' and frame.f_code is p2_inner.__code__:
                final.update({k:clean(frame.f_locals.get(k)) for k in ['cpt','gamma','relative_delta_fval','abs_delta_fval','numItermax','stopThr','stopThr2']})
        sys.setprofile(profiler)
        try:result=p2_inner(*a,**kw)
        finally:sys.setprofile(old_profile)
        P,log=result;loss=np.asarray(log.get('loss',[]),float);cap=int(values['numItermax'])
        stopped=(final.get('gamma')==0 or (final.get('relative_delta_fval') is not None and final['relative_delta_fval']<values['stopThr']) or (final.get('abs_delta_fval') is not None and final['abs_delta_fval']<values['stopThr2']))
        key='PASTE2_CG_'+str(len(list(E.glob('PASTE2_CG*.json')))+1)
        info={'parameters':{k:clean(values[k]) for k in ['alpha','m','armijo','numItermax','tol','stopThr','stopThr2']},'loss':loss.tolist(),'iterations':len(loss),'final_internal_state':final,'cap_reached':len(loss)>=cap and not stopped,'stopping_criterion_met':bool(stopped),'plan_hash':arrayhash(P),'objective':clean(log.get('partial_fgw_cost')),'status':'package stopping diagnostics; not global optimality'}
        if driver=='serial_runtime_audit.py':pending_external_evidence.append((key,info,np.asarray(P),loss))
        else:
            (E/(key+'.json')).write_text(json.dumps(info,indent=2),encoding='utf-8');np.savez_compressed(E/(key+'.npz'),P=np.asarray(P),loss=loss)
        return result
    p2.partial_fused_gromov_wasserstein=p2_inner_record
    if driver in ['run_dlpfc_baselines.py','synthetic_staligner_reconstruction.py','staligner_full_protocol.py','staligner_default_epochs.py']:
        import torch;torch.set_num_threads(2)
        import STAligner as sta
        train=sta.train_STAligner
        def train_record(ad,*a,**kw):
            t=time.time();r=train(ad,*a,**kw);x=np.asarray(ad.obsm['STAligner'])
            assert np.isfinite(x).all()
            key='STAligner_'+str(len(list(E.glob('STAligner*.json')))+1)
            np.savez_compressed(E/(key+'.npz'),embedding=x)
            (E/(key+'.json')).write_text(json.dumps({'method':'STAligner','parameters':clean(kw),'seconds':time.time()-t,'embedding_shape':list(x.shape),'embedding_hash':arrayhash(x),'finite':True,'status':'fixed training schedule; not OT convergence certification'},indent=2),encoding='utf-8')
            return r
        sta.train_STAligner=train_record
path=C/'experiments'/driver
faulthandler.cancel_dump_traceback_later()
sys.argv=[str(path)]+([stage] if stage else [])
started=time.time()
returned=[]
def audit_main_locals(values):
    selected={k:clean(v) for k,v in values.items() if k in ['raw','results','rows','all_rows','all_results','out','paste_curves','lrot_curves','summary','grid','ref','base_rows','data','curves','res','stats','sums','acc','ent','ablation_results','dlpfc','breast']}
    (O/'main_local_results.json').write_text(json.dumps(selected,ensure_ascii=False,indent=2),encoding='utf-8')
def capture_return(fn):
    def wrapped(*a,**kw):
        result=fn(*a,**kw)
        returned.append({'function':fn.__name__,'context':context(),'result':clean(result)})
        (O/'functional_return_records.json').write_text(json.dumps(returned,ensure_ascii=False,indent=2),encoding='utf-8')
        return result
    return wrapped
try:
    tree=ast.parse(path.read_text(encoding='utf-8-sig'))
    for n in tree.body:
        if isinstance(n,ast.FunctionDef) and n.name=='main':
            n.body.append(ast.Expr(value=ast.Call(func=ast.Name(id='audit_main_locals',ctx=ast.Load()),args=[ast.Call(func=ast.Name(id='locals',ctx=ast.Load()),args=[],keywords=[])],keywords=[])))
        if isinstance(n,ast.FunctionDef) and n.name in ['run_pair','analyze_pair','run_config','solve_full','run_method','evaluate','run','run_lrot','run_official_paste2','run_one','run_seed','run_dlpfc_calibration','run_breast_calibration']:
            n.decorator_list.insert(0,ast.Name(id='capture_return',ctx=ast.Load()))
    globals_={'__name__':'__main__','__file__':str(path),'__package__':None,'capture_return':capture_return,'audit_main_locals':audit_main_locals}
    exec(compile(ast.fix_missing_locations(tree),str(path),'exec'),globals_)
    outputs={k:clean(globals_[k]) for k in ['results','rows','all_rows','all_results','ablation_results'] if k in globals_}
    (O/'global_results.json').write_text(json.dumps(outputs,ensure_ascii=False,indent=2),encoding='utf-8')
    status={'driver':driver,'stage':stage,'completed':True,'seconds':time.time()-started,'solver_calls':len(calls)}
except SystemExit as ex:
    status={'driver':driver,'stage':stage,'completed':ex.code in (None,0),'seconds':time.time()-started,'solver_calls':len(calls),'explicit_script_exit':ex.code}
    if not status['completed']:status['error']='Script exited '+str(ex.code)
except Exception as ex:
    traceback.print_exc();status={'driver':driver,'stage':stage,'completed':False,'seconds':time.time()-started,'solver_calls':len(calls),'error':repr(ex)}
(O/'driver_status.json').write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf-8')
(O/'resolved_input_hashes.json').write_text(json.dumps(hashes,ensure_ascii=False,indent=2),encoding='utf-8')
print('DRIVER_STATUS',json.dumps(status),flush=True)
lock.unlink(missing_ok=True)
sys.exit(0 if status['completed'] else 1)
